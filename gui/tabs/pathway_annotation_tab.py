import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext, filedialog
import logging
import os
import time
import tempfile
import pandas as pd
import numpy as np
import threading
from datetime import datetime

# Import shared components
from gui.shared.base_tab import BaseTab, _setup_global_styles
from gui.shared.utils import split_pathways
from gui.shared.column_assignment import show_column_assignment_dialog

from main_script.metabolite_pathways_annotator import MetabolitePathwayMapper
from main_script.metabolite_ID_annotator import MetaboliteIDAnnotator
from main_script.metabolite_pathway_network import is_disease_pathway_global
import main_script.metabolite_pathway_network as metabolite_pathway_network
import main_script.lipid_search as lipid_search

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PathwayAnnotationTab(BaseTab):
    """Pathway Annotation Tab - Placeholder for pathway annotation
    
    This tab will handle pathway annotation of metabolites.
    """
    
    def __init__(self, parent, data_manager):
        """Initialize Pathway Annotation tab"""
        super().__init__(parent, data_manager)
        
        # Setup global styles (runs only once)
        _setup_global_styles()
        
        # Get root window for dialogs
        self.root = parent.winfo_toplevel()
        
        # Thread tracking to prevent duplicate operations
        self.pathway_calculation_thread = None
        self.pathway_annotated_excel_path = None
        
        # Timer tracking variables
        self.pathway_start_time = None
        self.pathway_current_step = 0
        self.pathway_total_steps = 100
        
        # Lipid search results storage
        self.lipid_search_results = None
        self.verified_assignments = {}
        self.pvalue_verification_complete = False  # Gating: Track if p-value column verified with validation
        self.verified_pvalue_col = None  # Store verified p-value column name
        self.log2fc_verification_complete = False  # Gating: Track if log2FC column verified
        self.verified_log2fc_col = None  # Store verified log2FC column name
        
        # Create UI
        self.setup_ui()
        print("[OK] Pathway Annotation Tab initialized")

    def setup_ui(self):
        """Setup the UI for pathway annotation tab"""
        self.create_pathway_annotation_interface(self.frame)
    
    def create_pathway_annotation_interface(self, parent):
        """Create the 2-column pathway annotation interface with progress log underneath"""
        # Create main frame without scrollbar (only left column will have scrollbar)
        main_canvas_frame = tk.Frame(parent, bg='#f0f0f0')
        main_canvas_frame.pack(fill='both', expand=True)
        
        # Create main frame directly without canvas scrolling
        scrollable_content = tk.Frame(main_canvas_frame, bg='#f0f0f0')
        scrollable_content.pack(fill='both', expand=True)
        
        # Create main frame inside scrollable content
        main_frame = tk.Frame(scrollable_content, bg='#f0f0f0')
        main_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Create main content frame for 2 columns
        content_frame = tk.Frame(main_frame, bg='#f0f0f0')
        content_frame.pack(fill='both', expand=True)
            
        # Configure grid: 2 columns with 1:1 ratio
        content_frame.grid_columnconfigure(0, weight=1, uniform="columns")  # Left column
        content_frame.grid_columnconfigure(1, weight=1, uniform="columns")  # Right column - equal weight
        # Reduce minimum height for better scrollbar usability
        content_frame.grid_rowconfigure(0, weight=1, minsize=400)
            
        # ============ LEFT COLUMN: Controls ============
        left_frame = tk.LabelFrame(content_frame, 
                    text="⚙️ Controls & Settings", 
                    font=('Arial', 12, 'bold'),
                    bg='#f0f0f0', 
                    fg='#2c3e50')
        # Make left frame expandable in all directions like tabs 1 & 2
        left_frame.grid(row=0, column=0, sticky='nsew', padx=(0, 5))
            
        # Add canvas and scrollbar for left frame content only (matching tabs 1 & 2)
        left_canvas = tk.Canvas(left_frame, bg='#f0f0f0')
        left_scrollbar_y = ttk.Scrollbar(left_frame, orient="vertical", command=left_canvas.yview)
        # Remove horizontal scrollbar to match tabs 1 and 2
        # left_scrollbar_x = ttk.Scrollbar(left_frame, orient="horizontal", command=left_canvas.xview)
        left_scrollable_frame = tk.Frame(left_canvas, bg='#f0f0f0')
        
        left_scrollable_frame.bind(
            "<Configure>",
            lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all"))
        )
        
        canvas_window = left_canvas.create_window((0, 0), window=left_scrollable_frame, anchor="nw")
        left_canvas.configure(yscrollcommand=left_scrollbar_y.set)
        
        # Configure canvas and scrollbars with pack layout (like tabs 1 and 2)
        left_canvas.pack(side="left", fill="both", padx=15, pady=10, expand=True)
        left_scrollbar_y.pack(side="right", fill="y")
        
        # Remove horizontal scrollbar to match tabs 1 and 2 behavior
        # left_scrollbar_x.grid(row=1, column=0, sticky='ew')
        
        # Configure grid weights (keeping these for the frame structure)
        left_frame.grid_rowconfigure(0, weight=1)
        left_frame.grid_columnconfigure(0, weight=1)
        
        # Add canvas width configuration function (like tabs 1 and 2)
        def configure_scroll_region(event):
            left_canvas.configure(scrollregion=left_canvas.bbox("all"))
            canvas_width = event.width
            left_canvas.itemconfig(canvas_window, width=canvas_width)
        
        left_canvas.bind('<Configure>', configure_scroll_region)
        
        def _on_left_mousewheel(event):
            left_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        left_canvas.bind("<MouseWheel>", _on_left_mousewheel)
        left_scrollable_frame.bind("<MouseWheel>", _on_left_mousewheel)
        
        # Data Mode Selection (Metabolite/Lipid) - Moved to top
        mode_frame = tk.LabelFrame(left_scrollable_frame, 
                                  text="🔬 Data Mode", 
                                  bg='#f0f0f0',
                                  font=('Arial', 10, 'bold'))
        mode_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(mode_frame, 
                text="Select data type:", 
                bg='#f0f0f0',
                font=('Arial', 9)).pack(anchor='w', padx=5, pady=(5, 0))
        
        mode_button_frame = tk.Frame(mode_frame, bg='#f0f0f0')
        mode_button_frame.pack(fill='x', padx=20, pady=5)
        
        self.pathway_data_mode = tk.StringVar(value='metabolite')
        tk.Radiobutton(mode_button_frame, text='Metabolite', variable=self.pathway_data_mode,
                       value='metabolite', bg='#f0f0f0', command=self.on_pathway_mode_change).pack(side='left', padx=5)
        tk.Radiobutton(mode_button_frame, text='Lipid', variable=self.pathway_data_mode,
                       value='lipid', bg='#f0f0f0', command=self.on_pathway_mode_change).pack(side='left', padx=5)
        
        self.pathway_mode_info_label = tk.Label(mode_frame,
                text="ℹ️ Metabolite mode: Feature ID required. Optional: HMDB, KEGG, PubChem, ChEBI, CAS, SMILES, InChI, InChIKey.",
                bg='#f0f0f0',
                font=('Arial', 8, 'italic'),
                fg='#7f8c8d')
        self.pathway_mode_info_label.pack(anchor='w', padx=5, pady=(0, 5))
        
        # ============ Project Setup (REQUIRED before annotation) ============
        project_frame = tk.LabelFrame(left_scrollable_frame,
                                     text="📁 Project Setup (Required)",
                                     font=('Arial', 10, 'bold'),
                                     bg='#fff3cd',  # Light yellow background to draw attention
                                     fg='#856404')
        project_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(project_frame,
                text="⚠️ Configure project before running annotation",
                bg='#fff3cd',
                font=('Arial', 9, 'bold'),
                fg='#856404').pack(anchor='w', padx=5, pady=(5, 2))
        
        # Working folder selection
        tk.Label(project_frame,
                text="Working Folder:",
                bg='#fff3cd',
                font=('Arial', 9)).pack(anchor='w', padx=5, pady=(5, 0))
        
        folder_frame = tk.Frame(project_frame, bg='#fff3cd')
        folder_frame.pack(fill='x', padx=5, pady=2)
        
        self.project_working_folder = tk.StringVar()
        tk.Entry(folder_frame,
                textvariable=self.project_working_folder,
                font=('Arial', 9),
                state='readonly').pack(side='left', fill='x', expand=True, padx=(0, 5))
        tk.Button(folder_frame,
                 text="📁 Browse",
                 command=self.browse_working_folder,
                 bg='#856404',
                 fg='white',
                 font=('Arial', 8, 'bold')).pack(side='right')
        
        # Project name input
        tk.Label(project_frame,
                text="Project Name:",
                bg='#fff3cd',
                font=('Arial', 9)).pack(anchor='w', padx=5, pady=(5, 0))
        
        project_name_frame = tk.Frame(project_frame, bg='#fff3cd')
        project_name_frame.pack(fill='x', padx=5, pady=2)
        
        self.project_name = tk.StringVar()
        self.project_name.trace('w', lambda *args: self._update_project_status())
        self.project_name_entry = tk.Entry(project_name_frame,
                                           textvariable=self.project_name,
                                           font=('Arial', 9))
        self.project_name_entry.pack(side='left', fill='x', expand=True)
        
        # Project setup status indicator
        self.project_status_label = tk.Label(project_frame,
                text="❌ Project not configured",
                bg='#fff3cd',
                font=('Arial', 8, 'bold'),
                fg='#d32f2f')
        self.project_status_label.pack(anchor='w', padx=5, pady=(2, 5))
        
        tk.Label(project_frame,
                text="ℹ️ Results will auto-save as: [project_name]_pathway_annotated_[timestamp].xlsx",
                bg='#fff3cd',
                font=('Arial', 7, 'italic'),
                fg='#856404').pack(anchor='w', padx=5, pady=(0, 5))
        
        # File selection frame
        file_frame = tk.LabelFrame(left_scrollable_frame, 
                                  text="File Selection", 
                                  bg='#f0f0f0')
        file_frame.pack(fill='x', padx=10, pady=10)
        
        # Input file selection
        tk.Label(file_frame, 
                text="Input Excel file (with IDs):", 
                bg='#f0f0f0').pack(anchor='w', padx=5, pady=2)
        
        input_frame = tk.Frame(file_frame, bg='#f0f0f0')
        input_frame.pack(fill='x', padx=5, pady=2)
        
        self.pathway_input_file = tk.StringVar()
        tk.Entry(input_frame, 
                textvariable=self.pathway_input_file, 
                font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
        tk.Button(input_frame, 
                 text="Browse", 
                 command=lambda: self.browse_file(self.pathway_input_file, "Excel file with metabolite IDs"),
                 bg='#3498db', 
                 fg='white', 
                 font=('Arial', 8, 'bold')).pack(side='right')
        
        # Sheet selection
        sheet_frame = tk.Frame(file_frame, bg='#f0f0f0')
        sheet_frame.pack(fill='x', padx=5, pady=2)
        
        tk.Label(sheet_frame, 
                text="Select Sheet:", 
                bg='#f0f0f0',
                font=('Arial', 9)).pack(side='left')
        
        self.pathway_sheet_var = tk.StringVar()
        self.pathway_sheet_combo = ttk.Combobox(sheet_frame, 
                                               textvariable=self.pathway_sheet_var,
                                               font=('Arial', 8),  # Smaller font
                                               state='readonly',
                                               width=15)  # Reduced from 30 to 15
        self.pathway_sheet_combo.pack(side='left', padx=(5, 0), fill='x', expand=True)
        
        tk.Button(sheet_frame, 
                 text="Auto-Detect", 
                 command=lambda: self.auto_detect_sheets(self.pathway_input_file, self.pathway_sheet_combo, self.pathway_sheet_var),
                 bg='#e67e22', 
                 fg='white', 
                 font=('Arial', 8, 'bold')).pack(side='right', padx=(5, 0))
        
        # Verify Columns button
        verify_frame = tk.Frame(file_frame, bg='#f0f0f0')
        verify_frame.pack(fill='x', padx=5, pady=5)
        
        self.verify_columns_button = tk.Button(
            verify_frame,
            text="🔍 Verify Columns",
            command=self.verify_pathway_columns,
            bg='#9b59b6',
            fg='white',
            font=('Arial', 9, 'bold'),
            relief='raised',
            borderwidth=2
        )
        self.verify_columns_button.pack(fill='x', pady=2)
        
        # Verification status indicator
        self.verify_status_label = tk.Label(verify_frame,
                text="⚠️ Column verification REQUIRED before running",
                bg='#ffe6e6',
                font=('Arial', 8, 'bold'),
                fg='#d32f2f')
        self.verify_status_label.pack(anchor='w', pady=(2, 0))
        
        tk.Label(verify_frame,
                text="ℹ️ Verify Feature ID and optional ID columns",
                bg='#f0f0f0',
                font=('Arial', 8, 'italic'),
                fg='#7f8c8d').pack(anchor='w', pady=(2, 0))
        
        # Organism Selection
        organism_frame = tk.LabelFrame(left_scrollable_frame, 
                                      text="🧬 Organism Selection", 
                                      bg='#f0f0f0',
                                      font=('Arial', 10, 'bold'))
        organism_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(organism_frame, 
                text="Select organism for pathway databases:", 
                bg='#f0f0f0',
                font=('Arial', 9)).pack(anchor='w', padx=5, pady=(5, 5))
        
        self.pathway_organism_var = tk.StringVar(value="Homo sapiens")
        organism_options = ["Homo sapiens", "Rattus norvegicus", "Mus musculus"]
        
        # Create 2 organisms per row
        org_row1 = tk.Frame(organism_frame, bg='#f0f0f0')
        org_row1.pack(fill='x', padx=20, pady=2)
        tk.Radiobutton(org_row1, text="Homo sapiens", variable=self.pathway_organism_var,
                      value="Homo sapiens", bg='#f0f0f0', font=('Arial', 9)).pack(side='left', padx=(0, 20))
        tk.Radiobutton(org_row1, text="Rattus norvegicus", variable=self.pathway_organism_var,
                      value="Rattus norvegicus", bg='#f0f0f0', font=('Arial', 9)).pack(side='left')
        tk.Radiobutton(org_row1, text="Mus musculus", variable=self.pathway_organism_var,
                      value="Mus musculus", bg='#f0f0f0', font=('Arial', 9)).pack(side='left')

        # ============ Parallel Processing Settings ============
        import multiprocessing
        max_workers = multiprocessing.cpu_count()
        
        workers_frame = tk.LabelFrame(left_scrollable_frame,
                                      text="⚡ Parallel Processing",
                                      font=('Arial', 10, 'bold'),
                                      bg='#f0f0f0',
                                      fg='#2c3e50')
        workers_frame.pack(fill='x', padx=10, pady=5)
        
        # System info
        system_info_label = tk.Label(workers_frame,
                                    text=f"System has {max_workers} CPU cores available",
                                    bg='#f0f0f0',
                                    font=('Arial', 8, 'italic'),
                                    fg='#7f8c8d')
        system_info_label.pack(anchor='w', padx=5, pady=2)
        
        # Workers setting
        workers_inner_frame = tk.Frame(workers_frame, bg='#f0f0f0')
        workers_inner_frame.pack(fill='x', padx=5, pady=5)
        
        tk.Label(workers_inner_frame,
                text="Parallel workers:",
                bg='#f0f0f0',
                font=('Arial', 9)).pack(side='left')
        
        self.pathway_workers = tk.StringVar(value=str(min(6, max_workers)))
        workers_spinbox = tk.Spinbox(workers_inner_frame,
                                    from_=1, to=max_workers,
                                    textvariable=self.pathway_workers,
                                    width=5,
                                    font=('Arial', 9))
        workers_spinbox.pack(side='left', padx=(10, 5))
        
        tk.Button(workers_inner_frame,
                 text="Auto",
                 command=lambda: self.pathway_workers.set(str(max_workers)),
                 bg='#9b59b6',
                 fg='white',
                 font=('Arial', 8, 'bold'),
                 width=6).pack(side='left', padx=2)
        
        # Note: Pathway filters removed - filtering now done in Network tab via data table selections
        # Initialize pathway filter variables with default values (used by backend)
        # Default: Z score activated +1, Z score inhibited -1
        self.pathway_pos_zscore = tk.StringVar(value="1.0")
        self.pathway_neg_zscore = tk.StringVar(value="-1.0")
        self.pathway_neutral_z_band = tk.StringVar(value="0.5")
        # NOTE: self.pathway_min_metabolites is defined in Fisher ORA settings section (line ~345)
        # Do NOT redefine it here to avoid overwriting the user's selection

        # Note: Upstream Regulator and Disease Filters moved to Pathway Network tab under Network Components
        # This keeps pathway annotation focused on pathway filtering only
        
        # ============ RIGHT COLUMN: Progress & Status ============
        right_frame = tk.LabelFrame(content_frame, 
                                    text="📋 Progress & Status", 
                                    font=('Arial', 12, 'bold'),
                                    bg='#f0f0f0', 
                                    fg='#2c3e50')
        # Make sure right frame expands with the window and maintains padding from edge
        right_frame.grid(row=0, column=1, sticky='nsew', padx=(5, 10))  # Added right padding

        # Action Buttons - at the top
        action_btn_frame = tk.Frame(right_frame, bg='#f0f0f0')
        action_btn_frame.pack(fill='x', padx=10, pady=(10, 5))
        
        # Run Pathway Annotation button
        self.pathway_annotation_button = tk.Button(action_btn_frame, 
                                                  text="🔄 Run Annotation", 
                                                  command=self.run_pathway_annotation_step,
                                                  bg='#95a5a6', 
                                                  fg='white', 
                                                  font=('Arial', 9, 'bold'),
                                                  relief='raised',
                                                  borderwidth=2,
                                                  state='disabled')  # Initially disabled - requires verification
        self.pathway_annotation_button.pack(side='left', padx=(0, 5), expand=True, fill='x')
        
        # Stop button
        self.pathway_stop_button = tk.Button(action_btn_frame, 
                                            text="⏹️ Stop", 
                                            command=self.stop_pathway_annotation,
                                            bg='#e74c3c', 
                                            fg='white', 
                                            font=('Arial', 9, 'bold'),
                                            relief='raised',
                                            borderwidth=2,
                                            state='disabled')  # Initially disabled
        self.pathway_stop_button.pack(side='left', padx=5)
        
        # Export Results button
        self.export_pathway_results_button = tk.Button(action_btn_frame, 
                                                      text="💾 Export", 
                                                      command=self.export_pathway_annotation_results,
                                                      bg='#16a085', 
                                                      fg='white', 
                                                      font=('Arial', 9, 'bold'),
                                                      state='disabled')  # Initially disabled
        self.export_pathway_results_button.pack(side='left', padx=(5, 0), expand=True, fill='x')

        # Progress controls in a horizontal layout
        progress_controls_frame = tk.Frame(right_frame, bg='#f0f0f0')
        progress_controls_frame.pack(fill='x', padx=10, pady=10)
        
        # Progress bar (left side)
        progress_left_frame = tk.Frame(progress_controls_frame, bg='#f0f0f0')
        progress_left_frame.pack(side='left', fill='x', expand=True)
        
        self.pathway_progress_bar = ttk.Progressbar(progress_left_frame,
                                                   mode='determinate',
                                                   length=300)
        self.pathway_progress_bar.pack(side='left', padx=(0, 10))
        
        # Progress percentage and status (right side)
        progress_right_frame = tk.Frame(progress_controls_frame, bg='#f0f0f0')
        progress_right_frame.pack(side='right')
        
        self.pathway_progress_percentage = tk.Label(progress_right_frame,
                                                   text="0%",
                                                   bg='#f0f0f0',
                                                   font=('Arial', 9, 'bold'),
                                                   fg='#3498db')
        self.pathway_progress_percentage.pack(side='left', padx=(0, 15))
        
        # Timer and ETA (middle section)
        self.pathway_timer_label = tk.Label(progress_right_frame,
                                           text="Time: 00:00",
                                           bg='#f0f0f0',
                                           font=('Arial', 9),
                                           fg='#7f8c8d')
        self.pathway_timer_label.pack(side='left', padx=(0, 10))
        
        self.pathway_eta_label = tk.Label(progress_right_frame,
                                         text="ETA: --:--",
                                         bg='#f0f0f0',
                                         font=('Arial', 9),
                                         fg='#7f8c8d')
        self.pathway_eta_label.pack(side='left', padx=(0, 15))
        
        # Progress text area with scrollbar
        self.pathway_progress_text = scrolledtext.ScrolledText(
            right_frame,
            font=('Courier', 9),
            wrap=tk.WORD
        )
        self.pathway_progress_text.pack(fill='both', expand=True, padx=10, pady=(5, 10))
        
        # Initialize data storage
        self.pathway_original_metabolites_data = None
        self.pathway_filtered_metabolites_data = None
        self.pathway_original_pathways_data = None
        self.pathway_filtered_pathways_data = None
        self.pathway_selections = {}  # Track pathway selections
        
        # Storage for upstream regulators and diseases (for filtering)
        self.original_upstream_data = None  # Stores {(name, type): [metabolites]}
        self.original_disease_data = None   # Stores {disease_name: {pathways, metabolite_count, pvalue, etc.}}
        
        # Initialize metabolite filter variables (for legacy mode compatibility)
        if not hasattr(self, 'metabolite_pvalue_threshold'):
            self.metabolite_pvalue_threshold = tk.StringVar(value="0.05")
        if not hasattr(self, 'metabolite_pos_fc'):
            self.metabolite_pos_fc = tk.StringVar(value="0")
        if not hasattr(self, 'metabolite_neg_fc'):
            self.metabolite_neg_fc = tk.StringVar(value="0")
    
    def update_log_safe(self, message):
        """Thread-safe logging to the progress text widget"""
        try:
            if hasattr(self, 'pathway_progress_text') and self.pathway_progress_text:
                self.root.after(0, lambda: (
                    self.pathway_progress_text.insert(tk.END, f"{message}\n"),
                    self.pathway_progress_text.see(tk.END)
                ))
        except Exception as e:
            print(f"Failed to update log: {str(e)}")
    
    # def handle_pathway_calculation_error(self, error_msg):
    #     """Handle pathway calculation error in main thread"""
    #     print(f"Error calculating pathway statistics: {error_msg}")
    #     import traceback
    #     traceback.print_exc()
    #     self.update_pathway_progress_with_percentage(0, f"❌ Error: {error_msg}")
    #     self.pathway_start_time = None  # Stop timer updates
    #     self.pathway_eta_label.config(text="ETA: --:--")
    #     messagebox.showerror("Error", f"Failed to calculate pathway statistics: {error_msg}")
    
    def _save_results_after_calculation(self, log_message):
        """Auto-save results to Excel after pathway statistics calculation"""
        try:
            if self.pathway_filtered_metabolites_data is None:
                log_message(f"⚠️ No data to save (annotation may have failed)\n")
                return
            
            # ========================================
            # AUTO-SAVE WITH PROJECT NAME
            # ========================================
            project_folder = self.project_working_folder.get().strip()
            project_name = self.project_name.get().strip()
            
            if not project_folder or not project_name:
                log_message(f"⚠️ Project setup incomplete - cannot auto-save\n")
                return
            
            # Generate filename with timestamp
            from datetime import datetime
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            filename = f"{project_name}_pathway_annotated_{timestamp}.xlsx"
            save_path = os.path.join(project_folder, filename)
            
            log_message(f"\n💾 Auto-saving results...\n")
            log_message(f"   Project: {project_name}\n")
            log_message(f"   Folder: {project_folder}\n")
            log_message(f"   File: {filename}\n")
            
            # Save to Excel
            annotated_df = self.pathway_filtered_metabolites_data.copy()
            annotated_df.to_excel(save_path, index=False, sheet_name='Annotated_Data')
            
            # Store the saved file path
            self.pathway_annotated_excel_path = save_path
            
            log_message(f"✅ Results auto-saved successfully!\n")
            log_message(f"   Full path: {save_path}\n\n")
            
            # ========================================
            # AUTO-LOAD TO NETWORK ANALYSIS TAB
            # ========================================
            log_message(f"🔄 Loading results into Network Analysis tab...\n")
            
            # Get reference to network tab through parent
            network_tab = getattr(self, 'network_tab_instance', None)
            
            if network_tab and hasattr(network_tab, 'load_annotated_data'):
                # Load from saved file path so auto-load and manual Browse import
                # follow the exact same parsing/normalization path.
                network_tab.load_annotated_data(save_path)
                log_message(f"✅ Data loaded into Network Analysis tab\n")
                
                # Switch to Network Analysis subtab (with small delay for UI update)
                self.root.after(500, self._switch_to_network_tab)
                log_message(f"🔄 Switching to Network Analysis tab...\n")
            else:
                log_message(f"⚠️ Network Analysis tab not found for auto-loading\n")
                
        except Exception as e:
            log_message(f"⚠️ Failed to auto-save/load results: {str(e)}\n")
            import traceback
            traceback.print_exc()
    
    def _switch_to_network_tab(self):
        """Switch to Network Analysis subtab.

        Prefer calling the parent container's helper if available.  The hardcoded
        index used previously became invalid after adding the Comparative tab,
        which pushed Network to position 1 (and Multi-Omics to 2).  Using the
        parent's `switch_to_network_tab` method guarantees we end up on the
        correct tab regardless of order changes.
        """
        try:
            # Try parent object method first (PathwayAnalysisParentTab)
            if hasattr(self, 'frame') and hasattr(self.frame, 'master'):
                parent = self.frame.master
                if hasattr(parent, 'switch_to_network_tab'):
                    parent.switch_to_network_tab()
                    print("[DEBUG] Switched to Network Analysis tab via parent helper")
                    return
                # Fallback: old logic using direct notebook indices
                if isinstance(parent, ttk.Notebook):
                    # previous assumption was index 2, but tabs may have been
                    # reordered; network should now be index 1 in current layout
                    # (annotation=0, network=1, multiomics=2, comparative=3)
                    try:
                        parent.select(1)
                        print("[DEBUG] Switched to Network Analysis tab (fallback index=1)")
                        return
                    except Exception:
                        pass
            print("[DEBUG] Could not switch to Network Analysis tab (no parent helper)")
        except Exception as e:
            print(f"[ERROR] Failed to switch tabs: {e}")
    
    def _populate_network_tab_after_annotation(self, log_message=None):
        """Safely load annotated results into the Network Analysis tab after async work."""
        try:
            network_tab = getattr(self, 'network_tab_instance', None)
            excel_path = getattr(self, 'pathway_annotated_excel_path', None)
            annotated_df = getattr(self, 'pathway_filtered_metabolites_data', None)
            if not (network_tab and hasattr(network_tab, 'load_annotated_data')):
                if log_message:
                    log_message("⚠️ Network Analysis tab unavailable for auto-load\n")
                return
            if not excel_path or annotated_df is None:
                if log_message:
                    log_message("⚠️ Annotated data not ready for Network Analysis tab\n")
                return
            # Use saved file path to ensure consistency with manual Browse import.
            network_tab.load_annotated_data(excel_path)
            if log_message:
                log_message("✅ Network Analysis tab refreshed with latest annotation\n")
            # Keep UI in sync by switching tabs (non-blocking)
            self.root.after(0, self._switch_to_network_tab)
        except Exception as exc:
            if log_message:
                log_message(f"⚠️ Failed to populate Network Analysis tab: {exc}\n")

    def _update_tree_table_after_calculation(self, log_message=None):
        """Refresh the pathway results tree after stats complete."""
        try:
            self.update_pathway_table()
            if log_message:
                log_message("✅ Pathway results table updated\n")
        except Exception as exc:
            if log_message:
                log_message(f"⚠️ Failed to update pathway table: {exc}\n")

    def update_pathway_progress_with_percentage(self, percentage, message):
        """Update pathway progress with percentage, timer, and message"""
        self.pathway_current_step = percentage
        self.pathway_progress_bar.config(value=percentage)
        self.pathway_progress_percentage.config(text=f"{percentage}%")
        
        # Update progress text
        self.pathway_progress_text.insert(tk.END, f"[{percentage}%] {message}\n")
        self.pathway_progress_text.see(tk.END)
        
        self.root.update_idletasks()
    
    def update_pathway_table(self):
        """Update the pathway table display"""
        # Clear existing items
        for item in self.pathway_tree.get_children():
            self.pathway_tree.delete(item)
        
        if self.pathway_filtered_pathways_data is None:
            self.pathway_stats_label.config(text="No pathway data")
            return
        
        # Import the centralized pathway filtering function
        try:
            from main_script.metabolite_pathway_network import should_filter_pathway
        except ImportError:
            # Fallback if import fails - define locally (minimal version)
            def should_filter_pathway(pathway_name):
                """Fallback filtering - basic disease detection"""
                if not pathway_name:
                    return True
                pathway_lower = pathway_name.lower()
                disease_keywords = ["disease", "disorder", "syndrome", "emia", "uria", "osis"]
                for keyword in disease_keywords:
                    if keyword in pathway_lower:
                        return True
                return False
        
        # Reset selection tracking
        # FIXED: Explicitly mark all pathways as selected by default (changed from implicit default-True to explicit True)
        self.pathway_selections = {}
        
        # Update column header based on FDR status and method
        selected_method = getattr(self, 'pathway_analysis_method', tk.StringVar(value="Fisher ORA")).get()
        fdr_enabled = getattr(self, 'pathway_fdr_enabled', True)  # Default to True if not set
        if selected_method == "IWPA":
            pvalue_header = 'P-Adj'
        elif fdr_enabled:
            pvalue_header = 'P-Adj'
        else:
            pvalue_header = 'P-Value'
        self.pathway_tree.heading('P-Value', text=pvalue_header)
        
        # Log table update
        pathway_count = len([k for k in self.pathway_filtered_pathways_data.keys() if k != '_validation_report'])
        self.pathway_progress_text.insert(tk.END, f"Updating pathway table with {pathway_count} pathways...\n")
        self.pathway_progress_text.insert(tk.END, f"📊 P-value column showing: {pvalue_header}\n")
        self.pathway_progress_text.insert(tk.END, f"🔬 Analysis method: {selected_method}\n")
        self.pathway_progress_text.see(tk.END)
        
        filtered_count = 0
        for pathway_name, stats in self.pathway_filtered_pathways_data.items():
            # Skip the validation report entry (it's a string, not a pathway stats dict)
            # Also skip any other special diagnostic keys (they start with an underscore)
            if isinstance(pathway_name, str) and pathway_name.startswith('_'):
                continue
            # Skip disease pathways and overly general pathways using centralized filter
            if should_filter_pathway(pathway_name):
                filtered_count += 1
                continue
            
            pvalue = stats.get('combined_pvalue', 0)
            # Original simpler z-score extraction
            zscore = stats.get('z_score', 0.0)
            status = stats.get('status', 'No Change')
            n_metabolites = stats.get('n_metabolites', 0)
            
            # FIXED: Explicitly mark pathway as selected by default
            self.pathway_selections[pathway_name] = True
            select_symbol = "☑️"
            
            # Format display values
            pvalue_display = f"{pvalue:.6f}"
            zscore_display = f"{zscore:.3f}" if zscore is not None else "0.000"
            
            # Insert into table (columns: Select, Pathway, Status, Z-Score, P-Value, # Metabolites)
            # Insert full pathway name (Treeview will handle visual truncation via column width)
            item_id = self.pathway_tree.insert('', 'end', values=(
                select_symbol, pathway_name, status, zscore_display, pvalue_display, n_metabolites
            ))
        
        # Update stats - count only DISPLAYED pathways (not filtered out by disease/threshold)
        # Get all displayed pathway names from the tree
        displayed_pathways = []
        for item_id in self.pathway_tree.get_children():
            values = self.pathway_tree.item(item_id)['values']
            if values and len(values) >= 2:
                pathway_name = values[1]  # Pathway name is second column
                displayed_pathways.append(pathway_name)
        
        # Count status from displayed pathways only
        total_pathways = len(displayed_pathways)
        activated = 0
        inhibited = 0
        no_change = 0
        
        for pw_name in displayed_pathways:
            if pw_name in self.pathway_filtered_pathways_data:
                status = self.pathway_filtered_pathways_data[pw_name].get('status', 'No Change')
                if status == 'Activated':
                    activated += 1
                elif status == 'Inhibited':
                    inhibited += 1
                else:
                    no_change += 1
        
        stats_text = f"Total: {total_pathways} | Activated: {activated} | Inhibited: {inhibited} | No Change: {no_change}"
        self.pathway_stats_label.config(text=stats_text)
        
        # Log completion
        self.pathway_progress_text.insert(tk.END, f"Pathway table updated: {stats_text}\n")
        self.pathway_progress_text.see(tk.END)
        
        # If diagnostic info about excluded pathways exists, show a compact hint
        try:
            excluded_preview = self.pathway_filtered_pathways_data.get('_excluded_debug')
            if excluded_preview:
                sample_names = ', '.join([e.get('pathway', '') for e in excluded_preview[:3] if isinstance(e, dict)])
                self.pathway_progress_text.insert(tk.END, f"ℹ️ Some pathways were excluded by thresholds (e.g., FDR/min_mets). Sample: {sample_names}\n")
                self.pathway_progress_text.see(tk.END)
        except Exception:
            pass
    
    def update_enzymes_table(self):
        """Update the enzymes table display"""
        # Clear existing items
        for item in self.enzymes_tree.get_children():
            self.enzymes_tree.delete(item)
        
        if self.pathway_filtered_metabolites_data is None:
            self.enzymes_stats_label.config(text="No enzyme data")
            return
        
        df = self.pathway_filtered_metabolites_data
        
        # Debug: Print available columns
        print("\n=== DEBUG: Enzyme Table Update ===")
        print(f"Available columns: {df.columns.tolist()}")
        
        # Reset selection tracking
        self.enzyme_selections = {}
        
        # Collect enzymes data - try multiple column name variations
        enzymes_dict = {}  # {enzyme_name: [metabolites]}
        
        # Possible column names for enzymes
        # Prefer gene symbol column for enzymes if available
        enzyme_column_names = ['Enzyme_Gene_name', 'Enzymes', 'Enzyme', 'Enzyme_Name', 'enzymes']
        enzyme_col = None
        for col_name in enzyme_column_names:
            if col_name in df.columns:
                enzyme_col = col_name
                print(f"Found enzyme column: {col_name}")
                break
        
        if enzyme_col is None:
            print(f"WARNING: No enzyme column found. Available columns: {df.columns.tolist()}")
            self.enzymes_stats_label.config(text="No enzyme column found in data")
            return
        
        for _, row in df.iterrows():
            metabolite_name = row.get('Name', 'Unknown')
            
            # Process Enzymes - try different separators
            enzymes = row.get(enzyme_col, '')
            
            if pd.notna(enzymes) and str(enzymes).strip():
                enzymes_str = str(enzymes)
                print(f"Processing enzymes for {metabolite_name}: {enzymes_str[:100]}")
                
                # Try pipe separator first, then comma
                if '|' in enzymes_str:
                    enzyme_list = [e.strip() for e in enzymes_str.split('|') if e.strip()]
                elif ';' in enzymes_str:
                    enzyme_list = [e.strip() for e in enzymes_str.split(';') if e.strip()]
                elif ',' in enzymes_str:
                    enzyme_list = [e.strip() for e in enzymes_str.split(',') if e.strip()]
                else:
                    enzyme_list = [enzymes_str.strip()] if enzymes_str.strip() else []
                
                print(f"  Found {len(enzyme_list)} enzymes: {enzyme_list[:3]}")
                
                for enzyme in enzyme_list:
                    if enzyme not in enzymes_dict:
                        enzymes_dict[enzyme] = []
                    enzymes_dict[enzyme].append(metabolite_name)
        
        print(f"Total unique enzymes found: {len(enzymes_dict)}")
        
        # Populate table with enzymes
        for enzyme_name, metabolites in sorted(enzymes_dict.items()):
            metabolite_list = ', '.join(sorted(set(metabolites))[:5])  # Show first 5
            if len(metabolites) > 5:
                metabolite_list += f'... (+{len(metabolites)-5} more)'
            
            # Default to selected
            key = f"Enzyme_{enzyme_name}"
            self.enzyme_selections[key] = True
            select_symbol = "☑️"
            
            self.enzymes_tree.insert('', 'end', values=(
                select_symbol, enzyme_name, '', metabolite_list, len(set(metabolites))
            ))
        
        # Update stats
        total_enzymes = len(enzymes_dict)
        stats_text = f"Total Enzymes: {total_enzymes}"
        self.enzymes_stats_label.config(text=stats_text)
        
        # Log completion
        self.pathway_progress_text.insert(tk.END, f"Enzymes table updated: {stats_text}\n")
        self.pathway_progress_text.see(tk.END)
    
    def update_transporters_table(self):
        """Update the transporters table display"""
        # Clear existing items
        for item in self.transporters_tree.get_children():
            self.transporters_tree.delete(item)
        
        if self.pathway_filtered_metabolites_data is None:
            self.transporters_stats_label.config(text="No transporter data")
            return
        
        df = self.pathway_filtered_metabolites_data
        
        # Debug: Print available columns
        print("\n=== DEBUG: Transporter Table Update ===")
        print(f"Available columns: {df.columns.tolist()}")
        
        # Reset selection tracking
        self.transporter_selections = {}
        
        # Collect transporters data - try multiple column name variations
        transporters_dict = {}  # {transporter_name: [metabolites]}
        
        # Possible column names for transporters
        # Prefer gene symbol column for transporters if available
        transporter_column_names = ['Transporter_Gene_name', 'Transporters', 'Transporter', 'Transporter_Name', 'transporters']
        transporter_col = None
        for col_name in transporter_column_names:
            if col_name in df.columns:
                transporter_col = col_name
                print(f"Found transporter column: {col_name}")
                break
        
        if transporter_col is None:
            print(f"WARNING: No transporter column found. Available columns: {df.columns.tolist()}")
            self.transporters_stats_label.config(text="No transporter column found in data")
            return
        
        for _, row in df.iterrows():
            metabolite_name = row.get('Name', 'Unknown')
            
            # Process Transporters - try different separators
            transporters = row.get(transporter_col, '')
            
            if pd.notna(transporters) and str(transporters).strip():
                transporters_str = str(transporters)
                print(f"Processing transporters for {metabolite_name}: {transporters_str[:100]}")
                
                # Try pipe separator first, then comma
                if '|' in transporters_str:
                    transporter_list = [t.strip() for t in transporters_str.split('|') if t.strip()]
                elif ';' in transporters_str:
                    transporter_list = [t.strip() for t in transporters_str.split(';') if t.strip()]
                elif ',' in transporters_str:
                    transporter_list = [t.strip() for t in transporters_str.split(',') if t.strip()]
                else:
                    transporter_list = [transporters_str.strip()] if transporters_str.strip() else []
                
                print(f"  Found {len(transporter_list)} transporters: {transporter_list[:3]}")
                
                for transporter in transporter_list:
                    if transporter not in transporters_dict:
                        transporters_dict[transporter] = []
                    transporters_dict[transporter].append(metabolite_name)
        
        print(f"Total unique transporters found: {len(transporters_dict)}")
        
        # Populate table with transporters
        for transporter_name, metabolites in sorted(transporters_dict.items()):
            metabolite_list = ', '.join(sorted(set(metabolites))[:5])  # Show first 5
            if len(metabolites) > 5:
                metabolite_list += f'... (+{len(metabolites)-5} more)'
            
            # Default to selected
            key = f"Transporter_{transporter_name}"
            self.transporter_selections[key] = True
            select_symbol = "☑️"
            
            self.transporters_tree.insert('', 'end', values=(
                select_symbol, transporter_name, '', metabolite_list, len(set(metabolites))
            ))
        
        # Update stats
        total_transporters = len(transporters_dict)
        stats_text = f"Total Transporters: {total_transporters}"
        self.transporters_stats_label.config(text=stats_text)
        
        # Log completion
        self.pathway_progress_text.insert(tk.END, f"Transporters table updated: {stats_text}\n")
        self.pathway_progress_text.see(tk.END)
    
    def apply_network_component_filters(self):
        """Apply global metabolite p-value filtering before updating all tables."""
        try:
            comp_thr = float(self.component_pvalue_threshold.get())
        except Exception:
            comp_thr = None
        
        if self.pathway_filtered_metabolites_data is not None and comp_thr is not None:
            # Filter to only significant metabolites globally
            original_count = len(self.pathway_filtered_metabolites_data)
            self.pathway_filtered_metabolites_data = self.pathway_filtered_metabolites_data[
                self.pathway_filtered_metabolites_data['pvalue'] < comp_thr  # Strict < for significance
            ].copy()
            filtered_count = len(self.pathway_filtered_metabolites_data)
            self.pathway_progress_text.insert(tk.END, f"Applied global metabolite filter: {original_count} → {filtered_count} metabolites (p < {comp_thr})\n")
            self.pathway_progress_text.see(tk.END)
        
        # Now update all tables with the filtered data
        self.update_metabolite_table()
        self.update_upstream_table()
        self.update_diseases_table()
        # Note: Pathway table doesn't need re-filtering as it's based on Fisher ORA results
            
    def update_diseases_table(self):
        """Update the diseases table display"""
        # Clear existing items
        for item in self.disease_tree.get_children():
            self.disease_tree.delete(item)
        
        if self.pathway_filtered_metabolites_data is None:
            self.diseases_stats_label.config(text="No disease data")
            return
        
        df = self.pathway_filtered_metabolites_data
        # Apply component significance threshold if available
        try:
            comp_thr = float(self.component_pvalue_threshold.get()) if hasattr(self, 'component_pvalue_threshold') else None
        except Exception:
            comp_thr = None
        if comp_thr is not None and 'pvalue' in df.columns:
            df = df[df['pvalue'] <= comp_thr]  # <-- FILTERING HAPPENS HERE
        
        # Preserve existing selections if available
        prev_disease_sel = getattr(self, 'disease_selections', {}).copy()
        new_disease_sel = {}
        
        # Collect diseases data from filtered df
        diseases_dict = {}  # {disease_name: [metabolites]}
        
        print(f"DEBUG: Looking for disease data. Available columns: {df.columns.tolist()}")
        
        for _, row in df.iterrows():
            metabolite_name = row.get('Name', 'Unknown')
            
            # Process Associated Diseases
            if 'Associated_Diseases' in df.columns:
                diseases = row.get('Associated_Diseases', '')
                #print(f"DEBUG: Found Associated_Diseases column for metabolite {metabolite_name}: {diseases}")
                
                if pd.notna(diseases) and str(diseases).strip():
                    disease_list = [d.strip() for d in str(diseases).split('|') if d.strip()]
                    
                    for disease in disease_list:
                        if disease not in diseases_dict:
                            diseases_dict[disease] = []
                        diseases_dict[disease].append(metabolite_name)
            else:
                print("DEBUG: Associated_Diseases column not found in data")
        
        print(f"DEBUG: Found {len(diseases_dict)} diseases total")
        
        # Store original disease data for filtering
        self.original_disease_data = diseases_dict.copy()
        
        # Get filter value
        try:
            min_metabolites = int(self.disease_min_metabolites.get())
        except (ValueError, AttributeError):
            min_metabolites = 1
        
        # Populate table with filtering
        displayed_count = 0
        filtered_out = 0
        
        for disease_name, metabolites in sorted(diseases_dict.items()):
            num_metabolites = len(set(metabolites))
            
            # Apply minimum metabolites filter
            if num_metabolites < min_metabolites:
                filtered_out += 1
                continue
            
            # Preserve selection state if previously toggled
            selected_state = prev_disease_sel.get(disease_name, True)
            new_disease_sel[disease_name] = selected_state
            select_symbol = "☑️" if selected_state else "☐"
            
            # FIXED: Insert values matching column order (Select, Disease, # Metabolites)
            self.disease_tree.insert('', 'end', values=(
                select_symbol, 
                disease_name, 
                num_metabolites
            ))
            displayed_count += 1
        
        # Commit updated selections
        self.disease_selections = new_disease_sel

        # Update stats
        stats_text = f"Showing: {displayed_count} diseases"
        if filtered_out > 0:
            stats_text += f" | Filtered: {filtered_out}"
        self.diseases_stats_label.config(text=stats_text)
        
        # Log completion
        self.pathway_progress_text.insert(tk.END, f"Diseases table updated: {stats_text}\n")
        self.pathway_progress_text.see(tk.END)
    
    def update_upstream_table(self):
        """Update the upstream regulators table (combining enzymes and transporters)"""
        # Clear existing items
        for item in self.upstream_tree.get_children():
            self.upstream_tree.delete(item)
        
        if self.pathway_filtered_metabolites_data is None:
            self.upstream_stats_label.config(text="No upstream regulator data")
            return
        
        df = self.pathway_filtered_metabolites_data
        # Apply component significance threshold if available
        try:
            comp_thr = float(self.component_pvalue_threshold.get()) if hasattr(self, 'component_pvalue_threshold') else None
        except Exception:
            comp_thr = None
        if comp_thr is not None and 'pvalue' in df.columns:
            df = df[df['pvalue'] <= comp_thr]
        
        # Debug: Print available columns
        print("\n=== DEBUG: Upstream Regulator Table Update ===")
        print(f"Available columns: {df.columns.tolist()}")
        
        # Preserve existing selections if available
        prev_upstream_sel = getattr(self, 'upstream_selections', {}).copy()
        new_upstream_sel = {}
            
        # Collect upstream regulators data - {(name, type): [metabolites]}
        upstream_dict = {}
        
        # Process ENZYMES (prefer gene symbol if available)
        enzyme_column_names = ['Enzyme_Gene_name', 'Enzymes', 'Enzyme', 'Enzyme_Name', 'enzymes']
        enzyme_col = None
        for col_name in enzyme_column_names:
            if col_name in df.columns:
                enzyme_col = col_name
                print(f"Found enzyme column: {col_name}")
                break
        
        if enzyme_col is None:
            print(f"DEBUG: No enzyme column found. Available columns: {df.columns.tolist()}")
        
        if enzyme_col:
            for _, row in df.iterrows():
                metabolite_name = row.get('Name', 'Unknown')
                enzymes = row.get(enzyme_col, '')
                
                if pd.notna(enzymes) and str(enzymes).strip() and str(enzymes).lower() != 'nan':
                    enzymes_str = str(enzymes)
                    # Try different separators
                    if '|' in enzymes_str:
                        enzyme_list = [e.strip() for e in enzymes_str.split('|') if e.strip()]
                    elif ';' in enzymes_str:
                        enzyme_list = [e.strip() for e in enzymes_str.split(';') if e.strip()]
                    elif ',' in enzymes_str:
                        enzyme_list = [e.strip() for e in enzymes_str.split(',') if e.strip()]
                    else:
                        enzyme_list = [enzymes_str.strip()] if enzymes_str.strip() else []
                    
                    for enzyme in enzyme_list:
                        key = (enzyme, 'Enzyme')
                        if key not in upstream_dict:
                            upstream_dict[key] = []
                        upstream_dict[key].append(metabolite_name)
            
        # Process TRANSPORTERS (prefer gene symbol if available)
        transporter_column_names = ['Transporter_Gene_name', 'Transporters', 'Transporter', 'Transporter_Name', 'transporters']
        transporter_col = None
        for col_name in transporter_column_names:
            if col_name in df.columns:
                transporter_col = col_name
                print(f"Found transporter column: {col_name}")
                break
        
        if transporter_col is None:
            print(f"DEBUG: No transporter column found. Available columns: {df.columns.tolist()}")
        
        if transporter_col:
            for _, row in df.iterrows():
                metabolite_name = row.get('Name', 'Unknown')
                transporters = row.get(transporter_col, '')
                
                if pd.notna(transporters) and str(transporters).strip() and str(transporters).lower() != 'nan':
                    transporters_str = str(transporters)
                    # Try different separators
                    if '|' in transporters_str:
                        transporter_list = [t.strip() for t in transporters_str.split('|') if t.strip()]
                    elif ';' in transporters_str:
                        transporter_list = [t.strip() for t in transporters_str.split(';') if t.strip()]
                    elif ',' in transporters_str:
                        transporter_list = [t.strip() for t in transporters_str.split(',') if t.strip()]
                    else:
                        transporter_list = [transporters_str.strip()] if transporters_str.strip() else []
                    
                    for transporter in transporter_list:
                        key = (transporter, 'Transporter')
                        if key not in upstream_dict:
                            upstream_dict[key] = []
                        upstream_dict[key].append(metabolite_name)
        
        print(f"Total unique upstream regulators found: {len(upstream_dict)}")
        
        # Store original upstream data for filtering
        self.original_upstream_data = upstream_dict.copy()
        
        # Populate table with upstream regulators
        for (name, reg_type), metabolites in sorted(upstream_dict.items()):
            # Preserve selection state if previously toggled
            key = f"{reg_type}_{name}"
            selected_state = prev_upstream_sel.get(key, True)
            new_upstream_sel[key] = selected_state
            select_symbol = "☑️" if selected_state else "☐"
            
            # FIXED: Insert values matching column order (Select, Name, Type, # Metabolites)
            self.upstream_tree.insert('', 'end', values=(
                select_symbol, 
                name, 
                reg_type, 
                len(set(metabolites))
            ))
        
        # Commit updated selections
        self.upstream_selections = new_upstream_sel

        # Update stats
        total_upstream = len(upstream_dict)
        stats_text = f"Total Upstream Regulators: {total_upstream}"
        self.upstream_stats_label.config(text=stats_text)
        
        # Log completion
        self.pathway_progress_text.insert(tk.END, f"Upstream regulators table updated: {stats_text}\n")
        self.pathway_progress_text.see(tk.END)
        
    def apply_pathway_filter(self):
        """Apply pathway filtering based on network tab parameters and update all tables"""
        try:
            if self.pathway_original_pathways_data is None:
                messagebox.showwarning("Warning", "No pathway data to filter. Please run pathway annotation first.")
                return
            
            # Get filter parameters from network tab
            pos_zscore = float(self.network_zscore_activation.get())
            neg_zscore = float(self.network_zscore_inhibition.get())
            min_metabolites = int(self.network_min_metabolites.get())
            max_pvalue = float(self.network_p_threshold.get())
            status_filter = self.network_status_filter.get()
            
            # Log filter application
            filter_msg = f"Applying pathway filters: Z-score Activation ≥{pos_zscore}, Inhibition ≤{neg_zscore}, Min {min_metabolites} metabolites, P-value <{max_pvalue}"
            if status_filter != "All":
                filter_msg += f", Status: {status_filter}"
            self._log_network(filter_msg + "\n")
            
            # CRITICAL: Always start from original data to avoid compounding filters
            source_data = self.pathway_original_pathways_data.copy()
            original_count = len([k for k in source_data.keys() if k != '_validation_report'])
            
            # Debug: Log sample pathway data structure and p-value distribution
            if original_count > 0:
                sample_name = next((k for k in source_data.keys() if k != '_validation_report'), None)
                if sample_name:
                    sample_data = source_data[sample_name]
                    print(f"DEBUG: Sample pathway '{sample_name}' keys: {list(sample_data.keys())}")
                    print(f"  z_score={sample_data.get('z_score')}, adjusted_pvalue={sample_data.get('adjusted_pvalue')}, n_metabolites={sample_data.get('n_metabolites')}, status={sample_data.get('status')}")
                
                # Analyze p-value distribution across all pathways
                pvalues = []
                for name, data in source_data.items():
                    if name == '_validation_report':
                        continue
                    pval = data.get('adjusted_pvalue', 
                           data.get('fdr_qvalue',
                           data.get('combined_pvalue', None)))
                    if pval is not None:
                        pvalues.append(pval)
                
                if pvalues:
                    pvalues.sort()
                    print(f"DEBUG: P-value distribution: min={min(pvalues):.6f}, max={max(pvalues):.6f}, median={pvalues[len(pvalues)//2]:.6f}")
                    print(f"DEBUG: P-values <= 0.05: {sum(1 for p in pvalues if p <= 0.05)}/{len(pvalues)}")
                    print(f"DEBUG: P-values <= 0.10: {sum(1 for p in pvalues if p <= 0.10)}/{len(pvalues)}")
                    print(f"DEBUG: P-values <= 1.00: {sum(1 for p in pvalues if p <= 1.00)}/{len(pvalues)}")
                    
                    # Log to GUI as well
                    self._log_network(f"📊 P-value analysis: {len(pvalues)} pathways, range [{min(pvalues):.6f}, {max(pvalues):.6f}]\n")
                    self._log_network(f"   ≤0.05: {sum(1 for p in pvalues if p <= 0.05)}, ≤0.10: {sum(1 for p in pvalues if p <= 0.10)}, ≤1.00: {sum(1 for p in pvalues if p <= 1.00)}\n")
            
            # Filter pathways based on criteria
            filtered_pathways = {}
            debug_filter_results = {'passed': 0, 'failed_zscore': 0, 'failed_pvalue': 0, 'failed_min_met': 0, 'failed_status': 0}
            
            for pathway_name, pathway_data in source_data.items():
                if pathway_name == '_validation_report':
                    continue
                
                # Apply filters - handle both Fisher ORA and Legacy method keys
                zscore = pathway_data.get('z_score', 0)
                # Try adjusted_pvalue first (standard), fall back to fdr_qvalue (Fisher ORA), then combined_pvalue
                pvalue = pathway_data.get('adjusted_pvalue', 
                          pathway_data.get('fdr_qvalue',
                          pathway_data.get('combined_pvalue', 1.0)))
                num_metabolites = pathway_data.get('n_metabolites', 0)
                status = pathway_data.get('status', 'No Change')
                
                # Check all filter criteria
                passes_zscore = True
                if status == 'Activated':
                    passes_zscore = zscore >= pos_zscore
                elif status == 'Inhibited':
                    passes_zscore = zscore <= neg_zscore
                # No Change pathways always pass z-score filter
                
                passes_pvalue = pvalue < max_pvalue  # Changed from < to <= to include exact threshold
                passes_min_metabolites = num_metabolites >= min_metabolites
                # Support combined status option from UI (e.g., "Activated + Inhibited")
                if status_filter == "All":
                    passes_status = True
                elif status_filter == "Activated + Inhibited":
                    passes_status = status in ("Activated", "Inhibited")
                else:
                    passes_status = (status == status_filter)
                
                # Debug tracking
                if not passes_zscore:
                    debug_filter_results['failed_zscore'] += 1
                if not passes_pvalue:
                    debug_filter_results['failed_pvalue'] += 1
                if not passes_min_metabolites:
                    debug_filter_results['failed_min_met'] += 1
                if not passes_status:
                    debug_filter_results['failed_status'] += 1
                
                if passes_zscore and passes_pvalue and passes_min_metabolites and passes_status:
                    filtered_pathways[pathway_name] = pathway_data
                    debug_filter_results['passed'] += 1
            
            # Log debug results
            print(f"DEBUG Filter Results: {debug_filter_results}")
            self._log_network(f"  Filter debug: Passed={debug_filter_results['passed']}, Failed(zscore={debug_filter_results['failed_zscore']}, pvalue={debug_filter_results['failed_pvalue']}, min_met={debug_filter_results['failed_min_met']}, status={debug_filter_results['failed_status']})\n")
            
            # Update the filtered pathways data
            self.pathway_filtered_pathways_data = filtered_pathways
            
            # Update pathway table
            self.update_pathway_table()
            
            # Filter metabolites to only those in current pathways
            if self.pathway_filtered_metabolites_data is not None and len(filtered_pathways) > 0:
                current_pathways = set(filtered_pathways.keys())
                
                pathway_cols = ['HMDB_Pathways', 'PathBank_Pathways', 'SMP_Pathways', 'WikiPathways', 'Metabolika_Pathways', 'All_Pathways']
                present_pathway_cols = [col for col in pathway_cols if col in self.pathway_filtered_metabolites_data.columns]
                
                def metabolite_has_current_pathway(row):
                    """Check if metabolite has any of the currently filtered pathways"""
                    for col in present_pathway_cols:
                        if pd.notna(row[col]):
                            if isinstance(row[col], str):
                                pathways = split_pathways(row[col])
                            elif isinstance(row[col], list):
                                pathways = row[col]
                            else:
                                continue
                            
                            if any(pathway in current_pathways for pathway in pathways):
                                return True
                    return False
                
                original_metabolite_count = len(self.pathway_filtered_metabolites_data)
                metabolite_mask = self.pathway_filtered_metabolites_data.apply(metabolite_has_current_pathway, axis=1)
                kept_df = self.pathway_filtered_metabolites_data[metabolite_mask].copy()
                removed_df = self.pathway_filtered_metabolites_data[~metabolite_mask].copy()
                self.pathway_filtered_metabolites_data = kept_df
                
                filtered_metabolite_count = len(self.pathway_filtered_metabolites_data)
                removed_count = original_metabolite_count - filtered_metabolite_count
                self._log_network(f"  • Metabolites: {original_metabolite_count} → {filtered_metabolite_count} (removed {removed_count})\n")
                # Additional debug on metabolite filtering logic
                self._log_network(f"    - Pathway columns used: {', '.join(present_pathway_cols) if present_pathway_cols else 'None found'}\n")
                self._log_network(f"    - Kept with ≥1 selected pathway: {filtered_metabolite_count}\n")
                self._log_network(f"    - Removed with 0 selected pathways: {removed_count}\n")
                # Show up to 5 examples of removed metabolite names for transparency
                try:
                    if removed_count > 0 and 'Name' in removed_df.columns:
                        sample_removed = removed_df['Name'].dropna().astype(str).head(5).tolist()
                        if sample_removed:
                            self._log_network(f"    - Examples removed (no selected pathway match): {', '.join(sample_removed)}\n")
                except Exception:
                    pass
                
                # Update metabolite table
                self.update_metabolite_table()
                
                # Update upstream and diseases tables
                try:
                    self.update_upstream_table()
                except Exception as e:
                    print(f"Error updating upstream table: {e}")
                
                try:
                    self.update_diseases_table()
                except Exception as e:
                    print(f"Error updating diseases table: {e}")
            
            # Report results
            new_count = len(filtered_pathways)
            removed_count = original_count - new_count
            
            activated = len([s for s in filtered_pathways.values() if s['status'] == 'Activated'])
            inhibited = len([s for s in filtered_pathways.values() if s['status'] == 'Inhibited'])
            no_change = len([s for s in filtered_pathways.values() if s['status'] == 'No Change'])
            
            filter_results = f"✅ Filter Results: {original_count} → {new_count} pathways (removed {removed_count})"
            filter_details = f"  • Activated: {activated}, Inhibited: {inhibited}, No Change: {no_change}"
            
            self._log_network(f"{filter_results}\n{filter_details}\n\n")
            
        except Exception as e:
            error_msg = f"Pathway filter error: {str(e)}"
            self._log_network(f"❌ ERROR: {error_msg}\n")
            messagebox.showerror("Filter Error", error_msg)
            import traceback
            traceback.print_exc()
    
    def reset_pathway_filters(self):
        """Reset pathway filters to original data and default values"""
        try:
            if self.pathway_original_pathways_data is None:
                messagebox.showwarning("Warning", "No original pathway data available")
                return
            
            self._log_network("🔄 Resetting pathway filters to original data...\n")
            
            # Reset filter values to defaults
            self.network_min_metabolites.set("2")
            self.network_p_threshold.set("1.0")  # Changed from 0.05 to 1.0
            self.network_zscore_activation.set("2.0")
            self.network_zscore_inhibition.set("-2.0")
            self.network_status_filter.set("All")
            
            # Restore original pathway data
            self.pathway_filtered_pathways_data = self.pathway_original_pathways_data.copy()
            
            # Restore original metabolite data
            if self.pathway_original_metabolites_data is not None:
                self.pathway_filtered_metabolites_data = self.pathway_original_metabolites_data.copy()
            
            # Update all tables
            self.update_pathway_table()
            self.update_metabolite_table()
            
            try:
                self.update_upstream_table()
            except Exception as e:
                print(f"Error updating upstream table: {e}")
            
            try:
                self.update_diseases_table()
            except Exception as e:
                print(f"Error updating diseases table: {e}")
            
            # Report reset
            pathway_count = len([k for k in self.pathway_filtered_pathways_data.keys() if k != '_validation_report'])
            metabolite_count = len(self.pathway_filtered_metabolites_data) if self.pathway_filtered_metabolites_data is not None else 0
            
            self._log_network(f"✅ Filters reset: {pathway_count} pathways, {metabolite_count} metabolites restored\n\n")
            
        except Exception as e:
            error_msg = f"Reset filter error: {str(e)}"
            self._log_network(f"❌ ERROR: {error_msg}\n")
            messagebox.showerror("Reset Error", error_msg)
    
    def show_pathway_statistics(self):
        """Show statistics about pathway p-values and other metrics"""
        if self.pathway_original_pathways_data is None:
            messagebox.showwarning("No Data", "No pathway data available. Please run pathway analysis first.")
            return
        
        # Analyze pathway data
        pathways = [name for name in self.pathway_original_pathways_data.keys() if name != '_validation_report']
        if not pathways:
            messagebox.showwarning("No Pathways", "No pathways found in the data.")
            return
        
        # Collect statistics
        pvalues = []
        zscores = []
        metabolites = []
        statuses = {'Activated': 0, 'Inhibited': 0, 'No Change': 0}
        
        for name in pathways:
            data = self.pathway_original_pathways_data[name]
            pval = data.get('adjusted_pvalue', 
                   data.get('fdr_qvalue',
                   data.get('combined_pvalue', None)))
            if pval is not None:
                pvalues.append(pval)
            
            zscore = data.get('z_score', 0)
            zscores.append(zscore)
            
            n_met = data.get('n_metabolites', 0)
            metabolites.append(n_met)
            
            status = data.get('status', 'No Change')
            if status in statuses:
                statuses[status] += 1
        
        # Create statistics message
        stats_msg = f"📊 Pathway Statistics ({len(pathways)} total pathways)\n\n"
        
        if pvalues:
            pvalues.sort()
            stats_msg += f"P-values:\n"
            stats_msg += f"  Range: {min(pvalues):.6f} - {max(pvalues):.6f}\n"
            stats_msg += f"  Median: {pvalues[len(pvalues)//2]:.6f}\n"
            stats_msg += f"  ≤ 0.05: {sum(1 for p in pvalues if p <= 0.05)} pathways\n"
            stats_msg += f"  ≤ 0.10: {sum(1 for p in pvalues if p <= 0.10)} pathways\n"
            stats_msg += f"  ≤ 0.20: {sum(1 for p in pvalues if p <= 0.20)} pathways\n"
            stats_msg += f"  ≤ 1.00: {sum(1 for p in pvalues if p <= 1.00)} pathways\n\n"
        
        if zscores:
            stats_msg += f"Z-scores:\n"
            stats_msg += f"  Range: {min(zscores):.2f} - {max(zscores):.2f}\n"
            stats_msg += f"  |Z| ≥ 2.0: {sum(1 for z in zscores if abs(z) >= 2.0)} pathways\n\n"
        
        if metabolites:
            stats_msg += f"Metabolites per pathway:\n"
            stats_msg += f"  Range: {min(metabolites)} - {max(metabolites)}\n"
            stats_msg += f"  ≥ 2 metabolites: {sum(1 for m in metabolites if m >= 2)} pathways\n\n"
        
        stats_msg += f"Status distribution:\n"
        for status, count in statuses.items():
            stats_msg += f"  {status}: {count} pathways\n"
        
        stats_msg += f"\n💡 Tip: If all pathways disappear when filtering, try increasing the P-value threshold (currently {self.network_p_threshold.get()})."
        
        # Show in a scrollable message box
        self.show_scrollable_message("Pathway Statistics", stats_msg)
    
    def show_scrollable_message(self, title, message):
        """Show a scrollable message dialog"""
        dialog = tk.Toplevel(self.root)
        dialog.title(title)
        dialog.geometry("600x400")
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Text widget with scrollbar
        frame = tk.Frame(dialog)
        frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        text = tk.Text(frame, wrap=tk.WORD, font=('Courier', 10))
        scrollbar = ttk.Scrollbar(frame, command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        
        text.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')
        
        text.insert('1.0', message)
        text.config(state='disabled')
        
        # OK button
        tk.Button(dialog, text="OK", command=dialog.destroy).pack(pady=10)
    
    def save_pathway_annotation_results(self):
        """Save the pathway annotation results to Excel in Pathway_annotation folder with auto-generated filename"""
        try:
            if self.pathway_filtered_metabolites_data is None:
                raise Exception("No metabolite data available")
            
            if self.pathway_filtered_pathways_data is None:
                raise Exception("No pathway data available")
            
            # Get output folder from project folder
            if hasattr(self, 'project_working_folder') and self.project_working_folder.get():
                output_folder = self.project_working_folder.get()
            else:
                # Fallback to current directory + Pathway_annotation
                output_folder = os.path.join(os.getcwd(), "Pathway_annotation")
            
            # Create folder if doesn't exist
            os.makedirs(output_folder, exist_ok=True)
            
            # Store the output folder for export function to access
            if not hasattr(self, 'pathways_networks_dir'):
                self.pathways_networks_dir = output_folder
            
            # Auto-generate filename from sheet name with timestamp
            sheet_name = getattr(self, 'loaded_sheet_name', 'Metabolites')
            # Clean sheet name for filename (remove invalid characters)
            clean_sheet_name = "".join(c for c in sheet_name if c.isalnum() or c in (' ', '_', '-')).strip()
            clean_sheet_name = clean_sheet_name.replace(' ', '_')
            
            # Add timestamp to prevent overwriting
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            filename = f"{clean_sheet_name}_pathways_{timestamp}.xlsx"
            output_file = os.path.join(output_folder, filename)
            
            print(f"DEBUG: Auto-generated output file: {output_file}")
            
            # Store the full path for network generation
            if not hasattr(self, 'pathway_output_file'):
                self.pathway_output_file = tk.StringVar()
            self.pathway_output_file.set(output_file)
            
            print(f"📁 Saving results to: {output_folder}")
            self.update_log_safe(f"📁 Saving results to: {output_folder}")
            self.update_log_safe(f"📄 Filename: {filename}")
            
            # Prepare metabolite data for export
            metabolite_data = self.pathway_filtered_metabolites_data.copy()
            print(f"DEBUG: Metabolite data shape: {metabolite_data.shape}")
            
            # Add pathway statistics - include ALL pathways (no selection filtering)
            # Also include role-scoring columns when present and de-duplicate pathway names
            import re
            def canonical_key(name: str) -> str:
                return re.sub(r"[^a-z0-9]+", "", str(name).lower())

            def display_name(name: str) -> str:
                n = str(name).strip().strip(" .;:,\t")
                if not n:
                    return ''
                small = {"and", "or", "of", "to", "in", "on", "the", "for"}
                tokens = re.split(r"(\s+|-)", n)
                rebuilt = []
                for t in tokens:
                    if t in {'-', ' ', '\t'} or t.isspace():
                        rebuilt.append(t)
                    else:
                        tt = t.lower()
                        rebuilt.append(tt if tt in small else tt.capitalize())
                return ''.join(rebuilt).strip()

            def valid_name(name: str) -> bool:
                base = str(name).strip()
                if not base or len(base) <= 3:
                    return False
                if re.fullmatch(r"\d+(?:\.\d+)?", base):
                    return False
                return True

            aggregated: dict = {}
            for pathway_name, stats in self.pathway_filtered_pathways_data.items():
                # Skip metadata entries
                if pathway_name.startswith('_'):
                    continue
                if not isinstance(stats, dict):
                    continue

                disp = display_name(pathway_name)
                if not valid_name(disp):
                    continue
                key = canonical_key(disp)

                # Extract Z-score from multiple possible field names
                z_score = stats.get('z_score')
                if z_score is None:
                    z_score = stats.get('Z')
                if z_score is None:
                    z_score = stats.get('direction_zscore')
                if z_score is None:
                    z_score = 0.0
                
                row = {
                    'Pathway': disp,
                    'Status': stats.get('status', 'No Change'),
                    'Z_Score': z_score,
                    'Combined_P_Value': stats.get('combined_pvalue', stats.get('pvalue', 1.0)),
                    'N_Metabolites': stats.get('n_metabolites', stats.get('N_Metabolites', 0)),
                    'Mean_Log2FC': stats.get('mean_log2fc', stats.get('Mean_Log2FC', 0.0)),
                    'Metabolites': ', '.join(stats.get('metabolites', stats.get('Metabolites', []))),
                    'ML_RWR_Score': stats.get('ml_rwr_score', None),
                    'ML_PageRank': stats.get('ml_pagerank', None),
                    'ML_Betweenness': stats.get('ml_betweenness', None),
                }
                # Reaction role/role-scoring columns removed per requirement

                # Aggregate duplicates: prefer lower p-value, then higher |Z|, union metabolites
                if key not in aggregated:
                    aggregated[key] = row
                    aggregated[key]['_met_set'] = set([m.strip() for m in row['Metabolites'].split(',') if m.strip()]) if row.get('Metabolites') else set()
                else:
                    existing = aggregated[key]
                    p_new = row['Combined_P_Value'] if row['Combined_P_Value'] is not None else 1.0
                    p_old = existing['Combined_P_Value'] if existing['Combined_P_Value'] is not None else 1.0
                    choose_new = False
                    if p_new < p_old:
                        choose_new = True
                    elif abs(row['Z_Score'] or 0.0) > abs(existing['Z_Score'] or 0.0):
                        choose_new = True
                    if choose_new:
                        # Keep union of metabolites when replacing
                        existing_set = existing.get('_met_set', set())
                        new_set = set([m.strip() for m in row['Metabolites'].split(',') if m.strip()]) if row.get('Metabolites') else set()
                        row['_met_set'] = existing_set.union(new_set)
                        aggregated[key] = row
                    else:
                        # Update metabolite union only
                        new_set = set([m.strip() for m in row['Metabolites'].split(',') if m.strip()]) if row.get('Metabolites') else set()
                        aggregated[key]['_met_set'] = aggregated[key].get('_met_set', set()).union(new_set)

            # Finalize rows
            pathway_summary = []
            for item in aggregated.values():
                if '_met_set' in item:
                    item['Metabolites'] = ', '.join(sorted(item['_met_set']))
                    item.pop('_met_set', None)
                pathway_summary.append(item)

            # Compute consensus rank (average of ranks by |Z| and ML_RWR_Score)
            pathway_df = pd.DataFrame(pathway_summary)
            if not pathway_df.empty:
                try:
                    # Coerce to numeric safely
                    z_vals = pathway_df['Z_Score'].apply(lambda v: abs(float(v)) if pd.notna(v) else 0.0)
                    ml_vals = pathway_df['ML_RWR_Score'].apply(lambda v: float(v) if pd.notna(v) else 0.0) if 'ML_RWR_Score' in pathway_df.columns else pd.Series([0.0]*len(pathway_df))
                    # Rank descending (highest gets rank 1)
                    rank_z = z_vals.rank(ascending=False, method='min')
                    rank_ml = ml_vals.rank(ascending=False, method='min')
                    consensus = ((rank_z + rank_ml) / 2.0).round().astype(int)
                    pathway_df['Consensus_Rank'] = consensus
                except Exception:
                    pathway_df['Consensus_Rank'] = None
            print(f"DEBUG: Pathway summary has {len(pathway_df)} pathways")
            
            # Save to Excel: ONLY the Metabolites sheet (per user request)
            print(f"DEBUG: Writing to Excel file: {output_file}")
            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                metabolite_data.to_excel(writer, sheet_name='Metabolites', index=False)
            # Note: Intentionally not writing 'Pathway_Summary' sheet
            
            # Verify file was created
            if os.path.exists(output_file):
                file_size = os.path.getsize(output_file)
                print(f"DEBUG: File created successfully, size: {file_size} bytes")
                self.update_log_safe(f"✅ Excel file created ({file_size} bytes)")
            else:
                print(f"DEBUG: WARNING - File not found after writing: {output_file}")
                self.update_log_safe(f"⚠️ WARNING - File not found after writing")
            
            print(f"✅ Pathway annotation results saved to: {output_file}")
            self.update_log_safe(f"✅ Results saved to: {os.path.basename(output_file)} (Metabolites sheet only)")
            
        except Exception as e:
            print(f"ERROR in save_pathway_annotation_results: {str(e)}")
            import traceback
            traceback.print_exc()
            raise Exception(f"Failed to save results: {str(e)}")
    
    def export_pathway_annotation_results(self):
        """Export pathway annotation results to the specified output directory"""
        try:
            if self.pathway_filtered_metabolites_data is None or self.pathway_filtered_pathways_data is None:
                messagebox.showwarning("No Data", "No pathway annotation results to export. Please run pathway annotation first.")
                return
            
            # Show progress
            self.update_log_safe("💾 Exporting pathway annotation results...")
            
            # Save the Excel file (this will create a new timestamped file)
            self.save_pathway_annotation_results()
            
            # Get the saved file path
            output_file = self.pathway_output_file.get() if hasattr(self, 'pathway_output_file') else None
            output_dir = self.pathways_networks_dir if hasattr(self, 'pathways_networks_dir') else os.path.dirname(output_file) if output_file else "Unknown"
            
            # Inform user
            if output_file and os.path.exists(output_file):
                messagebox.showinfo("Export Complete", 
                                  f"Pathway annotation results exported successfully!\n\n"
                                  f"Folder: {output_dir}\n"
                                  f"Excel file: {os.path.basename(output_file)}\n\n"
                                  f"• Metabolites with pathway annotations (Metabolites sheet only)\n\n"
                                  f"Note: Visualizations are generated from the Pathway Network tab.")
                
                self.update_log_safe(f"✅ Results exported to: {output_file}")
            else:
                raise Exception("File was not created successfully")
            
        except Exception as e:
            error_msg = f"Export failed: {str(e)}"
            self.update_log_safe(f"❌ {error_msg}")
            messagebox.showerror("Export Error", error_msg)
            self.metabolite_neg_fc = tk.StringVar(value="0")

    def on_pathway_mode_change(self):
        """Handle pathway data mode change (metabolite/lipid)"""
        mode = self.pathway_data_mode.get()
        if mode == 'lipid':
            self.pathway_mode_info_label.config(
                text="ℹ️ Lipid mode: Feature ID (LipidID) required, Class_name recommended. Performs ID annotation → Lipid pathway annotation."
            )
        else:
            self.pathway_mode_info_label.config(
                text="ℹ️ Metabolite mode: Feature ID required. Optional: HMDB, KEGG, PubChem, ChEBI, CAS, SMILES, InChI, InChIKey."
            )
    
    def browse_pathway_annotation_folder(self):
        """NOTE: browse_output_folder removed - using project_folder from project setup"""
        pass
    
    def browse_working_folder(self):
        """Browse for working folder where results will be saved"""
        folder = filedialog.askdirectory(
            title="Select Working Folder for Project Results"
        )
        if folder:
            self.project_working_folder.set(folder)
            self._update_project_status()
            print(f"[DEBUG] Selected working folder: {folder}")
    
    def _update_project_status(self):
        """Update project setup status indicator"""
        folder = self.project_working_folder.get().strip()
        name = self.project_name.get().strip()
        
        if folder and name:
            self.project_status_label.config(
                text="✅ Project configured",
                fg='#28a745'
            )
            return True
        else:
            self.project_status_label.config(
                text="❌ Project not configured",
                fg='#d32f2f'
            )
            return False
    
    def _validate_project_setup(self):
        """Validate project setup before running annotation"""
        folder = self.project_working_folder.get().strip()
        name = self.project_name.get().strip()
        
        if not folder:
            messagebox.showerror(
                "Project Setup Required",
                "Please select a working folder before running annotation.\n\n"
                "Results will be automatically saved to this folder."
            )
            return False
        
        if not name:
            messagebox.showerror(
                "Project Setup Required",
                "Please enter a project name before running annotation.\n\n"
                "The project name will be used in the output filename."
            )
            return False
        
        if not os.path.exists(folder):
            messagebox.showerror(
                "Invalid Folder",
                f"The selected folder does not exist:\n{folder}"
            )
            return False
        
        # Check write permission
        try:
            test_file = os.path.join(folder, f".test_write_{int(time.time())}.tmp")
            with open(test_file, 'w') as f:
                f.write("test")
            os.remove(test_file)
        except Exception as e:
            messagebox.showerror(
                "Permission Error",
                f"Cannot write to selected folder:\n{folder}\n\n"
                f"Error: {str(e)}"
            )
            return False
        
        return True
    
    def browse_file(self, file_var, title):
        """Browse for an Excel file and load it"""
        filename = filedialog.askopenfilename(
            title=f"Select {title}",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if filename:
            file_var.set(filename)
            print(f"[DEBUG] Selected file: {filename}")
            
            # Auto-detect sheets for the selected file
            if hasattr(self, 'pathway_sheet_combo') and hasattr(self, 'pathway_sheet_var'):
                self.auto_detect_sheets(file_var, self.pathway_sheet_combo, self.pathway_sheet_var)
    
    def auto_detect_sheets(self, file_var, combo_widget, sheet_var):
        """Auto-detect and populate sheets from Excel file"""
        try:
            file_path = file_var.get()
            if not file_path or not file_path.endswith(('.xlsx', '.xls')):
                messagebox.showerror("Error", "Please select a valid Excel file first")
                return
            
            if not os.path.exists(file_path):
                messagebox.showerror("Error", "Selected file does not exist")
                return
            
            # GATING: Reset verification when new file is selected
            self.pvalue_verification_complete = False
            self.verified_pvalue_col = None
            self.log2fc_verification_complete = False
            self.verified_log2fc_col = None
            self._disable_run_button()
            self.pathway_progress_text.insert(tk.END, "⚠️ File changed - Feature ID & other ID verification required again\n\n")
            self.pathway_progress_text.see(tk.END)
            
            # Read sheet names
            excel_file = pd.ExcelFile(file_path)
            sheet_names = excel_file.sheet_names
            
            if not sheet_names:
                messagebox.showerror("Error", "No sheets found in the Excel file")
                return
            
            # Populate the combobox
            combo_widget['values'] = sheet_names
            
            # Auto-select the first sheet or a sheet with relevant data
            default_sheet = sheet_names[0]
            
            # Look for sheets that might contain metabolite data
            priority_keywords = ['metabolite', 'compound', 'data', 'results', 'analysis']
            for sheet in sheet_names:
                sheet_lower = str(sheet).lower()
                if any(keyword in sheet_lower for keyword in priority_keywords):
                    default_sheet = sheet
                    break
            
            # Set the selected sheet
            sheet_var.set(default_sheet)
            combo_widget.set(default_sheet)
            
            # Auto-detect columns for the selected sheet
            self._auto_detect_pathway_columns(file_path, default_sheet)
            
            messagebox.showinfo("Success", f"Found {len(sheet_names)} sheet(s). Selected: {default_sheet}\n\nColumns auto-detected. Click 'Verify Columns' to confirm or modify.")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read Excel file: {str(e)}")
    
    def _auto_detect_pathway_columns(self, file_path, sheet_name):
        """Auto-detect columns from the selected sheet and store assignments"""
        try:
            from gui.shared.column_assignment import ColumnDetector
            
            df = pd.read_excel(file_path, sheet_name=sheet_name)
            
            # Initialize assignments dictionary
            self.verified_assignments = {}
            
            # Get current mode
            mode = self.pathway_data_mode.get()
            
            # For lipid mode, prioritize LipidID detection
            if mode == 'lipid':
                # First, detect LipidID specifically (before Feature ID)
                for col in df.columns:
                    if ColumnDetector.detect_column(col, 'LipidID'):
                        self.verified_assignments['Feature ID'] = col
                        break
                
                # Detect other required columns
                required_for_pathway = ['P-Value', 'Log2 Fold Change']
                for col_type in required_for_pathway:
                    for col in df.columns:
                        if ColumnDetector.detect_column(col, col_type):
                            self.verified_assignments[col_type] = col
                            break
                
                # Also detect Class and Class_name for lipid mode
                for col_type in ['Class', 'Class_name']:
                    for col in df.columns:
                        # Special handling for Class_name (not in standard patterns)
                        if col_type == 'Class_name' and str(col).lower() in ['class_name', 'classname', 'class name']:
                            self.verified_assignments[col_type] = col
                            break
                        elif ColumnDetector.detect_column(col, col_type):
                            self.verified_assignments[col_type] = col
                            break
            else:
                # Metabolite mode: Prioritize exact 'Name' match for Feature ID
                # First, look for exact 'Name' column (case-insensitive)
                name_col = None
                for col in df.columns:
                    if str(col).strip().lower() == 'name':
                        name_col = col
                        break
                
                if name_col:
                    self.verified_assignments['Feature ID'] = name_col
                else:
                    # Fall back to standard Feature ID detection
                    for col in df.columns:
                        if ColumnDetector.detect_column(col, 'Feature ID'):
                            self.verified_assignments['Feature ID'] = col
                            break
                
                # Detect other required columns
                for col_type in ['P-Value', 'Log2 Fold Change']:
                    for col in df.columns:
                        if ColumnDetector.detect_column(col, col_type):
                            self.verified_assignments[col_type] = col
                            break
            
        except Exception as e:
            # Silently fail - auto-detection is optional
            import traceback
            print(f"Auto-detection failed: {e}\n{traceback.format_exc()}")
    
    def verify_pathway_columns(self):
        """Verify and assign columns for pathway annotation using dialog"""
        import threading
        from gui.shared.column_assignment import show_column_assignment_dialog
        
        def _load_and_verify():
            """Background thread worker for loading and showing dialog"""
            try:
                # Check if file is selected
                file_path = self.pathway_input_file.get()
                if not file_path:
                    self.root.after(0, lambda: messagebox.showerror("Error", "Please select an input file first"))
                    return
                
                if not os.path.exists(file_path):
                    self.root.after(0, lambda: messagebox.showerror("Error", "Selected file does not exist"))
                    return
                
                # Check if sheet is selected
                selected_sheet = self.pathway_sheet_var.get()
                if not selected_sheet:
                    self.root.after(0, lambda: messagebox.showerror("Error", "Please select a sheet first (use Auto-Detect button)"))
                    return
                
                # Get current mode
                mode = self.pathway_data_mode.get()
                mode_text = "Lipid" if mode == 'lipid' else "Metabolite"
                
                # Update UI
                self.root.after(0, lambda: self.pathway_progress_text.insert(tk.END, f"\n📂 Opening {mode_text} file - sheet '{selected_sheet}'...\n"))
                self.root.after(0, lambda: self.pathway_progress_text.see(tk.END))
                
                # Load dataframe
                df = pd.read_excel(file_path, sheet_name=selected_sheet)
                
                self.root.after(0, lambda: self.pathway_progress_text.insert(tk.END, "✓ Analyzing columns...\n"))
                self.root.after(0, lambda: self.pathway_progress_text.see(tk.END))
                
                # Use correct tab type based on mode
                # Lipid mode uses lipid_pathway_annotation (Feature ID + Class_name)
                # Metabolite mode uses pathway_annotation (Feature ID + optional IDs)
                if mode == 'lipid':
                    tab_type = 'lipid_pathway_annotation'
                else:
                    tab_type = 'pathway_annotation'
                
                # Show column assignment dialog
                result = show_column_assignment_dialog(
                    parent=self.root,
                    df=df,
                    tab_type=tab_type,
                    auto_calculate=False,
                    excel_file_path=file_path
                )
                
                if result:
                    # Store assignments
                    self.verified_assignments = result['assignments']
                    
                    # Rename columns in dataframe based on assignments for downstream compatibility
                    # Convert spaces to underscores for standard column names (HMDB ID → HMDB_ID)
                    for standard_name, user_column in self.verified_assignments.items():
                        if not user_column or user_column not in df.columns:
                            continue
                        # Do not convert the feature ID column to 'Feature_ID'; leave as-is
                        # so that it remains identifiable when annotation runs later.
                        if standard_name == 'Feature ID':
                            continue
                        # Convert standard name to underscore version for downstream code
                        underscore_name = standard_name.replace(' ', '_')
                        if user_column != underscore_name:
                            # Rename user's column to underscore version
                            df.rename(columns={user_column: underscore_name}, inplace=True)
                            self.root.after(0, lambda sn=underscore_name, uc=user_column: 
                                self.pathway_progress_text.insert(tk.END, f"   📝 Renamed '{uc}' → '{sn}'\n"))
                    
                    # Store renamed dataframe
                    if result.get('dataframe') is not None:
                        self.verified_dataframe = df  # Use renamed df
                    
                    # Get Feature ID assignment
                    feature_id = result['assignments'].get('Feature ID')
                    
                    # Build summary based on mode
                    if mode == 'lipid':
                        # Lipid mode: show Class_name and Class
                        class_name = result['assignments'].get('Class_name')
                        class_col = result['assignments'].get('Class')
                        lipidmaps_id = result['assignments'].get('LipidMaps ID')
                        
                        summary_lines = [f"{mode_text} Pathway Annotation Verification:"]
                        summary_lines.append(f"  ✓ Feature ID (LipidID): {feature_id}")
                        if class_name:
                            summary_lines.append(f"  ✓ Class_name: {class_name} (recommended for KEGG lookup)")
                        else:
                            summary_lines.append(f"  ⚠️ Class_name: Not assigned (KEGG pathway lookup may be limited)")
                        if class_col:
                            summary_lines.append(f"  ✓ Class: {class_col}")
                        if lipidmaps_id:
                            summary_lines.append(f"  ✓ LipidMaps ID: {lipidmaps_id}")
                        summary_lines.append("")
                        summary_lines.append("  📋 Workflow: ID Annotation → Lipid Pathway Annotation")
                    else:
                        # Metabolite mode: show ID columns
                        id_column_types = ['HMDB ID', 'KEGG ID', 'PubChem CID', 'ChEBI ID', 
                                          'CAS', 'SMILES', 'InChI', 'InChIKey', 'LipidMaps ID', 'Formula']
                        assigned_ids = []
                        for id_type in id_column_types:
                            assigned_col = result['assignments'].get(id_type)
                            if assigned_col:  # User assigned this ID type
                                assigned_ids.append(id_type)
                        
                        summary_lines = [f"{mode_text} Pathway Annotation Verification:"]
                        summary_lines.append(f"  ✓ Feature ID: {feature_id}")
                        if assigned_ids:
                            summary_lines.append(f"  ✓ ID Columns: {', '.join(assigned_ids)}")
                        else:
                            summary_lines.append(f"  ⚠️ No ID Columns (will use name-matching only)")
                    
                    summary_text = "\n".join(summary_lines)
                    
                    # Mark verification as complete
                    self.pvalue_verification_complete = True
                    self.log2fc_verification_complete = True
                    
                    # Update UI
                    self.root.after(0, lambda: (
                        self.pathway_progress_text.insert(tk.END, f"✅ Column verification completed\n{summary_text}\n"),
                        self.pathway_progress_text.see(tk.END),
                        self._enable_run_button(),
                        messagebox.showinfo(
                            "✅ Column Verification Complete",
                            summary_text + "\n\nYou can now run pathway annotation."
                        )
                    ))
                else:
                    # User cancelled
                    self.root.after(0, lambda: (
                        self.pathway_progress_text.insert(tk.END, "❌ Column verification cancelled\n"),
                        self.pathway_progress_text.see(tk.END)
                    ))
                    
            except Exception as e:
                self.pvalue_verification_complete = False
                self.root.after(0, lambda: (
                    self._disable_run_button(),
                    messagebox.showerror("Error", f"Failed to verify columns: {str(e)}"),
                ))
                logger.error(f"Verification error: {e}", exc_info=True)
        
        # Run in background thread
        thread = threading.Thread(target=_load_and_verify, daemon=True)
        thread.start()
    
    def _validate_pvalue_column(self, df, pvalue_col):
        """Validate that a column contains valid p-value data
        
        Returns dict with:
          - valid: bool indicating if column has valid p-values
          - min, max: value range
          - type: 'p-value' (0-1), 'log10_p', or 'invalid'
          - suggestions: string of alternative columns
        """
        try:
            # Convert to numeric and drop NaN
            pval_data = pd.to_numeric(df[pvalue_col], errors='coerce').dropna()
            
            if len(pval_data) == 0:
                return {
                    'valid': False,
                    'min': 0,
                    'max': 0,
                    'type': 'invalid',
                    'suggestions': 'Column contains no numeric data'
                }
            
            min_val = pval_data.min()
            max_val = pval_data.max()
            median_val = pval_data.median()
            
            # Check if values look like p-values (0-1 range)
            looks_like_pvalue = max_val <= 1 and min_val >= 0
            
            # Check if values look like -log10(p) (typically 0-10 range)
            looks_like_log10p = max_val > 1 and max_val <= 20 and min_val >= 0
            
            if looks_like_pvalue:
                return {
                    'valid': True,
                    'min': min_val,
                    'max': max_val,
                    'type': 'p-value (0-1 range)',
                    'suggestions': ''
                }
            elif looks_like_log10p:
                return {
                    'valid': True,
                    'min': min_val,
                    'max': max_val,
                    'type': '-log10(p-value)',
                    'suggestions': ''
                }
            else:
                # Invalid - find alternative columns
                suggestions = self._find_pvalue_columns(df, exclude=pvalue_col)
                return {
                    'valid': False,
                    'min': min_val,
                    'max': max_val,
                    'type': 'invalid (likely fold-change or other metric)',
                    'suggestions': suggestions
                }
        
        except Exception as e:
            return {
                'valid': False,
                'min': 0,
                'max': 0,
                'type': 'error',
                'suggestions': f'Validation error: {str(e)}'
            }
    
    def _find_pvalue_columns(self, df, exclude=None):
        """Find columns that look like p-values"""
        suggestions = ""
        for col in df.columns:
            if exclude and col == exclude:
                continue
            try:
                col_data = pd.to_numeric(df[col], errors='coerce').dropna()
                if len(col_data) > 0:
                    col_min = col_data.min()
                    col_max = col_data.max()
                    
                    # Check if this looks like a p-value
                    looks_like_p = col_max <= 1 and col_min >= 0
                    looks_like_log = col_max > 1 and col_max <= 20 and col_min >= 0
                    
                    if looks_like_p or looks_like_log:
                        ptype = "p-value" if looks_like_p else "-log10(p)"
                        suggestions += f"  ✓ '{col}': [{col_min:.6f}, {col_max:.6f}] ({ptype})\n"
            except:
                pass
        
        if not suggestions:
            suggestions = "  No obvious p-value columns found"
        
        return suggestions
    
    def _validate_log2fc_column(self, df, log2fc_col):
        """Validate that a column contains valid log2 fold-change data
        
        Returns dict with:
          - valid: bool indicating if column has valid fold-change values
          - min, max: value range
          - type: 'log2_fold_change' or 'invalid'
          - suggestions: string of alternative columns
        """
        try:
            # Convert to numeric and drop NaN
            log2fc_data = pd.to_numeric(df[log2fc_col], errors='coerce').dropna()
            
            if len(log2fc_data) == 0:
                return {
                    'valid': False,
                    'min': 0,
                    'max': 0,
                    'type': 'invalid',
                    'suggestions': 'Column contains no numeric data'
                }
            
            min_val = log2fc_data.min()
            max_val = log2fc_data.max()
            median_val = log2fc_data.median()
            mean_val = log2fc_data.mean()
            std_val = log2fc_data.std()
            
            # Log2 fold-change values are typically in range [-10, 10]
            # But we allow wider range to accommodate different data
            # Key indicators: should include negative values (down-regulation)
            # and typically centered near 0
            
            is_likely_log2fc = (
                (min_val < 0 and max_val > 0) or  # Has both up and down regulation
                (abs(mean_val) < 5 and abs(median_val) < 5 and std_val < 10)  # Centered near 0
            )
            
            # Reject if values look like p-values (all positive, 0-1 range)
            is_pvalue_like = (max_val <= 1 and min_val >= 0 and max_val > 0)
            
            # Reject if values look like -log10(p) (typically 0-20, positive only)
            is_log10p_like = (min_val >= 0 and max_val > 1 and max_val <= 20 and 
                              (mean_val > 1 or median_val > 1))
            
            if is_pvalue_like:
                # Detected p-values in log2FC column
                suggestions = self._find_log2fc_columns(df, exclude=log2fc_col)
                return {
                    'valid': False,
                    'min': min_val,
                    'max': max_val,
                    'type': 'invalid (detected p-values, not fold-change)',
                    'suggestions': suggestions
                }
            elif is_log10p_like:
                # Detected -log10(p) values
                suggestions = self._find_log2fc_columns(df, exclude=log2fc_col)
                return {
                    'valid': False,
                    'min': min_val,
                    'max': max_val,
                    'type': 'invalid (detected -log10(p), not fold-change)',
                    'suggestions': suggestions
                }
            elif is_likely_log2fc or (max_val <= 20 and min_val >= -20):
                # Looks like valid log2 fold-change
                return {
                    'valid': True,
                    'min': min_val,
                    'max': max_val,
                    'type': f'log2_fold_change (mean: {mean_val:.4f}, std: {std_val:.4f})',
                    'suggestions': ''
                }
            else:
                # Invalid - find alternative columns
                suggestions = self._find_log2fc_columns(df, exclude=log2fc_col)
                return {
                    'valid': False,
                    'min': min_val,
                    'max': max_val,
                    'type': f'invalid (range [{min_val:.6f}, {max_val:.6f}])',
                    'suggestions': suggestions
                }
        
        except Exception as e:
            return {
                'valid': False,
                'min': 0,
                'max': 0,
                'type': 'error',
                'suggestions': f'Validation error: {str(e)}'
            }
    
    def _find_log2fc_columns(self, df, exclude=None):
        """Find columns that look like log2 fold-change data"""
        suggestions = ""
        for col in df.columns:
            if exclude and col == exclude:
                continue
            # Skip obvious metadata columns
            if any(skip in col.lower() for skip in ['id', 'name', 'gene', 'protein', 'accession', 'compound']):
                continue
            
            try:
                col_data = pd.to_numeric(df[col], errors='coerce').dropna()
                if len(col_data) > 10:  # Need enough data points
                    col_min = col_data.min()
                    col_max = col_data.max()
                    col_mean = col_data.mean()
                    col_median = col_data.median()
                    
                    # Check if this looks like fold-change
                    has_neg_pos = (col_min < 0 and col_max > 0)
                    is_centered = (abs(col_mean) < 5 and abs(col_median) < 5)
                    reasonable_range = (col_max <= 20 and col_min >= -20)
                    
                    if (has_neg_pos or is_centered) and reasonable_range:
                        suggestions += f"  ✓ '{col}': [{col_min:.6f}, {col_max:.6f}] (mean: {col_mean:.4f})\n"
            except:
                pass
        
        if not suggestions:
            suggestions = "  No obvious log2 fold-change columns found"
        
        return suggestions
    
    def _enable_run_button(self):
        """Enable Run Annotation button (gates on verified p-value and log2fc columns)"""
        self.pathway_annotation_button.config(
            state='normal',
            bg='#27ae60',
            relief='raised'
        )
        self.verify_status_label.config(
            text="✅ P-Value & Log2FC VERIFIED - Ready to run",
            bg='#e8f5e9',
            fg='#2e7d32'
        )
    
    def _disable_run_button(self):
        """Disable Run Annotation button (gates on verified p-value and log2fc columns)"""
        self.pathway_annotation_button.config(
            state='disabled',
            bg='#95a5a6',
            relief='sunken'
        )
        self.verify_status_label.config(
            text="⚠️ Column verification REQUIRED (Feature ID & Metabolite IDs)",
            bg='#ffe6e6',
            fg='#d32f2f'
        )
    
    def run_pathway_annotation_step(self):
        """Step 1: Run pathway annotation and populate tables"""
        
        # ========================================
        # VALIDATE PROJECT SETUP FIRST
        # ========================================
        if not self._validate_project_setup():
            return
        
        # GATING: Check if BOTH p-value and log2FC verification are complete
        if not self.pvalue_verification_complete:
            messagebox.showerror(
                "Verification Required",
                "⛔ You must verify the P-Value column before running annotation.\n\n"
                "Steps:\n"
                "1. Click '🔍 Verify Columns' button\n"
                "2. Select the correct P-Value column\n"
                "3. Wait for validation confirmation\n"
                "4. Then click 'Run Annotation'\n\n"
                "This ensures accurate pathway analysis."
            )
            return
        
        if not self.log2fc_verification_complete:
            messagebox.showerror(
                "Verification Required",
                "⛔ You must verify the Log2 Fold Change column before running annotation.\n\n"
                "Steps:\n"
                "1. Click '🔍 Verify Columns' button\n"
                "2. Select the correct Log2 Fold Change column\n"
                "3. Wait for validation confirmation\n"
                "4. Then click 'Run Annotation'\n\n"
                "This ensures accurate pathway analysis."
            )
            return
        
        if not self.pathway_input_file.get():
            messagebox.showerror("Error", "Please select an input file")
            return
        
        if not os.path.exists(self.pathway_input_file.get()):
            messagebox.showerror("Error", "Input file does not exist")
            return
        
        # Check if a sheet is selected
        if not self.pathway_sheet_var.get():
            messagebox.showerror("Error", "Please select a sheet first or click 'Auto-Detect Sheets'")
            return
        
        # Initialize cancel flag
        self.pathway_annotation_cancelled = False
        
        # Enable Stop button, disable Annotate and Export buttons during processing
        if hasattr(self, 'pathway_stop_button'):
            self.pathway_stop_button.config(state='normal')
        if hasattr(self, 'pathway_annotation_button'):
            self.pathway_annotation_button.config(state='disabled')
        if hasattr(self, 'export_pathway_results_button'):
            self.export_pathway_results_button.config(state='disabled')
        
        # Initialize progress tracking
        self.pathway_start_time = time.time()
        self.pathway_current_step = 0
        self.pathway_total_steps = 100
        
        # Setup progress bar
        self.pathway_progress_bar.config(mode='determinate', maximum=100, value=0)
        self.pathway_progress_percentage.config(text="0%")
        self.pathway_timer_label.config(text="Time: 00:00")
        self.pathway_eta_label.config(text="ETA: --:--")
        
        # Clear progress text
        self.pathway_progress_text.delete(1.0, tk.END)
        self.pathway_progress_text.insert(tk.END, "=== PATHWAY ANNOTATION PROCESS STARTED ===\n\n")
        self.pathway_progress_text.see(tk.END)
        
        # Start timer update
        self.update_pathway_timer()
        
        # Spawn background thread to run annotation
        self.annotation_thread = threading.Thread(target=self._run_annotation_process, daemon=True)
        self.annotation_thread.start()
    
    def update_pathway_timer(self):
        """Update pathway annotation timer"""
        if getattr(self, '_shutting_down', False):
            return
        if hasattr(self, 'pathway_start_time') and self.pathway_start_time is not None:
            elapsed = time.time() - self.pathway_start_time
            minutes = int(elapsed // 60)
            seconds = int(elapsed % 60)
            self.pathway_timer_label.config(text=f"Time: {minutes:02d}:{seconds:02d}")
            
            # Calculate ETA
            if hasattr(self, 'pathway_current_step') and self.pathway_current_step > 0:
                rate = elapsed / self.pathway_current_step
                remaining_steps = self.pathway_total_steps - self.pathway_current_step
                eta_seconds = remaining_steps * rate
                eta_minutes = int(eta_seconds // 60)
                eta_secs = int(eta_seconds % 60)
                self.pathway_eta_label.config(text=f"ETA: {eta_minutes:02d}:{eta_secs:02d}")
            
            # Schedule next update
            self.root.after(1000, self.update_pathway_timer)
    
    def _run_annotation_process(self):
        """Background thread method: Run pathway annotation without blocking GUI"""
        try:
            # Check if cancelled before starting
            if getattr(self, 'pathway_annotation_cancelled', False):
                self.root.after(0, lambda: self._handle_cancellation())
                return
            
            # Update progress
            self.root.after(0, lambda: self.update_pathway_progress_with_percentage(10, "Loading input data..."))
            
            # Read Excel file
            selected_sheet = self.pathway_sheet_var.get()
            df = pd.read_excel(self.pathway_input_file.get(), sheet_name=selected_sheet)
            
            # Check if cancelled after loading
            if getattr(self, 'pathway_annotation_cancelled', False):
                self.root.after(0, lambda: self._handle_cancellation())
                return
            
            # Log data info (thread-safe UI update)
            def log_message(msg):
                self.root.after(0, lambda: (
                    self.pathway_progress_text.insert(tk.END, msg),
                    self.pathway_progress_text.see(tk.END)
                ))
            
            log_message(f"\n📂 Loaded sheet: '{selected_sheet}'\n")
            log_message(f"📊 Loaded data: {len(df)} rows, {len(df.columns)} columns\n")
            log_message(f"📋 Columns: {list(df.columns)}\n\n")
            
            # Normalize column names
            df.columns = df.columns.map(lambda c: str(c).strip())

            # If user verified columns, normalize them to underscore names for downstream checks
            # NOTE: feature ID column is treated specially because the mapper requires a 'Name' column later.
            if hasattr(self, 'verified_assignments') and self.verified_assignments:
                try:
                    for standard_name, user_column in self.verified_assignments.items():
                        if not user_column or user_column not in df.columns:
                            continue
                        # Do NOT convert Feature ID to "Feature_ID"; leave it alone
                        # so that the later renaming logic can turn it into 'Name'.
                        if standard_name == 'Feature ID':
                            continue
                        underscore_name = standard_name.replace(' ', '_')
                        if user_column != underscore_name:
                            df.rename(columns={user_column: underscore_name}, inplace=True)
                except Exception:
                    pass
            
            # Get current mode
            mode = self.pathway_data_mode.get()
            
            # ===== STEP 1: SKIP P-Value/Log2FC for Annotation =====
            # Pathway annotation only needs Feature ID + optional ID columns
            # P-Value and Log2FC are for Network Analysis tab only
            self.root.after(0, lambda: self.update_pathway_progress_with_percentage(11, "Preparing metabolites for annotation..."))
            
            log_message(f"✅ All {len(df)} metabolites will be annotated\n")
            log_message(f"   (P-value/Log2FC not required for annotation)\n\n")
            
            # ===== EARLY CHECK: Detect existing pathway annotation to SKIP re-annotation =====
            # Check for a broad set of common pathway columns (variants included)
            pathway_column_variants = [
                'HMDB_Pathways', 'PathBank_Pathways', 'PathBank_Pathways', 'PathBank_Pathways',
                'SMPDB_Pathways', 'SMP_Pathways', 'SMPDB_Pathways', 'WikiPathways', 'Metabolika_Pathways',
                'Reactome_Pathways', 'KEGG_Pathways', 'All_Pathways', 'All_Pathways_Display', 'total_pathways'
            ]
            present_annot_cols = [c for c in pathway_column_variants if c in df.columns]
            skip_annotation_early = False
            if present_annot_cols:
                # Skip annotation if any aggregated column exists with data
                agg_cols = {'All_Pathways', 'All_Pathways_Display'}
                agg_present = [c for c in present_annot_cols if c in agg_cols]
                if agg_present:
                    try:
                        has_any = any(df[c].notna().any() for c in agg_present)
                    except Exception:
                        has_any = False
                    if has_any:
                        skip_annotation_early = True
                        log_message(f"✅ Detected aggregated pathway columns: {agg_present}\n")
                        log_message("   Will SKIP re-annotation and proceed directly to reclassification/statistics.\n\n")
                else:
                    log_message(f"ℹ️ Found pathway columns {present_annot_cols} but no aggregated field; annotation will run\n")

            # ===== LIPID MODE: Check if IDs already present, skip re-annotation if yes =====
            if mode == 'lipid' and not skip_annotation_early:
                # Check if lipid IDs are already present from ID annotation tab
                id_cols_to_check = ['LipidMaps_ID', 'KEGG_ID', 'HMDB_ID', 'PubChem_CID', 'ChEBI_ID']

                # Prefer user-verified assignments if available
                existing_id_cols = []
                if hasattr(self, 'verified_assignments') and self.verified_assignments:
                    assigned_types = ['LipidMaps ID', 'KEGG ID', 'HMDB ID', 'PubChem CID', 'ChEBI ID']
                    for col_type in assigned_types:
                        assigned_col = self.verified_assignments.get(col_type)
                        if assigned_col:
                            assigned_col = col_type.replace(' ', '_')
                            if assigned_col in df.columns:
                                existing_id_cols.append(assigned_col)

                # Fallback: Match columns by normalized header to accept spaces/underscores/case variants
                if not existing_id_cols:
                    def _norm_header(val: str) -> str:
                        return re.sub(r'[^a-z0-9]', '', str(val).lower())
                    normalized_cols = {_norm_header(c): c for c in df.columns}
                    for expected in id_cols_to_check:
                        key = _norm_header(expected)
                        if key in normalized_cols:
                            existing_id_cols.append(normalized_cols[key])
                
                has_existing_ids = False
                if existing_id_cols:
                    # Check if any of these ID columns have non-empty values
                    try:
                        def _count_non_empty(series: pd.Series) -> int:
                            # Treat NaN/None and blank strings as empty
                            s = series.dropna()
                            if s.empty:
                                return 0
                            if s.dtype == object:
                                s = s.astype(str).str.strip()
                                s = s[~s.str.lower().isin(['', 'nan', 'none'])]
                            return int(s.shape[0])

                        for id_col in existing_id_cols:
                            if id_col not in df.columns:
                                continue
                            non_empty = _count_non_empty(df[id_col])
                            if non_empty > 0:
                                has_existing_ids = True
                                log_message(f"✅ Found existing IDs in column '{id_col}': {non_empty} lipids\n")
                    except Exception:
                        pass
                
                if has_existing_ids:
                    log_message(f"\n🎯 LIPID IDs DETECTED - Proceeding with pathway annotation\n")
                    log_message(f"   ID columns present: {existing_id_cols}\n")
                    log_message(f"   Proceeding directly to pathway annotation...\n\n")
                    skip_annotation_early = True
                    self.root.after(0, lambda: self.update_pathway_progress_with_percentage(20, "Lipid IDs detected - proceeding with pathway annotation"))
                else:
                    # No existing IDs found - STOP and tell user to run ID annotation first
                    try:
                        log_message("⚠️ Lipid ID check failed. Column summary:\n")
                        for id_col in existing_id_cols or []:
                            if id_col in df.columns:
                                non_empty = df[id_col].dropna().shape[0]
                                log_message(f"   • {id_col}: {non_empty} non-null values\n")
                    except Exception:
                        pass
                    error_msg = (
                        "⚠️⚠️⚠️ LIPID IDs NOT FOUND ⚠️⚠️⚠️\n\n"
                        "The input file does not contain lipid ID annotations.\n"
                        "Expected ID columns: LipidMaps_ID, KEGG_ID, HMDB_ID, PubChem_CID, ChEBI_ID\n\n"
                        "📋 REQUIRED STEPS:\n"
                        "1. Go to the 'ID Annotation' tab\n"
                        "2. Load your cleaned lipid file\n"
                        "3. Run lipid ID annotation (LipidMAPS + KEGG + HMDB + PubChem)\n"
                        "4. Return to this tab and load the ID-annotated file\n\n"
                        "⛔ Pathway annotation cannot proceed without lipid IDs.\n"
                    )
                    log_message(error_msg)
                    
                    self.root.after(0, lambda: self.update_pathway_progress_with_percentage(0, "Error: Lipid IDs required"))
                    self.root.after(0, lambda: self.pathway_progress_text.config(state='normal'))
                    self.root.after(0, lambda: self.pathway_progress_text.insert('end', f"\n{error_msg}\n"))
                    self.root.after(0, lambda: self.pathway_progress_text.config(state='disabled'))
                    self.root.after(0, lambda: self.pathway_progress_text.see('end'))
                    self.root.after(0, lambda: self.pathway_annotation_button.config(state='normal'))
                    
                    # Show messagebox to user
                    self.root.after(0, lambda: messagebox.showerror(
                        "Lipid IDs Required",
                        "Your lipid file does not have ID annotations.\n\n"
                        "Please go to the ID Annotation tab first to:\n"
                        "• Run lipid ID annotation (LipidMAPS + KEGG + HMDB)\n"
                        "• Then return here with the ID-annotated file\n\n"
                        "Pathway annotation requires lipid IDs to map to pathways."
                    ))
                    
                    return  # Exit - cannot proceed without IDs
            
            # Before any auto-detection: if user verified p-value column, enforce it now.
            if hasattr(self, 'verified_pvalue_col') and self.verified_pvalue_col:
                pcol = self.verified_pvalue_col
                if pcol in df.columns:
                    # Re-apply canonical pvalue column mapping after any DataFrame replacement
                    canonical_pvals = pd.to_numeric(df[pcol], errors='coerce')
                    orig_na_mask = canonical_pvals.isna()
                    n_na = int(orig_na_mask.sum())
                    df['pvalue'] = canonical_pvals  # Keep NaN values as-is
                    if n_na > 0:
                        na_names = df.loc[orig_na_mask, 'Name'].tolist() if 'Name' in df.columns else []
                        log_message(f"⚠️ {n_na} metabolites have missing p-values (post-replacement)\n")
                        if na_names:
                            log_message(f"   Examples: {na_names[:5]}\n")
                        log_message(f"   These will be ANNOTATED but excluded from Fisher ORA / IWPA statistics later\n")
                        log_message(f"   Total metabolites for annotation: {len(df)} (including {n_na} without p-values)\n")
                    else:
                        log_message(f"✅ Verified p-value column present: '{pcol}' - {len(df)} metabolites retained\n")

            # Annotation tab should NOT detect or care about p-values or fold changes
            # Those are comparison columns that the network tab will handle
            # Just focus on adding pathway annotations
            
            # Check if pathway annotation is needed
            self.root.after(0, lambda: self.update_pathway_progress_with_percentage(15, "Checking existing pathway data..."))
            pathway_cols = ['HMDB_Pathways', 'PathBank_Pathways', 'SMP_Pathways', 'WikiPathways', 
                           'Metabolika_Pathways', 'All_Pathways', 'All_Pathways_Display']
            present_pathway_cols = [c for c in pathway_cols if c in df.columns]
            
            skip_annotation = False
            if present_pathway_cols:
                # Always skip if an aggregated column is present and contains data
                agg_cols = {'All_Pathways', 'All_Pathways_Display'}
                agg_present = [c for c in present_pathway_cols if c in agg_cols]
                if agg_present:
                    has_pathway_data = False
                    for col in agg_present:
                        if df[col].notna().any():
                            has_pathway_data = True
                            break
                else:
                    has_pathway_data = False
                
                if has_pathway_data:
                    skip_annotation = True
                    log_message(f"✅ Found existing pathway data!\n")
                    log_message(f"   Columns: {present_pathway_cols}\n")
                    log_message(f"   Skipping annotation step and proceeding directly to statistics...\n\n")
            
            if skip_annotation:
                # Use existing pathway data but enforce filtering and disease relocation
                log_message(f"📊 Using annotated data from loaded file\n")
                log_message(f"   Skipping database loading (not needed for statistics only)\n")
                
                # CRITICAL CHECK: If Associated_Diseases already exists with data, don't reclassify!
                # Reclassification is only needed for "dirty" data where diseases are still mixed in pathway columns
                has_diseases_column = 'Associated_Diseases' in df.columns
                diseases_populated = False
                if has_diseases_column:
                    non_empty = df['Associated_Diseases'].fillna('').astype(str).str.strip()
                    non_empty = non_empty[(non_empty != '') & (non_empty.str.lower() != 'nan')]
                    diseases_populated = len(non_empty) > 0
                
                if has_diseases_column and diseases_populated:
                    # Data already has diseases separated - use as-is
                    log_message(f"✓ Associated_Diseases column already populated ({len(non_empty)} metabolites)\n")
                    log_message(f"   Skipping reclassification (data already clean)\n\n")
                    # Still sanitize any pre-existing 'All_Pathways' semicolon fragmentation
                    try:
                        from main_script.metabolite_pathways_annotator import MetabolitePathwayMapper
                        if 'All_Pathways' in df.columns:
                            cleaned_count = 0
                            examples = []
                            for idx, val in df['All_Pathways'].fillna('').astype(str).items():
                                if not val.strip():
                                    continue
                                new_val = MetabolitePathwayMapper.group_semicolon_fragments(val)
                                if new_val != val:
                                    cleaned_count += 1
                                    if len(examples) < 5:
                                        examples.append((val, new_val))
                                    df.at[idx, 'All_Pathways'] = new_val
                            if cleaned_count > 0:
                                log_message(f"🔧 Cleaned {cleaned_count} pre-existing 'All_Pathways' entries by grouping semicolon fragments\n")
                                for old, new in examples:
                                    log_message(f"   • Example: '{old}' → '{new}'\n")
                    except Exception as _e:
                        log_message(f"⚠️ Warning: Could not run All_Pathways cleaning: {_e}\n")

                    annotated_df = df
                else:
                    # Need to reclassify diseases from pathway columns
                    log_message(f"⚠️ Associated_Diseases empty or missing - will reclassify from pathway columns\n")
                    
                    try:
                        # For skip_annotation path, we don't need full annotator initialization
                        # Just use the reclassify function directly
                        from main_script.metabolite_pathways_annotator import apply_annotation_filters, is_disease_pathway
                        
                        # Force rebuild All_Pathways columns to ensure disease filtering
                        # Drop old All_Pathways columns if they exist (may contain diseases)
                        df_clean = df.copy()

                        # Sanitize individual pathway columns for semicolon fragmentation
                        try:
                            from main_script.metabolite_pathways_annotator import MetabolitePathwayMapper
                            for pcol in pathway_cols:
                                if pcol in df_clean.columns:
                                    def _maybe_group(v):
                                        try:
                                            if pd.isna(v):
                                                return v
                                            s = str(v)
                                            if ';' in s and len(s) > 3:
                                                return MetabolitePathwayMapper.group_semicolon_fragments(s)
                                        except Exception:
                                            pass
                                        return v
                                    df_clean[pcol] = df_clean[pcol].apply(_maybe_group)
                        except Exception as _e:
                            log_message(f"⚠️ Warning: Could not sanitize pathway columns: {_e}\n")
                        
                        # Reclassify diseases using the pure function approach
                        # Build All_Pathways from individual pathway columns
                        pathway_cols = ['HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 
                                       'WikiPathways', 'Metabolika_Pathways', 'Reactome_Pathways', 'KEGG_Pathways']
                        
                        # Also check for Reactome_Disease column
                        disease_cols = ['Reactome_Disease']
                        
                        def reclassify_for_stats(row):
                            """Apply disease filtering to pathway columns"""
                            import re
                            all_pathways = []
                            all_diseases = []
                            
                            # DEBUG: Track input
                            debug_info = {'metabolite': row.get('Name', 'Unknown'), 'before': len(all_pathways)}
                            
                            # First, collect diseases from disease-specific columns
                            for col in disease_cols:
                                if col in row.index and pd.notna(row[col]):
                                    val = str(row[col])
                                    if val and val.lower() not in ['nan', 'none', '']:
                                        # Split by both ; and | delimiters
                                        diseases = [d.strip() for d in re.split(r'[;|]', val) if d.strip()]
                                        all_diseases.extend(diseases)
                            
                            # Then process pathway columns
                            for col in pathway_cols:
                                if col in row.index and pd.notna(row[col]):
                                    val = str(row[col])
                                    if val and val.lower() not in ['nan', 'none', '']:
                                        # Split by both ; and | delimiters (HMDB uses |, others use ;)
                                        pathways = [p.strip() for p in re.split(r'[;|]', val) if p.strip()]
                                        # DEBUG: Track before filtering
                                        before_filter = len(pathways)
                                        # Apply disease filtering
                                        clean_pathways, diseases = apply_annotation_filters(pathways)
                                        # DEBUG: Track after filtering
                                        after_filter = len(clean_pathways)
                                        if before_filter != after_filter:
                                            debug_info[col] = f"{before_filter}→{after_filter}"
                                        all_pathways.extend(clean_pathways)
                                        all_diseases.extend(diseases)
                            
                            # DEBUG: Track before deduplication
                            debug_info['total_collected'] = len(all_pathways)
                            
                            # Deduplicate
                            unique_pathways = []
                            seen_pathways = set()
                            duplicates_found = []
                            for p in all_pathways:
                                p_lower = p.lower()
                                if p_lower not in seen_pathways:
                                    seen_pathways.add(p_lower)
                                    unique_pathways.append(p)
                                else:
                                    duplicates_found.append(p)
                            
                            # DEBUG: Track deduplication
                            debug_info['after_dedup'] = len(unique_pathways)
                            debug_info['duplicates_removed'] = len(duplicates_found)
                            if duplicates_found:
                                debug_info['duplicate_examples'] = duplicates_found[:3]
                            
                            unique_diseases = []
                            seen_diseases = set()
                            for d in all_diseases:
                                d_lower = d.lower()
                                if d_lower not in seen_diseases:
                                    seen_diseases.add(d_lower)
                                    unique_diseases.append(d)
                            
                            # Use proper delimiters: ; for pathways, | for diseases
                            pathways_str = ';'.join(unique_pathways) if unique_pathways else ''
                            diseases_str = ' | '.join(unique_diseases) if unique_diseases else ''
                            
                            return pathways_str, diseases_str, debug_info
                        
                        # Apply reclassification
                        log_message(f"\n📋 Reclassifying diseases from pathway columns...\n")
                        
                        # DEBUG: Check BEFORE reclassification
                        if 'Associated_Diseases' in df_clean.columns:
                            before_count = (df_clean['Associated_Diseases'].fillna('').astype(str).str.strip() != '').sum()
                            before_sample = df_clean['Associated_Diseases'].dropna().head(2).tolist()
                            log_message(f"   BEFORE: {before_count} metabolites with Associated_Diseases\n")
                            log_message(f"   BEFORE sample: {before_sample}\n")
                        
                        result = df_clean.apply(reclassify_for_stats, axis=1, result_type='expand')
                        df_clean['All_Pathways'] = result[0]
                        df_clean['Associated_Diseases'] = result[1]
                        debug_info_list = result[2].tolist()
                        
                        # DEBUG: Analyze deduplication results
                        total_duplicates_removed = sum(d.get('duplicates_removed', 0) for d in debug_info_list)
                        metabolites_with_duplicates = sum(1 for d in debug_info_list if d.get('duplicates_removed', 0) > 0)
                        
                        log_message(f"\n🔍 DEDUPLICATION ANALYSIS:\n")
                        log_message(f"   • Total duplicate pathways removed: {total_duplicates_removed}\n")
                        log_message(f"   • Metabolites with duplicates: {metabolites_with_duplicates}/{len(debug_info_list)}\n")
                        
                        # Show examples of duplicates removed
                        examples_shown = 0
                        for d in debug_info_list:
                            if d.get('duplicates_removed', 0) > 0 and examples_shown < 5:
                                log_message(f"   • {d.get('metabolite', 'Unknown')}: {d.get('total_collected', 0)} collected → {d.get('after_dedup', 0)} unique\n")
                                if 'duplicate_examples' in d:
                                    log_message(f"     Duplicates: {', '.join(d['duplicate_examples'])}\n")
                                examples_shown += 1
                        
                        # DEBUG: Check AFTER reclassification
                        after_count = (df_clean['Associated_Diseases'].fillna('').astype(str).str.strip() != '').sum()
                        after_sample = df_clean['Associated_Diseases'].dropna().head(2).tolist()
                        log_message(f"\n   AFTER: {after_count} metabolites with Associated_Diseases\n")
                        log_message(f"   AFTER sample: {after_sample}\n")
                        
                        # Count diseases
                        diseases_count = (df_clean['Associated_Diseases'].fillna('').astype(str).str.strip() != '').sum()
                        pathways_count = (df_clean['All_Pathways'].fillna('').astype(str).str.strip() != '').sum()
                        log_message(f"\n   ✓ Reclassification complete:\n")
                        log_message(f"      - {diseases_count} metabolites with diseases\n")
                        log_message(f"      - {pathways_count} metabolites with pathways\n")
                            
                        annotated_df = df_clean
                            
                        # Note: All_Pathways is already a string from reclassify_for_stats, no conversion needed
                        log_message(f"✓ Applied disease relocation and pathway filtering to pre-annotated data\n")
                        log_message(f"   Pathways column ready ({len(annotated_df)} metabolites)\n")
                        
                    except Exception as _e:
                        log_message(f"⚠️ Could not post-process pre-annotated data: {_e}\n")
                        annotated_df = df
            else:
                # Check cancellation before running annotation
                if getattr(self, 'pathway_annotation_cancelled', False):
                    self.root.after(0, lambda: self._handle_cancellation())
                    return
                
                # Run pathway annotation
                self.root.after(0, lambda: self.update_pathway_progress_with_percentage(15, "Running pathway annotation..."))
                log_message(f"🛣️ Mapping metabolites to pathways...\n")
                
                # Import pathway mapper
                from main_script.metabolite_pathways_annotator import MetabolitePathwayMapper
                
                # Create a progress callback for the annotator (thread-safe)
                def progress_callback(current, total, message):
                    percentage = int(5 + (current / total) * 85)  # 5-90% range
                    self.root.after(0, lambda: self.update_pathway_progress_with_percentage(percentage, f"Pathway mapping: {current}/{total}"))
                
                # Get organism selection
                organism = self.pathway_organism_var.get() if hasattr(self, 'pathway_organism_var') else 'Homo sapiens'
                
                # Get parallel workers setting
                try:
                    workers = int(self.pathway_workers.get())
                except (ValueError, AttributeError):
                    workers = 6
                
                # NOTE: P-value threshold NOT needed for annotation
                # (only needed for Network Analysis tab)
                
                # Create temporary directories and files
                temp_dir = tempfile.gettempdir()
                # Ensure the DataFrame uses the user-verified feature ID column as 'Name'
                # so downstream annotators rely on the user's assignment rather than
                # attempting to auto-detect or use fallbacks.
                try:
                    feature_col = None
                    if hasattr(self, 'verified_assignments') and self.verified_assignments:
                        feature_col = self.verified_assignments.get('Feature ID') or self.verified_assignments.get('LipidID') or self.verified_assignments.get('KEGG_ID')

                    if feature_col and feature_col in df.columns and feature_col != 'Name':
                        # If 'Name' column already exists, drop it first to avoid duplicates
                        if 'Name' in df.columns:
                            df = df.drop(columns=['Name'])
                            log_message(f"🔄 Dropped existing 'Name' column (will use '{feature_col}' instead)\n")
                        df = df.rename(columns={feature_col: 'Name'})
                        log_message(f"🔁 Renamed verified feature column '{feature_col}' to 'Name' for annotation\n")
                    
                    # CRITICAL: Check for and fix duplicate column names
                    if df.columns.duplicated().any():
                        dup_cols = df.columns[df.columns.duplicated()].tolist()
                        log_message(f"⚠️ Found duplicate columns: {dup_cols} - removing duplicates\n")
                        df = df.loc[:, ~df.columns.duplicated()]
                    # DEBUG: ensure 'Name' exists before annotation
                    if 'Name' not in df.columns:
                        log_message("⚠️ WARNING: 'Name' column missing after renaming step – annotation may fail\n")
                except Exception as e:
                    # If any unexpected issue occurs, continue without renaming (mapper has fallbacks)
                    log_message(f"⚠️ Column rename issue: {e}\n")
                    pass

                temp_input_path = os.path.join(temp_dir, f'pathway_input_{os.getpid()}.xlsx')
                df.to_excel(temp_input_path, index=False, sheet_name='Sheet1')
                
                # Create output path using project folder and name
                output_folder = self.project_working_folder.get()
                if not os.path.exists(output_folder):
                    os.makedirs(output_folder, exist_ok=True)
                
                # Use project name for output file
                project_name = self.project_name.get() or 'pathway'
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                final_output_path = os.path.join(output_folder, f'{project_name}_pathway_annotated_{timestamp}.xlsx')
                temp_output_path = os.path.join(temp_dir, f'pathway_annotated_{os.getpid()}.xlsx')
                
                # Create a cancel check function for the annotator
                def check_cancelled():
                    return getattr(self, 'pathway_annotation_cancelled', False)
                
                # Run pathway annotation using MetabolitePathwayMapper
                try:
                    annotator = MetabolitePathwayMapper(
                        input_file=temp_input_path,
                        output_file=temp_output_path,
                        sheet_name='Sheet1',
                        progress_callback=progress_callback,
                        organism=organism,
                        max_workers=workers,
                        max_metabolite_workers=workers,  # Use GUI setting instead of hardcoded value
                        cancel_flag=check_cancelled  # Pass cancel check function
                    )
                    
                    log_message(f"🔄 Starting pathway mapping for all metabolites...\n")
                    log_message(f"   • Annotating all {len(df)} metabolites\n")
                    log_message(f"   • No p-value filtering (annotation only)\n\n")
                    
                    # Run the pathway mapping WITHOUT filtering (all metabolites annotated)
                    # P-value threshold NOT passed - annotation doesn't need it
                    # (Network Analysis tab will use p-values for statistics)
                    annotator.run_pathway_mapping_with_dataframe(df, pvalue_threshold=None)
                    
                    # Check if cancelled after annotation
                    if getattr(self, 'pathway_annotation_cancelled', False):
                        self.root.after(0, lambda: self._handle_cancellation())
                        return
                    
                    # Read the annotated results
                    self.root.after(0, lambda: self.update_pathway_progress_with_percentage(85, "Loading annotated results..."))
                    annotated_df = pd.read_excel(temp_output_path)
                    
                    # CRITICAL FIX: Convert All_Pathways list to semicolon-separated string
                    # This must happen BEFORE Fisher ORA to prevent rebuilding from dirty individual columns
                    if 'All_Pathways' in annotated_df.columns:
                        def _list_to_string(pathways):
                            import ast
                            if isinstance(pathways, list):
                                return ';'.join(pathways) if pathways else ''
                            elif isinstance(pathways, str):
                                # Check if it's a string representation of a list (from Excel)
                                if pathways.startswith('[') and pathways.endswith(']'):
                                    try:
                                        # Parse string representation back to list
                                        pathway_list = ast.literal_eval(pathways)
                                        if isinstance(pathway_list, list):
                                            return ';'.join(str(p) for p in pathway_list if p) if pathway_list else ''
                                    except:
                                        pass
                                # Already a semicolon-separated string
                                return pathways
                            else:
                                return ''
                        annotated_df['All_Pathways'] = annotated_df['All_Pathways'].apply(_list_to_string)
                        log_message(f"✓ Converted All_Pathways from list to string format ({len(annotated_df)} metabolites)\n")
                    else:
                        log_message(f"⚠️ WARNING: All_Pathways column not found in annotated data!\n")
                        log_message(f"   Available columns: {list(annotated_df.columns)}\n")
                    
                    # Skip writing the intermediate 'pathway_annotated_*.xlsx' file
                    # This file duplicates content that will be saved later as a multi-sheet
                    # '<sheet>_pathways_*.xlsx'. Keeping only the consolidated export avoids confusion.
                    log_message("ℹ️ Skipped saving intermediate pathway_annotated_*.xlsx (redundant).\n")

                    # Clean up temporary files
                    try:
                        os.unlink(temp_input_path)
                        os.unlink(temp_output_path)
                    except:
                        pass  # Ignore cleanup errors
                
                except Exception as e:
                    log_message(f"\n❌ Pathway mapping FAILED: {str(e)}\n")
                    log_message(f"\n⛔ ANNOTATION ABORTED - Cannot continue without successful pathway mapping\n")
                    log_message(f"\nPlease check:\n")
                    log_message(f"   1. Your file has valid p-value and fold change columns\n")
                    log_message(f"   2. Use 'Verify Columns' to confirm column assignments\n")
                    log_message(f"   3. Check that numeric columns contain valid numbers\n")
                    import traceback
                    log_message(f"\nTechnical details:\n{traceback.format_exc()}\n")
                    self.root.after(0, lambda: self.update_pathway_progress_with_percentage(100, "Failed - Pathway mapping error"))
                    self.root.after(0, lambda: self.pathway_annotation_button.config(state='normal'))
                    return
            
            # Check cancellation before final steps
            if getattr(self, 'pathway_annotation_cancelled', False):
                self.root.after(0, lambda: self._handle_cancellation())
                return
            
            # Store the annotated data (thread-safe)
            self.root.after(0, lambda: self._store_annotation_results(annotated_df, annotated_df.copy()))
            
            # Save results to Excel file (after a short delay to keep UI responsive)
            self.root.after(1000, lambda: self._save_results_after_calculation(log_message))
            
            # Schedule network tab population to happen after saving/export
            self.root.after(1500, lambda: self._populate_network_tab_after_annotation(log_message))
            
            # Mode-aware terminology
            entity_name = "lipids" if mode == 'lipid' else "metabolites"
            log_message(f"\n✅ Pathway annotation complete!\n")
            log_message(f"📊 Processed {len(annotated_df)} {entity_name}\n")
            
            # Stop the timer
            self.root.after(0, lambda: setattr(self, 'pathway_start_time', None))
            
            self.root.after(0, lambda: self.update_pathway_progress_with_percentage(100, "✅ Process complete!"))
            log_message(f"\n✅ Ready to use pathway network analysis\n")
            
            # Enable Export button, disable Stop button, enable Annotate button
            self.root.after(0, lambda: self._enable_export_button())
            
            # Show completion message (thread-safe)
            self.root.after(0, lambda: messagebox.showinfo(
                "Annotation Complete", 
                f"✅ Pathway annotation complete!\n\n"
                f"📊 Results auto-loaded to Network Analysis tab\n\n"
                f"Next: Go to 'Network Analysis' tab to:\n"
                f"  1. Verify comparison columns (P-value, Log2FC)\n"
                f"  2. Configure analysis settings\n"
                f"  3. Generate pathway network"
            ))
            
        except Exception as e:
            def log_error(err_msg):
                self.root.after(0, lambda: (
                    self.pathway_progress_text.insert(tk.END, f"\n❌ {err_msg}\n"),
                    self.pathway_progress_text.see(tk.END)
                ))
            error_msg = f"Error during pathway annotation: {str(e)}"
            log_error(error_msg)
            
            # Re-enable buttons after error to allow restart
            self.root.after(0, lambda: self.pathway_annotation_button.config(state='normal'))
            if hasattr(self, 'pathway_stop_button'):
                self.root.after(0, lambda: self.pathway_stop_button.config(state='disabled'))
            if hasattr(self, 'export_pathway_results_button'):
                self.root.after(0, lambda: self.export_pathway_results_button.config(state='normal'))
            
            # Update progress indicators
            self.root.after(0, lambda: self.pathway_progress_bar.config(value=0))
            self.root.after(0, lambda: self.pathway_progress_percentage.config(text="Error"))
            
            # Show error dialog
            self.root.after(0, lambda: messagebox.showerror("Error", error_msg))
    
    def _store_annotation_results(self, filtered_data, original_data):
        """Thread-safe storage of annotation results"""
        # CRITICAL FIX: Convert All_Pathways list column to semicolon-separated string
        # This ensures Fisher ORA can use it without rebuilding from (potentially dirty) individual columns
        if 'All_Pathways' in filtered_data.columns:
            def _list_to_string(pathways):
                if isinstance(pathways, list):
                    return ';'.join(pathways) if pathways else ''
                elif isinstance(pathways, str):
                    return pathways  # Already a string
                else:
                    return ''
            
            filtered_data = filtered_data.copy()
            filtered_data['All_Pathways'] = filtered_data['All_Pathways'].apply(_list_to_string)
            logger.info(f"[STORAGE FIX] Converted All_Pathways from list to semicolon-separated string")
        
        # DEBUG: Check disease data before storage
        logger.info(f"[STORAGE DEBUG] Storing annotated data...")
        logger.info(f"[STORAGE DEBUG] Shape: {filtered_data.shape}")
        logger.info(f"[STORAGE DEBUG] Has 'Associated_Diseases': {'Associated_Diseases' in filtered_data.columns}")
        logger.info(f"[STORAGE DEBUG] Has 'Reactome_Disease': {'Reactome_Disease' in filtered_data.columns}")
        
        if 'Associated_Diseases' in filtered_data.columns:
            non_empty = filtered_data['Associated_Diseases'].fillna('').astype(str).str.strip()
            non_empty = non_empty[non_empty != '']
            non_empty = non_empty[non_empty.str.lower() != 'nan']
            logger.info(f"[STORAGE DEBUG] Rows with non-empty Associated_Diseases: {len(non_empty)}/{len(filtered_data)}")
            if len(non_empty) > 0:
                logger.info(f"[STORAGE DEBUG] Sample diseases (first 2): {list(non_empty.head(2))}")
        
        self.pathway_filtered_metabolites_data = filtered_data
        self.pathway_original_metabolites_data = original_data
        
        # (Removed noisy debug logging of All_Pathways_Display samples)

    def _run_fisher_ora_calculation(self, progress_callback, log_message):
        """Run Fisher ORA pathway statistics calculation.
        
        Based on the reference implementation from metabolite_annotation_gui001_1.py lines 13284-13380
        """
        from main_script.fisher_ora_pathway_analysis import calculate_pathway_statistics_fisher_ora
        from main_script.metabolite_pathways_annotator import MetabolitePathwayMapper
        import pandas as pd
        import tempfile
        import os
        
        log_message("ℹ️ NOTE: Fisher ORA relies on universe_config.json. After downloading or refreshing databases, run calculate_universe.py to keep the background universe accurate.\n")

        # Get GUI settings for Fisher ORA
        fisher_min_metabolites = int(self.pathway_min_metabolites.get())
        fisher_pvalue_threshold = float(self.pathway_pvalue_threshold.get())
        fisher_z_threshold = float(self.pathway_z_threshold.get())
        apply_fdr = self.pathway_apply_fdr.get()
        
        # Get pathway limiting settings
        fisher_max_total = None if not self.pathway_limit_pathways.get() else int(self.pathway_max_total.get())
        
        # Get selected organism/species
        selected_species = self.pathway_organism_var.get() if hasattr(self, 'pathway_organism_var') else 'Homo sapiens'
        
        log_message(f"✓ Fisher ORA Parameters:\n")
        log_message(f"  • Species: {selected_species}\n")
        log_message(f"  • Min metabolites: {fisher_min_metabolites}\n")
        log_message(f"  • P-value threshold (metabolite): {fisher_pvalue_threshold}\n")
        log_message(f"  • Z-score threshold: {fisher_z_threshold}\n")
        log_message(f"  • Apply FDR correction: {apply_fdr}\n")
        log_message(f"  • Max total pathways: {fisher_max_total if fisher_max_total else 'Unlimited'}\n\n")
        
        # Prepare data - Create 'All_Pathways' column if missing
        df_for_fisher = self.pathway_filtered_metabolites_data.copy()
        
        # DEBUG: Log available columns to diagnose verification issues
        log_message(f"  🔍 DEBUG: Available columns in data ({len(df_for_fisher.columns)} total):\n")
        log_message(f"    {list(df_for_fisher.columns)}\n\n")
        
        # ========================================================================
        # CRITICAL: Column Standardization Layer
        # ALWAYS use verified column assignments from user - NO HARDCODING, NO DEFAULTS
        # ========================================================================
        
        log_message(f"  📋 Standardizing column names for internal processing...\n")
        
        # Get verified column assignments - REQUIRED
        verified = getattr(self, 'verified_assignments', {}) or {}
        log_message(f"    📝 Verified assignments: {verified}\n")
        
        if not verified:
            raise ValueError("No verified column assignments found. Please run 'Verify Columns' first.")
        
        # Get user's verified column names
        name_col = verified.get('Feature ID') or verified.get('LipidID')
        pval_col = verified.get('P-Value')
        fc_col = verified.get('Log2 Fold Change')
        
        # CRITICAL: The annotator may have already standardized column names during annotation
        # Check if columns were renamed to standardized names ('Name', 'pvalue', 'log2FC')
        # This happens at line 3176 where verified feature column is renamed to 'Name'
        if name_col and name_col not in df_for_fisher.columns and 'Name' in df_for_fisher.columns:
            log_message(f"    ℹ️ Verified Feature ID '{name_col}' was already standardized to 'Name'\n")
            name_col = 'Name'
        if pval_col and pval_col not in df_for_fisher.columns and 'pvalue' in df_for_fisher.columns:
            log_message(f"    ℹ️ Verified P-Value '{pval_col}' was already standardized to 'pvalue'\n")
            pval_col = 'pvalue'
        if fc_col and fc_col not in df_for_fisher.columns and 'log2FC' in df_for_fisher.columns:
            log_message(f"    ℹ️ Verified Log2FC '{fc_col}' was already standardized to 'log2FC'\n")
            fc_col = 'log2FC'
        
        # Validate all required columns are assigned
        if not name_col:
            raise ValueError("Feature ID / LipidID not assigned. Please run 'Verify Columns' and assign the ID column.")
        if not pval_col:
            raise ValueError("P-Value column not assigned. Please run 'Verify Columns' and assign the P-Value column.")
        if not fc_col:
            raise ValueError("Log2 Fold Change column not assigned. Please run 'Verify Columns' and assign the Log2 Fold Change column.")
        
        # Create standardized columns from user's verified assignments
        # 1. Name column
        if name_col in df_for_fisher.columns:
            if 'Name' not in df_for_fisher.columns or name_col != 'Name':
                df_for_fisher['Name'] = df_for_fisher[name_col]
            log_message(f"    ✓ Name: '{name_col}'\n")
        else:
            raise ValueError(f"Verified Name column '{name_col}' not found in data!")
        
        # 2. pvalue column - CRITICAL for pathway statistics
        if pval_col in df_for_fisher.columns:
            if 'pvalue' not in df_for_fisher.columns or pval_col != 'pvalue':
                df_for_fisher['pvalue'] = df_for_fisher[pval_col]
            log_message(f"    ✓ pvalue: '{pval_col}' → Range: [{df_for_fisher['pvalue'].min():.6f}, {df_for_fisher['pvalue'].max():.6f}]\n")
        else:
            raise ValueError(f"Verified P-Value column '{pval_col}' not found in data!")
        
        # 3. log2FC column - CRITICAL for pathway direction
        if fc_col in df_for_fisher.columns:
            if 'log2FC' not in df_for_fisher.columns or fc_col != 'log2FC':
                df_for_fisher['log2FC'] = df_for_fisher[fc_col]
            log_message(f"    ✓ log2FC: '{fc_col}' → Range: [{df_for_fisher['log2FC'].min():.4f}, {df_for_fisher['log2FC'].max():.4f}]\n")
        else:
            raise ValueError(f"Verified Log2 Fold Change column '{fc_col}' not found in data!")
        
        log_message(f"  ✅ Column standardization complete (using verified assignments)\n\n")
        
        # If pre-annotated 'All_Pathways' exists, clean semicolon-fragmented entries first
        try:
            from main_script.metabolite_pathways_annotator import MetabolitePathwayMapper
            if 'All_Pathways' in df_for_fisher.columns:
                cleaned_count = 0
                examples = []
                for idx, val in df_for_fisher['All_Pathways'].fillna('').astype(str).items():
                    if not val.strip():
                        continue
                    new_val = MetabolitePathwayMapper.group_semicolon_fragments(val)
                    if new_val != val:
                        cleaned_count += 1
                        if len(examples) < 5:
                            examples.append((val, new_val))
                        df_for_fisher.at[idx, 'All_Pathways'] = new_val
                if cleaned_count > 0:
                    log_message(f"🔧 Cleaned {cleaned_count} pre-existing 'All_Pathways' entries by grouping semicolon fragments\n")
                    for old, new in examples:
                        log_message(f"   • Example: '{old}' → '{new}'\n")
        except Exception as _e:
            log_message(f"⚠️ Warning: Could not run All_Pathways cleaning: {_e}\n")

        if 'All_Pathways' not in df_for_fisher.columns:
            log_message(f"  Creating 'All_Pathways' column from available pathway sources...\n")
            pathway_cols = [col for col in df_for_fisher.columns if 'pathway' in col.lower() or 'kegg' in col.lower()]

            # Sanitize individual pathway columns to group semicolon-fragmented names
            try:
                from main_script.metabolite_pathways_annotator import MetabolitePathwayMapper
                cleaned_count = 0
                examples = []
                for pcol in pathway_cols:
                    if pcol in df_for_fisher.columns:
                        def _maybe_group(v):
                            try:
                                if pd.isna(v):
                                    return v
                                s = str(v)
                                if ';' in s and len(s) > 3:
                                    new = MetabolitePathwayMapper.group_semicolon_fragments(s)
                                    return new
                            except Exception:
                                pass
                            return v
                        before_sample = df_for_fisher[pcol].astype(str).head(5).tolist()
                        df_for_fisher[pcol] = df_for_fisher[pcol].apply(_maybe_group)
                        # Quick detect changes in column (count non-equal samples)
                        after_sample = df_for_fisher[pcol].astype(str).head(5).tolist()
                        for b, a in zip(before_sample, after_sample):
                            if b != a:
                                cleaned_count += 1
                                if len(examples) < 5:
                                    examples.append((b, a))
                if cleaned_count > 0:
                    log_message(f"🔧 Sanitized {cleaned_count} pathway column entries by grouping semicolon fragments before building All_Pathways\n")
                    for old, new in examples:
                        log_message(f"   • Example: '{old}' → '{new}'\n")
            except Exception as _e:
                log_message(f"⚠️ Warning: Could not sanitize pathway columns before All_Pathways creation: {_e}\n")

            def combine_pathways(row):
                all_pathways = set()
                for col in pathway_cols:
                    if pd.notna(row[col]):
                        val = str(row[col])
                        if val and val.lower() not in ['nan', 'none', '']:
                            # Split by common delimiters
                            pathways = [p.strip() for p in val.replace(';', '|').replace(',', '|').split('|')]
                            all_pathways.update([p for p in pathways if p])
                # IMPORTANT: Use semicolon as delimiter (backend expects semicolons, not pipes!)
                return ';'.join(sorted(all_pathways)) if all_pathways else ''
            
            df_for_fisher['All_Pathways'] = df_for_fisher.apply(combine_pathways, axis=1)
            # Normalize and clean pathway strings: remove numeric-only/too-short, dedupe case variants
            import re
            def _normalize_pathway_name(name: str) -> str:
                n = str(name).strip().strip(" .;:,\t")
                if not n:
                    return ''
                # Remove too-short or numeric-like names (e.g., '5', '1.')
                if len(n) <= 3:
                    return ''
                if re.fullmatch(r"\d+(?:\.\d+)?", n):
                    return ''
                # Standardize case (Title Case) while keeping small words lowercase
                small = {"and", "or", "of", "to", "in", "on", "the", "for"}
                tokens = re.split(r"(\s+|-)", n)
                rebuilt = []
                for t in tokens:
                    if t in {'-', ' ', '\t'} or t.isspace():
                        rebuilt.append(t)
                    else:
                        tt = t.lower()
                        rebuilt.append(tt if tt in small else tt.capitalize())
                return ''.join(rebuilt).strip()

            def _clean_all_paths(val: str) -> str:
                if not isinstance(val, str) or not val.strip():
                    return ''
                parts = [p.strip() for p in re.split(r"[;|,]", val) if p and p.strip()]
                seen = set()
                out = []
                for p in parts:
                    norm = _normalize_pathway_name(p)
                    if not norm:
                        continue
                    key = re.sub(r"[^a-z0-9]+", "", norm.lower())
                    if key and key not in seen:
                        seen.add(key)
                        out.append(norm)
                return ';'.join(sorted(out))

            df_for_fisher['All_Pathways'] = df_for_fisher['All_Pathways'].apply(_clean_all_paths)
            log_message(f"  ✓ 'All_Pathways' column created and normalized\n\n")
        
        # ===== CRITICAL: Reclassify diseases BEFORE Fisher ORA =====
        log_message(f"\n🔬 Pre-processing pathways before Fisher ORA...\n")
        log_message(f"📊 Reclassifying disease pathways...\n")
        
        try:
            # Import disease filtering functions
            from main_script.metabolite_pathways_annotator import apply_annotation_filters, is_disease_pathway
            
            # Define reclassification function (same as skip_annotation path)
            def reclassify_for_stats(row):
                """Apply disease filtering to pathway columns"""
                import re
                all_pathways = []
                all_diseases = []
                
                pathway_cols = ['HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 
                               'WikiPathways', 'Metabolika_Pathways', 'Reactome_Pathways', 'KEGG_Pathways']
                disease_cols = ['Reactome_Disease']
                
                # First, collect diseases from disease-specific columns
                for col in disease_cols:
                    if col in row.index and pd.notna(row[col]):
                        val = str(row[col])
                        if val and val.lower() not in ['nan', 'none', '']:
                            # Split by both ; and | delimiters
                            diseases = [d.strip() for d in re.split(r'[;|]', val) if d.strip()]
                            all_diseases.extend(diseases)
                
                # Then process pathway columns
                for col in pathway_cols:
                    if col in row.index and pd.notna(row[col]):
                        val = str(row[col])
                        if val and val.lower() not in ['nan', 'none', '']:
                            # Split by both ; and | delimiters (HMDB uses |, others use ;)
                            pathways = [p.strip() for p in re.split(r'[;|]', val) if p.strip()]
                            # Apply disease filtering
                            clean_pathways, diseases = apply_annotation_filters(pathways)
                            all_pathways.extend(clean_pathways)
                            all_diseases.extend(diseases)
                
                # Deduplicate
                unique_pathways = []
                seen_pathways = set()
                duplicates_found = []
                for p in all_pathways:
                    p_lower = p.lower()
                    if p_lower not in seen_pathways:
                        seen_pathways.add(p_lower)
                        unique_pathways.append(p)
                    else:
                        duplicates_found.append(p)
                
                unique_diseases = []
                seen_diseases = set()
                for d in all_diseases:
                    d_lower = d.lower()
                    if d_lower not in seen_diseases:
                        seen_diseases.add(d_lower)
                        unique_diseases.append(d)
                
                # Use proper delimiters: ; for pathways, | for diseases
                pathways_str = ';'.join(unique_pathways) if unique_pathways else ''
                diseases_str = ' | '.join(unique_diseases) if unique_diseases else ''
                
                # Return debug info
                debug_info = {
                    'metabolite': row.get('Name', 'Unknown'),
                    'total_collected': len(all_pathways),
                    'after_dedup': len(unique_pathways),
                    'duplicates_removed': len(duplicates_found)
                }
                
                return pathways_str, diseases_str, debug_info
            
            # Apply reclassification
            result = df_for_fisher.apply(reclassify_for_stats, axis=1, result_type='expand')
            df_for_fisher['All_Pathways'] = result[0]
            df_for_fisher['Associated_Diseases'] = result[1]
            debug_info_list = result[2].tolist()
            
            # DEBUG: Report deduplication
            total_dups = sum(d.get('duplicates_removed', 0) for d in debug_info_list)
            if total_dups > 0:
                log_message(f"  🔍 Removed {total_dups} duplicate pathways during Fisher ORA prep\n")
            
            # Count results
            diseases_count = (df_for_fisher['Associated_Diseases'].fillna('').astype(str).str.strip() != '').sum()
            pathways_count = (df_for_fisher['All_Pathways'].fillna('').astype(str).str.strip() != '').sum()
            
            log_message(f"  ✅ Disease reclassification complete:\n")
            log_message(f"     • {diseases_count} metabolites with diseases\n")
            log_message(f"     • {pathways_count} metabolites with pathways\n\n")
                
        except Exception as e:
            log_message(f"  ⚠️ Warning: Disease reclassification failed: {e}\n")
            log_message(f"  Continuing with original pathways...\n\n")
        
        # ===== Apply min_metabolites filter to reduce workload =====
        log_message(f"📊 Pre-filtering pathways by size (min: {fisher_min_metabolites})...\n")
        
        # Build pathway database from All_Pathways and merge duplicates by normalized key
        pathway_database_norm = {}   # norm_key -> set(metabolite_names)
        orig_name_map = {}           # norm_key -> {orig_name: set(metabolite_names)}
        for idx, row in df_for_fisher.iterrows():
            metabolite_name = row.get('Name', '')
            pathways = row.get('All_Pathways', [])

            if isinstance(pathways, str):
                if pathways.strip():
                    pathways = [p.strip() for p in pathways.split(';') if p.strip()]
                else:
                    pathways = []
            elif not isinstance(pathways, list):
                pathways = []

            for pathway in pathways:
                # Normalized key: lowercase, remove non-alphanumeric
                norm = re.sub(r"[^a-z0-9]+", "", pathway.lower())
                if not norm:
                    continue
                pathway_database_norm.setdefault(norm, set()).add(metabolite_name)
                orig_name_map.setdefault(norm, {}).setdefault(pathway, set()).add(metabolite_name)

        # Create canonical pathway database by selecting the original name variant with most metabolites
        pathway_database = {}
        duplicates_merged = 0
        merge_examples = []
        for norm_key, metab_set in pathway_database_norm.items():
            names_map = orig_name_map.get(norm_key, {})
            if names_map:
                # Choose canonical original name as the variant with largest metabolite set
                canonical_name = max(names_map.items(), key=lambda kv: len(kv[1]))[0]
            else:
                canonical_name = norm_key

            # Union of all metabolites for this normalized key
            pathway_database[canonical_name] = list(metab_set)

            if len(names_map) > 1:
                duplicates_merged += 1
                if len(merge_examples) < 5:
                    merge_examples.append((list(names_map.keys()), canonical_name))
        
        # Count before filtering (after merging duplicates)
        total_pathways_before_filter = len(pathway_database)
        if duplicates_merged > 0:
            log_message(f"🔁 Merged {duplicates_merged} duplicate pathway name groups (normalized by removing punctuation and lowercasing)\n")
            for orig_list, chosen in merge_examples:
                log_message(f"   • Variants: {orig_list} → Kept: '{chosen}'\n")
        
        # Filter by min_metabolites
        pathway_database_filtered = {
            pathway: metabolites 
            for pathway, metabolites in pathway_database.items() 
            if len(metabolites) >= fisher_min_metabolites
        }
        
        total_pathways_after_filter = len(pathway_database_filtered)
        pathways_filtered_out = total_pathways_before_filter - total_pathways_after_filter
        
        log_message(f"  ✅ Pre-filtering complete:\n")
        log_message(f"     • Total pathways: {total_pathways_before_filter}\n")
        log_message(f"     • Pathways meeting size threshold: {total_pathways_after_filter}\n")
        log_message(f"     • Pathways filtered out: {pathways_filtered_out}\n")
        log_message(f"     • Workload reduced by {pathways_filtered_out} pathways\n\n")
        
        # Pass pre-filtered pathway database to Fisher ORA to avoid re-building
        
        # Call Fisher ORA calculation with pre-filtered database
        pathway_stats = calculate_pathway_statistics_fisher_ora(
            df_for_fisher,
            min_metabolites=fisher_min_metabolites,
            pvalue_threshold=fisher_pvalue_threshold,
            z_threshold=fisher_z_threshold,
            fdr_method='fdr_bh' if apply_fdr else None,
            max_total_pathways=fisher_max_total,
            species=selected_species,
            progress_callback=progress_callback,
            pathway_database=pathway_database_filtered  # Use pre-filtered database
        )
        
        # Field name mapping for compatibility (Fisher ORA uses different field names)
        if pathway_stats:
            for pathway_name, stats in pathway_stats.items():
                if isinstance(stats, dict):
                    # CRITICAL: Keep BOTH 'hits' (significant) and 'metabolites' (all measured)
                    # - 'hits': significant metabolites only (used for network cascade selection)
                    # - 'metabolites': all measured metabolites (used for stats display)
                    # Network tab uses 'hits' for cascade selection, stats display uses 'metabolites'
                    # DO NOT overwrite one with the other!
                    
                    # Map direction_zscore → z_score (always use direction_zscore if present)
                    if 'direction_zscore' in stats:
                        stats['z_score'] = stats['direction_zscore']
                    
                    # Map direction → status (always use direction if present)
                    if 'direction' in stats:
                        stats['status'] = stats['direction']
                    
                    # Map fisher_pvalue/fdr_qvalue → combined_pvalue/adjusted_pvalue
                    if apply_fdr and 'fdr_qvalue' in stats:
                        stats['combined_pvalue'] = stats['fdr_qvalue']
                        stats['adjusted_pvalue'] = stats['fdr_qvalue']
                    elif 'fisher_pvalue' in stats:
                        stats['combined_pvalue'] = stats['fisher_pvalue']
                        stats['adjusted_pvalue'] = stats['fisher_pvalue']
        
        # Optionally augment with ML network score
        try:
            if pathway_stats:
                from main_script.ml_pathway_models import augment_with_ml_scores
                augmented = augment_with_ml_scores(df_for_fisher, pathway_stats)
                pathway_stats = augmented
                log_message("  ✓ Added ML network scores (RWR, centrality) to results\n")
        except Exception as _e:
            log_message(f"  ⚠️ Skipped ML augmentation: {_e}\n")
        
        # RE-RANK pathways after ML scores are added using the full ranking criteria:
        # 1. |Z-score| (highest first)
        # 2. Number of metabolites (highest first)
        # 3. P-value (lowest first)
        # 4. ML score (lowest/best first) - consensus score
        # 5. Pathway name (alphabetical - deterministic tiebreaker)
        if pathway_stats and fisher_max_total:
            pathway_list = []
            for pathway_name, stats in pathway_stats.items():
                stats['_pathway_name'] = pathway_name
                pathway_list.append(stats)
            
            # Sort by the full criteria
            pathway_list.sort(key=lambda r: (
                -abs(r.get('direction_zscore', r.get('z_score', 0))),  # 1. |Z-score| highest first
                -r.get('n_metabolites', len(r.get('metabolites', []))),  # 2. Most metabolites first
                r.get('fisher_pvalue', r.get('combined_pvalue', 1.0)),  # 3. Lowest p-value first
                r.get('ml_consensus_score', r.get('ml_rwr_score', float('inf'))),  # 4. Lowest ML score first (best)
                r.get('_pathway_name', '')  # 5. Alphabetical tiebreaker
            ))
            
            # Apply limit
            if len(pathway_list) > fisher_max_total:
                pathway_list = pathway_list[:fisher_max_total]
                log_message(f"  ✓ Re-ranked and limited to top {fisher_max_total} pathways\n")
            
            # Convert back to dict
            pathway_stats = {}
            for stats in pathway_list:
                pathway_name = stats.pop('_pathway_name')
                pathway_stats[pathway_name] = stats

        # Store FDR status for table headers
        self.pathway_fdr_enabled = apply_fdr
        
        # Save settings to data_manager for network tab to use
        self.data_manager.set_value('pathway_settings', {
            'min_metabolites': fisher_min_metabolites,
            'pvalue_threshold': fisher_pvalue_threshold,
            'z_threshold': fisher_z_threshold,
            'apply_fdr': apply_fdr,
            'max_total': fisher_max_total,
            'species': selected_species
        })
        
        return pathway_stats
    
    def _run_iwpa_calculation(self, progress_callback, log_message):
        """Run IWPA pathway statistics calculation.
        
        Based on the reference implementation from metabolite_annotation_gui001_1.py lines 13404-13490
        """
        from main_script.fisher_ora_pathway_analysis import calculate_pathway_statistics_iwpa
        import pandas as pd
        
        # Get GUI settings for IWPA
        weight_mode = self.iwpa_weight_mode.get()
        z_threshold = float(self.pathway_z_threshold.get())
        min_metabolites = int(self.pathway_min_metabolites.get())
        apply_fdr = self.pathway_apply_fdr.get()
        
        # Get pathway limiting settings (must be defined before use)
        iwpa_max_total = None if not self.pathway_limit_pathways.get() else int(self.pathway_max_total.get())
        
        log_message(f"✓ IWPA Parameters:\n")
        log_message(f"  • Weight mode: {weight_mode}\n")
        log_message(f"  • Z-score threshold: {z_threshold}\n")
        log_message(f"  • Min metabolites: {min_metabolites}\n")
        log_message(f"  • Apply FDR correction: {apply_fdr}\n")
        log_message(f"  • Max pathways: {iwpa_max_total if iwpa_max_total else 'Unlimited'}\n\n")
        
        # Prepare data - Create 'All_Pathways' column if missing (like Fisher ORA)
        df_for_iwpa = self.pathway_filtered_metabolites_data.copy()

        # If pre-annotated 'All_Pathways' exists, clean semicolon-fragmented entries first
        try:
            from main_script.metabolite_pathways_annotator import MetabolitePathwayMapper
            if 'All_Pathways' in df_for_iwpa.columns:
                cleaned_count = 0
                examples = []
                for idx, val in df_for_iwpa['All_Pathways'].fillna('').astype(str).items():
                    if not val.strip():
                        continue
                    new_val = MetabolitePathwayMapper.group_semicolon_fragments(val)
                    if new_val != val:
                        cleaned_count += 1
                        if len(examples) < 5:
                            examples.append((val, new_val))
                        df_for_iwpa.at[idx, 'All_Pathways'] = new_val
                if cleaned_count > 0:
                    log_message(f"🔧 Cleaned {cleaned_count} pre-existing 'All_Pathways' entries by grouping semicolon fragments\n")
                    for old, new in examples:
                        log_message(f"   • Example: '{old}' → '{new}'\n")
        except Exception as _e:
            log_message(f"⚠️ Warning: Could not run All_Pathways cleaning: {_e}\n")
        
        # Map comparison-specific columns to generic names expected by IWPA
        # IWPA expects: 'log2FC', 'pvalue', 'adjusted_pvalue'
        # Data may have: 'Control_vs_TBI_log2FC', 'Control_vs_TBI_adj_p', etc.
        
        # REQUIREMENT: Use verified p-value / log2FC assignments only (no auto-detect).
        # Log2FC: prefer verified_log2fc_col; if not available but a FoldChange column exists, compute log2FC = log2(FoldChange).
        # P-value: must be the verified p-value column. If missing, abort.

        # LOG2FC handling
        if hasattr(self, 'verified_log2fc_col') and self.verified_log2fc_col:
            log2fc_col = self.verified_log2fc_col
            if log2fc_col in df_for_iwpa.columns:
                df_for_iwpa['log2FC'] = pd.to_numeric(df_for_iwpa[log2fc_col], errors='coerce')
                log_message(f"📊 ENFORCING verified Log2FC column: '{log2fc_col}'\n")
            else:
                log_message(f"❌ Verified Log2FC column '{log2fc_col}' not found in data.\n")
                log_message("⛔ IWPA ABORTED - please re-verify Log2 Fold Change column.\n")
                return
        else:
            # Try to compute from a FoldChange-like column if present (but only if user provided it)
            fc_candidates = [c for c in df_for_iwpa.columns if 'fold' in c.lower() and 'log' not in c.lower()]
            if fc_candidates:
                fc_col = fc_candidates[0]
                try:
                    fc_vals = pd.to_numeric(df_for_iwpa[fc_col], errors='coerce')
                    if fc_vals.isna().all() or (fc_vals <= 0).any():
                        raise ValueError('Invalid fold-change values')
                    df_for_iwpa['log2FC'] = np.log2(fc_vals)
                    log_message(f"📊 Computed 'log2FC' from fold-change column '{fc_col}'\n")
                except Exception:
                    log_message(f"❌ FoldChange column '{fc_col}' is present but invalid for log2 computation.\n")
                    log_message("⛔ IWPA ABORTED - please provide a valid Log2 Fold Change column.\n")
                    return
            else:
                log_message("❌ No verified Log2FC or FoldChange column available. Please verify columns.\n")
                log_message("⛔ IWPA ABORTED - missing log2FC.\n")
                return

        # PVALUE handling (must be verified)
        if hasattr(self, 'verified_pvalue_col') and self.verified_pvalue_col:
            pvalue_col = self.verified_pvalue_col
            if pvalue_col in df_for_iwpa.columns:
                df_for_iwpa['pvalue'] = pd.to_numeric(df_for_iwpa[pvalue_col], errors='coerce').fillna(1.0)
                log_message(f"📊 ENFORCING verified p-value column: '{pvalue_col}'\n")
            else:
                log_message(f"❌ Verified p-value column '{pvalue_col}' not found in data.\n")
                log_message("⛔ IWPA ABORTED - please re-verify p-value column.\n")
                return
        else:
            log_message("❌ No verified p-value column available for IWPA. Please verify p-value column.\n")
            log_message("⛔ IWPA ABORTED - missing p-value.\n")
            return
        
        if 'All_Pathways' not in df_for_iwpa.columns:
            log_message(f"  Creating 'All_Pathways' column from available pathway sources...\n")
            pathway_cols = [col for col in df_for_iwpa.columns if 'pathway' in col.lower() or 'kegg' in col.lower()]

            # Sanitize individual pathway columns to group semicolon-fragmented names
            try:
                from main_script.metabolite_pathways_annotator import MetabolitePathwayMapper
                cleaned_count = 0
                examples = []
                for pcol in pathway_cols:
                    if pcol in df_for_iwpa.columns:
                        def _maybe_group(v):
                            try:
                                if pd.isna(v):
                                    return v
                                s = str(v)
                                if ';' in s and len(s) > 3:
                                    new = MetabolitePathwayMapper.group_semicolon_fragments(s)
                                    return new
                            except Exception:
                                pass
                            return v
                        before_sample = df_for_iwpa[pcol].astype(str).head(5).tolist()
                        df_for_iwpa[pcol] = df_for_iwpa[pcol].apply(_maybe_group)
                        after_sample = df_for_iwpa[pcol].astype(str).head(5).tolist()
                        for b, a in zip(before_sample, after_sample):
                            if b != a:
                                cleaned_count += 1
                                if len(examples) < 5:
                                    examples.append((b, a))
                if cleaned_count > 0:
                    log_message(f"🔧 Sanitized {cleaned_count} pathway column entries by grouping semicolon fragments before building All_Pathways\n")
                    for old, new in examples:
                        log_message(f"   • Example: '{old}' → '{new}'\n")
            except Exception as _e:
                log_message(f"⚠️ Warning: Could not sanitize pathway columns before All_Pathways creation: {_e}\n")

            def combine_pathways(row):
                all_pathways = set()
                for col in pathway_cols:
                    if pd.notna(row[col]):
                        val = str(row[col])
                        if val and val.lower() not in ['nan', 'none', '']:
                            # Split by common delimiters
                            pathways = [p.strip() for p in val.replace(';', '|').replace(',', '|').split('|')]
                            all_pathways.update([p for p in pathways if p])
                # IMPORTANT: Use semicolon as delimiter (backend expects semicolons, not pipes!)
                return ';'.join(sorted(all_pathways)) if all_pathways else ''
            
            df_for_iwpa['All_Pathways'] = df_for_iwpa.apply(combine_pathways, axis=1)
            # Apply same normalization/cleanup as Fisher ORA
            import re
            def _normalize_pathway_name(name: str) -> str:
                n = str(name).strip().strip(" .;:,\t")
                if not n:
                    return ''
                if len(n) <= 3:
                    return ''
                if re.fullmatch(r"\d+(?:\.\d+)?", n):
                    return ''
                small = {"and", "or", "of", "to", "in", "on", "the", "for"}
                tokens = re.split(r"(\s+|-)", n)
                rebuilt = []
                for t in tokens:
                    if t in {'-', ' ', '\t'} or t.isspace():
                        rebuilt.append(t)
                    else:
                        tt = t.lower()
                        rebuilt.append(tt if tt in small else tt.capitalize())
                return ''.join(rebuilt).strip()

            def _clean_all_paths(val: str) -> str:
                if not isinstance(val, str) or not val.strip():
                    return ''
                parts = [p.strip() for p in re.split(r"[;|,]", val) if p and p.strip()]
                seen = set()
                out = []
                for p in parts:
                    norm = _normalize_pathway_name(p)
                    if not norm:
                        continue
                    key = re.sub(r"[^a-z0-9]+", "", norm.lower())
                    if key and key not in seen:
                        seen.add(key)
                        out.append(norm)
                return ';'.join(sorted(out))

            df_for_iwpa['All_Pathways'] = df_for_iwpa['All_Pathways'].apply(_clean_all_paths)
            log_message(f"  ✓ 'All_Pathways' column created and normalized from {len(pathway_cols)} sources\n\n")
        
        # ===== CRITICAL: Reclassify diseases BEFORE IWPA =====
        log_message(f"\n🔬 Pre-processing pathways before IWPA...\n")
        log_message(f"📊 Reclassifying disease pathways...\n")
        
        try:
            # Import disease filtering functions
            from main_script.metabolite_pathways_annotator import apply_annotation_filters, is_disease_pathway
            
            # Define reclassification function (same as skip_annotation path)
            def reclassify_for_stats(row):
                """Apply disease filtering to pathway columns"""
                import re
                all_pathways = []
                all_diseases = []
                
                pathway_cols = ['HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 
                               'WikiPathways', 'Metabolika_Pathways', 'Reactome_Pathways', 'KEGG_Pathways']
                disease_cols = ['Reactome_Disease']
                
                # First, collect diseases from disease-specific columns
                for col in disease_cols:
                    if col in row.index and pd.notna(row[col]):
                        val = str(row[col])
                        if val and val.lower() not in ['nan', 'none', '']:
                            # Split by both ; and | delimiters
                            diseases = [d.strip() for d in re.split(r'[;|]', val) if d.strip()]
                            all_diseases.extend(diseases)
                
                # Then process pathway columns
                for col in pathway_cols:
                    if col in row.index and pd.notna(row[col]):
                        val = str(row[col])
                        if val and val.lower() not in ['nan', 'none', '']:
                            # Split by both ; and | delimiters (HMDB uses |, others use ;)
                            pathways = [p.strip() for p in re.split(r'[;|]', val) if p.strip()]
                            # Apply disease filtering
                            clean_pathways, diseases = apply_annotation_filters(pathways)
                            all_pathways.extend(clean_pathways)
                            all_diseases.extend(diseases)
                
                # Deduplicate
                unique_pathways = []
                seen_pathways = set()
                duplicates_found = []
                for p in all_pathways:
                    p_lower = p.lower()
                    if p_lower not in seen_pathways:
                        seen_pathways.add(p_lower)
                        unique_pathways.append(p)
                    else:
                        duplicates_found.append(p)
                
                unique_diseases = []
                seen_diseases = set()
                for d in all_diseases:
                    d_lower = d.lower()
                    if d_lower not in seen_diseases:
                        seen_diseases.add(d_lower)
                        unique_diseases.append(d)
                
                # Use proper delimiters: ; for pathways, | for diseases
                pathways_str = ';'.join(unique_pathways) if unique_pathways else ''
                diseases_str = ' | '.join(unique_diseases) if unique_diseases else ''
                
                # Return debug info
                debug_info = {
                    'metabolite': row.get('Name', 'Unknown'),
                    'total_collected': len(all_pathways),
                    'after_dedup': len(unique_pathways),
                    'duplicates_removed': len(duplicates_found)
                }
                
                return pathways_str, diseases_str, debug_info
            
            # Apply reclassification
            result = df_for_iwpa.apply(reclassify_for_stats, axis=1, result_type='expand')
            df_for_iwpa['All_Pathways'] = result[0]
            df_for_iwpa['Associated_Diseases'] = result[1]
            debug_info_list = result[2].tolist()
            
            # DEBUG: Report deduplication
            total_dups = sum(d.get('duplicates_removed', 0) for d in debug_info_list)
            if total_dups > 0:
                log_message(f"  🔍 Removed {total_dups} duplicate pathways during IWPA prep\n")
            
            # Count results
            diseases_count = (df_for_iwpa['Associated_Diseases'].fillna('').astype(str).str.strip() != '').sum()
            pathways_count = (df_for_iwpa['All_Pathways'].fillna('').astype(str).str.strip() != '').sum()
            
            log_message(f"  ✅ Disease reclassification complete:\n")
            log_message(f"     • {diseases_count} metabolites with diseases\n")
            log_message(f"     • {pathways_count} metabolites with pathways\n\n")
                
        except Exception as e:
            log_message(f"  ⚠️ Warning: Disease reclassification failed: {e}\n")
            log_message(f"  Continuing with original pathways...\n\n")
        
        # ===== Apply min_metabolites filter to reduce workload =====
        log_message(f"📊 Pre-filtering pathways by size (min: {min_metabolites})...\n")
        
        # Build pathway database from All_Pathways and merge duplicates by normalized key
        pathway_database_norm = {}
        orig_name_map = {}
        for idx, row in df_for_iwpa.iterrows():
            metabolite_name = row.get('Name', '')
            pathways = row.get('All_Pathways', [])

            if isinstance(pathways, str):
                if pathways.strip():
                    pathways = [p.strip() for p in pathways.split(';') if p.strip()]
                else:
                    pathways = []
            elif not isinstance(pathways, list):
                pathways = []

            for pathway in pathways:
                norm = re.sub(r"[^a-z0-9]+", "", pathway.lower())
                if not norm:
                    continue
                pathway_database_norm.setdefault(norm, set()).add(metabolite_name)
                orig_name_map.setdefault(norm, {}).setdefault(pathway, set()).add(metabolite_name)

        pathway_database = {}
        duplicates_merged = 0
        merge_examples = []
        for norm_key, metab_set in pathway_database_norm.items():
            names_map = orig_name_map.get(norm_key, {})
            if names_map:
                canonical_name = max(names_map.items(), key=lambda kv: len(kv[1]))[0]
            else:
                canonical_name = norm_key
            pathway_database[canonical_name] = list(metab_set)
            if len(names_map) > 1:
                duplicates_merged += 1
                if len(merge_examples) < 5:
                    merge_examples.append((list(names_map.keys()), canonical_name))
        
        # Count before filtering (after merging duplicates)
        total_pathways_before_filter = len(pathway_database)
        if duplicates_merged > 0:
            log_message(f"🔁 Merged {duplicates_merged} duplicate pathway name groups (normalized by removing punctuation and lowercasing)\n")
            for orig_list, chosen in merge_examples:
                log_message(f"   • Variants: {orig_list} → Kept: '{chosen}'\n")
        
        # Filter by min_metabolites
        pathway_database_filtered = {
            pathway: metabolites 
            for pathway, metabolites in pathway_database.items() 
            if len(metabolites) >= min_metabolites
        }
        
        total_pathways_after_filter = len(pathway_database_filtered)
        pathways_filtered_out = total_pathways_before_filter - total_pathways_after_filter
        
        log_message(f"  ✅ Pre-filtering complete:\n")
        log_message(f"     • Total pathways: {total_pathways_before_filter}\n")
        log_message(f"     • Pathways meeting size threshold: {total_pathways_after_filter}\n")
        log_message(f"     • Pathways filtered out: {pathways_filtered_out}\n")
        log_message(f"     • Workload reduced by {pathways_filtered_out} pathways\n\n")
        
        # Pass pre-filtered pathway database to IWPA to avoid re-building
        
        # Call IWPA calculation with pre-filtered database
        pathway_stats = calculate_pathway_statistics_iwpa(
            df_for_iwpa,
            min_metabolites=min_metabolites,
            weight_mode=weight_mode,
            z_threshold=z_threshold,
            max_total_pathways=iwpa_max_total,
            fdr_method='fdr_bh' if apply_fdr else None,
            progress_callback=progress_callback,
            pathway_database=pathway_database_filtered  # Use pre-filtered database
        )
        
        # Field name mapping for compatibility (IWPA uses different field names)
        if pathway_stats:
            for pathway_name, stats in pathway_stats.items():
                if isinstance(stats, dict):
                    # Map Z → z_score
                    if 'Z' in stats and 'z_score' not in stats:
                        stats['z_score'] = stats['Z']
                    
                    # Map pvalue → combined_pvalue/adjusted_pvalue
                    if apply_fdr and 'adjusted_pvalue' in stats:
                        stats['combined_pvalue'] = stats['adjusted_pvalue']
                    elif 'pvalue' in stats:
                        stats['combined_pvalue'] = stats['pvalue']
                        if 'adjusted_pvalue' not in stats:
                            stats['adjusted_pvalue'] = stats['pvalue']
                    
                    # Map mean_weight → mean_log2fc (if using log2fc mode)
                    if 'mean_weight' in stats and 'mean_log2fc' not in stats:
                        stats['mean_log2fc'] = stats['mean_weight']
                    if 'std_weight' in stats and 'std_log2fc' not in stats:
                        stats['std_log2fc'] = stats['std_weight']
        
        # Optionally augment with ML network score
        try:
            if pathway_stats:
                from main_script.ml_pathway_models import augment_with_ml_scores
                augmented = augment_with_ml_scores(df_for_iwpa, pathway_stats)
                pathway_stats = augmented
                log_message("  ✓ Added ML network scores (RWR, centrality) to results\n")
        except Exception as _e:
            log_message(f"  ⚠️ Skipped ML augmentation: {_e}\n")
        
        # RE-RANK pathways after ML scores are added using the full ranking criteria:
        # 1. |Z-score| (highest first)
        # 2. Number of metabolites (highest first)
        # 3. P-value (lowest first)
        # 4. ML score (lowest/best first) - consensus score
        # 5. Pathway name (alphabetical - deterministic tiebreaker)
        if pathway_stats and iwpa_max_total:
            pathway_list = []
            for pathway_name, stats in pathway_stats.items():
                stats['_pathway_name'] = pathway_name
                pathway_list.append(stats)
            
            # Sort by the full criteria
            pathway_list.sort(key=lambda r: (
                -abs(r.get('z_score', r.get('Z', 0))),  # 1. |Z-score| highest first
                -r.get('n_metabolites', len(r.get('metabolites', []))),  # 2. Most metabolites first
                r.get('combined_pvalue', r.get('pvalue', 1.0)),  # 3. Lowest p-value first
                r.get('ml_consensus_score', r.get('ml_rwr_score', float('inf'))),  # 4. Lowest ML score first (best)
                r.get('_pathway_name', '')  # 5. Alphabetical tiebreaker
            ))
            
            # Apply limit
            if len(pathway_list) > iwpa_max_total:
                pathway_list = pathway_list[:iwpa_max_total]
                log_message(f"  ✓ Re-ranked and limited to top {iwpa_max_total} pathways\n")
            
            # Convert back to dict
            pathway_stats = {}
            for stats in pathway_list:
                pathway_name = stats.pop('_pathway_name')
                pathway_stats[pathway_name] = stats

        # Store FDR status for table headers
        self.pathway_fdr_enabled = apply_fdr
        
        # Save settings to data_manager for network tab to use
        # Get selected organism/species
        selected_species = self.pathway_organism_var.get() if hasattr(self, 'pathway_organism_var') else 'Homo sapiens'
        
        # Get p-value threshold (for metabolite significance)
        iwpa_pvalue_threshold = 0.05  # IWPA doesn't use this parameter directly, but save for consistency
        
        self.data_manager.set_value('pathway_settings', {
            'min_metabolites': min_metabolites,
            'pvalue_threshold': iwpa_pvalue_threshold,
            'z_threshold': z_threshold,
            'apply_fdr': apply_fdr,
            'max_total': iwpa_max_total,
            'species': selected_species
        })
        
        return pathway_stats

    # def _calculate_pathway_statistics(self, annotated_df, log_message):
    #     """Calculate and store pathway statistics using Fisher ORA or IWPA method.
        
    #     Checks the user's method selection and calls the appropriate analysis method:
    #     - Fisher ORA: Fisher's Exact Test with Overrepresentation Analysis
    #     - IWPA: Integrated Weighted Pathway Analysis with z-scores
        
    #     Based on the reference implementation from metabolite_annotation_gui001_1.py lines 13261-13264
    #     """
    #     try:
    #         # Check if we have the required columns for pathway statistics
    #         # Use the verified columns from column assignment dialog (stored in instance variables)
    #         has_pvalue = bool(getattr(self, 'verified_pvalue_col', None))
    #         has_fc = bool(getattr(self, 'verified_log2fc_col', None))
            
    #         # # If missing required columns, skip statistics calculation
    #         # if not has_pvalue or not has_fc:
    #         #     log_message(f"\n⚠️ Skipping pathway statistics calculation:\n")
    #         #     if not has_pvalue:
    #         #         log_message(f"   ❌ Missing p-value column\n")
    #         #         log_message(f"      Please verify columns using 'Verify Columns' button\n")
    #         #     if not has_fc:
    #         #         log_message(f"   ❌ Missing fold change column\n")
    #         #         log_message(f"      Please verify columns using 'Verify Columns' button\n")
    #         #     log_message(f"\n   For lipid mode without statistical comparison:\n")
    #         #     log_message(f"   • ID annotation completed successfully\n")
    #         #     log_message(f"   • Pathway statistics require p-values and fold changes from statistical tests\n")
    #         #     log_message(f"   • To enable statistics: perform statistical analysis first, then load that file here\n\n")
                
    #             # Store empty statistics
    #             # self.pathway_original_pathways_data = {}
    #             # self.pathway_filtered_pathways_data = {}
    #             # self.current_pathway_stats = {}
                
    #             # log_message(f"✅ Analysis complete (statistics skipped)\n")
    #             # return
            
    #         # Get selected analysis method from UI
    #         selected_method = self.pathway_analysis_method.get()
            
    #         log_message(f"\n📊 Running {selected_method} pathway statistics calculation...\n")
            
    #         # Create progress callback for UI updates
    #         def progress_callback(percentage, message):
    #             self.root.after(0, lambda: self.update_pathway_progress_with_percentage(percentage, message))
    #             log_message(f"  {message}\n")
            
    #         # Call the appropriate method based on user selection
    #         if selected_method == "Fisher ORA":
    #             pathway_stats = self._run_fisher_ora_calculation(progress_callback, log_message)
    #         elif selected_method == "IWPA":
    #             pathway_stats = self._run_iwpa_calculation(progress_callback, log_message)
    #         else:
    #             raise ValueError(f"Unknown pathway analysis method: {selected_method}")
            
    #         # Ensure pathway_stats is a dict
    #         if not isinstance(pathway_stats, dict):
    #             log_message(f"⚠️ Unexpected return type {type(pathway_stats)}, using empty dict\n")
    #             pathway_stats = {}
            
    #         # Store the calculated statistics
    #         self.pathway_original_pathways_data = pathway_stats.copy() if pathway_stats else {}
    #         self.pathway_filtered_pathways_data = pathway_stats.copy() if pathway_stats else {}
            
    #         log_message(f"\n✅ Pathway statistics calculated successfully!\n")
    #         log_message(f"   Total pathways found: {len([p for p in pathway_stats.keys() if p != '_validation_report'])}\n")
            
    #         # Log pathway status distribution
    #         activated = 0
    #         inhibited = 0
    #         no_change = 0
    #         for pathway_name, stats in pathway_stats.items():
    #             if pathway_name == '_validation_report' or pathway_name == '_excluded_debug':
    #                 continue
    #             status = stats.get('status', 'No Change') if isinstance(stats, dict) else 'No Change'
    #             if status == 'Activated':
    #                 activated += 1
    #             elif status == 'Inhibited':
    #                 inhibited += 1
    #             else:
    #                 no_change += 1
            
    #         log_message(f"\n📈 PATHWAY STATUS DISTRIBUTION:\n")
    #         log_message(f"   🔺 Activated pathways: {activated}\n")
    #         log_message(f"   🔻 Inhibited pathways: {inhibited}\n")
    #         log_message(f"   ⚪ No change pathways: {no_change}\n")
            
    #         # Log top significant pathways by z-score
    #         log_message(f"\n🏆 TOP SIGNIFICANT PATHWAYS:\n")
    #         sorted_pathways = sorted(
    #             [(name, stats) for name, stats in pathway_stats.items() if name not in ['_validation_report', '_excluded_debug']],
    #             key=lambda x: abs(x[1].get('z_score', 0)) if isinstance(x[1], dict) else 0,
    #             reverse=True
    #         )[:5]
            
    #         for rank, (pathway_name, stats) in enumerate(sorted_pathways, 1):
    #             if not isinstance(stats, dict):
    #                 continue
    #             z_score = stats.get('z_score', 0)
    #             p_value = stats.get('combined_pvalue', 1.0)
    #             n_met = stats.get('n_metabolites', 0)
    #             status = stats.get('status', 'No Change')
    #             log_message(f"   {rank}. {pathway_name}\n")
    #             log_message(f"      Z-score: {z_score:.3f} | P-value: {p_value:.2e} | Metabolites: {n_met} | Status: {status}\n")
            
    #         log_message(f"\n💾 Pathways ready for network visualization\n")
            
    #     except Exception as e:
    #         error_msg = f"Error calculating pathway statistics: {str(e)}"
    #         log_message(f"\n❌ {error_msg}\n")
    #         import traceback
    #         traceback.print_exc()
    #         self.pathway_original_pathways_data = {}
    #         self.pathway_filtered_pathways_data = {}

    def stop_pathway_annotation(self):
        """Stop the current pathway annotation process"""
        # Check if annotation thread is running
        if hasattr(self, 'annotation_thread') and self.annotation_thread and self.annotation_thread.is_alive():
            # Set cancel flag
            self.pathway_annotation_cancelled = True
            
            # Update UI to show cancellation is in progress
            self.pathway_progress_text.insert(tk.END, "\n🛑 Stopping annotation process...\n")
            self.pathway_progress_text.see(tk.END)
            self.pathway_progress_percentage.config(text="Stopping...")
            
            messagebox.showinfo("Process Stopped", "Pathway annotation process is being stopped. Please wait...")
        else:
            messagebox.showinfo("No Process Running", "No pathway annotation process is currently running.")
    
    def _handle_cancellation(self):
        """Handle cleanup after cancellation"""
        # Reset progress indicators
        self.pathway_progress_bar.config(value=0)
        self.pathway_progress_percentage.config(text="Cancelled")
        self.pathway_progress_text.insert(tk.END, "\n❌ Process cancelled by user\n")
        self.pathway_progress_text.see(tk.END)
        
        # Reset timer
        self.pathway_start_time = None
        
        # Re-enable buttons
        if hasattr(self, 'pathway_stop_button'):
            self.pathway_stop_button.config(state='disabled')
        if hasattr(self, 'pathway_annotation_button'):
            self.pathway_annotation_button.config(state='normal')
        
        messagebox.showwarning("Cancelled", "Pathway annotation process was cancelled.")
    
    def _enable_export_button(self):
        """Enable Export button and disable Stop button after successful completion"""
        if hasattr(self, 'export_pathway_results_button'):
            self.export_pathway_results_button.config(state='normal')
        if hasattr(self, 'pathway_stop_button'):
            self.pathway_stop_button.config(state='disabled')
        if hasattr(self, 'pathway_annotation_button'):
            self.pathway_annotation_button.config(state='normal')


