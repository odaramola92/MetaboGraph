import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import logging
import os
import threading
import pandas as pd
import time
import traceback


# Import backend and shared components
from gui.shared.base_tab import BaseTab, _setup_global_styles
from gui.shared.utils import TimerManager, SessionManager, ProgressUpdater
from gui.shared.column_assignment import show_column_assignment_dialog
from main_script.metabolite_ID_annotator import MetaboliteIDAnnotator

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)



class IDAnnotationTab(BaseTab):
    """ID Annotation Tab - Placeholder for metabolite ID annotation
    
    This tab handles metabolite ID annotation using multiple databases.
    Data flows from Data Cleaning tab -> ID Annotation tab
    """
    
    def __init__(self, parent, data_manager):
        """Initialize ID Annotation tab"""
        super().__init__(parent, data_manager)
        
        # Setup global styles (runs only once)
        _setup_global_styles()
        
        # Get root window for dialogs
        self.root = parent.winfo_toplevel()
        
        # Initialize variables
        self._initialize_variables()
        
        # Create UI
        self.setup_ui()
        print("[OK] ID Annotation Tab initialized")
    
    def _initialize_variables(self):
        """Initialize all variables before UI creation"""
        # Input file (auto-loaded from data cleaning)
        self.id_input_file = tk.StringVar()
        self.lipid_id_annotation_file = tk.StringVar()
        
        # Output settings
        self.id_output_file = tk.StringVar()
        self.id_output_directory = tk.StringVar(value=os.getcwd())
        self.id_output_filename = tk.StringVar(value="id_annotated_metabolites.xlsx")
        self.lipid_output_file = tk.StringVar()
        self.lipid_output_directory = tk.StringVar(value=os.getcwd())
        
        # Processing options
        self.id_workers = tk.StringVar(value="2")
        self.lipid_workers = tk.StringVar(value="2")
        self.id_mode = tk.StringVar(value="metabolite")
        self.use_database_hmdb = tk.BooleanVar(value=True)
        self.use_database_lipidmaps = tk.BooleanVar(value=True)
        self.use_database_pathbank = tk.BooleanVar(value=True)
        
        # Progress tracking
        self.id_progress_var = tk.DoubleVar(value=0.0)
        self.id_status_message = tk.StringVar(value="Ready")
        
        # Timer tracking
        self.id_start_time = None
        self.id_current_step = 0
        self.id_total_steps = 100
        
        # Annotation cancellation flag
        self._stop_annotation = False

        # Runtime annotation settings (safe defaults for first run)
        self._selected_id_cols_runtime = []
        self._effective_id_filter_mode_runtime = False
        self._skip_id_filtering_runtime = False
        self._custom_id_search_mode = False
        self._annotation_mode_runtime = None
        self._annotation_workers_runtime = None
        
        # Pathway annotation variables
        self.pathway_input_file = tk.StringVar()
        self.pathway_sheet_combo = None
        self.pathway_sheet_var = tk.StringVar()
    
    def setup_ui(self):
        """Create the tab interface"""
        # Main container
        main_container = tk.Frame(self.frame, bg='#f0f0f0')
        main_container.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Create the ID annotation interface
        self.create_id_annotation_interface(main_container)
     # ID Annotation Tab methods       
    
    def create_id_annotation_interface(self, parent):
        """Create the ID annotation interface with tabs for Metabolites and Lipids"""
        # Main scrollable container (Canvas + inner frame) so the whole
        # ID Annotation tab can be scrolled when window is small.
        main_canvas = tk.Canvas(parent, bg='#f0f0f0', highlightthickness=0)
        main_scrollbar = ttk.Scrollbar(parent, orient='vertical', command=main_canvas.yview)
        main_scrollable = tk.Frame(main_canvas, bg='#f0f0f0')

        # Keep canvas scrollregion in sync with inner frame size
        main_scrollable.bind(
            "<Configure>",
            lambda e: main_canvas.configure(scrollregion=main_canvas.bbox("all"))
        )

        canvas_window = main_canvas.create_window((0, 0), window=main_scrollable, anchor='nw')
        main_canvas.configure(yscrollcommand=main_scrollbar.set)

        def _configure_main_canvas(event):
            try:
                main_canvas.itemconfig(canvas_window, width=event.width)
            except Exception:
                try:
                    main_canvas.itemconfig(canvas_window, width=event.width)
                except Exception:
                    pass
            try:
                main_canvas.configure(scrollregion=main_canvas.bbox('all'))
            except Exception:
                pass

        main_canvas.bind('<Configure>', _configure_main_canvas)
        main_canvas.pack(side='left', fill='both', expand=True, padx=10, pady=10)
        main_scrollbar.pack(side='right', fill='y')

        # CRITICAL: Ensure main_scrollable expands to fill available canvas height
        # This allows content_frame.grid_rowconfigure(0, weight=1) to work properly
        def _on_main_canvas_configure(event):
            # Get the height of the canvas viewport
            canvas_height = main_canvas.winfo_height()
            # Get current scroll region height
            scroll_height = main_canvas.bbox('all')[3] if main_canvas.bbox('all') else 0
            # If content is smaller than canvas, expand the scrollable frame to fill
            if scroll_height < canvas_height:
                main_scrollable.configure(height=canvas_height)
                main_canvas.itemconfig(canvas_window, height=canvas_height)
        
        main_canvas.bind('<Configure>', _on_main_canvas_configure, add=True)

        # Mousewheel support (Windows). Bind to both canvas and inner frame so
        # wheel events over child widgets still scroll the container.
        def _on_main_mousewheel(event):
            try:
                main_canvas.yview_scroll(int(-1*(event.delta/120)), 'units')
            except Exception:
                pass
        main_canvas.bind('<MouseWheel>', _on_main_mousewheel)
        main_scrollable.bind('<MouseWheel>', _on_main_mousewheel)
        
        # Create three-column layout with proper weights
        # Reduced weight for left/middle columns (0,1) to make them narrower
        content_frame = tk.Frame(main_scrollable, bg='#f0f0f0')
        content_frame.pack(fill='both', expand=True)
        content_frame.grid_columnconfigure(0, weight=2, minsize=150)
        content_frame.grid_columnconfigure(1, weight=2, minsize=150)
        content_frame.grid_columnconfigure(2, weight=3, minsize=200)
        content_frame.grid_rowconfigure(0, weight=1, minsize=500)  # Increased minimum height for better visibility

        # Left/Middle columns - Combined frame with notebook for Metabolites and Lipids tabs
        left_middle_frame = tk.LabelFrame(content_frame, 
                                text="⚙️ Configuration & Settings", 
                                font=('Arial', 12, 'bold'),
                                bg='#f0f0f0', 
                                fg='#2c3e50')
        left_middle_frame.grid(row=0, column=0, columnspan=2, sticky='nsew', padx=(0, 3))

        # Sub-notebook for Metabolite vs Lipid ID annotation
        self.id_notebook = ttk.Notebook(left_middle_frame, style='IdAnnotation.TNotebook')
        self.id_notebook.pack(fill='both', expand=True, padx=5, pady=5)

        metabolite_tab = tk.Frame(self.id_notebook, bg='#f0f0f0')
        lipid_tab = tk.Frame(self.id_notebook, bg='#f0f0f0')

        self.id_notebook.add(metabolite_tab, text="🧪 Metabolites")
        self.id_notebook.add(lipid_tab, text="🧬 Lipids")

        # Track which mode is active based on selected subtab
        def _on_id_tab_change(event):
            try:
                tab_text = event.widget.tab(event.widget.select(), "text")
                if "Lipids" in tab_text:
                    if hasattr(self, 'id_mode'):
                        self.id_mode.set('lipid')
                else:
                    if hasattr(self, 'id_mode'):
                        self.id_mode.set('metabolite')
            except Exception:
                pass

        self.id_notebook.bind('<<NotebookTabChanged>>', _on_id_tab_change)

        # Create interfaces for each subtab
        self.create_metabolite_id_annotation_interface(metabolite_tab)
        # Lipid-specific ID annotation interface (new subtab)
        self.create_lipid_id_annotation_interface(lipid_tab)
        
        # ========== RIGHT COLUMN - Progress and results (shared) ==========
        right_frame = tk.LabelFrame(content_frame, 
                                text="📋 Progress & Results", 
                                font=('Arial', 12, 'bold'),
                                bg='#f0f0f0', 
                                fg='#2c3e50')
        right_frame.grid(row=0, column=2, sticky='nsew', padx=(3, 0))
        
        # Control buttons AT THE TOP
        buttons_frame = tk.Frame(right_frame, bg='#f0f0f0')
        buttons_frame.pack(fill='x', padx=10, pady=(10, 10))
        
        # Buttons container frame for side by side layout
        buttons_container = tk.Frame(buttons_frame, bg='#f0f0f0')
        buttons_container.pack(fill='x', pady=5)
        
        # Start button
        self.id_start_button = tk.Button(buttons_container, 
                        text="🚀 Start ID Annotation", 
                        command=self.on_start_annotation_clicked,
                                        bg='#27ae60', 
                                        fg='white', 
                                        font=('Arial', 10, 'bold'),
                                        height=2,
                                        relief='raised',
                                        borderwidth=3)
        self.id_start_button.pack(side='left', fill='x', expand=True, padx=(0, 5))
        
        # Stop button
        self.id_stop_button = tk.Button(buttons_container, 
                                    text="⏹️ Stop", 
                                    command=self.stop_id_annotation,
                                    bg='#e74c3c', 
                                    fg='white', 
                                    font=('Arial', 10, 'bold'),
                                    height=2,
                                    relief='raised',
                                    borderwidth=3,
                                    state='disabled')
        self.id_stop_button.pack(side='left', fill='x', expand=True, padx=(5, 0))
        
        # Open Output Folder button
        self.id_open_output_folder_button = tk.Button(buttons_container,
                                                     text="📂 Open Folder",
                                                     command=self.open_id_output_folder,
                                                     bg='#9b59b6',
                                                     fg='white',
                                                     font=('Arial', 10, 'bold'),
                                                     height=2,
                                                     relief='raised',
                                                     borderwidth=3)
        self.id_open_output_folder_button.pack(side='left', fill='x', expand=True, padx=(5, 0))
        
        # Progress bar section
        progress_section = tk.Frame(right_frame, bg='#f0f0f0')
        progress_section.pack(fill='x', padx=10, pady=(10, 0))
        
        # Main progress bar
        self.id_progress_bar = ttk.Progressbar(progress_section,
                                            mode='determinate',
                                            length=400)
        self.id_progress_bar.pack(fill='x', pady=2)
        
        # Progress info frame
        progress_info_frame = tk.Frame(progress_section, bg='#f0f0f0')
        progress_info_frame.pack(fill='x', pady=2)
        
        # Progress percentage
        self.id_progress_percentage = tk.Label(progress_info_frame,
                                            text="0%",
                                            bg='#f0f0f0',
                                            font=('Arial', 9, 'bold'),
                                            fg='#3498db')
        self.id_progress_percentage.pack(side='left')
        
        # Timer and ETA
        self.id_timer_label = tk.Label(progress_info_frame,
                                    text="Time: 00:00",
                                    bg='#f0f0f0',
                                    font=('Arial', 9),
                                    fg='#7f8c8d')
        self.id_timer_label.pack(side='right')
        
        self.id_eta_label = tk.Label(progress_info_frame,
                                    text="ETA: --:--",
                                    bg='#f0f0f0',
                                    font=('Arial', 9),
                                    fg='#7f8c8d')
        self.id_eta_label.pack(side='right', padx=(0, 20))
        
        # Progress text area - reduced height to 10 lines
        self.id_progress_text = scrolledtext.ScrolledText(right_frame, 
                                                        width=50,
                                                        font=('Courier', 9))
        self.id_progress_text.pack(fill='both', expand=True, padx=10, pady=(10, 5))
        
        # Status label
        self.id_process_status = tk.Label(right_frame,
                                        text="Ready to start ID annotation",
                                        bg='#f0f0f0',
                                        font=('Arial', 10),
                                        fg='#7f8c8d')
        self.id_process_status.pack(pady=2)

    def create_metabolite_id_annotation_interface(self, parent):
        """Create metabolite ID annotation interface (original functionality with re-filtering)"""
        # Create main frame with grid layout for two columns
        content_frame = tk.Frame(parent, bg='#f0f0f0')
        content_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Configure grid weights to distribute space (reduced weights for narrower columns)
        content_frame.grid_columnconfigure(0, weight=1, minsize=180)
        content_frame.grid_columnconfigure(1, weight=1, minsize=180)
        # Use minsize for row instead of weight to ensure content is visible
        content_frame.grid_rowconfigure(0, weight=1, minsize=500)  # Increased minimum height for better visibility
        
        # ========== LEFT COLUMN - Input/Output Settings ==========
        left_frame = tk.LabelFrame(content_frame, 
                                text="🔧 Input/Output Settings", 
                                font=('Arial', 10, 'bold'),
                                bg='#f0f0f0', 
                                fg='#2c3e50')
        left_frame.grid(row=0, column=0, sticky='nsew', padx=(0, 3))
        
        # Create container frame for canvas and scrollbars
        left_scroll_container = tk.Frame(left_frame, bg='#f0f0f0')
        left_scroll_container.pack(fill='both', expand=True, pady=5)
        
        # Create canvas and scrollbars for left frame
        left_canvas = tk.Canvas(left_scroll_container, bg='#f0f0f0', highlightthickness=0)
        left_v_scrollbar = ttk.Scrollbar(left_scroll_container, orient="vertical", command=left_canvas.yview)
        left_h_scrollbar = ttk.Scrollbar(left_scroll_container, orient="horizontal", command=left_canvas.xview)
        left_scrollable_frame = tk.Frame(left_canvas, bg='#f0f0f0')
        
        left_scrollable_frame.bind(
            "<Configure>",
            lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all"))
        )
        
        canvas_window = left_canvas.create_window((0, 0), window=left_scrollable_frame, anchor="nw")
        left_canvas.configure(yscrollcommand=left_v_scrollbar.set, xscrollcommand=left_h_scrollbar.set)
        
        def configure_scroll_region(event):
            left_canvas.configure(scrollregion=left_canvas.bbox("all"))
            # Only expand width if scrollable frame is smaller than canvas
            if left_scrollable_frame.winfo_reqwidth() < event.width:
                left_canvas.itemconfig(canvas_window, width=event.width)
        
        left_canvas.bind('<Configure>', configure_scroll_region)
        
        # Grid layout for canvas and scrollbars
        left_canvas.grid(row=0, column=0, sticky='nsew')
        left_v_scrollbar.grid(row=0, column=1, sticky='ns')
        left_h_scrollbar.grid(row=1, column=0, sticky='ew')
        left_scroll_container.grid_rowconfigure(0, weight=1)
        left_scroll_container.grid_columnconfigure(0, weight=1)
        
        def _on_left_mousewheel(event):
            left_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        def _on_left_shift_mousewheel(event):
            left_canvas.xview_scroll(int(-1*(event.delta/120)), "units")
        
        def _bind_mousewheel_recursive(widget):
            """Recursively bind mousewheel to all child widgets"""
            widget.bind("<MouseWheel>", _on_left_mousewheel)
            widget.bind("<Shift-MouseWheel>", _on_left_shift_mousewheel)
            for child in widget.winfo_children():
                _bind_mousewheel_recursive(child)
        
        left_canvas.bind("<MouseWheel>", _on_left_mousewheel)
        left_canvas.bind("<Shift-MouseWheel>", _on_left_shift_mousewheel)
        # Schedule recursive binding after all widgets are created
        left_scrollable_frame.after(100, lambda: _bind_mousewheel_recursive(left_scrollable_frame))
        
        tk.Label(
            left_scrollable_frame,
            text="Metabolomics ID Annotation",
            bg='#f0f0f0',
            fg='#2c3e50',
            font=('Arial', 16, 'bold')
        ).pack(anchor='w', padx=8, pady=(10, 0))

        # tk.Label(
        #     left_scrollable_frame,
        #     text="This subtab uses cleaned metabolite data from the Data Cleaning tab.",
        #     bg='#f0f0f0',
        #     fg='#7f8c8d',
        #     font=('Arial', 9, 'italic')
        # ).pack(anchor='w', padx=8, pady=(0, 10))

        # File selection frame
        file_frame = tk.LabelFrame(left_scrollable_frame, 
                                text="File Selection", 
                                bg='#f0f0f0')
        file_frame.pack(fill='x', padx=8, pady=8)
        
        # Auto-loaded file display (read-only)
        tk.Label(file_frame, 
                text="Auto-loaded from Data Cleaning:", 
                bg='#f0f0f0',
                font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=(5, 2))
        
        self.id_auto_file_display = tk.Entry(file_frame, 
                                            font=('Arial', 8),
                                            state='readonly',
                                            bg='#e8e8e8')
        self.id_auto_file_display.pack(fill='x', padx=5, pady=2)
        
        # Manual file selection
        tk.Label(file_frame, 
                text="Or select custom Excel file:", 
                bg='#f0f0f0',
                font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=(8, 2))
        
        input_frame = tk.Frame(file_frame, bg='#f0f0f0')
        input_frame.pack(fill='x', padx=5, pady=2)
        
        # Reuse existing self.id_input_file from previous tab creation (don't recreate)
        if not hasattr(self, 'id_input_file'):
            self.id_input_file = tk.StringVar()
        tk.Entry(input_frame, 
                textvariable=self.id_input_file, 
                font=('Arial', 8)).pack(side='left', fill='x', expand=True, padx=(0, 3))
        tk.Button(input_frame, 
                text="Browse", 
                command=lambda: self.browse_file(self.id_input_file, "Excel file for ID annotation"),
                bg='#3498db', 
                fg='white', 
                font=('Arial', 8, 'bold')).pack(side='right')
        
        # Output settings
        output_frame = tk.LabelFrame(left_scrollable_frame, 
                                    text="Output Settings", 
                                    bg='#f0f0f0')
        output_frame.pack(fill='x', padx=8, pady=8)
        
        tk.Label(output_frame, 
                text="Output file:", 
                bg='#f0f0f0',
                font=('Arial', 9)).pack(anchor='w', padx=5)
        
        output_file_frame = tk.Frame(output_frame, bg='#f0f0f0')
        output_file_frame.pack(fill='x', padx=5, pady=5)
        
        # Reuse existing self.id_output_file from previous tab creation (don't recreate)
        if not hasattr(self, 'id_output_file') or not self.id_output_file.get():
            if not hasattr(self, 'id_output_file'):
                self.id_output_file = tk.StringVar()
            self.id_output_file.set(os.path.join(os.getcwd(), "metabolite_ids_annotated.xlsx"))
        tk.Entry(output_file_frame, 
                textvariable=self.id_output_file, 
                font=('Arial', 8)).pack(side='left', fill='x', expand=True, padx=(0, 3))
        tk.Button(output_file_frame, 
                text="Browse", 
                command=lambda: self.browse_file(self.id_output_file, "Excel file for ID annotation results"),
                bg='#3498db', 
                fg='white', 
                font=('Arial', 8, 'bold')).pack(side='right')
        
        # Annotation Mode Selection
        mode_frame = tk.LabelFrame(left_scrollable_frame,
                                  text="Annotation Mode",
                                  bg='#f0f0f0')
        mode_frame.pack(fill='x', padx=8, pady=8)
        
        self.id_annotation_mode = tk.StringVar(value='standard')
        self.id_annotation_mode.trace('w', self._on_annotation_mode_change)
        
        tk.Radiobutton(mode_frame,
                      text="Standard Mode (Pos/Neg sheets required)",
                      variable=self.id_annotation_mode,
                      value='standard',
                      bg='#f0f0f0',
                      font=('Arial', 9)).pack(anchor='w', padx=5, pady=2)
        
        tk.Radiobutton(mode_frame,
                      text="Custom ID Search (any sheet, ID lookup only)",
                      variable=self.id_annotation_mode,
                      value='custom',
                      bg='#f0f0f0',
                      font=('Arial', 9)).pack(anchor='w', padx=5, pady=2)
        
        tk.Label(mode_frame,
                text="ℹ️ Custom mode: No Pos/Neg required, returns IDs only",
                bg='#f0f0f0',
                font=('Arial', 8, 'italic'),
                fg='#7f8c8d').pack(anchor='w', padx=5, pady=(0, 5))
        
        # Verify Columns button (only shown in Custom ID Search mode) - nested inside mode_frame
        self.id_verify_frame = tk.Frame(mode_frame, bg='#f0f0f0')
        
        self.id_verify_columns_button = tk.Button(
            self.id_verify_frame,
            text="🔍 Verify Columns",
            command=self.verify_id_columns,
            bg='#9b59b6',
            fg='white',
            font=('Arial', 9, 'bold'),
            relief='raised',
            borderwidth=2
        )
        
        self.id_verify_warning_label = tk.Label(self.id_verify_frame,
                text="⚠️ Click to verify all required columns before annotation",
                bg='#f0f0f0',
                font=('Arial', 8, 'italic'),
                fg='#7f8c8d')
        
        # Don't pack by default - will be shown when custom mode is selected
        
        # Processing settings
        settings_frame = tk.LabelFrame(left_scrollable_frame, 
            text="Processing Settings", 
            bg='#f0f0f0')
        settings_frame.pack(fill='x', padx=8, pady=8)
        
        # System info
        max_workers = os.cpu_count() or 4
        system_info_label = tk.Label(
            settings_frame,
            text=f"System has {max_workers} CPU cores available",
            bg='#f0f0f0',
            font=('Arial', 8, 'italic'),
            fg='#7f8c8d'
        )
        system_info_label.pack(anchor='w', padx=5, pady=2)
        
        # Workers setting
        workers_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        workers_frame.pack(fill='x', padx=5, pady=5)
        
        tk.Label(workers_frame, 
                text="Parallel workers:", 
                bg='#f0f0f0',
                font=('Arial', 9)).pack(side='left')
        
        # Use optimal worker count (up to 6 workers for best performance)
        self.id_workers = tk.StringVar(value=str(min(6, max_workers)))
        workers_spinbox = tk.Spinbox(workers_frame, 
                                    from_=1, to=max_workers,  # Allow full CPU count
                                    textvariable=self.id_workers, 
                                    width=5)
        workers_spinbox.pack(side='left', padx=(5, 0))
        
        tk.Button(workers_frame, 
                text="Auto", 
                command=lambda: self.id_workers.set(str(max_workers)),  # Use full CPU count
                bg='#9b59b6', 
                fg='white', 
                font=('Arial', 8, 'bold')).pack(side='left', padx=(8, 0))
        
        # ========== MIDDLE COLUMN - Re-filter & Filter Settings ==========
        middle_frame = tk.LabelFrame(content_frame, 
                                    text="🔄 Re-filter & Settings", 
                                    font=('Arial', 10, 'bold'),
                                    bg='#f0f0f0', 
                                    fg='#2c3e50')
        middle_frame.grid(row=0, column=1, sticky='nsew', padx=(3, 0))
        
        # Create container frame for canvas and scrollbars
        middle_scroll_container = tk.Frame(middle_frame, bg='#f0f0f0')
        middle_scroll_container.pack(fill='both', expand=True, pady=5)
        
        # Create canvas and scrollbars for middle frame
        middle_canvas = tk.Canvas(middle_scroll_container, bg='#f0f0f0', highlightthickness=0)
        middle_v_scrollbar = ttk.Scrollbar(middle_scroll_container, orient="vertical", command=middle_canvas.yview)
        middle_h_scrollbar = ttk.Scrollbar(middle_scroll_container, orient="horizontal", command=middle_canvas.xview)
        middle_scrollable_frame = tk.Frame(middle_canvas, bg='#f0f0f0')
        
        middle_scrollable_frame.bind(
            "<Configure>",
            lambda e: middle_canvas.configure(scrollregion=middle_canvas.bbox("all"))
        )
        
        middle_canvas_window = middle_canvas.create_window((0, 0), window=middle_scrollable_frame, anchor="nw")
        middle_canvas.configure(yscrollcommand=middle_v_scrollbar.set, xscrollcommand=middle_h_scrollbar.set)
        
        def configure_middle_scroll_region(event):
            middle_canvas.configure(scrollregion=middle_canvas.bbox("all"))
            # Only expand width if scrollable frame is smaller than canvas
            if middle_scrollable_frame.winfo_reqwidth() < event.width:
                middle_canvas.itemconfig(middle_canvas_window, width=event.width)
        
        middle_canvas.bind('<Configure>', configure_middle_scroll_region)
        
        # Grid layout for canvas and scrollbars
        middle_canvas.grid(row=0, column=0, sticky='nsew')
        middle_v_scrollbar.grid(row=0, column=1, sticky='ns')
        middle_h_scrollbar.grid(row=1, column=0, sticky='ew')
        middle_scroll_container.grid_rowconfigure(0, weight=1)
        middle_scroll_container.grid_columnconfigure(0, weight=1)
        
        def _on_middle_mousewheel(event):
            middle_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        def _on_middle_shift_mousewheel(event):
            middle_canvas.xview_scroll(int(-1*(event.delta/120)), "units")
        
        def _bind_middle_mousewheel_recursive(widget):
            """Recursively bind mousewheel to all child widgets"""
            widget.bind("<MouseWheel>", _on_middle_mousewheel)
            widget.bind("<Shift-MouseWheel>", _on_middle_shift_mousewheel)
            for child in widget.winfo_children():
                _bind_middle_mousewheel_recursive(child)
        
        middle_canvas.bind("<MouseWheel>", _on_middle_mousewheel)
        middle_canvas.bind("<Shift-MouseWheel>", _on_middle_shift_mousewheel)
        # Schedule recursive binding after all widgets are created
        middle_scrollable_frame.after(100, lambda: _bind_middle_mousewheel_recursive(middle_scrollable_frame))
        
        # ========== RE-FILTER EXISTING RESULTS SECTION ==========
        refilter_frame = tk.LabelFrame(middle_scrollable_frame,
                                    text="🔄 Re-filter Existing Annotation Results",
                                    bg='#e8f4f8',
                                    font=('Arial', 9, 'bold'),
                                    fg='#2980b9')
        refilter_frame.pack(fill='x', padx=8, pady=(5, 12))
        
        tk.Label(refilter_frame,
                text="Upload a previously annotated Excel file to apply new filters without re-running annotation.",
                wraplength=350,
                justify='left',
                bg='#e8f4f8',
                fg='#34495e',
                font=('Arial', 8)).pack(anchor='w', padx=6, pady=(5, 5))
        
        # File selection for re-filtering
        refilter_file_frame = tk.Frame(refilter_frame, bg='#e8f4f8')
        refilter_file_frame.pack(fill='x', padx=6, pady=5)
        
        tk.Label(refilter_file_frame,
                text="Annotated Excel file:",
                bg='#e8f4f8',
                font=('Arial', 8, 'bold')).pack(anchor='w')
        
        refilter_input_frame = tk.Frame(refilter_file_frame, bg='#e8f4f8')
        refilter_input_frame.pack(fill='x', pady=2)
        
        self.refilter_input_file = tk.StringVar()
        tk.Entry(refilter_input_frame,
                textvariable=self.refilter_input_file,
                font=('Arial', 8)).pack(side='left', fill='x', expand=True, padx=(0, 3))
        tk.Button(refilter_input_frame,
                text="Browse",
                command=lambda: self.browse_file(self.refilter_input_file, "Select annotated Excel file to re-filter"),
                bg='#3498db',
                fg='white',
                font=('Arial', 7, 'bold')).pack(side='right')
        
        # Output for re-filtered file
        tk.Label(refilter_file_frame,
                text="Save re-filtered results as:",
                bg='#e8f4f8',
                font=('Arial', 8, 'bold')).pack(anchor='w', pady=(6, 0))
        
        refilter_output_frame = tk.Frame(refilter_file_frame, bg='#e8f4f8')
        refilter_output_frame.pack(fill='x', pady=2)
        
        self.refilter_output_file = tk.StringVar(value=os.path.join(os.getcwd(), "refiltered_ids_annotated.xlsx"))
        tk.Entry(refilter_output_frame,
                textvariable=self.refilter_output_file,
                font=('Arial', 8)).pack(side='left', fill='x', expand=True, padx=(0, 3))
        tk.Button(refilter_output_frame,
                text="Browse",
                command=lambda: self.browse_file(self.refilter_output_file, "Save re-filtered results"),
                bg='#3498db',
                fg='white',
                font=('Arial', 7, 'bold')).pack(side='right')
        
        # Re-filter button
        self.refilter_button = tk.Button(refilter_frame,
                                        text="🔄 Apply Filters to Existing Results",
                                        command=self.refilter_existing_results,
                                        bg='#16a085',
                                        fg='white',
                                        font=('Arial', 9, 'bold'),
                                        height=2,
                                        relief='raised',
                                        borderwidth=2)
        self.refilter_button.pack(fill='x', padx=6, pady=(5, 6))
        
        tk.Label(refilter_frame,
                text="Note: Use the filter settings below to configure which filters to apply.",
                wraplength=350,
                justify='left',
                bg='#e8f4f8',
                fg='#7f8c8d',
                font=('Arial', 7, 'italic')).pack(anchor='w', padx=6, pady=(0, 5))
        
        # ========== FILTER SETTINGS ==========
        posneg_frame = tk.LabelFrame(middle_scrollable_frame, 
                                    text="Filter Settings (Applied to Both New Annotation & Re-filtering)", 
                                    bg='#f0f0f0', 
                                    font=('Arial', 9, 'bold'))
        posneg_frame.pack(fill='x', padx=5, pady=(8, 5))

        # ===== SKIP ID FILTERING CHECKBOX =====
        skip_filter_frame = tk.Frame(posneg_frame, bg='#f0f0f0')
        skip_filter_frame.pack(fill='x', padx=8, pady=(5, 8), side='top')
        
        # Load preferences first
        self._loaded_prefs = self.load_id_annotation_prefs()
        
        # Default to unchecked - only internally set True for custom mode
        self.skip_id_filtering_var = tk.BooleanVar(value=False)
        skip_cb = tk.Checkbutton(
            skip_filter_frame,
            text="☑ Skip ID Filtering",
            variable=self.skip_id_filtering_var,
            bg='#f0f0f0',
            font=('Arial', 8, 'bold'),
            fg='#2c3e50',
            command=self._toggle_id_filter_controls
        )
        skip_cb.pack(anchor='w', side='left', padx=(0, 10))
        
        # skip_info = tk.Label(
        #     skip_filter_frame,
        #     text="When checked: Save all metabolites (with or without IDs). ID filtering section below will be disabled.",
        #     wraplength=380,
        #     justify='left',
        #     bg='#f0f0f0',
        #     fg='#7f8c8d',
        #     font=('Arial', 7, 'italic')
        # )
        # skip_info.pack(anchor='w', side='left')

        # ID column checkbox selection panel
        self.available_id_columns = [
            'HMDB_ID', 'PubChem_CID', 'KEGG_ID', 'LipidMaps_ID',
            'ChEBI_ID', 'CAS', 'InChI', 'InChIKey'
        ]

        # Get preference selections
        pref_selected = set(self._loaded_prefs.get('id_selected_columns', self.available_id_columns))

        idcols_outer = tk.LabelFrame(posneg_frame, text="ID Columns to Use For Filtering", bg='#f0f0f0')
        idcols_outer.pack(fill='x', padx=8, pady=(0, 5))

        # Grid frame for 4-across layout
        id_grid = tk.Frame(idcols_outer, bg='#f0f0f0')
        id_grid.pack(fill='x', padx=5, pady=5)

        # Store checkbox variables (4 per row)
        self.id_filter_vars = {}
        cols_per_row = 4
        for idx, col in enumerate(self.available_id_columns):
            var = tk.BooleanVar(value=(col in pref_selected))
            self.id_filter_vars[col] = var
            r = idx // cols_per_row
            c = idx % cols_per_row
            cb = tk.Checkbutton(id_grid, text=col, variable=var, bg='#f0f0f0', 
                            font=('Arial', 8),
                            command=lambda: self.update_id_selected_count())
            cb.grid(row=r, column=c, sticky='w', padx=3, pady=2)
            # Add tooltips
            tooltip_map = {
                'HMDB_ID': 'Human Metabolome Database ID',
                'PubChem_CID': 'PubChem Compound Identifier',
                'KEGG_ID': 'KEGG Compound/Lipid ID',
                'LipidMaps_ID': 'LipidMaps Structural Database ID',
                'ChEBI_ID': 'ChEBI (Chemical Entities of Biological Interest) ID',
                'CAS': 'CAS Registry Number',
                'InChI': 'IUPAC International Chemical Identifier',
                'InChIKey': 'Hashed (27-char) InChIKey'
            }
            self._create_tooltip(cb, tooltip_map.get(col, col))
        
        # Add stretch columns
        for c in range(cols_per_row):
            id_grid.grid_columnconfigure(c, weight=1)

        # Control buttons row
        btn_row = tk.Frame(idcols_outer, bg='#f0f0f0')
        btn_row.pack(fill='x', padx=5, pady=(0, 5))
        tk.Button(btn_row, text="Select All", width=10,
            command=lambda: [v.set(True) or self.update_id_selected_count() for v in self.id_filter_vars.values()],
                    bg='#27ae60', fg='white', font=('Arial', 7, 'bold')).pack(side='left', padx=(0, 3))
        tk.Button(btn_row, text="Deselect All", width=10,
            command=lambda: [v.set(False) or self.update_id_selected_count() for v in self.id_filter_vars.values()],
                    bg='#e74c3c', fg='white', font=('Arial', 7, 'bold')).pack(side='left')

        # Live selection count label
        self.id_selected_count_label = tk.Label(idcols_outer, text='', bg='#f0f0f0', 
                                                font=('Arial', 7, 'italic'), fg='#7f8c8d')
        self.id_selected_count_label.pack(anchor='w', padx=8, pady=(0, 3))
        self.update_id_selected_count()
        
        # Apply initial toggle state based on loaded preferences
        self._toggle_id_filter_controls()

        # Info label
        tk.Label(idcols_outer,
            text="Rule: Selecting at least one ID column will automatically filter Pos/Neg sheets to rows with any selected ID present.",
            wraplength=380,
            justify='left',
            bg='#f0f0f0',
            fg='#566573',
            font=('Arial', 7, 'italic')).pack(anchor='w', padx=8, pady=(0, 4))

        # Backward compatibility string var
        self.selected_id_columns_str = tk.StringVar(value=','.join(self.available_id_columns))

        def _sync_id_selection_string():
            selected = [c for c, v in self.id_filter_vars.items() if v.get()]
            self.selected_id_columns_str.set(','.join(selected))
        self.sync_id_selection_string = _sync_id_selection_string

        # MS2 filter
        ms2_filter_frame = tk.Frame(posneg_frame, bg='#f0f0f0')
        ms2_filter_frame.pack(fill='x', padx=8, pady=(2, 2))
        tk.Label(ms2_filter_frame, text="MS2 Filter:", bg='#f0f0f0', font=('Arial', 8)).pack(side='left')
        self.id_ms2_filter_var = tk.StringVar(value='none')
        ms2_options_display = [
            'No MS2 filter',
            "Keep 'DDA for preferred ion' only",
            "Keep 'DDA for preferred ion' and 'DDA for other ion'"
        ]
        self._id_ms2_display_to_code = {
            'No MS2 filter': 'none',
            "Keep 'DDA for preferred ion' only": 'preferred_only',
            "Keep 'DDA for preferred ion' and 'DDA for other ion'": 'preferred_or_other'
        }
        ms2_combo = ttk.Combobox(ms2_filter_frame, values=ms2_options_display, state='readonly', width=38)
        # Restore preference
        pref_ms2 = self._loaded_prefs.get('ms2_filter', 'none')
        reverse_map = {v: k for k, v in self._id_ms2_display_to_code.items()}
        if pref_ms2 in reverse_map:
            ms2_combo.set(reverse_map[pref_ms2])
            self.id_ms2_filter_var.set(pref_ms2)
        else:
            ms2_combo.set(ms2_options_display[0])
        ms2_combo.pack(side='left', padx=5, fill='x', expand=True)
        self._create_tooltip(ms2_combo, "Filter polarity tables by MS2 acquisition type (DDA).")

        def _sync_ms2_choice(event=None):
            self.id_ms2_filter_var.set(self._id_ms2_display_to_code.get(ms2_combo.get(), 'none'))
        ms2_combo.bind('<<ComboboxSelected>>', _sync_ms2_choice)

        # Confidence filter
        confidence_filter_frame = tk.Frame(posneg_frame, bg='#f0f0f0')
        confidence_filter_frame.pack(fill='x', padx=8, pady=(2, 2))
        tk.Label(confidence_filter_frame, text="Confidence Filter:", bg='#f0f0f0', font=('Arial', 8)).pack(side='left')
        self.id_confidence_filter_var = tk.StringVar(value='exclude_low')
        confidence_options_display = [
            'Exclude low confidence (default)',
            'No confidence filter'
        ]
        self._id_conf_display_to_code = {
            'Exclude low confidence (default)': 'exclude_low',
            'No confidence filter': 'none'
        }
        conf_combo = ttk.Combobox(confidence_filter_frame, values=confidence_options_display, state='readonly', width=38)
        pref_conf = self._loaded_prefs.get('confidence_filter', 'exclude_low')
        conf_reverse_map = {v: k for k, v in self._id_conf_display_to_code.items()}
        if pref_conf in conf_reverse_map:
            conf_combo.set(conf_reverse_map[pref_conf])
            self.id_confidence_filter_var.set(pref_conf)
        else:
            conf_combo.set(confidence_options_display[0])
            self.id_confidence_filter_var.set('exclude_low')
        conf_combo.pack(side='left', padx=5, fill='x', expand=True)
        self._create_tooltip(conf_combo, "Filter rows by Metabograph_Confidence. Default removes low-confidence rows.")

        def _sync_conf_choice(event=None):
            self.id_confidence_filter_var.set(self._id_conf_display_to_code.get(conf_combo.get(), 'exclude_low'))
        conf_combo.bind('<<ComboboxSelected>>', _sync_conf_choice)

        # ID-stage deduplication RT window (minutes)
        rtwin_frame = tk.Frame(posneg_frame, bg='#f0f0f0')
        rtwin_frame.pack(fill='x', padx=8, pady=(2, 2))
        tk.Label(rtwin_frame, text="RT window for ID dedup (min):", bg='#f0f0f0', font=('Arial', 8)).pack(side='left')
        default_rtwin = str(self._loaded_prefs.get('dedup_rt_window_minutes', 2.0))
        self.id_dedup_rt_window = tk.StringVar(value=default_rtwin)
        tk.Entry(rtwin_frame, textvariable=self.id_dedup_rt_window, width=8).pack(side='left', padx=6)
        self._create_tooltip(rtwin_frame, "Applies only during ID annotation: when collapsing Formula+ID duplicates within each polarity, rows within ±this RT window of the group's mean RT will be combined (numeric summed).")

        # Toggle: require Endogenous_Source == Yes
        self.require_endogenous_yes = tk.BooleanVar(value=bool(self._loaded_prefs.get('require_endogenous_yes', False)))
        tk.Checkbutton(posneg_frame,
            text="Also require Endogenous_Source == Yes",
            variable=self.require_endogenous_yes,
            font=('Arial', 8),
            bg='#f0f0f0').pack(anchor='w', padx=8, pady=(2, 5))

    def create_lipid_id_annotation_interface(self, parent):
        """Create lipid ID annotation interface (new Lipids subtab)."""
        content_frame = tk.Frame(parent, bg='#f0f0f0')
        content_frame.pack(fill='both', expand=True, padx=5, pady=5)

        # Left column: input/output settings for lipids
        left_frame = tk.LabelFrame(content_frame,
                                   text="🔧 Lipid Input/Output Settings",
                                   font=('Arial', 10, 'bold'),
                                   bg='#f0f0f0',
                                   fg='#2c3e50')
        left_frame.grid(row=0, column=0, sticky='nsew', padx=(0, 3))

        content_frame.grid_columnconfigure(0, weight=1)
        content_frame.grid_rowconfigure(0, weight=1)

        # Create container frame for canvas and scrollbars
        left_scroll_container = tk.Frame(left_frame, bg='#f0f0f0')
        left_scroll_container.pack(fill='both', expand=True, padx=5, pady=5)

        # Scrollable wrapper for left column with both scrollbars
        left_canvas = tk.Canvas(left_scroll_container, bg='#f0f0f0', highlightthickness=0)
        left_v_scrollbar = ttk.Scrollbar(left_scroll_container, orient='vertical', command=left_canvas.yview)
        left_h_scrollbar = ttk.Scrollbar(left_scroll_container, orient='horizontal', command=left_canvas.xview)
        left_scrollable = tk.Frame(left_canvas, bg='#f0f0f0')

        left_scrollable.bind(
            '<Configure>',
            lambda e: left_canvas.configure(scrollregion=left_canvas.bbox('all'))
        )

        canvas_window = left_canvas.create_window((0, 0), window=left_scrollable, anchor='nw')
        left_canvas.configure(yscrollcommand=left_v_scrollbar.set, xscrollcommand=left_h_scrollbar.set)

        def _configure_left_canvas(event):
            left_canvas.configure(scrollregion=left_canvas.bbox('all'))
            try:
                # Only expand width if scrollable frame is smaller than canvas
                if left_scrollable.winfo_reqwidth() < event.width:
                    left_canvas.itemconfig(canvas_window, width=event.width)
            except Exception:
                pass

        left_canvas.bind('<Configure>', _configure_left_canvas)
        
        # Grid layout for canvas and scrollbars
        left_canvas.grid(row=0, column=0, sticky='nsew')
        left_v_scrollbar.grid(row=0, column=1, sticky='ns')
        left_h_scrollbar.grid(row=1, column=0, sticky='ew')
        left_scroll_container.grid_rowconfigure(0, weight=1)
        left_scroll_container.grid_columnconfigure(0, weight=1)

        def _on_left_wheel(event):
            left_canvas.yview_scroll(int(-1 * (event.delta / 120)), 'units')

        def _on_left_shift_wheel(event):
            left_canvas.xview_scroll(int(-1 * (event.delta / 120)), 'units')

        left_canvas.bind('<MouseWheel>', _on_left_wheel)
        left_canvas.bind('<Shift-MouseWheel>', _on_left_shift_wheel)
        left_scrollable.bind('<MouseWheel>', _on_left_wheel)
        left_scrollable.bind('<Shift-MouseWheel>', _on_left_shift_wheel)

        # Title and description
        tk.Label(
            left_scrollable,
            text="Lipid ID Annotation",
            bg='#f0f0f0',
            fg='#2c3e50',
            font=('Arial', 16, 'bold')
        ).pack(anchor='w', padx=8, pady=(10, 0))

        # tk.Label(
        #     left_scrollable,
        #     text="Uses cleaned lipid data from the Data Cleaning tab to perform ID annotation only (no pathway annotation).",
        #     bg='#f0f0f0',
        #     fg='#7f8c8d',
        #     wraplength=380,
        #     justify='left',
        #     font=('Arial', 9, 'italic')
        # ).pack(anchor='w', padx=8, pady=(0, 10))

        # File selection frame
        file_frame = tk.LabelFrame(left_scrollable,
                                   text="File Selection",
                                   bg='#f0f0f0')
        file_frame.pack(fill='x', padx=8, pady=8)

        # Auto-loaded file display (read-only)
        tk.Label(file_frame,
                 text="Auto-loaded from Data Cleaning:",
                 bg='#f0f0f0',
                 font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=(5, 2))

        self.lipid_auto_file_display = tk.Entry(file_frame,
                                                font=('Arial', 8),
                                                state='readonly',
                                                bg='#e8e8e8')
        self.lipid_auto_file_display.pack(fill='x', padx=5, pady=2)

        # Manual file selection
        tk.Label(file_frame,
                 text="Or select custom Excel file:",
                 bg='#f0f0f0',
                 font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=(8, 2))

        input_frame = tk.Frame(file_frame, bg='#f0f0f0')
        input_frame.pack(fill='x', padx=5, pady=2)

        if not hasattr(self, 'lipid_id_annotation_file'):
            self.lipid_id_annotation_file = tk.StringVar()

        tk.Entry(input_frame,
                 textvariable=self.lipid_id_annotation_file,
                 font=('Arial', 8)).pack(side='left', fill='x', expand=True, padx=(0, 3))

        tk.Button(input_frame,
                  text="Browse",
                  command=lambda: self.browse_file(self.lipid_id_annotation_file, "cleaned lipids Excel file"),
                  bg='#3498db',
                  fg='white',
                  font=('Arial', 8, 'bold')).pack(side='right')

        # Output settings
        output_frame = tk.LabelFrame(left_scrollable,
                                     text="Output Settings",
                                     bg='#f0f0f0')
        output_frame.pack(fill='x', padx=8, pady=8)

        tk.Label(output_frame,
                 text="Output file:",
                 bg='#f0f0f0',
                 font=('Arial', 9)).pack(anchor='w', padx=5)

        output_file_frame = tk.Frame(output_frame, bg='#f0f0f0')
        output_file_frame.pack(fill='x', padx=5, pady=5)

        if not hasattr(self, 'lipid_output_file'):
            self.lipid_output_file = tk.StringVar()
        if not self.lipid_output_file.get():
            self.lipid_output_file.set(os.path.join(os.getcwd(), "lipid_ids_annotated.xlsx"))

        tk.Entry(output_file_frame,
                 textvariable=self.lipid_output_file,
                 font=('Arial', 8)).pack(side='left', fill='x', expand=True, padx=(0, 3))

        tk.Button(output_file_frame,
                  text="Browse",
                  command=lambda: self.browse_file(self.lipid_output_file, "Excel file for lipid ID annotation results"),
                  bg='#3498db',
                  fg='white',
                  font=('Arial', 8, 'bold')).pack(side='right')

        # Workers for lipid mode
        settings_frame = tk.LabelFrame(left_scrollable,
                                       text="Processing Settings",
                                       bg='#f0f0f0')
        settings_frame.pack(fill='x', padx=8, pady=8)

        max_workers = os.cpu_count() or 4
        tk.Label(settings_frame,
                 text=f"System has {max_workers} CPU cores available",
                 bg='#f0f0f0',
                 font=('Arial', 8, 'italic'),
                 fg='#7f8c8d').pack(anchor='w', padx=5, pady=(5, 2))

        if not hasattr(self, '_loaded_prefs'):
            self._loaded_prefs = self.load_id_annotation_prefs()

        if not hasattr(self, 'lipid_keep_all_rows_var'):
            self.lipid_keep_all_rows_var = tk.BooleanVar(
                value=bool(self._loaded_prefs.get('lipid_keep_all_rows', False))
            )

        keep_rows_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        keep_rows_frame.pack(fill='x', padx=5, pady=(4, 2))

        tk.Checkbutton(
            keep_rows_frame,
            text="Keep rows with no ID (do not remove unmatched rows)",
            variable=self.lipid_keep_all_rows_var,
            bg='#f0f0f0',
            font=('Arial', 8, 'bold'),
            fg='#2c3e50'
        ).pack(anchor='w')

        tk.Label(
            keep_rows_frame,
            text="Default: unchecked (rows without IDs are removed)",
            bg='#f0f0f0',
            fg='#7f8c8d',
            font=('Arial', 7, 'italic')
        ).pack(anchor='w', padx=(20, 0))

        workers_frame = tk.Frame(settings_frame, bg='#f0f0f0')
        workers_frame.pack(fill='x', padx=5, pady=5)

        tk.Label(workers_frame,
                 text="Parallel workers:",
                 bg='#f0f0f0',
                 font=('Arial', 9)).pack(side='left')

        if not hasattr(self, 'lipid_workers'):
            self.lipid_workers = tk.StringVar(value=str(min(6, max_workers)))

        tk.Spinbox(workers_frame,
                   from_=1,
                   to=max_workers,
                   textvariable=self.lipid_workers,
                   width=5).pack(side='left', padx=(5, 0))

        tk.Button(workers_frame,
                  text="Auto",
                  command=lambda: self.lipid_workers.set(str(min(6, max_workers))),
                  bg='#9b59b6',
                  fg='white',
                  font=('Arial', 8, 'bold')).pack(side='left', padx=(8, 0))

    
    def _on_annotation_mode_change(self, *args):
        """Handle annotation mode change to show/hide Verify Columns button"""
        mode = self.id_annotation_mode.get()
        
        if mode == 'custom':
            # Show Verify Columns button and warning for custom mode (nested within mode_frame)
            self.id_verify_frame.pack(fill='x', padx=5, pady=(5, 10))
            self.id_verify_columns_button.pack(fill='x', pady=2)
            self.id_verify_warning_label.pack(anchor='w', pady=(2, 0))
            # Note: skip_id_filtering is enforced programmatically in start_id_annotation()
            # for custom mode - we don't change the checkbox here to avoid affecting standard mode
        else:
            # Hide for standard mode - no column verification needed
            self.id_verify_frame.pack_forget()
    

    def on_start_annotation_clicked(self):
        """Dispatch start based on current subtab (metabolite vs lipid)."""
        current_mode = 'metabolite'
        try:
            if hasattr(self, 'id_mode'):
                current_mode = self.id_mode.get() or 'metabolite'
        except Exception:
            current_mode = 'metabolite'

        if current_mode == 'lipid':
            self.start_lipid_id_annotation()
        else:
            self.start_id_annotation()


    def open_id_output_folder(self):
        """Open the folder containing ID annotation output files"""
        import subprocess
        import platform
        
        output_file = self.id_output_file.get()
        
        if not output_file:
            messagebox.showinfo("No Output File", "No output file path set. Please run ID annotation first or set an output path.")
            return
        
        folder_path = os.path.dirname(output_file)
        if not os.path.exists(folder_path):
            messagebox.showerror("Folder Not Found", f"Output folder does not exist:\n{folder_path}")
            return
        
        try:
            if platform.system() == 'Windows':
                os.startfile(folder_path)
            elif platform.system() == 'Darwin':  # macOS
                subprocess.Popen(['open', folder_path])
            else:  # Linux
                subprocess.Popen(['xdg-open', folder_path])
        except Exception as e:
            messagebox.showerror("Error", f"Could not open folder:\n{str(e)}")
    
    def verify_id_columns(self):
        """Verify and assign columns for Custom ID Search mode"""
        try:
            input_file = self.id_input_file.get() or getattr(self, 'cleaned_excel_path', None)
            
            if not input_file:
                messagebox.showerror(
                    "Error",
                    "No file selected.\n\n"
                    "Please select a file or run data cleaning first."
                )
                return
            
            if not os.path.exists(input_file):
                messagebox.showerror("Error", f"Selected file does not exist:\n{input_file}")
                return
            
            # Custom ID Search mode - load first available sheet
            self.id_progress_text.insert(tk.END, f"\n🔍 Opening column assignment dialog for Custom ID Search mode...\n")
            self.id_progress_text.see(tk.END)

            result = self._prompt_id_annotation_column_assignment(input_file, store_as_custom=True)
            
            if result:
                # Store column assignments and dataframe for use in annotation
                self._custom_id_column_mapping = result['assignments']
                self._custom_id_dataframe = result['dataframe']
                self._custom_id_selected_sheet = result.get('selected_sheet')
                
                # Debug: show what columns are in the dataframe
                print(f"🔍 DEBUG: Stored dataframe columns: {list(result['dataframe'].columns)}")
                print(f"🔍 DEBUG: Column mappings: {result['assignments']}")
                
                col_info = "\n".join([f"  • {k}: {v}" for k, v in result['assignments'].items() if v])
                selected_sheet = result.get('selected_sheet')
                sheet_info = f"\n  • Using sheet: {selected_sheet}" if selected_sheet else ""
                self.id_progress_text.insert(tk.END, f"✅ Custom ID search column assignment complete!{sheet_info}\n{col_info}\n")
                messagebox.showinfo(
                    "✅ Column Assignment Complete",
                    f"Column assignments confirmed:\n\n{col_info}{sheet_info}"
                )
            else:
                self.id_progress_text.insert(tk.END, f"❌ Column assignment cancelled.\n")
                # Clear any previous mappings
                self._custom_id_column_mapping = None
                self._custom_id_dataframe = None
                self._custom_id_selected_sheet = None
            
            self.id_progress_text.see(tk.END)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to verify columns: {str(e)}")
            logger.error(f"Verification error: {e}", exc_info=True)

    def start_id_annotation(self):
        """Start the ID annotation process"""
        # Get input file (metabolite mode)
        if self.id_input_file.get():
            input_file = self.id_input_file.get()
        elif hasattr(self, 'cleaned_excel_path') and self.cleaned_excel_path:
            input_file = self.cleaned_excel_path
        else:
            input_file = None
        
        output_file = self.id_output_file.get()
        workers = self._get_workers_count(self.id_workers)
        mode = 'metabolite'
        
        if not input_file:
            messagebox.showerror("Error", "Please select an input file or run data cleaning first")
            return
        
        if not os.path.exists(input_file):
            messagebox.showerror("Error", "Input file does not exist")
            return
        
        if not output_file:
            messagebox.showerror("Error", "Please specify an output file")
            return
        
        # Initialize progress tracking
        self.id_start_time = time.time()
        self.id_current_step = 0
        self.id_total_steps = 100
        
        # Setup progress bar
        self.id_progress_bar.config(mode='determinate', maximum=100, value=0)
        self.id_progress_percentage.config(text="0%")
        self.id_timer_label.config(text="Time: 00:00")
        self.id_eta_label.config(text="ETA: --:--")
        
        # Update UI
        self.id_start_button.config(state='disabled', text="Processing...")
        self.id_stop_button.config(state='normal')
        self.id_process_status.config(text="Initializing ID annotation...", fg='#f39c12')
        
        # Clear progress
        self.id_progress_text.delete(1.0, tk.END)
        self.update_id_progress_with_percentage(0, "Starting ID annotation...")
        
        # Initialize stop flag
        self.id_annotation_stop_flag = False
        
        # Sync selected ID columns from checkboxes to hidden string var
        if hasattr(self, 'sync_id_selection_string'):
            self.sync_id_selection_string()

        # Check if custom ID search mode is enabled FIRST
        self._custom_id_search_mode = self.id_annotation_mode.get() == 'custom' if hasattr(self, 'id_annotation_mode') else False
        
        # Derive selected ID columns list for passing to annotator
        selected_id_cols = [c for c, v in getattr(self, 'id_filter_vars', {}).items() if v.get()]
        
        # CRITICAL: Custom ID Search mode should NEVER filter rows - always skip filtering
        skip_id_filtering = getattr(self, 'skip_id_filtering_var', tk.BooleanVar(value=False)).get()
        if self._custom_id_search_mode:
            skip_id_filtering = True  # Force skip filtering for custom mode
            
        # Automatic rule: filtering active iff at least one column selected AND skip_id_filtering is False
        effective_filter_mode = (len(selected_id_cols) > 0) and (not skip_id_filtering)

        # Attach for thread consumption
        self._selected_id_cols_runtime = selected_id_cols if not skip_id_filtering else []
        self._effective_id_filter_mode_runtime = effective_filter_mode
        self._skip_id_filtering_runtime = skip_id_filtering
        self._annotation_mode_runtime = mode  # Store mode for run_id_annotation
        self._annotation_workers_runtime = workers  # Store workers count
        
        # Warn if custom mode is enabled but columns haven't been verified
        if self._custom_id_search_mode:
            column_mapping = getattr(self, '_custom_id_column_mapping', None)
            
            if not column_mapping:
                response = messagebox.askyesno(
                    "⚠️ Column Verification Recommended",
                    "You haven't verified column assignments for Custom ID Search mode.\n\n"
                    "Without verification, the column names in your file must exactly match:\n"
                    "  • MetID, Name\n\n"
                    "Click 'Verify Columns' button to map your columns to the required names.\n\n"
                    "Continue anyway?",
                    icon='warning'
                )
                if not response:
                    self.id_start_button.config(state='normal', text="🚀 Start ID Annotation")
                    self.id_stop_button.config(state='disabled')
                    self.id_process_status.config(text="Ready", fg='#27ae60')
                    return

        # Persist preferences immediately
        self.save_id_annotation_prefs()

        # Start ID annotation in separate thread
        self.id_annotation_thread = threading.Thread(
            target=self.run_id_annotation, 
            args=(input_file,)
        )
        self.id_annotation_thread.daemon = True
        self.id_annotation_thread.start()
        
        # Start timer updates
        self.update_id_timer()

    def start_lipid_id_annotation(self):
        """Start the lipid ID annotation process (Lipids subtab)."""
        # Determine input file for lipids
        input_file = None
        # 1) Explicit user selection
        if hasattr(self, 'lipid_id_annotation_file') and self.lipid_id_annotation_file.get():
            input_file = self.lipid_id_annotation_file.get()
        # 2) Cleaned lipid path stored on this tab (from Data Cleaning)
        elif hasattr(self, 'cleaned_lipid_excel_path') and getattr(self, 'cleaned_lipid_excel_path', None):
            input_file = self.cleaned_lipid_excel_path
        # 3) DataManager cleaned lipid path
        elif hasattr(self, 'data_manager') and getattr(self.data_manager, 'cleaned_lipid_excel_path', None):
            input_file = self.data_manager.cleaned_lipid_excel_path

        # Output file for lipid IDs - now a full path from the interface
        output_file = os.path.join(os.getcwd(), "lipid_ids_annotated.xlsx")  # Default
        try:
            if hasattr(self, 'lipid_output_file') and self.lipid_output_file.get():
                output_file = self.lipid_output_file.get()
                # Ensure it has .xlsx extension
                if not output_file.lower().endswith('.xlsx'):
                    output_file += '.xlsx'
        except Exception:
            pass

        # Workers
        workers = self._get_workers_count(getattr(self, 'lipid_workers', None) or self.id_workers)
        mode = 'lipid'

        if not input_file:
            messagebox.showerror("Error", "Please run lipid data cleaning first or select a cleaned lipids Excel file.")
            return
        if not os.path.exists(input_file):
            messagebox.showerror("Error", f"Lipid input file does not exist:\n{input_file}")
            return

        # Initialize progress tracking
        self.id_start_time = time.time()
        self.id_current_step = 0
        self.id_total_steps = 100
        
        # Setup progress bar
        self.id_progress_bar.config(mode='determinate', maximum=100, value=0)
        self.id_progress_percentage.config(text="0%")
        self.id_timer_label.config(text="Time: 00:00")
        self.id_eta_label.config(text="ETA: --:--")
        
        # Update UI
        self.id_start_button.config(state='disabled', text="Processing...")
        self.id_stop_button.config(state='normal')
        self.id_process_status.config(text="Initializing lipid ID annotation...", fg='#f39c12')
        
        # Clear progress
        self.id_progress_text.delete(1.0, tk.END)
        self.update_id_progress_with_percentage(0, "Starting lipid ID annotation...")
        
        # Initialize stop flag
        self.id_annotation_stop_flag = False

        # For lipid mode, row-removal can be controlled by the lipid checkbox:
        # unchecked (default) -> remove rows without IDs
        # checked -> keep all rows
        selected_id_cols = []
        effective_filter_mode = False
        skip_id_filtering = getattr(self, 'lipid_keep_all_rows_var', tk.BooleanVar(value=False)).get()
        self._selected_id_cols_runtime = selected_id_cols
        self._effective_id_filter_mode_runtime = effective_filter_mode
        self._skip_id_filtering_runtime = bool(skip_id_filtering)
        self._annotation_mode_runtime = mode
        self._annotation_workers_runtime = workers
        self._custom_id_search_mode = False

        # Persist lipid output path on the tab
        try:
            if hasattr(self, 'lipid_output_file'):
                self.lipid_output_file.set(output_file)
        except Exception:
            pass

        # Start ID annotation in separate thread (lipid mode)
        self.id_annotation_thread = threading.Thread(
            target=self.run_id_annotation,
            args=(input_file,)
        )
        self.id_annotation_thread.daemon = True
        self.id_annotation_thread.start()
        
        # Start timer updates
        self.update_id_timer()
    
    def stop_id_annotation(self):
        """Stop the ID annotation process"""
        # Set the stop flag
        self.id_annotation_stop_flag = True
        
        # Update UI
        self.id_start_button.config(state='normal', text="🚀 Start ID Annotation")
        self.id_stop_button.config(state='disabled')
        self.id_process_status.config(text="Stopping annotation process...", fg='#e74c3c')
        self.id_start_time = None
        
        # Log the stop request
        self.id_progress_text.insert(tk.END, "\n⛔ Stop requested by user - waiting for current operation to complete...\n")
        self.id_progress_text.see(tk.END)
    
    def refilter_existing_results(self):
        """Re-filter an existing annotated Excel file with new filter settings"""
        input_file = self.refilter_input_file.get()
        output_file = self.refilter_output_file.get()
        
        if not input_file:
            messagebox.showerror("Error", "Please select an annotated Excel file to re-filter")
            return
        
        if not os.path.exists(input_file):
            messagebox.showerror("Error", "Input file does not exist")
            return
        
        if not output_file:
            messagebox.showerror("Error", "Please specify an output file for re-filtered results")
            return
        
    # Get current filter settings
        self.sync_id_selection_string()  # Sync checkbox state to string var
        selected_id_cols = [c for c, v in self.id_filter_vars.items() if v.get()]
        skip_id_filtering = getattr(self, 'skip_id_filtering_var', tk.BooleanVar(value=False)).get()
        ms2_filter_mode = self.id_ms2_filter_var.get()
        confidence_filter_mode = getattr(self, 'id_confidence_filter_var', tk.StringVar(value='exclude_low')).get()
        require_endogenous = self.require_endogenous_yes.get()
        
        # Update status
        self.id_process_status.config(text="Re-filtering existing results...", fg='#f39c12')
        self.id_progress_text.delete('1.0', tk.END)
        self.id_progress_text.insert(tk.END, f"Starting re-filter process...\n")
        self.id_progress_text.insert(tk.END, f"Input: {input_file}\n")
        self.id_progress_text.insert(tk.END, f"Output: {output_file}\n")
        if skip_id_filtering:
            self.id_progress_text.insert(tk.END, f"Skip ID Filtering: YES (all metabolites will be saved)\n")
        else:
            self.id_progress_text.insert(tk.END, f"Selected ID columns: {', '.join(selected_id_cols) if selected_id_cols else 'None (no ID filtering)'}\n")
        self.id_progress_text.insert(tk.END, f"MS2 filter: {ms2_filter_mode}\n")
        self.id_progress_text.insert(tk.END, f"Confidence filter: {confidence_filter_mode}\n")
        self.id_progress_text.insert(tk.END, f"Require Endogenous_Source=Yes: {require_endogenous}\n")
        self.id_progress_text.insert(tk.END, f"\nLoading existing annotation file...\n")
        self.id_progress_text.see(tk.END)
        self.root.update()
        
        try:
            # Probe schema to detect lipid-mode, to inform logs (backend also enforces bypass)
            lipid_mode_hint = False
            try:
                xl_probe = pd.ExcelFile(input_file)
                probe_df = None
                for s in ['Pos_id', 'Neg_id', 'Merged_IDs']:
                    if s in xl_probe.sheet_names:
                        probe_df = pd.read_excel(input_file, sheet_name=s, nrows=1)
                        if probe_df is not None:
                            break
                if probe_df is not None and not probe_df.empty:
                    cols = {str(c).strip() for c in probe_df.columns}
                    if ('LipidID' in cols) or ({'ABBREVIATION', 'MAIN_CLASS'} & cols):
                        lipid_mode_hint = True
            except Exception:
                pass

            if lipid_mode_hint:
                self.id_progress_text.insert(tk.END, "Detected lipid-mode annotated file → will bypass ID/MS2/Endogenous filters.\n")
                # Reflect intended bypass in the local variables shown to user; backend already enforces
                selected_id_cols = []
                ms2_filter_mode = 'none'
                confidence_filter_mode = 'none'
                require_endogenous = False

            # Import the re-filter function from metabolite_ID_annotator
            from main_script.metabolite_ID_annotator import refilter_annotated_excel
            
            # Call the re-filter function
            result = refilter_annotated_excel(
                input_file=input_file,
                output_file=output_file,
                selected_id_columns=selected_id_cols if not skip_id_filtering else [],
                skip_id_filtering=skip_id_filtering,
                ms2_filter_mode=ms2_filter_mode,
                confidence_filter_mode=confidence_filter_mode,
                require_endogenous_yes=require_endogenous,
                progress_callback=lambda msg: self._update_refilter_progress(msg)
            )
            
            if result.get('success'):
                self.id_process_status.config(text="Re-filtering completed successfully!", fg='#27ae60')
                self.id_progress_text.insert(tk.END, f"\n✓ Re-filtering completed successfully!\n")
                self.id_progress_text.insert(tk.END, f"Output saved to: {output_file}\n")
                
                if result.get('summary'):
                    self.id_progress_text.insert(tk.END, f"\nSummary:\n{result['summary']}\n")
                
                self.id_progress_text.see(tk.END)
                
                # Store the output file path for auto-loading to statistics tab
                self.annotated_metabolites_excel_path = output_file
                self.id_annotated_excel_path = output_file
                self.annotated_ids_excel_path = output_file
                
                # Auto-load to Statistics tab (same as after ID annotation)
                # Only auto-load to Statistics tab if not in lipid mode
                # Don't try to auto-load to Statistics tab - it's in a separate tool
                self.id_progress_text.insert(tk.END, f"\n✓ Re-filtering completed successfully!\n")
                self.id_progress_text.insert(tk.END, f"\n📋 Next: Use Statistics tool to generate p-values and fold changes\n")
                self.id_progress_text.see(tk.END)
                messagebox.showinfo("Success", f"Re-filtering completed!\n\nOutput saved to:\n{output_file}\n\nNext: Use Statistics tool to add p-values/FC")
            else:
                error_msg = result.get('error', 'Unknown error')
                self.id_process_status.config(text="Re-filtering failed", fg='#e74c3c')
                self.id_progress_text.insert(tk.END, f"\n✗ Error: {error_msg}\n")
                self.id_progress_text.see(tk.END)
                messagebox.showerror("Error", f"Re-filtering failed:\n{error_msg}")
                
        except Exception as e:
            self.id_process_status.config(text="Re-filtering failed", fg='#e74c3c')
            self.id_progress_text.insert(tk.END, f"\n✗ Error: {str(e)}\n")
            self.id_progress_text.see(tk.END)
            messagebox.showerror("Error", f"Re-filtering failed:\n{str(e)}")
    
    def _update_refilter_progress(self, message):
        """Update progress text during re-filtering"""
        self.id_progress_text.insert(tk.END, f"{message}\n")
        self.id_progress_text.see(tk.END)
        self.root.update()
    
    def get_cleaned_file_path(self):
        """Get the path of cleaned file from Tab 1"""
        if hasattr(self, 'cleaned_data') and self.cleaned_data and self.cleaned_data.get('success'):
            files_created = self.cleaned_data.get('files_created', [])
            for file_path in files_created:
                if 'cleaned_metabolites' in os.path.basename(file_path):
                    return file_path
        return None
    
    def update_auto_loaded_file_display(self, file_path=None):
        """Update the auto-loaded file display"""
        # If file_path is provided (from auto-load), store it
        if file_path:
            self.cleaned_excel_path = file_path
        
        if hasattr(self, 'cleaned_excel_path') and self.cleaned_excel_path and os.path.exists(self.cleaned_excel_path):
            self.id_auto_file_display.config(state='normal')
            self.id_auto_file_display.delete(0, tk.END)
            self.id_auto_file_display.insert(0, self.cleaned_excel_path)
            self.id_auto_file_display.config(state='readonly')
            
            # Add log message
            if hasattr(self, 'id_progress_text'):
                self.id_progress_text.config(state='normal')
                self.id_progress_text.insert(tk.END, f"\n✅ Auto-loaded cleaned file: {self.cleaned_excel_path}\n")
                self.id_progress_text.see(tk.END)
                self.id_progress_text.config(state='disabled')
        else:
            self.id_auto_file_display.config(state='normal')
            self.id_auto_file_display.delete(0, tk.END)
            self.id_auto_file_display.insert(0, "No cleaned data available")
            self.id_auto_file_display.config(state='readonly')
        
        # Memory data is logged in console output instead of GUI display
    
    def run_id_annotation(self, input_file):
        """Run ID annotation in background thread"""
        try:
            # Install logging handler for GUI updates
            import logging
            
            class GUILogHandler(logging.Handler):
                def __init__(self, append_callback):
                    super().__init__()
                    self.append_callback = append_callback
                    
                def emit(self, record):
                    try:
                        msg = self.format(record)
                        # Schedule GUI update in main thread
                        self.append_callback(msg)
                    except Exception:
                        pass
            
            # Create handler that writes to ID annotation log widget
            log_widget = getattr(self, 'id_progress_text', None)
            if log_widget:
                gui_handler = GUILogHandler(lambda msg: self.root.after(0, lambda: self.append_to_id_log(msg)))
                gui_handler.setLevel(logging.INFO)
                gui_handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
                
                # Add handler to root logger
                root_logger = logging.getLogger()
                root_logger.addHandler(gui_handler)
            
            # Get mode and workers from runtime variables set by start_id_annotation
            mode = getattr(self, '_annotation_mode_runtime', 'metabolite')
            workers = getattr(self, '_annotation_workers_runtime', self._get_workers_count(self.id_workers))
            
            # Prepare configuration
            output_file = self.id_output_file.get() if mode == 'metabolite' else self.lipid_output_file.get()
            
            # Ensure output file has the correct extension
            if not output_file.endswith('.xlsx'):
                output_file += '.xlsx'
            
            # Create annotator
            self.root.after(0, lambda: self.update_id_progress_with_percentage(10, "Initializing ID annotator..."))
            
            # Define progress callback function with aggressive throttling and GUI event processing
            last_update_time = [0]  # Use list to allow modification in nested function
            update_counter = [0]
            last_percentage = [0]  # Track last percentage to ensure visual updates
            
            def progress_callback(i, total, msg):
                # Throttling: update when percentage changes or every 100ms minimum
                current_time = time.time()
                percentage = int((i / total) * 100) if total > 0 else 0
                
                # Always update if percentage changed, otherwise throttle by time
                percentage_changed = percentage != last_percentage[0]
                time_elapsed = current_time - last_update_time[0] >= 0.1  # 100ms = 10 updates/sec max
                
                if not percentage_changed and not time_elapsed and i < total:
                    return  # Skip this update
                
                last_update_time[0] = current_time
                last_percentage[0] = percentage
                update_counter[0] += 1
                
                # Extract metabolite name from message
                current_metabolite = "Unknown"
                if msg and "Completed:" in msg:
                    try:
                        current_metabolite = msg.split("Completed:")[1].strip()
                    except:
                        current_metabolite = "Unknown"
                
                # Build progress message - don't show i/total for percentage-based progress (Phase 1-3)
                # Instead, show the phase message directly
                if "Phase" in str(msg):
                    # Phase-based progress: show message as-is
                    detailed_msg = str(msg)
                else:
                    # Regular progress: show count/total with percentage
                    detailed_msg = f"🔬 {i}/{total} ({percentage}%)"
                    if current_metabolite != "Unknown" and len(current_metabolite) < 50:
                        detailed_msg += f" - {current_metabolite}"
                
                # Estimate time remaining (only show occasionally to reduce overhead)
                if hasattr(self, 'id_start_time') and self.id_start_time and update_counter[0] % 5 == 0:
                    elapsed = time.time() - self.id_start_time
                    if i > 0 and elapsed > 0:
                        rate = i / elapsed
                        remaining = (total - i) / rate if rate > 0 else 0
                        remaining_min = int(remaining // 60)
                        remaining_sec = int(remaining % 60)
                        detailed_msg += f" | ETA: {remaining_min:02d}:{remaining_sec:02d}"
                
                # Thread-safe GUI update using root.after() - CRITICAL for background thread updates
                try:
                    # Capture variables for lambda closure
                    pct = percentage
                    dmsg = detailed_msg
                    self.root.after(0, lambda p=pct, m=dmsg: self.update_id_progress_with_percentage(p, m))
                except Exception as e:
                    pass  # Fail silently if GUI update fails
            
            # Define stop check callback function
            def stop_check_callback():
                """Check if user requested stop"""
                return getattr(self, 'id_annotation_stop_flag', False)
            
            # Build annotator with selected mode and options
            lipid_mode = (mode == 'lipid')
            
            # Get DataFrames from memory store for merging
            cleaned_metabolites_df = None
            metabolite_ids_df = None
            pos_df = None
            neg_df = None
            
            if hasattr(self, 'memory_store') and self.memory_store:
                print("\n" + "="*80)
                print("🔍 DEBUG: LOADING DATA FROM MEMORY STORE FOR ID ANNOTATION")
                print("="*80)

                if lipid_mode:
                    # Lipid-specific keys populated by Data Cleaning (lipid workflow)
                    cleaned_metabolites_df = self.memory_store.get('lipid_combined_df')
                    if isinstance(cleaned_metabolites_df, pd.DataFrame) and not cleaned_metabolites_df.empty:
                        print(f"✅ LIPID COMBINED SHEET LOADED: {len(cleaned_metabolites_df)} lipids")
                        cleaned_metabolites_df = cleaned_metabolites_df.copy()
                    else:
                        print("❌ LIPID COMBINED SHEET: Not found or empty")
                        cleaned_metabolites_df = None

                    pos_df = self.memory_store.get('lipid_pos_df')
                    if isinstance(pos_df, pd.DataFrame) and not pos_df.empty:
                        print(f"✅ LIPID POSITIVE SHEET LOADED: {len(pos_df)} lipids")
                    else:
                        print("ℹ️  LIPID POSITIVE SHEET: Not found or empty")
                        pos_df = None

                    neg_df = self.memory_store.get('lipid_neg_df')
                    if isinstance(neg_df, pd.DataFrame) and not neg_df.empty:
                        print(f"✅ LIPID NEGATIVE SHEET LOADED: {len(neg_df)} lipids")
                    else:
                        print("ℹ️  LIPID NEGATIVE SHEET: Not found or empty")
                        neg_df = None

                    metabolite_ids_df = None  # not used in lipid-only workflow
                else:
                    # PRIORITY: Use Combined sheet (all metabolites) for ID annotation
                    cleaned_metabolites_df = self.memory_store.get('cleaned_metabolites_df')
                    if isinstance(cleaned_metabolites_df, pd.DataFrame) and not cleaned_metabolites_df.empty:
                        print(f"✅ COMBINED SHEET LOADED: {len(cleaned_metabolites_df)} metabolites")
                        print(f"   Columns: {list(cleaned_metabolites_df.columns)[:5]}... ({len(cleaned_metabolites_df.columns)} total)")
                        cleaned_metabolites_df = cleaned_metabolites_df.copy()
                    else:
                        print("❌ COMBINED SHEET: Not found or empty")
                        cleaned_metabolites_df = None

                    metabolite_ids_df = self.memory_store.get('metabolite_ids_df')
                    if isinstance(metabolite_ids_df, pd.DataFrame) and not metabolite_ids_df.empty:
                        print(f"✅ METABOLITE_IDS LOADED: {len(metabolite_ids_df)} entries")
                    else:
                        print("ℹ️  METABOLITE_IDS: Not found or empty")

                    # Get Positive and Negative DataFrames from memory store (for merging back later)
                    pos_df = self.memory_store.get('Positive')
                    if isinstance(pos_df, pd.DataFrame) and not pos_df.empty:
                        print(f"✅ POSITIVE SHEET LOADED: {len(pos_df)} metabolites (for merging results back)")
                        print(f"   Columns: {list(pos_df.columns)[:5]}... ({len(pos_df.columns)} total)")
                    else:
                        print("ℹ️  POSITIVE SHEET: Not found or empty")
                        pos_df = None

                    neg_df = self.memory_store.get('Negative')
                    if isinstance(neg_df, pd.DataFrame) and not neg_df.empty:
                        print(f"✅ NEGATIVE SHEET LOADED: {len(neg_df)} metabolites (for merging results back)")
                        print(f"   Columns: {list(neg_df.columns)[:5]}... ({len(neg_df.columns)} total)")
                    else:
                        print("ℹ️  NEGATIVE SHEET: Not found or empty")
                        neg_df = None

                print("\n🎯 PRIORITY FOR ID ANNOTATION:")
                if cleaned_metabolites_df is not None:
                    label = "lipids" if lipid_mode else "metabolites"
                    print(f"   1️⃣  COMBINED sheet will be used: {len(cleaned_metabolites_df)} {label}")
                    if not lipid_mode:
                        print(f"   2️⃣  Results will be merged back to Positive and Negative sheets")
                elif pos_df is not None:
                    label = "lipids" if lipid_mode else "metabolites"
                    print(f"   ⚠️  Fallback: Positive sheet will be used: {len(pos_df)} {label}")
                else:
                    print("   ⚠️  No data found - will try loading from file")
                print("="*80 + "\n")

                if cleaned_metabolites_df is not None:
                    print(f"🧠 Passing cleaned DataFrame to annotator: {len(cleaned_metabolites_df)} rows; columns: {list(cleaned_metabolites_df.columns)}")
                else:
                    print("🧠 No cleaned DataFrame found in memory store")

                if metabolite_ids_df is not None:
                    print(f"🧠 Passing metabolite IDs DataFrame to annotator: {len(metabolite_ids_df)} entries")
                else:
                    if not lipid_mode:
                        print("🧠 No metabolite IDs DataFrame found in memory store")

                if pos_df is not None and isinstance(pos_df, pd.DataFrame):
                    print(f"🧠 Passing Positive DataFrame to annotator: {len(pos_df)} rows, {len(pos_df.columns)} columns")
                else:
                    print("🧠 No Positive DataFrame found in memory store")

                if neg_df is not None and isinstance(neg_df, pd.DataFrame):
                    print(f"🧠 Passing Negative DataFrame to annotator: {len(neg_df)} rows, {len(neg_df.columns)} columns")
                else:
                    print("🧠 No Negative DataFrame found in memory store")
            else:
                print("🧠 No memory store available")
            
            # Check if custom ID search mode is enabled
            custom_mode = getattr(self, '_custom_id_search_mode', False)
            
            if custom_mode:
                # Custom ID Search Mode - no Pos/Neg sheets required
                print("🔍 Custom ID Search Mode enabled - loading first available sheet")
                self.root.after(0, lambda: self.update_id_progress_with_percentage(10, "Custom ID Search Mode..."))
                
                try:
                    # Check if column mapping was verified
                    mode_name = "lipid" if lipid_mode else "metabolite"
                    column_mapping = getattr(self, '_custom_id_column_mapping', None) if not lipid_mode else getattr(self, '_custom_lipid_column_mapping', None)
                    verified_df = getattr(self, '_custom_id_dataframe', None) if not lipid_mode else getattr(self, '_custom_lipid_dataframe', None)
                    selected_sheet = getattr(self, '_custom_id_selected_sheet', None) if not lipid_mode else getattr(self, '_custom_lipid_selected_sheet', None)
                    
                    print(f"🔍 Custom mode debug ({mode_name}):")
                    print(f"   - Column mapping found: {column_mapping is not None}")
                    print(f"   - Verified dataframe found: {verified_df is not None}")
                    if column_mapping:
                        print(f"   - Mappings: {column_mapping}")
                    
                    if verified_df is not None and column_mapping:
                        # Use the verified dataframe - columns are already renamed by the dialog
                        print(f"✅ Using verified dataframe (columns already renamed by dialog)")
                        custom_df = verified_df.copy()
                        print(f"✅ Dataframe columns: {list(custom_df.columns)}")
                    else:
                        # No verification - load from file
                        print("⚠️ No column verification - loading raw sheet (may fail if columns don't match expected names)")
                        xl = pd.ExcelFile(input_file)
                        if len(xl.sheet_names) == 0:
                            raise ValueError("No sheets found in Excel file")
                        
                        # Load first sheet only
                        first_sheet = selected_sheet if selected_sheet else xl.sheet_names[0]
                        custom_df = pd.read_excel(input_file, sheet_name=first_sheet)
                        print(f"✅ Loaded sheet '{first_sheet}': {len(custom_df)} rows, {len(custom_df.columns)} columns")
                    
                    # In custom mode, we treat this as both pos and neg for simplicity
                    # The annotator will just do ID lookup without filtering/merging
                    pos_df = custom_df
                    neg_df = None  # Single sheet mode
                    
                    # Disable filtering in custom mode
                    selected_id_cols = []
                    effective_filter_mode = False
                    id_filter_mode = 'none'
                    ms2_filter_mode = 'none'
                    confidence_filter_mode = getattr(self, 'id_confidence_filter_var', tk.StringVar(value='exclude_low')).get()
                    require_endogenous = False
                    
                    print("🔧 Custom mode settings:")
                    print("   - Single sheet ID lookup only")
                    print("   - No Pos/Neg sheet requirement")
                    print("   - No filtering/cleaning/merging")
                    print(f"   - Confidence filter: {confidence_filter_mode}")
                    print("   - Returns IDs only")
                    
                except Exception as e:
                    error_msg = f"Failed to load custom sheet: {str(e)}"
                    print(f"❌ {error_msg}")
                    self.root.after(0, lambda: self.handle_id_annotation_error(error_msg))
                    return
                    
            else:
                # Standard Mode - try to load Pos/Neg sheets
                # If Positive/Negative/Combined not in memory, try loading from input Excel file
                if (pos_df is None or neg_df is None or cleaned_metabolites_df is None) and input_file and os.path.exists(input_file):
                    try:
                        print(f"🔍 Loading sheets from input file: {input_file}")
                        xl = pd.ExcelFile(input_file)
                        print(f"   Available sheets: {xl.sheet_names}")
                        
                        # PRIORITY: Load Combined sheet first for ID annotation
                        if 'Combined' in xl.sheet_names and cleaned_metabolites_df is None:
                            cleaned_metabolites_df = pd.read_excel(input_file, sheet_name='Combined')
                            print(f"✅ Loaded Combined sheet from file: {len(cleaned_metabolites_df)} metabolites, {len(cleaned_metabolites_df.columns)} columns")
                            print(f"   This will be used for ID annotation (all unique metabolites)")
                        
                        if 'Positive' in xl.sheet_names and pos_df is None:
                            pos_df = pd.read_excel(input_file, sheet_name='Positive')
                            print(f"✅ Loaded Positive sheet from file: {len(pos_df)} metabolites (for merging results back)")
                        if 'Negative' in xl.sheet_names and neg_df is None:
                            neg_df = pd.read_excel(input_file, sheet_name='Negative')
                            print(f"✅ Loaded Negative sheet from file: {len(neg_df)} metabolites (for merging results back)")
                    except Exception as e:
                        print(f"⚠️ Warning: Could not load sheets from file: {e}")

                # Get filter settings from runtime variables
                selected_id_cols = getattr(self, '_selected_id_cols_runtime', [])
                effective_filter_mode = getattr(self, '_effective_id_filter_mode_runtime', False)
                id_filter_mode = 'any_selected' if effective_filter_mode else 'none'
                ms2_filter_mode = getattr(self, 'id_ms2_filter_var', tk.StringVar(value='none')).get()
                confidence_filter_mode = getattr(self, 'id_confidence_filter_var', tk.StringVar(value='exclude_low')).get()
                require_endogenous = getattr(self, 'require_endogenous_yes', tk.BooleanVar(value=False)).get()
                
                # Log filter settings
                print(f"🔧 Filter settings:")
                print(f"   ID filter mode: {id_filter_mode}")
                if selected_id_cols:
                    print(f"   Selected ID columns: {', '.join(selected_id_cols)}")
                else:
                    print(f"   Selected ID columns: None (no filtering)")
                print(f"   MS2 filter mode: {ms2_filter_mode}")
                print(f"   Confidence filter mode: {confidence_filter_mode}")
                print(f"   Require Endogenous=Yes: {require_endogenous}")

            print("\n" + "="*80)
            print("🚀 DEBUG: PASSING DATA TO ID ANNOTATOR")
            print("="*80)
            print(f"📄 Input file: {input_file}")
            print(f"📄 Output file: {output_file}")
            print(f"👷 Workers: {workers}")
            print(f"🧬 Lipid mode: {lipid_mode}")
            print(f"\n📊 DataFrames being passed:")
            if cleaned_metabolites_df is not None:
                print(f"   ✅ cleaned_metabolites_df (COMBINED): {len(cleaned_metabolites_df)} rows, {len(cleaned_metabolites_df.columns)} cols")
            else:
                print(f"   ❌ cleaned_metabolites_df: None")
            if metabolite_ids_df is not None:
                print(f"   ✅ metabolite_ids_df: {len(metabolite_ids_df)} rows")
            else:
                print(f"   ℹ️  metabolite_ids_df: None")
            if pos_df is not None:
                print(f"   ✅ pos_df (POSITIVE): {len(pos_df)} rows, {len(pos_df.columns)} cols")
            else:
                print(f"   ❌ pos_df: None")
            if neg_df is not None:
                print(f"   ✅ neg_df (NEGATIVE): {len(neg_df)} rows, {len(neg_df.columns)} cols")
            else:
                print(f"   ❌ neg_df: None")
            print("="*80 + "\n")

            mapped_override_df = getattr(self, '_id_annotation_input_df', None)
            mapped_override_sheet = str(getattr(self, '_id_annotation_selected_sheet', '') or '').strip().lower()
            use_override_df = bool(custom_mode) or (mapped_override_sheet == 'combined')

            # If user explicitly selected a custom upload file, force file input over pre-loaded memory data.
            manual_selected_path = ''
            auto_loaded_path = ''
            try:
                manual_selected_path = (
                    self.lipid_id_annotation_file.get().strip() if lipid_mode and hasattr(self, 'lipid_id_annotation_file')
                    else self.id_input_file.get().strip() if hasattr(self, 'id_input_file')
                    else ''
                )
            except Exception:
                manual_selected_path = ''

            try:
                auto_loaded_path = (
                    str(getattr(self, 'cleaned_lipid_excel_path', '') or '').strip() if lipid_mode
                    else str(getattr(self, 'cleaned_excel_path', '') or '').strip()
                )
            except Exception:
                auto_loaded_path = ''

            norm_manual = os.path.normcase(os.path.abspath(manual_selected_path)) if manual_selected_path else ''
            norm_auto = os.path.normcase(os.path.abspath(auto_loaded_path)) if auto_loaded_path else ''
            force_file_input = bool(norm_manual and (not norm_auto or norm_manual != norm_auto))

            if force_file_input:
                print(f"📁 Custom upload selected -> force_file_input=True")
                print(f"   Manual file: {manual_selected_path}")
                if auto_loaded_path:
                    print(f"   Auto-loaded file: {auto_loaded_path}")
            else:
                print("🧠 Using pre-loaded memory dataframe priority (no custom override file selected)")

            if not custom_mode and mapped_override_df is not None and not use_override_df:
                print("ℹ️  Column mapping was done on a polarity sheet; merged output will still be built from the Pos/Neg workflow.")
                if hasattr(self, 'id_progress_text'):
                    self.root.after(0, lambda: self.id_progress_text.insert(
                        tk.END,
                        "ℹ️ Pos/Neg mappings were captured and will be used in the merged output workflow.\n"
                    ))
            
            annotator = MetaboliteIDAnnotator(
                input_file=input_file,
                output_file=output_file,
                progress_callback=progress_callback,
                stop_check_callback=stop_check_callback,        # Pass stop check callback
                lipid_mode=lipid_mode,
                skip_pubchem=False,          # PubChem searches not skipped
                skip_kegg=lipid_mode,        # Skip KEGG searches only for lipid mode
                skip_hmdb=lipid_mode,        # Skip HMDB searches only for lipid mode
                cleaned_metabolites_df=cleaned_metabolites_df,  # Pass DataFrame from memory
                input_df_override=(mapped_override_df if use_override_df else None),  # Preserve Combined-first workflow unless Combined/custom override is used
                metabolite_ids_df=metabolite_ids_df,            # Pass metabolite IDs for merging
                force_file_input=force_file_input,              # Custom upload must override pre-loaded memory
                pos_df=pos_df,                                  # Pass Positive DataFrame
                neg_df=neg_df,                                  # Pass Negative DataFrame
                id_filter_mode=id_filter_mode,                  # Pass ID filter mode
                selected_id_columns=selected_id_cols if not self._skip_id_filtering_runtime else [],  # Pass selected ID columns
                skip_id_filtering=self._skip_id_filtering_runtime,  # Pass skip ID filtering flag
                ms2_filter_mode=ms2_filter_mode,                # Pass MS2 filter setting
                confidence_filter_mode=confidence_filter_mode,  # Pass confidence filter setting
                require_endogenous_yes=require_endogenous       # Pass endogenous filter setting
            )
            
            # Run annotation
            self.root.after(0, lambda: self.update_id_progress_with_percentage(5, "Starting ID annotation..."))
            
            if lipid_mode:
                annotator.run_lipid_id_annotation(max_workers=workers)
            else:
                annotator.run_id_annotation(max_workers=workers)
            
            # If custom mode, save only Merged_IDs sheet
            custom_mode = getattr(self, '_custom_id_search_mode', False)
            if custom_mode and os.path.exists(output_file):
                try:
                    print("🔧 Custom mode: Extracting only Merged_IDs sheet...")
                    xl = pd.ExcelFile(output_file)
                    if 'Merged_IDs' in xl.sheet_names:
                        merged_df = xl.parse('Merged_IDs')
                        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                            merged_df.to_excel(writer, sheet_name='Merged_IDs', index=False)
                        print("✅ Saved only Merged_IDs sheet")
                    else:
                        print("⚠️ Merged_IDs sheet not found, keeping all sheets")
                except Exception as e:
                    print(f"⚠️ Warning: Could not filter to Merged_IDs only: {e}")
            
            # Check if stopped by user
            if getattr(self, 'id_annotation_stop_flag', False):
                # Annotation was stopped
                results = {
                    'success': False,
                    'stopped': True,
                    'message': 'ID annotation stopped by user'
                }
                self.root.after(0, lambda: self.handle_id_annotation_stopped())
                return
            
            # Store annotator reference for direct DataFrame access
            self._last_annotator_results = annotator
            
            # Success
            results = {
                'success': True,
                'output_file': output_file,
                'input_file': input_file
            }
            
            self.root.after(0, lambda: self.handle_id_annotation_results(results))
            
        except Exception as e:
            error_msg = str(e)
            self.root.after(0, lambda: self.handle_id_annotation_error(error_msg))
    
    def update_id_timer(self):
        """Update ID annotation timer"""
        if getattr(self, '_shutting_down', False):
            return
        if hasattr(self, 'id_start_time') and self.id_start_time is not None:
            elapsed = time.time() - self.id_start_time
            minutes = int(elapsed // 60)
            seconds = int(elapsed % 60)
            self.id_timer_label.config(text=f"Time: {minutes:02d}:{seconds:02d}")
            
            # Calculate ETA
            if hasattr(self, 'id_current_step') and self.id_current_step > 0:
                rate = elapsed / self.id_current_step
                remaining_steps = self.id_total_steps - self.id_current_step
                eta_seconds = remaining_steps * rate
                eta_minutes = int(eta_seconds // 60)
                eta_seconds = int(eta_seconds % 60)
                self.id_eta_label.config(text=f"ETA: {eta_minutes:02d}:{eta_seconds:02d}")
            
            # Schedule next update
            self.root.after(1000, self.update_id_timer)
    
    def append_to_id_log(self, message):
        """Append message to ID annotation progress log - thread-safe"""
        if hasattr(self, 'id_progress_text'):
            try:
                self.id_progress_text.config(state='normal')
                self.id_progress_text.insert(tk.END, message + "\n")
                self.id_progress_text.see(tk.END)
                self.id_progress_text.config(state='disabled')
            except Exception as e:
                print(f"Error appending to ID log: {e}")
    
    def update_id_progress_with_percentage(self, percentage, message):
        """Update ID annotation progress with percentage"""
        self.id_current_step = percentage
        self.id_progress_bar['value'] = percentage  # Force update progress bar value
        self.id_progress_percentage.config(text=f"{percentage}%")
        
        # Update progress text
        self.id_progress_text.insert(tk.END, f"[{percentage}%] {message}\n")
        self.id_progress_text.see(tk.END)
        
        # Only update status for major milestones (not every progress update)
        if percentage in [0, 5, 100] or "Error" in message or "completed" in message.lower():
            status_msg = message[:50] + "..." if len(message) > 50 else message
            self.id_process_status.config(text=status_msg, fg='#3498db')
        
        # Force GUI update for progress bar visibility
        try:
            self.id_progress_bar.update()
        except Exception:
            pass
    
    def update_id_progress(self, message):
        """Update ID annotation progress text"""
        self.id_progress_text.insert(tk.END, message + "\n")
        self.id_progress_text.see(tk.END)
        
        # Only update status for important messages (not routine progress updates)
        if any(keyword in message.lower() for keyword in ["error", "completed", "failed", "success", "starting", "finished"]):
            status_msg = message[:50] + "..." if len(message) > 50 else message
            self.id_process_status.config(text=status_msg, fg='#3498db')
        
        # Removed update_idletasks() - let natural event loop handle updates
    
    def handle_id_annotation_stopped(self):
        """Handle ID annotation stopped by user"""
        # Reset UI
        self.id_start_button.config(state='normal', text="🚀 Start ID Annotation")
        self.id_stop_button.config(state='disabled')
        self.id_start_time = None
        
        self.id_process_status.config(text="⛔ ID annotation stopped by user", fg='#e74c3c')
        self.id_progress_text.insert(tk.END, "\n⛔ ID annotation process stopped by user\n")
        self.id_progress_text.see(tk.END)
        
        # Reset ETA
        self.id_eta_label.config(text="ETA: 00:00")
        self.id_timer_label.config(text="Time: 00:00")
        
        messagebox.showinfo("Stopped", "ID annotation process was stopped by user")
    
    def handle_id_annotation_results(self, results):
        """Handle ID annotation results"""
        # Reset UI
        self.id_start_button.config(state='normal', text="🚀 Start ID Annotation")
        self.id_stop_button.config(state='disabled')
        self.id_start_time = None
        
        if results['success']:
            self.id_annotated_data = results
            self.id_process_status.config(text="✅ ID annotation completed!", fg='#27ae60')
            self.update_id_progress_with_percentage(100, "✅ ID annotation completed!")
            
            # Reset ETA
            self.id_eta_label.config(text="ETA: 00:00")
            
            # Store data for pathway annotation (Tab 3)
            self.annotated_data = results['output_file']
            self.id_annotated_excel_path = results['output_file']
            # Ensure legacy/debug attribute exists so other parts of the GUI
            # (and the debug report) can detect annotated metabolites data.
            # Try to load the annotated Excel into a DataFrame; fall back to
            # storing the path if loading fails.
            try:
                if isinstance(results.get('output_file'), str) and os.path.exists(results['output_file']):
                    df_ann = pd.read_excel(results['output_file'])
                    self.annotated_metabolites_data = df_ann
                    # also keep a path reference for compatibility
                    self.annotated_metabolites_excel_path = results['output_file']
                else:
                    # Store the path (or whatever output_file holds) so hasattr checks pass
                    self.annotated_metabolites_data = results['output_file']
                    self.annotated_metabolites_excel_path = results.get('output_file')
            except Exception as e:
                # Non-fatal: keep the path so debug shows available, but warn in console
                print(f"Warning: could not load annotated metabolites file: {e}")
                self.annotated_metabolites_data = results.get('output_file')
                self.annotated_metabolites_excel_path = results.get('output_file')
            
            # Capture polarity-specific ID DataFrames directly from results if available
            try:
                # Check if annotator exposed posneg_results in the results
                if hasattr(self, '_last_annotator_results') and self._last_annotator_results:
                    posneg_results = getattr(self._last_annotator_results, 'posneg_results', None)
                    print(f"🔍 DEBUG: posneg_results type: {type(posneg_results)}")
                    if isinstance(posneg_results, dict):
                        print(f"🔍 DEBUG: posneg_results keys: {list(posneg_results.keys())}")
                    pos_df_out = None
                    neg_df_out = None
                    if isinstance(posneg_results, dict):
                        pos_df_out = posneg_results.get('pos_id_df')
                        neg_df_out = posneg_results.get('neg_id_df')
                        print(f"🔍 DEBUG: pos_df_out type: {type(pos_df_out)}, is DataFrame: {isinstance(pos_df_out, pd.DataFrame)}")
                        if isinstance(pos_df_out, pd.DataFrame):
                            print(f"🔍 DEBUG: pos_df_out shape: {pos_df_out.shape}, empty: {pos_df_out.empty}")
                        print(f"🔍 DEBUG: neg_df_out type: {type(neg_df_out)}, is DataFrame: {isinstance(neg_df_out, pd.DataFrame)}")
                        if isinstance(neg_df_out, pd.DataFrame):
                            print(f"🔍 DEBUG: neg_df_out shape: {neg_df_out.shape}, empty: {neg_df_out.empty}")
                        
                        if isinstance(pos_df_out, pd.DataFrame) and not pos_df_out.empty and hasattr(self, 'memory_store'):
                            self.memory_store['pos_id_df'] = pos_df_out
                            print(f"✅ Captured pos_id_df directly: {len(pos_df_out)} rows, {len(pos_df_out.columns)} columns")
                        elif isinstance(pos_df_out, pd.DataFrame) and pos_df_out.empty:
                            print(f"⚠️ pos_df_out is EMPTY - not storing in memory")
                        
                        if isinstance(neg_df_out, pd.DataFrame) and not neg_df_out.empty and hasattr(self, 'memory_store'):
                            self.memory_store['neg_id_df'] = neg_df_out
                            print(f"✅ Captured neg_id_df directly: {len(neg_df_out)} rows, {len(neg_df_out.columns)} columns")
                        elif isinstance(neg_df_out, pd.DataFrame) and neg_df_out.empty:
                            print(f"⚠️ neg_df_out is EMPTY - not storing in memory")

                    # If direct capture failed (either missing or empty), fallback to Excel sheets
                    if (not isinstance(pos_df_out, pd.DataFrame) or pos_df_out.empty) or (not isinstance(neg_df_out, pd.DataFrame) or neg_df_out.empty):
                        print("ℹ️ Falling back to reading Pos_id/Neg_id from Excel output or splitting merged sheet.")
                        if isinstance(results.get('output_file'), str) and os.path.exists(results['output_file']):
                            import pandas as _pd_capture
                            try:
                                xl_capture = _pd_capture.ExcelFile(results['output_file'])
                                found_any = False
                                if 'Pos_id' in xl_capture.sheet_names:
                                    pos_id_df = xl_capture.parse('Pos_id')
                                    if hasattr(self, 'memory_store') and isinstance(self.memory_store, dict):
                                        self.memory_store['pos_id_df'] = pos_id_df
                                        found_any = True
                                        print(f"📁 Loaded pos_id_df from Excel: {len(pos_id_df)} rows, {len(pos_id_df.columns)} columns")
                                if 'Neg_id' in xl_capture.sheet_names:
                                    neg_id_df = xl_capture.parse('Neg_id')
                                    if hasattr(self, 'memory_store') and isinstance(self.memory_store, dict):
                                        self.memory_store['neg_id_df'] = neg_id_df
                                        found_any = True
                                        print(f"📁 Loaded neg_id_df from Excel: {len(neg_id_df)} rows, {len(neg_id_df.columns)} columns")
                                # If Pos_id/Neg_id not present, try Positive/Negative
                                if not found_any:
                                    if 'Positive' in xl_capture.sheet_names:
                                        pos_df = xl_capture.parse('Positive')
                                        self.memory_store['pos_id_df'] = pos_df
                                        found_any = True
                                        print(f"📁 Loaded Positive sheet: {len(pos_df)} rows, {len(pos_df.columns)} columns")
                                    if 'Negative' in xl_capture.sheet_names:
                                        neg_df = xl_capture.parse('Negative')
                                        self.memory_store['neg_id_df'] = neg_df
                                        found_any = True
                                        print(f"📁 Loaded Negative sheet: {len(neg_df)} rows, {len(neg_df.columns)} columns")
                                # If still not found, split the merged sheet by Polarity if available
                                if not found_any:
                                    try:
                                        merged_sheet = 'Merged_IDs'
                                        if merged_sheet in xl_capture.sheet_names:
                                            merged_df = xl_capture.parse(merged_sheet)
                                            if 'Polarity' in merged_df.columns:
                                                pos_df = merged_df[merged_df['Polarity'].astype(str).str.lower().isin(['+', 'positive'])]
                                                neg_df = merged_df[merged_df['Polarity'].astype(str).str.lower().isin(['-', 'negative'])]
                                                if not pos_df.empty:
                                                    self.memory_store['pos_id_df'] = pos_df
                                                    found_any = True
                                                    print(f"🔀 Split merged sheet into Positive: {len(pos_df)} rows")
                                                if not neg_df.empty:
                                                    self.memory_store['neg_id_df'] = neg_df
                                                    found_any = True
                                                    print(f"🔀 Split merged sheet into Negative: {len(neg_df)} rows")
                                    except Exception as _split_err:
                                        print(f"Warning: could not split merged sheet: {_split_err}")
                            except Exception as _cap_err:
                                print(f"Warning: could not load polarity ID DataFrames from Excel: {_cap_err}")
            except Exception as _outer_cap_err:
                print(f"Warning: polarity capture block failed: {_outer_cap_err}")
            
            # Auto-load data into Tab 3 (Pathway Annotation)
            try:
                self.pathway_input_file.set(results['output_file'])
                # Auto-detect sheets for pathway annotation
                self.auto_detect_sheets(self.pathway_input_file, self.pathway_sheet_combo, self.pathway_sheet_var)
            except Exception as e:
                print(f"Could not auto-load data into Tab 3: {str(e)}")
            
            # Show completion message with next steps guidance
            annotation_mode = getattr(self, '_annotation_mode_runtime', 'metabolite')
            custom_mode = getattr(self, '_custom_id_search_mode', False)
            
            # Add helpful guidance to progress log
            self.id_progress_text.insert(tk.END, "\n" + "="*60 + "\n")
            self.id_progress_text.insert(tk.END, "📋 NEXT STEPS FOR PATHWAY ANALYSIS\n")
            self.id_progress_text.insert(tk.END, "="*60 + "\n\n")
            self.id_progress_text.insert(tk.END, "To proceed with pathway analysis, you need to:\n\n")
            self.id_progress_text.insert(tk.END, "1️⃣  Run STATISTICS ANALYSIS (separate Statistics tool)\n")
            self.id_progress_text.insert(tk.END, "    • Load your ID-annotated file\n")
            self.id_progress_text.insert(tk.END, "    • Generate statistics (p-values, fold changes, log2FC)\n")
            self.id_progress_text.insert(tk.END, "    • This creates the required statistical output\n\n")
            self.id_progress_text.insert(tk.END, "2️⃣  Return to MetaboGraph - Pathway Analysis tab\n")
            self.id_progress_text.insert(tk.END, "    • Load the statistics output file\n")
            self.id_progress_text.insert(tk.END, "    • File must contain: p-value, FC, log2FC columns\n")
            self.id_progress_text.insert(tk.END, "    • Perform pathway enrichment analysis\n\n")
            self.id_progress_text.insert(tk.END, "💡 TIP: The statistics file links your IDs to biological\n")
            self.id_progress_text.insert(tk.END, "   significance (p-values) and magnitude (fold changes)\n")
            self.id_progress_text.insert(tk.END, "="*60 + "\n")
            self.id_progress_text.see(tk.END)
            
            messagebox.showinfo("✅ ID Annotation Complete", 
                              f"ID annotation completed successfully!\n\n"
                              f"Output: {results['output_file']}\n\n"
                              f"📋 NEXT STEPS:\n"
                              f"1. Run Statistics Analysis (separate tool)\n"
                              f"   to generate p-values, FC, and log2FC\n\n"
                              f"2. Load the statistics output in\n"
                              f"   Pathway Analysis tab for enrichment\n\n"
                              f"See Progress Log for detailed instructions.")
        else:
            error_msg = results.get('error', 'Unknown error occurred')
            self.id_process_status.config(text="❌ Error occurred", fg='#e74c3c')
            self.update_id_progress_with_percentage(0, f"❌ Error: {error_msg}")
    
    def browse_file(self, file_var, title):
        """Browse for a file - supports both open and save dialogs"""
        try:
            # Determine if this is for input (open) or output (save)
            is_output = 'result' in title.lower() or 'output' in title.lower()
            
            if is_output:
                # For output files, use asksaveasfilename to allow new files
                filename = filedialog.asksaveasfilename(
                    title=f"Select {title}",
                    filetypes=[("Excel files", "*.xlsx"), ("Excel 97-2003", "*.xls"), ("All files", "*.*")],
                    defaultextension=".xlsx"
                )
            else:
                # For input files, use askopenfilename to select existing files
                filename = filedialog.askopenfilename(
                    title=f"Select {title}",
                    filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
                )
            
            if filename:
                file_var.set(filename)
                self.root.update_idletasks()

                if file_var is getattr(self, 'id_input_file', None):
                    self._prompt_id_annotation_column_assignment(filename)
        except Exception as e:
            messagebox.showerror("Browse Error", f"Failed to open file dialog:\n{str(e)}")

    def _prompt_id_annotation_column_assignment(self, input_file, *, store_as_custom=False):
        """Open the ID annotation column assignment dialog for a selected workbook."""
        try:
            if not input_file or not os.path.exists(input_file):
                return None

            self._id_annotation_column_mapping = None
            self._id_annotation_input_df = None
            self._id_annotation_selected_sheet = None

            if store_as_custom:
                self._custom_id_column_mapping = None
                self._custom_id_dataframe = None
                self._custom_id_selected_sheet = None

            xl = pd.ExcelFile(input_file)
            if len(xl.sheet_names) == 0:
                messagebox.showerror("Error", "No sheets found in the selected Excel file")
                return None

            preview_sheet = xl.sheet_names[0]
            df = pd.read_excel(input_file, sheet_name=preview_sheet)

            self.id_progress_text.insert(tk.END, f"\n🔎 Opening column assignment dialog for: {os.path.basename(input_file)}\n")
            self.id_progress_text.insert(tk.END, f"   Preview sheet: {preview_sheet}\n")
            if len(xl.sheet_names) > 1:
                self.id_progress_text.insert(tk.END, f"   Available sheets: {', '.join(xl.sheet_names)}\n")
            self.id_progress_text.see(tk.END)

            result = show_column_assignment_dialog(
                parent=self.root,
                df=df,
                tab_type='id_annotation',
                auto_calculate=False,
                excel_file_path=input_file,
            )

            if not result:
                self.id_progress_text.insert(tk.END, "⚠️ Column assignment was cancelled. The annotation will continue using available column names.\n")
                self.id_progress_text.see(tk.END)
                return None

            assignments = result['assignments']
            assigned_df = result['dataframe']
            selected_sheet = result.get('selected_sheet')
            mapped_sheets = result.get('mapped_sheets', [])

            self._id_annotation_column_mapping = assignments
            self._id_annotation_input_df = assigned_df
            self._id_annotation_selected_sheet = selected_sheet

            self._custom_id_column_mapping = assignments
            self._custom_id_dataframe = assigned_df
            self._custom_id_selected_sheet = selected_sheet

            missing_roles = [
                role for role in ['MS2', 'MS2 Purity [%]']
                if not assignments.get(role)
            ]

            col_info = "\n".join([f"  • {k}: {v}" for k, v in assignments.items() if v])
            mapped_info = f"\n  • Mapped sheets: {', '.join(mapped_sheets)}" if mapped_sheets else ""
            sheet_info = f"\n  • Active run sheet: {selected_sheet}" if selected_sheet else ""
            self.id_progress_text.insert(tk.END, f"✅ Column assignment complete!{mapped_info}{sheet_info}\n{col_info}\n")

            if missing_roles:
                warning_text = (
                    "Some optional confidence columns were not assigned:\n"
                    f"{', '.join(missing_roles)}\n\n"
                    "Confidence logic that depends on the missing columns will be skipped.\n"
                    "The analysis will continue with the columns that are available."
                )
                messagebox.showwarning("Incomplete Column Assignment", warning_text)
                self.id_progress_text.insert(tk.END, f"⚠️ Missing assignments: {', '.join(missing_roles)}\n")

            self.id_progress_text.see(tk.END)
            return result

        except Exception as e:
            messagebox.showerror("Error", f"Failed to open column assignment dialog:\n{str(e)}")
            logger.error(f"ID annotation column assignment error: {e}", exc_info=True)
            return None
    
    def browse_output_directory(self):
        """Browse for output directory"""
        directory = filedialog.askdirectory(title="Select Output Directory")
        if directory:
            self.id_output_directory.set(directory)
    
    def browse_output_directory_lipid(self):
        """Browse for lipid output directory"""
        directory = filedialog.askdirectory(title="Select Output Directory")
        if directory:
            self.lipid_output_directory.set(directory)
    
    def handle_id_annotation_error(self, error_msg):
        """Handle ID annotation error"""
        # Reset UI
        self.id_progress_bar.config(value=0)
        self.id_progress_percentage.config(text="0%")
        self.id_start_button.config(state='normal', text="🚀 Start ID Annotation")
        self.id_stop_button.config(state='disabled')
        self.id_process_status.config(text="❌ Error occurred", fg='#e74c3c')
        self.id_start_time = None
        
        # Reset ETA
        self.id_eta_label.config(text="ETA: --:--")
        
        self.update_id_progress_with_percentage(0, f"❌ Error: {error_msg}")
        messagebox.showerror("Error", error_msg)
    
    def load_id_annotation_prefs(self):
        """Load ID annotation preferences from file"""
        import json
        prefs_file = os.path.join(os.getcwd(), 'id_annotation_prefs.json')
        default_prefs = {
            'id_selected_columns': [
                'HMDB_ID', 'PubChem_CID', 'KEGG_ID', 'LipidMaps_ID',
                'ChEBI_ID', 'CAS', 'InChI', 'InChIKey'
            ],
            'skip_id_filtering': False,
            'lipid_keep_all_rows': False,
            'ms2_filter': 'none',
            'confidence_filter': 'exclude_low',
            'require_endogenous_yes': False,
            'dedup_rt_window_minutes': 2.0
        }
        
        try:
            if os.path.exists(prefs_file):
                with open(prefs_file, 'r') as f:
                    loaded = json.load(f)
                    # Merge with defaults to handle missing keys
                    default_prefs.update(loaded)
                return default_prefs
        except Exception as e:
            print(f"Warning: Could not load preferences from {prefs_file}: {e}")
        
        return default_prefs
    
    def save_id_annotation_prefs(self):
        """Save ID annotation preferences to file"""
        import json
        prefs_file = os.path.join(os.getcwd(), 'id_annotation_prefs.json')
        
        # Collect current preferences from UI
        prefs = {
            'id_selected_columns': [col for col, var in self.id_filter_vars.items() if var.get()],
            'skip_id_filtering': self.skip_id_filtering_var.get() if hasattr(self, 'skip_id_filtering_var') else False,
            'lipid_keep_all_rows': self.lipid_keep_all_rows_var.get() if hasattr(self, 'lipid_keep_all_rows_var') else False,
            'ms2_filter': self.id_ms2_filter_var.get() if hasattr(self, 'id_ms2_filter_var') else 'none',
            'confidence_filter': self.id_confidence_filter_var.get() if hasattr(self, 'id_confidence_filter_var') else 'exclude_low',
            'require_endogenous_yes': self.require_endogenous_yes.get() if hasattr(self, 'require_endogenous_yes') else False,
            'dedup_rt_window_minutes': float(self.dedup_rt_window.get()) if hasattr(self, 'dedup_rt_window') else 2.0
        }
        
        try:
            with open(prefs_file, 'w') as f:
                json.dump(prefs, f, indent=2)
            print(f"✓ Saved preferences to {prefs_file}")
        except Exception as e:
            print(f"Warning: Could not save preferences to {prefs_file}: {e}")
    
    def update_id_selected_count(self):
        """Update the label showing how many ID columns are selected and reflect implicit filtering rule."""
        try:
            selected = [c for c, v in getattr(self, 'id_filter_vars', {}).items() if v.get()]
            total = len(getattr(self, 'id_filter_vars', {}))
            if hasattr(self, 'id_selected_count_label'):
                if len(selected) == 0:
                    self.id_selected_count_label.config(
                        text=f"0 of {total} ID columns selected (no ID-based row filtering)",
                        fg='#c0392b'
                    )
                else:
                    self.id_selected_count_label.config(
                        text=f"{len(selected)} of {total} ID columns selected → rows will be filtered to those having at least one.",
                        fg='#2c3e50'
                    )
        except Exception as e:
            # Fail silently but log to console for debugging
            print(f"update_id_selected_count error: {e}")
    
    def _toggle_id_filter_controls(self):
        """Enable/disable ID filter controls based on skip_id_filtering_var state."""
        skip_filtering = getattr(self, 'skip_id_filtering_var', tk.BooleanVar(value=False)).get()
        
        # Helper function to get all descendants of a widget
        def get_all_descendants(widget):
            descendants = []
            for child in widget.winfo_children():
                descendants.append(child)
                descendants.extend(get_all_descendants(child))
            return descendants
        
        # Find the frame and disable/enable all its children
        try:
            for widget in get_all_descendants(self.frame):
                # Disable ID column checkboxes and buttons
                if isinstance(widget, tk.Checkbutton) and hasattr(widget, 'cget'):
                    text = widget.cget('text')
                    if text in self.available_id_columns:
                        widget.config(state='disabled' if skip_filtering else 'normal')
                elif isinstance(widget, tk.Button) and hasattr(widget, 'cget'):
                    text = widget.cget('text')
                    if text in ['Select All', 'Deselect All']:
                        widget.config(state='disabled' if skip_filtering else 'normal')
            
            # Update label
            if hasattr(self, 'id_selected_count_label'):
                if skip_filtering:
                    self.id_selected_count_label.config(
                        text="ID filtering is disabled (all metabolites will be saved)",
                        fg='#27ae60'
                    )
                else:
                    self.update_id_selected_count()
        except Exception as e:
            print(f"Error toggling ID filter controls: {e}")
    
    def auto_detect_sheets(self, file_var, combo_widget, sheet_var):
        """Auto-detect and populate sheet names in a combobox from an Excel file"""
        try:
            file_path = file_var.get() if hasattr(file_var, 'get') else str(file_var)
            if file_path and os.path.exists(file_path):
                xl = pd.ExcelFile(file_path)
                sheet_names = xl.sheet_names
                if combo_widget and sheet_names:
                    combo_widget['values'] = sheet_names
                    if sheet_names:
                        combo_widget.current(0)
                        sheet_var.set(sheet_names[0]) if hasattr(sheet_var, 'set') else None
        except Exception as e:
            print(f"Warning: Could not auto-detect sheets: {e}")
    
    def auto_prepare_statistics_tab(self):
        """Prepare the statistics tab with cleaned/annotated data and auto-switch"""
        try:
            # Notify Statistics tab that annotated data is ready
            self.notify_data_ready("📈 Statistics", "annotated_metabolites")
            
            # Get Statistics tab and auto-load the annotated data
            stats_tab = self.get_tab_by_name("📈 Statistics")
            if stats_tab:
                # Store the annotated file path in Statistics tab's memory
                if hasattr(self, 'id_annotated_excel_path'):
                    stats_tab.memory_store['annotated_excel_path'] = self.id_annotated_excel_path
                    print(f"✅ Stored annotated file path in Statistics tab: {self.id_annotated_excel_path}")
                
                # If Statistics tab has auto-load method, call it
                if hasattr(stats_tab, 'auto_load_annotated_file'):
                    stats_tab.auto_load_annotated_file(self.id_annotated_excel_path)
                    print(f"✅ Auto-loaded annotated data to Statistics tab")
                
                # Switch to Statistics tab after a short delay
                self.root.after(500, lambda: self.switch_to_tab("📈 Statistics"))
                print(f"✅ Switched to Statistics tab for analysis")
            else:
                print(f"⚠️ Statistics tab not found")
        except Exception as e:
            print(f"Warning: Could not auto-prepare statistics tab: {e}")

