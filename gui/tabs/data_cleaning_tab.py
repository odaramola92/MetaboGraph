import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import logging
import sys
import os
import stat
import pandas as pd
import threading
import time
import traceback

# Import backend and shared components
from gui.shared.base_tab import BaseTab
from gui.shared.utils import TimerManager, SessionManager, ProgressUpdater
from gui.shared.column_assignment import show_column_assignment_dialog, TAB_REQUIREMENTS
from main_script.metabolite_data_cleaner import CleanColumnDataCleaner as DataCleaner

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DataCleaningTab(BaseTab):
    """Data Cleaning Tab - Complete Standalone Implementation
    
    Provides full metabolite and lipid data cleaning functionality.
    """
    
    def __init__(self, parent, data_manager):
        """Initialize the Data Cleaning tab"""
        super().__init__(parent, data_manager)
        
        # Get root window for dialogs
        self.root = parent.winfo_toplevel()
        
        # Initialize all variables BEFORE calling setup_ui
        self._initialize_variables()
        
        # Thread tracking
        self.cleaning_thread = None
        self.lipid_cleaning_thread = None
        
        # Timer and progress tracking using shared utilities
        self.timer_manager = TimerManager()
        self.session_manager = SessionManager()
        
        # Create UI
        self.setup_ui()
        print("[OK] Data Cleaning Tab initialized (STANDALONE)")
    
    def _initialize_variables(self):
        """Initialize all StringVar and BooleanVar before UI creation"""
        # Metabolite file variables
        self.cleaning_mode = tk.StringVar(value="separate")
        self.combined_file = tk.StringVar()
        self.positive_file = tk.StringVar()
        self.negative_file = tk.StringVar()
        self.metabolika_file = tk.StringVar()
        self.mzcloud_pos_file = tk.StringVar()
        self.mzcloud_neg_file = tk.StringVar()
        self.use_negative_mode = tk.BooleanVar(value=True)
        
        # Output settings - single file path for output
        self.metabolite_output_file = tk.StringVar(value=os.path.join(os.getcwd(), "cleaned_metabolites.xlsx"))
        
        # Enhanced cleaning options
        self.use_enhanced_cleaning = tk.BooleanVar(value=True)
        self.ppm_threshold = tk.StringVar(value="10.0")
        self.rt_threshold = tk.StringVar(value="4.0")
        self.ms2_filter_mode = tk.StringVar(value='none')
        
        # Column cleaning patterns
        self.column_clean_patterns = tk.StringVar(
            value="Group Area:, _pos, _neg, .raw, Area[, ], [, ], pos, pos_, neg, neg_, positive, negative"
        )
        
        # Lipid file variables
        self.lipid_cleaning_mode = tk.StringVar(value="separate")  # combined or separate
        self.lipid_combined_file = tk.StringVar()  # For combined mode
        self.lipid_positive_file = tk.StringVar()  # For separate mode
        self.lipid_negative_file = tk.StringVar()  # For separate mode
        self.lipid_input_file = tk.StringVar()  # Keep for backward compatibility
        self.lipid_sample_details_file = tk.StringVar()  # Optional sample details mapping (for combined mode)
        self.lipid_positive_sample_details_file = tk.StringVar()  # File ID for positive mode
        self.lipid_negative_sample_details_file = tk.StringVar()  # File ID for negative mode
        self.lipid_output_file = tk.StringVar(value=os.path.join(os.getcwd(), "cleaned_lipids.xlsx"))
        self.lipid_column_clean_patterns = tk.StringVar(value="Area[, ], [, ], _pos, _neg, .raw, _Area, Area_")
        
        # Reference ion checkboxes
        self.reference_ion_vars = {}
        self.lipid_adduct_vars = {}
        
        # Initialize reference ions (all 10 ions, checked by default)
        REFERENCE_IONS_ALL = [
            '[M+H]+1',      # Positive
            '[M-H]-1',      # Negative
            '[M+2H]+2',     # Positive
            '[M-2H]-2',     # Negative
            '[M+H-H2O]+1',  # Positive
            '[M-H-H2O]-1',  # Negative
            '[M+H+MeOH]+1', # Positive
            '[M+FA+H]+1',   # Positive
            '[M-H-MeOH]-1', # Negative
            '[M+FA-H]-1'    # Negative
        ]
        for ion in REFERENCE_IONS_ALL:
            self.reference_ion_vars[ion] = tk.BooleanVar(value=True)
        
        # MS2 filter options
        self.ms2_filter_options = {
            'No MS2 filter': 'none',
            "Keep 'DDA for preferred ion' only": 'preferred_only',
            "Keep 'DDA for preferred ion' or 'DDA for non-preferred ion'": 'preferred_or_nonpreferred'
        }
        
        # Progress tracking
        self.current_step = 0
        self.cleaned_data = None
        self.start_time = None
    
    def setup_ui(self):
        """Create the complete tab interface"""
        # Main container - pack into self.frame (inherited from BaseTab)
        main_container = tk.Frame(self.frame, bg='#f0f0f0')
        main_container.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Create notebook for Metabolites and Lipids tabs
        self.cleaning_notebook = ttk.Notebook(main_container)
        self.cleaning_notebook.pack(fill='both', expand=True)
        
        # Create tabs
        metabolites_tab = tk.Frame(self.cleaning_notebook, bg='#f0f0f0')
        lipids_tab = tk.Frame(self.cleaning_notebook, bg='#f0f0f0')
        
        self.cleaning_notebook.add(metabolites_tab, text="🧪 Metabolites")
        self.cleaning_notebook.add(lipids_tab, text="🧬 Lipids")
        
        # Create interfaces
        self.create_metabolite_cleaning_interface(metabolites_tab)
        self.create_lipid_cleaning_interface(lipids_tab)

    def create_metabolite_cleaning_interface(self, parent):
        """Create metabolite data cleaning interface with two-column layout:
        Column 1 (Left): Basic/Advanced settings tabs
        Column 2 (Right): Progress section with log and controls
        """
        # Create main container with grid layout (TWO COLUMNS)
        main_container = tk.Frame(parent, bg='#f0f0f0')
        main_container.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Configure grid weights for equal column distribution
        main_container.grid_columnconfigure(0, weight=1)  # Left column (settings)
        main_container.grid_columnconfigure(1, weight=1)  # Right column (progress)
        main_container.grid_rowconfigure(0, weight=1)
        
        # ===== COLUMN 1 (LEFT): Basic/Advanced Settings Notebook =====
        settings_frame = tk.LabelFrame(main_container, 
                                      text="⚙️ Settings", 
                                      font=('Arial', 12, 'bold'),
                                      bg='#f0f0f0', 
                                      fg='#2c3e50')
        settings_frame.grid(row=0, column=0, sticky='nsew', padx=(0, 3))
        
        # Create sub-notebook for Basic and Advanced tabs INSIDE Column 1
        self.metabolite_notebook = ttk.Notebook(settings_frame)
        self.metabolite_notebook.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Basic Settings Tab
        basic_tab = tk.Frame(self.metabolite_notebook, bg='#f0f0f0')
        self.metabolite_notebook.add(basic_tab, text="📁 Basic")
        
        # Advanced Settings Tab
        advanced_tab = tk.Frame(self.metabolite_notebook, bg='#f0f0f0')
        self.metabolite_notebook.add(advanced_tab, text="⚙️ Advanced")
        
        # Create content for each tab
        self.create_metabolite_basic_tab(basic_tab)
        self.create_metabolite_advanced_tab(advanced_tab)
        
        # ===== COLUMN 2 (RIGHT): Progress Section with Scrollbar =====
        progress_container = tk.LabelFrame(main_container, 
                                          text="📊 Progress & Control", 
                                          font=('Arial', 12, 'bold'),
                                          bg='#f0f0f0', 
                                          fg='#2c3e50')
        progress_container.grid(row=0, column=1, sticky='nsew', padx=(3, 0))
        
        # Create scrollable container for progress section
        right_canvas = tk.Canvas(progress_container, bg='#f0f0f0', highlightthickness=0)
        right_scrollbar = ttk.Scrollbar(progress_container, orient='vertical', command=right_canvas.yview)
        right_scrollable_frame = tk.Frame(right_canvas, bg='#f0f0f0')
        
        right_scrollable_frame.bind('<Configure>', 
                                    lambda e: right_canvas.configure(scrollregion=right_canvas.bbox('all')))
        
        canvas_window = right_canvas.create_window((0, 0), window=right_scrollable_frame, anchor='nw')
        
        def configure_right_scroll_region(event):
            right_canvas.configure(scrollregion=right_canvas.bbox("all"))
            canvas_width = event.width
            right_canvas.itemconfig(canvas_window, width=canvas_width)
        
        right_canvas.bind('<Configure>', configure_right_scroll_region)
        right_canvas.configure(yscrollcommand=right_scrollbar.set)
        
        right_canvas.pack(side="left", fill="both", expand=True, padx=5, pady=5)
        right_scrollbar.pack(side="right", fill="y")
        
        # Bind mouse wheel for right column scrolling
        def _on_right_mousewheel(event):
            right_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        right_canvas.bind("<MouseWheel>", _on_right_mousewheel)
        right_scrollable_frame.bind("<MouseWheel>", _on_right_mousewheel)
        
        # Create progress content inside scrollable frame
        self.create_metabolite_progress_section(right_scrollable_frame)
    
    def create_metabolite_basic_tab(self, parent):
        """Create basic settings tab for metabolite cleaning - Column 1 content only"""
        # Create scrollable frame for content
        canvas = tk.Canvas(parent, bg='#f0f0f0', highlightthickness=0)
        scrollbar = ttk.Scrollbar(parent, orient='vertical', command=canvas.yview)
        scrollable_frame = tk.Frame(canvas, bg='#f0f0f0')
        
        scrollable_frame.bind('<Configure>', 
                            lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        
        canvas_window = canvas.create_window((0, 0), window=scrollable_frame, anchor='nw')
        
        def configure_scroll_region(event):
            canvas.configure(scrollregion=canvas.bbox("all"))
            canvas_width = event.width
            canvas.itemconfig(canvas_window, width=canvas_width)
        
        canvas.bind('<Configure>', configure_scroll_region)
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True, padx=5, pady=5)
        scrollbar.pack(side="right", fill="y")
        
        # Bind mouse wheel
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind("<MouseWheel>", _on_mousewheel)
        scrollable_frame.bind("<MouseWheel>", _on_mousewheel)
        
        # File selection frame
        file_frame = tk.LabelFrame(scrollable_frame, 
                                  text="📁 File Selection", 
                                  font=('Arial', 11, 'bold'),
                                  bg='#f0f0f0')
        file_frame.pack(fill='x', padx=10, pady=10)
        
        # Processing mode selection
        mode_frame = tk.LabelFrame(file_frame, 
                                  text="Processing Mode", 
                                  bg='#f0f0f0')
        mode_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Radiobutton(mode_frame, 
                      text="Combined file (split by polarity)", 
                      variable=self.cleaning_mode, 
                      value="combined",
                      bg='#f0f0f0',
                      command=self.update_cleaning_interface).pack(anchor='w', padx=5, pady=2)
        tk.Radiobutton(mode_frame, 
                      text="Separate positive/negative files", 
                      variable=self.cleaning_mode, 
                      value="separate",
                      bg='#f0f0f0',
                      command=self.update_cleaning_interface).pack(anchor='w', padx=5, pady=2)
        
        # File selection frame
        self.file_frame = tk.LabelFrame(file_frame, 
                                       text="Input Files", 
                                       bg='#f0f0f0')
        self.file_frame.pack(fill='both', expand=True, padx=10, pady=10)
            
        # Create file selection widgets (initially for combined mode)
        self.create_file_selection_widgets()
            
        # Output settings
        output_frame = tk.LabelFrame(scrollable_frame, 
                        text="💾 Output Settings", 
                        font=('Arial', 11, 'bold'),
                        bg='#f0f0f0')
        output_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(output_frame, 
            text="Select output file path:", 
            bg='#f0f0f0',
            font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=(5, 8))
            
        # Single output file path
        output_path_frame = tk.Frame(output_frame, bg='#f0f0f0')
        output_path_frame.pack(fill='x', padx=5, pady=(0, 8))
            
        self.metabolite_output_file = tk.StringVar(value=os.path.join(os.getcwd(), "cleaned_metabolites.xlsx"))
        tk.Entry(output_path_frame, 
            textvariable=self.metabolite_output_file, 
            font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
        tk.Button(output_path_frame, 
            text="Browse", 
            command=self.browse_metabolite_output_file,
            bg='#3498db', 
            fg='white', 
            font=('Arial', 8, 'bold')).pack(side='right')
        
        tk.Label(output_frame,
                text="Note: Cleaned metabolites saved with Positive, Negative, and Combined sheets",
                bg='#f0f0f0',
                font=('Arial', 7, 'italic'),
                fg='#7f8c8d').pack(anchor='w', padx=5, pady=(2, 5))
    
    def create_metabolite_progress_section(self, parent):
        """Create progress section for Column 2 (Right) - Progress log, bar, and controls"""
        # Create container frame
        container = tk.Frame(parent, bg='#f0f0f0')
        container.pack(fill='both', expand=True, padx=3, pady=3)
        
        # Button frame for action buttons
        button_frame = tk.Frame(container, bg='#f0f0f0')
        button_frame.pack(fill='x', padx=5, pady=5)
        
        # Start button - DISABLED until column verification is done
        self.start_button = tk.Button(button_frame, 
                                     text="🚀 Start Data Cleaning", 
                                     command=self.start_data_cleaning,
                                     bg='#27ae60', 
                                     fg='white', 
                                     font=('Arial', 11, 'bold'),
                                     height=1,
                                     relief='raised',
                                     borderwidth=1,
                                     padx=10, 
                                     pady=5,
                                     state='disabled')  # DISABLED until verification
        self.start_button.pack(side='left', fill='x', expand=True, padx=(0, 3))
        
        # Verify Columns button
        self.verify_button = tk.Button(button_frame, 
                                      text="🔍 Verify Columns", 
                                      command=self.verify_data_cleaning_columns,
                                      bg='#3498db', 
                                      fg='white', 
                                      font=('Arial', 11, 'bold'),
                                      height=1,
                                      relief='raised',
                                      borderwidth=1,
                                      padx=10, 
                                      pady=5)
        self.verify_button.pack(side='left', fill='x', expand=True, padx=(0, 3))
        
        # Open Output Folder button
        self.open_output_folder_button = tk.Button(button_frame,
                                                  text="📂 Open Output Folder",
                                                  command=self.open_metabolite_output_folder,
                                                  bg='#9b59b6',
                                                  fg='white',
                                                  font=('Arial', 11, 'bold'),
                                                  height=1,
                                                  relief='raised',
                                                  borderwidth=1,
                                                  padx=10,
                                                  pady=5)
        self.open_output_folder_button.pack(side='left', fill='x', expand=True)
        
        # Status label
        self.process_status = tk.Label(container,
                                      text="Ready to start data cleaning",
                                      bg='#f0f0f0',
                                      font=('Arial', 9, 'bold'),
                                      fg='#7f8c8d',
                                      wraplength=350)
        self.process_status.pack(pady=(0, 5), padx=5, fill='x')
        
        # Progress bar section
        progress_bar_frame = tk.LabelFrame(container,
                                          text="📊 Progress",
                                          font=('Arial', 10, 'bold'),
                                          bg='#f0f0f0')
        progress_bar_frame.pack(fill='x', padx=5, pady=3)
        
        # Progress bar
        self.main_progress_bar = ttk.Progressbar(progress_bar_frame,
                                                mode='determinate',
                                                length=300)
        self.main_progress_bar.pack(fill='x', padx=5, pady=5)
        
        # Progress info row
        info_frame = tk.Frame(progress_bar_frame, bg='#f0f0f0')
        info_frame.pack(fill='x', padx=5, pady=(0, 5))
        
        self.progress_percentage = tk.Label(info_frame,
                                           text="0%",
                                           bg='#f0f0f0',
                                           font=('Arial', 9, 'bold'),
                                           fg='#3498db')
        self.progress_percentage.pack(side='left')
        
        self.timer_label = tk.Label(info_frame,
                                   text="Time: 00:00",
                                   bg='#f0f0f0',
                                   font=('Arial', 8),
                                   fg='#7f8c8d')
        self.timer_label.pack(side='right')
        
        self.eta_label = tk.Label(info_frame,
                                 text="ETA: --:--",
                                 bg='#f0f0f0',
                                 font=('Arial', 8),
                                 fg='#7f8c8d')
        self.eta_label.pack(side='right', padx=(0, 15))
        
        # Progress log section
        from tkinter import scrolledtext
        log_frame = tk.LabelFrame(container,
                                 text="📝 Progress Log",
                                 font=('Arial', 10, 'bold'),
                                 bg='#f0f0f0')
        log_frame.pack(fill='both', expand=True, padx=8, pady=5)
        
        # Progress text widget
        self.progress_text = scrolledtext.ScrolledText(log_frame,
                                                       height=25,
                                                       width=40,
                                                       font=('Consolas', 8),
                                                       bg='#ffffff',
                                                       fg='#2c3e50',
                                                       wrap=tk.WORD,
                                                       relief='sunken',
                                                       borderwidth=1)
        self.progress_text.pack(fill='both', expand=True, padx=3, pady=3)
        
        # Add initial message
        self.progress_text.insert('1.0', "Waiting to start data cleaning...\n Before you begin:\n Import Excels and Verify Required columns for Data Cleaning\n")
        self.progress_text.config(state='disabled')  # Make read-only initially
        
    def create_file_selection_widgets(self):
        """Create file selection widgets based on current mode"""
        # Clear existing widgets
        for widget in self.file_frame.winfo_children():
            widget.destroy()
            
        if self.cleaning_mode.get() == "combined":
            self.create_combined_mode_widgets()
        else:
            self.create_separate_mode_widgets()
            
    def create_combined_mode_widgets(self):
        """Create widgets for combined file mode"""
        # Combined file
        tk.Label(self.file_frame, 
                text="Combined Excel file:", 
                bg='#f0f0f0').pack(anchor='w', padx=5, pady=2)
        
        combined_frame = tk.Frame(self.file_frame, bg='#f0f0f0')
        combined_frame.pack(fill='x', padx=5, pady=2)
        
        tk.Entry(combined_frame, 
                textvariable=self.combined_file, 
                font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
        tk.Button(combined_frame, 
                 text="Browse", 
                 command=lambda: self.browse_file(self.combined_file, "Combined Excel file"),
                 bg='#3498db', 
                 fg='white', 
                 font=('Arial', 8, 'bold')).pack(side='right')
        
        # mzcloud file (optional)
        tk.Label(self.file_frame, 
                text="MzCloud file (optional):", 
                bg='#f0f0f0').pack(anchor='w', padx=5, pady=(10, 2))

        mzcloud_frame = tk.Frame(self.file_frame, bg='#f0f0f0')
        mzcloud_frame.pack(fill='x', padx=5, pady=2)
        
        tk.Entry(mzcloud_frame, 
                textvariable=self.metabolika_file, 
                font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
        tk.Button(mzcloud_frame, 
                 text="Browse", 
                 command=lambda: self.browse_file(self.metabolika_file, "Metabolika Excel file"),
                 bg='#3498db', 
                 fg='white', 
                 font=('Arial', 8, 'bold')).pack(side='right')
        
    def create_separate_mode_widgets(self):
        """Create widgets for separate files mode"""
        # Positive file
        tk.Label(self.file_frame, 
                text="Positive mode Excel file:", 
                bg='#f0f0f0').pack(anchor='w', padx=5, pady=2)
        
        pos_frame = tk.Frame(self.file_frame, bg='#f0f0f0')
        pos_frame.pack(fill='x', padx=5, pady=2)
        
        tk.Entry(pos_frame, 
                textvariable=self.positive_file, 
                font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
        tk.Button(pos_frame, 
                 text="Browse", 
                 command=lambda: self.browse_file(self.positive_file, "Positive mode Excel file"),
                 bg='#3498db', 
                 fg='white', 
                 font=('Arial', 8, 'bold')).pack(side='right')
        
        # MzCloud files (optional)
        tk.Label(self.file_frame, 
                text="MzCloud Positive (optional):", 
                bg='#f0f0f0').pack(anchor='w', padx=5, pady=(10, 2))
        
        mzcloud_pos_frame = tk.Frame(self.file_frame, bg='#f0f0f0')
        mzcloud_pos_frame.pack(fill='x', padx=5, pady=2)
        
        tk.Entry(mzcloud_pos_frame, 
                textvariable=self.mzcloud_pos_file, 
                font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
        tk.Button(mzcloud_pos_frame, 
                 text="Browse", 
                 command=lambda: self.browse_file(self.mzcloud_pos_file, "MzCloud Positive Excel file"),
                 bg='#3498db', 
                 fg='white', 
                 font=('Arial', 8, 'bold')).pack(side='right')
        
        # Negative mode option (specific to separate mode) - moved here
        neg_option_frame = tk.Frame(self.file_frame, bg='#f0f0f0')
        neg_option_frame.pack(fill='x', padx=5, pady=(10, 5))
        
        tk.Checkbutton(neg_option_frame, 
                      text="Include negative mode processing", 
                      variable=self.use_negative_mode,
                      bg='#f0f0f0',
                      command=self.toggle_negative_file_widgets).pack(anchor='w')
        
        # Negative file frame (will be shown/hidden based on checkbox)
        self.negative_file_frame = tk.Frame(self.file_frame, bg='#f0f0f0')
        self.negative_file_frame.pack(fill='x', padx=5, pady=2)
        
        tk.Label(self.negative_file_frame, 
                text="Negative mode Excel file:", 
                bg='#f0f0f0').pack(anchor='w', pady=2)
        
        neg_frame = tk.Frame(self.negative_file_frame, bg='#f0f0f0')
        neg_frame.pack(fill='x', pady=2)
        
        tk.Entry(neg_frame, 
                textvariable=self.negative_file, 
                font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
        tk.Button(neg_frame, 
                 text="Browse", 
                 command=lambda: self.browse_file(self.negative_file, "Negative mode Excel file"),
                 bg='#3498db', 
                 fg='white', 
                 font=('Arial', 8, 'bold')).pack(side='right')
        
        # MzCloud Negative (moved inside negative file frame)
        tk.Label(self.negative_file_frame, 
                text="MzCloud Negative (optional):", 
                bg='#f0f0f0').pack(anchor='w', padx=5, pady=(10, 2))
        
        self.mzcloud_neg_frame = tk.Frame(self.negative_file_frame, bg='#f0f0f0')
        self.mzcloud_neg_frame.pack(fill='x', padx=5, pady=2)
        
        tk.Entry(self.mzcloud_neg_frame, 
                textvariable=self.mzcloud_neg_file, 
                font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
        tk.Button(self.mzcloud_neg_frame, 
                 text="Browse", 
                 command=lambda: self.browse_file(self.mzcloud_neg_file, "MzCloud Negative Excel file"),
                 bg='#3498db', 
                 fg='white', 
                 font=('Arial', 8, 'bold')).pack(side='right')
        
        # Initialize the toggle state
        self.toggle_negative_file_widgets()
    
    def create_metabolite_advanced_tab(self, parent):
        """Create advanced settings tab for metabolite cleaning - Single column only"""
        # Create scrollable frame for content
        canvas = tk.Canvas(parent, bg='#f0f0f0', highlightthickness=0)
        scrollbar = ttk.Scrollbar(parent, orient='vertical', command=canvas.yview)
        scrollable_frame = tk.Frame(canvas, bg='#f0f0f0')
        
        # Configure scrolling
        scrollable_frame.bind('<Configure>', 
                            lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        
        # Create window in canvas
        canvas_window = canvas.create_window((0, 0), window=scrollable_frame, anchor='nw')
        
        def configure_scroll_region(event):
            canvas.configure(scrollregion=canvas.bbox("all"))
            canvas_width = event.width
            canvas.itemconfig(canvas_window, width=canvas_width)
        
        canvas.bind('<Configure>', configure_scroll_region)
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Pack canvas and scrollbar
        canvas.pack(side="left", fill="both", expand=True, padx=5, pady=5)
        scrollbar.pack(side="right", fill="y")
        
        # Bind mouse wheel to canvas for scrolling
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind("<MouseWheel>", _on_mousewheel)
        scrollable_frame.bind("<MouseWheel>", _on_mousewheel)
        
        # 1. Filtering Thresholds Section
        threshold_frame = tk.LabelFrame(scrollable_frame, text="📊 Filtering Thresholds", 
                                       bg='#f0f0f0', font=('Arial', 11, 'bold'))
        threshold_frame.pack(fill='x', padx=10, pady=10)
        
        # PPM Threshold
        ppm_frame = tk.Frame(threshold_frame, bg='#f0f0f0')
        ppm_frame.pack(fill='x', padx=10, pady=5)
        tk.Label(ppm_frame, text="PPM Threshold:", bg='#f0f0f0', width=15, anchor='w').pack(side='left')
        self.ppm_threshold = tk.StringVar(value="10.0")
        tk.Entry(ppm_frame, textvariable=self.ppm_threshold, width=10,
                font=('Arial', 9)).pack(side='left', padx=5)
        tk.Label(ppm_frame, text="(max allowed ppm difference)", bg='#f0f0f0',
                font=('Arial', 8), fg='#7f8c8d').pack(side='left')
        
        # RT Threshold
        rt_frame = tk.Frame(threshold_frame, bg='#f0f0f0')
        rt_frame.pack(fill='x', padx=10, pady=5)
        tk.Label(rt_frame, text="RT Threshold:", bg='#f0f0f0', width=15, anchor='w').pack(side='left')
        self.rt_threshold = tk.StringVar(value="4.0")
        tk.Entry(rt_frame, textvariable=self.rt_threshold, width=10,
                font=('Arial', 9)).pack(side='left', padx=5)
        tk.Label(rt_frame, text="(RT deviation for deduplication)", bg='#f0f0f0',
                font=('Arial', 8), fg='#7f8c8d').pack(side='left')
        
        # MS2 Filter
        ms2_frame = tk.Frame(threshold_frame, bg='#f0f0f0')
        ms2_frame.pack(fill='x', padx=10, pady=5)
        tk.Label(ms2_frame, text="MS2 filter:", bg='#f0f0f0', width=15, anchor='w').pack(side='left')
        self.ms2_filter_mode = tk.StringVar(value='No MS2 filter')
        self.ms2_filter_options = {
            'No MS2 filter': 'none',
            "Keep 'DDA for preferred ion' only": 'preferred_only',
            "Keep 'DDA for preferred ion' or 'DDA for non-preferred ion'": 'preferred_or_nonpreferred'
        }
        ms2_combo = ttk.Combobox(ms2_frame, values=list(self.ms2_filter_options.keys()),
                                textvariable=tk.StringVar(value=list(self.ms2_filter_options.keys())[0]),
                                state='readonly', width=50)
        ms2_combo.pack(side='left', padx=5)
        
        def _sync_ms2_mode(event=None):
            label = ms2_combo.get()
            self.ms2_filter_mode.set(self.ms2_filter_options.get(label, 'none'))
        ms2_combo.bind('<<ComboboxSelected>>', _sync_ms2_mode)
        ms2_frame.bind('<Map>', _sync_ms2_mode)
        
        # 2. Reference Ion Filtering Section
        ion_filter_frame = tk.LabelFrame(scrollable_frame, text="⚡ Reference Ion Filtering", 
                                        bg='#f0f0f0', font=('Arial', 11, 'bold'))
        ion_filter_frame.pack(fill='x', padx=10, pady=10)
        
        # Description
        tk.Label(ion_filter_frame, 
                text="Select which reference ions to include in the analysis. All ions are checked by default.",
                bg='#f0f0f0', wraplength=600, justify='left',
                font=('Arial', 9)).pack(anchor='w', padx=10, pady=5)
        
        # Button row
        button_row = tk.Frame(ion_filter_frame, bg='#f0f0f0')
        button_row.pack(fill='x', padx=10, pady=5)
        tk.Button(button_row, text="✓ Select All", command=self.select_all_reference_ions,
                 bg='#27ae60', fg='white', font=('Arial', 9, 'bold')).pack(side='left', padx=2)
        tk.Button(button_row, text="✗ Deselect All", command=self.deselect_all_reference_ions,
                 bg='#e74c3c', fg='white', font=('Arial', 9, 'bold')).pack(side='left', padx=2)
        tk.Button(button_row, text="↻ Reset to Default", command=self.reset_reference_ions_to_default,
                 bg='#3498db', fg='white', font=('Arial', 9, 'bold')).pack(side='left', padx=2)
        
        # Create checkboxes in two columns - USE EXISTING reference_ion_vars
        checkbox_container = tk.Frame(ion_filter_frame, bg='#f0f0f0')
        checkbox_container.pack(fill='x', padx=20, pady=5)
        
        left_col = tk.Frame(checkbox_container, bg='#f0f0f0')
        left_col.pack(side='left', fill='both', expand=True)
        right_col = tk.Frame(checkbox_container, bg='#f0f0f0')
        right_col.pack(side='left', fill='both', expand=True)
        
        # Populate checkboxes from existing reference_ion_vars
        for idx, (ion, var) in enumerate(self.reference_ion_vars.items()):
            polarity = "Positive" if '+' in ion else "Negative"
            color = '#27ae60' if '+' in ion else '#e74c3c'
            
            # Alternate between columns
            parent_col = left_col if idx % 2 == 0 else right_col
            
            cb_frame = tk.Frame(parent_col, bg='#f0f0f0')
            cb_frame.pack(anchor='w', padx=5, pady=2)
            tk.Checkbutton(cb_frame, text=f"{ion}", variable=var, bg='#f0f0f0',
                          font=('Arial', 9)).pack(side='left')
            tk.Label(cb_frame, text=f"({polarity})", bg='#f0f0f0',
                    font=('Arial', 8), fg=color).pack(side='left', padx=(5, 0))
        
        # 3. Column Name Cleaning Section
        column_clean_frame = tk.LabelFrame(scrollable_frame, text="🧹 Column Name Cleaning", 
                                          bg='#f0f0f0', font=('Arial', 11, 'bold'))
        column_clean_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(column_clean_frame, 
                text="Enter text patterns to remove from numeric column names (comma-separated):",
                bg='#f0f0f0', wraplength=600, justify='left',
                font=('Arial', 9)).pack(anchor='w', padx=10, pady=5)
        
        # Text input with default patterns (use existing self.column_clean_patterns)
        pattern_entry = tk.Entry(column_clean_frame, textvariable=self.column_clean_patterns, 
                                font=('Arial', 9))
        pattern_entry.pack(fill='x', padx=10, pady=5)
        
        # Examples
        tk.Label(column_clean_frame, 
                text="Examples: 'Group Area:' → removes from 'Group Area: Sample1' → 'Sample1'",
                bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w', padx=10)
        tk.Label(column_clean_frame, 
                text="          '_pos' → removes from 'Control_pos' → 'Control'",
                bg='#f0f0f0', font=('Arial', 8, 'italic'), fg='#7f8c8d').pack(anchor='w', padx=10)
        
        # Info box
        info_frame = tk.Frame(scrollable_frame, bg='#e8f4f8', relief='solid', borderwidth=1)
        info_frame.pack(fill='x', padx=10, pady=20)
        tk.Label(info_frame, text="ℹ️ Note:", bg='#e8f4f8', 
                font=('Arial', 10, 'bold'), fg='#2980b9').pack(anchor='w', padx=10, pady=(10, 5))
        tk.Label(info_frame, 
                text="• Reference ion filtering applies to both Stage 1 and Stage 2 cleaning",
                bg='#e8f4f8', font=('Arial', 9), fg='#34495e', justify='left').pack(anchor='w', padx=20, pady=2)
        tk.Label(info_frame, 
                text="• Column name cleaning only affects numeric sample columns, not metadata",
                bg='#e8f4f8', font=('Arial', 9), fg='#34495e', justify='left').pack(anchor='w', padx=20, pady=2)
        tk.Label(info_frame, 
                text="• Settings are saved automatically and restored when you reopen the application",
                bg='#e8f4f8', font=('Arial', 9), fg='#34495e', justify='left').pack(anchor='w', padx=20, pady=(2, 10))
    
    def select_all_reference_ions(self):
        """Select all reference ion checkboxes"""
        for var in self.reference_ion_vars.values():
            var.set(True)
    
    def deselect_all_reference_ions(self):
        """Deselect all reference ion checkboxes"""
        for var in self.reference_ion_vars.values():
            var.set(False)
    
    def reset_reference_ions_to_default(self):
        """Reset all reference ions to default (all checked)"""
        for var in self.reference_ion_vars.values():
            var.set(True)
    
    def toggle_negative_file_widgets(self):
        """Toggle visibility of negative file widgets based on checkbox"""
        if hasattr(self, 'negative_file_frame'):
            if self.use_negative_mode.get():
                self.negative_file_frame.pack(fill='x', padx=5, pady=2)
                # Also show MzCloud negative if it exists
                if hasattr(self, 'mzcloud_neg_frame'):
                    self.mzcloud_neg_frame.pack(fill='x', padx=5, pady=2)
            else:
                self.negative_file_frame.pack_forget()
                # Also hide MzCloud negative if it exists
                if hasattr(self, 'mzcloud_neg_frame'):
                    self.mzcloud_neg_frame.pack_forget()
    
    def create_lipid_cleaning_interface(self, parent):
        """Create lipid data cleaning interface"""
        # Create main content frame with grid layout (TWO-COLUMN LAYOUT LIKE ID ANNOTATION)
        content_frame = tk.Frame(parent, bg='#f0f0f0')
        content_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Configure grid weights with equal distribution
        content_frame.grid_columnconfigure(0, weight=1)  # Left column
        content_frame.grid_columnconfigure(1, weight=1)  # Right column
        # Equal row weight so both columns get equal space, no minsize to avoid overflow
        content_frame.grid_rowconfigure(0, weight=1)
        
        # Left column - Settings
        left_frame = tk.LabelFrame(content_frame, 
                                  text="⚙️ Lipid Processing Settings", 
                                  font=('Arial', 12, 'bold'),
                                  bg='#f0f0f0', 
                                  fg='#2c3e50')
        left_frame.grid(row=0, column=0, sticky='nsew', padx=(0, 3))
        left_frame.grid_propagate(False)  # Prevent frame from expanding beyond grid cell
        
        # Create scrollable frame
        left_canvas = tk.Canvas(left_frame, bg='#f0f0f0', highlightthickness=0, width=400)
        left_scrollbar = ttk.Scrollbar(left_frame, orient='vertical', command=left_canvas.yview)
        left_scrollable_frame = tk.Frame(left_canvas, bg='#f0f0f0')
        
        left_scrollable_frame.bind('<Configure>', 
                                  lambda e: left_canvas.configure(scrollregion=left_canvas.bbox('all')))
        
        canvas_window = left_canvas.create_window((0, 0), window=left_scrollable_frame, anchor='nw')
        
        def configure_scroll_region(event):
            left_canvas.configure(scrollregion=left_canvas.bbox("all"))
            canvas_width = event.width
            left_canvas.itemconfig(canvas_window, width=canvas_width)
        
        left_canvas.bind('<Configure>', configure_scroll_region)
        left_canvas.configure(yscrollcommand=left_scrollbar.set)
        left_canvas.pack(side="left", fill="both", expand=True, padx=5, pady=5)
        left_scrollbar.pack(side="right", fill="y", padx=(0, 2))
        
        # Bind mouse wheel
        def _on_left_mousewheel(event):
            left_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        left_canvas.bind("<MouseWheel>", _on_left_mousewheel)
        left_scrollable_frame.bind("<MouseWheel>", _on_left_mousewheel)
        
        # File selection frame
        file_frame = tk.LabelFrame(left_scrollable_frame, 
                                  text="📁 File Selection", 
                                  bg='#f0f0f0',
                                  font=('Arial', 10, 'bold'))
        file_frame.pack(fill='x', padx=10, pady=10)
        
        # Processing mode selection
        mode_frame = tk.LabelFrame(file_frame, 
                                  text="Processing Mode", 
                                  bg='#f0f0f0')
        mode_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Radiobutton(mode_frame, 
                      text="Combined file (split by polarity)", 
                      variable=self.lipid_cleaning_mode, 
                      value="combined",
                      bg='#f0f0f0',
                      command=self.update_lipid_file_interface).pack(anchor='w', padx=5, pady=2)
        tk.Radiobutton(mode_frame, 
                      text="Separate positive/negative files", 
                      variable=self.lipid_cleaning_mode, 
                      value="separate",
                      bg='#f0f0f0',
                      command=self.update_lipid_file_interface).pack(anchor='w', padx=5, pady=2)
        
        # Dynamic file input frame
        self.lipid_file_input_frame = tk.LabelFrame(file_frame, 
                                                    text="Input Files", 
                                                    bg='#f0f0f0')
        self.lipid_file_input_frame.pack(fill='x', padx=10, pady=10)
        
        # Create initial file selection widgets (for combined mode)
        self.create_lipid_file_selection_widgets()
        
        # Dynamic File ID imports frame
        self.lipid_file_id_frame = tk.LabelFrame(file_frame,
                                                 text="File ID Mapping (Optional)",
                                                 bg='#f0f0f0')
        self.lipid_file_id_frame.pack(fill='x', padx=10, pady=10)
        
        # Create initial file ID widgets (for combined mode)
        self.create_lipid_file_id_widgets()
 
        # Output settings
        output_frame = tk.LabelFrame(left_scrollable_frame, 
                                    text="💾 Output Settings", 
                                    font=('Arial', 10, 'bold'),
                                    bg='#f0f0f0')
        output_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(output_frame, 
                text="Select output file path:", 
                bg='#f0f0f0',
                font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=(5, 8))
        
        # Single output file path
        output_path_frame = tk.Frame(output_frame, bg='#f0f0f0')
        output_path_frame.pack(fill='x', padx=5, pady=(0, 8))
        
        self.lipid_output_file = tk.StringVar(value=os.path.join(os.getcwd(), "cleaned_lipids.xlsx"))
        tk.Entry(output_path_frame, 
                textvariable=self.lipid_output_file, 
                font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
        tk.Button(output_path_frame, 
                 text="Browse", 
                 command=self.browse_lipid_output_file,
                 bg='#3498db', 
                 fg='white', 
                 font=('Arial', 8, 'bold')).pack(side='right')
        
        tk.Label(output_frame,
                text="Note: Positive, Negative, and Class tables saved as separate sheets",
                bg='#f0f0f0',
                font=('Arial', 7, 'italic'),
                fg='#7f8c8d').pack(anchor='w', padx=5, pady=(2, 5))
        
        # Adduct filtering options
        adduct_frame = tk.LabelFrame(left_scrollable_frame, 
                                    text="Adduct Ion Filtering", 
                                    bg='#f0f0f0')
        adduct_frame.pack(fill='x', padx=10, pady=10)
        
        # Available adducts list
        available_adducts = [
            'M+H', 'M-H', 'M+2H', 'M-2H', 'M+NH4', 'M+H-H2O', 'M+Na',
            'M+HCOO', 'M+K', 'M+Li', 'M+CH3COO', 'M-CH3'
        ]
        
        # Default selected adducts
        default_adducts = ['M+H', 'M-H', 'M+2H', 'M-2H', 'M+NH4', 'M+H-H2O']
        
        # Create checkboxes for each adduct
        self.lipid_adduct_vars = {}
        adduct_checkbox_frame = tk.Frame(adduct_frame, bg='#f0f0f0')
        adduct_checkbox_frame.pack(fill='x', padx=5, pady=5)
        
        # Arrange in 3 columns
        for i, adduct in enumerate(available_adducts):
            row = i // 3
            col = i % 3
            is_default = adduct in default_adducts
            var = tk.BooleanVar(value=is_default)
            self.lipid_adduct_vars[adduct] = var
            cb = tk.Checkbutton(adduct_checkbox_frame, 
                               text=adduct, 
                               variable=var,
                               bg='#f0f0f0')
            cb.grid(row=row, column=col, sticky='w', padx=5, pady=2)
        
        # Grade Threshold Section (NEW)
        grade_frame = tk.LabelFrame(left_scrollable_frame,
                                   text="📊 Grade-Based Quality Filtering",
                                   bg='#f0f0f0')
        grade_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(grade_frame,
                text="Specific grades to retain (comma-separated, e.g., A,B,C or A,B,C,D or A,B,C,P):",
                bg='#f0f0f0', wraplength=400, justify='left',
                font=('Arial', 9)).pack(anchor='w', padx=5, pady=5)
        
        self.lipid_grade_threshold = tk.StringVar(value="A,B,C")
        tk.Entry(grade_frame, textvariable=self.lipid_grade_threshold,
                font=('Arial', 9), width=20).pack(anchor='w', padx=5, pady=2)
        
        tk.Label(grade_frame,
                text="Area values will be set to 0 for grades NOT in this list (D, E, F, NP, ND, etc.)",
                bg='#f0f0f0', font=('Arial', 8, 'italic'),
                fg='#7f8c8d', wraplength=400, justify='left').pack(anchor='w', padx=5, pady=(0, 5))
        
        tk.Label(grade_frame,
                text="⚠️ Note: Grade columns must match Area columns (same count and sample names)",
                bg='#f0f0f0', font=('Arial', 8, 'bold'),
                fg='#e74c3c', wraplength=400, justify='left').pack(anchor='w', padx=5, pady=(0, 5))
        
        # Column Name Cleaning Section (NEW)
        column_clean_frame = tk.LabelFrame(left_scrollable_frame, 
                                          text="🧹 Column Name Cleaning", 
                                          bg='#f0f0f0')
        column_clean_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(column_clean_frame, 
                text="Enter text patterns to remove from sample column names (comma-separated):",
                bg='#f0f0f0', wraplength=400, justify='left',
                font=('Arial', 9)).pack(anchor='w', padx=5, pady=5)
        
        # Lipid-specific default patterns
        self.lipid_column_clean_patterns = tk.StringVar(
            value="Area[, ], [, ], _pos, _neg, .raw, _Area, Area_"
        )
        tk.Entry(column_clean_frame, textvariable=self.lipid_column_clean_patterns, 
                font=('Arial', 9)).pack(fill='x', padx=5, pady=5)
        
        # Examples
        tk.Label(column_clean_frame, 
                text="Examples: 'Area[' removes from 'Area[Sample1]', '_pos' removes from 'Control_pos'",
                bg='#f0f0f0', font=('Arial', 8, 'italic'), 
                fg='#7f8c8d').pack(anchor='w', padx=5, pady=(0, 5))
        
        # Right column - Progress and results with scrollbar
        right_frame = tk.LabelFrame(content_frame, 
                                   text="📋 Progress & Results", 
                                   font=('Arial', 12, 'bold'),
                                   bg='#f0f0f0', 
                                   fg='#2c3e50')
        right_frame.grid(row=0, column=1, sticky='nsew', padx=(3, 0))
        right_frame.grid_propagate(False)  # Prevent frame from expanding beyond grid cell
        
        # Create scrollable container for right column
        right_canvas = tk.Canvas(right_frame, bg='#f0f0f0', highlightthickness=0, width=400)
        right_scrollbar = ttk.Scrollbar(right_frame, orient='vertical', command=right_canvas.yview)
        right_scrollable_frame = tk.Frame(right_canvas, bg='#f0f0f0')
        
        right_scrollable_frame.bind('<Configure>', 
                                    lambda e: right_canvas.configure(scrollregion=right_canvas.bbox('all')))
        
        right_canvas_window = right_canvas.create_window((0, 0), window=right_scrollable_frame, anchor='nw')
        
        def configure_lipid_right_scroll_region(event):
            right_canvas.configure(scrollregion=right_canvas.bbox("all"))
            canvas_width = event.width
            right_canvas.itemconfig(right_canvas_window, width=canvas_width)
        
        right_canvas.bind('<Configure>', configure_lipid_right_scroll_region)
        right_canvas.configure(yscrollcommand=right_scrollbar.set)
        
        right_canvas.pack(side="left", fill="both", expand=True, padx=5, pady=5)
        right_scrollbar.pack(side="right", fill="y", padx=(0, 2))
        
        # Bind mouse wheel for right column scrolling
        def _on_lipid_right_mousewheel(event):
            right_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        right_canvas.bind("<MouseWheel>", _on_lipid_right_mousewheel)
        right_scrollable_frame.bind("<MouseWheel>", _on_lipid_right_mousewheel)
        
        # Start button and status AT THE TOP
        action_frame = tk.Frame(right_scrollable_frame, bg='#f0f0f0')
        action_frame.pack(fill='x', padx=10, pady=(10, 10))
        
        # Button frame for action buttons (side by side like metabolite)
        lipid_button_frame = tk.Frame(action_frame, bg='#f0f0f0')
        lipid_button_frame.pack(fill='x', padx=0, pady=(0, 5))
        
        # Start button - DISABLED until column verification is done
        self.lipid_start_button = tk.Button(lipid_button_frame, 
                                           text="🚀 Start Lipid Cleaning", 
                                           command=self.start_lipid_cleaning,
                                           bg='#27ae60', 
                                           fg='white', 
                                           font=('Arial', 11, 'bold'),
                                           height=1,
                                           relief='raised',
                                           borderwidth=1,
                                           padx=10, 
                                           pady=5,
                                           state='disabled')  # DISABLED until verification
        self.lipid_start_button.pack(side='left', fill='x', expand=True, padx=(0, 3))
        
        # Verify Columns button
        self.lipid_verify_button = tk.Button(lipid_button_frame, 
                                            text="🔍 Verify Columns", 
                                            command=self.verify_lipid_cleaning_columns,
                                            bg='#3498db', 
                                            fg='white', 
                                            font=('Arial', 11, 'bold'),
                                            height=1,
                                            relief='raised',
                                            borderwidth=1,
                                            padx=10, 
                                            pady=5)
        self.lipid_verify_button.pack(side='left', fill='x', expand=True, padx=(0, 3))
        
        # Open Output Folder button for lipids
        self.lipid_open_output_folder_button = tk.Button(lipid_button_frame,
                                                         text="📂 Open Folder",
                                                         command=self.open_lipid_output_folder,
                                                         bg='#9b59b6',
                                                         fg='white',
                                                         font=('Arial', 11, 'bold'),
                                                         height=1,
                                                         relief='raised',
                                                         borderwidth=1,
                                                         padx=10,
                                                         pady=5)
        self.lipid_open_output_folder_button.pack(side='left', fill='x', expand=True)
        
        # Status label
        self.lipid_process_status = tk.Label(action_frame,
                                            text="Ready to start lipid cleaning",
                                            bg='#f0f0f0',
                                            font=('Arial', 10),
                                            fg='#7f8c8d')
        self.lipid_process_status.pack(pady=2)
        # Class table generation status label
        self.lipid_class_status = tk.Label(action_frame,
                                           text="Class tables: (pending)",
                                           bg='#f0f0f0',
                                           font=('Arial', 9),
                                           fg='#7f8c8d')
        self.lipid_class_status.pack(pady=(0,4))
        
        # Progress bar section
        progress_bar_frame = tk.LabelFrame(right_scrollable_frame,
                                          text="📊 Progress",
                                          font=('Arial', 11, 'bold'),
                                          bg='#f0f0f0')
        progress_bar_frame.pack(fill='x', padx=10, pady=10)
        
        # Progress bar
        self.lipid_progress_bar = ttk.Progressbar(progress_bar_frame,
                                                  mode='determinate',
                                                  length=300)
        self.lipid_progress_bar.pack(fill='x', padx=10, pady=10)
        
        # Progress info row (percentage, timer, ETA)
        lipid_info_frame = tk.Frame(progress_bar_frame, bg='#f0f0f0')
        lipid_info_frame.pack(fill='x', padx=10, pady=(0, 10))
        
        self.lipid_progress_percentage = tk.Label(lipid_info_frame,
                                                  text="0%",
                                                  bg='#f0f0f0',
                                                  font=('Arial', 10, 'bold'),
                                                  fg='#3498db')
        self.lipid_progress_percentage.pack(side='left')
        
        self.lipid_timer_label = tk.Label(lipid_info_frame,
                                         text="Time: 00:00",
                                         bg='#f0f0f0',
                                         font=('Arial', 9),
                                         fg='#7f8c8d')
        self.lipid_timer_label.pack(side='right')
        
        self.lipid_eta_label = tk.Label(lipid_info_frame,
                                       text="ETA: --:--",
                                       bg='#f0f0f0',
                                       font=('Arial', 9),
                                       fg='#7f8c8d')
        self.lipid_eta_label.pack(side='right', padx=(0, 15))
        
        # Progress text area - no fixed height, expands with window
        self.lipid_progress_text = scrolledtext.ScrolledText(right_scrollable_frame, 
                                                             width=50,
                                                             font=('Courier', 9))
        self.lipid_progress_text.pack(fill='both', expand=True, padx=10, pady=10)
    
    def create_lipid_file_selection_widgets(self):
        """Create file selection widgets based on lipid cleaning mode"""
        # Clear existing widgets
        for widget in self.lipid_file_input_frame.winfo_children():
            widget.destroy()
        
        mode = self.lipid_cleaning_mode.get()
        
        if mode == "combined":
            # Combined file mode
            tk.Label(self.lipid_file_input_frame, 
                    text="LipidSearch Combined File:", 
                    bg='#f0f0f0',
                    font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=5)
            
            combined_frame = tk.Frame(self.lipid_file_input_frame, bg='#f0f0f0')
            combined_frame.pack(fill='x', padx=5, pady=2)
            
            tk.Entry(combined_frame, 
                    textvariable=self.lipid_combined_file, 
                    font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
            tk.Button(combined_frame, 
                     text="Browse", 
                     command=lambda: self.browse_file(self.lipid_combined_file, "LipidSearch Combined File"),
                     bg='#3498db', 
                     fg='white', 
                     font=('Arial', 8, 'bold')).pack(side='right')
            
        else:  # separate mode
            # Positive file
            tk.Label(self.lipid_file_input_frame, 
                    text="LipidSearch Positive Mode File (Optional):", 
                    bg='#f0f0f0',
                    font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=5)
            
            pos_frame = tk.Frame(self.lipid_file_input_frame, bg='#f0f0f0')
            pos_frame.pack(fill='x', padx=5, pady=2)
            
            tk.Entry(pos_frame, 
                    textvariable=self.lipid_positive_file, 
                    font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
            tk.Button(pos_frame, 
                     text="Browse", 
                     command=lambda: self.browse_file(self.lipid_positive_file, "LipidSearch Positive File"),
                     bg='#3498db', 
                     fg='white', 
                     font=('Arial', 8, 'bold')).pack(side='right')
            
            # Negative file
            tk.Label(self.lipid_file_input_frame, 
                    text="LipidSearch Negative Mode File (Optional):", 
                    bg='#f0f0f0',
                    font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=5)
            
            neg_frame = tk.Frame(self.lipid_file_input_frame, bg='#f0f0f0')
            neg_frame.pack(fill='x', padx=5, pady=2)
            
            tk.Entry(neg_frame, 
                    textvariable=self.lipid_negative_file, 
                    font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
            tk.Button(neg_frame, 
                     text="Browse", 
                     command=lambda: self.browse_file(self.lipid_negative_file, "LipidSearch Negative File"),
                     bg='#3498db', 
                     fg='white', 
                     font=('Arial', 8, 'bold')).pack(side='right')
            
            # Add note about optional files
            tk.Label(self.lipid_file_input_frame,
                    text="Note: You can import only positive, only negative, or both files",
                    bg='#f0f0f0',
                    font=('Arial', 8, 'italic'),
                    fg='#7f8c8d').pack(anchor='w', padx=5, pady=(5, 0))
    
    def create_lipid_file_id_widgets(self):
        """Create File ID mapping widgets based on lipid cleaning mode"""
        # Clear existing widgets
        for widget in self.lipid_file_id_frame.winfo_children():
            widget.destroy()
        
        tk.Label(self.lipid_file_id_frame,
                text="Excel files with DataID and File Name columns:",
                bg='#f0f0f0',
                font=('Arial', 9, 'bold')).pack(anchor='w', padx=5, pady=(5, 2))
        
        mode = self.lipid_cleaning_mode.get()
        
        if mode == "combined":
            # For combined mode, only show one File ID field
            tk.Label(self.lipid_file_id_frame,
                    text="File ID Mapping (Optional):",
                    bg='#f0f0f0',
                    font=('Arial', 8)).pack(anchor='w', padx=5, pady=(8, 2))
            
            file_id_frame = tk.Frame(self.lipid_file_id_frame, bg='#f0f0f0')
            file_id_frame.pack(fill='x', padx=5, pady=2)
            
            tk.Entry(file_id_frame,
                    textvariable=self.lipid_positive_sample_details_file,
                    font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
            tk.Button(file_id_frame,
                     text="Browse",
                     command=lambda: self.browse_file(self.lipid_positive_sample_details_file, "File ID Mapping"),
                     bg='#3498db',
                     fg='white',
                     font=('Arial', 8, 'bold')).pack(side='right')
            
            note_text = "Note: Will replace DataIDs in Area columns with File Names from File ID file"
        else:  # separate mode
            # For separate mode, show both Positive and Negative File ID fields
            # Positive File ID
            tk.Label(self.lipid_file_id_frame,
                    text="Positive Mode File ID (Optional):",
                    bg='#f0f0f0',
                    font=('Arial', 8)).pack(anchor='w', padx=5, pady=(8, 2))
            
            pos_file_id_frame = tk.Frame(self.lipid_file_id_frame, bg='#f0f0f0')
            pos_file_id_frame.pack(fill='x', padx=5, pady=2)
            
            tk.Entry(pos_file_id_frame,
                    textvariable=self.lipid_positive_sample_details_file,
                    font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
            tk.Button(pos_file_id_frame,
                     text="Browse",
                     command=lambda: self.browse_file(self.lipid_positive_sample_details_file, "Positive Mode File ID"),
                     bg='#3498db',
                     fg='white',
                     font=('Arial', 8, 'bold')).pack(side='right')
            
            # Negative File ID
            tk.Label(self.lipid_file_id_frame,
                    text="Negative Mode File ID (Optional):",
                    bg='#f0f0f0',
                    font=('Arial', 8)).pack(anchor='w', padx=5, pady=(8, 2))
            
            neg_file_id_frame = tk.Frame(self.lipid_file_id_frame, bg='#f0f0f0')
            neg_file_id_frame.pack(fill='x', padx=5, pady=2)
            
            tk.Entry(neg_file_id_frame,
                    textvariable=self.lipid_negative_sample_details_file,
                    font=('Arial', 9)).pack(side='left', fill='x', expand=True, padx=(0, 5))
            tk.Button(neg_file_id_frame,
                     text="Browse",
                     command=lambda: self.browse_file(self.lipid_negative_sample_details_file, "Negative Mode File ID"),
                     bg='#3498db',
                     fg='white',
                     font=('Arial', 8, 'bold')).pack(side='right')
            
            note_text = "Note: Will replace DataIDs in Area columns with File Names from respective File ID files"
        
        tk.Label(self.lipid_file_id_frame,
                text=note_text,
                bg='#f0f0f0',
                font=('Arial', 7, 'italic'),
                fg='#7f8c8d').pack(anchor='w', padx=5, pady=(5, 5))
    
    def update_lipid_file_interface(self):
        """Update file selection interface when lipid mode changes"""
        self.create_lipid_file_selection_widgets()
        self.create_lipid_file_id_widgets()
    
    def browse_lipid_output_file(self):
        """Browse for lipid output file path"""
        file_path = filedialog.asksaveasfilename(
            title="Select Output File",
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")],
            initialfile="cleaned_lipids.xlsx"
        )
        if file_path:
            self.lipid_output_file.set(file_path)
    
    def start_lipid_cleaning(self):
        """Start the lipid cleaning process"""
        try:
            # Validate inputs based on mode
            mode = self.lipid_cleaning_mode.get()
            
            if mode == "combined":
                if not self.lipid_combined_file.get():
                    messagebox.showerror("Error", "Please select a combined lipid file")
                    return
                if not os.path.exists(self.lipid_combined_file.get()):
                    messagebox.showerror("Error", "Combined file does not exist")
                    return
                # Set for backward compatibility
                self.lipid_input_file.set(self.lipid_combined_file.get())
            else:  # separate mode
                pos_file = self.lipid_positive_file.get()
                neg_file = self.lipid_negative_file.get()
                
                # At least one file must be provided
                if not pos_file and not neg_file:
                    messagebox.showerror("Error", "Please select at least one file (positive or negative mode)")
                    return
                
                # Check if provided files exist
                if pos_file and not os.path.exists(pos_file):
                    messagebox.showerror("Error", "Positive mode file does not exist")
                    return
                if neg_file and not os.path.exists(neg_file):
                    messagebox.showerror("Error", "Negative mode file does not exist")
                    return
            
            # Initialize progress tracking
            self.lipid_start_time = time.time()
            self.lipid_current_step = 0
            self.lipid_total_steps = 100
            
            # Setup progress bar for determinate mode
            self.lipid_progress_bar.config(mode='determinate', maximum=100, value=0)
            self.lipid_progress_percentage.config(text="0%")
            self.lipid_timer_label.config(text="Time: 00:00")
            self.lipid_eta_label.config(text="ETA: --:--")
            
            # Disable start button
            self.lipid_start_button.config(state='disabled', text="Processing...")
            self.lipid_process_status.config(text="Initializing lipid cleaning...", fg='#f39c12')
            
            # Clear progress text
            self.lipid_progress_text.config(state='normal')
            self.lipid_progress_text.delete('1.0', tk.END)
            self.lipid_progress_text.insert('1.0', "🚀 Starting lipid cleaning process...\n" + "="*50 + "\n")
            self.lipid_progress_text.config(state='disabled')
            
            # Start timer update
            self.update_lipid_timer()
            
            # Prepare configuration
            output_file = self.lipid_output_file.get() or os.path.join(os.getcwd(), "cleaned_lipids.xlsx")
            
            # Get selected adducts
            selected_adducts = [adduct for adduct, var in self.lipid_adduct_vars.items() if var.get()]
            
            # Get column cleaning patterns
            patterns_str = self.lipid_column_clean_patterns.get()
            column_clean_patterns = [p.strip() for p in patterns_str.split(',') if p.strip()]
            
            # Process sample details mapping if provided
            # For separate mode, load separate File ID files for positive and negative
            sample_mapping = None
            positive_sample_mapping = None
            negative_sample_mapping = None
            
            # Load positive File ID if provided
            pos_sample_file = self.lipid_positive_sample_details_file.get()
            if pos_sample_file and os.path.exists(pos_sample_file):
                try:
                    sample_df = pd.read_excel(pos_sample_file)
                    if 'DataID' in sample_df.columns and 'File Name' in sample_df.columns:
                        positive_sample_mapping = {
                            str(data_id).strip(): str(file_name).strip() 
                            for data_id, file_name in zip(sample_df['DataID'], sample_df['File Name'])
                        }
                        self.append_to_lipid_log(f"✓ Loaded POSITIVE File ID mapping with {len(positive_sample_mapping)} entries")
                        for data_id, file_name in list(positive_sample_mapping.items())[:3]:
                            self.append_to_lipid_log(f"  → {data_id} → {file_name}")
                        # For combined mode, use positive mapping as the single mapping
                        if mode == "combined":
                            sample_mapping = positive_sample_mapping
                    else:
                        self.append_to_lipid_log("⚠ Positive File ID missing 'DataID' or 'File Name' columns")
                except Exception as e:
                    self.append_to_lipid_log(f"⚠ Failed to load positive File ID: {str(e)}")
            
            # Load negative File ID if provided
            neg_sample_file = self.lipid_negative_sample_details_file.get()
            if neg_sample_file and os.path.exists(neg_sample_file):
                try:
                    sample_df = pd.read_excel(neg_sample_file)
                    if 'DataID' in sample_df.columns and 'File Name' in sample_df.columns:
                        negative_sample_mapping = {
                            str(data_id).strip(): str(file_name).strip() 
                            for data_id, file_name in zip(sample_df['DataID'], sample_df['File Name'])
                        }
                        self.append_to_lipid_log(f"✓ Loaded NEGATIVE File ID mapping with {len(negative_sample_mapping)} entries")
                        for data_id, file_name in list(negative_sample_mapping.items())[:3]:
                            self.append_to_lipid_log(f"  → {data_id} → {file_name}")
                    else:
                        self.append_to_lipid_log("⚠ Negative File ID missing 'DataID' or 'File Name' columns")
                except Exception as e:
                    self.append_to_lipid_log(f"⚠ Failed to load negative File ID: {str(e)}")
            
            # Build config based on mode
            if mode == "combined":
                config = {
                    'mode': mode,
                    'input_file': self.lipid_combined_file.get(),
                    'output_file': output_file,
                    'adduct_filter': selected_adducts if selected_adducts else None,
                    'column_clean_patterns': column_clean_patterns if column_clean_patterns else None,
                    'grade_threshold': self.lipid_grade_threshold.get() if hasattr(self, 'lipid_grade_threshold') else '<= C',
                    'column_assignments': getattr(self, 'lipid_verified_assignments', {}),
                    'sample_mapping': sample_mapping,
                    'verified_dataframes': {},
                }
                
                # Log column assignments for verification
                assignments = getattr(self, 'lipid_verified_assignments', {})
                self.append_to_lipid_log(f"\n✅ Column Assignments Verified:")
                self.append_to_lipid_log(f"  • LipidID: {assignments.get('LipidID', '(not assigned)')}")
                self.append_to_lipid_log(f"  • Class: {assignments.get('Class', '(not assigned)')} [REQUIRED FOR CLASS ANALYSIS]")
                self.append_to_lipid_log(f"  • AdductIon: {assignments.get('AdductIon', '(not assigned)')}")
                self.append_to_lipid_log(f"  • Sample columns: {len(sample_mapping)} samples\n")
                
                # Add verified dataframe if available (from verify columns)
                if hasattr(self, 'lipid_verified_dataframe') and self.lipid_verified_dataframe is not None:
                    config['verified_dataframes']['combined'] = self.lipid_verified_dataframe
                    
            else:  # separate mode - need to process as two separate files
                pos_file = self.lipid_positive_file.get()
                neg_file = self.lipid_negative_file.get()
                
                # For separate mode, we'll process positive and negative separately
                # Backend clean_lipid_data expects single input_file, so we process in sequence
                config = {
                    'mode': 'separate',
                    'positive_file': pos_file,
                    'negative_file': neg_file,
                    'output_file': output_file,
                    'adduct_filter': selected_adducts if selected_adducts else None,
                    'column_clean_patterns': column_clean_patterns if column_clean_patterns else None,
                    'grade_threshold': self.lipid_grade_threshold.get() if hasattr(self, 'lipid_grade_threshold') else '<= C',
                    'positive_assignments': getattr(self, 'lipid_verified_pos_assignments', {}),
                    'negative_assignments': getattr(self, 'lipid_verified_neg_assignments', {}),
                    'positive_sample_mapping': positive_sample_mapping,  # Separate mapping for positive
                    'negative_sample_mapping': negative_sample_mapping,  # Separate mapping for negative
                    'verified_dataframes': {},
                }
                
                # Add verified dataframes if available (from verify columns)
                if hasattr(self, 'lipid_verified_dataframe') and self.lipid_verified_dataframe is not None:
                    config['verified_dataframes']['positive'] = self.lipid_verified_dataframe
                if hasattr(self, 'lipid_verified_neg_dataframe') and self.lipid_verified_neg_dataframe is not None:
                    config['verified_dataframes']['negative'] = self.lipid_verified_neg_dataframe
            
            # Run in background thread
            def run_cleaning():
                try:
                    # Create progress callback with logging handler
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
                    
                    # Create handler that writes to lipid log widget
                    gui_handler = GUILogHandler(lambda msg: self.root.after(0, self.append_to_lipid_log, msg))
                    gui_handler.setLevel(logging.INFO)
                    gui_handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
                    
                    # Add handler to root logger
                    root_logger = logging.getLogger()
                    root_logger.addHandler(gui_handler)
                    
                    try:
                        # Create progress callback
                        def progress_callback(message):
                            # Estimate progress based on message keywords
                            if "Starting" in message:
                                progress = 10
                            elif "Loading" in message:
                                progress = 20
                            elif "Processing" in message or "Cleaning" in message:
                                progress = 40
                            elif "Analyzing" in message:
                                progress = 60
                            elif "Saving" in message or "Writing" in message:
                                progress = 80
                            elif "completed" in message or "finished" in message:
                                progress = 100
                            else:
                                # Increment current progress
                                progress = min(self.lipid_current_step + 5, 95)
                            
                            # Update progress in main thread AND log the message
                            self.root.after(0, self.update_lipid_progress_with_percentage, progress, message)
                            self.root.after(0, self.append_to_lipid_log, message)
                        
                        # Create cleaner with progress callback
                        cleaner = DataCleaner(progress_callback=progress_callback)
                        
                        # Run cleaning
                        results = cleaner.clean_lipid_data(config)
                        
                        # Update UI with results
                        self.root.after(0, lambda: self.lipid_cleaning_complete(results))
                        
                    finally:
                        # Remove the handler to avoid duplicate logging
                        root_logger.removeHandler(gui_handler)
                    
                except Exception as e:
                    error_msg = f"Error during lipid cleaning: {str(e)}\n{traceback.format_exc()}"
                    self.root.after(0, lambda: self.lipid_cleaning_error(error_msg))
            
            # Start background thread
            thread = threading.Thread(target=run_cleaning, daemon=True)
            thread.start()
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to start lipid cleaning: {str(e)}")
            self.lipid_start_button.config(state='normal', text="🚀 Start Lipid Cleaning")
            self.lipid_process_status.config(text="Ready to start lipid cleaning", fg='#7f8c8d')
    
    def lipid_cleaning_complete(self, results):
        """Handle lipid cleaning completion"""
        self.lipid_start_button.config(state='normal', text="🚀 Start Lipid Cleaning")
        
        # Update progress bar to 100%
        self.lipid_progress_bar.config(value=100)
        self.lipid_progress_percentage.config(text="100%")
        
        # Reset timer and ETA
        self.lipid_start_time = None  # Stop timer updates
        self.lipid_eta_label.config(text="ETA: 00:00")  # Show 00:00 when complete
        
        if results['success']:
            self.lipid_process_status.config(text="Lipid cleaning completed successfully!", fg='#27ae60')
            
            # Log completion to progress text
            self.append_to_lipid_log("="*50)
            self.append_to_lipid_log("✅ LIPID CLEANING COMPLETED SUCCESSFULLY!")
            self.append_to_lipid_log("="*50)
            
            # Show summary
            stats = results.get('stats', {})
            summary = "Lipid Cleaning Summary:\n"
            
            # Separate mode shows positive and negative counts directly
            pos_count = stats.get('positive_after_cleanup', 0)
            neg_count = stats.get('negative_after_cleanup', 0)
            
            summary += f"  Positive lipids: {pos_count}\n"
            summary += f"  Negative lipids: {neg_count}\n"
            summary += f"  Total lipids: {pos_count + neg_count}\n"
            
            # Only show combined mode stats if they exist
            if 'original_count' in stats:
                summary += f"\n  Original lipids (combined): {stats.get('original_count', 0)}\n"
                summary += f"  After filtering: {stats.get('after_filtering', 0)}\n"
            
            summary += f"\nOutput file: {results['files_created'][0]}"
            
            # Check if the workflow logged a critical error about missing Area columns
            # The metabolite_data_cleaner.py logs "CRITICAL ERROR" and "No 'Area[' columns found!" if Area columns are missing
            area_error_found = results.get('area_columns_missing', False)
            
            if area_error_found:
                self.append_to_lipid_log("\n" + "="*50)
                self.append_to_lipid_log("⚠️⚠️⚠️ CRITICAL WARNING ⚠️⚠️⚠️")
                self.append_to_lipid_log("="*50)
                self.append_to_lipid_log("\nNo 'Area[' columns were detected!")
                self.append_to_lipid_log("\nThis indicates the Area column assignment may be missing.")
                self.append_to_lipid_log("\n📋 RECOMMENDATION:")
                self.append_to_lipid_log("   Click the 'Verify Columns' button and explicitly assign")
                self.append_to_lipid_log("   the Area columns before running lipid cleaning again.")
                self.append_to_lipid_log("\nThank you!")
                self.append_to_lipid_log("="*50 + "\n")
            
            # Store lipid DataFrames in shared memory_store similar to metabolite workflow
            if hasattr(self, 'memory_store') and isinstance(self.memory_store, dict):
                self.memory_store.clear()
            else:
                # Fallback: re-link to shared DataManager store
                self.memory_store = self.data_manager.memory_store
            # Combined sheet (for annotation) preferred for loading
            combined_df = results.get('combined_lipids')
            if isinstance(combined_df, pd.DataFrame) and not combined_df.empty:
                self.memory_store['lipid_combined_df'] = combined_df
                print(f"🧠 Stored combined lipid dataframe: {len(combined_df)} rows, cols: {list(combined_df.columns)}")
            pos_df = results.get('pos_lipids')
            if isinstance(pos_df, pd.DataFrame) and not pos_df.empty:
                self.memory_store['lipid_pos_df'] = pos_df
                print(f"🧠 Stored positive lipid dataframe: {len(pos_df)} rows")
            neg_df = results.get('neg_lipids')
            if isinstance(neg_df, pd.DataFrame) and not neg_df.empty:
                self.memory_store['lipid_neg_df'] = neg_df
                print(f"🧠 Stored negative lipid dataframe: {len(neg_df)} rows")
            if results.get('pos_class_df') is not None:
                self.memory_store['lipid_pos_class_df'] = results.get('pos_class_df')
            if results.get('neg_class_df') is not None:
                self.memory_store['lipid_neg_class_df'] = results.get('neg_class_df')

            # Update class status label
            pos_cls = results.get('pos_class_df')
            neg_cls = results.get('neg_class_df')
            if hasattr(self, 'lipid_class_status'):
                if (isinstance(pos_cls, pd.DataFrame) and not pos_cls.empty) or (isinstance(neg_cls, pd.DataFrame) and not neg_cls.empty):
                    parts = []
                    if isinstance(pos_cls, pd.DataFrame) and not pos_cls.empty:
                        parts.append(f"+{len(pos_cls)}")
                    if isinstance(neg_cls, pd.DataFrame) and not neg_cls.empty:
                        parts.append(f"-{len(neg_cls)}")
                    self.lipid_class_status.config(text=f"Class tables: generated ({', '.join(parts)})", fg='#27ae60')
                else:
                    self.lipid_class_status.config(text="Class tables: none generated", fg='#e67e22')

            # Auto-load combined lipids for ID annotation: store path and notify ID Annotation tab
            if results.get('files_created'):
                lipid_output_file_path = results['files_created'][0]
                self.cleaned_lipid_excel_path = lipid_output_file_path
                
                # Store separate paths for positive and negative sheets
                # These are in the same Excel file but different sheets
                self.cleaned_lipid_positive_path = lipid_output_file_path  # Sheet: Positive_Lipids
                self.cleaned_lipid_negative_path = lipid_output_file_path  # Sheet: Negative_Lipids
                
                # Note: Do NOT update lipid_input_file field to avoid re-loading into Data Cleaning tab
                # The cleaned data is available in memory_store for ID Annotation
                print(f"✅ Cleaned lipid file created: {lipid_output_file_path}")
                
                # Persist cleaned lipid path in DataManager for other tabs
                try:
                    if hasattr(self.data_manager, 'set_cleaned_lipid_excel_path'):
                        self.data_manager.set_cleaned_lipid_excel_path(self.cleaned_lipid_excel_path)
                    else:
                        self.data_manager.cleaned_lipid_excel_path = self.cleaned_lipid_excel_path
                except Exception:
                    self.data_manager.cleaned_lipid_excel_path = self.cleaned_lipid_excel_path

                # Notify ID Annotation tab that lipid data is ready
                self.notify_data_ready("🔬 ID Annotation", "lipid_combined_df")
                
                # Update the lipid auto-file displays in ID annotation tab (Lipids sub-tab)
                id_annotation_tab = self.get_tab_by_name("🔬 ID Annotation")
                if id_annotation_tab:
                    try:
                        # Store the cleaned lipid Excel path on the ID Annotation tab
                        setattr(id_annotation_tab, 'cleaned_lipid_excel_path', lipid_output_file_path)

                        # Check which file(s) were actually processed (could be only positive or only negative)
                        has_positive = results.get('pos_lipids') is not None and not results['pos_lipids'].empty
                        has_negative = results.get('neg_lipids') is not None and not results['neg_lipids'].empty
                        
                        # Update single lipid auto-display field
                        if hasattr(id_annotation_tab, 'lipid_auto_file_display'):
                            id_annotation_tab.lipid_auto_file_display.config(state='normal')
                            id_annotation_tab.lipid_auto_file_display.delete(0, tk.END)
                            id_annotation_tab.lipid_auto_file_display.insert(0, lipid_output_file_path)
                            id_annotation_tab.lipid_auto_file_display.config(state='readonly')
                            print(f"✅ Updated lipid auto-display: {lipid_output_file_path}")
                    except Exception as e:
                        print(f"Warning: could not update lipid auto-file display: {e}")

            # Build detailed summary message
            detailed_summary = ""
            if isinstance(combined_df, pd.DataFrame) and not combined_df.empty:
                detailed_summary += f"🧠 Stored combined lipid dataframe: {len(combined_df)} rows, cols: {list(combined_df.columns)}\n"
            if isinstance(pos_df, pd.DataFrame) and not pos_df.empty:
                detailed_summary += f"🧠 Stored positive lipid dataframe: {len(pos_df)} rows\n"
            if isinstance(neg_df, pd.DataFrame) and not neg_df.empty:
                detailed_summary += f"🧠 Stored negative lipid dataframe: {len(neg_df)} rows\n"
            if results.get('files_created'):
                detailed_summary += f"✅ Cleaned lipid file created: {results['files_created'][0]}\n"
            
            # ===== CLASS SHEET VERIFICATION =====
            # If class sheets were generated, verify their column structure
            pos_cls = results.get('pos_class_df')
            neg_cls = results.get('neg_class_df')
            
            if (isinstance(pos_cls, pd.DataFrame) and not pos_cls.empty) or (isinstance(neg_cls, pd.DataFrame) and not neg_cls.empty):
                # Show class sheet verification results
                class_verify_msg = "🏷️ CLASS SHEET VERIFICATION:\n\n"
                
                # Helper to count sample columns even after column cleaning (Area[ ] may be stripped)
                def _count_sample_cols(df_cls: pd.DataFrame) -> int:
                    feature_cols = {
                        'LipidID','Class','Class_name','LipidGroup','Charge','CalcMz','BaseRt','SubClass',
                        'AdductIon','IonFormula','MolStructure','ObsMz','ObsRt','Polarity'
                    }
                    area_like = [c for c in df_cls.columns if isinstance(c, str) and c.startswith('Area[') and c.endswith(']')]
                    if area_like:
                        return len(area_like)
                    numeric_cols = [c for c in df_cls.columns if pd.api.types.is_numeric_dtype(df_cls[c])]
                    fallback = [c for c in numeric_cols if c not in feature_cols and not (isinstance(c, str) and c.startswith('Grade['))]
                    return len(fallback)

                # Check positive class
                if isinstance(pos_cls, pd.DataFrame) and not pos_cls.empty:
                    has_class_col = 'Class' in pos_cls.columns
                    has_class_name = 'Class_name' in pos_cls.columns
                    sample_count = _count_sample_cols(pos_cls)
                    
                    class_verify_msg += f"✓ Positive_Lipid_Class: {len(pos_cls)} rows\n"
                    class_verify_msg += f"  {'✓' if has_class_col else '❌'} Class column: {'Found' if has_class_col else 'MISSING - REQUIRED'}\n"
                    class_verify_msg += f"  {'✓' if has_class_name else '✗'} Class_name column: {len(pos_cls.columns)} total columns\n"
                    class_verify_msg += f"  {'✓' if sample_count > 0 else '❌'} Area/sample columns: {sample_count} column(s)\n\n"
                
                # Check negative class
                if isinstance(neg_cls, pd.DataFrame) and not neg_cls.empty:
                    has_class_col = 'Class' in neg_cls.columns
                    has_class_name = 'Class_name' in neg_cls.columns
                    sample_count = _count_sample_cols(neg_cls)
                    
                    class_verify_msg += f"✓ Negative_Lipid_Class: {len(neg_cls)} rows\n"
                    class_verify_msg += f"  {'✓' if has_class_col else '❌'} Class column: {'Found' if has_class_col else 'MISSING - REQUIRED'}\n"
                    class_verify_msg += f"  {'✓' if has_class_name else '✗'} Class_name column: {len(neg_cls.columns)} total columns\n"
                    class_verify_msg += f"  {'✓' if sample_count > 0 else '❌'} Area/sample columns: {sample_count} column(s)\n\n"
                
                class_verify_msg += "✅ Class sheets are ready for class-level statistics analysis"
                
                # Store verified class dataframes for later use
                if isinstance(pos_cls, pd.DataFrame) and not pos_cls.empty:
                    self.lipid_verified_class_pos_df = pos_cls
                if isinstance(neg_cls, pd.DataFrame) and not neg_cls.empty:
                    self.lipid_verified_class_neg_df = neg_cls
                
                # Show verification info
                messagebox.showinfo("✅ Lipid Cleaning Complete - Class Sheets Verified", 
                                  detailed_summary + "\n\n" + class_verify_msg)
            else:
                # No class sheets generated - still a success but no class analysis possible
                messagebox.showinfo("Success", detailed_summary)
            
            # Auto-load to statistics tab in lipid mode
            stats_tab = self.get_tab_by_name("📊 Statistics")
            if stats_tab and hasattr(stats_tab, 'load_lipid_data_from_memory'):
                try:
                    stats_tab.load_lipid_data_from_memory()
                    # Message already shown above
                except Exception as e:
                    print(f"Warning: Could not auto-load to statistics tab: {e}")
        else:
            self.lipid_process_status.config(text="Lipid cleaning failed", fg='#e74c3c')
            if hasattr(self, 'lipid_class_status'):
                self.lipid_class_status.config(text="Class tables: not generated", fg='#e74c3c')
            messagebox.showerror("Error", f"Lipid cleaning failed: {results['message']}")
    
    def lipid_cleaning_error(self, error_msg):
        """Handle lipid cleaning error"""
        self.lipid_start_button.config(state='normal', text="🚀 Start Lipid Cleaning")
        self.lipid_process_status.config(text="Error during processing", fg='#e74c3c')
        if hasattr(self, 'lipid_class_status'):
            self.lipid_class_status.config(text="Class tables: error", fg='#e74c3c')
        self.lipid_progress_text.insert(tk.END, f"\n❌ Error: {error_msg}\n")
        self.lipid_progress_text.see(tk.END)
        messagebox.showerror("Error", f"Lipid cleaning failed:\n{error_msg}")
        
    def append_to_lipid_log(self, message):
        """Append message to lipid progress log text widget"""
        try:
            # Ensure the text widget is in a writable state
            current_state = self.lipid_progress_text.cget('state')
            self.lipid_progress_text.config(state='normal')
            
            # Insert the message
            self.lipid_progress_text.insert(tk.END, f"{message}\n")
            
            # Scroll to the end
            self.lipid_progress_text.see(tk.END)
            
            # Restore the disabled state
            self.lipid_progress_text.config(state='disabled')
            
        except Exception as e:
            print(f"[APPEND_TO_LIPID_LOG] ERROR: {e}")
    
    def update_lipid_progress_with_percentage(self, percentage, message):
        """Update lipid progress with percentage, timer, and message"""
        # Update current step for tracking
        self.lipid_current_step = percentage
        
        self.lipid_progress_bar.config(value=percentage)
        self.lipid_progress_percentage.config(text=f"{percentage}%")
        
        # Update process status with current message (shortened for status)
        status_msg = message[:50] + "..." if len(message) > 50 else message
        self.lipid_process_status.config(text=status_msg, fg='#3498db')
        
        self.root.update_idletasks()
    
    def update_lipid_timer(self):
        """Update the lipid timer display - Timer starts on button click and continues until 100% or error"""
        if getattr(self, '_shutting_down', False):
            return
        
        # Stop timer if lipid_start_time is None (cleaned up after analysis completed)
        if not hasattr(self, 'lipid_start_time') or self.lipid_start_time is None:
            return
        
        # Update elapsed time
        elapsed = time.time() - self.lipid_start_time
        minutes = int(elapsed // 60)
        seconds = int(elapsed % 60)
        self.lipid_timer_label.config(text=f"Time: {minutes:02d}:{seconds:02d}")
        
        # Calculate ETA only if progress is between 1-99%
        current_percentage = self.lipid_current_step
        if 0 < current_percentage < 100:
            # Estimate total time based on rate of progress
            estimated_total = elapsed / (current_percentage / 100.0)
            remaining = estimated_total - elapsed
            eta_minutes = int(remaining // 60)
            eta_seconds = int(remaining % 60)
            self.lipid_eta_label.config(text=f"ETA: {eta_minutes:02d}:{eta_seconds:02d}")
        elif current_percentage == 0:
            # No progress yet, show --:--
            self.lipid_eta_label.config(text="ETA: --:--")
        # If at 100%, ETA will be set to 00:00 in lipid_cleaning_complete
        
        # Schedule next update ONLY if timer is still running
        if hasattr(self, 'lipid_start_time') and self.lipid_start_time is not None:
            self.root.after(1000, self.update_lipid_timer)
        
    def update_cleaning_interface(self):
        """Update the interface when cleaning mode changes"""
        self.create_file_selection_widgets()
        
    def browse_file(self, file_var, title):
        """Browse for a file"""
        try:
            filename = filedialog.askopenfilename(
                title=f"Select {title}",
                filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
            )
            if filename:
                file_var.set(filename)
                
                # Force GUI update
                # Also mirror into common read-only auto-file displays if present
                try:
                    # lipid tab manual vs auto display
                    if hasattr(self, 'lipid_input_file') and file_var is getattr(self, 'lipid_input_file') and hasattr(self, 'lipid_auto_file_display'):
                        try:
                            self.lipid_auto_file_display.config(state='normal')
                            self.lipid_auto_file_display.delete(0, tk.END)
                            self.lipid_auto_file_display.insert(0, filename)
                            self.lipid_auto_file_display.config(state='readonly')
                        except Exception:
                            pass
                    # generic metabolite auto display
                    if hasattr(self, 'id_input_file') and file_var is getattr(self, 'id_input_file') and hasattr(self, 'id_auto_file_display'):
                        try:
                            self.id_auto_file_display.config(state='normal')
                            self.id_auto_file_display.delete(0, tk.END)
                            self.id_auto_file_display.insert(0, filename)
                            self.id_auto_file_display.config(state='readonly')
                        except Exception:
                            pass
                except Exception:
                    pass

                self.root.update_idletasks()
        except Exception as e:
            messagebox.showerror("Browse Error", f"Failed to open file dialog:\n{str(e)}")
    
    def browse_save_file(self, file_var, title):
        """Browse for a save file location"""
        try:
            filename = filedialog.asksaveasfilename(
                title=f"Save {title}",
                defaultextension=".xlsx",
                filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")]
            )
            if filename:
                file_var.set(filename)
                # Force GUI update
                self.winfo_toplevel().update_idletasks()
        except Exception as e:
            messagebox.showerror("Browse Error", f"Failed to open save dialog:\n{str(e)}")
            
    def browse_metabolite_output_file(self):
        """Browse for metabolite output file path"""
        file_path = filedialog.asksaveasfilename(
            title="Select Output File",
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")],
            initialfile="cleaned_metabolites.xlsx"
        )
        if file_path:
            self.metabolite_output_file.set(file_path)
    
    @staticmethod
    def _remove_readonly_attribute(filepath):
        """Remove read-only attribute from a file (Windows-specific fix)"""
        try:
            if os.path.exists(filepath):
                # Remove read-only flag on Windows
                os.chmod(filepath, stat.S_IWRITE | stat.S_IREAD)
        except Exception as e:
            # Non-critical, log but don't fail
            print(f"Warning: Could not remove read-only attribute from {filepath}: {e}")
    
    def update_timer(self):
        """Update the timer display - Timer starts on button click and continues until 100% or error"""
        if getattr(self, '_shutting_down', False):
            return
        
        # Stop timer if start_time is None (cleaned up after analysis completed)
        if self.timer_manager.start_time is None:
            return
        
        # Update elapsed time
        elapsed = time.time() - self.timer_manager.start_time
        minutes = int(elapsed // 60)
        seconds = int(elapsed % 60)
        self.timer_label.config(text=f"Time: {minutes:02d}:{seconds:02d}")
        
        # Calculate ETA only if progress is between 1-99%
        # ETA = (elapsed / progress%) - elapsed
        current_percentage = self.current_step
        if 0 < current_percentage < 100:
            # Estimate total time based on rate of progress
            estimated_total = elapsed / (current_percentage / 100.0)
            remaining = estimated_total - elapsed
            eta_minutes = int(remaining // 60)
            eta_seconds = int(remaining % 60)
            self.eta_label.config(text=f"ETA: {eta_minutes:02d}:{eta_seconds:02d}")
        elif current_percentage == 0:
            # No progress yet, show --:--
            self.eta_label.config(text="ETA: --:--")
        # If at 100%, ETA will be set to 00:00 in handle_cleaning_results
        
        # Schedule next update ONLY if timer is still running
        if self.timer_manager.start_time is not None:
            self.root.after(1000, self.update_timer)
    
    def update_progress_with_percentage(self, percentage, message):
        """Update progress with percentage, timer, and message"""
        # Update current step for tracking
        self.current_step = percentage
        
        # Update timer manager
        self.timer_manager.update_progress(percentage)
        
        self.main_progress_bar.config(value=percentage)
        self.progress_percentage.config(text=f"{percentage}%")
        
        # Update process status with current message (shortened for status)
        status_msg = message[:50] + "..." if len(message) > 50 else message
        self.process_status.config(text=status_msg, fg='#3498db')
        
        self.root.update_idletasks()
    
    def append_to_log(self, message):
        """Append message to progress log text widget"""
        try:
            # Ensure the text widget is in a writable state
            current_state = self.progress_text.cget('state')
            self.progress_text.config(state='normal')
            
            # Insert the message
            self.progress_text.insert(tk.END, f"{message}\n")
            
            # Scroll to the end
            self.progress_text.see(tk.END)
            
            # Restore the disabled state
            self.progress_text.config(state='disabled')
            
        except Exception as e:
            print(f"[APPEND_TO_LOG] ERROR: {e}")
            
    def update_progress(self, message):
        """Update progress text area and status"""
        self.append_to_log(message)
        
        # Update process status with current message (shortened for status)
        status_msg = message[:50] + "..." if len(message) > 50 else message
        self.process_status.config(text=status_msg, fg='#3498db')
        
        self.root.update_idletasks()
    
    def _validate_required_columns(self):
        """Validate that required columns (Group Area and Feature ID) are identified.
        
        Returns:
            bool: True if validation passes, False otherwise
        """
        mode = self.cleaning_mode.get()
        
        # Check if verify columns was used
        if mode == "combined":
            assignments = getattr(self, 'verified_assignments', {})
            sample_cols = getattr(self, 'verified_sample_cols', [])
        else:
            assignments = getattr(self, 'verified_pos_assignments', {})
            sample_cols = getattr(self, 'verified_pos_sample_cols', [])
        
        # Debug: Print assignments to console
        print("\n[DEBUG] Validation - Current assignments:")
        print(f"  Mode: {mode}")
        print(f"  Assignments keys: {list(assignments.keys())}")
        for key, value in assignments.items():
            if isinstance(value, list):
                print(f"  {key}: {len(value)} items")
            else:
                print(f"  {key}: {value}")
        print(f"  Sample columns: {len(sample_cols)} columns")
        
        # Check for Feature ID (required)
        feature_id = assignments.get('Feature ID')
        
        # Check for Group Area columns (required)
        # Group Area columns are stored in sample_cols from the dialog result
        has_group_area = isinstance(sample_cols, list) and len(sample_cols) > 0
        
        # Check for Formula (warning only)
        formula = assignments.get('Formula')
        has_formula = bool(formula)
        
        missing_required = []
        if not feature_id:
            missing_required.append("Feature ID")
            print(f"  ✗ Feature ID missing")
        else:
            print(f"  ✓ Feature ID: {feature_id}")
            
        if not has_group_area:
            missing_required.append("Group Area")
            print(f"  ✗ Group Area: No sample columns detected")
        else:
            print(f"  ✓ Group Area: {len(sample_cols)} sample columns detected")
        
        # If required columns are missing, show error
        if missing_required:
            error_msg = (
                f"⚠️ Required columns not identified: {', '.join(missing_required)}\n\n"
                "Please use the 'Verify Columns' button to assign these columns.\n\n"
                "Required columns:\n"
                "  • Feature ID: Main identifier for metabolites\n"
                "  • Group Area: Sample data columns\n"
            )
            if not has_formula:
                error_msg += "\nNote: 'Formula' column is recommended but not required."
            
            messagebox.showerror("Missing Required Columns", error_msg)
            return False
        
        # If Formula is missing, show warning but allow to continue
        if not has_formula:
            warning_msg = (
                "⚠️ Warning: 'Formula' column not identified.\n\n"
                "While not required, the Formula column is recommended for better analysis.\n\n"
                "Do you want to continue without it?"
            )
            if not messagebox.askyesno("Formula Column Missing", warning_msg):
                return False
        
        print("\n  ✅ All required columns validated successfully!\n")
        return True
    
    def _log_column_verification_status(self, pos_assignments, neg_assignments=None):
        """Log column verification status to console.
        
        Args:
            pos_assignments: Assignments for positive mode (or combined mode)
            neg_assignments: Optional assignments for negative mode
        """
        print("\n" + "="*60)
        print("COLUMN VERIFICATION STATUS")
        print("="*60)
        
        # Check positive/combined assignments
        feature_id = pos_assignments.get('Feature ID')
        group_area = pos_assignments.get('Group Area')
        formula = pos_assignments.get('Formula')
        
        print("\n✓ Positive/Combined Mode:")
        if feature_id:
            print(f"  ✅ Feature ID: '{feature_id}' - IDENTIFIED")
        else:
            print("  ❌ Feature ID: NOT IDENTIFIED - User must assign in Verify Columns")
        
        if isinstance(group_area, list) and len(group_area) > 0:
            print(f"  ✅ Group Area: {len(group_area)} columns - IDENTIFIED")
            print(f"      Columns: {', '.join(group_area[:3])}{'...' if len(group_area) > 3 else ''}")
        elif isinstance(group_area, str) and group_area:
            print(f"  ✅ Group Area: '{group_area}' pattern - IDENTIFIED")
        else:
            print("  ❌ Group Area: NOT IDENTIFIED - User must assign in Verify Columns")
        
        if formula:
            print(f"  ✅ Formula: '{formula}' - IDENTIFIED")
        else:
            print("  ⚠️  Formula: NOT IDENTIFIED - Recommended but not required")
        
        # Check negative assignments if in separate mode
        if neg_assignments:
            print("\n✓ Negative Mode:")
            neg_feature_id = neg_assignments.get('Feature ID')
            neg_group_area = neg_assignments.get('Group Area')
            neg_formula = neg_assignments.get('Formula')
            
            if neg_feature_id:
                print(f"  ✅ Feature ID: '{neg_feature_id}' - IDENTIFIED")
            else:
                print("  ❌ Feature ID: NOT IDENTIFIED - User must assign in Verify Columns")
            
            if isinstance(neg_group_area, list) and len(neg_group_area) > 0:
                print(f"  ✅ Group Area: {len(neg_group_area)} columns - IDENTIFIED")
            elif isinstance(neg_group_area, str) and neg_group_area:
                print(f"  ✅ Group Area: '{neg_group_area}' pattern - IDENTIFIED")
            else:
                print("  ❌ Group Area: NOT IDENTIFIED - User must assign in Verify Columns")
            
            if neg_formula:
                print(f"  ✅ Formula: '{neg_formula}' - IDENTIFIED")
            else:
                print("  ⚠️  Formula: NOT IDENTIFIED - Recommended but not required")
        
        # Summary
        print("\n" + "-"*60)
        has_all_required = bool(feature_id) and (
            (isinstance(group_area, list) and len(group_area) > 0) or 
            (isinstance(group_area, str) and group_area)
        )
        
        if has_all_required:
            print("✅ ALL REQUIRED COLUMNS IDENTIFIED - Ready to start cleaning")
        else:
            print("❌ MISSING REQUIRED COLUMNS - Use Verify Columns to assign them")
        
        print("="*60 + "\n")
    
    def open_metabolite_output_folder(self):
        """Open the folder containing metabolite output files"""
        import subprocess
        import platform
        
        # First check if we have a stored cleaned file path from processing
        output_file = getattr(self, 'last_cleaned_file_path', None)
        if not output_file:
            # Fall back to metabolite_output_file entry
            output_file = self.metabolite_output_file.get()
        if not output_file:
            messagebox.showinfo("No Output File", "No output file path set. Please run data cleaning first or set an output path.")
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
    
    def open_lipid_output_folder(self):
        """Open the folder containing lipid output files"""
        import subprocess
        import platform
        
        # Use the stored cleaned lipid excel path
        output_file = getattr(self, 'cleaned_lipid_excel_path', None)
        if not output_file:
            messagebox.showinfo("No Output File", "No output file path set. Please run lipid cleaning first or set an output path.")
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
    
    def verify_data_cleaning_columns(self):
        """Verify and assign columns for data cleaning using unified dialog"""
        import threading
        
        def _load_and_verify():
            """Background thread worker for loading files and showing dialogs"""
            try:
                mode = self.cleaning_mode.get()
                
                if mode == "combined":
                    # Combined mode - verify single file
                    file_path = self.combined_file.get()
                    
                    if not file_path or not os.path.exists(file_path):
                        self.root.after(0, lambda: messagebox.showerror("Error", "Please select a combined file first"))
                        return
                    
                    # Update status
                    self.root.after(0, lambda: self.process_status.config(text="Opening combined file..."))
                    self.root.after(0, lambda: self.progress_text.config(state='normal'))
                    self.root.after(0, lambda: self.progress_text.insert('end', "\n📂 Opening combined file...\n"))
                    self.root.after(0, lambda: self.progress_text.config(state='disabled'))
                    self.root.after(0, lambda: self.progress_text.see('end'))
                    
                    # Load the dataframe
                    if file_path.endswith(('.xlsx', '.xls')):
                        df = pd.read_excel(file_path)
                    elif file_path.endswith('.csv'):
                        df = pd.read_csv(file_path)
                    else:
                        self.root.after(0, lambda: messagebox.showerror("Error", "Unsupported file format"))
                        return
                    
                    # Update status
                    self.root.after(0, lambda: self.process_status.config(text="Analyzing columns..."))
                    self.root.after(0, lambda: self.progress_text.config(state='normal'))
                    self.root.after(0, lambda: self.progress_text.insert('end', "✓ Analyzing columns...\n"))
                    self.root.after(0, lambda: self.progress_text.config(state='disabled'))
                    self.root.after(0, lambda: self.progress_text.see('end'))
                    
                    # Show column assignment dialog
                    result = show_column_assignment_dialog(
                        parent=self.root,
                        df=df,
                        tab_type='data_cleaning',
                        auto_calculate=False,
                    )
                    
                    if result:
                        self.verified_assignments = result['assignments']
                        self.verified_sample_cols = result.get('sample_cols', [])  # Store sample columns
                        if result.get('dataframe') is not None:
                            self.verified_dataframe = result['dataframe']
                        
                        # Log column status to console
                        self._log_column_verification_status(result['assignments'])
                        
                        # Build detailed verification summary (required + essentials)
                        def _build_dc_summary(assignments, sample_cols, df, mode_label="Combined"):
                            lines = [f"{mode_label} Verification:"]
                            fid = assignments.get('Feature ID')
                            lines.append(f"  {'✓' if fid else '✗'} Feature ID: {fid if fid else 'Not assigned'}")
                            ga_ok = isinstance(sample_cols, list) and len(sample_cols) > 0
                            ga_detail = f"{len(sample_cols)} column(s)" if ga_ok else "No sample columns"
                            lines.append(f"  {'✓' if ga_ok else '✗'} Group Area: {ga_detail}")
                            essentials = TAB_REQUIREMENTS.get('data_cleaning', {}).get('essential', [])
                            if essentials:
                                lines.append("  Essentials:")
                                for e in essentials:
                                    col = assignments.get(e)
                                    ok = bool(col and col in df.columns)
                                    detail = col if col else 'Not assigned'
                                    lines.append(f"    {'✓' if ok else '✗'} {e}: {detail}")
                            return "\n".join(lines)
                        summary_text = _build_dc_summary(result['assignments'], self.verified_sample_cols, df, mode_label="Combined")
                        self.root.after(0, lambda: messagebox.showinfo(
                            "✅ Column Verification Complete",
                            summary_text
                        ))
                        self.root.after(0, lambda: (
                            self.progress_text.config(state='normal'),
                            self.progress_text.insert('end', f"✅ Column verification completed\n"),
                            self.progress_text.config(state='disabled'),
                            self.progress_text.see('end'),
                            self.process_status.config(text="✅ Column verification complete"),
                            self.start_button.config(state='normal')  # ENABLE start button
                        ))
                else:
                    # Separate mode - verify positive and optionally negative
                    pos_file = self.positive_file.get()
                    
                    if not pos_file or not os.path.exists(pos_file):
                        self.root.after(0, lambda: messagebox.showerror("Error", "Please select a positive file first"))
                        return
                    
                    # Update status - opening positive file
                    self.root.after(0, lambda: self.process_status.config(text="Opening Positive file..."))
                    self.root.after(0, lambda: self.progress_text.config(state='normal'))
                    self.root.after(0, lambda: self.progress_text.insert('end', "\n📂 Opening Positive file...\n"))
                    self.root.after(0, lambda: self.progress_text.config(state='disabled'))
                    self.root.after(0, lambda: self.progress_text.see('end'))
                    
                    if pos_file.endswith(('.xlsx', '.xls')):
                        pos_df = pd.read_excel(pos_file)
                    elif pos_file.endswith('.csv'):
                        pos_df = pd.read_csv(pos_file)
                    else:
                        self.root.after(0, lambda: messagebox.showerror("Error", "Unsupported file format for positive file"))
                        return
                    
                    # Update status - analyzing positive
                    self.root.after(0, lambda: self.process_status.config(text="Analyzing Positive file columns..."))
                    self.root.after(0, lambda: self.progress_text.config(state='normal'))
                    self.root.after(0, lambda: self.progress_text.insert('end', "✓ Analyzing Positive file columns...\n"))
                    self.root.after(0, lambda: self.progress_text.config(state='disabled'))
                    self.root.after(0, lambda: self.progress_text.see('end'))
                    
                    # Show dialog for positive
                    pos_result = show_column_assignment_dialog(
                        parent=self.root,
                        df=pos_df,
                        tab_type='data_cleaning',
                        auto_calculate=False,
                        dialog_title="Positive File – Data Cleaning Column Assignment",
                    )
                    
                    if not pos_result:
                        self.root.after(0, lambda: self.process_status.config(text="❌ Verification cancelled"))
                        return  # User cancelled
                    
                    # Store positive assignments AND dataframe
                    self.verified_pos_assignments = pos_result['assignments']
                    self.verified_pos_sample_cols = pos_result.get('sample_cols', [])  # Store sample columns
                    self.verified_assignments = pos_result['assignments']  # Keep for backward compatibility
                    if pos_result.get('dataframe') is not None:
                        self.verified_dataframe = pos_result['dataframe']  # Store positive dataframe
                    
                    # Verify NEGATIVE file if present
                    if self.use_negative_mode.get():
                        neg_file = self.negative_file.get()
                        
                        if neg_file and os.path.exists(neg_file):
                            # Update status - opening negative file
                            self.root.after(0, lambda: self.process_status.config(text="Opening Negative file..."))
                            self.root.after(0, lambda: self.progress_text.config(state='normal'))
                            self.root.after(0, lambda: self.progress_text.insert('end', "\n📂 Opening Negative file...\n"))
                            self.root.after(0, lambda: self.progress_text.config(state='disabled'))
                            self.root.after(0, lambda: self.progress_text.see('end'))
                            
                            if neg_file.endswith(('.xlsx', '.xls')):
                                neg_df = pd.read_excel(neg_file)
                            elif neg_file.endswith('.csv'):
                                neg_df = pd.read_csv(neg_file)
                            else:
                                self.root.after(0, lambda: messagebox.showerror("Error", "Unsupported file format for negative file"))
                                return
                            
                            # Update status - analyzing negative
                            self.root.after(0, lambda: self.process_status.config(text="Analyzing Negative file columns..."))
                            self.root.after(0, lambda: self.progress_text.config(state='normal'))
                            self.root.after(0, lambda: self.progress_text.insert('end', "✓ Analyzing Negative file columns...\n"))
                            self.root.after(0, lambda: self.progress_text.config(state='disabled'))
                            self.root.after(0, lambda: self.progress_text.see('end'))
                            
                            # Show dialog for negative
                            neg_result = show_column_assignment_dialog(
                                parent=self.root,
                                df=neg_df,
                                tab_type='data_cleaning',
                                auto_calculate=False,
                                dialog_title="Negative File – Data Cleaning Column Assignment",
                            )
                            
                            if neg_result:
                                self.verified_neg_assignments = neg_result['assignments']
                                self.verified_neg_sample_cols = neg_result.get('sample_cols', [])  # Store sample columns
                                if neg_result.get('dataframe') is not None:
                                    self.verified_neg_dataframe = neg_result['dataframe']  # Store negative dataframe
                    
                    # Show detailed confirmation per polarity
                    def _build_dc_summary(assignments, sample_cols, df, mode_label):
                        lines = [f"{mode_label} Verification:"]
                        fid = assignments.get('Feature ID')
                        lines.append(f"  {'✓' if fid else '✗'} Feature ID: {fid if fid else 'Not assigned'}")
                        ga_ok = isinstance(sample_cols, list) and len(sample_cols) > 0
                        ga_detail = f"{len(sample_cols)} column(s)" if ga_ok else "No sample columns"
                        lines.append(f"  {'✓' if ga_ok else '✗'} Group Area: {ga_detail}")
                        essentials = TAB_REQUIREMENTS.get('data_cleaning', {}).get('essential', [])
                        if essentials:
                            lines.append("  Essentials:")
                            for e in essentials:
                                col = assignments.get(e)
                                ok = bool(col and col in df.columns)
                                detail = col if col else 'Not assigned'
                                lines.append(f"    {'✓' if ok else '✗'} {e}: {detail}")
                        return "\n".join(lines)
                    pos_summary = _build_dc_summary(self.verified_pos_assignments, self.verified_pos_sample_cols, pos_df, "Positive")
                    neg_summary = None
                    if hasattr(self, 'verified_neg_assignments'):
                        neg_summary = _build_dc_summary(self.verified_neg_assignments, getattr(self, 'verified_neg_sample_cols', []), neg_df, "Negative")
                    full_msg = pos_summary + ("\n\n" + neg_summary if neg_summary else "")
                    self.root.after(0, lambda: messagebox.showinfo(
                        "✅ Column Verification Complete",
                        full_msg
                    ))
                    
                    self.root.after(0, lambda: (
                        self.progress_text.config(state='normal'),
                        self.progress_text.insert('end', f"✅ Column verification completed for separate mode\n"),
                        self.progress_text.config(state='disabled'),
                        self.progress_text.see('end'),
                        self.process_status.config(text="✅ Column verification complete"),
                        self.start_button.config(state='normal')  # ENABLE start button
                    ))
                    
            except Exception as e:
                logger.error(f"Error verifying columns: {e}")
                self.root.after(0, lambda: messagebox.showerror("Error", f"Failed to verify columns: {str(e)}"))
                self.root.after(0, lambda: self.process_status.config(text="❌ Error during column verification"))
        
        # Start verification in background thread
        self.process_status.config(text="Initializing verification...")
        verification_thread = threading.Thread(target=_load_and_verify, daemon=True)
        verification_thread.start()
    
    def verify_lipid_cleaning_columns(self):
        """Verify and assign columns for lipid data cleaning using unified dialog"""
        import threading
        
        def _load_and_verify():
            """Background thread worker for loading file and showing dialog"""
            try:
                # Get file path based on mode
                mode = self.lipid_cleaning_mode.get()
                
                if mode == "combined":
                    file_path = self.lipid_combined_file.get()
                    if not file_path or not os.path.exists(file_path):
                        self.root.after(0, lambda: messagebox.showerror("Error", "Please select a lipid file first"))
                        return
                    # Update status
                    self.root.after(0, lambda: self.lipid_process_status.config(text="Opening lipid file..."))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='normal'))
                    self.root.after(0, lambda: self.lipid_progress_text.insert('end', "\n📂 Opening lipid file...\n"))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='disabled'))
                    self.root.after(0, lambda: self.lipid_progress_text.see('end'))
                    
                    # Load the dataframe
                    if file_path.endswith(('.xlsx', '.xls')):
                        df = pd.read_excel(file_path)
                    elif file_path.endswith('.csv'):
                        df = pd.read_csv(file_path)
                    else:
                        self.root.after(0, lambda: messagebox.showerror("Error", "Unsupported file format"))
                        return
                    
                    # Update status
                    self.root.after(0, lambda: self.lipid_process_status.config(text="Analyzing columns..."))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='normal'))
                    self.root.after(0, lambda: self.lipid_progress_text.insert('end', "✓ Analyzing columns...\n"))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='disabled'))
                    self.root.after(0, lambda: self.lipid_progress_text.see('end'))
                    
                    # Show column assignment dialog for lipid cleaning (use lipid_cleaning tab type)
                    result = show_column_assignment_dialog(
                        parent=self.root,
                        df=df,
                        tab_type='lipid_cleaning',  # Use lipid-specific requirements
                        auto_calculate=False,  # Lipid cleaning doesn't need calculated columns
                    )
                    
                    if result:
                        # Store the assignments and modified dataframe for lipid cleaning
                        self.lipid_verified_assignments = result['assignments']
                        self.lipid_verified_sample_cols = result.get('sample_cols', [])  # Store Area columns
                        if result.get('dataframe') is not None:
                            self.lipid_verified_dataframe = result['dataframe']
                        
                        # Show detailed verification summary for lipids
                        def _build_lipid_summary(assignments, sample_cols, df, mode_label="Combined"):
                            lines = [f"{mode_label} Verification:"]
                            lid = assignments.get('LipidID')
                            lines.append(f"  {'✓' if lid else '✗'} LipidID: {lid if lid else 'Not assigned'}")
                            # Class column (CRITICAL for lipid class analysis)
                            class_col = assignments.get('Class')
                            class_ok = bool(class_col and class_col in df.columns)
                            lines.append(f"  {'✓' if class_ok else '❌'} Class (REQUIRED for class analysis): {class_col if class_col else 'Not assigned'}")

                            # Area columns: trust user assignments first, then fall back to detected sample_cols
                            assigned_area = assignments.get('Area')
                            if isinstance(assigned_area, list):
                                assigned_area_cols = assigned_area
                            elif assigned_area:
                                assigned_area_cols = [assigned_area]
                            else:
                                assigned_area_cols = []

                            # Grade columns: prefer assignments, otherwise infer by pattern
                            assigned_grade = assignments.get('Grade')
                            if isinstance(assigned_grade, list):
                                assigned_grade_cols = assigned_grade
                            elif assigned_grade:
                                assigned_grade_cols = [assigned_grade]
                            else:
                                assigned_grade_cols = [col for col in df.columns if isinstance(col, str) and col.startswith('Grade[') and col.endswith(']')]

                            # If user provided assignments use them; otherwise use sample_cols fallback
                            area_columns_display = assigned_area_cols if assigned_area_cols else (sample_cols if isinstance(sample_cols, list) else [])
                            area_ok = len(area_columns_display) > 0
                            area_detail = f"{len(area_columns_display)} Area column(s)" if area_ok else "No Area columns"
                            lines.append(f"  {'✓' if area_ok else '✗'} Area columns: {area_detail}")

                            if assigned_grade_cols:
                                lines.append(f"  ✓ Grade columns: {len(assigned_grade_cols)} column(s) (for quality filtering)")

                            # Recommended essentials for lipids
                            essentials = ['AdductIon', 'ObsMz', 'ObsRt']
                            lines.append("  Recommended:")
                            for e in essentials:
                                col = assignments.get(e)
                                ok = bool(col and col in df.columns)
                                detail = col if col else 'Not assigned'
                                lines.append(f"    {'✓' if ok else '✗'} {e}: {detail}")
                            return "\n".join(lines)
                        summary_text = _build_lipid_summary(result['assignments'], self.lipid_verified_sample_cols, df, mode_label="Combined")
                        self.root.after(0, lambda: messagebox.showinfo(
                            "✅ Column Verification Complete",
                            summary_text
                        ))
                        self.root.after(0, lambda: (
                            self.lipid_progress_text.config(state='normal'),
                            self.lipid_progress_text.insert('end', f"✅ Column verification completed\n"),
                            self.lipid_progress_text.config(state='disabled'),
                            self.lipid_progress_text.see('end'),
                            self.lipid_process_status.config(text="✅ Column verification complete"),
                            self.lipid_start_button.config(state='normal')  # ENABLE start button
                        ))
                else:  # separate mode
                    pos_file = self.lipid_positive_file.get()
                    neg_file = self.lipid_negative_file.get()
                    
                    if not pos_file or not os.path.exists(pos_file):
                        self.root.after(0, lambda: messagebox.showerror("Error", "Please select a positive lipid file first"))
                        return
                    
                    # Update status - opening positive file
                    self.root.after(0, lambda: self.lipid_process_status.config(text="Opening Positive lipid file..."))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='normal'))
                    self.root.after(0, lambda: self.lipid_progress_text.insert('end', "\n📂 Opening Positive lipid file...\n"))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='disabled'))
                    self.root.after(0, lambda: self.lipid_progress_text.see('end'))
                    
                    if pos_file.endswith(('.xlsx', '.xls')):
                        pos_df = pd.read_excel(pos_file)
                    elif pos_file.endswith('.csv'):
                        pos_df = pd.read_csv(pos_file)
                    else:
                        self.root.after(0, lambda: messagebox.showerror("Error", "Unsupported file format for positive file"))
                        return
                    
                    # Update status - analyzing positive
                    self.root.after(0, lambda: self.lipid_process_status.config(text="Analyzing Positive file columns..."))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='normal'))
                    self.root.after(0, lambda: self.lipid_progress_text.insert('end', "✓ Analyzing Positive file columns...\n"))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='disabled'))
                    self.root.after(0, lambda: self.lipid_progress_text.see('end'))
                    
                    # Show dialog for positive
                    pos_result = show_column_assignment_dialog(
                        parent=self.root,
                        df=pos_df,
                        tab_type='lipid_cleaning',
                        auto_calculate=False,
                        dialog_title="Positive Lipid File – Column Assignment",
                    )
                    
                    if not pos_result:
                        self.root.after(0, lambda: self.lipid_process_status.config(text="❌ Verification cancelled"))
                        return  # User cancelled
                    
                    # Store positive assignments AND dataframe
                    self.lipid_verified_pos_assignments = pos_result['assignments']
                    self.lipid_verified_pos_sample_cols = pos_result.get('sample_cols', [])
                    self.lipid_verified_assignments = pos_result['assignments']  # Keep for backward compatibility
                    if pos_result.get('dataframe') is not None:
                        self.lipid_verified_dataframe = pos_result['dataframe']
                    
                    # Verify NEGATIVE file
                    if not neg_file or not os.path.exists(neg_file):
                        self.root.after(0, lambda: messagebox.showerror("Error", "Please select a negative lipid file"))
                        return
                    
                    # Update status - opening negative file
                    self.root.after(0, lambda: self.lipid_process_status.config(text="Opening Negative lipid file..."))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='normal'))
                    self.root.after(0, lambda: self.lipid_progress_text.insert('end', "\n📂 Opening Negative lipid file...\n"))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='disabled'))
                    self.root.after(0, lambda: self.lipid_progress_text.see('end'))
                    
                    if neg_file.endswith(('.xlsx', '.xls')):
                        neg_df = pd.read_excel(neg_file)
                    elif neg_file.endswith('.csv'):
                        neg_df = pd.read_csv(neg_file)
                    else:
                        self.root.after(0, lambda: messagebox.showerror("Error", "Unsupported file format for negative file"))
                        return
                    
                    # Update status - analyzing negative
                    self.root.after(0, lambda: self.lipid_process_status.config(text="Analyzing Negative file columns..."))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='normal'))
                    self.root.after(0, lambda: self.lipid_progress_text.insert('end', "✓ Analyzing Negative file columns...\n"))
                    self.root.after(0, lambda: self.lipid_progress_text.config(state='disabled'))
                    self.root.after(0, lambda: self.lipid_progress_text.see('end'))
                    
                    # Show dialog for negative
                    neg_result = show_column_assignment_dialog(
                        parent=self.root,
                        df=neg_df,
                        tab_type='lipid_cleaning',
                        auto_calculate=False,
                        dialog_title="Negative Lipid File – Column Assignment",
                    )
                    
                    if neg_result:
                        self.lipid_verified_neg_assignments = neg_result['assignments']
                        self.lipid_verified_neg_sample_cols = neg_result.get('sample_cols', [])
                        if neg_result.get('dataframe') is not None:
                            self.lipid_verified_neg_dataframe = neg_result['dataframe']
                    
                    # Show detailed confirmation per polarity
                    def _build_lipid_summary(assignments, sample_cols, df, mode_label):
                        lines = [f"{mode_label} Verification:"]
                        lid = assignments.get('LipidID')
                        lines.append(f"  {'✓' if lid else '✗'} LipidID: {lid if lid else 'Not assigned'}")
                        # Class column (CRITICAL for lipid class analysis)
                        class_col = assignments.get('Class')
                        class_ok = bool(class_col and class_col in df.columns)
                        lines.append(f"  {'✓' if class_ok else '❌'} Class (REQUIRED for class analysis): {class_col if class_col else 'Not assigned'}")
                        # Area columns
                        area_ok = isinstance(sample_cols, list) and len(sample_cols) > 0
                        area_detail = f"{len(sample_cols)} Area column(s)" if area_ok else "No Area columns"
                        lines.append(f"  {'✓' if area_ok else '✗'} Area[ ] columns: {area_detail}")
                        # Recommended essentials
                        essentials = ['AdductIon', 'ObsMz', 'ObsRt']
                        lines.append("  Recommended:")
                        for e in essentials:
                            col = assignments.get(e)
                            ok = bool(col and col in df.columns)
                            detail = col if col else 'Not assigned'
                            lines.append(f"    {'✓' if ok else '✗'} {e}: {detail}")
                        return "\n".join(lines)
                    pos_summary = _build_lipid_summary(self.lipid_verified_pos_assignments, self.lipid_verified_pos_sample_cols, pos_df, "Positive")
                    neg_summary = None
                    if hasattr(self, 'lipid_verified_neg_assignments'):
                        neg_summary = _build_lipid_summary(self.lipid_verified_neg_assignments, getattr(self, 'lipid_verified_neg_sample_cols', []), neg_df, "Negative")
                    full_msg = pos_summary + ("\n\n" + neg_summary if neg_summary else "")
                    self.root.after(0, lambda: messagebox.showinfo(
                        "✅ Column Verification Complete",
                        full_msg
                    ))
                    
                    self.root.after(0, lambda: (
                        self.lipid_progress_text.config(state='normal'),
                        self.lipid_progress_text.insert('end', f"✅ Column verification completed for separate mode\n"),
                        self.lipid_progress_text.config(state='disabled'),
                        self.lipid_progress_text.see('end'),
                        self.lipid_process_status.config(text="✅ Column verification complete"),
                        self.lipid_start_button.config(state='normal')  # ENABLE start button
                    ))
                    
            except Exception as e:
                logger.error(f"Error verifying lipid columns: {e}")
                self.root.after(0, lambda: messagebox.showerror("Error", f"Failed to verify columns: {str(e)}"))
                self.root.after(0, lambda: self.lipid_process_status.config(text="❌ Error during column verification"))
        
        # Start verification in background thread
        self.lipid_process_status.config(text="Initializing verification...")
        verification_thread = threading.Thread(target=_load_and_verify, daemon=True)
        verification_thread.start()
   
    def start_data_cleaning(self):
        """Start the data cleaning process"""
        # Validate inputs
        if not self.validate_cleaning_inputs():
            return
        
        # Initialize progress tracking using TimerManager
        self.timer_manager.start(total_steps=100)
        
        # Setup progress bar for determinate mode
        self.main_progress_bar.config(mode='determinate', maximum=100, value=0)
        self.progress_percentage.config(text="0%")
        self.timer_label.config(text="Time: 00:00")
        self.eta_label.config(text="ETA: --:--")
        
        # Clear and prepare progress log
        self.progress_text.config(state='normal')
        self.progress_text.delete('1.0', tk.END)
        self.progress_text.insert('1.0', "🚀 Starting data cleaning process...\n" + "="*50 + "\n")
        
        # Check if verify columns was used and add warning if not
        verify_columns_used = hasattr(self, 'verified_assignments') and self.verified_assignments
        if not verify_columns_used:
            warning_msg = "⚠️  WARNING: Column verification not performed. Verifying that all required columns are correctly assigned allows confidence in analysis.\n\n"
            self.progress_text.insert('end', warning_msg)
        
        self.progress_text.config(state='disabled')
        
        # Disable start button and update status
        self.start_button.config(state='disabled', text="🔄 Processing...")
        self.process_status.config(text="Initializing data cleaning...", fg='#f39c12')
        
        # Start timer update
        self.update_timer()
        
        # Prepare configuration
        config = self.prepare_cleaning_config()
        
        # Start cleaning in a separate thread
        self.cleaning_thread = threading.Thread(target=self.run_data_cleaning, args=(config,))
        self.cleaning_thread.daemon = True
        self.cleaning_thread.start()
        
    def validate_cleaning_inputs(self):
        """Validate inputs for data cleaning"""
        if self.cleaning_mode.get() == "combined":
            if not self.combined_file.get():
                messagebox.showerror("Error", "Please select a combined Excel file")
                return False
        else:
            if not self.positive_file.get():
                messagebox.showerror("Error", "Please select a positive mode Excel file")
                return False
            if self.use_negative_mode.get() and not self.negative_file.get():
                messagebox.showerror("Error", "Please select a negative mode Excel file")
                return False
                
        if not self.metabolite_output_file.get():
            messagebox.showerror("Error", "Please select an output file path")
            return False
        
        # Validate enhanced cleaning parameters if enabled
        if self.use_enhanced_cleaning.get():
            try:
                ppm_val = float(self.ppm_threshold.get())
                if ppm_val <= 0:
                    messagebox.showerror("Error", "PPM threshold must be greater than 0")
                    return False
            except ValueError:
                messagebox.showerror("Error", "PPM threshold must be a valid number")
                return False
            
            try:
                rt_val = float(self.rt_threshold.get())
                if rt_val <= 0:
                    messagebox.showerror("Error", "RT threshold must be greater than 0")
                    return False
            except ValueError:
                messagebox.showerror("Error", "RT threshold must be a valid number")
                return False
            
        return True
    
    def prepare_cleaning_config(self):
        """Prepare configuration for data cleaning"""
        # Determine column assignments based on mode
        if self.cleaning_mode.get() == "combined":
            column_assignments = getattr(self, 'verified_assignments', {})
            pos_assignments = {}
            neg_assignments = {}
        else:
            # Separate mode - use specific pos/neg assignments
            pos_assignments = getattr(self, 'verified_pos_assignments', {})
            neg_assignments = getattr(self, 'verified_neg_assignments', {})
            column_assignments = pos_assignments  # Default to pos for backward compatibility
        
        # Get output file path and extract directory
        output_file_path = self.metabolite_output_file.get()
        output_dir = os.path.dirname(output_file_path) if output_file_path else os.getcwd()
        
        config = {
            'cleaning_mode': self.cleaning_mode.get(),
            'use_negative_mode': self.use_negative_mode.get(),
            'files': {},
            'output': {
                'directory': output_dir,
                'prefix': 'cleaned_metabolites'
            },
            'column_assignments': column_assignments,
            'pos_column_assignments': pos_assignments,
            'neg_column_assignments': neg_assignments,
            'verified_dataframes': {},  # Store pre-loaded dataframes from verification step
            'verify_columns_used': False,  # Track if verify columns was used
        }
        
        # Get selected reference ions (NEW)
        selected_ions = [ion for ion, var in self.reference_ion_vars.items() if var.get()]
        
        # Get column cleaning patterns (NEW)
        patterns_str = self.column_clean_patterns.get()
        column_clean_patterns = [p.strip() for p in patterns_str.split(',') if p.strip()]

        config['stage1_options'] = {}
        
        # Add enhanced cleaning parameters if enabled
        if self.use_enhanced_cleaning.get():
            config['enhanced_cleaning'] = {
                'enabled': True,
                'ppm_threshold': float(self.ppm_threshold.get()),
                'rt_threshold': float(self.rt_threshold.get()),
                'ms2_filter_mode': self.ms2_filter_mode.get(),
                'selected_ions': selected_ions if selected_ions else None,  # NEW
                'column_clean_patterns': column_clean_patterns if column_clean_patterns else None  # NEW
            }
        else:
            config['enhanced_cleaning'] = {
                'enabled': False,
                'selected_ions': selected_ions if selected_ions else None,  # NEW
                'column_clean_patterns': column_clean_patterns if column_clean_patterns else None  # NEW
            }
        
        if self.cleaning_mode.get() == "combined":
            config['files']['combined_file'] = self.combined_file.get()
            if self.metabolika_file.get():
                config['files']['combined_metabolika_file'] = self.metabolika_file.get()
            # Extract filename from full path without extension
            filename = os.path.basename(output_file_path)
            if isinstance(filename, str) and filename.lower().endswith('.xlsx'):
                filename = filename[:-5]
            config['output']['combined_filename'] = filename
            
            # If verify columns was used, pass the pre-loaded dataframe (avoids re-loading)
            if hasattr(self, 'verified_dataframe') and self.verified_dataframe is not None:
                config['verified_dataframes']['combined'] = self.verified_dataframe
                config['verify_columns_used'] = True
        else:
            config['files']['pos_raw_file'] = self.positive_file.get()
            if self.use_negative_mode.get() and self.negative_file.get():
                config['files']['neg_raw_file'] = self.negative_file.get()
            if self.mzcloud_pos_file.get():
                config['files']['pos_metabolika_file'] = self.mzcloud_pos_file.get()
            if self.mzcloud_neg_file.get():
                config['files']['neg_metabolika_file'] = self.mzcloud_neg_file.get()

            # Extract filename from full path without extension
            combined_filename = os.path.basename(output_file_path)
            if isinstance(combined_filename, str) and combined_filename.lower().endswith('.xlsx'):
                combined_filename = combined_filename[:-5]
            config['output']['combined_filename'] = combined_filename
            
            # If verify columns was used, pass the pre-loaded dataframes (avoids re-loading)
            if hasattr(self, 'verified_dataframe') and self.verified_dataframe is not None:
                config['verified_dataframes']['positive'] = self.verified_dataframe
                config['verify_columns_used'] = True
            if hasattr(self, 'verified_neg_dataframe') and self.verified_neg_dataframe is not None:
                config['verified_dataframes']['negative'] = self.verified_neg_dataframe
                config['verify_columns_used'] = True
                
        return config
        
    def run_data_cleaning(self, config):
        """Run data cleaning in background thread - WITH PROGRESS LOGGING"""
        # Create a logging handler to capture logging.info() messages
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
        
        # Create handler that writes to our log widget
        gui_handler = GUILogHandler(lambda msg: self.root.after(0, self.append_to_log, msg))
        gui_handler.setLevel(logging.INFO)
        gui_handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
        
        # Add handler to root logger
        root_logger = logging.getLogger()
        root_logger.addHandler(gui_handler)
        
        try:
            # Create custom progress callback that updates percentage AND logs messages
            def progress_callback(message):
                # Simulate progress steps based on message content
                if "Starting" in message:
                    progress = 10
                elif "Loading" in message:
                    progress = 20
                elif "Processing" in message or "Cleaning" in message:
                    progress = 40
                elif "Analyzing" in message:
                    progress = 60
                elif "Saving" in message or "Writing" in message:
                    progress = 80
                elif "completed" in message or "finished" in message:
                    progress = 100
                else:
                    # Increment current progress
                    progress = min(self.current_step + 5, 95)
                
                # Update progress in main thread AND log the message
                self.root.after(0, self.update_progress_with_percentage, progress, message)
                self.root.after(0, self.append_to_log, message)
            
            # Create data cleaner with progress callback
            cleaner = DataCleaner(progress_callback=progress_callback)
            
            # Run cleaning - use enhanced cleaning if enabled
            if config.get('enhanced_cleaning', {}).get('enabled', False):
                results = self.run_enhanced_cleaning(cleaner, config, progress_callback)
            else:
                results = cleaner.clean_data(config)
            
            # Update UI with results
            self.root.after(0, self.handle_cleaning_results, results)
            
        except Exception as e:
            error_msg = f"Error during data cleaning: {str(e)}\n{traceback.format_exc()}"
            self.root.after(0, self.handle_cleaning_error, error_msg)
        finally:
            # Remove the handler to avoid duplicate logging
            root_logger.removeHandler(gui_handler)
    
    def run_enhanced_cleaning(self, cleaner, config, progress_callback):
        """Run enhanced cleaning with new parameters - Two stage process"""
        try:
            # Initialize skipped steps tracking
            self._last_pos_skipped_steps = []
            self._last_neg_skipped_steps = []
            self._last_combined_skipped_steps = []
            
            enhanced_config = config['enhanced_cleaning']
            
            # STAGE 1: Run original cleaning first
            progress_callback("🔄 Stage 1: Running original data cleaning...")
            original_results = cleaner.clean_data(config)
            
            if not original_results['success']:
                return original_results
            
            # STAGE 2: Apply enhanced cleaning to the already-cleaned data (in-memory preferred)
            progress_callback("🔬 Stage 2: Applying enhanced cleaning to cleaned data...")

            results = {
                'success': False,
                'files_created': original_results['files_created'].copy(),
                'pos_result': None,
                'neg_result': None,
                'combined_result': None,
                'message': '',
                'error': None
            }

            output_dir = config['output']['directory']

            # Helper: derive default pos/neg filenames if not provided (combined mode)
            def _derive_pos_neg_basenames():
                base = config['output'].get('combined_filename', 'cleaned_metabolites_combined')
                if 'combined' in base.lower():
                    pos_base = base.lower().replace('combined', 'positive')
                    neg_base = base.lower().replace('combined', 'negative')
                else:
                    pos_base = f"{base}_positive"
                    neg_base = f"{base}_negative"
                # Preserve original case from base prefix for aesthetics
                prefix = config['output'].get('combined_filename', 'cleaned_metabolites_combined')
                pos_name = prefix.replace('combined', 'positive') if 'combined' in prefix else f"{prefix}_positive"
                neg_name = prefix.replace('combined', 'negative') if 'combined' in prefix else f"{prefix}_negative"
                return pos_name, neg_name

            # 1) Build in-memory pos/neg DataFrames to enhance
            pos_src_df, neg_src_df = None, None
            if config['cleaning_mode'] == 'separate':
                pos_src_df = original_results.get('pos_result')
                neg_src_df = original_results.get('neg_result') if config.get('use_negative_mode', False) else None
            else:
                # combined mode: split by Reference Ion into pos/neg
                combined_df = original_results.get('combined_result')
                if combined_df is None:
                    # Fallback: attempt to load from created file (Cleaned_Data sheet)
                    try:
                        combined_files = [f for f in original_results['files_created'] if f.lower().endswith('.xlsx')]
                        if combined_files:
                            combined_df = pd.read_excel(combined_files[0], sheet_name='Cleaned_Data')
                    except Exception:
                        combined_df = None
                if combined_df is not None and 'Reference Ion' in combined_df.columns:
                    tmp = combined_df.copy()
                    tmp['Reference Ion'] = tmp['Reference Ion'].astype(str).str.strip().apply(cleaner._normalize_ion)
                    pos_src_df = tmp[tmp['Reference Ion'].isin(getattr(cleaner, 'positive_ion_order', []))].copy()
                    neg_src_df = tmp[tmp['Reference Ion'].isin(getattr(cleaner, 'negative_ion_order', []))].copy()

            # 2) Enhance Positive (in-memory only; will be written inside Enhanced Combined)
            enhanced_pos_output = None
            pos_skipped_steps = []
            if pos_src_df is not None and len(pos_src_df) > 0:
                progress_callback("⚙️ Enhancing positive polarity data in-memory")
                pos_result = cleaner.enhanced_stage2_cleaning(
                    pos_src_df,
                    ppm_threshold=enhanced_config['ppm_threshold'],
                    rt_threshold=enhanced_config['rt_threshold'],
                    ion_mode="positive",
                    ms2_filter_mode=enhanced_config.get('ms2_filter_mode', 'none'),
                    progress_callback=progress_callback,
                    selected_ions=enhanced_config.get('selected_ions'),  # NEW
                    column_clean_patterns=enhanced_config.get('column_clean_patterns')  # NEW
                )
                pos_enhanced_df = pos_result['data']
                pos_skipped_steps = pos_result.get('skipped_steps', [])
                self._last_pos_skipped_steps = pos_skipped_steps  # Store for later display
                results['pos_result'] = pos_enhanced_df
                progress_callback("📦 Queued Enhanced Positive for inclusion in Enhanced Combined workbook")

            # 3) Enhance Negative (in-memory only; will be written inside Enhanced Combined)
            enhanced_neg_output = None
            neg_skipped_steps = []
            if neg_src_df is not None and len(neg_src_df) > 0:
                progress_callback("⚙️ Enhancing negative polarity data in-memory")
                neg_result = cleaner.enhanced_stage2_cleaning(
                    neg_src_df,
                    ppm_threshold=enhanced_config['ppm_threshold'],
                    rt_threshold=enhanced_config['rt_threshold'],
                    ion_mode="negative",
                    ms2_filter_mode=enhanced_config.get('ms2_filter_mode', 'none'),
                    progress_callback=progress_callback,
                    selected_ions=enhanced_config.get('selected_ions'),  # NEW
                    column_clean_patterns=enhanced_config.get('column_clean_patterns')  # NEW
                )
                neg_enhanced_df = neg_result['data']
                neg_skipped_steps = neg_result.get('skipped_steps', [])
                self._last_neg_skipped_steps = neg_skipped_steps  # Store for later display
                results['neg_result'] = neg_enhanced_df
                progress_callback("📦 Queued Enhanced Negative for inclusion in Enhanced Combined workbook")

            # 4) Build and save Enhanced Combined workbook (only file we save in enhanced mode)
            enhanced_parts = [df for df in [results.get('pos_result'), results.get('neg_result')] if df is not None]
            # Use user-provided enhanced filename if set, else prefix the stage-1 base
            # Extract filename from output path if available
            enhanced_name = os.path.basename(self.metabolite_output_file.get()) if hasattr(self, 'metabolite_output_file') else None
            if enhanced_name and str(enhanced_name).strip():
                # Keep only the filename part if a full path was accidentally pasted
                enhanced_name = os.path.basename(str(enhanced_name).strip())
                # Ensure .xlsx extension
                if not enhanced_name.lower().endswith('.xlsx'):
                    enhanced_name = enhanced_name + '.xlsx'
                enhanced_combined_output = os.path.join(output_dir, enhanced_name)
            else:
                combined_filename_base = config['output'].get('combined_filename', 'cleaned_metabolites_combined')
                enhanced_combined_output = os.path.join(output_dir, f"enhanced_{combined_filename_base}.xlsx")

            if enhanced_parts:
                enhanced_combined_df = pd.concat(enhanced_parts, ignore_index=True)
                # Deduplicate combined data by Name to avoid redundant ID annotation searches
                original_count = len(enhanced_combined_df)
                if 'Name' in enhanced_combined_df.columns:
                    enhanced_combined_df = enhanced_combined_df.drop_duplicates(subset=['Name'], keep='first').reset_index(drop=True)
                    dedup_count = original_count - len(enhanced_combined_df)
                    if dedup_count > 0:
                        progress_callback(f"🔍 Deduplicated enhanced combined data: removed {dedup_count} duplicate metabolite(s)")
                        progress_callback(f"   Unique metabolites for ID annotation: {len(enhanced_combined_df)}")
                results['combined_result'] = enhanced_combined_df
                # Build Combined subset (Name, Formula)
                cols = [c for c in ['Name', 'Formula'] if c in enhanced_combined_df.columns]
                combined_subset = enhanced_combined_df[cols].drop_duplicates(subset=['Name'], keep='first').reset_index(drop=True) if cols else pd.DataFrame(columns=['Name', 'Formula'])
                # Expose the 2-column Combined subset in results for downstream consumers (annotator)
                results['combined_subset'] = combined_subset

                # Build Metabolites_IDS
                ids_combined = None
                try:
                    mz_parts = []
                    if config['cleaning_mode'] == 'separate':
                        if config['files'].get('pos_metabolika_file'):
                            mz_parts.append(pd.read_excel(config['files']['pos_metabolika_file']))
                        if config['files'].get('neg_metabolika_file'):
                            mz_parts.append(pd.read_excel(config['files']['neg_metabolika_file']))
                    else:
                        if config['files'].get('combined_metabolika_file'):
                            mz_parts.append(pd.read_excel(config['files']['combined_metabolika_file']))
                    if mz_parts:
                        mz_all = pd.concat(mz_parts, ignore_index=True)
                        ids_combined = cleaner._clean_metabolika_ids_df(mz_all)
                        results['metabolite_ids_result'] = ids_combined  # Store in results for memory
                except Exception:
                    ids_combined = None

                progress_callback(f"💾 Saving Enhanced Combined to: {os.path.basename(enhanced_combined_output)}")
                with pd.ExcelWriter(enhanced_combined_output, engine='openpyxl') as writer:
                    combined_subset.to_excel(writer, sheet_name='Combined', index=False)
                    if results.get('pos_result') is not None:
                        results['pos_result'].to_excel(writer, sheet_name='Positive', index=False)
                    if results.get('neg_result') is not None:
                        results['neg_result'].to_excel(writer, sheet_name='Negative', index=False)
                    if ids_combined is not None and not ids_combined.empty:
                        ids_combined.to_excel(writer, sheet_name='Metabolites_IDS', index=False)
                # Remove read-only attribute
                self._remove_readonly_attribute(enhanced_combined_output)
                # Only save the enhanced combined workbook
                results['files_created'] = [enhanced_combined_output]
                results['success'] = True
                results['message'] = "Two-stage enhanced cleaning completed successfully"
            else:
                # Fallback: enhance the combined dataframe as a single sheet if split failed
                combined_df_fallback = original_results.get('combined_result')
                if combined_df_fallback is not None:
                    progress_callback("⚙️ Enhancing combined data (fallback)")
                    # Choose a neutral ion_mode; positive is acceptable since stage2 will still apply other filters
                    combined_result = cleaner.enhanced_stage2_cleaning(
                        combined_df_fallback,
                        ppm_threshold=enhanced_config['ppm_threshold'],
                        rt_threshold=enhanced_config['rt_threshold'],
                        ion_mode="positive",
                        ms2_filter_mode=enhanced_config.get('ms2_filter_mode', 'none'),
                        progress_callback=progress_callback
                    )
                    combined_enhanced = combined_result['data']
                    combined_skipped_steps = combined_result.get('skipped_steps', [])
                    self._last_combined_skipped_steps = combined_skipped_steps  # Store for later display
                    # Deduplicate combined data by Name to avoid redundant ID annotation searches
                    original_count = len(combined_enhanced)
                    if 'Name' in combined_enhanced.columns:
                        combined_enhanced = combined_enhanced.drop_duplicates(subset=['Name'], keep='first').reset_index(drop=True)
                        dedup_count = original_count - len(combined_enhanced)
                        if dedup_count > 0:
                            progress_callback(f"🔍 Deduplicated enhanced combined data: removed {dedup_count} duplicate metabolite(s)")
                            progress_callback(f"   Unique metabolites for ID annotation: {len(combined_enhanced)}")
                    progress_callback(f"💾 Saving Enhanced Combined to: {os.path.basename(enhanced_combined_output)}")
                    with pd.ExcelWriter(enhanced_combined_output, engine='openpyxl') as writer:
                        combined_enhanced.to_excel(writer, sheet_name='Cleaned_Data', index=False)
                    # Remove read-only attribute
                    self._remove_readonly_attribute(enhanced_combined_output)
                    results['combined_result'] = combined_enhanced
                    # Only save the enhanced combined workbook
                    results['files_created'] = [enhanced_combined_output]
                    results['success'] = True
                    results['message'] = "Two-stage enhanced cleaning completed successfully"

            return results
            
        except Exception as e:
            return {
                'success': False,
                'files_created': [],
                'message': f"Enhanced cleaning failed: {str(e)}",
                'error': str(e)
            }
            
    def handle_cleaning_results(self, results):
        """Handle cleaning results in main thread"""
        # Reset UI elements
        self.main_progress_bar.config(value=100)
        self.progress_percentage.config(text="100%")
        self.start_button.config(state='normal', text="🚀 Start Data Cleaning")
        self.start_time = None  # Stop timer updates
        self.timer_manager.start_time = None  # Stop timer manager
        
        if results['success']:
            self.cleaned_data = results
            self.process_status.config(text="Data cleaning completed successfully!", fg='#27ae60')
            
            # Log completion to progress text
            self.append_to_log("="*50)
            self.append_to_log("✅ DATA CLEANING COMPLETED SUCCESSFULLY!")
            self.append_to_log("="*50)
            
            # Reset ETA to 00:00 when completed
            self.eta_label.config(text="ETA: 00:00")
            
            # Log summary information
            self.append_to_log(f"📁 Files created: {len(results['files_created'])}")
            
            if results.get('pos_result') is not None:
                self.append_to_log(f"➕ Positive metabolites: {len(results['pos_result'])}")
            if results.get('neg_result') is not None:
                self.append_to_log(f"➖ Negative metabolites: {len(results['neg_result'])}")
            if results.get('combined_result') is not None:
                self.append_to_log(f"🔗 Combined metabolites: {len(results['combined_result'])}")
            
            # Log created files and store the first one for folder opening
            if results.get('files_created'):
                self.append_to_log("\n📄 Output files:")
                for file_path in results['files_created']:
                    self.append_to_log(f"   • {os.path.basename(file_path)}")
                # Store the first created file path for the "Open Folder" button
                if results['files_created']:
                    self.last_cleaned_file_path = results['files_created'][0]
            
            # Display warnings about skipped cleaning steps due to missing columns
            all_skipped_steps = []
            if hasattr(self, '_last_pos_skipped_steps') and self._last_pos_skipped_steps:
                all_skipped_steps.extend([(step, reason, "Positive") for step, reason in self._last_pos_skipped_steps])
            if hasattr(self, '_last_neg_skipped_steps') and self._last_neg_skipped_steps:
                all_skipped_steps.extend([(step, reason, "Negative") for step, reason in self._last_neg_skipped_steps])
            if hasattr(self, '_last_combined_skipped_steps') and self._last_combined_skipped_steps:
                all_skipped_steps.extend([(step, reason, "Combined") for step, reason in self._last_combined_skipped_steps])
            
            if all_skipped_steps:
                self.append_to_log("\n" + "="*50)
                self.append_to_log("⚠️  OPTIONAL CLEANING STEPS SKIPPED")
                self.append_to_log("="*50)
                self.append_to_log("\nSome optional cleaning steps were skipped due to missing columns:")
                for step, reason, mode in all_skipped_steps:
                    self.append_to_log(f"  • [{mode}] {step}: {reason}")
                self.append_to_log("\n💡 To enable these steps:")
                self.append_to_log("   1. Click 'Verify Columns' before cleaning")
                self.append_to_log("   2. Assign the missing columns to their correct names")
                self.append_to_log("   3. Re-run the cleaning process")
                self.append_to_log("="*50)
            
            # Check if Group Area columns were missing during cleaning
            group_area_missing = False
            if isinstance(results.get('combined_result'), dict) and results['combined_result'].get('group_area_missing'):
                group_area_missing = True
            elif isinstance(results.get('pos_result'), dict) and results['pos_result'].get('group_area_missing'):
                group_area_missing = True
            elif isinstance(results.get('neg_result'), dict) and results['neg_result'].get('group_area_missing'):
                group_area_missing = True
            
            # Only show critical warning if Group Area was actually missing
            if group_area_missing:
                self.append_to_log("\n" + "="*50)
                self.append_to_log("⚠️⚠️⚠️ CRITICAL WARNING ⚠️⚠️⚠️")
                self.append_to_log("="*50)
                self.append_to_log("\nNo 'Group Area:' columns were detected!")
                self.append_to_log("\nThis indicates the 'Group Area' column assignment may be missing.")
                self.append_to_log("\n📋 RECOMMENDATION:")
                self.append_to_log("   Click the 'Verify Columns' button and explicitly assign")
                self.append_to_log("   the Group Area column before running data cleaning again.")
                self.append_to_log("\nThank you!")
                self.append_to_log("="*50 + "\n")
                
            self.process_status.config(text="✅ Cleaning completed successfully!", fg='green')
            
            # Store current data for next steps
            if results.get('combined_result') is not None:
                self.current_data = results['combined_result']
            elif results.get('pos_result') is not None:
                self.current_data = results['pos_result']
            
            # Set up memory store for DataFrame persistence (shared with DataManager)
            if hasattr(self, 'memory_store') and isinstance(self.memory_store, dict):
                self.memory_store.clear()
            else:
                # Fallback: re-link to shared DataManager store
                self.memory_store = self.data_manager.memory_store
            
            # Store enhanced cleaning results (from enhanced mode)
            # Store the 2-column Combined subset (Name, Formula) for ID annotation
            if results.get('combined_subset') is not None and isinstance(results['combined_subset'], pd.DataFrame):
                self.memory_store['cleaned_metabolites_df'] = results['combined_subset']
                self.append_to_log(f"🧠 Stored Combined sheet (Name+Formula) for ID annotation: {len(results['combined_subset'])} unique metabolites")
          
            if results.get('pos_result') is not None:
                self.memory_store['pos_enhanced_df'] = results['pos_result']
                self.append_to_log(f"🧠 Stored positive enhanced data: {len(results['pos_result'])} metabolites")
            if results.get('neg_result') is not None:
                self.memory_store['neg_enhanced_df'] = results['neg_result']
                self.append_to_log(f"🧠 Stored negative enhanced data: {len(results['neg_result'])} metabolites")
            
            # Store metabolite IDs data (from both enhanced and regular cleaning)
            if results.get('metabolite_ids_result') is not None:
                self.memory_store['metabolite_ids_df'] = results['metabolite_ids_result']
                self.append_to_log(f"🧠 Stored metabolite IDs data: {len(results['metabolite_ids_result'])} entries")
            elif results.get('metabolite_ids_df') is not None:
                self.memory_store['metabolite_ids_df'] = results['metabolite_ids_df']
                self.append_to_log(f"🧠 Stored metabolite IDs data: {len(results['metabolite_ids_df'])} entries")
            
            # Store regular cleaning results (from regular mode) if no enhanced data
            if results.get('combined_result') is None:
                if results.get('combined_df') is not None:
                    self.memory_store['cleaned_metabolites_df'] = results['combined_df']
                    self.append_to_log(f"🧠 Stored combined cleaned data: {len(results['combined_df'])} metabolites")
                elif results.get('pos_df') is not None:
                    self.memory_store['cleaned_metabolites_df'] = results['pos_df']
                    self.append_to_log(f"🧠 Stored positive cleaned data: {len(results['pos_df'])} metabolites")
            
            # Store cleaned file path for auto-loading to Tab 2
            files_created = results.get('files_created', [])
            # Prioritize enhanced combined, then cleaned combined, then enhanced positive, then any cleaned file
            self.cleaned_excel_path = None
            combined_file = None
            positive_file = None
            any_cleaned_file = None
            
            for file_path in files_created:
                filename = os.path.basename(file_path)
                name_lower = filename.lower()
                # Prefer enhanced combined first
                if name_lower.startswith('enhanced_') and 'combined' in name_lower:
                    combined_file = file_path
                elif 'cleaned_metabolites' in name_lower and 'combined' in name_lower and combined_file is None:
                    combined_file = file_path
                elif name_lower.startswith('enhanced_') and 'positive' in name_lower and positive_file is None:
                    positive_file = file_path
                elif ('cleaned_metabolites' in name_lower or name_lower.startswith('enhanced_')) and any_cleaned_file is None:
                    any_cleaned_file = file_path
            
            # Set priority: combined > positive > any cleaned file
            self.cleaned_excel_path = combined_file or positive_file or any_cleaned_file

            # Fallback: if we did not detect any file from results, use the
            # configured metabolite output path (if set).
            if not self.cleaned_excel_path:
                try:
                    if hasattr(self, 'metabolite_output_file') and self.metabolite_output_file.get():
                        self.cleaned_excel_path = self.metabolite_output_file.get()
                except Exception:
                    pass

            # Store in DataManager for other tabs to access
            if self.cleaned_excel_path:
                try:
                    if hasattr(self.data_manager, 'set_cleaned_excel_path'):
                        self.data_manager.set_cleaned_excel_path(self.cleaned_excel_path)
                    else:
                        self.data_manager.cleaned_excel_path = self.cleaned_excel_path
                except Exception:
                    self.data_manager.cleaned_excel_path = self.cleaned_excel_path

            # Log successful update
            if self.cleaned_excel_path:
                print(f"✅ Updated cleaned file path: {self.cleaned_excel_path}")
            
            # Auto-load to ID Annotation tab with shared memory_store (but don't auto-open the tab)
            try:
                # Notify ID Annotation tab that cleaned data is ready
                self.notify_data_ready("🔬 ID Annotation", "cleaned_metabolites_df")
                
                # Get ID Annotation tab and update its display
                id_annotation_tab = self.get_tab_by_name("🔬 ID Annotation")
                if id_annotation_tab and hasattr(id_annotation_tab, 'update_auto_loaded_file_display'):
                    id_annotation_tab.update_auto_loaded_file_display(self.cleaned_excel_path)
                    print(f"✅ Auto-loaded file to ID Annotation tab: {self.cleaned_excel_path}")
                    self.append_to_log(f"\n📚 Data has been automatically loaded to ID Annotation tab!")
                    self.append_to_log(f"    You can begin annotation whenever ready.")
                else:
                    print(f"Note: Could not access ID Annotation tab")
            except Exception as e:
                print(f"Note: Could not auto-load to ID Annotation tab: {e}")
            
            messagebox.showinfo("Success", "Data cleaning completed successfully!\n" + 
                              f"Files created: {', '.join(results['files_created'])}\n" +
                              "Switched to ID Annotation tab for next step.")
                              
        else:
            error_msg = results.get('error', 'Unknown error occurred')
            self.process_status.config(text="Error occurred during processing", fg='#e74c3c')
            self.update_progress_with_percentage(0, f"ERROR: {error_msg}")
            messagebox.showerror("Error", f"Data cleaning failed: {error_msg}")
            
    def handle_cleaning_error(self, error_msg):
        """Handle cleaning error in main thread"""
        # Reset UI elements
        self.main_progress_bar.config(value=0)
        self.progress_percentage.config(text="0%")
        self.start_button.config(state='normal', text="🚀 Start Data Cleaning")
        self.process_status.config(text="Error occurred during processing", fg='#e74c3c')
        self.start_time = None  # Stop timer updates
        self.timer_manager.start_time = None  # Stop timer manager
        
        # Reset ETA to --:-- on error
        self.eta_label.config(text="ETA: --:--")
        
        # Log error to progress text
        self.append_to_log("="*50)
        self.append_to_log("❌ ERROR OCCURRED")
        self.append_to_log("="*50)
        self.append_to_log(error_msg)
        
        messagebox.showerror("Error", error_msg)


