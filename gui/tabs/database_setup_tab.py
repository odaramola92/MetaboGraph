"""
Database Setup Tab - Configure and manage metabolite databases
"""
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import logging
import os
import sys
from pathlib import Path
import pandas as pd
import threading
import time

from gui.shared.base_tab import BaseTab

logger = logging.getLogger(__name__)


class DatabaseSetupTab(BaseTab):
    """Tab for database setup and configuration"""
    
    def __init__(self, parent, data_manager):
        """Initialize the database setup tab"""
        print("[DEBUG] DatabaseSetupTab.__init__: Starting initialization...")
        super().__init__(parent, data_manager)
        # Get reference to root window for thread-safe GUI updates
        self.root = self.frame.winfo_toplevel()
        print(f"[DEBUG] DatabaseSetupTab.__init__: self.root = {self.root}")
        self.setup_ui()
        print("[DEBUG] DatabaseSetupTab.__init__: Initialization complete")
    
    def setup_ui(self):
        """Setup the Database Setup tab for database configuration."""
        print("[DEBUG] DatabaseSetupTab.setup_ui: Starting UI setup...")
        # Title (Full Width)
        title_frame = tk.Frame(self.frame, bg='#2c3e50', height=60)
        title_frame.pack(fill='x', pady=(0, 5))
        title_frame.pack_propagate(False)
        
        tk.Label(
            title_frame,
            text="🗄️ Database Setup & Configuration",
            font=('Arial', 16, 'bold'),
            bg='#2c3e50',
            fg='white'
        ).pack(pady=15)
        
        # Create main canvas with scrollbar for entire tab
        main_canvas_frame = tk.Frame(self.frame, bg='#ecf0f1')
        main_canvas_frame.pack(fill='both', expand=True, padx=0, pady=0)
        
        main_canvas = tk.Canvas(main_canvas_frame, bg='#ecf0f1', highlightthickness=0)
        main_scrollbar = ttk.Scrollbar(main_canvas_frame, orient="vertical", command=main_canvas.yview)
        
        # Scrollable content frame
        scrollable_content = tk.Frame(main_canvas, bg='#ecf0f1')
        
        scrollable_content.bind(
            "<Configure>",
            lambda e: main_canvas.configure(scrollregion=main_canvas.bbox("all"))
        )
        
        main_canvas_window = main_canvas.create_window((0, 0), window=scrollable_content, anchor="nw")
        main_canvas.configure(yscrollcommand=main_scrollbar.set)
        
        # Update canvas window width when canvas resizes
        def _configure_main_canvas(event):
            main_canvas.itemconfig(main_canvas_window, width=event.width)
        
        main_canvas.bind('<Configure>', _configure_main_canvas)
        
        # Pack canvas and scrollbar
        main_canvas.pack(side='left', fill='both', expand=True)
        main_scrollbar.pack(side='right', fill='y')
        
        # Mouse wheel binding for main canvas
        def _on_main_mousewheel(event):
            main_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        main_canvas.bind("<MouseWheel>", _on_main_mousewheel)
        scrollable_content.bind("<MouseWheel>", _on_main_mousewheel)
        
        # Main container for two columns - maximize screen usage
        main_container = tk.Frame(scrollable_content, bg='#ecf0f1')
        main_container.pack(fill='both', expand=True, padx=2, pady=2)
        
        # Configure grid for two equal columns - allow free expansion with minimum height
        main_container.grid_columnconfigure(0, weight=1, uniform="column")
        main_container.grid_columnconfigure(1, weight=1, uniform="column")
        # Set minimum height to ensure content is visible
        main_container.grid_rowconfigure(0, weight=1, minsize=800)
        
        # ============ LEFT COLUMN: ID ANNOTATION DATABASES ============
        left_frame = tk.LabelFrame(main_container, 
                                   text="🆔 ID Annotation Databases", 
                                   font=('Arial', 12, 'bold'),
                                   bg='#ecf0f1', 
                                   fg='#2c3e50',
                                   relief='raised',
                                   borderwidth=2)
        left_frame.grid(row=0, column=0, sticky='nsew', padx=(2, 1), pady=2)
        
        # Configure left frame grid
        left_frame.grid_rowconfigure(0, weight=1)
        left_frame.grid_columnconfigure(0, weight=1)
        
        # Left column canvas and scrollbar
        left_canvas = tk.Canvas(left_frame, bg='#ecf0f1', highlightthickness=0)
        left_scrollbar = ttk.Scrollbar(left_frame, orient="vertical", command=left_canvas.yview)
        left_scrollable = ttk.Frame(left_canvas, style='Custom.TFrame')
        
        left_scrollable.bind(
            "<Configure>",
            lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all"))
        )
        
        left_canvas_window = left_canvas.create_window((0, 0), window=left_scrollable, anchor="nw")
        left_canvas.configure(yscrollcommand=left_scrollbar.set)
        
        # Function to update canvas window width when canvas resizes
        def _configure_left_canvas(event):
            canvas_width = event.width
            left_canvas.itemconfig(left_canvas_window, width=canvas_width)
        
        left_canvas.bind('<Configure>', _configure_left_canvas)
        
        left_canvas.grid(row=0, column=0, sticky='nsew')
        left_scrollbar.grid(row=0, column=1, sticky='ns')
        
        # Mouse wheel binding for left column
        def _on_left_mousewheel(event):
            left_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        left_canvas.bind("<MouseWheel>", _on_left_mousewheel)
        left_scrollable.bind("<MouseWheel>", _on_left_mousewheel)
        
        # ============ RIGHT COLUMN: PATHWAY ANNOTATION DATABASES ============
        right_frame = tk.LabelFrame(main_container, 
                                    text="🛤️ Pathway Annotation Databases", 
                                    font=('Arial', 12, 'bold'),
                                    bg='#ecf0f1', 
                                    fg='#2c3e50',
                                    relief='raised',
                                    borderwidth=2)
        right_frame.grid(row=0, column=1, sticky='nsew', padx=(1, 2), pady=2)
        
        # Configure right frame grid
        right_frame.grid_rowconfigure(0, weight=1)
        right_frame.grid_columnconfigure(0, weight=1)
        
        # Right column canvas and scrollbar
        right_canvas = tk.Canvas(right_frame, bg='#ecf0f1', highlightthickness=0)
        right_scrollbar = ttk.Scrollbar(right_frame, orient="vertical", command=right_canvas.yview)
        right_scrollable = ttk.Frame(right_canvas, style='Custom.TFrame')
        
        right_scrollable.bind(
            "<Configure>",
            lambda e: right_canvas.configure(scrollregion=right_canvas.bbox("all"))
        )
        
        right_canvas_window = right_canvas.create_window((0, 0), window=right_scrollable, anchor="nw")
        right_canvas.configure(yscrollcommand=right_scrollbar.set)
        
        # Function to update canvas window width when canvas resizes
        def _configure_right_canvas(event):
            canvas_width = event.width
            right_canvas.itemconfig(right_canvas_window, width=canvas_width)
        
        right_canvas.bind('<Configure>', _configure_right_canvas)
        
        right_canvas.grid(row=0, column=0, sticky='nsew')
        right_scrollbar.grid(row=0, column=1, sticky='ns')
        
        # Mouse wheel binding for right column
        def _on_right_mousewheel(event):
            right_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        right_canvas.bind("<MouseWheel>", _on_right_mousewheel)
        right_scrollable.bind("<MouseWheel>", _on_right_mousewheel)
        
        # ============ LEFT COLUMN CONTENT ============
        
        # HMDB Database Section
        hmdb_frame = ttk.LabelFrame(left_scrollable, text="🧬 HMDB Database", padding=15)
        hmdb_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(
            hmdb_frame,
            text="Step 1: Download HMDB XML",
            font=('Arial', 10, 'bold'),
            bg='#ecf0f1'
        ).pack(anchor='w', pady=(0, 5))
        
        tk.Label(
            hmdb_frame,
            text="Visit: https://hmdb.ca/downloads\nDownload: hmdb_metabolites.xml",
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#555',
            justify='left'
        ).pack(anchor='w', pady=(0, 10))
        
        tk.Label(
            hmdb_frame,
            text="Step 2: Process the file",
            font=('Arial', 10, 'bold'),
            bg='#ecf0f1'
        ).pack(anchor='w', pady=(0, 5))
        
        hmdb_file_frame = ttk.Frame(hmdb_frame)
        hmdb_file_frame.pack(fill='x', pady=5)
        
        self.hmdb_xml_path = tk.StringVar()
        ttk.Entry(hmdb_file_frame, textvariable=self.hmdb_xml_path).pack(side='left', fill='x', expand=True, padx=(0, 5))
        ttk.Button(hmdb_file_frame, text="Browse", command=self.browse_hmdb_xml).pack(side='right')
        
        ttk.Button(hmdb_frame, text="🔄 Process HMDB", command=self.process_hmdb_database).pack(pady=10, fill='x')
        
        # Progress bar for HMDB
        self.hmdb_progress = ttk.Progressbar(hmdb_frame, mode='determinate', maximum=100)
        self.hmdb_progress.pack(fill='x', pady=5)
        
        self.hmdb_status = tk.StringVar(value="No file selected")
        tk.Label(
            hmdb_frame,
            textvariable=self.hmdb_status,
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#555',
            wraplength=400,
            justify='left'
        ).pack(anchor='w')
        
        # LipidMaps Database Section
        lipid_frame = ttk.LabelFrame(left_scrollable, text="💧 LipidMaps Database", padding=15)
        lipid_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(
            lipid_frame,
            text="Step 1: Download LipidMaps SDF",
            font=('Arial', 10, 'bold'),
            bg='#ecf0f1'
        ).pack(anchor='w', pady=(0, 5))
        
        tk.Label(
            lipid_frame,
            text="Visit: https://www.lipidmaps.org/databases/lmsd/download\nDownload: structures.sdf",
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#555',
            justify='left'
        ).pack(anchor='w', pady=(0, 10))
        
        tk.Label(
            lipid_frame,
            text="Step 2: Process the file",
            font=('Arial', 10, 'bold'),
            bg='#ecf0f1'
        ).pack(anchor='w', pady=(0, 5))
        
        lipid_file_frame = ttk.Frame(lipid_frame)
        lipid_file_frame.pack(fill='x', pady=5)
        
        self.lipid_sdf_path = tk.StringVar()
        ttk.Entry(lipid_file_frame, textvariable=self.lipid_sdf_path).pack(side='left', fill='x', expand=True, padx=(0, 5))
        ttk.Button(lipid_file_frame, text="Browse", command=self.browse_lipid_sdf).pack(side='right')
        
        ttk.Button(lipid_frame, text="🔄 Process LipidMaps", command=self.process_lipid_database).pack(pady=10, fill='x')
        
        # Progress bar for LipidMaps
        self.lipid_progress = ttk.Progressbar(lipid_frame, mode='determinate', maximum=100)
        self.lipid_progress.pack(fill='x', pady=5)
        
        self.lipid_status = tk.StringVar(value="No file selected")
        tk.Label(
            lipid_frame,
            textvariable=self.lipid_status,
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#555',
            wraplength=400,
            justify='left'
        ).pack(anchor='w')
        
        # ============ RIGHT COLUMN CONTENT ============
        
        # PathBank Database Section
        pathbank_frame = ttk.LabelFrame(right_scrollable, text="🧬 PathBank Database", padding=15)
        pathbank_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(
            pathbank_frame,
            text="Step 1: Download PathBank CSV",
            font=('Arial', 10, 'bold'),
            bg='#ecf0f1'
        ).pack(anchor='w', pady=(0, 5))
        
        tk.Label(
            pathbank_frame,
            text="Visit: https://pathbank.org/downloads\nDownload: pathbank_all_metabolites.csv\n⚠️ Filtered for Human/Rat/Mouse (Species column retained)",
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#555',
            justify='left'
        ).pack(anchor='w', pady=(0, 10))
        
        tk.Label(
            pathbank_frame,
            text="Step 2: Process the file",
            font=('Arial', 10, 'bold'),
            bg='#ecf0f1'
        ).pack(anchor='w', pady=(0, 5))
        
        pathbank_file_frame = ttk.Frame(pathbank_frame)
        pathbank_file_frame.pack(fill='x', pady=5)
        
        self.pathbank_csv_path = tk.StringVar()
        ttk.Entry(pathbank_file_frame, textvariable=self.pathbank_csv_path).pack(side='left', fill='x', expand=True, padx=(0, 5))
        ttk.Button(pathbank_file_frame, text="Browse", command=self.browse_pathbank_csv).pack(side='right')
        
        ttk.Button(pathbank_frame, text="🔄 Process PathBank", command=self.process_pathbank_database).pack(pady=10, fill='x')
        
        # Progress bar for PathBank
        self.pathbank_progress = ttk.Progressbar(pathbank_frame, mode='determinate', maximum=100)
        self.pathbank_progress.pack(fill='x', pady=5)
        
        self.pathbank_status = tk.StringVar(value="No file selected")
        tk.Label(
            pathbank_frame,
            textvariable=self.pathbank_status,
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#555',
            wraplength=400,
            justify='left'
        ).pack(anchor='w')
        
        # SMPDB Database Section
        smpdb_frame = ttk.LabelFrame(right_scrollable, text="🔬 SMPDB Database", padding=15)
        smpdb_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(
            smpdb_frame,
            text="Step 1: Download SMPDB CSV files",
            font=('Arial', 10, 'bold'),
            bg='#ecf0f1'
        ).pack(anchor='w', pady=(0, 5))
        
        tk.Label(
            smpdb_frame,
            text="Visit: https://smpdb.ca/downloads\nDownload: Metabolites CSV (ZIP), extract to folder\n⚠️ Merges ~48,687 pathway CSV files",
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#555',
            justify='left'
        ).pack(anchor='w', pady=(0, 10))
        
        tk.Label(
            smpdb_frame,
            text="Step 2: Select the folder",
            font=('Arial', 10, 'bold'),
            bg='#ecf0f1'
        ).pack(anchor='w', pady=(0, 5))
        
        smpdb_file_frame = ttk.Frame(smpdb_frame)
        smpdb_file_frame.pack(fill='x', pady=5)
        
        self.smpdb_folder_path = tk.StringVar()
        ttk.Entry(smpdb_file_frame, textvariable=self.smpdb_folder_path).pack(side='left', fill='x', expand=True, padx=(0, 5))
        ttk.Button(smpdb_file_frame, text="Browse Folder", command=self.browse_smpdb_folder).pack(side='right')
        
        ttk.Button(smpdb_frame, text="🔄 Process SMPDB", command=self.process_smpdb_database).pack(pady=10, fill='x')
        
        # Progress bar for SMPDB
        self.smpdb_progress = ttk.Progressbar(smpdb_frame, mode='determinate', maximum=100)
        self.smpdb_progress.pack(fill='x', pady=5)
        
        self.smpdb_status = tk.StringVar(value="No folder selected")
        tk.Label(
            smpdb_frame,
            textvariable=self.smpdb_status,
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#555',
            wraplength=400,
            justify='left'
        ).pack(anchor='w')
        
        # WikiPathways Database Section
        wikipathways_frame = ttk.LabelFrame(right_scrollable, text="🌐 WikiPathways Database", padding=15)
        wikipathways_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(
            wikipathways_frame,
            text="Step 1: Download WikiPathways GPML files",
            font=('Arial', 10, 'bold'),
            bg='#ecf0f1'
        ).pack(anchor='w', pady=(0, 5))
        
        tk.Label(
            wikipathways_frame,
            text="Visit: https://zenodo.org/communities/wikipathways\nDownload GPML archives for each organism:\n• wikipathways-*-gpml-Homo_sapiens.zip\n• wikipathways-*-gpml-Rattus_norvegicus.zip\n• wikipathways-*-gpml-Mus_musculus.zip\nExtract each ZIP to separate folders",
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#555',
            justify='left'
        ).pack(anchor='w', pady=(0, 10))
        
        tk.Label(
            wikipathways_frame,
            text="Step 2: Add organism folders",
            font=('Arial', 10, 'bold'),
            bg='#ecf0f1'
        ).pack(anchor='w', pady=(0, 5))
        
        # Storage for WikiPathways organism folders
        self.wikipathways_folders = {
            'Homo sapiens': tk.StringVar(),
            'Rattus norvegicus': tk.StringVar(),
            'Mus musculus': tk.StringVar()
        }
        
        # Create folder selection for each organism
        for organism in ['Homo sapiens', 'Rattus norvegicus', 'Mus musculus']:
            org_frame = ttk.Frame(wikipathways_frame)
            org_frame.pack(fill='x', pady=3)
            
            # Organism label
            org_label = tk.Label(
                org_frame,
                text=f"{organism}:",
                font=('Arial', 9, 'bold'),
                bg='#ecf0f1',
                width=18,
                anchor='w'
            )
            org_label.pack(side='left', padx=(0, 5))
            
            # Path entry
            org_entry = ttk.Entry(org_frame, textvariable=self.wikipathways_folders[organism])
            org_entry.pack(side='left', fill='x', expand=True, padx=(0, 5))
            
            # Browse button with organism context
            browse_btn = ttk.Button(
                org_frame,
                text="Browse",
                command=lambda org=organism: self.browse_wikipathways_folder(org)
            )
            browse_btn.pack(side='right')
        
        # Process button
        ttk.Button(
            wikipathways_frame,
            text="🔄 Process WikiPathways",
            command=self.process_wikipathways_database
        ).pack(pady=10, fill='x')
        
        # Progress bar for WikiPathways
        self.wikipathways_progress = ttk.Progressbar(wikipathways_frame, mode='determinate', maximum=100)
        self.wikipathways_progress.pack(fill='x', pady=5)
        
        self.wikipathways_status = tk.StringVar(value="No folders selected")
        tk.Label(
            wikipathways_frame,
            textvariable=self.wikipathways_status,
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#555',
            wraplength=400,
            justify='left'
        ).pack(anchor='w')
        
        # ============ DATABASE STATUS LOG (Full Width Bottom) ============
        status_frame = ttk.LabelFrame(main_container, text="📊 Database Status", padding=15)
        status_frame.grid(row=1, column=0, columnspan=2, sticky='nsew', padx=2, pady=(10, 2))
        
        # Configure row to allow status frame to expand freely
        main_container.grid_rowconfigure(1, weight=0)
        
        # Auto-refresh button and last updated label
        status_header_frame = tk.Frame(status_frame, bg='#ecf0f1')
        status_header_frame.pack(fill='x', pady=(0, 5))
        
        tk.Label(
            status_header_frame,
            text="Real-time database availability monitor",
            font=('Arial', 9, 'italic'),
            bg='#ecf0f1',
            fg='#555'
        ).pack(side='left')
        
        ttk.Button(
            status_header_frame,
            text="🔄 Refresh Status",
            command=self.update_database_status_log
        ).pack(side='right', padx=5)
        
        self.db_status_last_updated = tk.StringVar(value="Never")
        tk.Label(
            status_header_frame,
            textvariable=self.db_status_last_updated,
            font=('Arial', 8),
            bg='#ecf0f1',
            fg='#7f8c8d'
        ).pack(side='right')
        
        # Progress bar for status checking
        self.db_status_progress = ttk.Progressbar(
            status_frame,
            mode='determinate',
            maximum=100,
            length=200
        )
        self.db_status_progress.pack(fill='x', padx=5, pady=(0, 5))
        
        # Status log text area with increased height
        self.database_status_log = scrolledtext.ScrolledText(
            status_frame, 
            height=12, 
            font=('Courier', 9),
            bg='#ffffff',
            fg='#2c3e50',
            wrap=tk.WORD
        )
        self.database_status_log.pack(fill='both', expand=True, pady=5)
        
        # Configure tags for colored output
        self.database_status_log.tag_config('available', foreground='#27ae60', font=('Courier', 9, 'bold'))
        self.database_status_log.tag_config('missing', foreground='#e74c3c', font=('Courier', 9, 'bold'))
        self.database_status_log.tag_config('warning', foreground='#f39c12', font=('Courier', 9))
        self.database_status_log.tag_config('info', foreground='#3498db', font=('Courier', 9))
        self.database_status_log.tag_config('header', foreground='#2c3e50', font=('Courier', 10, 'bold'))
        
        # ============ FISHER ORA UNIVERSE MANAGEMENT ============
        # COMMENTED OUT: Universe is now calculated from the user's metabolite list
        # universe_frame = ttk.LabelFrame(main_container, text="🌍 Fisher ORA Universe (U) Management", padding=15)
        # universe_frame.grid(row=2, column=0, columnspan=2, sticky='nsew', padx=2, pady=(10, 2))
        # 
        # # Configure row to allow universe frame to expand
        # main_container.grid_rowconfigure(2, weight=0)
        # 
        # # Info text
        # universe_info = tk.Frame(universe_frame, bg='#ecf0f1')
        # universe_info.pack(fill='x', pady=(0, 10))
        # 
        # tk.Label(
        #     universe_info,
        #     text="The universe (U) is the total number of unique metabolites across all databases with pathway annotations.",
        #     font=('Arial', 9),
        #     bg='#ecf0f1',
        #     fg='#555',
        #     wraplength=700,
        #     justify='left'
        # ).pack(anchor='w')
        # 
        # tk.Label(
        #     universe_info,
        #     text="⚠️  Recalculate ONLY when you update or re-download databases. Otherwise, use the current cached values.",
        #     font=('Arial', 8, 'italic'),
        #     bg='#ecf0f1',
        #     fg='#e74c3c',
        #     wraplength=700,
        #     justify='left'
        # ).pack(anchor='w', pady=(5, 0))
        # 
        # # Universe status display
        # universe_status_frame = tk.Frame(universe_frame, bg='#ecf0f1')
        # universe_status_frame.pack(fill='x', pady=10)
        # 
        # # Load and display current universe values
        # self.universe_status_label = tk.Label(
        #     universe_status_frame,
        #     text="Loading universe status...",
        #     font=('Arial', 9),
        #     bg='#ecf0f1',
        #     fg='#2c3e50',
        #     justify='left'
        # )
        # self.universe_status_label.pack(anchor='w')
        # 
        # # Recalculate button
        # recalc_button_frame = tk.Frame(universe_frame, bg='#ecf0f1')
        # recalc_button_frame.pack(fill='x', pady=(0, 10))
        # 
        # ttk.Button(
        #     recalc_button_frame,
        #     text="🔄 Recalculate Universe (U)",
        #     command=self.recalculate_universe
        # ).pack(side='left', padx=5)
        # 
        # tk.Label(
        #     recalc_button_frame,
        #     text="(Use only when databases have been updated)",
        #     font=('Arial', 8, 'italic'),
        #     bg='#ecf0f1',
        #     fg='#7f8c8d'
        # ).pack(side='left', padx=5)
        # 
        # # Progress bar for universe calculation
        # self.universe_progress = ttk.Progressbar(
        #     universe_frame,
        #     mode='indeterminate',
        #     length=400
        # )
        # self.universe_progress.pack(fill='x', padx=5, pady=5)
        
        # Initial status check in background (non-blocking)
        print("[DEBUG] DatabaseSetupTab.setup_ui: Updating initial database status...")
        self.update_database_status_log()
        # print("[DEBUG] DatabaseSetupTab.setup_ui: Loading universe status...")
        # self.load_universe_status()  # Commented out - universe now calculated from user's metabolite list
        print("[DEBUG] DatabaseSetupTab.setup_ui: UI setup complete")

    # COMMENTED OUT: Universe is now calculated from the user's metabolite list
    # def load_universe_status(self):
    #     """Load and display current universe status from config file."""
    #     import json
    #     from pathlib import Path
    #     from datetime import datetime
    #     
    #     print("[DEBUG] load_universe_status: Starting...")
    #     
    #     try:
    #         # Check both locations: root and main_script
    #         config_paths = [
    #             Path(__file__).parent.parent.parent / "universe_config.json",  # Root
    #             Path(__file__).parent.parent.parent / "main_script" / "universe_config.json"  # main_script
    #         ]
    #         
    #         print(f"[DEBUG] Checking config paths:")
    #         for i, path in enumerate(config_paths):
    #             exists = path.exists()
    #             print(f"[DEBUG]   Path {i+1}: {path}")
    #             print(f"[DEBUG]   Exists: {exists}")
    #         
    #         config_path = None
    #         for path in config_paths:
    #             if path.exists():
    #                 config_path = path
    #                 print(f"[DEBUG] Using config from: {config_path}")
    #                 break
    #         
    #         if config_path:
    #             print(f"[DEBUG] Loading config from: {config_path}")
    #             with open(config_path, 'r') as f:
    #                 config = json.load(f)
    #             
    #             print(f"[DEBUG] Config loaded: {config.keys()}")
    #             universes = config.get('universes', {})
    #             calc_date = config.get('calculated_date', 'Unknown')
    #             print(f"[DEBUG] Universes found: {universes}")
    #             
    #             # Parse datetime
    #             try:
    #                 calc_dt = datetime.fromisoformat(calc_date)
    #                 calc_date_str = calc_dt.strftime('%Y-%m-%d %H:%M:%S')
    #             except:
    #                 calc_date_str = calc_date
    #             
    #             # Format status text
    #             status_text = "Current Universe (U) Values:\n"
    #             for species, u_value in universes.items():
    #                 status_text += f"  • {species:25s}: U = {u_value:,} metabolites\n"
    #             status_text += f"\nLast Calculated: {calc_date_str}"
    #             
    #             print(f"[DEBUG] Setting status label to show universe values")
    #             self.universe_status_label.config(text=status_text, fg='#27ae60')
    #         else:
    #             # Config not found - auto-generate it
    #             print("[DEBUG] Config not found - triggering auto-generation")
    #             self.universe_status_label.config(
    #                 text="⏳ Universe config not found. Auto-generating...",
    #                 fg='#f39c12'
    #             )
    #             # Auto-calculate universe in background
    #             self._auto_calculate_universe()
    #     except Exception as e:
    #         error_msg = f"Error loading universe status: {str(e)}"
    #         print(f"[DEBUG ERROR] {error_msg}")
    #         import traceback
    #         traceback.print_exc()
    #         self.universe_status_label.config(
    #             text=error_msg,
    #             fg='#e74c3c'
    #         )
    
    # def _auto_calculate_universe(self):
    #     """Automatically calculate universe in the background without user confirmation."""
    #     import subprocess
    #     from pathlib import Path
    #     
    #     print("[DEBUG] _auto_calculate_universe: Starting automatic universe calculation...")
    #     
    #     def run_calculation():
    #         try:
    #             # Get the path to calculate_universe.py in main_script folder
    #             script_path = Path(__file__).parent.parent.parent / "main_script" / "calculate_universe.py"
    #             print(f"[DEBUG] Looking for calculate_universe.py at: {script_path}")
    #             print(f"[DEBUG] Script exists: {script_path.exists()}")
    #             
    #             if not script_path.exists():
    #                 error_msg = f"calculate_universe.py not found at {script_path}"
    #                 print(f"[DEBUG ERROR] {error_msg}")
    #                 raise FileNotFoundError(error_msg)
    #             
    #             # Get Python executable
    #             python_exe = Path(os.sys.executable)
    #             print(f"[DEBUG] Using Python executable: {python_exe}")
    #             
    #             # Run the script
    #             print(f"[DEBUG] Running command: {python_exe} {script_path}")
    #             result = subprocess.run(
    #                 [str(python_exe), str(script_path)],
    #                 capture_output=True,
    #                 text=True,
    #                 timeout=600  # 10 minute timeout
    #             )
    #             
    #             print(f"[DEBUG] Subprocess return code: {result.returncode}")
    #             print(f"[DEBUG] STDOUT:\n{result.stdout}")
    #             if result.stderr:
    #                 print(f"[DEBUG] STDERR:\n{result.stderr}")
    #             
    #             # Update status based on result
    #             if result.returncode == 0:
    #                 print("[DEBUG] Universe calculation successful - reloading status")
    #                 self.root.after(0, lambda: self._on_auto_universe_complete(True, None))
    #             else:
    #                 error_msg = result.stderr or result.stdout or "Unknown error"
    #                 print(f"[DEBUG ERROR] Universe calculation failed with exit code {result.returncode}")
    #                 print(f"[DEBUG ERROR] Error output: {error_msg[:200]}")
    #                 self.root.after(0, lambda: self._on_auto_universe_complete(False, error_msg))
    #                 
    #         except subprocess.TimeoutExpired:
    #             error_msg = "Auto-calculation timed out (>10 minutes)"
    #             print(f"[DEBUG ERROR] {error_msg}")
    #             self.root.after(0, lambda: self._on_auto_universe_complete(False, error_msg))
    #         except Exception as e:
    #             error_msg = str(e)
    #             print(f"[DEBUG ERROR] Exception in auto-calculation: {error_msg}")
    #             import traceback
    #             traceback.print_exc()
    #             self.root.after(0, lambda: self._on_auto_universe_complete(False, error_msg))
    #     
    #     # Run in background thread
    #     print("[DEBUG] Starting calculation thread...")
    #     calculation_thread = threading.Thread(target=run_calculation, daemon=True)
    #     calculation_thread.start()
    #     print("[DEBUG] Calculation thread started")
    
    # def _on_auto_universe_complete(self, success, error):
    #     """Handle completion of automatic universe calculation."""
    #     print(f"[DEBUG] _on_auto_universe_complete called: success={success}, error={error}")
    #     
    #     if success:
    #         print("[DEBUG] Auto-calculation succeeded - reloading universe status")
    #         self.universe_status_label.config(
    #             text="✅ Universe auto-generated successfully! Loading...",
    #             fg='#27ae60'
    #         )
    #         # Reload the status to show the new values
    #         self.load_universe_status()
    #         print("[DEBUG] Universe status reloaded")
    #     else:
    #         error_short = str(error)[:100] if error else "Unknown error"
    #         print(f"[DEBUG] Auto-calculation failed: {error_short}")
    #         self.universe_status_label.config(
    #             text=f"⚠️ Auto-generation failed. Click 'Recalculate Universe' button.\n{error_short}",
    #             fg='#e74c3c'
    #         )
    
    # def recalculate_universe(self):
    #     """Recalculate Fisher ORA universe from databases."""
    #     import subprocess
    #     from pathlib import Path
    #     
    #     # Confirm with user
    #     result = messagebox.askyesno(
    #         "Recalculate Universe",
    #         "This will scan all database files and recalculate the universe size for all species.\n\n"
    #         "This process may take a few minutes.\n\n"
    #         "Do you want to proceed?"
    #     )
    #     
    #     if not result:
    #         return
    #     
    #     # Disable button and show progress
    #     self.universe_progress.start()
    #     self.universe_status_label.config(text="Calculating universe... Please wait...", fg='#f39c12')
    #     
    #     def run_calculation():
    #         try:
    #             # Get the path to calculate_universe.py in main_script folder
    #             script_path = Path(__file__).parent.parent.parent / "main_script" / "calculate_universe.py"
    #             
    #             if not script_path.exists():
    #                 raise FileNotFoundError(f"calculate_universe.py not found at {script_path}")
    #             
    #             # Run the script
    #             result = subprocess.run(
    #                 [str(Path(os.sys.executable)), str(script_path)],
    #                 capture_output=True,
    #                 text=True,
    #                 timeout=600  # 10 minute timeout
    #             )
    #             
    #             # Update status based on result
    #             if result.returncode == 0:
    #                 self.root.after(0, lambda: self._on_universe_calculation_complete(True, None))
    #             else:
    #                 error_msg = result.stderr or "Unknown error"
    #                 self.root.after(0, lambda: self._on_universe_calculation_complete(False, error_msg))
    #                 
    #         except subprocess.TimeoutExpired:
    #             self.root.after(0, lambda: self._on_universe_calculation_complete(False, "Calculation timed out"))
    #         except Exception as e:
    #             self.root.after(0, lambda: self._on_universe_calculation_complete(False, str(e)))
    #     
    #     # Run in background thread
    #     calculation_thread = threading.Thread(target=run_calculation, daemon=True)
    #     calculation_thread.start()
    
    # def _on_universe_calculation_complete(self, success, error):
    #     """Handle completion of universe calculation."""
    #     self.universe_progress.stop()
    #     
    #     if success:
    #         messagebox.showinfo(
    #             "Success",
    #             "Universe (U) has been successfully recalculated!\n\n"
    #             "Fisher ORA will use the updated values in the next analysis."
    #         )
    #         self.load_universe_status()
    #     else:
    #         messagebox.showerror(
    #             "Error",
    #             f"Failed to recalculate universe:\n\n{error}\n\n"
    #             "Please check the console output for details."
    #         )
    #         self.universe_status_label.config(
    #             text="Error during calculation. Check console output.",
    #             fg='#e74c3c'
    #         )
    
    def browse_hmdb_xml(self):
        """Browse for HMDB XML file."""
        filename = filedialog.askopenfilename(
            title="Select HMDB XML File",
            filetypes=[("XML files", "*.xml"), ("All files", "*.*")]
        )
        if filename:
            self.hmdb_xml_path.set(filename)
            self.hmdb_status.set(f"Selected: {os.path.basename(filename)}")
    
    def browse_lipid_sdf(self):
        """Browse for LipidMaps SDF file."""
        filename = filedialog.askopenfilename(
            title="Select LipidMaps SDF File",
            filetypes=[("SDF files", "*.sdf"), ("All files", "*.*")]
        )
        if filename:
            self.lipid_sdf_path.set(filename)
            self.lipid_status.set(f"Selected: {os.path.basename(filename)}")
    
    def browse_pathbank_csv(self):
        """Browse for PathBank CSV file."""
        filename = filedialog.askopenfilename(
            title="Select PathBank CSV File",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if filename:
            self.pathbank_csv_path.set(filename)
            self.pathbank_status.set(f"Selected: {os.path.basename(filename)}")
    
    def browse_smpdb_folder(self):
        """Browse for SMPDB folder containing CSV files."""
        foldername = filedialog.askdirectory(
            title="Select SMPDB Folder (containing SMP*.csv files)"
        )
        if foldername:
            self.smpdb_folder_path.set(foldername)
            csv_count = len([f for f in os.listdir(foldername) if f.endswith('.csv')])
            self.smpdb_status.set(f"Selected: {csv_count} CSV files in folder")
    
    def browse_wikipathways_folder(self, organism):
        """Browse for WikiPathways GPML folder for specific organism."""
        foldername = filedialog.askdirectory(
            title=f"Select WikiPathways GPML Folder for {organism}"
        )
        if foldername:
            self.wikipathways_folders[organism].set(foldername)
            gpml_count = len([f for f in os.listdir(foldername) if f.endswith('.gpml')])
            self.wikipathways_status.set(f"✓ {organism}: {gpml_count} GPML files")
    
    def log_database_message(self, message):
        """Log a message to the database status log."""
        if hasattr(self, 'database_status_log'):
            self.database_status_log.insert(tk.END, f"{message}\n")
            self.database_status_log.see(tk.END)
            self.database_status_log.update()
    
    def process_hmdb_database(self):
        """Process HMDB database"""
        xml_file = self.hmdb_xml_path.get().strip()
        
        if not xml_file or not os.path.exists(xml_file):
            messagebox.showerror("Error", "Please select a valid HMDB XML file first")
            return
        
        # Run in background thread
        def process():
            try:
                # Start progress bar
                self.root.after(0, lambda: self.hmdb_progress.config(mode='indeterminate'))
                self.root.after(0, lambda: self.hmdb_progress.start(10))
                
                self.log_database_message("🔄 Processing HMDB database...")
                self.log_database_message("   This may take 5-10 minutes depending on file size...")
                self.root.after(0, lambda: self.hmdb_status.set("Processing... Please wait (5-10 min estimate)"))
                
                from main_script.database_builder import HMDBDatabaseBuilder
                
                builder = HMDBDatabaseBuilder(xml_file)
                
                # Determine output directory - use Databases folder
                output_dir = os.path.join(os.getcwd(), 'Databases')
                os.makedirs(output_dir, exist_ok=True)
                
                hmdb_feather, synonyms_feather = builder.build(output_dir)
                
                # Stop progress bar
                self.root.after(0, lambda: self.hmdb_progress.stop())
                self.root.after(0, lambda: self.hmdb_progress.config(mode='determinate', value=100))
                
                self.log_database_message(f"✅ HMDB processing complete!")
                self.log_database_message(f"   Created: {os.path.basename(hmdb_feather)}")
                self.log_database_message(f"   Created: {os.path.basename(synonyms_feather)}")
                self.log_database_message(f"   Location: {output_dir}")
                
                self.root.after(0, lambda: self.hmdb_status.set(f"✅ Complete! Files saved to Databases folder"))
                self.root.after(0, lambda: messagebox.showinfo("Success", f"HMDB database processed successfully!\n\nFiles created in: {output_dir}\n- {os.path.basename(hmdb_feather)}\n- {os.path.basename(synonyms_feather)}"))
            except Exception as e:
                self.root.after(0, lambda: self.hmdb_progress.stop())
                self.root.after(0, lambda: self.hmdb_progress.config(value=0))
                self.log_database_message(f"❌ HMDB processing failed: {e}")
                self.root.after(0, lambda: self.hmdb_status.set(f"❌ Failed: {str(e)[:50]}"))
                self.root.after(0, lambda: messagebox.showerror("Error", f"Failed to process HMDB database:\n{e}"))
        
        threading.Thread(target=process, daemon=True).start()
    
    def process_lipid_database(self):
        """Process LipidMaps database"""
        sdf_file = self.lipid_sdf_path.get().strip()
        
        if not sdf_file or not os.path.exists(sdf_file):
            messagebox.showerror("Error", "Please select a valid LipidMaps SDF file first")
            return
        
        # Run in background thread
        def process():
            try:
                # Start progress bar
                self.root.after(0, lambda: self.lipid_progress.config(mode='indeterminate'))
                self.root.after(0, lambda: self.lipid_progress.start(10))
                
                self.log_database_message("🔄 Processing LipidMaps database...")
                self.log_database_message("   This may take 3-5 minutes...")
                self.root.after(0, lambda: self.lipid_status.set("Processing... Please wait (3-5 min estimate)"))
                
                from main_script.database_builder import LipidMapsDatabaseBuilder
                
                builder = LipidMapsDatabaseBuilder(sdf_file)
                
                # Determine output directory - use Databases folder
                output_dir = os.path.join(os.getcwd(), 'Databases')
                os.makedirs(output_dir, exist_ok=True)
                
                lipid_feather = builder.build(output_dir)
                
                # Stop progress bar
                self.root.after(0, lambda: self.lipid_progress.stop())
                self.root.after(0, lambda: self.lipid_progress.config(mode='determinate', value=100))
                
                self.log_database_message(f"✅ LipidMaps processing complete!")
                self.log_database_message(f"   Created: {os.path.basename(lipid_feather)}")
                self.log_database_message(f"   Location: {output_dir}")
                
                self.root.after(0, lambda: self.lipid_status.set(f"✅ Complete! Files saved to Databases folder"))
                self.root.after(0, lambda: messagebox.showinfo("Success", f"LipidMaps database processed successfully!\n\nFile created in: {output_dir}\n- {os.path.basename(lipid_feather)}"))
            except Exception as e:
                self.root.after(0, lambda: self.lipid_progress.stop())
                self.root.after(0, lambda: self.lipid_progress.config(value=0))
                self.log_database_message(f"❌ LipidMaps processing failed: {e}")
                self.root.after(0, lambda: self.lipid_status.set(f"❌ Failed: {str(e)[:50]}"))
                self.root.after(0, lambda: messagebox.showerror("Error", f"Failed to process LipidMaps database:\n{e}"))
        
        threading.Thread(target=process, daemon=True).start()
    
    def process_pathbank_database(self):
        """Process PathBank database"""
        csv_file = self.pathbank_csv_path.get().strip()
        
        if not csv_file or not os.path.exists(csv_file):
            messagebox.showerror("Error", "Please select a valid PathBank CSV file first")
            return
        
        # Run in background thread
        def process():
            try:
                # Start progress bar
                self.root.after(0, lambda: self.pathbank_progress.config(mode='indeterminate'))
                self.root.after(0, lambda: self.pathbank_progress.start(10))
                
                self.log_database_message("🔄 Processing PathBank database...")
                self.log_database_message("   This may take 2-3 minutes...")
                self.root.after(0, lambda: self.pathbank_status.set("Processing... Please wait (2-3 min estimate)"))
                
                from main_script.database_builder import PathBankDatabaseBuilder
                
                builder = PathBankDatabaseBuilder(csv_file)
                
                # Determine output directory - use Databases folder
                output_dir = os.path.join(os.getcwd(), 'Databases')
                os.makedirs(output_dir, exist_ok=True)
                
                pathbank_feather = builder.build(output_dir)
                
                # Stop progress bar
                self.root.after(0, lambda: self.pathbank_progress.stop())
                self.root.after(0, lambda: self.pathbank_progress.config(mode='determinate', value=100))
                
                self.log_database_message(f"✅ PathBank processing complete!")
                self.log_database_message(f"   Created: {os.path.basename(pathbank_feather)}")
                self.log_database_message(f"   Location: {output_dir}")
                
                self.root.after(0, lambda: self.pathbank_status.set(f"✅ Complete! Files saved to Databases folder"))
                self.root.after(0, lambda: messagebox.showinfo("Success", f"PathBank database processed successfully!\n\nFile created in: {output_dir}\n- {os.path.basename(pathbank_feather)}"))
            except Exception as e:
                self.root.after(0, lambda: self.pathbank_progress.stop())
                self.root.after(0, lambda: self.pathbank_progress.config(value=0))
                self.log_database_message(f"❌ PathBank processing failed: {e}")
                self.root.after(0, lambda: self.pathbank_status.set(f"❌ Failed: {str(e)[:50]}"))
                self.root.after(0, lambda: messagebox.showerror("Error", f"Failed to process PathBank database:\n{e}"))
        
        threading.Thread(target=process, daemon=True).start()
    
    def process_smpdb_database(self):
        """Process SMPDB database"""
        folder = self.smpdb_folder_path.get().strip()
        
        if not folder or not os.path.exists(folder):
            messagebox.showerror("Error", "Please select a valid SMPDB folder first")
            return
        
        # Run in background thread
        def process():
            try:
                # Start progress bar
                self.root.after(0, lambda: self.smpdb_progress.config(mode='indeterminate'))
                self.root.after(0, lambda: self.smpdb_progress.start(10))
                
                self.log_database_message("🔄 Processing SMPDB database...")
                self.log_database_message("   This may take 10-15 minutes (merging ~48k files)...")
                self.root.after(0, lambda: self.smpdb_status.set("Processing... Please wait (10-15 min estimate)"))
                
                from main_script.database_builder import SMPDBDatabaseBuilder
                
                builder = SMPDBDatabaseBuilder(folder)
                
                # Determine output directory - use Databases folder
                output_dir = os.path.join(os.getcwd(), 'Databases')
                os.makedirs(output_dir, exist_ok=True)
                
                smpdb_feather = builder.build(output_dir)
                
                # Stop progress bar
                self.root.after(0, lambda: self.smpdb_progress.stop())
                self.root.after(0, lambda: self.smpdb_progress.config(mode='determinate', value=100))
                
                self.log_database_message(f"✅ SMPDB processing complete!")
                self.log_database_message(f"   Created: {os.path.basename(smpdb_feather)}")
                self.log_database_message(f"   Location: {output_dir}")
                
                self.root.after(0, lambda: self.smpdb_status.set(f"✅ Complete! Files saved to Databases folder"))
                self.root.after(0, lambda: messagebox.showinfo("Success", f"SMPDB database processed successfully!\n\nFile created in: {output_dir}\n- {os.path.basename(smpdb_feather)}"))
            except Exception as e:
                self.root.after(0, lambda: self.smpdb_progress.stop())
                self.root.after(0, lambda: self.smpdb_progress.config(value=0))
                self.log_database_message(f"❌ SMPDB processing failed: {e}")
                self.root.after(0, lambda: self.smpdb_status.set(f"❌ Failed: {str(e)[:50]}"))
                self.root.after(0, lambda: messagebox.showerror("Error", f"Failed to process SMPDB database:\n{e}"))
        
        threading.Thread(target=process, daemon=True).start()
    
    def process_wikipathways_database(self):
        """Process WikiPathways database"""
        # Check if at least one organism folder is selected
        selected_organisms = {}
        for organism, path_var in self.wikipathways_folders.items():
            folder = path_var.get().strip()
            if folder and os.path.exists(folder):
                selected_organisms[organism] = folder
        
        if not selected_organisms:
            messagebox.showerror("Error", "Please select at least one WikiPathways folder first")
            return
        
        # Run in background thread
        def process():
            try:
                # Start progress bar
                self.root.after(0, lambda: self.wikipathways_progress.config(mode='indeterminate'))
                self.root.after(0, lambda: self.wikipathways_progress.start(10))
                
                self.log_database_message("🔄 Processing WikiPathways database...")
                self.log_database_message("   This may take 5-8 minutes...")
                self.root.after(0, lambda: self.wikipathways_status.set("Processing... Please wait (5-8 min estimate)"))
                
                from main_script.database_builder import WikiPathwaysDatabaseBuilder
                
                # Determine output directory - use Databases folder
                output_dir = os.path.join(os.getcwd(), 'Databases')
                os.makedirs(output_dir, exist_ok=True)
                
                builder = WikiPathwaysDatabaseBuilder(output_dir)
                created_files = []
                
                # Process each organism
                for organism, folder in selected_organisms.items():
                    self.log_database_message(f"   Processing {organism}...")
                    # Build expects organism_name parameter
                    organism_files = builder.build(output_dir)
                    if isinstance(organism_files, list):
                        created_files.extend(organism_files)
                    else:
                        created_files.append(organism_files)
                
                # Stop progress bar
                self.root.after(0, lambda: self.wikipathways_progress.stop())
                self.root.after(0, lambda: self.wikipathways_progress.config(mode='determinate', value=100))
                
                self.log_database_message(f"✅ WikiPathways processing complete!")
                for file in created_files:
                    self.log_database_message(f"   Created: {os.path.basename(file)}")
                
                files_list = "\n- ".join([os.path.basename(f) for f in created_files])
                self.root.after(0, lambda: self.wikipathways_status.set(f"✅ Complete! Files saved to Databases folder"))
                self.root.after(0, lambda: messagebox.showinfo("Success", f"WikiPathways database processed successfully!\n\nFiles created:\n- {files_list}"))
            except Exception as e:
                self.root.after(0, lambda: self.wikipathways_progress.stop())
                self.root.after(0, lambda: self.wikipathways_progress.config(value=0))
                self.log_database_message(f"❌ WikiPathways processing failed: {e}")
                self.root.after(0, lambda: self.wikipathways_status.set(f"❌ Failed: {str(e)[:50]}"))
                self.root.after(0, lambda: messagebox.showerror("Error", f"Failed to process WikiPathways database:\n{e}"))
        
        threading.Thread(target=process, daemon=True).start()
    
    def update_database_status_log(self):
        """Check and display database file status"""
        from gui.shared.utils import find_database_paths
        
        self.database_status_log.config(state='normal')
        self.database_status_log.delete('1.0', tk.END)
        
        # Get the prioritized list of database search paths
        possible_paths = find_database_paths()
        
        self.database_status_log.insert(tk.END, "🔍 Scanning database files...\n\n", 'header')
        
        databases_path = None
        for path in possible_paths:
            # Check if any .feather files exist in this path
            if os.path.isdir(path):
                try:
                    feather_files = [f for f in os.listdir(path) if f.endswith('.feather')]
                    if feather_files:
                        databases_path = path
                        self.database_status_log.insert(tk.END, f"✓ Found {len(feather_files)} database files in: {path}\n", 'info')
                        break
                except Exception as e:
                    # Skip paths that can't be read
                    continue
        
        if not databases_path:
            self.database_status_log.insert(tk.END, "❌ No database files (.feather) found\n\n", 'missing')
            self.database_status_log.insert(tk.END, f"Searched in:\n", 'info')
            for path in possible_paths:
                exists = "✓ exists" if os.path.isdir(path) else "✗ not found"
                self.database_status_log.insert(tk.END, f"  • {path} ({exists})\n", 'info')
            self.database_status_log.insert(tk.END, "\n💡 Tip: Place .feather database files in the same folder as the executable\n", 'info')
            self.database_status_log.insert(tk.END, "✨ Database status check complete\n", 'info')
            self.database_status_log.config(state='disabled')
            self.db_status_last_updated.set(f"Last updated: {time.strftime('%H:%M:%S')}")
            return
        
        self.database_status_log.insert(tk.END, f"\n📁 Database location: {databases_path}\n\n", 'info')
        
        # Check core databases
        core_databases = [
            ("hmdb_database.feather", "HMDB Core Database"),
            ("lipidmap.feather", "LipidMaps Core Database"),
            ("merged_SMP_metabolites.feather", "Merged SMP Metabolites"),
            ("pathbank_selected.feather", "PathBank Database"),
            ("wikipathways_homo_sapiens.feather", "WikiPathways - Homo Sapiens"),
            ("wikipathways_mus_musculus.feather", "WikiPathways - Mus Musculus"),
            ("wikipathways_rattus_norvegicus.feather", "WikiPathways - Rattus Norvegicus"),
        ]
        
        for db_file, description in core_databases:
            full_path = os.path.join(databases_path, db_file)
            if os.path.exists(full_path):
                size_mb = os.path.getsize(full_path) / (1024**2)
                self.database_status_log.insert(tk.END, f"✅ {db_file} ({description})\n   Size: {size_mb:.1f} MB\n\n", 'available')
            else:
                self.database_status_log.insert(tk.END, f"❌ {db_file} ({description}) - NOT FOUND\n\n", 'missing')
        
        self.database_status_log.insert(tk.END, "✨ Database status check complete\n", 'info')
        self.database_status_log.config(state='disabled')
        
        # Update last checked time
        self.db_status_last_updated.set(f"Last updated: {time.strftime('%H:%M:%S')}")
