"""
MetaboGraph GUI - Main Application Entry Point

This module coordinates all MetaboGraph tabs:
- Data Cleaning
- ID Annotation  
- Pathway Analysis
- Database Setup
- Help
"""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import logging
import sys
import os
import json

# Add current directory to path so we can import gui module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

# Import shared components
from gui.shared.data_manager import DataManager
from gui.shared.base_tab import BaseTab

# Import tabs that are ready
from gui.tabs.data_cleaning_tab import DataCleaningTab
from gui.tabs.id_annotation_tab import IDAnnotationTab
from gui.tabs.pathway_analysis_parent_tab import PathwayAnalysisParentTab
from gui.tabs.database_setup_tab import DatabaseSetupTab
from gui.tabs.help_tab import HelpTab


class MetaboGraphGUI:
    """
    Main MetaboGraph GUI application using modular tab architecture.
    
    This class:
    1. Creates the main window
    2. Initializes the shared DataManager
    3. Creates tab instances
    4. Adds tabs to the notebook (tabbed interface)
    
    MetaboGraph includes:
    - Data Cleaning
    - ID Annotation
    - Pathway Analysis
    """
    
    def __init__(self):
        """Initialize the main GUI application"""
        logger.info("Initializing MetaboGraph GUI")
        
        # Create root window
        self.root = tk.Tk()
        self.root.title("MetaboGraph - Metabolite Annotation & Pathway Analysis")
        self.root.geometry("1000x600")
        self.root.configure(bg='#f0f0f0')
        
        # Configure ttk styles for better appearance
        self._setup_styles()
        
        # Initialize shared data manager (all tabs will use this)
        self.data_manager = DataManager()
        logger.info("DataManager initialized")
        
        # Create title bar with settings controls
        self._setup_title_bar()
        
        # Create notebook (tab container)
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Store reference to this GUI instance on the notebook for inter-tab communication
        self.notebook.gui_instance = self
        
        # Store tab instances for inter-tab communication
        self.tab_instances = {}
        
        # Shutdown flag
        self._shutting_down = False
        self.root.protocol("WM_DELETE_WINDOW", self._on_closing)
        
        # Setup all tabs
        self._setup_tabs()
        
        logger.info("GUI initialization complete")
    
    def _setup_tabs(self):
        """
        Setup all tabs by instantiating tab classes.
        
        Each tab is a separate class that handles its own UI.
        This method just instantiates them and adds them to the notebook.
        """
        logger.info("Setting up tabs...")
        
        # Define all tabs: (display_name, TabClass)
        tabs_to_create = [
            ("📊 Data Cleaning", DataCleaningTab),
            ("🔬 ID Annotation", IDAnnotationTab),
            ("🧬 Pathway Analysis", PathwayAnalysisParentTab),
            ("🗄️ Database Setup", DatabaseSetupTab),
            ("❓ Help", HelpTab),
        ]
        
        # Try to create each tab
        for tab_name, TabClass in tabs_to_create:
            self._add_tab(tab_name, TabClass)
        
        logger.info(f"Loaded {len(tabs_to_create)} tabs")
    
    def _add_tab(self, tab_name: str, TabClass):
        """
        Add a single tab to the notebook.
        
        Args:
            tab_name: Display name for the tab
            TabClass: The tab class to instantiate
        """
        try:
            logger.debug(f"Loading tab: {tab_name}")
            
            # Instantiate the tab
            tab_instance = TabClass(self.notebook, self.data_manager)
            
            # Store tab instance for inter-tab communication
            self.tab_instances[tab_name] = tab_instance
            
            # Tabs in this project build their UI during __init__,
            # so do not call setup_ui() again to avoid duplicating widgets.
            # If a future tab requires an explicit setup call, handle it there.
            
            # Get the frame and add it to the notebook
            frame = tab_instance.get_frame()
            self.notebook.add(frame, text=tab_name)
            
            logger.debug(f"Successfully loaded tab: {tab_name}")
        
        except Exception as e:
            logger.error(f"Error loading tab '{tab_name}': {e}", exc_info=True)
            # Add error tab instead
            error_frame = ttk.Frame(self.notebook)
            error_label = ttk.Label(
                error_frame,
                text=f"Error loading {tab_name}:\n{str(e)}",
                foreground="red"
            )
            error_label.pack(pady=20, padx=20)
            self.notebook.add(error_frame, text=tab_name)
    
    def _setup_styles(self):
        """Configure ttk styles for better appearance"""
        style = ttk.Style()
        style.theme_use('clam')  # Use clam theme for better customization
        
        # Configure tab font - smaller size
        style.configure('TNotebook.Tab', font=('Arial', 10), padding=[15, 8])
        
        # Configure notebook background to match old GUI
        style.configure('TNotebook', background="#bfbfbf", borderwidth=1)
        
        # Configure tab colors - match old GUI gray scheme
        # Unselected tabs: lighter gray, Selected tabs: medium gray with white text
        style.configure('TNotebook.Tab', background="#bfbfbf", foreground="#000000")
        style.map('TNotebook.Tab',
                  background=[('selected', '#bfbfbf'), ('!selected', '#bfbfbf')],
                  foreground=[('selected', '#000000'), ('!selected', '#000000')],
                  lightcolor=[('selected', '#bfbfbf')],
                  darkcolor=[('selected', '#666666')])
        
        # Configure active/hover state
        style.configure('TNotebook.Tab', borderwidth=2, relief='raised')
        style.map('TNotebook.Tab', relief=[('selected', 'sunken'), ('!selected', 'raised')])
    
    def _setup_title_bar(self):
        """Create title bar with settings controls"""
        # Title frame (dark blue background)
        title_frame = tk.Frame(self.root, bg='#2c3e50', height=60)
        title_frame.pack(fill='x', pady=(0, 5))
        title_frame.pack_propagate(False)
        
        # Title label (left side)
        title_label = tk.Label(
            title_frame,
            text="🧬 MetaboGraph - Metabolite Annotation & Pathway Analysis",
            font=('Arial', 16, 'bold'),
            fg='white',
            bg='#2c3e50'
        )
        title_label.pack(side='left', padx=12, pady=8)
        
        # Settings controls (right side)
        settings_frame = tk.Frame(title_frame, bg='#2c3e50')
        settings_frame.pack(side='right', padx=10, pady=8)
        
        # Auto-load checkbox
        self.auto_load_settings = tk.BooleanVar(value=True)
        auto_load_check = tk.Checkbutton(
            settings_frame,
            text='Auto-load settings on start',
            variable=self.auto_load_settings,
            bg='#2c3e50',
            fg='white',
            selectcolor='#2c3e50',
            font=('Arial', 9)
        )
        auto_load_check.pack(side='left', padx=(0, 10))
        
        # Load Settings button
        load_btn = tk.Button(
            settings_frame,
            text='Load Settings',
            command=self.load_settings_dialog,
            bg='#3498db',
            fg='white',
            font=('Arial', 9, 'bold'),
            padx=8,
            pady=4,
            relief='raised',
            cursor='hand2'
        )
        load_btn.pack(side='left', padx=3)
        
        # Save Settings button
        save_btn = tk.Button(
            settings_frame,
            text='Save Settings',
            command=self.save_settings_dialog,
            bg='#27ae60',
            fg='white',
            font=('Arial', 9, 'bold'),
            padx=8,
            pady=4,
            relief='raised',
            cursor='hand2'
        )
        save_btn.pack(side='left', padx=3)
    
    def save_settings_dialog(self):
        """Save current settings to JSON file"""
        try:
            path = filedialog.asksaveasfilename(
                title='Save GUI Settings',
                defaultextension='.json',
                filetypes=[('JSON files', '*.json'), ('All files', '*.*')]
            )
            if not path:
                return
            
            # Collect settings from all tabs
            settings = {'version': '1.0', 'tabs': {}}
            
            # Get settings from each tab if they have the method
            for tab_name, tab_instance in self.tab_instances.items():
                if hasattr(tab_instance, 'get_settings'):
                    try:
                        settings['tabs'][tab_name] = tab_instance.get_settings()
                    except Exception as e:
                        logger.warning(f"Could not get settings from {tab_name}: {e}")
            
            # Save to file
            with open(path, 'w', encoding='utf-8') as f:
                json.dump(settings, f, indent=2, default=str)
            
            logger.info(f"Settings saved to {path}")
            messagebox.showinfo("Success", f"Settings saved successfully to:\n{path}")
        except Exception as e:
            logger.error(f"Error saving settings: {e}")
            messagebox.showerror("Error", f"Failed to save settings:\n{str(e)}")
    
    def load_settings_dialog(self):
        """Load settings from JSON file"""
        try:
            path = filedialog.askopenfilename(
                title='Load GUI Settings',
                filetypes=[('JSON files', '*.json'), ('All files', '*.*')]
            )
            if not path:
                return
            
            # Load from file
            with open(path, 'r', encoding='utf-8') as f:
                settings = json.load(f)
            
            # Apply settings to each tab if they have the method
            for tab_name, tab_instance in self.tab_instances.items():
                if hasattr(tab_instance, 'set_settings'):
                    tab_settings = settings.get('tabs', {}).get(tab_name, {})
                    try:
                        tab_instance.set_settings(tab_settings)
                    except Exception as e:
                        logger.warning(f"Could not apply settings to {tab_name}: {e}")
            
            logger.info(f"Settings loaded from {path}")
            messagebox.showinfo("Success", f"Settings loaded successfully from:\n{path}")
        except Exception as e:
            logger.error(f"Error loading settings: {e}")
            messagebox.showerror("Error", f"Failed to load settings:\n{str(e)}")
    
    def _on_closing(self):
        """Handle window closing"""
        logger.info("Application closing")
        self._shutting_down = True
        self.root.destroy()
    
    def run(self):
        """Start the application event loop"""
        logger.info("Starting application event loop")
        self.root.mainloop()
    
    def get_data_manager(self):
        """Get the shared data manager"""
        return self.data_manager


def main():
    """Entry point for the MetaboGraph application"""
    try:
        app = MetaboGraphGUI()
        app.run()
    except Exception as e:
        logger.critical(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
