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
from datetime import datetime

# Add current directory to path so we can import gui module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Setup file-based logging for diagnostics (especially useful for frozen executables)
try:
    import tempfile
    log_dir = tempfile.gettempdir()
    log_file = os.path.join(log_dir, f"metabograph_gui_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    
    # Create both console and file handlers
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('[%(asctime)s] %(name)s - %(levelname)s: %(message)s', datefmt='%H:%M:%S')
    file_handler.setFormatter(file_formatter)
except Exception:
    file_handler = None
    log_file = None

# Setup root logger
logging.basicConfig(
    level=logging.DEBUG,  # More verbose logging for diagnostics
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)

if file_handler:
    logger.addHandler(file_handler)
    logger.debug(f"File logging enabled: {log_file}")

# Import shared components
from gui.shared.data_manager import DataManager
from gui.shared.base_tab import BaseTab

logger.debug("Imported gui.shared components")

# Import tabs that are ready
from gui.tabs.data_cleaning_tab import DataCleaningTab
from gui.tabs.id_annotation_tab import IDAnnotationTab
from gui.tabs.pathway_analysis_parent_tab import PathwayAnalysisParentTab
from gui.tabs.database_setup_tab import DatabaseSetupTab
from gui.tabs.help_tab import HelpTab

logger.debug("Imported all tab modules")


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
        try:
            logger.info("=" * 80)
            logger.info("Initializing MetaboGraph GUI")
            logger.info("=" * 80)
            
            # Create root window
            logger.debug("Creating root Tk window...")
            self.root = tk.Tk()
            logger.debug("[OK] Root window created successfully")
            
            self.root.title("MetaboGraph - Metabolite Annotation & Pathway Analysis")
            self.root.geometry("1000x600")
            self.root.configure(bg='#f0f0f0')
            logger.debug("[OK] Root window configured")
            
            # Configure ttk styles for better appearance
            logger.debug("Setting up ttk styles...")
            self._setup_styles()
            logger.debug("[OK] ttk styles configured")
            
            # Initialize shared data manager (all tabs will use this)
            logger.debug("Initializing DataManager...")
            self.data_manager = DataManager()
            logger.info("[OK] DataManager initialized successfully")
            
            # Create title bar
            logger.debug("Creating title bar...")
            self._setup_title_bar()
            logger.debug("[OK] Title bar created")
            
            # Create notebook (tab container)
            logger.debug("Creating notebook (tab container)...")
            self.notebook = ttk.Notebook(self.root)
            self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
            logger.debug("[OK] Notebook created and packed")
            
            # Store tab instances for inter-tab communication
            self.tab_instances = {}
            
            # Shutdown flag
            self._shutting_down = False
            self.root.protocol("WM_DELETE_WINDOW", self._on_closing)
            logger.debug("[OK] Shutdown protocol registered")
            
            # Setup all tabs
            logger.debug("Setting up tabs...")
            self._setup_tabs()
            logger.info("[OK] All tabs initialized successfully")
            
            logger.info("=" * 80)
            logger.info("GUI initialization complete - ready for event loop")
            logger.info("=" * 80)
            
        except Exception as e:
            logger.exception(f"[ERROR] Fatal error during GUI initialization: {e}")
            raise
    
    def _setup_tabs(self):
        """
        Setup all tabs by instantiating tab classes.
        
        Each tab is a separate class that handles its own UI.
        This method just instantiates them and adds them to the notebook.
        """
        logger.info("Setting up tabs...")
        
        # Define all tabs: (display_name, TabClass)
        tabs_to_create = [
            ("Data Cleaning", DataCleaningTab),
            ("ID Annotation", IDAnnotationTab),
            ("Pathway Analysis", PathwayAnalysisParentTab),
            ("Database Setup", DatabaseSetupTab),
            ("Help", HelpTab),
        ]
        
        # Try to create each tab
        successful_tabs = 0
        for tab_name, TabClass in tabs_to_create:
            if self._add_tab(tab_name, TabClass):
                successful_tabs += 1
        
        logger.info(f"[OK] Successfully loaded {successful_tabs}/{len(tabs_to_create)} tabs")
        
        if successful_tabs == 0:
            logger.error("[ERROR] No tabs loaded successfully - application may not be functional")
    
    def _add_tab(self, tab_name: str, TabClass):
        """
        Add a single tab to the notebook.
        
        Args:
            tab_name: Display name for the tab
            TabClass: The tab class to instantiate
            
        Returns:
            True if tab loaded successfully, False otherwise
        """
        try:
            logger.debug(f"Loading tab: {tab_name}")
            
            # Instantiate the tab
            logger.debug(f"  [-] Instantiating {TabClass.__name__}...")
            tab_instance = TabClass(self.notebook, self.data_manager)
            logger.debug(f"  [OK] {TabClass.__name__} instantiated")
            
            # Store tab instance for inter-tab communication
            self.tab_instances[tab_name] = tab_instance
            
            # Tabs in this project build their UI during __init__,
            # so do not call setup_ui() again to avoid duplicating widgets.
            # If a future tab requires an explicit setup call, handle it there.
            
            # Get the frame and add it to the notebook
            logger.debug(f"  [-] Getting frame from tab...")
            frame = tab_instance.get_frame()
            logger.debug(f"  [OK] Frame retrieved")
            
            logger.debug(f"  [-] Adding frame to notebook...")
            self.notebook.add(frame, text=tab_name)
            logger.debug(f"  [OK] Frame added to notebook")
            
            logger.info(f"[OK] Successfully loaded tab: {tab_name}")
            return True
        
        except Exception as e:
            logger.error(f"[ERROR] Error loading tab '{tab_name}': {e}", exc_info=True)
            # Add error tab instead
            try:
                error_frame = ttk.Frame(self.notebook)
                error_label = ttk.Label(
                    error_frame,
                    text=f"Error loading {tab_name}:\n{str(e)}",
                    foreground="red"
                )
                error_label.pack(pady=20, padx=20)
                self.notebook.add(error_frame, text=tab_name)
                logger.debug(f"  [OK] Error tab placeholder added")
            except Exception as frame_error:
                logger.error(f"[ERROR] Failed to add error placeholder for {tab_name}: {frame_error}")
            
            return False
    
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
        """Create title bar"""
        # Title frame (dark blue background)
        title_frame = tk.Frame(self.root, bg='#2c3e50', height=60)
        title_frame.pack(fill='x', pady=(0, 5))
        title_frame.pack_propagate(False)
        
        # Title label (centered)
        title_label = tk.Label(
            title_frame,
            text="MetaboGraph - Metabolite Annotation & Pathway Analysis",
            font=('Arial', 16, 'bold'),
            fg='white',
            bg='#2c3e50'
        )
        title_label.pack(side='left', padx=12, pady=8)
    
    def _on_closing(self):
        """Handle window closing"""
        logger.info("Application closing")
        self._shutting_down = True
        
        # Close all matplotlib figures before destroying tkinter window
        # This prevents "main thread is not in main loop" error at exit
        try:
            import matplotlib.pyplot as plt
            plt.close('all')
        except Exception:
            pass  # Matplotlib might not be imported
        
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
    logger.info("=" * 80)
    logger.info("MetaboGraph Main Application Entry Point")
    logger.info("=" * 80)
    
    try:
        logger.info("Creating MetaboGraphGUI instance...")
        app = MetaboGraphGUI()
        logger.info("[OK] MetaboGraphGUI instance created successfully")
        
        logger.info("Starting event loop...")
        app.run()
        logger.info("[OK] Application exited normally")
        
    except Exception as e:
        logger.critical(f"[ERROR] Fatal error: {e}", exc_info=True)
        
        # Show error dialog if running as frozen executable (PyInstaller)
        try:
            import tkinter as tk
            from tkinter import messagebox
            
            logger.info("Attempting to show error dialog...")
            root = tk.Tk()
            root.withdraw()  # Hide the root window
            
            error_msg = f"MetaboGraph Failed to Start\n\n{str(e)}"
            if log_file:
                error_msg += f"\n\nDiagnostic log:\n{log_file}"
            
            messagebox.showerror("MetaboGraph Error", error_msg)
            root.destroy()
            logger.info("[OK] Error dialog shown")
        except Exception as dialog_error:
            logger.error(f"Failed to show error dialog: {dialog_error}", exc_info=True)
        
        sys.exit(1)


if __name__ == "__main__":
    main()
