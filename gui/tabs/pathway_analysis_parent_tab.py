"""
Pathway Analysis Parent Tab - Container for Pathway Annotation, Network Analysis, and Multi-Omics subtabs.

This parent tab creates a nested notebook structure similar to the Data Cleaning tab,
containing three subtabs in this order:
1. Pathway Annotation - Annotate metabolites/lipids with pathways (one-time process)
2. Network Analysis - Analyze and visualize networks for specific comparisons
3. Multi-Omics Integration - Merge multiple pathway-ready datasets before loading into the network tab
"""
import tkinter as tk
from tkinter import ttk
import logging

from gui.shared.base_tab import BaseTab, _setup_global_styles
from gui.tabs.pathway_annotation_tab import PathwayAnnotationTab
from gui.tabs.pathway_network_tab import PathwayNetworkTab
from gui.tabs.multiomics_analysis_tab import MultiOmicsAnalysisTab
from gui.tabs.comparative_analysis_tab import ComparativeAnalysisTab

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PathwayAnalysisParentTab(BaseTab):
    """
    Parent tab containing Pathway Annotation, Network Analysis, and Multi-Omics subtabs.
    
    Structure:
        Pathway Analysis (parent)
        ├── Pathway Annotation (subtab) - One-time annotation with auto-save
        ├── Network Analysis (subtab) - Comparison-specific visualization & export tools
        └── Multi-Omics Integration (subtab) - Merge multiple datasets before network loading
    """
    
    def __init__(self, parent, data_manager):
        """Initialize the parent tab with nested subtabs"""
        super().__init__(parent, data_manager)
        
        # Setup global styles
        _setup_global_styles()
        
        # Get root window for dialogs
        self.root = parent.winfo_toplevel()
        
        # Setup the UI (required pattern for BaseTab subclasses)
        self.setup_ui()
        
        logger.info("Pathway Analysis Parent Tab initialized")
    
    def setup_ui(self):
        """Setup the UI for this parent tab (required by BaseTab)"""
        # Create nested notebook for subtabs
        self.subtab_notebook = ttk.Notebook(self.frame)
        self.subtab_notebook.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Create subtab instances
        self._create_subtabs()
        
        logger.info("Pathway Analysis subtabs created")
    
    def _create_subtabs(self):
        """Create and add Pathway Annotation, Network Analysis, and Multi-Omics subtabs"""
        try:
            # Create Pathway Annotation subtab
            logger.debug("Creating Pathway Annotation subtab...")
            self.annotation_tab = PathwayAnnotationTab(self.subtab_notebook, self.data_manager)
            annotation_frame = self.annotation_tab.get_frame()
            self.subtab_notebook.add(annotation_frame, text="📝 Pathway Annotation")

            # Create Network Analysis subtab next (immediately after annotation)
            logger.debug("Creating Network Analysis subtab...")
            self.network_tab = PathwayNetworkTab(self.subtab_notebook, self.data_manager)
            network_frame = self.network_tab.get_frame()
            self.subtab_notebook.add(network_frame, text="🕸️ Network Analysis")

            # Create Multi-Omics Integration subtab last
            logger.debug("Creating Multi-Omics Integration subtab...")
            self.multiomics_tab = MultiOmicsAnalysisTab(self.subtab_notebook, self.data_manager)
            multiomics_frame = self.multiomics_tab.get_frame()
            self.subtab_notebook.add(multiomics_frame, text="🧬 Multi-Omics Integration")

            # Create Comparative Analysis subtab after multi-omics
            logger.debug("Creating Comparative Analysis subtab...")
            self.comparative_tab = ComparativeAnalysisTab(self.subtab_notebook, self.data_manager)
            comparative_frame = self.comparative_tab.get_frame()
            self.subtab_notebook.add(comparative_frame, text="📈 Comparative Analysis")

            # Hold on to frames for tab switching helpers
            self.annotation_frame = annotation_frame
            self.network_frame = network_frame
            self.multiomics_frame = multiomics_frame
            self.comparative_frame = comparative_frame

            # Wire tab references for data passing
            self.annotation_tab.network_tab_instance = self.network_tab
            self.network_tab.annotation_tab_instance = self.annotation_tab
            self.multiomics_tab.network_tab_instance = self.network_tab
            
            logger.info("Successfully created Pathway Analysis subtabs")
            
        except Exception as e:
            logger.error(f"Error creating Pathway Analysis subtabs: {e}", exc_info=True)
            raise
    
    def get_annotation_tab(self):
        """Get reference to Pathway Annotation subtab"""
        return self.annotation_tab
    
    def get_multiomics_tab(self):
        """Get reference to Multi-Omics Analysis subtab"""
        return self.multiomics_tab
    
    def get_network_tab(self):
        """Get reference to Network Analysis subtab"""
        return self.network_tab
    
    def switch_to_network_tab(self):
        """Switch to Network Analysis subtab (called after annotation completes)"""
        try:
            target = getattr(self, 'network_frame', None) or self.network_tab.get_frame()
            self.subtab_notebook.select(target)
            logger.debug("Switched to Network Analysis subtab")
        except Exception as e:
            logger.error(f"Error switching to Network Analysis tab: {e}")
    
    def switch_to_multiomics_tab(self):
        """Switch to Multi-Omics Analysis subtab"""
        try:
            target = getattr(self, 'multiomics_frame', None) or self.multiomics_tab.get_frame()
            self.subtab_notebook.select(target)
            logger.debug("Switched to Multi-Omics Analysis subtab")
        except Exception as e:
            logger.error(f"Error switching to Multi-Omics Analysis tab: {e}")

