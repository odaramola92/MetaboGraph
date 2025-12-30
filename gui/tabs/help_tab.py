"""
Help Tab - User guide and documentation
"""
import tkinter as tk
from tkinter import ttk, scrolledtext
import logging

from gui.shared.base_tab import BaseTab

logger = logging.getLogger(__name__)


class HelpTab(BaseTab):
    """Tab for help and documentation"""
    
    def __init__(self, parent, data_manager):
        """Initialize the help tab"""
        super().__init__(parent, data_manager)
        self.setup_ui()
        print("[OK] Help Tab initialized")
    
    def setup_ui(self):
        """Setup the comprehensive Help tab with sub-tabs for each main tab"""
        # Title banner
        title_frame = tk.Frame(self.frame, bg='#27ae60', height=60)
        title_frame.pack(fill='x', pady=(0, 5))
        title_frame.pack_propagate(False)
        
        tk.Label(
            title_frame,
            text="User Guide & Help",
            font=('Arial', 16, 'bold'),
            bg='#27ae60',
            fg='white'
        ).pack(pady=15)
        
        # Create main content frame
        content_frame = tk.Frame(self.frame)
        content_frame.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Create notebook for sub-tabs
        help_notebook = ttk.Notebook(content_frame)
        help_notebook.pack(fill='both', expand=True)
        
        # Get help content
        help_contents = self._get_default_help_content()
        
        # Create sub-tabs for each main tab
        for tab_name, content in help_contents.items():
            self._create_help_subtab(help_notebook, tab_name, content)
        
        # Additional Resources section
        resources_frame = ttk.LabelFrame(self.frame, text="Additional Resources", padding=15)
        resources_frame.pack(fill='x', padx=5, pady=(0, 5))
        
        resources_text = "Documentation files available in the License/ folder"
        
        tk.Label(
            resources_frame,
            text=resources_text,
            font=('Arial', 9),
            bg='#ecf0f1',
            justify='left'
        ).pack(anchor='w')
    
    def _create_help_subtab(self, notebook, tab_name, content):
        """Create a help sub-tab with formatted content."""
        tab_frame = ttk.Frame(notebook)
        notebook.add(tab_frame, text=" " + tab_name + " ")
        
        # Create scrollable text area
        text_area = scrolledtext.ScrolledText(
            tab_frame,
            wrap=tk.WORD,
            font=('Arial', 10),
            bg='#ffffff',
            fg='#2c3e50',
            padx=15,
            pady=15
        )
        text_area.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Configure tags for formatting
        text_area.tag_config('title', font=('Arial', 14, 'bold'), foreground='#2c3e50')
        text_area.tag_config('heading', font=('Arial', 12, 'bold'), foreground='#3498db')
        text_area.tag_config('subheading', font=('Arial', 10, 'bold'), foreground='#27ae60')
        text_area.tag_config('bullet', font=('Arial', 10), foreground='#555')
        text_area.tag_config('warning', font=('Arial', 10), foreground='#e74c3c')
        text_area.tag_config('tip', font=('Arial', 10), foreground='#f39c12')
        
        # Insert content
        text_area.insert('1.0', content)
        text_area.config(state='disabled')
    
    def _get_default_help_content(self):
        """Provide default help content if help.py is not available."""
        return {
            "Data Cleaning": """DATA CLEANING TAB

This tab allows you to clean and prepare your metabolite data for analysis.

OVERVIEW:
The Data Cleaning tab processes raw metabolite data files and prepares them for ID annotation and pathway analysis.

KEY FEATURES:
- Excel file import
- Automatic data validation
- Missing value handling
- Duplicate removal
- Column standardization

HOW TO USE:
1. Click "Browse" to select your Excel file
2. Select the sheet containing your data
3. Click "Start Cleaning"
4. Review the cleaning report

TIPS:
- Ensure your Excel file has headers in the first row
- Remove any empty rows or columns before upload
""",
            
            "ID Annotation": """ID ANNOTATION TAB

This tab annotates metabolite IDs using HMDB and LipidMaps databases.

OVERVIEW:
The ID Annotation tab matches your metabolite names/IDs against comprehensive databases to provide standardized identifiers.

KEY FEATURES:
- HMDB database integration
- LipidMaps database for lipids
- Multiple ID format support
- Synonym matching

DATABASE REQUIREMENTS:
OK - hmdb_database.feather (required)
OK - lipidmap.feather (for lipid annotation)

TIPS:
- Process databases first in Database Setup tab
- Use multiple workers for faster processing
""",
            
            "Pathway Annotation": """PATHWAY ANNOTATION TAB

Annotate metabolites with pathway information.

OVERVIEW:
The Pathway Annotation tab maps metabolites to biological pathways using PathBank, SMPDB, and WikiPathways databases.

KEY FEATURES:
- Pathway enrichment analysis
- Network visualization
- Z-score calculation
- Statistical significance testing

DATABASE REQUIREMENTS:
OK - pathbank_selected.feather
OK - merged_SMP_metabolites.feather
OK - wikipathways databases

TIPS:
- Process pathway databases in Database Setup tab first
- Apply metabolite filters before pathway analysis
""",
            
            "Pathway Network": """PATHWAY NETWORK TAB

Analyze and visualize metabolite network interactions.

OVERVIEW:
Explore metabolite-pathway relationships, upstream regulators, and disease associations.

KEY FEATURES:
- Network graph generation
- Pathway selection and filtering
- Upstream regulator analysis
- Disease association exploration

TABS:
- Pathways: Select pathways for network
- Metabolites: View annotated metabolites
- Upstream: Analyze regulators
- Diseases: Explore associations

TIPS:
- Start with moderate filters
- Include key upstream regulators
- Export network in vector format
""",
            
            "Database Setup": """DATABASE SETUP TAB

Configure and manage metabolite databases.

OVERVIEW:
Handle downloading, processing, and verifying all required databases.

ID ANNOTATION:
- HMDB Database (https://hmdb.ca/downloads)
- LipidMaps (https://www.lipidmaps.org)

PATHWAY ANNOTATION:
- PathBank (https://pathbank.org/downloads)
- SMPDB (https://smpdb.ca/downloads)
- WikiPathways (https://zenodo.org/communities/wikipathways)

HOW TO SETUP:
1. Download raw database file
2. Click "Browse" and select file
3. Click "Process" button
4. Wait for processing

TIPS:
- Processing runs in background
- Databases only need processing once
- Click "Refresh Status" to check state
"""
        }
