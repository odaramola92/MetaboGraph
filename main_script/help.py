#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metabolite Annotation Tool - Help System
Standalone help script that can be run independently or imported by the GUI.

Usage:
    python help.py                 # Show help menu
    python help.py --tab cleaning  # Show specific tab help
    from help import get_help_content  # Import in GUI
"""

import sys
import io

# Set UTF-8 encoding for console output
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

def get_help_content():
    """
    Returns comprehensive help content for all tabs in the Metabolite Annotation Tool.
    
    Returns:
        dict: Dictionary with tab names as keys and help content as values
    """
    return {
        "Data Cleaning": """📋 DATA CLEANING TAB

This tab allows you to clean and prepare your metabolite data for analysis.

OVERVIEW:
The Data Cleaning tab processes raw metabolite data files and prepares them for ID annotation and pathway analysis.

KEY FEATURES:
• Excel file import (.xlsx, .xls)
• Automatic data validation
• Missing value handling
• Duplicate removal
• Column standardization
• Data quality reporting

HOW TO USE:
1. Click "Browse" to select your Excel file
2. Select the sheet containing your data
3. Configure cleaning options:
   - Remove rows with missing values
   - Handle duplicates
   - Standardize column names
4. Click "Start Cleaning"
5. Review the cleaning report
6. Cleaned file is automatically saved with "_cleaned" suffix

REQUIRED COLUMNS:
• Metabolite name or ID (any column with metabolite identifiers)
• P-value (for statistical filtering)
• Fold change (log2FC recommended)
• Optional: Additional metadata columns

OUTPUT:
• Cleaned Excel file saved in same directory
• Summary report showing:
  - Number of rows processed
  - Missing values handled
  - Duplicates removed
  - Columns standardized

TIPS:
⭐ Ensure your Excel file has headers in the first row
⭐ Remove any empty rows or columns before upload
⭐ Use consistent metabolite naming conventions
⭐ Check for special characters in metabolite names
⭐ Save original file as backup before cleaning
        """,
        
        "ID Annotation": """🔬 ID ANNOTATION TAB

This tab annotates metabolite IDs using HMDB and LipidMaps databases.

OVERVIEW:
The ID Annotation tab matches your metabolite names/IDs against comprehensive databases to provide standardized identifiers and additional information.

KEY FEATURES:
• HMDB database integration (>220,000 metabolites)
• LipidMaps database for lipids (>47,000 lipids)
• Multiple ID format support
• Synonym matching with fuzzy search
• Automatic database selection
• Multi-threaded processing for speed
• Confidence scoring for matches

HOW TO USE:
1. Load cleaned data (auto-loaded from Data Cleaning tab)
   OR manually select an Excel file
2. Configure annotation settings:
   - Number of workers (4-8 recommended)
   - Match threshold (0.7-0.9)
   - Database priority
3. Click "Start ID Annotation"
4. Wait for processing to complete (progress bar shown)
5. Review annotated results in tabs:
   - All Results: Complete annotation data
   - Filtered Results: High-confidence matches
   - Summary: Statistics and quality metrics

SUPPORTED ID FORMATS:
• HMDB IDs (HMDB0000001, HMDB00001)
• ChEBI IDs (CHEBI:15377)
• PubChem CIDs (CID5280450)
• KEGG IDs (C00031)
• Common metabolite names
• Systematic chemical names
• Synonyms and trivial names

DATABASE REQUIREMENTS:
✅ hmdb_database.feather (required for all metabolites)
✅ lipidmap.feather (required for lipid annotation)

MATCHING ALGORITHM:
1. Exact ID match (highest priority)
2. Exact name match
3. Synonym match
4. Fuzzy name match (similarity ≥ threshold)
5. Chemical formula match (if available)

OUTPUT COLUMNS:
• Original metabolite name
• Matched HMDB ID
• Matched name
• Match confidence score
• Database source
• Chemical formula
• Molecular weight
• Classification

TIPS:
⭐ Process databases first in Database Setup tab
⭐ Use 4-8 workers for optimal performance
⭐ Higher match threshold (0.9) = fewer but higher quality matches
⭐ Lower match threshold (0.7) = more matches but check quality
⭐ Check the filtered results tab for quality control
⭐ Review low-confidence matches manually
        """,
        
        "Statistics": """📊 STATISTICS TAB

Perform statistical analysis on annotated metabolite data.

OVERVIEW:
The Statistics tab provides comprehensive tools for analyzing metabolite abundance, fold changes, and statistical significance.

KEY FEATURES:
• Differential expression analysis
• Volcano plots
• PCA (Principal Component Analysis)
• Heatmap generation
• Statistical summaries and reporting
• Multiple testing correction
• Data normalization options

HOW TO USE:
1. Load ID-annotated data (from ID Annotation tab)
2. Select analysis type:
   - Differential Expression
   - PCA
   - Clustering
3. Configure statistical parameters:
   - Significance threshold (p-value)
   - Fold change cutoff
   - Correction method
4. Run analysis
5. Review results and plots
6. Export results and visualizations

STATISTICAL METHODS:
• T-test (paired/unpaired)
• ANOVA (one-way, two-way)
• Multiple testing correction:
  - False Discovery Rate (FDR/Benjamini-Hochberg)
  - Bonferroni correction
  - Holm-Bonferroni
• Fold change calculations (linear or log2)

FILTERING OPTIONS:
• P-value threshold: 0.05, 0.01, 0.001, or custom
• Log2FC threshold: ±0.5, ±1.0, ±1.5, or custom
• Minimum detection threshold
• Sample size requirements

OUTPUTS:
• Statistical summary tables
• Significant metabolites list
• Volcano plot (log2FC vs -log10(p-value))
• PCA score plots
• Loading plots
• Heatmaps with hierarchical clustering

DATA REQUIREMENTS:
• Minimum 3 replicates per group (recommended)
• Numerical abundance/intensity values
• Missing value handling options
• Normalization recommended

TIPS:
⭐ Ensure sufficient sample size for statistical power
⭐ Apply appropriate multiple testing correction (FDR recommended)
⭐ Check data distribution before analysis (normality tests)
⭐ Use log-transformation for abundance data
⭐ Visualize data before and after normalization
⭐ Document analysis parameters for reproducibility
        """,
        
        "Visualization": """📈 VISUALIZATION TAB

Create publication-quality plots and charts.

OVERVIEW:
The Visualization tab generates various plots to visualize metabolite data and analysis results with professional styling and export options.

KEY FEATURES:
• Interactive plots with zoom/pan
• Multiple chart types
• Customizable styling and themes
• Export options (PNG, PDF, SVG, EPS)
• Color scheme selection
• Font and size customization
• Multi-plot layouts

AVAILABLE PLOTS:
• Volcano plots (differential expression)
• Bar charts (metabolite comparisons)
• Heatmaps (hierarchical clustering)
• PCA plots (2D and 3D)
• Box plots (distribution comparisons)
• Scatter plots (correlation analysis)
• Network graphs (pathway visualization)
• Venn diagrams (overlap analysis)

HOW TO USE:
1. Load analyzed data
2. Select plot type from dropdown
3. Configure plot parameters:
   - X/Y axis variables
   - Color coding scheme
   - Point/bar sizes
   - Labels and titles
4. Preview the visualization
5. Adjust styling as needed
6. Export or save plot

CUSTOMIZATION OPTIONS:
• Color schemes: Viridis, Plasma, Set1, Set2, Custom
• Point markers: Circle, Square, Triangle, Diamond
• Line styles: Solid, Dashed, Dotted
• Font families: Arial, Times, Helvetica, Courier
• Font sizes: Title, axis labels, tick labels
• Figure dimensions: Width, height, DPI

EXPORT FORMATS:
• PNG: High-resolution raster (300+ DPI)
• PDF: Vector format for publications
• SVG: Scalable vector graphics (web/editing)
• EPS: Encapsulated PostScript (legacy journals)

INTERACTIVE FEATURES:
• Zoom: Box select to zoom in
• Pan: Drag to move view
• Reset: Double-click to reset view
• Hover: Show data point information
• Select: Click points to highlight

TIPS:
⭐ Use consistent color schemes across all plots
⭐ Add clear axis labels with units
⭐ Include informative titles
⭐ Export in vector format (PDF/SVG) for publications
⭐ Use 300+ DPI for high-quality PNG images
⭐ Save plot configurations for reproducibility
⭐ Test colorblind-friendly palettes
        """,
        
        "Pathway Annotation": """🛤️ PATHWAY ANNOTATION TAB

Annotate metabolites with pathway information and perform enrichment analysis.

OVERVIEW:
The Pathway Annotation tab maps metabolites to biological pathways using PathBank, SMPDB, and WikiPathways databases, then performs enrichment analysis to identify significantly affected pathways.

KEY FEATURES:
• Multiple pathway database integration
  - PathBank: Human/Rat/Mouse pathways
  - SMPDB: Small molecule pathways
  - WikiPathways: Community-curated pathways ✨ (organism-specific)
• Pathway enrichment analysis
• Z-score calculation for pathway ranking
• Network visualization with Cytoscape-style layouts
• Statistical significance testing
• Interactive pathway filtering
• Enzyme and transporter integration

HOW TO USE:
1. Load ID-annotated data (from ID Annotation tab)
2. Select organism (Human, Rat, or Mouse) ✨
3. Click "Run Pathway Annotation"
4. Wait for processing (progress shown)
5. Review metabolite-pathway mappings in table
6. Apply filters:
   - Metabolite filters: p-value, log2FC
   - Pathway filters: Z-score, min metabolites, max p-value
7. Select pathways of interest
8. Click "Generate Network" for visualization
9. Export results and networks

DATABASE REQUIREMENTS:
✅ pathbank_selected.feather (Human/Rat/Mouse pathways)
✅ merged_SMP_metabolites.feather (SMPDB metabolic pathways)
✅ wikipathways_homo_sapiens.feather (Human pathways) ✨ NEW!
✅ wikipathways_rattus_norvegicus.feather (Rat pathways) ✨ NEW!
✅ wikipathways_mus_musculus.feather (Mouse pathways) ✨ NEW!

Note: WikiPathways database automatically selected based on chosen organism

PATHWAY ENRICHMENT METRICS:
• Z-score: Pathway enrichment score
  Formula: (hits - expected) / sqrt(expected)
  - Positive Z-score: Pathway enriched
  - Z > 2.0: Statistically significant
  
• P-value: Hypergeometric test significance
  - p < 0.05: Significant enrichment
  
• Metabolite Coverage: % of pathway metabolites detected
• Fold Enrichment: Observed/Expected ratio

FILTERING OPTIONS:

Metabolite Filters:
• P-value threshold: 0.05 (default), 0.01, 0.001
• Log2FC threshold: ±0.5, ±1.0, ±1.5 (default: ±1.0)
• Direction: Up-regulated, Down-regulated, or Both

Pathway Filters:
• Z-score cutoff: 1.5, 2.0 (default), 2.5, 3.0
• Minimum metabolites: 2, 3 (default), 5, 10
• Max p-value: 0.05 (default), 0.01
• Min pathway coverage: 10%, 20%, 30%

NETWORK VISUALIZATION OPTIONS:
☑ Include enzymes (recommended)
☐ Include transporters
☐ Include associated diseases
☐ Show metabolite names as labels
☑ Color by fold change
☑ Size by -log10(p-value)

NETWORK LAYOUTS:
• Force-directed (default): Natural clustering
• Hierarchical: Top-down organization
• Circular: Radial arrangement
• Grid: Regular grid layout

NETWORK EXPORT:
• XGMML: Cytoscape import format
• GraphML: General graph format
• PNG/SVG: Static images
• Interactive HTML: Web-based viewer

OUTPUT FILES:
• pathway_annotation_results.xlsx: Complete results table
• enriched_pathways.csv: Filtered significant pathways
• pathway_network.xgmml: Cytoscape network file
• pathway_report.html: Interactive report

PATHWAY SOURCES IN OUTPUT:
• HMDB_Pathways: From HMDB database
• PathBank_Pathways: From PathBank (organism-filtered)
• SMP_Pathways: From SMPDB
• WikiPathways: From WikiPathways (organism-specific) ✨ NEW!
• All_Pathways: Combined from all sources

TIPS:
⭐ Process pathway databases in Database Setup tab first
⭐ Select correct organism to match your study species ✨
⭐ WikiPathways adds unique pathways not in other databases ✨
⭐ Apply metabolite filters BEFORE pathway analysis to reduce noise
⭐ Start with Z-score > 2.0, then relax to 1.5 if too few pathways
⭐ Include enzymes for comprehensive pathway networks
⭐ Use min metabolites = 3 to avoid spurious single-metabolite pathways
⭐ Export networks to Cytoscape for advanced visualization
⭐ Compare pathways across different conditions/groups
⭐ Document filter thresholds for reproducibility
⭐ For multi-species studies, process WikiPathways for all organisms ✨
        """,
        
        "Database Setup": """🗄️ DATABASE SETUP & CONFIGURATION

Configure and manage metabolite databases required for annotation and pathway analysis.

OVERVIEW:
The Database Setup tab handles downloading, processing, and verifying all required databases. Processing is required only ONCE - databases are saved as optimized Feather files for fast loading.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DATABASE CATALOG
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

ID ANNOTATION DATABASES (Left Column):

1. HMDB Database (Human Metabolome Database)
   Source: https://hmdb.ca/downloads
   File: hmdb_metabolites.xml (All Metabolites)
   Size: ~1.5 GB download, ~500 MB processed
   Content: 220,000+ metabolites
   Fields: HMDB ID, Name, Formula, MW, Synonyms, Pathways
   Processing Time: 5-15 minutes
   Button: "Process HMDB"
   Output: hmdb_database.feather
   
2. LipidMaps Database
   Source: https://www.lipidmaps.org/databases/lmsd/download
   File: structures.sdf (All Structures SDF file)
   Size: ~100 MB download, ~50 MB processed
   Content: 47,000+ lipids
   Fields: LM ID, Name, Formula, MW, Category, Class
   Processing Time: 2-5 minutes
   Button: "Process LipidMaps"
   Output: lipidmap.feather

PATHWAY ANNOTATION DATABASES (Right Column):

3. PathBank Database
   Source: https://pathbank.org/downloads
   File: pathbank_all_metabolites.csv
   Size: ~50 MB download, ~20 MB processed
   Organisms: Human, Rat, Mouse (auto-filtered)
   Content: 110,000+ pathway-metabolite associations
   Fields: Pathway ID, Name, Metabolite, HMDB ID, Species
   Processing Time: 1-3 minutes
   Button: "Process PathBank"
   Output: pathbank_selected.feather
   
4. SMPDB Database (Small Molecule Pathway Database)
   Source: https://smpdb.ca/downloads
   File: Metabolites CSV (ZIP archive)
   Size: ~200 MB download, ~80 MB processed
   Content: 48,687 pathway files → merged single database
   Fields: Pathway ID, Name, Metabolite, HMDB ID, KEGG ID
   Processing Time: 5-15 minutes (merges all files)
   Button: "Process SMPDB"
   Output: merged_SMP_metabolites.feather
   
   IMPORTANT: Extract the ZIP file first!
   You should have a folder with ~48,687 CSV files like:
   - SMP0000001_metabolites.csv
   - SMP0000002_metabolites.csv
   - ... (many more)
   
   Then select the FOLDER, not individual files.

5. WikiPathways Database (NEW!) ✨
   Source: https://zenodo.org/communities/wikipathways
   Files: GPML archives for each organism
   - wikipathways-*-gpml-Homo_sapiens.zip
   - wikipathways-*-gpml-Rattus_norvegicus.zip
   - wikipathways-*-gpml-Mus_musculus.zip
   Size: ~50-100 MB per organism (extracted)
   Organisms: Human, Rat, Mouse (organism-specific)
   Content: ~500 pathways (Human), ~200 each (Rat/Mouse)
   Fields: Metabolite, Pathways, Pathway IDs, HMDB, KEGG, ChEBI
   Processing Time: 2-5 minutes per organism
   Button: "Process WikiPathways"
   Output: 
   - wikipathways_homo_sapiens.feather
   - wikipathways_rattus_norvegicus.feather
   - wikipathways_mus_musculus.feather
   
   IMPORTANT: Extract each ZIP, then select folders!
   
   ORGANISM-SPECIFIC SELECTION:
   • Each organism has its own folder selector
   • Use dropdown to select correct organism-folder pairing
   • Can process one, two, or all three organisms
   • Prevents data mixing between species
   
   Workflow:
   1. Download GPML archive(s) from Zenodo
   2. Extract each ZIP to separate folders
   3. Click "Browse" next to organism name
   4. Select the extracted GPML folder
   5. Verify "X GPML files found" message
   6. Repeat for other organisms (if needed)
   7. Click "Process WikiPathways"
   8. Wait for completion (all selected organisms processed)
   
   Verification:
   • Status shows: "Configured: Homo sapiens, Rattus..."
   • File count displayed: "534 GPML files found"
   • Success message shows metabolite counts per organism

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SETUP INSTRUCTIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

STEP-BY-STEP GUIDE:

1. Download Raw Database File
   - Click the "Download from:" link below each database
   - Save to a convenient location
   - For SMPDB: Extract ZIP after download

2. Process Database
   - Click "Browse" button
   - Select the downloaded file (or folder for SMPDB)
   - Click "Process" button
   - Wait for completion (progress bar shown)
   
3. Verify Processing
   - Check Database Status log for green "✓ Available"
   - Verify file exists in program folder
   - Note: Processing runs in background (non-blocking)

4. Repeat for All Databases
   - Process databases as needed
   - Not all databases required for all features
   - See requirements below

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FEATURE REQUIREMENTS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Which databases do you need?

For ID Annotation:
✅ HMDB (required)
☐ LipidMaps (only if annotating lipids)

For Pathway Annotation:
✅ PathBank (required)
✅ SMPDB (required)
☑ WikiPathways (recommended) ✨ NEW!
  - Provides organism-specific pathway data
  - Adds unique pathways not in PathBank/SMPDB
  - Community-curated, regularly updated

Minimum Setup:
• HMDB + PathBank + SMPDB = Basic functionality

Recommended Setup:
• HMDB + PathBank + SMPDB + WikiPathways = Enhanced coverage

Full Setup:
• All databases = Complete functionality with maximum pathway coverage

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DATABASE STATUS LOG
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The Database Status section shows:
• Real-time availability of all databases
• Auto-updates after processing
• Download links for missing databases
• Last check timestamp

Status Indicators:
✓ Available (green): Database processed and ready
✗ Missing (red): Database not found
⚠ Warning (orange): Outdated or needs update
ℹ Info (blue): Processing in progress

Auto-Check:
• Status checked on application startup
• Runs in background (non-blocking)
• Progress bar shown during check
• Click "Refresh Status" to manually update

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
TROUBLESHOOTING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Common Issues:

1. Processing Fails
   ➜ Re-download source file (may be corrupted)
   ➜ Check file format matches requirements
   ➜ Ensure sufficient disk space (~2 GB free)
   ➜ Check console output for error details

2. Database Shows as Missing
   ➜ Click "Refresh Status" button
   ➜ Check file exists in program folder
   ➜ Verify .feather extension
   ➜ Re-process if needed

3. SMPDB Processing Slow
   ➜ Normal! Merging 48,687 files takes time
   ➜ Progress bar shows current file
   ➜ Don't close application during processing
   ➜ Background processing allows other work

4. WikiPathways - No GPML Files Found ✨ NEW!
   ➜ Verify you extracted the ZIP files
   ➜ Select the folder containing .gpml files (not ZIP)
   ➜ Check for files like WP####.gpml
   ➜ Download correct GPML archives (not GMT format)

5. WikiPathways - Wrong Organism Data ✨ NEW!
   ➜ Ensure correct folder selected for each organism
   ➜ Check organism name matches folder content
   ➜ Don't mix GPML files from different species
   ➜ Use organism-specific dropdowns correctly

6. Out of Memory Error
   ➜ Close other applications
   ➜ Process databases one at a time
   ➜ Increase system virtual memory
   ➜ Use 64-bit Python

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
TIPS & BEST PRACTICES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

⭐ Process databases only once - they're reusable
⭐ Processing runs in background - won't freeze GUI
⭐ Databases saved in same folder as program
⭐ Click "Refresh Status" to check current state
⭐ If processing fails, try re-downloading source
⭐ Keep source files for future re-processing
⭐ Feather files are optimized for fast loading
⭐ No internet required after initial download
⭐ Databases update periodically - check sources quarterly
⭐ WikiPathways: Use organism-specific folders to prevent data mixing ✨
⭐ WikiPathways: Download latest from Zenodo for most current pathways ✨
⭐ WikiPathways: Process all 3 organisms for multi-species studies ✨

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
FILE REFERENCE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Source Files → Processed Files:
• hmdb_metabolites.xml → hmdb_database.feather
• structures.sdf → lipidmap.feather
• pathbank_all_metabolites.csv → pathbank_selected.feather
• smpdb_metabolites.csv/ (folder) → merged_SMP_metabolites.feather
• wikipathways-*-gpml-Homo_sapiens/ → wikipathways_homo_sapiens.feather ✨
• wikipathways-*-gpml-Rattus_norvegicus/ → wikipathways_rattus_norvegicus.feather ✨
• wikipathways-*-gpml-Mus_musculus/ → wikipathways_mus_musculus.feather ✨

Cache Files (auto-generated):
• hmdb_cache.pkl
• lipidmaps_cache.pkl
• kegg_cache.pkl
• pubchem_cache.pkl

Configuration Files:
• id_annotation_prefs.json (ID annotation settings)
• database_structure_inspector.py (verification tool)
        """
    }


def display_help_menu():
    """Display interactive help menu in console."""
    help_content = get_help_content()
    
    print("\n" + "="*70)
    print("  METABOLITE ANNOTATION TOOL - HELP SYSTEM")
    print("="*70 + "\n")
    
    print("Available Help Topics:\n")
    topics = list(help_content.keys())
    for i, topic in enumerate(topics, 1):
        print(f"  {i}. {topic}")
    
    print("\n  0. Exit")
    print("\n" + "-"*70)
    
    while True:
        try:
            choice = input("\nSelect a topic (0-{}): ".format(len(topics)))
            
            if choice == '0':
                print("\nExiting help system. Goodbye!\n")
                break
            
            choice_num = int(choice)
            if 1 <= choice_num <= len(topics):
                topic = topics[choice_num - 1]
                print("\n" + "="*70)
                print(help_content[topic])
                print("="*70)
                input("\nPress Enter to continue...")
                print("\n" + "-"*70)
                print("Select another topic or 0 to exit:")
                for i, topic in enumerate(topics, 1):
                    print(f"  {i}. {topic}")
                print("  0. Exit")
                print("-"*70)
            else:
                print("Invalid choice. Please enter a number between 0 and {}.".format(len(topics)))
        
        except ValueError:
            print("Invalid input. Please enter a number.")
        except KeyboardInterrupt:
            print("\n\nExiting help system. Goodbye!\n")
            break


def main():
    """Main function for standalone execution."""
    import sys
    
    if len(sys.argv) > 1:
        # Command-line argument provided
        arg = sys.argv[1].lower()
        
        if arg in ['--help', '-h']:
            print(__doc__)
            return
        
        # Check if requesting specific tab
        if arg.startswith('--tab='):
            tab_name = arg.split('=')[1]
        elif arg == '--tab' and len(sys.argv) > 2:
            tab_name = sys.argv[2]
        else:
            tab_name = arg
        
        # Try to find matching tab
        help_content = get_help_content()
        matched = None
        
        for key in help_content.keys():
            if tab_name.lower() in key.lower().replace(' ', ''):
                matched = key
                break
        
        if matched:
            print("\n" + "="*70)
            print(help_content[matched])
            print("="*70 + "\n")
        else:
            print(f"\nError: Unknown tab '{tab_name}'")
            print("\nAvailable tabs:")
            for key in help_content.keys():
                print(f"  - {key}")
            print("\nUsage: python help.py --tab=<name>")
            print("       python help.py <name>\n")
    else:
        # No arguments - show interactive menu
        display_help_menu()


if __name__ == "__main__":
    main()
