import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext
import logging
import threading
import os
import sys
import time
import pandas as pd
import numpy as np
import subprocess
import platform
import re
from typing import Optional

# Import shared components
from gui.shared.base_tab import BaseTab, _setup_global_styles

# Import interactive network viewer (MetaboAnalyst-like)
try:
    from main_script.interactive_network_viewer import PathwayNetworkViewer, DASH_AVAILABLE
    INTERACTIVE_VIEWER_AVAILABLE = DASH_AVAILABLE
except ImportError:
    INTERACTIVE_VIEWER_AVAILABLE = False

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PathwayNetworkTab(BaseTab):
    CHECKED_SYMBOLS = {'☑️', '☑', '✓', '1', 'True'}
    """Pathway Network Tab - Pathway Network Visualization and Analysis"""
    
    def __init__(self, parent, data_manager):
        """Initialize Pathway Network tab"""
        super().__init__(parent, data_manager)
        
        # Setup global styles (runs only once)
        _setup_global_styles()
        
        # Get root window for dialogs
        self.root = parent.winfo_toplevel()
        
        # Initialize UI components (set in setup_ui, marked as Any for Pylance)
        self.pathway_progress_text: tk.scrolledtext.ScrolledText = None  # type: ignore
        self.pathway_stats_label: tk.Label = None  # type: ignore
        self.metabolite_stats_label: tk.Label = None  # type: ignore
        self.upstream_stats_label: tk.Label = None  # type: ignore
        self.diseases_stats_label: tk.Label = None  # type: ignore
        self.pathway_tree: ttk.Treeview = None  # type: ignore
        self.metabolite_tree: ttk.Treeview = None  # type: ignore
        self.upstream_tree: ttk.Treeview = None  # type: ignore
        self.disease_tree: ttk.Treeview = None  # type: ignore
        self.generate_network_button: tk.Button = None  # type: ignore
        self.open_folder_button: tk.Button = None  # type: ignore
        
        # Data storage
        self.pathway_original_metabolites_data = None
        self.pathway_filtered_metabolites_data = None
        self.pathway_original_pathways_data = None
        self.pathway_filtered_pathways_data = None
        self.pathway_selections = {}
        self.upstream_selections = {}
        self.disease_selections = {}
        self._last_viewer_selection = {'pathways': set(), 'upstream': set(), 'diseases': set()}
        self.verified_network_dataframe_raw = None
        self.verified_network_dataframe = None
        self.display_network_dataframe = None
        self.metabolite_class_lookup = {}
        self.class_compression_active = False
        self.class_level_mapping_table = None
        # CRITICAL: Store original uploaded raw data (never modified, always source for enrichment)
        self.original_uploaded_raw_data = None
        
        # Auto-connect checkboxes (default True = auto-connect enabled)
        self.pathway_auto_connect = tk.BooleanVar(value=True)
        self.upstream_auto_connect = tk.BooleanVar(value=True)
        self.disease_auto_connect = tk.BooleanVar(value=True)
        # Track whether we're rendering metabolite or lipid results
        self.current_data_mode: str = 'metabolite'
        self.assigned_columns: dict = {}
        self.feature_col: str = 'Name'
        self.pvalue_col: str = 'pvalue'
        self.log2fc_col: str = 'log2FC'
        self._last_multiomics_source_info = {
            'source_col': None,
            'n_metabolites_with_source': 0,
            'unique_sources': [],
        }
        
        # Remember last browsed directory
        self._last_browse_dir = None
        
        # Initialize UI state variables (will be set by setup_ui)
        self.network_log2fc_threshold = tk.DoubleVar(value=0.0)
        self.network_input_file = tk.StringVar(value='')
        self.verified_network_rename_map = {}
        self.pathway_rank_mode = tk.StringVar(value='Default')
        self.upstream_rank_mode = tk.StringVar(value='Default')
        self.disease_rank_mode = tk.StringVar(value='Default')
        self.pathway_organism_var = tk.StringVar(value='Homo sapiens')
        self.pathway_font_size_scale = tk.DoubleVar(value=1.0)
        self.chord_keep_metabolite_case = tk.BooleanVar(value=False)
        
        # Network visualization parameters (stored as instance variables)
        self.network_metabolite_font_size = tk.IntVar(value=10)
        self.network_pathway_font_size = tk.IntVar(value=9)
        self.network_enzyme_font_size = tk.IntVar(value=8)
        self.network_metabolite_max_chars = tk.IntVar(value=10)
        self.network_metabolite_max_lines = tk.IntVar(value=3)
        self.network_pathway_max_chars = tk.IntVar(value=20)
        self.network_pathway_max_lines = tk.IntVar(value=3)
        self.bargraph_pathway_max_chars = tk.IntVar(value=40)
        self.bargraph_pathway_max_lines = tk.IntVar(value=3)
        self.bargraph_title_font_size = tk.IntVar(value=22)
        self.bargraph_axis_title_font_size = tk.IntVar(value=20)
        self.bargraph_tick_font_size = tk.IntVar(value=15)
        self.bargraph_pathway_label_font_size = tk.IntVar(value=20)
        self.bargraph_width_px = tk.IntVar(value=1200)
        self.bargraph_height_px = tk.IntVar(value=700)
        self.bargraph_legend_title_font_size = tk.IntVar(value=14)
        self.bargraph_legend_tick_font_size = tk.IntVar(value=12)
        self.bargraph_legend_thickness = tk.IntVar(value=10)
        self.bargraph_legend_length_percent = tk.IntVar(value=70)
        self.chord_metabolite_max_chars = tk.IntVar(value=40)
        self.chord_metabolite_max_lines = tk.IntVar(value=3)
        self.chord_pathway_max_chars = tk.IntVar(value=40)
        self.chord_pathway_max_lines = tk.IntVar(value=3)
        self.chord_metabolite_font_size = tk.IntVar(value=12)
        self.chord_pathway_font_size = tk.IntVar(value=12)
        self.chord_figure_width_inches = tk.IntVar(value=12)
        self.chord_figure_height_inches = tk.IntVar(value=12)
        
        # Create UI
        self.setup_ui()
        print("[OK] Pathway Network Tab initialized")

    # ------------------------------
    # Data loading from annotation
    # ------------------------------
    def load_annotated_data(self, file_path, df=None):
        """
        Load pathway-annotated data from Pathway Annotation tab.
        
        This method is called automatically after annotation completes.
        
        Expected columns in annotated data:
        - Feature ID / Name
        - ID columns (HMDB, KEGG, etc.)
        - All_Pathways (semicolon-separated list)
        - Pathway columns from annotation
        - User's comparison columns (with p-values, log2FC, etc.)
        
        Args:
            file_path: Path to annotated Excel file
            df: Optional pre-loaded DataFrame
        """
        try:
            # Load data if not provided
            if df is None:
                df = pd.read_excel(file_path)

            # Strict source normalization: only an exact canonical 'Source' column is accepted.
            # No aliases or partial matches are promoted here.
            
            # Clear any previous enrichment data to force recalculation with new file
            self.pathway_filtered_pathways_data = None
            self.pathway_original_pathways_data = None
            self.verified_network_assignments = {}
            
            # Store data
            self.pathway_original_metabolites_data = df.copy()
            # CRITICAL: Store original uploaded raw data (used as source for all enrichment runs)
            self.original_uploaded_raw_data = df.copy()
            self.annotated_file_path = file_path
            
            # Show success FIRST
            filename = os.path.basename(file_path)
            self._log_network(f"✅ Annotated data loaded successfully!\n")
            self._log_network(f"📂 File: {filename}\n")
            self._log_network(f"📊 Rows: {len(df)}, Columns: {len(df.columns)}\n")
            
            # Auto-detect comparison columns
            self._log_network(f"\n🔍 Auto-detecting comparison columns...\n")
            self._auto_detect_network_columns(df)
            
            # Show info dialog
            messagebox.showinfo(
                "Data Loaded",
                f"✅ Annotated data loaded successfully!\n\n"
                f"File: {filename}\n"
                f"Rows: {len(df)}\n"
                f"Columns: {len(df.columns)}\n\n"
                f"Next: Click 'Verify Columns' to confirm/modify comparison columns"
            )
            
        except Exception as e:
            error_msg = f"Failed to load annotated data: {str(e)}"
            self._log_network(f"❌ {error_msg}\n")
            messagebox.showerror("Error", error_msg)
            logger.error(error_msg, exc_info=True)

    def import_annotated_data(self):
        """Import already-processed pathway annotation file"""
        try:
            # Use current output folder as initial directory if set, otherwise use last browsed location
            initial_dir = None
            if hasattr(self, 'network_output_folder') and self.network_output_folder.get():
                initial_dir = self.network_output_folder.get()
            elif hasattr(self, '_last_browse_dir') and self._last_browse_dir:
                initial_dir = self._last_browse_dir
            else:
                initial_dir = os.getcwd()
            
            # Browse for Excel file
            file_path = filedialog.askopenfilename(
                title="Select Pathway Annotated Data File",
                filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")],
                initialdir=initial_dir
            )
            
            if not file_path:
                return  # User cancelled
            
            # Remember the directory for next time
            self._last_browse_dir = os.path.dirname(file_path)
            
            # Also set output folder to the same directory if not already set
            if not self.network_output_folder.get():
                self.network_output_folder.set(self._last_browse_dir)
            
            # Load the file
            self.load_annotated_data(file_path)
            
        except Exception as e:
            error_msg = f"Failed to import annotated data: {str(e)}"
            messagebox.showerror("Error", error_msg)
            logger.error(error_msg, exc_info=True)
    
    def _auto_detect_network_columns(self, df: pd.DataFrame):
        """Auto-detect comparison columns (P-Value, Log2FC, Feature ID, Class, Class_name) for network analysis"""
        from gui.shared.column_assignment import ColumnDetector
        
        try:
            # Initialize assignments if not exists
            if not hasattr(self, 'verified_network_assignments'):
                self.verified_network_assignments = {}
            
            feature_id = None
            pvalue_candidates = []
            log2fc = None
            class_col = None
            class_name_col = None

            for col in df.columns:
                if feature_id is None and ColumnDetector.detect_column(col, 'Feature ID'):
                    feature_id = col
                if ColumnDetector.detect_column(col, 'P-Value'):
                    pvalue_candidates.append(col)
                if log2fc is None and ColumnDetector.detect_column(col, 'Log2 Fold Change'):
                    log2fc = col
                if class_col is None and ColumnDetector.detect_column(col, 'Class'):
                    class_col = col
                if class_name_col is None and ColumnDetector.detect_column(col, 'Class_name'):
                    class_name_col = col

            preferred_pvalue = ColumnDetector.select_best_pvalue_column(pvalue_candidates)

            if feature_id:
                self.verified_network_assignments['Feature ID'] = feature_id
                self._log_network(f"  ✓ Feature ID: {feature_id}\n")

            if preferred_pvalue:
                self.verified_network_assignments['P-Value'] = preferred_pvalue
                priority_note = "(adjusted)" if any(keyword in re.sub(r'[^a-z0-9]', '', preferred_pvalue.lower()) for keyword in ['adjp', 'padj', 'fdr', 'qvalue', 'adjustedp']) else ''
                priority_note = f" {priority_note}" if priority_note else ''
                self._log_network(f"  ✓ P-Value: {preferred_pvalue}{priority_note}\n")

            if log2fc:
                self.verified_network_assignments['Log2 Fold Change'] = log2fc
                self._log_network(f"  ✓ Log2 Fold Change: {log2fc}\n")
            
            if class_col:
                self.verified_network_assignments['Class'] = class_col
                self._log_network(f"  ✓ Class: {class_col}\n")
            
            if class_name_col:
                self.verified_network_assignments['Class_name'] = class_name_col
                self._log_network(f"  ✓ Class_name: {class_name_col}\n")
            
            # Check if all required columns detected
            required = ['Feature ID', 'P-Value', 'Log2 Fold Change']
            all_detected = all(col in self.verified_network_assignments for col in required)
            
            if all_detected:
                self._log_network(f"\n✅ All required columns auto-detected\n")
                self._log_network(f"💡 Click 'Verify Columns' to confirm or modify\n\n")
            else:
                missing = [col for col in required if col not in self.verified_network_assignments]
                self._log_network(f"\n⚠️ Could not auto-detect: {', '.join(missing)}\n")
                self._log_network(f"💡 Click 'Verify Columns' to manually assign\n\n")
            
        except Exception as e:
            self._log_network(f"⚠️ Auto-detection failed: {str(e)}\n")
            logger.error(f"Auto-detection error: {e}", exc_info=True)
    
    def verify_network_columns(self):
        """Verify and assign comparison columns for network analysis using dialog"""
        import threading
        from gui.shared.column_assignment import show_column_assignment_dialog
        
        compression_preference = self._is_class_compression_requested()

        def _load_and_verify():
            """Background thread worker for loading and showing dialog"""
            try:
                # Check if we have loaded data - try multiple sources
                df = None
                file_path = None
                
                # Priority 1: Check if data was loaded via load_annotated_data
                if hasattr(self, 'pathway_original_metabolites_data') and self.pathway_original_metabolites_data is not None:
                    df = self.pathway_original_metabolites_data
                    file_path = getattr(self, 'annotated_file_path', None)
                
                # Priority 2: Check data_manager
                if df is None:
                    df = self.data_manager.get_value('pathway_annotated_data')
                
                if df is None:
                    self.root.after(0, lambda: messagebox.showerror(
                        "No Data Loaded",
                        "Please load annotated pathway data first:\n\n"
                        "Option 1: Run annotation in 'Pathway Annotation' subtab\n"
                        "Option 2: Import pre-processed file using 'Browse & Import' button"
                    ))
                    return
                
                # Update UI
                self.root.after(0, lambda: self._log_network("📂 Analyzing columns for comparison selection...\n"))
                
                # Show column assignment dialog for pathway analysis
                result = show_column_assignment_dialog(
                    parent=self.root,
                    df=df,
                    tab_type='pathway',  # Requires Feature ID, P-Value, Log2FC
                    auto_calculate=False,
                    excel_file_path=file_path
                )
                
                if result:
                    # Store assignments
                    self.verified_network_assignments = result['assignments']
                    
                    # Rename columns for downstream compatibility
                    df_renamed = df.copy()
                    rename_map = {}
                    
                    for standard_name, user_column in result['assignments'].items():
                        if user_column and user_column in df_renamed.columns:
                            if standard_name == 'Feature ID':
                                rename_map[user_column] = 'Name'
                            elif standard_name == 'P-Value':
                                rename_map[user_column] = 'pvalue'
                            elif standard_name == 'Log2 Fold Change':
                                rename_map[user_column] = 'log2FC'
                            else:
                                rename_map[user_column] = standard_name.replace(' ', '_')
                    
                    # Store rename map for use in enrichment (maps original names to standard names)
                    self.verified_network_rename_map = rename_map
                    
                    if rename_map:
                        df_renamed.rename(columns=rename_map, inplace=True)
                        self.root.after(0, lambda: self._log_network("📝 Renamed columns for compatibility:\n"))
                        for old, new in rename_map.items():
                            self.root.after(0, lambda o=old, n=new: self._log_network(f"   {o} → {n}\n"))
                    
                    if df_renamed.columns.duplicated().any():
                        duplicate_cols = df_renamed.columns[df_renamed.columns.duplicated()].unique().tolist()
                        self.root.after(0, lambda: self._log_network(f"\n⚠️ Removing duplicate columns: {duplicate_cols}\n"))
                        df_renamed = df_renamed.loc[:, ~df_renamed.columns.duplicated(keep='last')]
                        self.root.after(0, lambda: self._log_network(f"✓ Duplicates removed, {len(df_renamed.columns)} columns remain\n\n"))

                    # Preserve only an exact canonical 'Source' column if it already exists.
                    if 'Source' not in df_renamed.columns and 'Source' in df.columns and len(df) == len(df_renamed):
                        try:
                            df_renamed['Source'] = df['Source'].values
                        except Exception:
                            pass
                    
                    # ═══════════════════════════════════════════════════════════════
                    # WORKFLOW BRANCH LOGIC
                    # ═══════════════════════════════════════════════════════════════
                    compress_requested = compression_preference
                    
                    # Three explicit workflow branches:
                    # 1. Compression OFF → individual nodes only
                    # 2. Compression ON + no Class_name column → fallback to individual nodes  
                    # 3. Compression ON + Class_name exists → aggregate classes
                    
                    df_for_network = None
                    class_mapping_df = None
                    
                    if not compress_requested:
                        # ─────────────────────────────────────────────────────────
                        # BRANCH 1: NORMAL MODE (Compression checkbox UNCHECKED)
                        # ─────────────────────────────────────────────────────────
                        self._log_network("\n🔹 WORKFLOW MODE: Individual Compounds (No Compression)\n")
                        self._log_network("   Using all individual lipid/metabolite nodes\n")
                        self._log_network("   Statistics will be calculated per compound\n\n")
                        
                        df_for_network = df_renamed.copy()
                        class_mapping_df = None
                        
                    elif compress_requested and 'Class_name' not in df_renamed.columns:
                        # ─────────────────────────────────────────────────────────
                        # BRANCH 2b: COMPRESSION REQUESTED BUT NO Class_name
                        # Fallback to individual mode (same as Branch 1)
                        # ─────────────────────────────────────────────────────────
                        self._log_network("\n⚠️ WORKFLOW MODE: Individual Compounds (Class_name column missing)\n")
                        self._log_network("   Compression requested but Class_name not found\n")
                        self._log_network("   Falling back to individual compound mode\n\n")
                        
                        df_for_network = df_renamed.copy()
                        class_mapping_df = None
                        
                    elif compress_requested and 'Class_name' in df_renamed.columns:
                        try:
                            # Helper to split pathways consistently
                            def _parse_pathways(value):
                                if value is None or (isinstance(value, float) and pd.isna(value)):
                                    return []
                                text = str(value).strip()
                                if not text or text.lower() == 'nan':
                                    return []
                                cleaned = text.replace(';', '|')
                                parts = [p.strip() for p in cleaned.split('|') if p.strip()]
                                return parts

                            # Step 1: Get user thresholds for significance
                            try:
                                p_threshold = float(self.network_p_threshold.get()) if hasattr(self, 'network_p_threshold') else 0.05
                            except (ValueError, AttributeError):
                                p_threshold = 0.05

                            try:
                                log2fc_threshold = float(self.network_log2fc_threshold.get()) if hasattr(self, 'network_log2fc_threshold') else 0.0
                            except (ValueError, AttributeError):
                                log2fc_threshold = 0.0

                            # Step 2: Filter to significant lipids (stats already calculated)
                            if 'pvalue' in df_renamed.columns and 'log2FC' in df_renamed.columns:
                                significant_mask = (
                                    (pd.to_numeric(df_renamed['pvalue'], errors='coerce') <= p_threshold) &
                                    (pd.to_numeric(df_renamed['log2FC'], errors='coerce').abs() >= log2fc_threshold)
                                )
                                df_significant = df_renamed[significant_mask].copy()
                                n_significant = len(df_significant)

                                self._log_network(f"🎯 Applying significance thresholds: p≤{p_threshold}, |log2FC|≥{log2fc_threshold}\n")
                                self._log_network(f"   • Total lipids: {len(df_renamed)}\n")
                                self._log_network(f"   • Significant lipids: {n_significant}\n\n")

                                if n_significant == 0:
                                    self._log_network("⚠️ No significant lipids found with current thresholds. Using all lipids.\n")
                                    df_significant = df_renamed.copy()
                            else:
                                self._log_network("⚠️ P-value or log2FC column missing. Using all lipids for compression.\n")
                                df_significant = df_renamed.copy()

                            # Step 3: Separate lipids with/without Class_name
                            class_name_values = df_significant['Class_name'].astype(str).str.strip()
                            has_class_mask = (
                                df_significant['Class_name'].notna() &
                                (class_name_values != '') &
                                (class_name_values.str.lower() != 'nan') &
                                (class_name_values.str.lower() != 'na')
                            )

                            df_with_class = df_significant[has_class_mask].copy()
                            df_without_class = df_significant[~has_class_mask].copy()

                            self._log_network("📊 Class assignment analysis:\n")
                            self._log_network(f"   • Significant lipids WITH Class_name: {len(df_with_class)}\n")
                            self._log_network(f"   • Significant lipids WITHOUT Class_name: {len(df_without_class)}\n")

                            # Step 4: Aggregate lipids WITH Class_name into class nodes while preserving unique members
                            preserved_individual_rows = []
                            preserved_details = []
                            unique_pathway_ratio = getattr(self, 'class_unique_pathway_ratio', 0.5)

                            if len(df_with_class) > 0:
                                class_rows = []
                                mapping_rows = []
                                numeric_cols = df_with_class.select_dtypes(include=[np.number]).columns.tolist()
                                pathway_cols = ['All_Pathways', 'All_Pathways_Display', 'HMDB_Pathways',
                                               'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways',
                                               'Metabolika_Pathways', 'Reactome_Pathways', 'KEGG_Pathways']

                                for class_name, block in df_with_class.groupby('Class_name'):
                                    block = block.copy()
                                    class_size = len(block)
                                    member_names = block['Name'].astype(str).tolist()

                                    pathway_col = next((c for c in pathway_cols if c in block.columns), None)
                                    member_pathways = {}
                                    pathway_counts = {}
                                    if pathway_col:
                                        for _, brow in block.iterrows():
                                            member = str(brow['Name']).strip()
                                            p_list = list(set(_parse_pathways(brow.get(pathway_col, ''))))
                                            member_pathways[member] = p_list
                                            for p in p_list:
                                                pathway_counts[p] = pathway_counts.get(p, 0) + 1

                                    members_to_preserve = set()
                                    preservation_reasons = {}
                                    if pathway_counts:
                                        for member, p_list in member_pathways.items():
                                            unique_paths = []
                                            for pathway in p_list:
                                                share = pathway_counts.get(pathway, 0) / class_size
                                                if share < unique_pathway_ratio:
                                                    unique_paths.append(pathway)
                                            if unique_paths:
                                                members_to_preserve.add(member)
                                                preservation_reasons[member] = unique_paths

                                    if members_to_preserve:
                                        preserved_rows = block[block['Name'].isin(members_to_preserve)].copy()
                                        preserved_individual_rows.append(preserved_rows)
                                        preserved_details.append((class_name, preservation_reasons))
                                        block = block[~block['Name'].isin(members_to_preserve)].copy()
                                        member_names = block['Name'].astype(str).tolist()
                                        class_size = len(block)

                                    if class_size == 0:
                                        continue  # Nothing left to aggregate

                                    agg_row = block.iloc[0].copy()
                                    # Preserve original class names unless they collide with a member metabolite name.
                                    # Use an explicit marker instead of pluralization to keep names semantically correct.
                                    class_display_name = str(class_name).strip()
                                    member_name_set = {str(m).strip().lower() for m in member_names if str(m).strip()}
                                    if class_display_name.lower() in member_name_set:
                                        class_display_name = f"{class_display_name} (class)"
                                    agg_row['Name'] = class_display_name
                                    agg_row['Class_name'] = class_display_name
                                    agg_row['Class_members'] = '|'.join(member_names)
                                    agg_row['Class_member_count'] = len(member_names)

                                    for col in numeric_cols:
                                        agg_row[col] = pd.to_numeric(block[col], errors='coerce').mean(skipna=True)

                                    for pcol in pathway_cols:
                                        if pcol in block.columns:
                                            all_pathways = set()
                                            for val in block[pcol].dropna():
                                                all_pathways.update(_parse_pathways(val))
                                            if all_pathways:
                                                agg_row[pcol] = '; '.join(sorted(all_pathways))

                                    disease_cols = ['Associated_Diseases', 'Reactome_Disease']
                                    for dcol in disease_cols:
                                        if dcol in block.columns:
                                            all_diseases = set()
                                            for val in block[dcol].dropna():
                                                all_diseases.update(_parse_pathways(val))
                                            if all_diseases:
                                                agg_row[dcol] = '|'.join(sorted(all_diseases))

                                    upstream_cols = ['Enzyme_Gene_name', 'Transporter_Gene_name']
                                    for ucol in upstream_cols:
                                        if ucol in block.columns:
                                            all_upstream = set()
                                            for val in block[ucol].dropna():
                                                tokens = [t.strip() for t in str(val).split('|') if t.strip()]
                                                all_upstream.update(tokens)
                                            if all_upstream:
                                                agg_row[ucol] = '|'.join(sorted(all_upstream))

                                    class_rows.append(agg_row)

                                    for member_name in member_names:
                                        mapping_rows.append({
                                            'Class': class_display_name,
                                            'Member': member_name,
                                            'Member_Count': len(member_names)
                                        })

                                df_class_nodes = pd.DataFrame(class_rows)
                                if df_class_nodes.empty:
                                    df_class_nodes = pd.DataFrame(columns=df_renamed.columns)
                                class_mapping_df = pd.DataFrame(mapping_rows) if mapping_rows else None

                                preserved_df = pd.concat(preserved_individual_rows, ignore_index=True) if preserved_individual_rows else pd.DataFrame(columns=df_renamed.columns)

                                df_for_network = pd.concat([df_class_nodes, preserved_df, df_without_class], ignore_index=True)

                                n_class_nodes = len(df_class_nodes)
                                n_preserved = len(preserved_df)
                                n_without_class = len(df_without_class)

                                self._log_network("\n✅ Class compression complete:\n")
                                self._log_network(f"   • Class nodes created: {n_class_nodes}\n")
                                self._log_network(f"   • Unique lipids preserved from their classes: {n_preserved}\n")
                                self._log_network(f"   • Lipids without Class_name: {n_without_class}\n")
                                self._log_network(f"   • Total compounds for analysis: {len(df_for_network)}\n")

                                if preserved_details:
                                    self._log_network("\n🔍 Preserved individual lipids (unique pathways <50% coverage):\n")
                                    for class_name, reasons in preserved_details:
                                        for member, paths in reasons.items():
                                            preview = ', '.join(paths[:5]) if paths else 'unique signal'
                                            # Show original class name in logs for clarity
                                            self._log_network(f"   • {member} (class {class_name}) → {preview}\n")
                                    self._log_network("\n")
                            else:
                                df_for_network = df_without_class.copy()
                                self._log_network("ℹ️ No significant lipids with Class_name found. Using individual lipids.\n")

                        except Exception as compression_error:
                            self._log_network(
                                f"⚠️ Class compression failed: {compression_error}\n"
                                "   Falling back to individual lipid nodes.\n"
                            )
                            df_for_network = df_renamed.copy()
                            class_mapping_df = None
                            import traceback
                            traceback.print_exc()
                    
                    # Verify df_for_network was set
                    if df_for_network is None:
                        raise ValueError("Internal error: df_for_network was not set by any workflow branch")
                    
                    # ═══════════════════════════════════════════════════════════════
                    # Store results and set compression state flags
                    # ═══════════════════════════════════════════════════════════════
                    self.class_level_mapping_table = class_mapping_df
                    if class_mapping_df is not None:
                        self.data_manager.set_value('class_level_mapping_table', class_mapping_df)
                    else:
                        self.data_manager.set_value('class_level_mapping_table', None)

                    self.verified_network_dataframe_raw = df_renamed
                    self.verified_network_dataframe = df_for_network
                    self.display_network_dataframe = df_for_network
                    # CRITICAL: Always update original raw data with renamed columns when verifying
                    # This ensures enrichment always works with properly named columns
                    self.original_uploaded_raw_data = df_renamed.copy()
                    
                    # Build metabolite→class lookup from mapping table
                    self.metabolite_class_lookup = {}
                    if class_mapping_df is not None and not class_mapping_df.empty:
                        for _, row in class_mapping_df.iterrows():
                            member = row.get('Member', '')
                            class_name = row.get('Class', '')
                            if member and class_name:
                                self.metabolite_class_lookup[str(member)] = str(class_name)
                    # Set compression active flag based on whether mapping was successfully created
                    self.class_compression_active = bool(
                        compress_requested and 
                        class_mapping_df is not None and 
                        not class_mapping_df.empty and
                        len(self.metabolite_class_lookup) > 0
                    )
                    
                    # Explicit state cleanup when compression is inactive
                    if not self.class_compression_active:
                        self.metabolite_class_lookup = {}
                        self.class_compression_active = False
                        self._log_network("✓ Compression state: INACTIVE (individual compound mode)\n\n")
                    else:
                        self._log_network(
                            "✓ Compression state: ACTIVE\n"
                            "   Statistics: calculated on individuals → averaged for classes\n"
                            "   Unique lipids with <50% pathway overlap preserved\n\n"
                        )

                    self.data_manager.set_value('pathway_annotated_data', df_for_network)
                    
                    # Get assignments for summary
                    feature_id = result['assignments'].get('Feature ID')
                    pvalue = result['assignments'].get('P-Value')
                    log2fc = result['assignments'].get('Log2 Fold Change')
                    
                    summary_lines = ["Network Analysis Column Verification:"]
                    summary_lines.append(f"  ✓ Feature ID: {feature_id} → Name")
                    summary_lines.append(f"  ✓ P-Value: {pvalue} → pvalue")
                    summary_lines.append(f"  ✓ Log2 Fold Change: {log2fc} → log2FC")
                    summary_text = "\n".join(summary_lines)
                    
                    def enable_generate():
                        if hasattr(self, 'generate_network_button'):
                            self.generate_network_button.config(state='disabled', bg='#95a5a6')
                        self._log_network(f"✅ Column verification completed\n{summary_text}\n")
                        self._log_network("\n🔄 Running pathway enrichment analysis...\n\n")
                        messagebox.showinfo(
                            "✅ Column Verification Complete",
                            summary_text + "\n\nRunning pathway enrichment analysis..."
                        )
                        self._log_network("🔄 Forcing fresh calculation with current settings...\n\n")
                        self.pathway_filtered_pathways_data = None
                        self.frame.after(100, self.run_pathway_enrichment)
                    
                    self.root.after(0, enable_generate)
                else:
                    self.root.after(0, lambda: self._log_network("❌ Column verification cancelled\n"))
                    
            except Exception as e:
                self.root.after(0, lambda: (
                    messagebox.showerror("Error", f"Failed to verify columns: {str(e)}"),
                ))
                logger.error(f"Network column verification error: {e}", exc_info=True)
        
        thread = threading.Thread(target=_load_and_verify, daemon=True)
        thread.start()

    def _build_class_lookup_from_mapping(self, mapping_df):
        """Create metabolite→class mapping for remapping hits after compression."""
        if mapping_df is None or mapping_df.empty:
            return {}

        lookup = {}
        for _, row in mapping_df.iterrows():
            class_name = str(row.get('Class', '')).strip()
            members = str(row.get('Member_Lipids', '')).split('|') if 'Member_Lipids' in row.index else []
            for member in members:
                member_clean = member.strip()
                if not member_clean or member_clean.lower() == 'nan':
                    continue
                # Preserve first mapping to avoid oscillations when lipids appear in many pathways
                lookup.setdefault(member_clean, class_name or member_clean)
        return lookup

    def _map_pathway_hits_to_display_nodes(self, pathways_dict):
        """Replace metabolite names in pathway stats with class nodes when compressed."""
        if not self.class_compression_active:
            return
        lookup = self.metabolite_class_lookup or {}
        if not lookup:
            return

        for pathway_name, stats in pathways_dict.items():
            if not isinstance(stats, dict):
                continue
            for key, count_key in (('hits', 'k_hits'), ('metabolites', 'n_metabolites')):
                values = stats.get(key)
                if not isinstance(values, list):
                    continue
                original_values = [str(v).strip() for v in values if str(v).strip()]
                mapped = []
                seen = set()
                for val in original_values:
                    display_val = lookup.get(val, val)
                    if not display_val:
                        continue
                    if display_val not in seen:
                        mapped.append(display_val)
                        seen.add(display_val)
                if key == 'hits' and 'raw_hits' not in stats:
                    stats['raw_hits'] = original_values
                if key == 'metabolites' and 'raw_metabolites' not in stats:
                    stats['raw_metabolites'] = original_values
                stats[key] = mapped
                stats[count_key] = len(mapped)

    # ------------------------------
    # Tooltip helper
    # ------------------------------
    def _create_tooltip(self, widget, text):
        """Create a tooltip for a widget that shows on hover."""
        tooltip_window = None
        
        def on_enter(event):
            nonlocal tooltip_window
            try:
                x = widget.winfo_rootx() + 25
                y = widget.winfo_rooty() + 25
            except Exception:
                x, y = 0, 0
            
            tooltip_window = tk.Toplevel(widget)
            tooltip_window.wm_overrideredirect(True)
            tooltip_window.wm_geometry(f"+{x}+{y}")
            
            label = tk.Label(tooltip_window, text=text, justify='left',
                           background="#ffffe0", relief='solid', borderwidth=1,
                           font=("Arial", 8, "normal"), padx=6, pady=4)
            label.pack()
        
        def on_leave(event):
            nonlocal tooltip_window
            if tooltip_window:
                tooltip_window.destroy()
                tooltip_window = None
        
        widget.bind("<Enter>", on_enter)
        widget.bind("<Leave>", on_leave)

    def _get_shared_metabolite_pvalue_var(self):
        """Return the StringVar backing all metabolite p-value controls."""
        if not hasattr(self, '_shared_metabolite_pvalue_var'):
            self._shared_metabolite_pvalue_var = tk.StringVar(value="0.05")
        # Keep both legacy attribute names pointing to the same variable
        self.metabolite_pvalue_filter = self._shared_metabolite_pvalue_var
        self.component_pvalue_threshold = self._shared_metabolite_pvalue_var
        return self._shared_metabolite_pvalue_var

    def _update_network_suboptions(self):
        """Enable/disable network sub-options based on Generate Network Graph checkbox state."""
        gen_network = self.network_gen_network.get()
        state = 'normal' if gen_network else 'disabled'
        
        # Update sub-option checkbuttons
        if hasattr(self, 'auto_launch_cb'):
            self.auto_launch_cb.config(state=state)
        if hasattr(self, 'export_cytoscape_cb'):
            self.export_cytoscape_cb.config(state=state)

    # ------------------------------
    # Mode helpers
    # ------------------------------
    def _is_lipid_mode(self) -> bool:
        """Return True when the incoming dataset represents lipids."""
        return str(getattr(self, 'current_data_mode', 'metabolite')).lower() == 'lipid'

    def _is_class_compression_requested(self) -> bool:
        """Safely read the Compress by Class checkbox (handles stringy Tk values)."""
        var = getattr(self, 'compress_by_class', None)
        if var is None:
            return False
        try:
            value = var.get()
        except Exception:
            return False

        if isinstance(value, str):
            lowered = value.strip().lower()
            if lowered in ('', '0', 'false', 'off', 'no'):
                return False
            if lowered in ('1', 'true', 'on', 'yes'):
                return True
            try:
                return bool(int(lowered))
            except ValueError:
                return False
        return bool(value)

    def _component_label(self, plural: bool = True) -> str:
        """Human-readable label for the active component type."""
        if self._is_lipid_mode():
            return 'Lipids' if plural else 'Lipid'
        return 'Compounds' if plural else 'Compound'

    def _update_component_labels(self):
        """Refresh UI labels/headings that differ between metabolite and lipid data."""
        label = self._component_label(plural=True)

        # Notebook tab label
        tab = getattr(self, 'metabolite_tab', None)
        notebook = getattr(self, 'table_notebook', None)
        if tab is not None and notebook is not None:
            try:
                notebook.tab(tab, text=f"🧪 {label}")
            except Exception:
                pass

        # Stats placeholder when no data available
        stats_label = getattr(self, 'metabolite_stats_label', None)
        if stats_label is not None:
            try:
                if (self.pathway_filtered_metabolites_data is None or
                        len(self.pathway_filtered_metabolites_data) == 0):
                    stats_label.config(text=f"No {label.lower()} data loaded")
            except Exception:
                pass

        # Tree heading for identifier column
        tree = getattr(self, 'metabolite_tree', None)
        if tree is not None:
            heading = 'Lipid ID' if self._is_lipid_mode() else 'Name'
            try:
                tree.heading('Name', text=heading)
            except Exception:
                pass

        # Update pathway table column heading for # Metabolites/Lipids
        pathway_tree = getattr(self, 'pathway_tree', None)
        if pathway_tree is not None:
            component_count_label = f"# {label}"
            try:
                pathway_tree.heading('# Compounds', text=component_count_label)
            except Exception:
                pass

    # ------------------------------
    # Column assignment helpers
    # ------------------------------
    def _resolve_assigned_column(self, logical_key, df=None, fallbacks=None):
        """Return the column name for a logical key using verified assignments with fallbacks."""
        df_columns = set(df.columns) if df is not None else set()
        assigned = None
        if getattr(self, 'assigned_columns', None):
            assigned = self.assigned_columns.get(logical_key)
            if assigned and (not df_columns or assigned in df_columns):
                return assigned
        if fallbacks:
            for candidate in fallbacks:
                if not df_columns or candidate in df_columns:
                    return candidate
            # Last resort: return first fallback even if missing so downstream warnings can trigger
            return fallbacks[0] if isinstance(fallbacks, (list, tuple)) and fallbacks else fallbacks
        return assigned

    def _cache_assigned_columns(self, df=None):
        """Derive cached column names for feature id, p-value, and log2FC.
        
        IMPORTANT: If DataFrame has standardized columns (Name, pvalue, log2FC) from
        annotation processing, use those. Otherwise detect from verified assignments.
        """
        # First check if we already have standardized columns in the DataFrame
        # (created by annotation tab's column mapping layer)
        if df is not None and not df.empty:
            has_standard_name = 'Name' in df.columns
            has_standard_pval = 'pvalue' in df.columns
            has_standard_log2 = 'log2FC' in df.columns
            
            if has_standard_name and has_standard_pval and has_standard_log2:
                # Use standardized columns directly
                self.feature_col = 'Name'
                self.pvalue_col = 'pvalue'
                self.log2fc_col = 'log2FC'
                logger.info(f"[ASSIGNMENTS] Using standardized columns: Name, pvalue, log2FC")
                return
        
        # Otherwise, detect from verified assignments with fallbacks
        fallback_name = ['Name', 'LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite']
        fallback_pval = ['pvalue', 'p_value', 'P-Value', 'p-value', 'adj_p']
        fallback_log2 = ['log2FC', 'log2_fc', 'log2_fold_change', 'log2(Fold Change)', 'FC']

        self.feature_col = self._resolve_assigned_column('Feature ID', df, fallback_name) or fallback_name[0]  # type: ignore
        self.pvalue_col = self._resolve_assigned_column('P-Value', df, fallback_pval) or fallback_pval[0]  # type: ignore
        self.log2fc_col = self._resolve_assigned_column('Log2 Fold Change', df, fallback_log2) or fallback_log2[0]  # type: ignore

        logger.info(f"[ASSIGNMENTS] Feature column: {self.feature_col}")
        logger.info(f"[ASSIGNMENTS] P-value column: {self.pvalue_col}")
        logger.info(f"[ASSIGNMENTS] Log2FC column: {self.log2fc_col}")

    def _get_metabolite_value(self, row, logical_key, default=None):
        """Fetch a metabolite value using cached assignments with graceful fallbacks."""
        column_map = {
            'name': self.feature_col,
            'pvalue': self.pvalue_col,
            'log2fc': self.log2fc_col,
        }
        fallback_map = {
            'name': ['Name', 'Compound_Name', 'Compound Name', 'Feature_ID', 'feature_id'],
            'pvalue': ['pvalue', 'P-Value', 'p_value', 'p-value', 'adj_p', 'adj.p'],
            'log2fc': ['log2FC', 'log2_fc', 'log2 fold change', 'log2_fold_change', 'log2(Fold Change)']
        }

        col = column_map.get(logical_key)
        if col and col in row.index:
            val = row.get(col)
            if pd.notna(val):
                return val
        for candidate in fallback_map.get(logical_key, []):
            if candidate in row.index:
                val = row.get(candidate)
                if pd.notna(val):
                    return val
        return default

    def _as_float(self, value, default=0.0):
        """Convert arbitrary values to float, honoring NaNs and blanks."""
        try:
            if value is None:
                return default
            if isinstance(value, str):
                if not value.strip():
                    return default
            if pd.isna(value):
                return default
            return float(value)
        except (TypeError, ValueError):
            return default

    # ------------------------------
    # Logging helpers
    # ------------------------------
    def _silence_third_party_logs(self, enable: bool = True):
        """Temporarily silence very verbose third-party logs (py4cytoscape, kaleido, etc.).

        When enable=True, stores current levels and sets WARNING for known noisy loggers.
        When enable=False, restores the stored levels if present.
        """
        import logging
        if not hasattr(self, "_logger_levels_backup"):
            self._logger_levels_backup = {}
        noisy_loggers = [
            'py4cytoscape', 'py4cytoscape.cyrest', 'py4cytoscape.tuning',
            'kaleido', 'kaleido.scopes.plotly',
            'urllib3', 'PIL', 'matplotlib',
        ]
        if enable:
            for name in noisy_loggers:
                lg = logging.getLogger(name)
                # backup only once
                if name not in self._logger_levels_backup:
                    self._logger_levels_backup[name] = lg.level
                lg.setLevel(logging.WARNING)
            # Also reduce root INFO noise to WARNING without affecting our prints
            self._logger_levels_backup['root'] = logging.getLogger().level
            logging.getLogger().setLevel(logging.WARNING)
            # Globally disable DEBUG/INFO logs from any library
            self._logger_levels_backup['disable_level'] = getattr(logging, 'current_disable_level', logging.NOTSET)
            logging.disable(logging.INFO)
        else:
            # restore
            for name, lvl in self._logger_levels_backup.items():
                if name == 'root':
                    logging.getLogger().setLevel(lvl)
                else:
                    logging.getLogger(name).setLevel(lvl)
            # Re-enable logging globally
            logging.disable(logging.NOTSET)
            self._logger_levels_backup = {}
    
    def set_pathway_data(self, pathways_data, metabolites_data, z_threshold=None, min_metabolites=None,
                         pvalue_threshold=None, data_mode=None, column_assignments=None, output_folder=None):
        """
        Receive pathway and metabolite data from Pathway Annotation tab.
        This method is called when annotation completes to populate the network tab.
        
        Parameters:
        -----------
        pathways_data : dict
            Dictionary of pathway statistics from annotation
        metabolites_data : pd.DataFrame
            DataFrame of annotated metabolites
        z_threshold : float, optional
            Z-score threshold used in pathway annotation (if None, keeps current value)
        min_metabolites : int, optional
            Minimum metabolites per pathway from annotation (if None, keeps current value)
        pvalue_threshold : float, optional
            P-value threshold from annotation (if None, keeps current value)
        data_mode : str, optional
            Either 'metabolite' or 'lipid' so UI labels and exports match the dataset type
        column_assignments : dict, optional
            Mapping of verified column roles (e.g., Feature ID, P-Value, Log2 Fold Change)
        output_folder : str, optional
            Output folder path from annotation tab to use as default for network outputs
        """
        try:
            mode_value = str(data_mode).lower() if data_mode else 'metabolite'
            self.current_data_mode = 'lipid' if mode_value == 'lipid' else 'metabolite'
            self.assigned_columns = (column_assignments or {}).copy()
            component_label = self._component_label(plural=True)
            
            # Set output folder from annotation if provided and not already set
            if output_folder and os.path.isdir(output_folder):
                if not self.network_output_folder.get():
                    self.network_output_folder.set(output_folder)
                    logger.info(f"📁 Output folder set from annotation: {output_folder}")
                self._last_browse_dir = output_folder
            
            logger.info(f"📊 Receiving pathway data from annotation tab...")
            logger.info(f"   Pathways: {len([k for k in pathways_data.keys() if not k.startswith('_')])}")
            logger.info(f"   {component_label}: {len(metabolites_data) if metabolites_data is not None else 0}")
            
            # DEBUG: Check disease data on receipt
            if metabolites_data is not None:
                logger.info(f"[RECEIVE DEBUG] {component_label} shape: {metabolites_data.shape}")
                logger.info(f"[RECEIVE DEBUG] Has 'Associated_Diseases': {'Associated_Diseases' in metabolites_data.columns}")
                logger.info(f"[RECEIVE DEBUG] Has 'Associated': {'Associated' in metabolites_data.columns}")
                logger.info(f"[RECEIVE DEBUG] Has 'Reactome_Disease': {'Reactome_Disease' in metabolites_data.columns}")

                disease_col = None
                if 'Associated_Diseases' in metabolites_data.columns:
                    disease_col = 'Associated_Diseases'
                elif 'Associated' in metabolites_data.columns:
                    disease_col = 'Associated'

                if disease_col:
                    non_empty = metabolites_data[disease_col].fillna('').astype(str).str.strip()
                    non_empty = non_empty[non_empty != '']
                    non_empty = non_empty[non_empty.str.lower() != 'nan']
                    logger.info(f"[RECEIVE DEBUG] Non-empty {disease_col}: {len(non_empty)}/{len(metabolites_data)}")
                    if len(non_empty) > 0:
                        logger.info(f"[RECEIVE DEBUG] Sample (first 2): {list(non_empty.head(2))}")
            
            # CRITICAL: Update ALL filter settings to match annotation settings
            if z_threshold is not None:
                logger.info(f"   Z-score threshold from annotation: {z_threshold}")
                self.network_zscore_activation.set(str(z_threshold))
                self.network_zscore_inhibition.set(str(-z_threshold))
                logger.info(f"✓ Network z-score thresholds updated: activation={z_threshold}, inhibition={-z_threshold}")
            
            if min_metabolites is not None:
                logger.info(f"   Min metabolites from annotation: {min_metabolites}")
                logger.info(f"   Current network min metabolites BEFORE update: {self.network_min_metabolites.get()}")
                self.network_min_metabolites.set(str(min_metabolites))
                logger.info(f"   Current network min metabolites AFTER update: {self.network_min_metabolites.get()}")
                logger.info(f"✓ Network min metabolites updated: {min_metabolites}")
            
            if pvalue_threshold is not None:
                logger.info(f"   P-value threshold from annotation: {pvalue_threshold}")
                self.network_p_threshold.set(str(pvalue_threshold))
                logger.info(f"✓ Network p-value threshold updated: {pvalue_threshold}")
            
            # DEBUG: Track pathway count before cleaning
            pathways_before_clean = len([k for k in (pathways_data or {}).keys() if not k.startswith('_')])
            logger.info(f"[CLEAN DEBUG] Pathways before _clean_and_merge_pathways_data: {pathways_before_clean}")
            
            # Store the data (CLEAN like old GUI before using)
            cleaned = self._clean_and_merge_pathways_data(pathways_data or {})
            
            # DEBUG: Track pathway count after cleaning
            pathways_after_clean = len([k for k in cleaned.keys() if not k.startswith('_')])
            logger.info(f"[CLEAN DEBUG] Pathways after _clean_and_merge_pathways_data: {pathways_after_clean}")
            if pathways_after_clean < pathways_before_clean:
                removed = pathways_before_clean - pathways_after_clean
                logger.warning(f"[CLEAN DEBUG] {removed} pathway(s) removed by cleaning function")
            
            self.pathway_original_pathways_data = cleaned.copy()
            self.pathway_filtered_pathways_data = cleaned.copy()
            
            # ========================================================
            # CRITICAL: ENSURE STANDARDIZED COLUMNS IN METABOLITES DATA
            # Downstream code expects: 'Name', 'pvalue', 'log2FC'
            # ========================================================
            if metabolites_data is not None:
                metabolites_data = metabolites_data.copy()
                
                # Check if standardization is needed (data from annotation should already have these)
                has_name = 'Name' in metabolites_data.columns
                has_pvalue = 'pvalue' in metabolites_data.columns
                has_log2fc = 'log2FC' in metabolites_data.columns
                
                if not (has_name and has_pvalue and has_log2fc):
                    logger.info(f"[STANDARDIZE] Standardizing column names for downstream processing...")
                    
                    # Standardize Name column (from LipidID, Lipid_ID, etc.)
                    if not has_name:
                        for candidate in ['LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite', 'name', 'Metabolite_Name']:
                            if candidate in metabolites_data.columns:
                                metabolites_data['Name'] = metabolites_data[candidate]
                                logger.info(f"[STANDARDIZE] Created 'Name' from '{candidate}'")
                                break
                    
                    # Standardize pvalue column
                    if not has_pvalue:
                        for candidate in ['p_value', 'P-Value', 'p-value', 'adj_p', 'P_value', 'pValue']:
                            if candidate in metabolites_data.columns:
                                metabolites_data['pvalue'] = metabolites_data[candidate]
                                logger.info(f"[STANDARDIZE] Created 'pvalue' from '{candidate}'")
                                break
                    
                    # Standardize log2FC column
                    if not has_log2fc:
                        for candidate in ['log2_FC', 'Log2FC', 'log2FoldChange', 'log2_fold_change', 'FC']:
                            if candidate in metabolites_data.columns:
                                metabolites_data['log2FC'] = metabolites_data[candidate]
                                logger.info(f"[STANDARDIZE] Created 'log2FC' from '{candidate}'")
                                break
                else:
                    logger.info(f"[STANDARDIZE] Standard columns already present (Name, pvalue, log2FC)")
            
            self.pathway_original_metabolites_data = metabolites_data.copy() if metabolites_data is not None else None
            self.pathway_filtered_metabolites_data = metabolites_data.copy() if metabolites_data is not None else None
            self._cache_assigned_columns(self.pathway_filtered_metabolites_data)

            # Update UI labels once we know which component type is active
            self._update_component_labels()
            
            # BUILD RELATIONSHIP MAPPINGS for smart cascading selections
            self._build_relationship_mappings(metabolites_data, cleaned)
            
            # Apply filters once so the component p-value threshold affects
            # metabolites, upstream, and diseases immediately
            self.apply_pathway_filter()
            
            logger.info(f"✅ Pathway Network tab populated successfully!")
            
        except Exception as e:
            logger.error(f"❌ Failed to set pathway data: {str(e)}")
            import traceback
            traceback.print_exc()

    # ===== Helpers copied/adapted from old GUI behavior =====
    def _split_pathways(self, pathway_str: str) -> list:
        """Split pathway strings by both ';' and '|' (old GUI compatible)."""
        import re
        if not pathway_str or not isinstance(pathway_str, str):
            return []
        parts = re.split(r'[;|]', pathway_str)
        return [p.strip() for p in parts if isinstance(p, str) and p.strip()]

    @staticmethod
    def _normalize_space_case(name: str) -> str:
        try:
            name = ' '.join(str(name).strip().split())
            return name
        except Exception:
            return str(name)

    def _excluded_pathway_keywords(self):
        """Keywords to exclude"""
        return [
            "action", "disease", "disorder", "syndrome", "deficiency", "toxicity",
            "drug", "medication", "therapy", "treatment", "pharmacology",
        ]

    def _is_database_label(self, name: str) -> bool:
        if not isinstance(name, str):
            return False
        s = name.strip().lower()
        if not s:
            return False
        db_tokens = {
            'hmdb', 'hmdb_id', 'kegg', 'kegg_id', 'smpdb', 'pathbank', 'pathbank_pathways',
            'wikipathways', 'metabolika', 'metabolika_pathways', 'reactome'
        }
        if s in db_tokens:
            return True
        if s.endswith('_id'):
            return True
        return False

    def _should_exclude_pathway_name(self, pathway_name: str) -> bool:
        """Return True if pathway should be excluded based on old-GUI keyword rules.
        
        NOTE: Disease/keyword filtering is now done BEFORE Fisher ORA statistics calculation.
        This function only removes obvious garbage (database labels, numeric junk, very short names).
        """
        if not isinstance(pathway_name, str):
            return True
        p = pathway_name.strip()
        if p == '':
            return True
        # remove single digits / numeric garbage
        if len(p) < 3:
            return True
        if p.replace(' ', '').isdigit():
            return True
        # explicit database labels
        if self._is_database_label(p):
            return True
        
        # COMMENTED OUT: Keyword-based exclusions now happen BEFORE statistics calculation
        # Disease pathways and unwanted keywords are filtered in pathway_annotation_tab.py
        # before Fisher ORA runs, so we don't need to filter again here
        # 
        # import re
        # pl = p.lower()
        # for kw in self._excluded_pathway_keywords():
        #     if re.search(r'\b' + re.escape(kw) + r'\b', pl):
        #         return True
        
        return False

    def _build_relationship_mappings(self, metabolites_df, pathways_dict):
        """
        Build relationship mappings for smart cascading selections.
        Creates bidirectional lookup dictionaries between pathways, metabolites,
        diseases, and upstream regulators.
        
        When class compression is active, metabolites_df already contains class nodes
        with aggregated pathways/diseases/upstream, so we build relationships directly
        from the class-level data.
        """
        import re
        
        # Initialize mapping dictionaries
        self.pathway_to_metabolites = {}
        self.metabolite_to_pathways = {}
        self.metabolite_to_diseases = {}
        self.metabolite_to_upstream = {}
        self.disease_to_metabolites = {}
        self.upstream_to_metabolites = {}
        
        if metabolites_df is None or metabolites_df.empty:
            logger.info("No metabolite data for relationship mapping")
            return
        
        # When class compression is active, metabolites_df is already at class level
        # No need to use raw data and remap - just build from the compressed data directly
        class_compression_active = getattr(self, 'class_compression_active', False)
        
        # DEBUG: Log available columns
        logger.info(f"📋 Available columns in metabolites data: {metabolites_df.columns.tolist()}")
        
        # DEBUG: Track pathway processing
        logger.info(f"[MAPPING DEBUG] Starting pathway relationship mapping...")
        logger.info(f"[MAPPING DEBUG] Input pathways_dict has {len([k for k in pathways_dict.keys() if not k.startswith('_')])} pathways")
        pathways_processed = 0
        pathways_skipped = []
        total_hits = 0
        total_metabolites = 0
        sample_pathway_shown = False
        
        # First, extract pathway-metabolite relationships from pathway stats
        # CRITICAL: Use 'hits' (significant metabolites) for cascade selection, NOT 'metabolites' (all measured)
        for pathway_name, stats in pathways_dict.items():
            if pathway_name.startswith('_'):  # Skip metadata
                continue
            if isinstance(stats, dict):
                # Prefer 'hits' (significant metabolites) over 'metabolites' (all measured)
                hits_list = stats.get('hits', [])
                metabolites_list = stats.get('metabolites', [])
                
                # DEBUG: Show sample pathway to verify data structure
                if not sample_pathway_shown:
                    logger.info(f"[MAPPING DEBUG] Sample pathway '{pathway_name}':")
                    logger.info(f"   - 'hits' field: {len(hits_list) if hits_list else 'MISSING/EMPTY'}")
                    logger.info(f"   - 'metabolites' field: {len(metabolites_list) if metabolites_list else 'MISSING/EMPTY'}")
                    logger.info(f"   - First 3 hits: {hits_list[:3] if hits_list else 'N/A'}")
                    sample_pathway_shown = True
                
                # Use hits if available, otherwise fall back to metabolites
                metabolite_list = hits_list if hits_list else metabolites_list
                total_hits += len(hits_list) if hits_list else 0
                total_metabolites += len(metabolites_list) if metabolites_list else 0
                
                if isinstance(metabolite_list, list) and len(metabolite_list) > 0:
                    pathways_processed += 1
                    self.pathway_to_metabolites[pathway_name] = set(metabolite_list)
                    for met in metabolite_list:
                        if met not in self.metabolite_to_pathways:
                            self.metabolite_to_pathways[met] = set()
                        self.metabolite_to_pathways[met].add(pathway_name)
                else:
                    pathways_skipped.append(f"{pathway_name} (no hits/metabolites or empty)")
            else:
                pathways_skipped.append(f"{pathway_name} (not a dict)")
        
        logger.info(f"[MAPPING DEBUG] Pathways processed: {pathways_processed}")
        logger.info(f"[MAPPING DEBUG] Pathways in pathway_to_metabolites: {len(self.pathway_to_metabolites)}")
        logger.info(f"[MAPPING DEBUG] Total 'hits' across all pathways: {total_hits}")
        logger.info(f"[MAPPING DEBUG] Total 'metabolites' across all pathways: {total_metabolites}")
        logger.info(f"[MAPPING DEBUG] Using 'hits' field (significant metabolites) for cascade selection")
        if pathways_skipped:
            logger.warning(f"[MAPPING DEBUG] Pathways skipped ({len(pathways_skipped)}): {pathways_skipped[:5]}...")
        
        # NOTE: Removed duplicate loop that was overwriting with 'metabolites' instead of 'hits'
        
        # Collect ALL metabolites that appear in any pathway's hits (for pathway cascade)
        all_hits_metabolites = set()
        for pathway_name, mets in self.pathway_to_metabolites.items():
            all_hits_metabolites.update(mets)
        logger.info(f"[MAPPING DEBUG] Unique metabolites in pathway hits: {len(all_hits_metabolites)}")
        
        # Build metabolite-disease and metabolite-upstream relationships
        # from the metabolites dataframe (already at class level if compression is active)
        all_significant_metabolites = set()
        for idx, row in metabolites_df.iterrows():
            metabolite_name = self._get_metabolite_value(row, 'name', '')
            if metabolite_name and not pd.isna(metabolite_name):
                all_significant_metabolites.add(metabolite_name)
        logger.info(f"[MAPPING DEBUG] Total compounds in data: {len(all_significant_metabolites)}")
        
        # Now build metabolite-disease and metabolite-upstream relationships
        # Directly from metabolites_df (which is class-level if compression is active)
        for idx, row in metabolites_df.iterrows():
            metabolite_name = self._get_metabolite_value(row, 'name', '')
            if not metabolite_name or pd.isna(metabolite_name):
                continue
            
            # Extract diseases - for ALL significant metabolites
            diseases_str = row.get('Associated_Diseases', '')
            if diseases_str and pd.notna(diseases_str) and str(diseases_str).strip():
                # Split by common separators
                disease_list = re.split(r'[;|,]', str(diseases_str))
                disease_list = [d.strip() for d in disease_list if d.strip()]
                
                for disease in disease_list:
                    if not self.metabolite_to_diseases.get(metabolite_name):
                        self.metabolite_to_diseases[metabolite_name] = set()
                    self.metabolite_to_diseases[metabolite_name].add(disease)
                    
                    if not self.disease_to_metabolites.get(disease):
                        self.disease_to_metabolites[disease] = set()
                    self.disease_to_metabolites[disease].add(metabolite_name)
            
            # Extract upstream regulators (enzymes and transporters)
            # Use the same columns as _populate_upstream_tree for consistency
            upstream_list = []
            
            # Enzymes - use 'Enzyme_Gene_name' column
            enzymes_str = row.get('Enzyme_Gene_name', '')
            if enzymes_str and pd.notna(enzymes_str) and str(enzymes_str).strip():
                enzyme_parts = re.split(r'[|;,]', str(enzymes_str))
                upstream_list.extend([e.strip() for e in enzyme_parts if e.strip() and str(e).lower() != 'nan'])
            
            # Transporters - use 'Transporter_Gene_name' if available, else fallback
            transporters_str = row.get('Transporter_Gene_name', '') or row.get('Transporters', '')
            if transporters_str and pd.notna(transporters_str) and str(transporters_str).strip():
                transporter_parts = re.split(r'[|;,]', str(transporters_str))
                upstream_list.extend([t.strip() for t in transporter_parts if t.strip() and str(t).lower() != 'nan'])
            
            # Build upstream mappings - for ALL significant metabolites
            for upstream in upstream_list:
                if not upstream:
                    continue
                if not self.metabolite_to_upstream.get(metabolite_name):
                    self.metabolite_to_upstream[metabolite_name] = set()
                self.metabolite_to_upstream[metabolite_name].add(upstream)
                
                if not self.upstream_to_metabolites.get(upstream):
                    self.upstream_to_metabolites[upstream] = set()
                self.upstream_to_metabolites[upstream].add(metabolite_name)
        
        # Log relationship stats
        comp_label = self._component_label(plural=True)
        logger.info(f"📊 Relationship Mappings Built:")
        logger.info(f"   Pathways → {comp_label}: {len(self.pathway_to_metabolites)} pathways")
        logger.info(f"   {comp_label} → Pathways: {len(self.metabolite_to_pathways)} {comp_label.lower()}")
        logger.info(f"   {comp_label} → Diseases: {len(self.metabolite_to_diseases)} {comp_label.lower()}")
        logger.info(f"   Diseases → {comp_label}: {len(self.disease_to_metabolites)} diseases")
        logger.info(f"   {comp_label} → Upstream: {len(self.metabolite_to_upstream)} {comp_label.lower()}")
        logger.info(f"   Upstream → {comp_label}: {len(self.upstream_to_metabolites)} upstream regulators")
        
        # Debug: show what upstream regulators we found
        if self.upstream_to_metabolites:
            logger.info(f"   Upstream regulators found: {sorted(list(self.upstream_to_metabolites.keys())[:10])}")
        else:
            logger.warning(f"   ⚠️ NO upstream regulators found in {comp_label.lower()} data!")

    def _compute_metabolite_sources(self, df: Optional[pd.DataFrame] = None):
        """Build a mapping from metabolite name to contributing sources (omics).

        Uses the 'Source' column created by the Multi-Omics tab.
        Falls back gracefully to original metabolites data if the
        filtered table is empty. If no usable 'Source' information is
        found, returns an empty dict and has no effect on ranking.
        """
        def _find_source_column(frame: Optional[pd.DataFrame]) -> Optional[str]:
            """Return best matching source column name from a dataframe."""
            if frame is None or frame.empty:
                return None
            cols = list(getattr(frame, 'columns', []))
            if not cols:
                return None

            # Strict mode: only the exact canonical name is accepted.
            if 'Source' in cols:
                return 'Source'

            return None

        try:
            if df is None:
                df = getattr(self, 'pathway_filtered_metabolites_data', None)
        except Exception:
            df = None

        source_col = _find_source_column(df)

        if source_col is not None:
            try:
                logger.info(
                    f"[MULTIOMICS] Source column detected: {source_col} | "
                    f"rows={len(df) if df is not None else 0} | "
                    f"columns={len(df.columns) if df is not None else 0}"
                )
            except Exception:
                pass

        # If filtered data is unavailable or has no usable source column,
        # fall back to the original metabolites.
        if df is None or df.empty or source_col is None:
            try:
                orig_df = getattr(self, 'pathway_original_metabolites_data', None)
            except Exception:
                orig_df = None
            source_col = _find_source_column(orig_df)
            if orig_df is None or orig_df.empty or source_col is None:
                logger.info("[MULTIOMICS] No usable 'Source' column found in metabolites data; multi-omics rank disabled.")
                return {}
            df = orig_df

        metabolite_sources = {}

        import re as _re

        if df is None:  # Should not reach here due to guards above, but add type guard
            return metabolite_sources

        for _, row in df.iterrows():
            # Prefer canonical Name to avoid mismatches from stale/manual assignments.
            name = row.get('Name', None) if 'Name' in row.index else None
            if name is None or (isinstance(name, float) and pd.isna(name)):
                name = self._get_metabolite_value(row, 'name', None)
            if name is None or (isinstance(name, float) and pd.isna(name)):
                continue

            src_val = row.get(source_col, None)
            if src_val is None or pd.isna(src_val):
                continue

            src_str = str(src_val).strip()
            if not src_str or src_str.lower() == 'nan':
                continue

            # Split on common separators between source labels
            parts = _re.split(r'[+;|,]', src_str)
            tokens = [p.strip() for p in parts if p.strip()]
            if not tokens:
                continue

            if name not in metabolite_sources:
                metabolite_sources[name] = set()
            metabolite_sources[name].update(tokens)

        # DEBUG: log basic stats so users can verify multi-omics info
        try:
            unique_sources = sorted({s for v in metabolite_sources.values() for s in v})
            self._last_multiomics_source_info = {
                'source_col': source_col,
                'n_metabolites_with_source': len(metabolite_sources),
                'unique_sources': unique_sources,
            }
            logger.info(f"[MULTIOMICS] Metabolites with Source info: {len(metabolite_sources)}")
            logger.info(f"[MULTIOMICS] Unique sources detected: {unique_sources}")
            if metabolite_sources:
                preview = []
                for met_name in list(metabolite_sources.keys())[:5]:
                    preview.append(f"{met_name}={sorted(metabolite_sources[met_name])}")
                logger.info(f"[MULTIOMICS] Source preview: {preview}")
        except Exception:
            pass

        return metabolite_sources

    def _compute_multiomics_score(self, metabolites, metabolite_sources, entity_name=None, debug=False):
        """Compute a multi-omics integration score for a set of metabolites.

        The score rewards entities that are supported by multiple sources
        (e.g. metabolites and lipids) in a *balanced* way.

        Heuristic:
        - Count how many times each source contributes across the metabolite set
        - If fewer than 2 sources contribute, score is 0 (single-omics)
        - Otherwise, score = min(counts) + 0.001 * total_counts
          so 5+5 (min=5,total=10) outranks 7+2 (min=2,total=9).
        Returns (score, n_sources, total_hits).
        
        Args:
            metabolites: List/set of metabolite names
            metabolite_sources: Dict mapping metabolite names to source sets
            entity_name: Name of pathway/upstream/disease for debug logging
            debug: Enable detailed logging
        """
        if not metabolites or not metabolite_sources:
            if debug and entity_name:
                logger.info(
                    f"[MULTIOMICS DEBUG] {entity_name}: skipped before scoring because "
                    f"metabolites_present={bool(metabolites)} metabolite_sources_present={bool(metabolite_sources)}"
                )
            return 0.0, 0, 0

        if debug and entity_name:
            logger.info(
                f"[MULTIOMICS DEBUG] {entity_name}: starting score calculation with "
                f"{len(metabolites)} metabolites and {len(metabolite_sources)} metabolite-source mappings"
            )

        source_counts = {}
        contributing_metabolites = {}  # Track which metabolites contribute to each source
        
        for met in metabolites:
            src_set = metabolite_sources.get(met, set())
            if not src_set:
                continue
            for src in src_set:
                source_counts[src] = source_counts.get(src, 0) + 1
                if src not in contributing_metabolites:
                    contributing_metabolites[src] = []
                contributing_metabolites[src].append(met)

        n_sources = len(source_counts)
        total_hits = sum(source_counts.values())

        if debug and entity_name:
            logger.info(
                f"[MULTIOMICS DEBUG] {entity_name}: source_counts={source_counts}, total_hits={total_hits}, n_sources={n_sources}"
            )
        
        # DEBUG: Log detailed breakdown for entities with multi-omics support
        if debug and n_sources > 1 and entity_name:
            logger.info(f"[MULTIOMICS DEBUG] {entity_name}:")
            logger.info(f"   Total metabolites: {len(metabolites)}")
            logger.info(f"   Sources detected: {n_sources}")
            for src, count in sorted(source_counts.items()):
                mets = contributing_metabolites[src]
                logger.info(f"   • {src}: {count} metabolites - {mets[:3]}{'...' if len(mets) > 3 else ''}")
        
        if n_sources <= 1 or total_hits == 0:
            if debug and entity_name:
                logger.info(
                    f"[MULTIOMICS DEBUG] {entity_name}: no multi-omics score because "
                    f"n_sources={n_sources} total_hits={total_hits}"
                )
            return 0.0, n_sources, total_hits

        min_count = min(source_counts.values())
        score = float(min_count) + 0.001 * float(total_hits)

        if debug and entity_name:
            logger.info(
                f"[MULTIOMICS DEBUG] {entity_name}: score computed as min_count={min_count} + "
                f"0.001*total_hits={total_hits} => {score:.3f}"
            )
        return score, n_sources, total_hits

    def _normalize_pathway_key(self, pathway_name: str) -> str:
        """Normalize pathway name for deduplication (case-insensitive, collapsed spaces)."""
        if not isinstance(pathway_name, str):
            return str(pathway_name)
        import re
        name = self._normalize_space_case(pathway_name)
        # Remove all non-alphanumeric characters and collapse to lowercase to match
        # normalization used in analysis modules. This ensures punctuation/case
        # variants deduplicate consistently across backend and GUI.
        return re.sub(r"[^a-z0-9]+", "", name.lower())

    def _merge_pathway_stats(self, a: dict, b: dict) -> dict:
        """Merge two pathway stats dicts conservatively (old GUI spirit: reconcile as one).
        - n_metabolites: max of the two (or union of metabolite lists when present)
        - combined_pvalue: min (more significant)
        - z_score: choose by absolute magnitude
        - status: derived from chosen z_score if missing
        - mean_log2fc: average if both present
        - metabolites: union unique (if lists/strings provided)
        - hits: union unique (CRITICAL for network - significant metabolites only)
        """
        res = dict(a) if isinstance(a, dict) else {}
        other = b if isinstance(b, dict) else {}

        # combine metabolites list if available
        def _extract_mets(d, key='metabolites'):
            mets = d.get(key) if isinstance(d, dict) else None
            if isinstance(mets, list):
                return [str(x).strip() for x in mets if str(x).strip()]
            if isinstance(mets, str):
                # try split by comma
                return [m.strip() for m in mets.split(',') if m.strip()]
            return []

        # Merge 'metabolites' (all measured)
        set_a = set(_extract_mets(res, 'metabolites'))
        set_b = set(_extract_mets(other, 'metabolites'))
        if set_a or set_b:
            merged = sorted(set_a.union(set_b))
            res['metabolites'] = merged
            res['n_metabolites'] = len(merged)
        else:
            # fallback to max of provided counts
            na = res.get('n_metabolites') if isinstance(res.get('n_metabolites'), (int, float)) else 0
            nb = other.get('n_metabolites') if isinstance(other.get('n_metabolites'), (int, float)) else 0
            res['n_metabolites'] = int(max(na or 0, nb or 0))

        # CRITICAL: Also merge 'hits' (significant metabolites only)
        hits_a = set(_extract_mets(res, 'hits'))
        hits_b = set(_extract_mets(other, 'hits'))
        if hits_a or hits_b:
            merged_hits = sorted(hits_a.union(hits_b))
            res['hits'] = merged_hits
            res['k_hits'] = len(merged_hits)
        elif 'k_hits' in res or 'k_hits' in other:
            # Fallback to max of k_hits counts
            ka = res.get('k_hits', 0) if isinstance(res.get('k_hits'), (int, float)) else 0
            kb = other.get('k_hits', 0) if isinstance(other.get('k_hits'), (int, float)) else 0
            res['k_hits'] = int(max(ka, kb))

        # p-value: take the most significant (min)
        pa = res.get('combined_pvalue', 1.0)
        pb = other.get('combined_pvalue', 1.0)
        try:
            res['combined_pvalue'] = min(float(pa), float(pb))
        except Exception:
            res['combined_pvalue'] = pa if isinstance(pa, (int, float)) else pb

        # z-score: choose by absolute magnitude when both present
        za = res.get('z_score')
        zb = other.get('z_score')
        try:
            if za is None:
                res['z_score'] = zb
            elif zb is None:
                res['z_score'] = za
            else:
                res['z_score'] = za if abs(float(za)) >= abs(float(zb)) else zb
        except Exception:
            pass

        # mean log2fc: simple average if both present
        la = res.get('mean_log2fc')
        lb = other.get('mean_log2fc')
        try:
            if la is not None and lb is not None:
                res['mean_log2fc'] = (float(la) + float(lb)) / 2.0
            elif lb is not None and la is None:
                res['mean_log2fc'] = lb
        except Exception:
            pass

        # status: PRESERVE the original status from pathway annotation
        # Do NOT recalculate from z-score because the original calculation uses
        # thresholds and neutral bands, not just positive/negative
        # Prefer the status from the entry with larger |z_score|
        za = res.get('z_score')
        zb = other.get('z_score')
        if za is not None and zb is not None:
            # Use status from the entry with the largest absolute z-score
            try:
                if abs(float(za)) >= abs(float(zb)):
                    res['status'] = res.get('status', 'No Change')
                else:
                    res['status'] = other.get('status', 'No Change')
            except Exception:
                res['status'] = res.get('status', other.get('status', 'No Change'))
        else:
            # Use whichever has a status defined
            res['status'] = res.get('status', other.get('status', 'No Change'))

        return res

    def _clean_and_merge_pathways_data(self, pathways_data: dict) -> dict:
        """Clean pathway dict like old GUI:
        - Exclude records with unwanted keywords (Action, Disease, etc.)
        - Drop database labels (HMDB, KEGG_ID, etc.)
        - Remove numeric/short garbage like '5', '1'
        - Case-insensitive dedup of names (beta-Alanine metabolism variants)
        - Merge duplicates conservatively
        Returns a new cleaned dict.
        """
        try:
            import re
            if not isinstance(pathways_data, dict):
                return {}

            # DEBUG: Track what gets filtered
            excluded_pathways = []
            merged_pathways = []

            # Heuristic keywords that often appear at the end of pathway names
            suffix_keywords = [
                'metabolism', 'degradation', 'biosynthesis', 'synthesis', 'pathway',
                'oxidation', 'reduction', 'transport', 'signaling', 'signalling',
                'metabolic', 'catabolism', 'anabolism', 'cycle'
            ]

            def _split_concatenated(pathname: str) -> list:
                """DISABLED: Previously tried to split concatenated pathway labels.
                This was causing pathway count inflation (25 -> 28) because it incorrectly
                split valid pathway names that happened to contain multiple suffix keywords.
                Now we simply return the original pathway name as-is.
                """
                if not pathname or not isinstance(pathname, str):
                    return []
                # DISABLED: Return original name without splitting
                # The splitting logic was incorrectly inflating pathway counts
                return [pathname.strip()] if pathname.strip() else []
                
                # OLD CODE (disabled):
                # If obvious separators present, rely on existing split
                if False and re.search(r'[;|]', pathname):
                    return [p.strip() for p in re.split(r'[;|]', pathname) if p.strip()]

                low = pathname.lower()
                for kw in suffix_keywords:
                    # find all occurrences of the keyword
                    matches = list(re.finditer(r'\b' + re.escape(kw) + r'\b', low))
                    if len(matches) >= 2:
                        # split after the end of the first occurrence
                        split_idx = matches[0].end()
                        left = pathname[:split_idx].strip()
                        right = pathname[split_idx:].strip()
                        # Validate both parts roughly before accepting split
                        if left and right and self._is_valid_pathway(left) and self._is_valid_pathway(right):
                            return [left, right]

                # No good split found; return original as single candidate
                return [pathname]

            merged = {}
            display_name_for = {}

            for k, stats in pathways_data.items():
                if not isinstance(k, str):
                    excluded_pathways.append((k, "not a string"))
                    continue
                if k.startswith('_'):
                    # ignore metadata entries
                    continue

                # Break concatenated entries into candidate names first
                candidates = _split_concatenated(k)
                for cand in candidates:
                    name = self._normalize_space_case(cand)
                    if self._should_exclude_pathway_name(name):
                        excluded_pathways.append((name, "excluded by keyword filter"))
                        continue
                    # final validity check (reuse existing validator)
                    if not self._is_valid_pathway(name):
                        excluded_pathways.append((name, "failed validity check"))
                        continue
                    norm = self._normalize_pathway_key(name)
                    if norm in merged:
                        merged_pathways.append((name, f"merged with {display_name_for[norm]}"))
                        merged[norm] = self._merge_pathway_stats(merged[norm], stats if isinstance(stats, dict) else {})
                    else:
                        merged[norm] = stats if isinstance(stats, dict) else {}
                        display_name_for[norm] = name

            # Rebuild dict with display names preserved
            cleaned = {}
            for norm, stats in merged.items():
                disp = display_name_for.get(norm, norm)
                cleaned[disp] = stats
            # Preserve any metadata keys from original (prefixed with '_')
            for k, v in pathways_data.items():
                if isinstance(k, str) and k.startswith('_'):
                    cleaned[k] = v

            # DEBUG: Log filtering details
            if excluded_pathways:
                logger.warning(f"[CLEAN DEBUG] Pathways excluded ({len(excluded_pathways)}):")
                for pathway, reason in excluded_pathways:
                    logger.warning(f"  - '{pathway}': {reason}")
            if merged_pathways:
                logger.info(f"[CLEAN DEBUG] Pathways merged ({len(merged_pathways)}):")
                for pathway, reason in merged_pathways[:3]:  # Show first 3
                    logger.info(f"  - '{pathway}': {reason}")

            return cleaned
        except Exception as e:
            logger.error(f"Failed to clean/merge pathways: {e}")
            import traceback
            traceback.print_exc()
            return pathways_data
    
    def _populate_pathway_tree(self):
        """Populate the pathway tree table with pathway data"""
        try:
            if not hasattr(self, 'pathway_tree') or self.pathway_tree is None:
                logger.warning("Pathway tree widget not initialized yet")
                return
            
            # Clear existing data
            for item in self.pathway_tree.get_children():
                self.pathway_tree.delete(item)
            
            if not self.pathway_filtered_pathways_data:
                return
            
            # DEBUG: Track pathway filtering
            logger.info(f"[POPULATE DEBUG] Starting pathway tree population...")
            logger.info(f"[POPULATE DEBUG] pathway_filtered_pathways_data has {len([k for k in self.pathway_filtered_pathways_data.keys() if not k.startswith('_')])} pathways")
            
            # Import centralized pathway filtering function
            try:
                from main_script.metabolite_pathway_network import should_filter_pathway
            except ImportError:
                # Fallback if import fails
                def should_filter_pathway(pathway_name):
                    """Minimal fallback filtering"""
                    if not pathway_name:
                        return True
                    pathway_lower = pathway_name.lower()
                    disease_keywords = ["disease", "disorder", "syndrome", "emia", "uria", "osis"]
                    for keyword in disease_keywords:
                        if keyword in pathway_lower:
                            return True
                    return False
            
            # Reset selections
            self.pathway_selections = {}
            
            # CRITICAL: Additional deduplication check before displaying
            # Ensure no duplicate pathways (case-insensitive) appear in tree
            # Prioritize highest |z-score| over lowest p-value (biological effect > statistical significance)
            normalized_to_original = {}
            deduplicated_data = {}
            for pathway_name, stats in self.pathway_filtered_pathways_data.items():
                if not isinstance(pathway_name, str) or pathway_name.startswith('_'):
                    continue
                # Normalize name for comparison
                normalized_name = ' '.join(pathway_name.lower().strip().split())
                if normalized_name in normalized_to_original:
                    # Duplicate found - keep the one with higher |z-score| (stronger biological effect)
                    existing_name = normalized_to_original[normalized_name]
                    existing_z = abs(deduplicated_data[existing_name].get('z_score', 0.0))
                    current_z = abs(stats.get('z_score', 0.0))
                    if current_z > existing_z:
                        # Current has stronger effect, replace
                        del deduplicated_data[existing_name]
                        deduplicated_data[pathway_name] = stats
                        normalized_to_original[normalized_name] = pathway_name
                        logger.warning(f"[DEDUP] Found duplicate in tree: kept '{pathway_name}' (|z|={current_z:.3f}), removed '{existing_name}' (|z|={existing_z:.3f})")
                else:
                    deduplicated_data[pathway_name] = stats
                    normalized_to_original[normalized_name] = pathway_name
            
            # Use deduplicated data
            self.pathway_filtered_pathways_data = deduplicated_data
            
            # ALWAYS apply ML augmentation to ensure columns populated (idempotent)
            try:
                if self.pathway_filtered_metabolites_data is not None:
                    from main_script.ml_pathway_models import augment_with_ml_scores
                    self.pathway_filtered_pathways_data = augment_with_ml_scores(
                        self.pathway_filtered_metabolites_data, 
                        self.pathway_filtered_pathways_data
                    )
                    logger.info("✓ ML augmentation applied for pathways (network tab).")
            except Exception as _e:
                logger.warning(f"ML augmentation skipped: {_e}")

            # First pass: collect displayed pathways ONLY (after disease/general filtering)
            items = []
            disease_pathways_filtered = []
            general_pathways_filtered = []
            invalid_stats_filtered = []
            
            for pathway_name, stats in self.pathway_filtered_pathways_data.items():
                if not isinstance(pathway_name, str) or pathway_name.startswith('_'):
                    continue
                # Filter out disease pathways and overly general pathways using centralized function
                if should_filter_pathway(pathway_name):
                    disease_pathways_filtered.append(pathway_name)
                    continue
                if not isinstance(stats, dict):
                    invalid_stats_filtered.append(f"{pathway_name} (stats not dict)")
                    continue

                try:
                    z_score = float(stats.get('z_score', 0.0) or 0.0)
                except Exception:
                    z_score = 0.0
                # Extract BOTH raw p-value and FDR for display
                try:
                    # Try to get raw p-value (fisher_pvalue or pvalue)
                    raw_pvalue = float(stats.get('fisher_pvalue', stats.get('pvalue', 1.0)) or 1.0)
                except Exception:
                    raw_pvalue = 1.0
                
                # Try to get FDR-adjusted p-value (fdr_qvalue or adjusted_pvalue)
                # Keep as None if not present (don't fall back to raw_pvalue)
                fdr_value = stats.get('fdr_qvalue', stats.get('adjusted_pvalue', None))
                if fdr_value is not None:
                    try:
                        fdr_pvalue = float(fdr_value)
                    except Exception:
                        fdr_pvalue = None
                else:
                    fdr_pvalue = None
                
                # For backward compatibility: combined_pvalue is what gets used for filtering
                try:
                    p_value = float(stats.get('combined_pvalue', raw_pvalue) or raw_pvalue)
                except Exception:
                    p_value = raw_pvalue
                n_metabolites = stats.get('n_metabolites', 0) or 0
                # CRITICAL: Use k_hits (significant metabolites) for display in network tab
                # k_hits = number of significant metabolites that passed p-value threshold
                # n_metabolites = ALL measured metabolites in pathway (for Fisher ORA stats)
                k_hits = stats.get('k_hits', n_metabolites) or n_metabolites
                status = stats.get('status', 'No Change')
                try:
                    ml_score = float(stats.get('ml_rwr_score', 0.0) or 0.0)
                except Exception:
                    ml_score = 0.0
                try:
                    pr = float(stats.get('ml_pagerank', 0.0) or 0.0)
                except Exception:
                    pr = 0.0
                try:
                    btw = float(stats.get('ml_betweenness', 0.0) or 0.0)
                except Exception:
                    btw = 0.0

                items.append({
                    'name': pathway_name,
                    'status': status,
                    'z': z_score,
                    'p': p_value,
                    'raw_p': raw_pvalue,  # Store raw p-value separately
                    'fdr_p': fdr_pvalue,  # Store FDR p-value separately
                    'n': k_hits,  # Use k_hits (significant metabolites) instead of n_metabolites (all measured)
                    'n_sources': 0,  # Will be populated after multi-omics calculation
                    'mo_score': 0.0,  # Will be populated after multi-omics calculation
                    'ml': ml_score,
                    'pr': pr,
                    'btw': btw,
                })
                
                # DEBUG: Log k_hits vs n_metabolites for first few pathways
                if len(items) <= 5:
                    logger.info(f"[PATHWAY DEBUG] {pathway_name}: k_hits={k_hits}, n_metabolites={n_metabolites}")
            
            # DEBUG: Log filtering results
            logger.info(f"[POPULATE DEBUG] Pathways after disease filter: {len(items)}")
            if disease_pathways_filtered:
                logger.warning(f"[POPULATE DEBUG] Disease pathways filtered out ({len(disease_pathways_filtered)}): {disease_pathways_filtered[:3]}...")
            if general_pathways_filtered:
                logger.warning(f"[POPULATE DEBUG] Overly general pathways filtered out ({len(general_pathways_filtered)}): {general_pathways_filtered}")
            if invalid_stats_filtered:
                logger.warning(f"[POPULATE DEBUG] Invalid stats filtered out ({len(invalid_stats_filtered)}): {invalid_stats_filtered}")

            # CRITICAL: Compute multi-omics scores ONLY for displayed pathways (after disease/general filtering)
            # This ensures multi-omics ranks are relative to the actual displayed set
            logger.info(
                f"[MULTIOMICS] About to score displayed pathways: {len(items)} visible pathways, "
                f"{len(self.pathway_to_metabolites) if hasattr(self, 'pathway_to_metabolites') else 0} pathways with metabolite mappings"
            )
            metabolite_sources = self._compute_metabolite_sources()
            try:
                info = getattr(self, '_last_multiomics_source_info', {}) or {}
                logger.info(
                    f"[MULTIOMICS] Source summary before scoring: source_col={info.get('source_col')} | "
                    f"n_metabolites_with_source={info.get('n_metabolites_with_source', 0)} | "
                    f"unique_sources={info.get('unique_sources', [])}"
                )
            except Exception:
                pass
            pathway_multiomics = {}
            if metabolite_sources and hasattr(self, 'pathway_to_metabolites'):
                # Only compute for pathways in items (displayed pathways)
                displayed_pathway_names = {it['name'] for it in items}
                for idx, p_name in enumerate(sorted(displayed_pathway_names)):
                    if p_name in self.pathway_to_metabolites:
                        mets = self.pathway_to_metabolites[p_name]
                        # Enable debug for top 10 pathways to track multi-omics breakdown
                        debug_enabled = (idx < 10)
                        if debug_enabled:
                            logger.info(
                                f"[MULTIOMICS] Scoring pathway '{p_name}' with {len(mets)} mapped metabolites"
                            )
                        score, n_sources, total_hits = self._compute_multiomics_score(
                            mets, metabolite_sources, entity_name=p_name, debug=debug_enabled
                        )
                        pathway_multiomics[p_name] = {
                            'score': score,
                            'n_sources': n_sources,
                            'total_hits': total_hits,
                        }
            
            # Update items with multi-omics scores
            for it in items:
                mo_info = pathway_multiomics.get(it['name'], {})
                it['mo_score'] = float(mo_info.get('score', 0.0) or 0.0)
                it['n_sources'] = int(mo_info.get('n_sources', 0) or 0)
                it['mo_total_hits'] = int(mo_info.get('total_hits', 0) or 0)

            # Compute ranks for consensus (rank by |z| desc and ML desc)
            # Handle ties by stable ordering; missing values already coerced to 0
            z_sorted = sorted(items, key=lambda x: abs(x['z']), reverse=True)
            ml_sorted = sorted(items, key=lambda x: x['ml'], reverse=True)
            rank_z = {it['name']: idx + 1 for idx, it in enumerate(z_sorted)}
            rank_ml = {it['name']: idx + 1 for idx, it in enumerate(ml_sorted)}

            # Compute multi-omics ranks (only pathways with positive score)
            mo_sorted = sorted(items, key=lambda x: x['mo_score'], reverse=True)
            mo_rank = {}
            rank_counter = 1
            mo_positive = 0
            for it in mo_sorted:
                if it['mo_score'] <= 0:
                    continue
                mo_positive += 1
                mo_rank[it['name']] = rank_counter
                rank_counter += 1

            # DEBUG: log how many pathways have multi-omics support
            try:
                logger.info(f"[MULTIOMICS] Pathways with multi-omics support (score>0): {mo_positive}")
                if mo_positive:
                    top_preview = mo_sorted[: min(5, len(mo_sorted))]
                    logger.info("[MULTIOMICS] Top pathways by multi-omics score:")
                    for it in top_preview:
                        if it['mo_score'] > 0:
                            logger.info(
                                f"   {it['name']}: score={it['mo_score']:.3f}, sources={it['n_sources']}, "
                                f"hits={it['n']}, total_source_hits={it.get('mo_total_hits', 0)}"
                            )
            except Exception:
                pass

            # Order items according to Rank by selection
            try:
                mode = self.pathway_rank_mode.get() if self.pathway_rank_mode else 'Default'
            except Exception:
                mode = 'Default'

            # CRITICAL: Always prioritize Activated/Inhibited over No Change
            # CRITICAL: Pathways with z-score = 0 (No Change) should ALWAYS rank below those with activation/inhibition
            # Then sort by: |z-score| desc, #metabolites desc, consensus asc, p-value asc
            def get_status_priority(status, z_score):
                """Return priority: Activated/Inhibited with z≠0 = 0 (top), No Change or z=0 = 1 (bottom)"""
                if status in ('Activated', 'Inhibited') and abs(z_score) > 0:
                    return 0  # Real activation/inhibition with non-zero z-score
                return 1  # No Change or zero z-score goes to bottom
            
            # Add consensus to items for sorting
            for it in items:
                it['consensus'] = int(round((rank_z.get(it['name'], 0) + rank_ml.get(it['name'], 0)) / 2.0)) if items else 0
                it['mo_rank'] = mo_rank.get(it['name'], 0)

            if mode == 'ML':
                # ML mode: status priority, then multi-omics score, then # compounds, then ML score, then |z|, then p-value
                items_sorted = sorted(items, key=lambda x: (
                    get_status_priority(x['status'], x['z']),  # 1. Activated/Inhibited (z≠0) first, No Change/z=0 last
                    -x['mo_score'],                            # 2. Stronger multi-omics support first
                    -x['n'],                                   # 3. More metabolites first (# significant metabolites)
                    -(x['ml']),                                # 4. Higher ML score first
                    -abs(x['z']),                              # 5. Higher |z| first
                    x['p']                                     # 6. Lower p-value first
                ))
            elif mode == 'Consensus':
                items_sorted = sorted(items, key=lambda x: (
                    get_status_priority(x['status'], x['z']),  # 1. Activated/Inhibited (z≠0) first, No Change/z=0 last
                    -x['mo_score'],                            # 2. Stronger multi-omics support first
                    -x['n'],                                   # 3. More metabolites first (# significant metabolites)
                    x['consensus'],                            # 4. Lower consensus rank first
                    -(x['ml']),                                # 5. Higher ML score first
                    -abs(x['z']),                              # 6. Higher |z| first
                    x['p']                                     # 7. Lower p-value first
                ))
            else:
                # Default: ALWAYS prioritize Activated/Inhibited over No Change
                # CRITICAL: Pathways with z-score = 0 should NEVER rank above those with z ≠ 0
                # Then: multi-omics score desc, # compounds desc, |z-score| desc, consensus asc, p-value asc
                items_sorted = sorted(items, key=lambda x: (
                    get_status_priority(x['status'], x['z']),  # 1. Activated/Inhibited (z≠0) first, No Change/z=0 last
                    -x['mo_score'],                            # 2. Stronger multi-omics support first
                    -x['n'],                                   # 3. More metabolites first (# significant metabolites)
                    -abs(x['z']),                              # 4. Higher |z-score| first
                    x['consensus'],                            # 5. Lower consensus rank first
                    x['p']                                     # 6. Lower p-value first
                ))
            
            logger.info(f"[SORT DEBUG] Sorting mode: {mode}")
            
            # Count status distribution
            activated_count = sum(1 for it in items_sorted if it['status'] == 'Activated')
            inhibited_count = sum(1 for it in items_sorted if it['status'] == 'Inhibited')
            no_change_count = sum(1 for it in items_sorted if it['status'] == 'No Change')
            logger.info(f"[SORT DEBUG] Status distribution: Activated={activated_count}, Inhibited={inhibited_count}, No Change={no_change_count}")
            
            logger.info(f"[SORT DEBUG] First 10 pathways after sorting (prioritized by Activated/Inhibited):")
            for i, it in enumerate(items_sorted[:10]):
                logger.info(f"   {i+1}. {it['name']}: status={it['status']}, z={it['z']:.3f}, n={it['n']}, consensus={it['consensus']}, p={it['p']:.2e}")

            # Apply P-value threshold filter
            try:
                pvalue_threshold = float(self.pathway_pvalue_threshold.get())
            except Exception:
                pvalue_threshold = 0.05
            
            logger.info(f"[POPULATE DEBUG] P-value threshold: {pvalue_threshold}")
            
            # Apply ML filters only if explicitly set (not default values)
            try:
                ml_score_min_str = self.pathway_ml_score_min.get().strip()
                consensus_rank_max_str = self.pathway_consensus_rank_max.get().strip()
                
                # Only apply if user changed from defaults
                apply_ml_filter = (ml_score_min_str and ml_score_min_str != "0.0" and float(ml_score_min_str) > 0)
                apply_consensus_filter = (consensus_rank_max_str and consensus_rank_max_str != "999999" and int(consensus_rank_max_str) < 999999)
                
                ml_score_min = float(ml_score_min_str) if apply_ml_filter else -999999
                consensus_rank_max = int(consensus_rank_max_str) if apply_consensus_filter else 999999
            except Exception:
                apply_ml_filter = False
                apply_consensus_filter = False
                ml_score_min = -999999
                consensus_rank_max = 999999
            
            # Get max pathways limit
            try:
                limit_enabled = hasattr(self, 'pathway_limit_pathways') and self.pathway_limit_pathways.get()
                max_pathways = self.pathway_max_total.get() if hasattr(self, 'pathway_max_total') else 25
            except Exception:
                limit_enabled = True
                max_pathways = 25
            
            logger.info(f"[POPULATE DEBUG] Max pathways limit: {max_pathways} (enabled: {limit_enabled})")
            
            # Populate tree using computed ranks
            ml_filtered = 0
            pvalue_filtered = 0
            pathway_count = 0
            displayed_items = []
            for it in items_sorted:
                # FIRST: Apply p-value threshold filter (ALWAYS applied)
                if it['p'] > pvalue_threshold:
                    pvalue_filtered += 1
                    continue
                
                # Check max pathways limit SECOND
                if limit_enabled and pathway_count >= max_pathways:
                    logger.info(f"[POPULATE DEBUG] Reached max pathways limit ({max_pathways}), stopping population")
                    break
                    
                name = it['name']
                consensus = it.get('consensus', 0)  # Use pre-computed consensus
                
                # Apply ML filters only if user set them
                if (apply_ml_filter and it['ml'] < ml_score_min) or (apply_consensus_filter and consensus > consensus_rank_max):
                    ml_filtered += 1
                    continue

                # Mark as UNSELECTED by default
                self.pathway_selections[name] = False
                select_symbol = "☐"
                
                # Determine FDR display value
                apply_fdr = self.pathway_apply_fdr.get() if hasattr(self, 'pathway_apply_fdr') else True
                if not apply_fdr or it['fdr_p'] is None:
                    fdr_display = "-"  # Show "-" if FDR not applied or no FDR value available
                else:
                    fdr_display = f"{it['fdr_p']:.2e}"

                # Only display a rank when multi-omics support exists; otherwise leave blank
                mo_rank_display = it.get('mo_rank', 0)
                mo_rank_display = mo_rank_display if mo_rank_display > 0 else ''

                self.pathway_tree.insert('', 'end', values=(
                    select_symbol,
                    name,
                    it['status'],
                    f"{it['z']:.3f}",
                    f"{it['raw_p']:.2e}",  # Raw p-value
                    fdr_display,           # FDR-adjusted p-value (or "-")
                    it['n'],
                    mo_rank_display,
                    f"{it['ml']:.4f}",
                    f"{it['pr']:.4f}",
                    f"{it['btw']:.4f}",
                    consensus
                ))
                displayed_items.append(it)
                pathway_count += 1
            
            if pvalue_filtered > 0:
                logger.info(f"   P-value filter removed {pvalue_filtered} pathways (p ≥ {pvalue_threshold})")
            if ml_filtered > 0:
                logger.info(f"   ML filters removed {ml_filtered} pathways")
            logger.info(f"✅ Pathway tree populated with {len(self.pathway_selections)} pathways")
            
            # Count status distribution of DISPLAYED pathways only (for stats label)
            displayed_activated = sum(1 for name in self.pathway_selections.keys() 
                                     if self.pathway_filtered_pathways_data.get(name, {}).get('status') == 'Activated')
            displayed_inhibited = sum(1 for name in self.pathway_selections.keys() 
                                     if self.pathway_filtered_pathways_data.get(name, {}).get('status') == 'Inhibited')
            displayed_no_change = sum(1 for name in self.pathway_selections.keys() 
                                     if self.pathway_filtered_pathways_data.get(name, {}).get('status') == 'No Change')
            
            # Store displayed counts for stats label
            self._displayed_pathway_counts = {
                'total': len(self.pathway_selections),
                'activated': displayed_activated,
                'inhibited': displayed_inhibited,
                'no_change': displayed_no_change
            }
            logger.info(f"   📊 Displayed pathway counts: Total={len(self.pathway_selections)}, Activated={displayed_activated}, Inhibited={displayed_inhibited}, No Change={displayed_no_change}")

            # Log a concise summary that matches the current table exactly
            try:
                self._log_network(f"\n{'='*70}\n")
                self._log_network(f"✅ ENRICHMENT RESULTS (current display settings)\n")
                self._log_network(f"{'='*70}\n")
                self._log_network(f"📊 Pathways shown: {len(self.pathway_selections)}\n")
                self._log_network(f"   • Activated: {displayed_activated}\n")
                self._log_network(f"   • Inhibited: {displayed_inhibited}\n")
                self._log_network(f"   • No change: {displayed_no_change}\n")

                # Multi-omics transparency summary (based on pathways currently shown in the table)
                shown_count = len(displayed_items)
                shown_multi_source = sum(1 for it in displayed_items if int(it.get('n_sources', 0) or 0) > 1)
                shown_positive_mo = sum(1 for it in displayed_items if float(it.get('mo_score', 0.0) or 0.0) > 0.0)
                single_source_failures = [
                    it for it in displayed_items
                    if int(it.get('n_sources', 0) or 0) == 1 and float(it.get('mo_score', 0.0) or 0.0) <= 0.0
                ]
                single_source_failures = sorted(
                    single_source_failures,
                    key=lambda x: (x.get('p', 1.0), -int(x.get('n', 0) or 0), str(x.get('name', '')))
                )

                self._log_network(f"\n🔎 MULTI-OMICS TRANSPARENCY\n")
                self._log_network(f"   • pathways shown: {shown_count}\n")
                self._log_network(f"   • pathways with >1 sources: {shown_multi_source}\n")
                self._log_network(f"   • pathways with positive mo_score: {shown_positive_mo}\n")

                if single_source_failures:
                    self._log_network("   • top pathways that failed multi-omics due to single-source only:\n")
                    for it in single_source_failures[:10]:
                        self._log_network(
                            f"      - {it.get('name', '')} | n_sources=1, mo_score={float(it.get('mo_score', 0.0) or 0.0):.3f}, "
                            f"#compounds={int(it.get('n', 0) or 0)}, p={float(it.get('p', 1.0) or 1.0):.2e}\n"
                        )
                else:
                    self._log_network("   • top pathways that failed multi-omics due to single-source only: none\n")

                # Compact run fingerprint for side-by-side reproducibility checks.
                try:
                    src_info = getattr(self, '_last_multiomics_source_info', {}) or {}
                    source_col_used = src_info.get('source_col') or 'None'
                    unique_sources = src_info.get('unique_sources', []) or []
                    n_sources_labels = len(unique_sources)
                    source_labels_preview = ', '.join([str(s) for s in unique_sources[:6]])
                    if len(unique_sources) > 6:
                        source_labels_preview += ', ...'

                    self._log_network("\n🧾 RUN FINGERPRINT\n")
                    self._log_network(
                        f"   source_col={source_col_used} | source_labels={n_sources_labels} [{source_labels_preview}] | "
                        f"pathway_p={pvalue_threshold} | max_pathways={max_pathways if limit_enabled else 'unlimited'} | "
                        f"shown={shown_count} | gt1_sources={shown_multi_source} | mo_positive={shown_positive_mo}\n"
                    )
                except Exception:
                    pass

                self._log_network(f"{'='*70}\n\n")
            except Exception:
                pass
            
        except Exception as e:
            logger.error(f"Failed to populate pathway tree: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def _is_valid_pathway(self, pathway_name):
        """Validate that a pathway name is legitimate (not a fake entry like '5 1')
        
        A valid pathway should:
        - Not be purely numeric
        - Not be a malformed entry from improper splitting
        - Have reasonable content (letters/meaningful chars)
        """
        if not pathway_name or not isinstance(pathway_name, str):
            return False
        
        pathway_name = pathway_name.strip()
        # Old-GUI compatible exclusions
        if self._should_exclude_pathway_name(pathway_name):
            return False
        
        # Skip if too short or just numbers/spaces
        if len(pathway_name) < 3:
            return False
        
        # Skip if only numbers (like '5', '1', '5 1', etc.)
        if pathway_name.replace(' ', '').isdigit():
            return False
        
        # Skip if it looks like garbage (mostly numbers)
        digit_ratio = sum(1 for c in pathway_name if c.isdigit()) / len(pathway_name)
        if digit_ratio > 0.7:
            return False
        
        return True
    
    def _filter_valid_pathways(self, pathways_input):
        """Split pathways and filter out invalid/unwanted entries (old GUI rules)
        
        Parameters:
        -----------
        pathways_input : str | list
            Either a semicolon/pipe-separated pathway string or a list of names
            
        Returns:
        --------
        list
            List of valid pathway names
        """
        if pathways_input is None:
            return []

        # Build a flat list of candidate names
        candidates = []
        if isinstance(pathways_input, str):
            candidates = self._split_pathways(pathways_input)
        elif isinstance(pathways_input, (list, tuple, set)):
            for item in pathways_input:
                if isinstance(item, str):
                    candidates.extend(self._split_pathways(item))
                else:
                    # best-effort string conversion
                    s = str(item)
                    candidates.extend(self._split_pathways(s))
        else:
            # Fallback: treat as string
            candidates = self._split_pathways(str(pathways_input))

        valid_pathways = []
        for pathway in candidates:
            p = self._normalize_space_case(pathway)
            if self._is_valid_pathway(p):
                valid_pathways.append(p)
        
        # Deduplicate case-insensitively
        seen = set()
        result = []
        for p in valid_pathways:
            key = p.lower()
            if key not in seen:
                seen.add(key)
                result.append(p)
        return result
    
    def _populate_metabolite_tree(self):
        """Populate the metabolite tree table with metabolite data"""
        try:
            if not hasattr(self, 'metabolite_tree') or self.metabolite_tree is None:
                logger.warning("Metabolite tree widget not initialized yet")
                return
            
            # Clear existing data
            for item in self.metabolite_tree.get_children():
                self.metabolite_tree.delete(item)
            
            if self.pathway_filtered_metabolites_data is None or len(self.pathway_filtered_metabolites_data) == 0:
                self._update_component_labels()
                return
            
            df = self.pathway_filtered_metabolites_data
            
            # Get p-value threshold from network settings
            try:
                p_threshold = float(self.network_p_threshold.get()) if hasattr(self, 'network_p_threshold') else 0.05
            except (ValueError, AttributeError):
                p_threshold = 0.05
            
            # CRITICAL FIX: When class compression is active, don't filter by p-value on class nodes
            # Class nodes have averaged p-values which are usually >0.05
            # Instead, show class nodes that appear in relationship mappings (contain significant metabolites)
            class_compression_active = getattr(self, 'class_compression_active', False)
            
            # Populate metabolite/compounds table
            # Columns: ('Compounds', 'P-Value', 'log2FC', 'Pathways', 'Upstream', 'Diseases', 'Regulation')
            total_metabolites = 0  # Count of compounds displayed (significant + connections)
            filtered_not_significant = 0  # Count filtered due to p-value > threshold
            filtered_no_connections = 0   # Count filtered due to zero connections
            total_in_df = len(df)  # Total compounds in source data
            
            for _, row in df.iterrows():
                name = self._get_metabolite_value(row, 'name', 'Unknown')
                pvalue = self._get_metabolite_value(row, 'pvalue', 1.0)
                log2fc = self._get_metabolite_value(row, 'log2fc', 0.0)

                try:
                    pvalue = float(pvalue)  # type: ignore
                except (TypeError, ValueError):
                    pvalue = 1.0

                try:
                    log2fc = float(log2fc)  # type: ignore
                except (TypeError, ValueError):
                    log2fc = 0.0
                
                # FILTER 1: Apply p-value filter ONLY when class compression is NOT active
                # When compression is active, class nodes represent groups of metabolites
                # and we should show them if they have connections
                if not class_compression_active:
                    if pvalue > p_threshold:
                        filtered_not_significant += 1
                        continue  # Skip non-significant metabolites
                
                # Count pathways using relationship mapping (NOT filtered table)
                pathways_count = 0
                if hasattr(self, 'metabolite_to_pathways') and name in self.metabolite_to_pathways:
                    pathways_count = len(self.metabolite_to_pathways.get(name, set()))
                
                # Count upstream using relationship mapping (NOT displayed tree)
                upstream_count = 0
                if hasattr(self, 'metabolite_to_upstream') and name in self.metabolite_to_upstream:
                    upstream_count = len(self.metabolite_to_upstream.get(name, set()))
                
                # Count diseases using relationship mapping (NOT displayed tree)
                diseases_count = 0
                if hasattr(self, 'metabolite_to_diseases') and name in self.metabolite_to_diseases:
                    diseases_count = len(self.metabolite_to_diseases.get(name, set()))
                
                # FILTER 2: Must have at least one connection (pathway/upstream/disease)
                total_connections = pathways_count + upstream_count + diseases_count
                if total_connections == 0:
                    filtered_no_connections += 1
                    continue  # Skip compounds with no connections
                
                # Determine regulation
                if log2fc > 0:
                    regulation = 'Up-regulated'
                elif log2fc < 0:
                    regulation = 'Down-regulated'
                else:
                    regulation = 'No Change'
                
                self.metabolite_tree.insert('', 'end', values=(
                    name,
                    f"{pvalue:.2e}",
                    f"{log2fc:.3f}",
                    pathways_count,
                    upstream_count,
                    diseases_count,
                    regulation
                ))
                total_metabolites += 1
            
            # Verify actual tree row count matches our counter
            actual_tree_rows = len(self.metabolite_tree.get_children())
            
            comp_label = self._component_label(plural=True)
            total_filtered = filtered_not_significant + filtered_no_connections
            logger.info(f"✅ {comp_label} table populated: {total_metabolites} displayed")
            if class_compression_active:
                logger.info(f"   📊 Filtered out: {filtered_no_connections} (no connections) = {total_filtered} hidden (p-value filter skipped for class nodes)")
            else:
                logger.info(f"   📊 Filtered out: {filtered_not_significant} (p>{p_threshold}) + {filtered_no_connections} (no connections) = {total_filtered} hidden")
            logger.info(f"   🔍 Verification: Tree widget contains {actual_tree_rows} actual rows, {total_in_df} total in source")
            
            if actual_tree_rows != total_metabolites:
                logger.warning(f"   ⚠️ MISMATCH: Expected {total_metabolites} rows but tree has {actual_tree_rows} rows!")
            
            # Update stats label if it exists
            if hasattr(self, 'metabolite_stats_label') and self.metabolite_stats_label:
                if class_compression_active:
                    new_text = f"Shown: {total_metabolites} | Hidden: {total_filtered} (no links: {filtered_no_connections}) | Total: {total_in_df}"
                else:
                    new_text = f"Shown: {total_metabolites} | Hidden: {total_filtered} (p>{p_threshold}: {filtered_not_significant}, no links: {filtered_no_connections}) | Total: {total_in_df}"
                self.metabolite_stats_label.config(text=new_text)
                self.metabolite_stats_label.update_idletasks()  # Force immediate UI refresh
                logger.info(f"   📊 Stats label updated to: {new_text}")
            
        except Exception as e:
            logger.error(f"Failed to populate metabolite tree: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def _populate_upstream_tree(self):
        """Populate the upstream regulators tree (enzymes + transporters)"""
        try:
            if not hasattr(self, 'upstream_tree') or self.upstream_tree is None:
                logger.warning("Upstream tree widget not initialized yet")
                return
            
            # Clear existing data
            for item in self.upstream_tree.get_children():
                self.upstream_tree.delete(item)
            
            if self.pathway_filtered_metabolites_data is None or len(self.pathway_filtered_metabolites_data) == 0:
                return
            
            # Use the filtered metabolites data directly (already at class level if compression is active)
            df = self.pathway_filtered_metabolites_data

            # IMPORTANT: keep upstream table consistent with the Metabolites tab
            # Only metabolites that pass the user-defined p-value threshold
            # should contribute to upstream counts.
            try:
                p_threshold = float(self.network_p_threshold.get()) if hasattr(self, 'network_p_threshold') else 0.05
            except (ValueError, AttributeError):
                p_threshold = 0.05

            # When class compression is active, don't filter by p-value (already filtered during compression)
            class_compression_active = getattr(self, 'class_compression_active', False)
            if not class_compression_active and 'pvalue' in df.columns:
                try:
                    numeric_p = pd.to_numeric(df['pvalue'], errors='coerce')
                    df = df.loc[numeric_p <= p_threshold].copy()
                except Exception:
                    # If anything goes wrong, fall back to unfiltered df
                    pass
            
            self.upstream_selections = {}
            
            # Collect enzymes and transporters
            upstream_dict = {}  # {name: {'type': type, 'metabolites': [list]}}
            
            # Helper function to extract gene name from enzyme/transporter entries
            def extract_gene_name(entry_str):
                """
                Extract gene name from enzyme/transporter entries.
                Format examples:
                  - "Name;Gene;ID" -> extract "Gene" (middle field)
                  - "Gene" -> use as-is
                  - Multiple entries separated by | or semicolon at top level
                """
                if not entry_str or not isinstance(entry_str, str):
                    return []
                
                entry_str = entry_str.strip()
                if not entry_str or entry_str.lower() == 'nan':
                    return []
                
                # First split by pipe (|) for multiple entries
                entries = entry_str.split('|')
                genes = []
                
                for entry in entries:
                    entry = entry.strip()
                    if not entry or entry.lower() == 'nan':
                        continue
                    
                    # Check if entry contains semicolons (Name;Gene;ID format)
                    if ';' in entry:
                        parts = entry.split(';')
                        if len(parts) >= 2:
                            # Extract gene name (second field)
                            gene = parts[1].strip()
                            if gene and gene.lower() != 'nan':
                                genes.append(gene)
                    else:
                        # No semicolons, use the whole entry as gene name
                        if entry:
                            genes.append(entry)
                
                return genes
            
            # Process enzymes
            if 'Enzyme_Gene_name' in df.columns:
                for _, row in df.iterrows():
                    metabolite_name = self._get_metabolite_value(row, 'name', 'Unknown')
                    enzymes_str = str(row.get('Enzyme_Gene_name', ''))
                    
                    gene_names = extract_gene_name(enzymes_str)
                    for gene in gene_names:
                        if gene not in upstream_dict:
                            upstream_dict[gene] = {'type': 'Enzyme', 'metabolites': []}
                        upstream_dict[gene]['metabolites'].append(metabolite_name)
            
            # Process transporters (use Gene name as display key, old GUI behavior)
            # Prefer 'Transporter_Gene_name'. Fallback to 'Transporters' or 'Transporter' only if gene is missing.
            transporter_gene_available = 'Transporter_Gene_name' in df.columns
            if transporter_gene_available:
                for _, row in df.iterrows():
                    metabolite_name = self._get_metabolite_value(row, 'name', 'Unknown')
                    transporters_str = str(row.get('Transporter_Gene_name', ''))
                    
                    gene_names = extract_gene_name(transporters_str)
                    for gene in gene_names:
                        if gene not in upstream_dict:
                            upstream_dict[gene] = {'type': 'Transporter', 'metabolites': []}
                        else:
                            # If enzyme already exists, mark as both
                            if upstream_dict[gene]['type'] == 'Enzyme':
                                upstream_dict[gene]['type'] = 'Enzyme & Transporter'
                        upstream_dict[gene]['metabolites'].append(metabolite_name)
            else:
                # Fallback to names if no gene column exists
                for tcol in ['Transporters', 'Transporter']:
                    if tcol in df.columns:
                        for _, row in df.iterrows():
                            metabolite_name = self._get_metabolite_value(row, 'name', 'Unknown')
                            tval = row.get(tcol, '')
                            if isinstance(tval, str) and tval and tval.strip() and str(tval).lower() != 'nan':
                                for transporter in str(tval).split('|'):
                                    transporter = transporter.strip()
                                    if transporter and transporter.lower() != 'nan':
                                        if transporter not in upstream_dict:
                                            upstream_dict[transporter] = {'type': 'Transporter', 'metabolites': []}
                                        upstream_dict[transporter]['metabolites'].append(metabolite_name)
            
            # Store original data for filtering
            self.original_upstream_data = upstream_dict.copy()
            
            # Get min_metabolites filter from upstream filter controls
            try:
                min_metabolites_upstream = int(self.upstream_min_metabolites.get())
            except (ValueError, AttributeError):
                min_metabolites_upstream = 1  # Default
            
            # Apply min_metabolites filter
            filtered_upstream_dict = {}
            for name, info in upstream_dict.items():
                n_metabolites = len(set(info['metabolites']))
                if n_metabolites >= min_metabolites_upstream:
                    filtered_upstream_dict[name] = info
            
            logger.info(f"🔍 Upstream filter: min_metabolites={min_metabolites_upstream}")
            logger.info(f"   Before filter: {len(upstream_dict)} regulators")
            logger.info(f"   After filter: {len(filtered_upstream_dict)} regulators")
            
            # Use filtered dict for display
            upstream_dict = filtered_upstream_dict
            
            # CRITICAL FIX: Also build upstream_to_metabolites mapping from the local upstream_dict
            # This ensures cascading works even if enzyme columns are missing from metabolites_df
            # Clear the mapping first
            self.upstream_to_metabolites = {}
            
            # Build reverse mapping: upstream -> metabolites
            for upstream_name, info in upstream_dict.items():
                metabolite_set = set(info['metabolites'])
                self.upstream_to_metabolites[upstream_name] = metabolite_set
            
            # Also update metabolite_to_upstream mapping to match
            # Clear and rebuild from upstream_dict to ensure consistency
            for metabolite_name in set():
                if metabolite_name in self.metabolite_to_upstream:
                    del self.metabolite_to_upstream[metabolite_name]
            
            # Rebuild metabolite_to_upstream from upstream_dict
            for upstream_name, info in upstream_dict.items():
                for metabolite_name in info['metabolites']:
                    if metabolite_name not in self.metabolite_to_upstream:
                        self.metabolite_to_upstream[metabolite_name] = set()
                    self.metabolite_to_upstream[metabolite_name].add(upstream_name)
            
            # Debug: log what we found
            if upstream_dict:
                logger.info(f"📍 _populate_upstream_tree found {len(upstream_dict)} upstream regulators")
                logger.info(f"   Sample: {sorted(list(upstream_dict.keys())[:5])}")
                logger.info(f"   Updated upstream_to_metabolites with {len(self.upstream_to_metabolites)} entries")
            else:
                logger.warning("   ⚠️ _populate_upstream_tree found NO upstream regulators!")
            
            # Compute ML metrics for upstream regulators
            try:
                from main_script.ml_pathway_models import compute_upstream_diffusion_scores, compute_upstream_centrality

                # Prefer filtered metabolites if available; fall back to df
                ml_df = self.pathway_filtered_metabolites_data
                if ml_df is None or getattr(ml_df, 'empty', True):
                    ml_df = df

                ml_series = compute_upstream_diffusion_scores(ml_df, self.upstream_to_metabolites)
                cent_df = compute_upstream_centrality(self.upstream_to_metabolites)
                # Fix: Check if cent_df is actually a DataFrame before using it
                if not isinstance(cent_df, pd.DataFrame) or cent_df.empty:
                    raise ValueError("Centrality computation returned invalid data")
            except Exception as _e:
                logger.warning(f"Upstream ML augmentation skipped: {_e}")
                import pandas as _pd
                ml_series = _pd.Series(dtype=float)
                cent_df = _pd.DataFrame(columns=['pagerank', 'betweenness'])

            # Compute metabolite→source mapping once for multi-omics scoring
            metabolite_df = self.pathway_filtered_metabolites_data if self.pathway_filtered_metabolites_data is not None else df
            if isinstance(metabolite_df, pd.DataFrame) and not metabolite_df.empty:
                metabolite_sources = self._compute_metabolite_sources(metabolite_df)
            else:
                metabolite_sources = {}

            items = []
            for name, info in upstream_dict.items():
                mets = set(info['metabolites'])
                n_met = len(mets)
                # Multi-omics score for this upstream regulator based on contributing metabolites
                mo_score, n_sources, total_hits = self._compute_multiomics_score(mets, metabolite_sources)
                ml = float(ml_series.get(name, 0.0)) if not ml_series.empty else 0.0
                pr = float(cent_df.loc[name, 'pagerank']) if name in getattr(cent_df, 'index', []) else 0.0  # type: ignore
                btw = float(cent_df.loc[name, 'betweenness']) if name in getattr(cent_df, 'index', []) else 0.0  # type: ignore
                items.append({
                    'name': name,
                    'type': info['type'],
                    'n': n_met,
                    'mo_score': mo_score,
                    'n_sources': n_sources,
                    'ml': ml,
                    'pr': pr,
                    'btw': btw,
                })

            # Ranking: compute ranks for consensus
            count_sorted = sorted(items, key=lambda x: x['n'], reverse=True)
            ml_sorted = sorted(items, key=lambda x: x['ml'], reverse=True)
            rank_n = {it['name']: idx + 1 for idx, it in enumerate(count_sorted)}
            rank_ml = {it['name']: idx + 1 for idx, it in enumerate(ml_sorted)}

            # Get ML filter thresholds for pre-filtering
            try:
                ml_score_min = float(self.upstream_ml_score_min.get())
            except (ValueError, AttributeError):
                ml_score_min = 0.0
            
            try:
                max_items = int(self.upstream_consensus_rank_max.get())
            except (ValueError, AttributeError):
                max_items = 999999
            
            # Apply ML score filter first
            ml_filtered_items = [it for it in items if it['ml'] >= ml_score_min]

            # Determine sorting mode
            try:
                mode = self.upstream_rank_mode.get() if self.upstream_rank_mode else 'Default'
            except Exception:
                mode = 'Default'

            if mode == 'ML':
                # Prioritize stronger multi-omics support, then ML score and count
                items_sorted = sorted(ml_filtered_items, key=lambda x: (-x['mo_score'], -x['ml'], -x['n'], x['name'].lower()))
            elif mode == 'Consensus':
                tmp = []
                for it in ml_filtered_items:
                    c = (rank_n.get(it['name'], 0) + rank_ml.get(it['name'], 0)) / 2.0 if ml_filtered_items else 0
                    jt = dict(it)
                    jt['consensus'] = c
                    tmp.append(jt)
                items_sorted = sorted(tmp, key=lambda x: (x['consensus'], -x['mo_score'], -x['ml'], -x['n'], x['name'].lower()))
            else:
                # Default: count first, but still nudge by multi-omics score
                items_sorted = sorted(ml_filtered_items, key=lambda x: (-x['n'], -x['mo_score'], -x['ml'], x['name'].lower()))
            
            # Check if Apply Auto Filter is checked
            apply_auto_filter = self.upstream_apply_filter.get() if hasattr(self, 'upstream_apply_filter') else False
            
            # If Auto Filter is ON: limit to top 25
            # If Auto Filter is OFF: use GUI settings with no top limit
            if apply_auto_filter:
                items_sorted = items_sorted[:25]
                logger.info(f"   Auto filter ON: Limiting to top 25 upstream regulators")
            else:
                logger.info(f"   Auto filter OFF: Using GUI settings (showing all {len(items_sorted)} items passing filters)")
            
            # Recompute multi-omics ranks based on displayed items only
            mo_sorted = sorted(items_sorted, key=lambda x: x['mo_score'], reverse=True)
            mo_rank = {}
            rank_counter = 1
            for it in mo_sorted:
                if it['mo_score'] <= 0:
                    continue
                mo_rank[it['name']] = rank_counter
                rank_counter += 1
            
            # Calculate consensus for all items (needed for tree display)
            for it in items_sorted:
                if 'consensus' not in it:
                    # Calculate consensus rank as average of count and ML ranks
                    c = (rank_n.get(it['name'], 0) + rank_ml.get(it['name'], 0)) / 2.0
                    it['consensus'] = c
            
            # Populate tree
            for it in items_sorted:
                
                composite_key = f"{it['type']}_{it['name']}"
                self.upstream_selections[composite_key] = False
                select_symbol = "☐"

                # Only display a multi-omics rank if present
                mo_rank_display = mo_rank.get(it['name'], 0)
                mo_rank_display = mo_rank_display if mo_rank_display > 0 else ''

                self.upstream_tree.insert('', 'end', values=(
                    select_symbol,
                    it['name'],
                    it['type'],
                    it['n'],
                    mo_rank_display,
                    f"{it['ml']:.4f}",
                    f"{it['pr']:.4f}",
                    f"{it['btw']:.4f}",
                    f"{it.get('consensus', 0):.2f}"
                ))
            
            displayed = len(self.upstream_tree.get_children())
            total = len(items)
            filtered = total - displayed
            logger.info(f"✅ Upstream tree populated with {displayed} items shown, {filtered} filtered (all UNSELECTED)")
            if hasattr(self, 'upstream_stats_label') and self.upstream_stats_label:
                self.upstream_stats_label.config(text=f"Showing: {displayed} | Total: {total}")
            
        except Exception as e:
            logger.error(f"Failed to populate upstream tree: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def _populate_disease_tree(self):
        """Populate the diseases tree table with disease data"""
        try:
            if not hasattr(self, 'disease_tree') or self.disease_tree is None:
                logger.warning("Disease tree widget not initialized yet")
                return
            
            # Clear existing data
            for item in self.disease_tree.get_children():
                self.disease_tree.delete(item)
            
            # Use filtered metabolite data (significant compounds only)
            # CRITICAL: Disease table should show counts based on SIGNIFICANT metabolites only
            if self.pathway_filtered_metabolites_data is None or len(self.pathway_filtered_metabolites_data) == 0:
                logger.warning("Disease table: No filtered metabolite data available")
                return
            
            df = self.pathway_filtered_metabolites_data
            self.disease_selections = {}
            
            # DEBUG: Check if Associated_Diseases column exists
            logger.info(f"[DISEASE DEBUG] Checking disease data (using FILTERED significant metabolites)...")
            logger.info(f"[DISEASE DEBUG] DataFrame shape: {df.shape}")
            logger.info(f"[DISEASE DEBUG] Has 'Associated_Diseases' column: {'Associated_Diseases' in df.columns}")
            
            # Determine which column holds disease annotations
            disease_col = None
            if 'Associated_Diseases' in df.columns:
                disease_col = 'Associated_Diseases'
            elif 'Associated' in df.columns:
                disease_col = 'Associated'

            if disease_col:
                # Check how many non-empty disease entries we have
                non_empty = df[disease_col].fillna('').astype(str).str.strip()
                non_empty = non_empty[non_empty != '']
                non_empty = non_empty[non_empty.str.lower() != 'nan']
                logger.info(f"[DISEASE DEBUG] Rows with non-empty Associated_Diseases: {len(non_empty)}/{len(df)}")
                if len(non_empty) > 0:
                    logger.info(f"[DISEASE DEBUG] Sample entries (first 3):")
                    for i, val in enumerate(non_empty.head(3), 1):
                        logger.info(f"  {i}. {val[:100]}")
            else:
                logger.error(f"[DISEASE DEBUG] Disease annotation column is missing ('Associated_Diseases' or 'Associated')!")
                logger.info(f"[DISEASE DEBUG] Available columns: {list(df.columns)}")
            
            # Collect diseases from metabolites
            # NOTE: Disease reclassification happens in annotation tab BEFORE Fisher ORA
            # Associated_Diseases column already contains ONLY true diseases (cleaned)
            # CRITICAL: Using pathway_filtered_metabolites_data which contains ONLY significant metabolites
            # No need for additional p-value filtering here since df is already filtered
            
            diseases_dict = {}  # {disease_name: [metabolites]}
            
            if disease_col:
                for _, row in df.iterrows():
                    metabolite_name = self._get_metabolite_value(row, 'name', 'Unknown')
                    diseases = str(row.get(disease_col, ''))
                    
                    if diseases and diseases.strip() and diseases.lower() != 'nan':
                        for disease in diseases.split('|'):
                            disease = disease.strip()
                            # Skip empty, 'nan', and whitespace-only entries
                            if disease and disease.lower() != 'nan':
                                if disease not in diseases_dict:
                                    diseases_dict[disease] = []
                                diseases_dict[disease].append(metabolite_name)
                
                # DEBUG: Show compound count for a few diseases
                logger.info(f"[DISEASE DEBUG] Total significant metabolites in filtered data: {len(df)}")
                if diseases_dict:
                    logger.info(f"[DISEASE DEBUG] Sample disease compound counts (from SIGNIFICANT metabolites only):")
                    for i, (disease_name, mets) in enumerate(list(diseases_dict.items())[:5], 1):
                        unique_mets = len(set(mets))
                        logger.info(f"   {i}. {disease_name}: {unique_mets} significant compound(s)")
            
            # Store original data for filtering
            self.original_disease_data = diseases_dict.copy()
            
            # Get min_metabolites filter from disease filter controls
            try:
                min_metabolites_disease = int(self.disease_min_metabolites.get())
            except (ValueError, AttributeError):
                min_metabolites_disease = 1  # Default
            
            # Apply min_metabolites filter
            filtered_diseases_dict = {}
            for disease_name, metabolites_list in diseases_dict.items():
                n_metabolites = len(set(metabolites_list))
                if n_metabolites >= min_metabolites_disease:
                    filtered_diseases_dict[disease_name] = metabolites_list
            
            logger.info(f"🔍 Disease filter: min_metabolites={min_metabolites_disease}")
            logger.info(f"   Before filter: {len(diseases_dict)} diseases")
            logger.info(f"   After filter: {len(filtered_diseases_dict)} diseases")
            
            # Use filtered dict for display
            diseases_dict = filtered_diseases_dict
            
            # CRITICAL FIX: Also build disease_to_metabolites and metabolite_to_diseases mappings
            # This ensures cascading works properly
            self.disease_to_metabolites = {}
            self.metabolite_to_diseases = {}
            
            # Build mappings from diseases_dict
            for disease_name, metabolites_list in diseases_dict.items():
                metabolite_set = set(metabolites_list)
                self.disease_to_metabolites[disease_name] = metabolite_set
                
                # Also update metabolite_to_diseases
                for metabolite_name in metabolite_set:
                    if metabolite_name not in self.metabolite_to_diseases:
                        self.metabolite_to_diseases[metabolite_name] = set()
                    self.metabolite_to_diseases[metabolite_name].add(disease_name)
            
            # Compute ML metrics and consensus for diseases
            try:
                from main_script.ml_pathway_models import compute_disease_diffusion_scores, compute_disease_centrality
                ml_series = compute_disease_diffusion_scores(df, diseases_dict)
                cent_df = compute_disease_centrality(diseases_dict)
            except Exception as _e:
                logger.warning(f"Disease ML augmentation skipped: {_e}")
                import pandas as _pd
                ml_series = _pd.Series(dtype=float)
                cent_df = _pd.DataFrame(columns=['pagerank', 'betweenness'])

            # Build items and ranks
            # Also compute multi-omics scores per disease when Source information exists
            # Prefer filtered metabolites if available, otherwise fall back to full df
            metabolite_df = self.pathway_filtered_metabolites_data
            if metabolite_df is None or getattr(metabolite_df, 'empty', True):
                metabolite_df = df
            metabolite_sources = self._compute_metabolite_sources(metabolite_df)

            items = []
            for disease_name, metabolites in diseases_dict.items():
                n_metabolites = len(set(metabolites))
                mo_score, n_sources, total_hits = self._compute_multiomics_score(metabolites, metabolite_sources)
                ml = float(ml_series.get(disease_name, 0.0)) if not ml_series.empty else 0.0
                pr = float(cent_df.loc[disease_name, 'pagerank']) if disease_name in getattr(cent_df, 'index', []) else 0.0  # type: ignore
                btw = float(cent_df.loc[disease_name, 'betweenness']) if disease_name in getattr(cent_df, 'index', []) else 0.0  # type: ignore
                items.append({
                    'name': disease_name,
                    'n': n_metabolites,
                    'mo_score': mo_score,
                    'n_sources': n_sources,
                    'ml': ml,
                    'pr': pr,
                    'btw': btw,
                })

            # Ranks: by count desc and ML desc
            count_sorted = sorted(items, key=lambda x: x['n'], reverse=True)
            ml_sorted = sorted(items, key=lambda x: x['ml'], reverse=True)
            rank_n = {it['name']: idx + 1 for idx, it in enumerate(count_sorted)}
            rank_ml = {it['name']: idx + 1 for idx, it in enumerate(ml_sorted)}

            # Get ML filter thresholds for pre-filtering
            try:
                ml_score_min = float(self.disease_ml_score_min.get())
            except (ValueError, AttributeError):
                ml_score_min = 0.0
            
            try:
                max_items = int(self.disease_consensus_rank_max.get())
            except (ValueError, AttributeError):
                max_items = 999999
            
            # Apply ML score filter first
            ml_filtered_items = [it for it in items if it['ml'] >= ml_score_min]

            # Ranking mode selection
            try:
                mode = self.disease_rank_mode.get() if self.disease_rank_mode else 'Default'
            except Exception:
                mode = 'Default'

            if mode == 'ML':
                # Prioritize stronger multi-omics support, then ML score and count
                items_sorted = sorted(ml_filtered_items, key=lambda x: (-x['mo_score'], -x['ml'], -x['n'], x['name'].lower()))
            elif mode == 'Consensus':
                tmp = []
                for it in ml_filtered_items:
                    name = it['name']
                    c = (rank_n.get(name, 0) + rank_ml.get(name, 0)) / 2.0 if ml_filtered_items else 0
                    jt = dict(it)
                    jt['consensus'] = c
                    tmp.append(jt)
                items_sorted = sorted(tmp, key=lambda x: (x['consensus'], -x['mo_score'], -x['ml'], -x['n'], x['name'].lower()))
            else:
                # Default order: count desc, then ML desc
                items_sorted = sorted(ml_filtered_items, key=lambda x: (-x['n'], -x['mo_score'], -x['ml'], x['name'].lower()))

            # Check if Apply Auto Filter is checked
            apply_auto_filter = self.disease_apply_filter.get() if hasattr(self, 'disease_apply_filter') else False
            
            # If Auto Filter is ON: limit to top 25
            # If Auto Filter is OFF: use GUI settings with no top limit
            if apply_auto_filter:
                items_sorted = items_sorted[:25]
                logger.info(f"   Auto filter ON: Limiting to top 25 diseases")
            else:
                logger.info(f"   Auto filter OFF: Using GUI settings (showing all {len(items_sorted)} items passing filters)")
            
            # Recompute multi-omics ranks based on displayed items only
            mo_sorted = sorted(items_sorted, key=lambda x: x['mo_score'], reverse=True)
            mo_rank = {}
            rank_counter = 1
            for it in mo_sorted:
                if it['mo_score'] <= 0:
                    continue
                mo_rank[it['name']] = rank_counter
                rank_counter += 1

            # Calculate consensus for all items (needed for tree display)
            for it in items_sorted:
                if 'consensus' not in it:
                    # Calculate consensus rank as average of count and ML ranks
                    c = (rank_n.get(it['name'], 0) + rank_ml.get(it['name'], 0)) / 2.0
                    it['consensus'] = c

            # Populate tree
            for it in items_sorted:
                name = it['name']
                
                # Mark as UNSELECTED by default
                self.disease_selections[name] = False
                select_symbol = "☐"

                # Only display a multi-omics rank if present
                mo_rank_display = mo_rank.get(name, 0)
                mo_rank_display = mo_rank_display if mo_rank_display > 0 else ''

                self.disease_tree.insert('', 'end', values=(
                    select_symbol,
                    name,
                    it['n'],
                    mo_rank_display,
                    f"{it['ml']:.4f}",
                    f"{it['pr']:.4f}",
                    f"{it['btw']:.4f}",
                    f"{it.get('consensus', 0):.2f}"
                ))
            
            displayed = len(self.disease_tree.get_children())
            total = len(items)
            filtered = total - displayed
            logger.info(f"✅ Disease tree populated with {displayed} diseases shown, {filtered} filtered (all UNSELECTED)")
            logger.info(f"   Updated disease_to_metabolites with {len(self.disease_to_metabolites)} entries")
            logger.info(f"   Updated metabolite_to_diseases for {len(self.metabolite_to_diseases)} metabolites")
            if hasattr(self, 'diseases_stats_label') and self.diseases_stats_label:
                self.diseases_stats_label.config(text=f"Showing: {displayed} | Total: {total}")
            
        except Exception as e:
            logger.error(f"Failed to populate disease tree: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def _update_stats_labels(self):
        """Update statistics labels with current data counts"""
        try:
            # Update pathway stats - show BOTH displayed counts and total significant pathways
            if hasattr(self, 'pathway_stats_label') and self.pathway_stats_label:
                # Use pre-computed displayed counts if available
                if hasattr(self, '_displayed_pathway_counts') and self._displayed_pathway_counts:
                    counts = self._displayed_pathway_counts
                    pathway_count = counts.get('total', 0)
                    activated = counts.get('activated', 0)
                    inhibited = counts.get('inhibited', 0)
                    no_change = counts.get('no_change', 0)
                else:
                    # Fallback: count from pathway_selections (displayed pathways)
                    pathway_count = len(getattr(self, 'pathway_selections', {}))
                    activated = sum(1 for name in getattr(self, 'pathway_selections', {}).keys() 
                                   if (getattr(self, 'pathway_filtered_pathways_data', {}) or {}).get(name, {}).get('status') == 'Activated')
                    inhibited = sum(1 for name in getattr(self, 'pathway_selections', {}).keys() 
                                   if (getattr(self, 'pathway_filtered_pathways_data', {}) or {}).get(name, {}).get('status') == 'Inhibited')
                    no_change = sum(1 for name in getattr(self, 'pathway_selections', {}).keys() 
                                   if (getattr(self, 'pathway_filtered_pathways_data', {}) or {}).get(name, {}).get('status') == 'No Change')

                # Totals from enrichment results (should match "Results Summary" in log)
                all_total = 0
                all_activated = 0
                all_inhibited = 0
                all_no_change = 0
                pathways_dict = getattr(self, 'pathway_filtered_pathways_data', None)
                if isinstance(pathways_dict, dict) and pathways_dict:
                    all_total = len(pathways_dict)
                    for stats in pathways_dict.values():
                        if not isinstance(stats, dict):
                            continue
                        status = stats.get('status', 'No Change')
                        if status == 'Activated':
                            all_activated += 1
                        elif status == 'Inhibited':
                            all_inhibited += 1
                        elif status == 'No Change':
                            all_no_change += 1

                # Stats label: emphasize what is currently visible
                # - "Shown" counts are the entries actually present in the table
                # - "Total" counts reflect all significant pathways from the enrichment step
                if all_total > 0:
                    stats_text = (
                        f"Shown: {pathway_count} / {all_total}  | "
                        f"Activated: {activated} / {all_activated}  | "
                        f"Inhibited: {inhibited} / {all_inhibited}  | "
                        f"No Change: {no_change} / {all_no_change}"
                    )
                else:
                    stats_text = (
                        f"Shown: {pathway_count} | "
                        f"Activated: {activated} | Inhibited: {inhibited} | No Change: {no_change}"
                    )
                self.pathway_stats_label.config(text=stats_text)
            
            # Update metabolite stats
            # NOTE: Metabolite stats are updated by _populate_metabolite_tree() with detailed breakdown
            # We don't overwrite it here to preserve "Shown: X | Hidden: Y | Total: Z" format
            # if hasattr(self, 'metabolite_stats_label') and self.metabolite_stats_label:
            #     metabolite_count = len(self.pathway_filtered_metabolites_data) if self.pathway_filtered_metabolites_data is not None else 0
            #     self.metabolite_stats_label.config(text=f"Metabolites: {metabolite_count}")
            
            # Update upstream stats
            if hasattr(self, 'upstream_stats_label') and self.upstream_stats_label:
                upstream_count = len(getattr(self, 'upstream_selections', {}))
                self.upstream_stats_label.config(text=f"Upstream: {upstream_count}")
            
            # Update disease stats
            if hasattr(self, 'diseases_stats_label') and self.diseases_stats_label:
                disease_count = len(getattr(self, 'disease_selections', {}))
                self.diseases_stats_label.config(text=f"Diseases: {disease_count}")
                
        except Exception as e:
            logger.error(f"Failed to update stats labels: {str(e)}")

    def setup_ui(self):
        """Setup the UI for pathway network tab"""
        # Create the pathway network interface
        self.create_pathway_network_interface(self.frame)
    
    def create_pathway_network_interface(self, parent):
        """Create the Pathway Network Analysis interface with 2-column layout"""
        # Create main frame
        main_frame = tk.Frame(parent, bg='#f0f0f0')
        main_frame.pack(fill='both', expand=True)
        
        # Configure grid: left column minimal fixed width, right expands fully
        main_frame.grid_columnconfigure(0, weight=0)  # Left column - no expansion, natural width
        main_frame.grid_columnconfigure(1, weight=1)  # Right column expands to fill ALL remaining space
        main_frame.grid_rowconfigure(0, weight=1)
        
        # ============ LEFT COLUMN: CONTROLS ============
        left_frame = tk.Frame(main_frame, bg='#ecf0f1', relief='groove', borderwidth=2, width=350)
        left_frame.grid(row=0, column=0, sticky='ns', padx=(5, 0), pady=5)
        left_frame.grid_propagate(False)  # Keep fixed width
        
        # Add canvas and scrollbars for left column
        left_canvas = tk.Canvas(left_frame, bg='#ecf0f1', highlightthickness=0)
        left_scrollbar_y = ttk.Scrollbar(left_frame, orient='vertical', command=left_canvas.yview)
        left_scrollbar_x = ttk.Scrollbar(left_frame, orient='horizontal', command=left_canvas.xview)
        left_scrollable = tk.Frame(left_canvas, bg='#ecf0f1')
        
        left_scrollable.bind("<Configure>", lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all")))
        left_canvas.create_window((0, 0), window=left_scrollable, anchor='nw')
        left_canvas.configure(yscrollcommand=left_scrollbar_y.set, xscrollcommand=left_scrollbar_x.set)
        
        left_canvas.grid(row=0, column=0, sticky='nsew')
        left_scrollbar_y.grid(row=0, column=1, sticky='ns')
        left_scrollbar_x.grid(row=1, column=0, sticky='ew')
        
        left_frame.grid_rowconfigure(0, weight=1)
        left_frame.grid_columnconfigure(0, weight=1)
        
        # Enable mousewheel scrolling for left column
        def _on_mousewheel_left(event):
            try:
                left_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
            except:
                pass
        left_canvas.bind_all("<MouseWheel>", _on_mousewheel_left)
        
        # Left column content
        left_content = tk.Frame(left_scrollable, bg='#ecf0f1')
        left_content.pack(fill='both', expand=True, padx=3, pady=3)
        
        # Debug function to track widget widths after rendering
        def debug_widths():
            try:
                left_width = left_frame.winfo_width()
                left_content_width = left_content.winfo_width()
                canvas_width = left_canvas.winfo_width()
                scrollable_width = left_scrollable.winfo_width()
                
                print(f"\n=== LEFT PANEL WIDTH DEBUG ===")
                print(f"Left frame actual width: {left_width}px")
                print(f"Left canvas width: {canvas_width}px")
                print(f"Left scrollable width: {scrollable_width}px")
                print(f"Left content width: {left_content_width}px")
                print(f"Grid minsize setting: 420px")
                print(f"Overflow: {left_content_width - 420}px beyond minsize")
                
                # Check some button widths
                if hasattr(self, 'generate_network_button'):
                    btn_width = self.generate_network_button.winfo_width()
                    btn_req_width = self.generate_network_button.winfo_reqwidth()
                    print(f"\nGenerate Network button:")
                    print(f"  Actual width: {btn_width}px")
                    print(f"  Requested width: {btn_req_width}px")
                
                # Check entry widths
                children = output_frame.winfo_children()
                for child in children:
                    if isinstance(child, tk.Entry):
                        entry_width = child.winfo_width()
                        entry_req_width = child.winfo_reqwidth()
                        print(f"\nOutput folder entry:")
                        print(f"  Actual width: {entry_width}px")
                        print(f"  Requested width: {entry_req_width}px")
                        break
                
                # Check LabelFrame widths
                if hasattr(self, 'network_progress_text'):
                    text_width = self.network_progress_text.winfo_width()
                    text_req_width = self.network_progress_text.winfo_reqwidth()
                    print(f"\nProgress text widget:")
                    print(f"  Actual width: {text_width}px")
                    print(f"  Requested width: {text_req_width}px")
                
                print(f"==============================\n")
            except Exception as e:
                print(f"Width debug failed: {e}")
        
        # Schedule debug after widgets are rendered
        left_content.after(1000, debug_widths)
        
        # ============================================
        # IMPORT ANNOTATED DATA (Browse button - moved to top)
        # ============================================
        import_frame = tk.LabelFrame(left_content, text="📥 Import Annotated Data",
                                     font=('Arial', 9, 'bold'),
                                     bg='#fff3cd',  # Light yellow background
                                     fg='#856404', padx=5, pady=5)
        import_frame.pack(fill='x', anchor='w', pady=(0, 10))
        
        tk.Label(import_frame, text="Load previously annotated pathway data:",
                font=('Arial', 8), bg='#fff3cd', fg='#856404').pack(anchor='w', pady=(0, 5))
        
        import_btn_frame = tk.Frame(import_frame, bg='#fff3cd')
        import_btn_frame.pack(fill='x')
        
        tk.Button(import_btn_frame, text="📁 Browse & Import",
                 command=self.import_annotated_data,
                 bg='#856404', fg='white', font=('Arial', 8, 'bold'),
                 padx=6, pady=4).pack(fill='x')
        
        # ============================================
        # OUTPUT FOLDER (moved to second position)
        # ============================================
        output_frame = tk.LabelFrame(left_content, text="📂 Output Folder", font=('Arial', 9, 'bold'),
                                     bg='#ecf0f1', fg='#2c3e50', padx=5, pady=5)
        output_frame.pack(fill='x', anchor='w', pady=(0, 10))

        if not hasattr(self, 'network_output_folder'):
            self.network_output_folder = tk.StringVar()

        tk.Entry(output_frame, textvariable=self.network_output_folder, font=('Arial', 8),
                 state='readonly', relief='groove').pack(fill='x', pady=(0, 8))

        folder_buttons = tk.Frame(output_frame, bg='#ecf0f1')
        folder_buttons.pack(fill='x')

        tk.Button(folder_buttons, text="Select", command=self.browse_network_output_folder,
                  bg='#9b59b6', fg='white', font=('Arial', 8, 'bold'),
                  padx=4, pady=3).pack(side='left', fill='x', expand=True, padx=(0, 3))
        self.open_folder_button = tk.Button(folder_buttons, text="Open",
                                            command=self.open_results_folder,
                                            bg='#2ecc71', fg='white', font=('Arial', 8, 'bold'),
                                            padx=4, pady=3)
        self.open_folder_button.pack(side='left', fill='x', expand=True, padx=(3, 0))

        # ============================================
        # PATHWAY ANALYSIS SETTINGS (moved to third position)
        # ============================================
        settings_frame = tk.LabelFrame(left_content, text="⚙️ Pathway Analysis Settings",
                                      font=('Arial', 9, 'bold'),
                                      bg='#d5f4e6',  # Light green background
                                      fg='#2c3e50', padx=5, pady=5)
        settings_frame.pack(fill='x', anchor='w', pady=(0, 10))
        
        # Analysis Method and FDR Correction
        row1 = tk.Frame(settings_frame, bg='#d5f4e6')
        row1.pack(fill='x', pady=(0, 5))
        
        tk.Label(row1, text="Analysis Method:", font=('Arial', 8), bg='#d5f4e6').pack(side='left')
        if not hasattr(self, 'pathway_analysis_method'):
            self.pathway_analysis_method = tk.StringVar(value="Fisher ORA")
        method_combo = ttk.Combobox(row1, textvariable=self.pathway_analysis_method,
                                   values=["Fisher ORA", "IWPA"], state='readonly', font=('Arial', 8), width=12)
        method_combo.pack(side='left', padx=(5, 10))
        
        if not hasattr(self, 'pathway_apply_fdr'):
            self.pathway_apply_fdr = tk.BooleanVar(value=True)
        tk.Checkbutton(row1, text="FDR (BH)", variable=self.pathway_apply_fdr,
                      font=('Arial', 8), bg='#d5f4e6').pack(side='left')
        
        # Min Metabolites
        row2 = tk.Frame(settings_frame, bg='#d5f4e6')
        row2.pack(fill='x', pady=(0, 5))
        
        tk.Label(row2, text="Min Compounds:", font=('Arial', 8), bg='#d5f4e6').pack(side='left')
        if not hasattr(self, 'pathway_min_metabolites'):
            self.pathway_min_metabolites = tk.IntVar(value=2)
        tk.Spinbox(row2, from_=1, to=20, textvariable=self.pathway_min_metabolites,
                  font=('Arial', 8), width=8).pack(side='left', padx=(5, 0))
        
        # Compound P-value Filter
        row3 = tk.Frame(settings_frame, bg='#d5f4e6')
        row3.pack(fill='x', pady=(0, 5))
        
        tk.Label(row3, text="Compound P-value Filter:", font=('Arial', 8), bg='#d5f4e6').pack(side='left')
        shared_metabolite_pvalue_var = self._get_shared_metabolite_pvalue_var()
        tk.Entry(row3, textvariable=shared_metabolite_pvalue_var,
            font=('Arial', 8), width=10).pack(side='left', padx=(5, 0))
        
        # Compound |log2FC| Filter (absolute fold change threshold)
        row3b = tk.Frame(settings_frame, bg='#d5f4e6')
        row3b.pack(fill='x', pady=(0, 5))
        
        tk.Label(row3b, text="Compound |log2FC| Filter:", font=('Arial', 8), bg='#d5f4e6').pack(side='left')
        if not hasattr(self, 'metabolite_log2fc_filter'):
            self.metabolite_log2fc_filter = tk.StringVar(value="0")
        tk.Entry(row3b, textvariable=self.metabolite_log2fc_filter,
                font=('Arial', 8), width=10).pack(side='left', padx=(5, 0))
        
        # Pathway P-value Threshold
        row4 = tk.Frame(settings_frame, bg='#d5f4e6')
        row4.pack(fill='x', pady=(0, 5))
        
        tk.Label(row4, text="Pathway P-value Threshold:", font=('Arial', 8), bg='#d5f4e6').pack(side='left')
        if not hasattr(self, 'pathway_pvalue_threshold'):
            self.pathway_pvalue_threshold = tk.StringVar(value="0.05")
        tk.Entry(row4, textvariable=self.pathway_pvalue_threshold,
                font=('Arial', 8), width=10).pack(side='left', padx=(5, 0))
        
        # Z-score Threshold
        row5 = tk.Frame(settings_frame, bg='#d5f4e6')
        row5.pack(fill='x', pady=(0, 5))
        
        tk.Label(row5, text="Z-score Threshold (±):", font=('Arial', 8), bg='#d5f4e6').pack(side='left')
        if not hasattr(self, 'pathway_z_threshold'):
            self.pathway_z_threshold = tk.StringVar(value="0")
        tk.Entry(row5, textvariable=self.pathway_z_threshold,
                font=('Arial', 8), width=10).pack(side='left', padx=(5, 0))
        
        # Max Total Pathways
        row6 = tk.Frame(settings_frame, bg='#d5f4e6')
        row6.pack(fill='x', pady=(0, 5))
        
        tk.Label(row6, text="Max Total Pathways:", font=('Arial', 8), bg='#d5f4e6').pack(side='left')
        if not hasattr(self, 'pathway_limit_pathways'):
            self.pathway_limit_pathways = tk.BooleanVar(value=True)
        tk.Checkbutton(row6, text="Limit", variable=self.pathway_limit_pathways,
                      font=('Arial', 8), bg='#d5f4e6').pack(side='left', padx=(5, 0))
        if not hasattr(self, 'pathway_max_total'):
            self.pathway_max_total = tk.IntVar(value=25)
        tk.Spinbox(row6, from_=5, to=10000, increment=5, textvariable=self.pathway_max_total,
                  font=('Arial', 8), width=8).pack(side='left', padx=(5, 0))
        
        # IWPA Weight Mode
        row7 = tk.Frame(settings_frame, bg='#d5f4e6')
        row7.pack(fill='x', pady=(0, 5))
        
        tk.Label(row7, text="IWPA Weight Mode:", font=('Arial', 8), bg='#d5f4e6').pack(side='left')
        if not hasattr(self, 'iwpa_weight_mode'):
            self.iwpa_weight_mode = tk.StringVar(value="signed_p")
        weight_combo = ttk.Combobox(row7, textvariable=self.iwpa_weight_mode,
                                   values=["signed_p", "combined", "log2fc"],
                                   state='readonly', font=('Arial', 8), width=12)
        weight_combo.pack(side='left', padx=(5, 0))
        
        # Font Size Scale for Statistics Plots
        row8 = tk.Frame(settings_frame, bg='#d5f4e6')
        row8.pack(fill='x', pady=(0, 5))
        
        tk.Label(row8, text="Stats Plot Font Size:", font=('Arial', 8), bg='#d5f4e6').pack(side='left')
        if not hasattr(self, 'pathway_font_size_scale'):
            self.pathway_font_size_scale = tk.DoubleVar(value=1.0)
        tk.Spinbox(row8, from_=0.5, to=3.0, increment=0.1, textvariable=self.pathway_font_size_scale,
                  font=('Arial', 8), width=8).pack(side='left', padx=(5, 0))
        tk.Label(row8, text="(0.5x - 3.0x)", font=('Arial', 7), bg='#d5f4e6', fg='#555').pack(side='left', padx=(3, 0))
        
        # Reset button
        def reset_analysis_defaults():
            self.pathway_analysis_method.set("Fisher ORA")
            self.pathway_apply_fdr.set(True)
            self.pathway_min_metabolites.set(2)
            self.metabolite_pvalue_filter.set("0.05")
            self.metabolite_log2fc_filter.set("0")
            self.pathway_pvalue_threshold.set("0.05")
            self.pathway_z_threshold.set("1.0")
            self.pathway_limit_pathways.set(True)
            self.pathway_max_total.set(25)
            self.iwpa_weight_mode.set("signed_p")
            self.pathway_font_size_scale.set(1.0)
        
        tk.Button(settings_frame, text="↺ Reset to Defaults", command=reset_analysis_defaults,
                 bg='#95a5a6', fg='white', font=('Arial', 7, 'bold'), padx=4, pady=2).pack(pady=(5, 0))
        
        # ============================================
        # COLUMN VERIFICATION (Populate - moved to fourth position)
        # ============================================
        column_frame = tk.LabelFrame(left_content, text="✅ Column Verification",
                                     font=('Arial', 9, 'bold'),
                                     bg='#d1ecf1',  # Light blue background
                                     fg='#0c5460', padx=5, pady=5)
        column_frame.pack(fill='x', anchor='w', pady=(0, 10))
        
        tk.Label(column_frame, text="Verify comparison columns for network analysis:",
                font=('Arial', 8), bg='#d1ecf1', fg='#0c5460').pack(anchor='w', pady=(0, 5))
        
        verify_btn_frame = tk.Frame(column_frame, bg='#d1ecf1')
        verify_btn_frame.pack(fill='x')
        
        tk.Button(verify_btn_frame, text="🔄 Populate Pathways",
                 command=self.verify_network_columns,
                 bg='#0c5460', fg='white', font=('Arial', 8, 'bold'),
                 padx=6, pady=4).pack(fill='x')
        
        # ============================================
        # ACTIONS (Generate Network - moved to fifth position)
        # ============================================
        action_frame = tk.LabelFrame(left_content, text="🚀 Actions", font=('Arial', 9, 'bold'),
                                     bg='#ecf0f1', fg='#2c3e50', padx=5, pady=5)
        action_frame.pack(fill='x', anchor='w', pady=(0, 10))

        # Generate button
        action_buttons = tk.Frame(action_frame, bg='#ecf0f1')
        action_buttons.pack(fill='x')

        self.generate_network_button = tk.Button(action_buttons, text="🌐 Generate Network",
                                                 command=self.start_generate_network,
                                                 bg='#95a5a6', fg='white', font=('Arial', 8, 'bold'),
                                                 state='disabled',  # Disabled until columns verified
                                                 padx=6, pady=4)
        self.generate_network_button.pack(fill='x', expand=True)
        self._create_tooltip(self.generate_network_button, 
                             "Generate network (click 'Verify Columns' first)")
        
        # Output options section
        output_opts_label = tk.Label(action_frame, text="Output Options:", font=('Arial', 8, 'bold'), 
                                     bg='#ecf0f1', fg='#2c3e50')
        output_opts_label.pack(anchor='w', pady=(8, 2))
        
        output_opts = tk.Frame(action_frame, bg='#ecf0f1')
        output_opts.pack(fill='x')
        
        # ═══════════════════════════════════════════════════════════════════
        # TOP LEVEL OPTIONS (at least one must be selected)
        # ═══════════════════════════════════════════════════════════════════
        
        # Generate Network Graph (checked by default) - ALWAYS create fresh BooleanVar
        self.network_gen_network = tk.BooleanVar(value=True)
        self.network_gen_network_cb = tk.Checkbutton(output_opts, text="📊 Generate Network Graph", 
                      variable=self.network_gen_network,
                      command=self._update_network_suboptions,
                      font=('Arial', 8, 'bold'), bg='#ecf0f1', fg='#2c3e50',
                      activebackground='#ecf0f1', selectcolor='#ecf0f1')
        self.network_gen_network_cb.pack(anchor='w')
        
        # Sub-options frame for Network Graph (indented)
        self.network_subopts_frame = tk.Frame(output_opts, bg='#ecf0f1')
        self.network_subopts_frame.pack(fill='x', padx=(20, 0))  # Indent sub-options
        
        # Sub-option 1: Open Interactive Viewer (checked by default) - ALWAYS create fresh
        self.auto_launch_interactive = tk.BooleanVar(value=True)
        self.auto_launch_cb = tk.Checkbutton(self.network_subopts_frame, text="🌐 Open Interactive Viewer", 
                      variable=self.auto_launch_interactive,
                      font=('Arial', 8), bg='#ecf0f1', fg='#2c3e50',
                      activebackground='#ecf0f1', selectcolor='#ecf0f1')
        self.auto_launch_cb.pack(anchor='w')
        
        # Sub-option 2: Export Cytoscape file (UNchecked by default) - ALWAYS create fresh
        self.export_cytoscape = tk.BooleanVar(value=False)  # Default OFF per user request
        self.export_cytoscape_cb = tk.Checkbutton(self.network_subopts_frame, text="💾 Export Cytoscape File (.graphml)", 
                      variable=self.export_cytoscape,
                      font=('Arial', 8), bg='#ecf0f1', fg='#2c3e50',
                      activebackground='#ecf0f1', selectcolor='#ecf0f1')
        self.export_cytoscape_cb.pack(anchor='w')
        
        # Generate Statistics (checked by default) - ALWAYS create fresh
        self.network_gen_plots = tk.BooleanVar(value=True)
        tk.Checkbutton(output_opts, text="📈 Generate Statistics",
                      variable=self.network_gen_plots,
                      font=('Arial', 8, 'bold'), bg='#ecf0f1', fg='#2c3e50',
                      activebackground='#ecf0f1', selectcolor='#ecf0f1').pack(anchor='w', pady=(4, 0))
        
        # Generate Statistics Plots button
        tk.Button(output_opts, text="⚙️ Generate Statistics Plots",
                 command=self.open_visualization_settings,
                 bg='#3498db', fg='white', font=('Arial', 8, 'bold'),
                 padx=8, pady=4, cursor='hand2').pack(anchor='w', padx=(20, 0), pady=(3, 0))
        
        # Initialize plot selections (all enabled by default)
        self.network_plot_selections = {
            'pathway_states': tk.BooleanVar(value=True),
            'pathway_sig_v': tk.BooleanVar(value=True),
            'pathway_sig_h': tk.BooleanVar(value=True),
            'pathway_enrichment': tk.BooleanVar(value=True),
            'pathway_enrichment_selected': tk.BooleanVar(value=True),
            'chord_diagram': tk.BooleanVar(value=True),
        }
        
        # ═══════════════════════════════════════════════════════════════════
        # MODE OPTIONS
        # ═══════════════════════════════════════════════════════════════════
        mode_opts_label = tk.Label(action_frame, text="Mode:", font=('Arial', 8, 'bold'), 
                                   bg='#ecf0f1', fg='#2c3e50')
        mode_opts_label.pack(anchor='w', pady=(8, 2))
        
        mode_opts = tk.Frame(action_frame, bg='#ecf0f1')
        mode_opts.pack(fill='x')
        
        # Focused Mode checkbox (checked by default) - ALWAYS create fresh
        self.focused_mode = tk.BooleanVar(value=True)
        focused_cb = tk.Checkbutton(mode_opts, text="🎯 Focused Mode", 
                      variable=self.focused_mode,
                      font=('Arial', 8, 'bold'), bg='#ecf0f1', fg='#2c3e50',
                      activebackground='#ecf0f1', selectcolor='#ecf0f1')
        focused_cb.pack(anchor='w')
        self._create_tooltip(focused_cb, 
                             "Focused: Show only pathway-linked metabolites\n"
                             "Unchecked: Show all metabolites including secondary connections")
        
        # Class Compression checkbox (checked by default)
        self.compress_by_class = tk.BooleanVar(value=True)
        compress_cb = tk.Checkbutton(mode_opts, text="📦 Compress Lipid Classes", 
                      variable=self.compress_by_class,
                      font=('Arial', 8), bg='#ecf0f1', fg='#2c3e50',
                      activebackground='#ecf0f1', selectcolor='#ecf0f1')
        compress_cb.pack(anchor='w', padx=(20, 0))
        self._create_tooltip(compress_cb, 
                             "Compress: Merge class members sharing ALL same pathways\n"
                             "Individual nodes shown for members with unique pathways\n"
                             "Requires Class_name column in data")

        # PATHWAY ENRICHMENT FILTERS (used during calculation)
        enrichment_filter_frame = tk.LabelFrame(left_content, text="🔬 Pathway Enrichment Filters", font=('Arial', 9, 'bold'),
                                               bg='#e8f5e9', fg='#2c3e50', padx=5, pady=5)
        enrichment_filter_frame.pack(fill='x', anchor='w', pady=(0, 10))

        # Component significance control (compound p-value threshold)
        comp_sig_frame = tk.Frame(enrichment_filter_frame, bg='#e8f5e9')
        comp_sig_frame.pack(fill='x', pady=(0, 6))
        tk.Label(comp_sig_frame, text="Compound p <", font=('Arial', 8), bg='#e8f5e9').pack(side='left')
        shared_metabolite_pvalue_var = self._get_shared_metabolite_pvalue_var()
        tk.Entry(comp_sig_frame, textvariable=shared_metabolite_pvalue_var, font=('Arial', 8)).pack(side='left', fill='x', expand=True, padx=(4, 0))

        # Load saved settings from pathway annotation if available
        saved_settings = self.data_manager.get_value('pathway_settings', {})
        min_met_default = str(saved_settings.get('min_metabolites', 1))
        pval_default = str(saved_settings.get('p_threshold', 0.05))

        # Min Metabolites and P-value
        row1 = tk.Frame(enrichment_filter_frame, bg='#e8f5e9')
        row1.pack(fill='x', pady=(0, 6))
        left_col1 = tk.Frame(row1, bg='#e8f5e9')
        left_col1.pack(side='left', fill='x', expand=True, padx=(0, 4))
        tk.Label(left_col1, text="Min Mets:", font=('Arial', 8), bg='#e8f5e9').pack(anchor='w')
        self.network_min_metabolites = tk.StringVar(value=min_met_default)
        tk.Entry(left_col1, textvariable=self.network_min_metabolites, font=('Arial', 8)).pack(fill='x')

        right_col1 = tk.Frame(row1, bg='#e8f5e9')
        right_col1.pack(side='left', fill='x', expand=True, padx=(4, 0))
        tk.Label(right_col1, text="P-val ≤:", font=('Arial', 8), bg='#e8f5e9').pack(anchor='w')
        self.network_p_threshold = tk.StringVar(value=pval_default)
        tk.Entry(right_col1, textvariable=self.network_p_threshold, font=('Arial', 8)).pack(fill='x')

        # NETWORK VISUALIZATION FILTERS (used for display selection)
        viz_filter_frame = tk.LabelFrame(left_content, text="🎨 Network Visualization Filters", font=('Arial', 9, 'bold'),
                                        bg='#fff3e0', fg='#2c3e50', padx=5, pady=5)
        viz_filter_frame.pack(fill='x', anchor='w', pady=(0, 10))

        # Z-score thresholds
        z_threshold_default = saved_settings.get('z_threshold', 1.0)
        row2 = tk.Frame(viz_filter_frame, bg='#fff3e0')
        row2.pack(fill='x', pady=(0, 6))
        left_col2 = tk.Frame(row2, bg='#fff3e0')
        left_col2.pack(side='left', fill='x', expand=True, padx=(0, 4))
        tk.Label(left_col2, text="Z+ Act:", font=('Arial', 8), bg='#fff3e0').pack(anchor='w')
        self.network_zscore_activation = tk.StringVar(value=str(z_threshold_default))
        tk.Entry(left_col2, textvariable=self.network_zscore_activation, font=('Arial', 8)).pack(fill='x')

        right_col2 = tk.Frame(row2, bg='#fff3e0')
        right_col2.pack(side='left', fill='x', expand=True, padx=(4, 0))
        tk.Label(right_col2, text="Z- Inh:", font=('Arial', 8), bg='#fff3e0').pack(anchor='w')
        self.network_zscore_inhibition = tk.StringVar(value=str(-z_threshold_default))
        tk.Entry(right_col2, textvariable=self.network_zscore_inhibition, font=('Arial', 8)).pack(fill='x')

        # Status filter
        tk.Label(viz_filter_frame, text="Pathway Status:", font=('Arial', 8), bg='#fff3e0').pack(anchor='w')
        self.network_status_filter = tk.StringVar(value="All")
        status_combo = ttk.Combobox(viz_filter_frame, textvariable=self.network_status_filter,
                                    values=["All", "Activated", "Inhibited", "Activated + Inhibited", "No Change"],
                                    state='readonly', font=('Arial', 8))
        status_combo.pack(fill='x', pady=(0, 6))

        # Pathway ML filters
        tk.Label(viz_filter_frame, text="Pathway ML Filters:", font=('Arial', 8, 'bold'), bg='#fff3e0').pack(anchor='w', pady=(0, 2))
        pathway_ml_row = tk.Frame(viz_filter_frame, bg='#fff3e0')
        pathway_ml_row.pack(fill='x', pady=(0, 6))
        pathway_ml_left = tk.Frame(pathway_ml_row, bg='#fff3e0')
        pathway_ml_left.pack(side='left', fill='x', expand=True, padx=(0, 4))
        tk.Label(pathway_ml_left, text="ML Score ≥:", font=('Arial', 8), bg='#fff3e0').pack(anchor='w')
        self.pathway_ml_score_min = tk.StringVar(value="0.0")
        tk.Entry(pathway_ml_left, textvariable=self.pathway_ml_score_min, font=('Arial', 8)).pack(fill='x')

        pathway_ml_right = tk.Frame(pathway_ml_row, bg='#fff3e0')
        pathway_ml_right.pack(side='left', fill='x', expand=True, padx=(4, 0))
        tk.Label(pathway_ml_right, text="Consensus Rank ≤:", font=('Arial', 8), bg='#fff3e0').pack(anchor='w')
        self.pathway_consensus_rank_max = tk.StringVar(value="999999")
        tk.Entry(pathway_ml_right, textvariable=self.pathway_consensus_rank_max, font=('Arial', 8)).pack(fill='x')

        # Filter buttons
        filter_buttons_frame = tk.Frame(viz_filter_frame, bg='#fff3e0')
        filter_buttons_frame.pack(fill='x')
        tk.Button(filter_buttons_frame, text="Apply", command=self.apply_pathway_filter,
                  bg='#27ae60', fg='white', font=('Arial', 8, 'bold'),
                  padx=6, pady=3).pack(side='left', fill='x', expand=True, padx=(0, 3))
        tk.Button(filter_buttons_frame, text="Reset", command=self.reset_pathway_filters,
                  bg='#95a5a6', fg='white', font=('Arial', 8, 'bold'),
                  padx=6, pady=3).pack(side='left', fill='x', expand=True, padx=(3, 0))
        
        # NETWORK COMPONENTS (Upstream and Diseases filters)
        components_frame = tk.LabelFrame(left_content, text="🧬 Network Components", font=('Arial', 11, 'bold'),
                                        bg='#ecf0f1', fg='#2c3e50', padx=5, pady=5)
        components_frame.pack(fill='x', pady=(0, 10))
        
        # Upstream Section
        upstream_section = tk.Frame(components_frame, bg='#ecf0f1')
        upstream_section.pack(fill='x', pady=(0, 8))
        
        self.network_include_upstream = tk.BooleanVar(value=True)
        tk.Checkbutton(upstream_section, text="Include Upstream Regulators", variable=self.network_include_upstream,
                      bg='#ecf0f1', font=('Arial', 9, 'bold'), anchor='w').pack(fill='x')
        
        # Apply Filter checkbox at top
        upstream_filter_checkbox = tk.Frame(upstream_section, bg='#ecf0f1')
        upstream_filter_checkbox.pack(fill='x', padx=(20, 0), pady=(5, 0))
        
        self.upstream_apply_filter = tk.BooleanVar(value=True)
        tk.Checkbutton(upstream_filter_checkbox, text="Apply Auto Filter", variable=self.upstream_apply_filter,
                      command=self.toggle_upstream_filter, bg='#ecf0f1', font=('Arial', 8, 'bold'),
                      anchor='w').pack(side='left')
        
        # Upstream filters
        self.upstream_filter_subframe = tk.Frame(upstream_section, bg='#ecf0f1')
        self.upstream_filter_subframe.pack(fill='x', padx=(20, 0), pady=(5, 0))
        
        self.upstream_min_label = tk.Label(self.upstream_filter_subframe, text="Min metabolites:", bg='#ecf0f1', font=('Arial', 8), width=12, anchor='w', state='normal')
        self.upstream_min_label.pack(side='left')
        self.upstream_min_metabolites = tk.StringVar(value="3")
        self.upstream_min_entry = tk.Entry(self.upstream_filter_subframe, textvariable=self.upstream_min_metabolites, width=6, font=('Arial', 8), state='normal')
        self.upstream_min_entry.pack(side='left', padx=2)
        
        # Upstream ML filters
        self.upstream_ml_subframe = tk.Frame(upstream_section, bg='#ecf0f1')
        self.upstream_ml_subframe.pack(fill='x', padx=(20, 0), pady=(5, 0))
        
        upstream_ml_left = tk.Frame(self.upstream_ml_subframe, bg='#ecf0f1')
        upstream_ml_left.pack(side='left', fill='x', expand=True, padx=(0, 5))
        self.upstream_ml_label = tk.Label(upstream_ml_left, text="ML Score ≥:", bg='#ecf0f1', font=('Arial', 8), state='normal')
        self.upstream_ml_label.pack(side='left', padx=(0, 2))
        self.upstream_ml_score_min = tk.StringVar(value="0.0")
        self.upstream_ml_entry = tk.Entry(upstream_ml_left, textvariable=self.upstream_ml_score_min, width=6, font=('Arial', 8), state='normal')
        self.upstream_ml_entry.pack(side='left')
        
        upstream_ml_right = tk.Frame(self.upstream_ml_subframe, bg='#ecf0f1')
        upstream_ml_right.pack(side='left', fill='x', expand=True, padx=(5, 0))
        self.upstream_consensus_label = tk.Label(upstream_ml_right, text="Consensus Rank ≤:", bg='#ecf0f1', font=('Arial', 8), state='normal')
        self.upstream_consensus_label.pack(side='left', padx=(0, 2))
        self.upstream_consensus_rank_max = tk.StringVar(value="25")
        self.upstream_consensus_entry = tk.Entry(upstream_ml_right, textvariable=self.upstream_consensus_rank_max, width=6, font=('Arial', 8), state='normal')
        self.upstream_consensus_entry.pack(side='left')
        
        # Upstream filter buttons
        self.upstream_filter_buttons = tk.Frame(upstream_section, bg='#ecf0f1')
        self.upstream_filter_buttons.pack(fill='x', padx=(20, 0), pady=(5, 0))
        
        self.upstream_apply_button = tk.Button(self.upstream_filter_buttons, text="Apply", command=self.apply_upstream_filter,
                 bg='#27ae60', fg='white', font=('Arial', 8, 'bold'),
                 padx=8, pady=3, cursor='hand2', state='normal')
        self.upstream_apply_button.pack(side='left', padx=(0, 5))
        
        self.upstream_reset_button = tk.Button(self.upstream_filter_buttons, text="Reset", command=self.reset_upstream_filter,
                 bg='#95a5a6', fg='white', font=('Arial', 8, 'bold'),
                 padx=8, pady=3, cursor='hand2', state='normal')
        self.upstream_reset_button.pack(side='left')
        
        # Diseases Section
        disease_section = tk.Frame(components_frame, bg='#ecf0f1')
        disease_section.pack(fill='x', pady=(8, 8))
        
        self.network_include_diseases = tk.BooleanVar(value=True)
        tk.Checkbutton(disease_section, text="Include Diseases", variable=self.network_include_diseases,
                      bg='#ecf0f1', font=('Arial', 9, 'bold'), anchor='w').pack(fill='x')
        
        # Apply Filter checkbox at top
        disease_filter_checkbox = tk.Frame(disease_section, bg='#ecf0f1')
        disease_filter_checkbox.pack(fill='x', padx=(20, 0), pady=(5, 0))
        
        self.disease_apply_filter = tk.BooleanVar(value=True)
        tk.Checkbutton(disease_filter_checkbox, text="Apply Auto Filter", variable=self.disease_apply_filter,
                      command=self.toggle_disease_filter, bg='#ecf0f1', font=('Arial', 8, 'bold'),
                      anchor='w').pack(side='left')
        
        # Disease filters
        self.disease_filter_subframe = tk.Frame(disease_section, bg='#ecf0f1')
        self.disease_filter_subframe.pack(fill='x', padx=(20, 0), pady=(5, 0))
        
        self.disease_min_label = tk.Label(self.disease_filter_subframe, text="Min metabolites:", bg='#ecf0f1', font=('Arial', 8), width=12, anchor='w', state='normal')
        self.disease_min_label.pack(side='left')
        self.disease_min_metabolites = tk.StringVar(value="1")
        self.disease_min_entry = tk.Entry(self.disease_filter_subframe, textvariable=self.disease_min_metabolites, width=6, font=('Arial', 8), state='normal')
        self.disease_min_entry.pack(side='left', padx=2)
        
        # Disease ML filters
        self.disease_ml_subframe = tk.Frame(disease_section, bg='#ecf0f1')
        self.disease_ml_subframe.pack(fill='x', padx=(20, 0), pady=(5, 0))
        
        disease_ml_left = tk.Frame(self.disease_ml_subframe, bg='#ecf0f1')
        disease_ml_left.pack(side='left', fill='x', expand=True, padx=(0, 5))
        self.disease_ml_label = tk.Label(disease_ml_left, text="ML Score ≥:", bg='#ecf0f1', font=('Arial', 8), state='normal')
        self.disease_ml_label.pack(side='left', padx=(0, 2))
        self.disease_ml_score_min = tk.StringVar(value="0.0")
        self.disease_ml_entry = tk.Entry(disease_ml_left, textvariable=self.disease_ml_score_min, width=6, font=('Arial', 8), state='normal')
        self.disease_ml_entry.pack(side='left')
        
        disease_ml_right = tk.Frame(self.disease_ml_subframe, bg='#ecf0f1')
        disease_ml_right.pack(side='left', fill='x', expand=True, padx=(5, 0))
        self.disease_consensus_label = tk.Label(disease_ml_right, text="Consensus Rank ≤:", bg='#ecf0f1', font=('Arial', 8), state='normal')
        self.disease_consensus_label.pack(side='left', padx=(0, 2))
        self.disease_consensus_rank_max = tk.StringVar(value="25")
        self.disease_consensus_entry = tk.Entry(disease_ml_right, textvariable=self.disease_consensus_rank_max, width=6, font=('Arial', 8), state='normal')
        self.disease_consensus_entry.pack(side='left')
        
        # Disease filter buttons
        self.disease_filter_buttons = tk.Frame(disease_section, bg='#ecf0f1')
        self.disease_filter_buttons.pack(fill='x', padx=(20, 0), pady=(5, 0))
        
        self.disease_apply_button = tk.Button(self.disease_filter_buttons, text="Apply", command=self.apply_disease_filter,
                 bg='#27ae60', fg='white', font=('Arial', 8, 'bold'),
                 padx=8, pady=3, cursor='hand2', state='normal')
        self.disease_apply_button.pack(side='left', padx=(0, 5))
        
        self.disease_reset_button = tk.Button(self.disease_filter_buttons, text="Reset", command=self.reset_disease_filter,
                 bg='#95a5a6', fg='white', font=('Arial', 8, 'bold'),
                 padx=8, pady=3, cursor='hand2', state='normal')
        self.disease_reset_button.pack(side='left')
        
        # ============ RIGHT COLUMN: TABLES + PROGRESS LOG ============
        # Wrap in canvas with scrollbar to scroll entire right column
        right_outer = tk.Frame(main_frame, bg='white', relief='groove', borderwidth=2)
        right_outer.grid(row=0, column=1, sticky='nsew', padx=(5, 5), pady=5)
        right_outer.grid_rowconfigure(0, weight=1)
        right_outer.grid_columnconfigure(0, weight=1)
        
        # Create canvas and scrollbar for entire right column
        right_canvas = tk.Canvas(right_outer, bg='white', highlightthickness=0)
        right_scrollbar = ttk.Scrollbar(right_outer, orient='vertical', command=right_canvas.yview)
        right_frame = tk.Frame(right_canvas, bg='white')
        
        # Bind canvas to update scroll region
        right_frame.bind('<Configure>', lambda e: right_canvas.configure(scrollregion=right_canvas.bbox('all')))
        
        # Create window in canvas
        canvas_window = right_canvas.create_window((0, 0), window=right_frame, anchor='nw')
        
        # Configure canvas
        right_canvas.configure(yscrollcommand=right_scrollbar.set)
        
        # Bind canvas resize to update window width
        def _on_canvas_configure(event):
            right_canvas.itemconfig(canvas_window, width=event.width)
        right_canvas.bind('<Configure>', _on_canvas_configure)
        
        # Grid canvas and scrollbar
        right_canvas.grid(row=0, column=0, sticky='nsew')
        right_scrollbar.grid(row=0, column=1, sticky='ns')
        
        # Mouse wheel scrolling
        def _on_mousewheel(event):
            right_canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        right_canvas.bind_all("<MouseWheel>", _on_mousewheel)
        
        # ===== TABLE NOTEBOOK (fixed tall height) =====
        table_container = tk.Frame(right_frame, bg='white', height=600)
        table_container.pack(fill='both', expand=False, padx=5, pady=5)
        table_container.pack_propagate(False)  # Maintain fixed height
        
        # Notebook for tabs
        table_notebook = ttk.Notebook(table_container)
        self.table_notebook = table_notebook
        table_notebook.pack(fill='both', expand=True)
        
        # ===== PROGRESS LOG (fixed height, scrollable within right column) =====
        progress_frame = tk.LabelFrame(right_frame, text="📊 Progress Log", font=('Arial', 9, 'bold'),
                                      bg='white', fg='#2c3e50', padx=5, pady=5, height=300)
        progress_frame.pack(fill='both', expand=False, padx=5, pady=(0, 5))
        progress_frame.pack_propagate(False)  # Maintain fixed height
        
        # Progress text
        self.network_progress_text = scrolledtext.ScrolledText(progress_frame, height=15, 
                                                              font=('Consolas', 9),
                                                              bg='#2c3e50', fg='#ecf0f1', wrap='word')
        self.network_progress_text.pack(fill='both', expand=True)
        self.network_progress_text.insert(tk.END, "📋 Ready to generate network.\n")
        self.network_progress_text.insert(tk.END, "⚠️ TIP: Launch Cytoscape before generating network for best experience.\n")
        self.network_progress_text.config(state='disabled')
        
        # ===== TAB 1: PATHWAYS =====
        pathway_tab = tk.Frame(table_notebook, bg='white')
        table_notebook.add(pathway_tab, text='🛤️ Pathways')
        
        # Pathway controls - use grid for better control
        pathway_controls = tk.Frame(pathway_tab, bg='white')
        pathway_controls.pack(fill='x', padx=5, pady=5)
        
        # Left side buttons
        left_buttons = tk.Frame(pathway_controls, bg='white')
        left_buttons.pack(side='left')
        
        tk.Button(left_buttons, text="✅ Select All", command=self.select_all_pathways,
                 bg='#27ae60', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='left', padx=2)
        tk.Button(left_buttons, text="❌ Deselect All", command=self.deselect_all_pathways,
                 bg='#e74c3c', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='left', padx=2)
        tk.Button(left_buttons, text="🗑️ Remove Selected", command=self.remove_selected_pathways,
                 bg='#f39c12', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='left', padx=2)
        
        # Auto-connect checkbox
        tk.Checkbutton(left_buttons, text="Auto-connect", 
                      variable=self.pathway_auto_connect,
                      bg='white', font=('Arial', 8), fg='#2c3e50').pack(side='left', padx=8)
        
        # Right side buttons
        right_buttons = tk.Frame(pathway_controls, bg='white')
        right_buttons.pack(side='right')
        
        # Export button on the right side
        tk.Button(right_buttons, text="📥 Export All Tables", command=self.export_pathway_table,
                 bg='#3498db', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='right', padx=2)
        
        # Rank-by selector for default ordering of the pathway table
        rank_frame = tk.Frame(right_buttons, bg='white')
        rank_frame.pack(side='right', padx=5)
        tk.Label(rank_frame, text="Rank by:", bg='white', font=('Arial', 8, 'bold')).pack(side='left', padx=(0, 3))
        saved_pathway_rank = self.data_manager.get_value('pathway_rank_mode', 'Default')
        self.pathway_rank_mode = tk.StringVar(value=saved_pathway_rank)
        rank_combo = ttk.Combobox(rank_frame, textvariable=self.pathway_rank_mode,
                                  values=["Default", "ML", "Consensus"], state='readonly', width=10)
        rank_combo.pack(side='left')
        def _on_pathway_rank_change(e):
            self.data_manager.set_value('pathway_rank_mode', self.pathway_rank_mode.get())
            self._populate_pathway_tree()
        rank_combo.bind("<<ComboboxSelected>>", _on_pathway_rank_change)
        # Tooltip
        self._create_tooltip(rank_combo, "Default: P-value → |Z|\nML: ML Score (diffusion) → |Z|\nConsensus: Avg rank of |Z| and ML")
        
        # Stats label (between left and right buttons)
        self.pathway_stats_label = tk.Label(pathway_controls, text="",
                                           bg='white', font=('Arial', 8), fg='#7f8c8d')
        self.pathway_stats_label.pack(side='left', padx=5)
        
        # Pathway table
        pathway_table_frame = tk.Frame(pathway_tab, bg='white')
        pathway_table_frame.pack(fill='both', expand=True, padx=5, pady=(0, 5))
        
        pathway_cols = (
            'Select', 'Pathway', 'Status', 'Z-Score', 'P-Value', 'FDR', '# Compounds',
            'Multi-Omics Rank', 'ML Score', 'PageRank', 'Betweenness', 'Consensus Rank'
        )
        self.pathway_tree = ttk.Treeview(pathway_table_frame, columns=pathway_cols, show='headings', height=35, selectmode='browse')
        
        self.pathway_tree.heading('Select', text='✓')
        self.pathway_tree.heading('Pathway', text='Pathway Name')
        self.pathway_tree.heading('Status', text='Status')
        self.pathway_tree.heading('Z-Score', text='Z-Score')
        self.pathway_tree.heading('P-Value', text='P-Value (raw)')
        self.pathway_tree.heading('FDR', text='FDR (adj)')
        self.pathway_tree.heading('# Compounds', text='# Compounds')
        self.pathway_tree.heading('Multi-Omics Rank', text='Multi-Omics Rank')
        self.pathway_tree.heading('ML Score', text='ML Score')
        self.pathway_tree.heading('PageRank', text='PageRank')
        self.pathway_tree.heading('Betweenness', text='Betweenness')
        self.pathway_tree.heading('Consensus Rank', text='Consensus Rank')
        
        # Tooltips for ML columns
        # self._create_tooltip(self.pathway_tree, 
        #     "ML Score: Diffusion signal from metabolites (higher=stronger)\n"
        #     "PageRank: Structural centrality (higher=hub)\n"
        #     "Betweenness: Bridge metric (higher=connector)\n"
        #     "Consensus Rank: Avg of |Z| rank & ML rank (lower=better)")
        
        self.pathway_tree.column('Select', width=50, anchor='center')
        self.pathway_tree.column('Pathway', width=340, anchor='w')
        self.pathway_tree.column('Status', width=100, anchor='center')
        self.pathway_tree.column('Z-Score', width=90, anchor='center')
        self.pathway_tree.column('P-Value', width=100, anchor='center')
        self.pathway_tree.column('FDR', width=100, anchor='center')
        self.pathway_tree.column('# Compounds', width=110, anchor='center')
        self.pathway_tree.column('Multi-Omics Rank', width=120, anchor='center')
        self.pathway_tree.column('ML Score', width=100, anchor='center')
        self.pathway_tree.column('PageRank', width=100, anchor='center')
        self.pathway_tree.column('Betweenness', width=110, anchor='center')
        self.pathway_tree.column('Consensus Rank', width=120, anchor='center')
        
        pathway_scrollbar_y = ttk.Scrollbar(pathway_table_frame, orient='vertical', command=self.pathway_tree.yview)
        pathway_scrollbar_x = ttk.Scrollbar(pathway_table_frame, orient='horizontal', command=self.pathway_tree.xview)
        self.pathway_tree.configure(yscrollcommand=pathway_scrollbar_y.set, xscrollcommand=pathway_scrollbar_x.set)
        
        self.pathway_tree.grid(row=0, column=0, sticky='nsew')
        pathway_scrollbar_y.grid(row=0, column=1, sticky='ns')
        pathway_scrollbar_x.grid(row=1, column=0, sticky='ew')
        
        pathway_table_frame.grid_rowconfigure(0, weight=1)
        pathway_table_frame.grid_columnconfigure(0, weight=1)
        
        # Enable double-click toggling of pathway selection (checkbox column)
        self.pathway_tree.bind('<Double-1>', self.toggle_pathway_selection)

        # Enable column sorting for pathway table
        self._enable_treeview_sorting(self.pathway_tree, [
            'check', 'str', 'str', 'float', 'float', 'float', 'int', 'int', 'float', 'float', 'float', 'int'
        ])
        
        # ===== TAB 2: COMPOUNDS =====
        metabolite_tab = tk.Frame(table_notebook, bg='white')
        self.metabolite_tab = metabolite_tab
        table_notebook.add(metabolite_tab, text='🧪 Compounds')
        
        metabolite_controls = tk.Frame(metabolite_tab, bg='white')
        metabolite_controls.pack(fill='x', padx=5, pady=5)
        
        self.metabolite_stats_label = tk.Label(metabolite_controls, text="",
                                              bg='white', font=('Arial', 8), fg='#7f8c8d')
        self.metabolite_stats_label.pack(side='left', padx=5)
        
        # Export button on the right side
        tk.Button(metabolite_controls, text="📥 Export All Tables", command=self.export_pathway_table,
                 bg='#3498db', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='right', padx=5)
        
        metabolite_table_frame = tk.Frame(metabolite_tab, bg='white')
        metabolite_table_frame.pack(fill='both', expand=True, padx=5, pady=(0, 5))
        
        metabolite_cols = ('Compounds', 'P-Value', 'log2FC', 'Pathways', 'Upstream', 'Diseases', 'Regulation')
        self.metabolite_tree = ttk.Treeview(metabolite_table_frame, columns=metabolite_cols, show='headings', height=35, selectmode='browse')
        
        for col in metabolite_cols:
            self.metabolite_tree.heading(col, text=col)
        
        self.metabolite_tree.column('Compounds', width=200, anchor='w')
        self.metabolite_tree.column('P-Value', width=100, anchor='center')
        self.metabolite_tree.column('log2FC', width=100, anchor='center')
        self.metabolite_tree.column('Pathways', width=80, anchor='center')
        self.metabolite_tree.column('Upstream', width=80, anchor='center')
        self.metabolite_tree.column('Diseases', width=80, anchor='center')
        self.metabolite_tree.column('Regulation', width=120, anchor='center')
        
        metabolite_scrollbar_y = ttk.Scrollbar(metabolite_table_frame, orient='vertical', command=self.metabolite_tree.yview)
        metabolite_scrollbar_x = ttk.Scrollbar(metabolite_table_frame, orient='horizontal', command=self.metabolite_tree.xview)
        self.metabolite_tree.configure(yscrollcommand=metabolite_scrollbar_y.set, xscrollcommand=metabolite_scrollbar_x.set)
        
        self.metabolite_tree.grid(row=0, column=0, sticky='nsew')
        metabolite_scrollbar_y.grid(row=0, column=1, sticky='ns')
        metabolite_scrollbar_x.grid(row=1, column=0, sticky='ew')
        
        metabolite_table_frame.grid_rowconfigure(0, weight=1)
        metabolite_table_frame.grid_columnconfigure(0, weight=1)
        # Enable column sorting for metabolite table
        self._enable_treeview_sorting(self.metabolite_tree, [
            'str', 'float', 'float', 'int', 'int', 'int', 'str'
        ])
        
        # ===== TAB 3: UPSTREAM =====
        upstream_tab = tk.Frame(table_notebook, bg='white')
        table_notebook.add(upstream_tab, text='⬆️ Upstream')
        
        upstream_controls = tk.Frame(upstream_tab, bg='white')
        upstream_controls.pack(fill='x', padx=5, pady=5)
        
        # Left side buttons
        upstream_left = tk.Frame(upstream_controls, bg='white')
        upstream_left.pack(side='left')
        
        tk.Button(upstream_left, text="✅ Select All", command=self.select_all_upstream,
                 bg='#27ae60', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='left', padx=2)
        tk.Button(upstream_left, text="❌ Deselect All", command=self.deselect_all_upstream,
                 bg='#e74c3c', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='left', padx=2)
        tk.Button(upstream_left, text="🗑️ Remove Selected", command=self.remove_selected_upstream,
                 bg='#f39c12', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='left', padx=2)
        
        # Auto-connect checkbox
        tk.Checkbutton(upstream_left, text="Auto-connect", 
                      variable=self.upstream_auto_connect,
                      bg='white', font=('Arial', 8), fg='#2c3e50').pack(side='left', padx=5)
        
        # Right side
        upstream_right = tk.Frame(upstream_controls, bg='white')
        upstream_right.pack(side='right')
        
        # Export button on the right side
        tk.Button(upstream_right, text="📥 Export All Tables", command=self.export_pathway_table,
                 bg='#3498db', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='right', padx=2)
        
        # Rank by selector for upstream regulators
        upstream_rank_frame = tk.Frame(upstream_right, bg='white')
        upstream_rank_frame.pack(side='right', padx=5)
        tk.Label(upstream_rank_frame, text="Rank by:", bg='white', font=('Arial', 8, 'bold')).pack(side='left', padx=(0, 3))
        saved_upstream_rank = self.data_manager.get_value('upstream_rank_mode', 'Default')
        self.upstream_rank_mode = tk.StringVar(value=saved_upstream_rank)
        upstream_rank_combo = ttk.Combobox(upstream_rank_frame, textvariable=self.upstream_rank_mode,
                                           values=["Default", "ML", "Consensus"], state='readonly', width=10)
        upstream_rank_combo.pack(side='left')
        def _on_upstream_rank_change(e):
            self.data_manager.set_value('upstream_rank_mode', self.upstream_rank_mode.get())
            self._populate_upstream_tree()
        upstream_rank_combo.bind("<<ComboboxSelected>>", _on_upstream_rank_change)
        
        self.upstream_stats_label = tk.Label(upstream_controls, text="",
                                            bg='white', font=('Arial', 8), fg='#7f8c8d')
        self.upstream_stats_label.pack(side='left', padx=5)

        upstream_table_frame = tk.Frame(upstream_tab, bg='white')
        upstream_table_frame.pack(fill='both', expand=True, padx=5, pady=(0, 5))
        
        upstream_cols = ('Select', 'Name', 'Type', '# Compounds', 'Multi-Omics Rank', 'ML Score', 'PageRank', 'Betweenness', 'Consensus Rank')
        self.upstream_tree = ttk.Treeview(upstream_table_frame, columns=upstream_cols, show='headings', height=35, selectmode='browse')
        
        for col in upstream_cols:
            self.upstream_tree.heading(col, text=col)
        
        self.upstream_tree.column('Select', width=50, anchor='center')
        self.upstream_tree.column('Name', width=400, anchor='w')
        self.upstream_tree.column('Type', width=150, anchor='center')
        self.upstream_tree.column('# Compounds', width=110, anchor='center')
        self.upstream_tree.column('Multi-Omics Rank', width=120, anchor='center')
        self.upstream_tree.column('ML Score', width=110, anchor='center')
        self.upstream_tree.column('PageRank', width=110, anchor='center')
        self.upstream_tree.column('Betweenness', width=120, anchor='center')
        self.upstream_tree.column('Consensus Rank', width=130, anchor='center')
        
        upstream_scrollbar_y = ttk.Scrollbar(upstream_table_frame, orient='vertical', command=self.upstream_tree.yview)
        upstream_scrollbar_x = ttk.Scrollbar(upstream_table_frame, orient='horizontal', command=self.upstream_tree.xview)
        self.upstream_tree.configure(yscrollcommand=upstream_scrollbar_y.set, xscrollcommand=upstream_scrollbar_x.set)
        
        self.upstream_tree.grid(row=0, column=0, sticky='nsew')
        upstream_scrollbar_y.grid(row=0, column=1, sticky='ns')
        upstream_scrollbar_x.grid(row=1, column=0, sticky='ew')
        
        upstream_table_frame.grid_rowconfigure(0, weight=1)
        upstream_table_frame.grid_columnconfigure(0, weight=1)
        self.upstream_tree.bind('<Double-1>', self.toggle_upstream_selection)
        # Enable column sorting for upstream table
        self._enable_treeview_sorting(self.upstream_tree, [
            'check', 'str', 'str', 'int', 'int', 'float', 'float', 'float', 'int'
        ])
        
        # ===== TAB 4: DISEASES =====
        disease_tab = tk.Frame(table_notebook, bg='white')
        table_notebook.add(disease_tab, text='🏥 Diseases')
        
        disease_controls = tk.Frame(disease_tab, bg='white')
        disease_controls.pack(fill='x', padx=5, pady=5)
        
        # Left side buttons
        disease_left = tk.Frame(disease_controls, bg='white')
        disease_left.pack(side='left')
        
        tk.Button(disease_left, text="✅ Select All", command=self.select_all_diseases,
                 bg='#27ae60', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='left', padx=2)
        tk.Button(disease_left, text="❌ Deselect All", command=self.deselect_all_diseases,
                 bg='#e74c3c', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='left', padx=2)
        tk.Button(disease_left, text="🗑️ Remove Selected", command=self.remove_selected_diseases,
                 bg='#f39c12', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='left', padx=2)
        
        # Auto-connect checkbox
        tk.Checkbutton(disease_left, text="Auto-connect", 
                      variable=self.disease_auto_connect,
                      bg='white', font=('Arial', 8), fg='#2c3e50').pack(side='left', padx=5)
        
        # Right side
        disease_right = tk.Frame(disease_controls, bg='white')
        disease_right.pack(side='right')
        
        # Export button on the right side
        tk.Button(disease_right, text="📥 Export All Tables", command=self.export_pathway_table,
                 bg='#3498db', fg='white', font=('Arial', 9, 'bold'),
                 padx=8, pady=3).pack(side='right', padx=2)

        # Rank by selector for diseases
        disease_rank_frame = tk.Frame(disease_right, bg='white')
        disease_rank_frame.pack(side='right', padx=5)
        tk.Label(disease_rank_frame, text="Rank by:", bg='white', font=('Arial', 8, 'bold')).pack(side='left', padx=(0, 3))
        saved_disease_rank = self.data_manager.get_value('disease_rank_mode', 'Default')
        self.disease_rank_mode = tk.StringVar(value=saved_disease_rank)
        disease_rank_combo = ttk.Combobox(disease_rank_frame, textvariable=self.disease_rank_mode,
                                          values=["Default", "ML", "Consensus"], state='readonly', width=10)
        disease_rank_combo.pack(side='left')
        def _on_disease_rank_change(e):
            self.data_manager.set_value('disease_rank_mode', self.disease_rank_mode.get())
            self._populate_disease_tree()
        disease_rank_combo.bind("<<ComboboxSelected>>", _on_disease_rank_change)
        
        self.diseases_stats_label = tk.Label(disease_controls, text="",
                                            bg='white', font=('Arial', 8), fg='#7f8c8d')
        self.diseases_stats_label.pack(side='left', padx=5)
        
        disease_table_frame = tk.Frame(disease_tab, bg='white')
        disease_table_frame.pack(fill='both', expand=True, padx=5, pady=(0, 5))
        
        disease_cols = ('Select', 'Disease', '# Compounds', 'Multi-Omics Rank', 'ML Score', 'PageRank', 'Betweenness', 'Consensus Rank')
        self.disease_tree = ttk.Treeview(disease_table_frame, columns=disease_cols, show='headings', height=35, selectmode='browse')
        
        for col in disease_cols:
            self.disease_tree.heading(col, text=col)
        
        self.disease_tree.column('Select', width=60, anchor='center')
        self.disease_tree.column('Disease', width=440, anchor='w')
        self.disease_tree.column('# Compounds', width=110, anchor='center')
        self.disease_tree.column('Multi-Omics Rank', width=120, anchor='center')
        self.disease_tree.column('ML Score', width=110, anchor='center')
        self.disease_tree.column('PageRank', width=110, anchor='center')
        self.disease_tree.column('Betweenness', width=120, anchor='center')
        self.disease_tree.column('Consensus Rank', width=130, anchor='center')
        
        disease_scrollbar_y = ttk.Scrollbar(disease_table_frame, orient='vertical', command=self.disease_tree.yview)
        disease_scrollbar_x = ttk.Scrollbar(disease_table_frame, orient='horizontal', command=self.disease_tree.xview)
        self.disease_tree.configure(yscrollcommand=disease_scrollbar_y.set, xscrollcommand=disease_scrollbar_x.set)
        
        self.disease_tree.grid(row=0, column=0, sticky='nsew')
        disease_scrollbar_y.grid(row=0, column=1, sticky='ns')
        disease_scrollbar_x.grid(row=1, column=0, sticky='ew')
        
        disease_table_frame.grid_rowconfigure(0, weight=1)
        disease_table_frame.grid_columnconfigure(0, weight=1)
        self.disease_tree.bind('<Double-1>', self.toggle_disease_selection)
        # Enable column sorting for disease table
        self._enable_treeview_sorting(self.disease_tree, [
            'check', 'str', 'int', 'int', 'float', 'float', 'float', 'int'
        ])

        # Ensure default labels reflect the current data mode
        self._update_component_labels()

    # ============ SORTING HELPERS ============

    def _enable_treeview_sorting(self, tree: ttk.Treeview, col_types):
        """Attach clickable sorting to each heading.

        col_types: list of type hints per column position: 'str' | 'float' | 'int' | 'check'
        """
        try:
            if not hasattr(self, '_tree_sort_state'):
                self._tree_sort_state = {}
            cols = list(tree['columns'])
            for idx, col in enumerate(cols):
                tree.heading(col, command=lambda cindex=idx: self._sort_tree_column(tree, cindex, col_types))
        except Exception as e:
            logger.warning(f"Failed to enable sorting: {e}")

    def _coerce_for_sort(self, val, typ):
        if typ == 'check':
            # map checkbox symbol to 1/0 to keep selected on top when ascending
            return 1 if str(val).strip() in ('☑️', '✓', '1', 'True') else 0
        if typ == 'int':
            try:
                return int(str(val).replace(',', '').strip())
            except Exception:
                try:
                    return int(float(str(val)))
                except Exception:
                    return 0
        if typ == 'float':
            try:
                return float(str(val).replace(',', '').strip())
            except Exception:
                return float('inf')
        # default string
        return str(val).lower()

    def _sort_tree_column(self, tree: ttk.Treeview, col_index: int, col_types):
        try:
            cols = list(tree['columns'])
            col = cols[col_index]
            typ = col_types[col_index] if col_index < len(col_types) else 'str'
            # Retrieve rows
            data = []
            for iid in tree.get_children(''):
                values = tree.item(iid, 'values')
                key_val = values[col_index] if col_index < len(values) else ''
                data.append((self._coerce_for_sort(key_val, typ), values, iid))

            # Toggle sort order
            state_key = id(tree)
            sort_state = self._tree_sort_state.get(state_key, {'col': None, 'reverse': False})
            reverse = sort_state['reverse'] if sort_state['col'] == col_index else False
            data.sort(key=lambda x: x[0], reverse=reverse)

            # Apply new order
            for index, (_, _values, iid) in enumerate(data):
                tree.move(iid, '', index)

            # Update state (toggle for next click)
            self._tree_sort_state[state_key] = {'col': col_index, 'reverse': not reverse}
        except Exception as e:
            logger.warning(f"Failed to sort column {col_index}: {e}")
    
    # ============ LOGGING HELPERS ============
    
    def _log_network(self, message):
        """Helper to log messages to network progress text"""
        try:
            if hasattr(self, 'network_progress_text'):
                self.network_progress_text.config(state='normal')
                self.network_progress_text.insert(tk.END, message)
                self.network_progress_text.see(tk.END)
                self.network_progress_text.update_idletasks()
                self.network_progress_text.config(state='disabled')
            print(f"[NETWORK LOG] {message.strip()}")
        except Exception as e:
            print(f"[ERROR] Failed to log to network text: {e}")
            print(f"[MESSAGE] {message}")
    
    # ============ FILE OPERATIONS ============
    
    def browse_network_output_folder(self):
        """Browse for network output folder"""
        # Use last browsed directory, or current folder value, or current working directory
        if self._last_browse_dir:
            initial_dir = self._last_browse_dir
        elif hasattr(self, 'network_output_folder') and self.network_output_folder.get():
            initial_dir = self.network_output_folder.get()
        else:
            initial_dir = os.getcwd()
        
        folder_path = filedialog.askdirectory(
            title="Select Output Folder for Network Results",
            initialdir=initial_dir
        )
        if folder_path:
            self.network_output_folder.set(folder_path)
            self._last_browse_dir = folder_path  # Remember this directory
            self._log_network(f"\n📁 Output folder set to: {folder_path}\n")
    
    def browse_network_file(self):
        """Browse for pathway-annotated Excel file to manually load into network tab"""
        # Use last browsed directory if available
        if self._last_browse_dir:
            initial_dir = self._last_browse_dir
        else:
            initial_dir = os.getcwd()
        
        file_path = filedialog.askopenfilename(
            title="Select Pathway-Annotated Excel File",
            filetypes=[
                ("Excel files", "*.xlsx"),
                ("All files", "*.*")
            ],
            initialdir=initial_dir
        )
        if file_path:
            self.network_input_file.set(file_path)
            self._last_browse_dir = os.path.dirname(file_path)  # Remember this directory
            self._log_network(f"\n📁 Manually selected file: {os.path.basename(file_path)}\n")
            self._log_network(f"   Path: {file_path}\n")
            
            # Load the data from the Excel file
            self._load_network_data_from_excel(file_path)
    
    def _load_network_data_from_excel(self, file_path):
        """Load pathway network data from Excel file and populate all tables"""
        try:
            import pandas as pd
            
            self._log_network(f"\n📊 Loading data from Excel file...\n")
            
            # Read all sheets
            excel_file = pd.ExcelFile(file_path)
            available_sheets = excel_file.sheet_names
            self._log_network(f"   Available sheets: {', '.join(str(s) for s in available_sheets)}\n")
            
            # Load Metabolites sheet
            if 'Metabolites' in available_sheets:
                df_metabolites = pd.read_excel(file_path, sheet_name='Metabolites')
                self._log_network(f"   ✓ Loaded {len(df_metabolites)} metabolites\n")
                
                # Store as filtered metabolites data
                self.pathway_filtered_metabolites_data = df_metabolites.copy()
                self.pathway_original_metabolites_data = df_metabolites.copy()
                self._cache_assigned_columns(self.pathway_filtered_metabolites_data)
            else:
                self._log_network(f"   ⚠️ No 'Metabolites' sheet found\n")
                df_metabolites = None
            
            # Load Selected_Pathways sheet and convert to pathways_data dict
            if 'Selected_Pathways' in available_sheets:
                df_pathways = pd.read_excel(file_path, sheet_name='Selected_Pathways')
                self._log_network(f"   ✓ Loaded {len(df_pathways)} pathways\n")
                
                # Convert DataFrame to pathways_data dictionary format
                pathways_data = {}
                for _, row in df_pathways.iterrows():
                    pathway_name = row.get('Pathway', '')
                    pathways_data[pathway_name] = {
                        'status': row.get('Status', 'No Change'),
                        'z_score': row.get('z_score', 0.0),
                        'combined_pvalue': row.get('combined_pvalue', 1.0),
                        'n_metabolites': int(row.get('n_metabolites', 0)),
                        'mean_log2fc': row.get('mean_log2fc', 0.0),
                    }
                
                # Store the pathways data
                self.pathway_original_pathways_data = pathways_data.copy()
                
                # Filter out "Other" and other generic pathways
                filtered_pathways = {}
                excluded_pathway_keywords = ['other', 'unknown', 'unclassified', 'miscellaneous']
                for pathway_name, stats in pathways_data.items():
                    # Check if pathway name contains any excluded keywords
                    pathway_lower = str(pathway_name).lower()
                    should_exclude = any(keyword in pathway_lower for keyword in excluded_pathway_keywords)
                    
                    if not should_exclude:
                        filtered_pathways[pathway_name] = stats
                    else:
                        self._log_network(f"   ⚠️ Filtered out pathway: '{pathway_name}'\n")
                
                self.pathway_filtered_pathways_data = filtered_pathways
                
                self._log_network(f"   ✓ Loaded {len(filtered_pathways)} pathways (filtered out {len(pathways_data) - len(filtered_pathways)} generic/invalid pathways)\n")
                
                # Update the pathway table
                self.update_pathway_table()
            else:
                self._log_network(f"   ⚠️ No 'Selected_Pathways' sheet found\n")
            
            # Load Upstream sheet
            if 'Upstream' in available_sheets:
                df_upstream = pd.read_excel(file_path, sheet_name='Upstream')
                if len(df_upstream) > 0:
                    self._log_network(f"   ✓ Loaded {len(df_upstream)} upstream regulators\n")
                    
                    # Convert DataFrame to upstream_data dict format
                    upstream_data = {}
                    for _, row in df_upstream.iterrows():
                        # Use cached feature column or fallback 'Name'
                        name = row.get(self.feature_col, '') if self.feature_col in row.index else row.get('Name', '')
                        reg_type = row.get('Type', '')
                        metabolites_str = row.get('Associated_Metabolites', '')
                        
                        # Parse metabolites from pipe-separated string
                        metabolites = []
                        if isinstance(metabolites_str, str) and metabolites_str:
                            metabolites = [m.strip() for m in metabolites_str.split('|') if m.strip()]
                        
                        upstream_data[name] = {
                            'type': reg_type,
                            'metabolites': metabolites,
                            'count': len(metabolites)
                        }
                    
                    self.original_upstream_data = upstream_data
                    
                    # Rebuild mappings to keep identifier alignment
                    self.upstream_to_metabolites = {}
                    self.metabolite_to_upstream = {}
                    for uname, info in upstream_data.items():
                        mets = set(info['metabolites'])
                        self.upstream_to_metabolites[uname] = mets
                        for m in mets:
                            if m not in self.metabolite_to_upstream:
                                self.metabolite_to_upstream[m] = set()
                            self.metabolite_to_upstream[m].add(uname)
                    
                    # Update upstream table
                    self._populate_upstream_table(upstream_data)
                else:
                    self._log_network(f"   ℹ️ Upstream sheet is empty\n")
            else:
                self._log_network(f"   ℹ️ No 'Upstream' sheet found\n")
            
            # Load Diseases sheet
            if 'Diseases' in available_sheets:
                df_diseases = pd.read_excel(file_path, sheet_name='Diseases')
                if len(df_diseases) > 0:
                    self._log_network(f"   ✓ Loaded {len(df_diseases)} diseases\n")
                    
                    # Convert DataFrame to disease_data dict format
                    disease_data = {}
                    for _, row in df_diseases.iterrows():
                        disease_name = row.get('Disease Name', '')
                        metabolites_str = row.get('Associated_Metabolites', '')
                        
                        # Parse metabolites from pipe-separated string
                        metabolites = []
                        if isinstance(metabolites_str, str) and metabolites_str:
                            metabolites = [m.strip() for m in metabolites_str.split('|') if m.strip()]
                        
                        disease_data[disease_name] = metabolites
                    
                    self.original_disease_data = disease_data
                    
                    # Update disease table
                    self._populate_disease_table(disease_data)
                else:
                    self._log_network(f"   ℹ️ Diseases sheet is empty\n")
            else:
                self._log_network(f"   ℹ️ No 'Diseases' sheet found\n")
            
            # Update metabolite table if we have metabolites
            if df_metabolites is not None and len(df_metabolites) > 0:
                self._populate_metabolite_table(df_metabolites)
            
            self._log_network(f"\n✅ Data loaded successfully from Excel file!\n")
            self.frame.after(0, lambda: messagebox.showinfo("Success", "Data loaded successfully from Excel file!"))
            
        except Exception as e:
            import traceback
            error_msg = str(e)
            self._log_network(f"\n❌ Error loading data from Excel: {error_msg}\n")
            self._log_network(f"{traceback.format_exc()}\n")
            self.frame.after(0, lambda msg=error_msg: messagebox.showerror("Error", f"Failed to load data from Excel:\n{msg}"))
    
    def _populate_upstream_table(self, upstream_data):
        """Populate upstream regulators table from data dict"""
        # Clear existing items
        for item in self.upstream_tree.get_children():
            self.upstream_tree.delete(item)
        
        if not upstream_data:
            return
        
        # Populate table
        for name, info in upstream_data.items():
            reg_type = info.get('type', '')
            count = info.get('count', 0)
            
            self.upstream_tree.insert('', 'end', values=(
                "☑️",  # Selected by default
                name,
                reg_type,
                count
            ))
        
        self.upstream_stats_label.config(text=f"Total: {len(upstream_data)}")
    
    def _populate_disease_table(self, disease_data):
        """Populate disease table from data dict"""
        # Clear existing items
        for item in self.disease_tree.get_children():
            self.disease_tree.delete(item)
        
        if not disease_data:
            return
        
        # Populate table
        for disease_name, metabolites in disease_data.items():
            count = len(metabolites) if isinstance(metabolites, list) else 0
            
            self.disease_tree.insert('', 'end', values=(
                "☑️",  # Selected by default
                disease_name,
                count
            ))
        
        self.diseases_stats_label.config(text=f"Total: {len(disease_data)}")
    
    def _populate_metabolite_table(self, df_metabolites):
        """Populate metabolite table from DataFrame"""
        # Clear existing items
        for item in self.metabolite_tree.get_children():
            self.metabolite_tree.delete(item)
        
        if df_metabolites is None or len(df_metabolites) == 0:
            return
        
        # Cache column assignments from imported DataFrame
        self._cache_assigned_columns(df_metabolites)
        
        # Populate table
        for _, row in df_metabolites.iterrows():
            compound_name = self._get_metabolite_value(row, 'name', '')
            log2fc = self._as_float(self._get_metabolite_value(row, 'log2fc', 0.0), 0.0)
            pvalue = self._as_float(self._get_metabolite_value(row, 'pvalue', 1.0), 1.0)
            
            # Determine regulation status
            if log2fc > 0:
                regulation = "Upregulated"
            elif log2fc < 0:
                regulation = "Downregulated"
            else:
                regulation = "No Change"
            
            # Get pathways
            pathways_str = ''
            if 'All_Pathways_Display' in row:
                pathways_str = str(row['All_Pathways_Display'])
            elif 'All_Pathways' in row:
                pathways_str = str(row['All_Pathways'])
            
            # Count pathways
            if pathways_str and pathways_str != 'nan':
                pathway_list = self._split_pathways(pathways_str)
                n_pathways = len(pathway_list)
            else:
                n_pathways = 0
            
            self.metabolite_tree.insert('', 'end', values=(
                "☑️",  # Selected by default
                compound_name,
                regulation,
                f"{log2fc:.3f}",
                f"{pvalue:.6f}",
                n_pathways
            ))
        
        # Update stats label (legacy format for Excel import)
        self.metabolite_stats_label.config(text=f"Shown: {len(df_metabolites)} | Hidden: 0 | Total: {len(df_metabolites)}")
        self.metabolite_stats_label.update_idletasks()
    
    def open_results_folder(self):
        """Open the results folder in file explorer"""
        try:
            if hasattr(self, 'network_output_folder') and self.network_output_folder.get():
                folder_path = self.network_output_folder.get()
                if os.path.exists(folder_path):
                    if platform.system() == 'Windows':
                        os.startfile(folder_path)
                    elif platform.system() == 'Darwin':
                        subprocess.Popen(['open', folder_path])
                    else:
                        subprocess.Popen(['xdg-open', folder_path])
                    self._log_network(f"📁 Opened results folder: {folder_path}\n")
                else:
                    messagebox.showinfo("Folder Not Found", f"Results folder doesn't exist yet.\n\nExpected location:\n{folder_path}")
            else:
                messagebox.showinfo("No Output Path", "Please set an output folder first.")
        except Exception as e:
            messagebox.showerror("Error", f"Could not open results folder:\n{str(e)}")
    
    # ============ PATHWAY SELECTION ============
    # ============ CASCADE SELECTION METHODS ============
    
    def _cascade_pathway_selection(self, pathway_name, selected_state):
        """
        When a pathway is selected/deselected, cascade to related metabolites/upstream/diseases.
        This creates a comprehensive network context automatically.
        Only cascades if the 'Auto-connect to other categories' checkbox is enabled.
        """
        # Check if auto-connect is enabled for pathways
        if not self.pathway_auto_connect.get():
            # Auto-connect is disabled, only select the pathway itself (no cascading)
            logger.info(f"ℹ️ Pathway auto-connect disabled: {pathway_name} selected independently")
            return
        
        # Get all metabolites in this pathway (from 'hits' - significant metabolites only)
        metabolites = self.pathway_to_metabolites.get(pathway_name, set())
        
        if selected_state:
            # SELECTING: cascade selection to related items
            # Track unique selections to avoid double counting
            diseases_selected = set()
            upstream_selected = set()
            
            for met in metabolites:
                # Select all diseases for this metabolite
                diseases = self.metabolite_to_diseases.get(met, set())
                for disease in diseases:
                    if disease in self.disease_selections:
                        self.disease_selections[disease] = True
                        diseases_selected.add(disease)
                
                # Select all upstream for this metabolite
                upstreams = self.metabolite_to_upstream.get(met, set())
                for upstream in upstreams:
                    # The upstream name from metabolite data is plain (e.g., "ALDOA")
                    # But upstream_selections uses composite keys (e.g., "Enzyme_ALDOA")
                    # We need to find and select the matching composite key
                    found = False
                    for sel_key in self.upstream_selections.keys():
                        # Extract the name part from composite key
                        # Composite keys are like "Enzyme_ALDOA" or "Transporter_SLC25A1"
                        if '_' in sel_key:
                            key_parts = sel_key.split('_', 1)
                            if len(key_parts) == 2:
                                key_type, key_name = key_parts
                                # Match if the name part matches the upstream name
                                if key_name == upstream and key_type in ['Enzyme', 'Transporter']:
                                    self.upstream_selections[sel_key] = True
                                    upstream_selected.add(sel_key)
                                    found = True
                        # Also try exact match
                        elif sel_key == upstream:
                            self.upstream_selections[sel_key] = True
                            upstream_selected.add(sel_key)
                            found = True
                    
                    if not found:
                        logger.debug(f"Could not find upstream selection key for: {upstream}")
            
            logger.info(f"✨ Cascaded selection for pathway: {pathway_name}")
            logger.info(f"   → Selected {len(metabolites)} metabolites (significant hits in pathway)")
            logger.info(f"   → Selected {len(diseases_selected)} unique diseases")
            logger.info(f"   → Selected {len(upstream_selected)} unique upstream regulators")
        else:
            # DESELECTING: cascade deselection to items ONLY used by this pathway
            # Collect metabolites from all other selected pathways, diseases, and upstream
            all_other_connected_metabolites = set()
            
            # Compounds from other selected pathways
            pathway_protected_mets = set()
            for other_pathway, other_mets in self.pathway_to_metabolites.items():
                if other_pathway != pathway_name and self.pathway_selections.get(other_pathway, False):
                    pathway_protected_mets.update(other_mets)
                    all_other_connected_metabolites.update(other_mets)
            
            # Compounds from selected diseases
            disease_protected_mets = set()
            for disease_name_key, disease_mets in self.disease_to_metabolites.items():
                if self.disease_selections.get(disease_name_key, False):
                    disease_protected_mets.update(disease_mets)
                    all_other_connected_metabolites.update(disease_mets)
                    logger.debug(f"   Disease '{disease_name_key}' protects: {disease_mets}")
            
            # Compounds from selected upstream regulators
            upstream_protected_mets = set()
            for upstream_key, upstream_mets in self.upstream_to_metabolites.items():
                # Check if this upstream is selected (may have composite keys)
                is_selected = False
                for sel_key, sel_value in self.upstream_selections.items():
                    if sel_value and (sel_key == upstream_key or sel_key.endswith(f"_{upstream_key}")):
                        is_selected = True
                        break
                if is_selected:
                    upstream_protected_mets.update(upstream_mets)
                    all_other_connected_metabolites.update(upstream_mets)
            
            logger.info(f"🔍 Deselection analysis for pathway: {pathway_name}")
            logger.info(f"   Pathway has {len(metabolites)} metabolites")
            logger.info(f"   Protected by other pathways: {len(pathway_protected_mets)} metabolites")
            logger.info(f"   Protected by diseases: {len(disease_protected_mets)} metabolites")
            logger.info(f"   Protected by upstream: {len(upstream_protected_mets)} metabolites")
            logger.info(f"   Total protected: {len(all_other_connected_metabolites)} metabolites")
            
            deselected_count = 0
            kept_count = 0
            
            # Only deselect items if metabolite has NO other connections
            for met in metabolites:
                if met not in all_other_connected_metabolites:
                    logger.info(f"   ❌ Removing '{met}' (no other connections)")
                    deselected_count += 1
                    
                    # DO NOT deselect diseases or upstream when deselecting pathways!
                    # Those are user-selected items and should only be deselected by the user
                    # or when deselecting that specific disease/upstream
                    
                    # The metabolite is no longer connected, but we don't touch
                    # user's explicit disease/upstream selections
                else:
                    reasons = []
                    if met in pathway_protected_mets:
                        reasons.append("other pathways")
                    if met in disease_protected_mets:
                        reasons.append("diseases")
                    if met in upstream_protected_mets:
                        reasons.append("upstream")
                    logger.info(f"   ✅ Keeping '{met}' (connected to: {', '.join(reasons)})")
                    kept_count += 1
            
            logger.info(f"   Summary: Kept {kept_count}, Removed {deselected_count}")
        
        # Refresh all trees to show selections
        self._refresh_tree_selections()
        
        if selected_state:
            logger.info(f"✨ Cascaded selection for pathway: {pathway_name}")
            logger.info(f"   → Selected {len(metabolites)} metabolites' connections")
        else:
            logger.info(f"✨ Cascaded deselection for pathway: {pathway_name}")
            logger.info(f"   → Deselected {deselected_count} metabolites' connections (Kept {kept_count})")
    
    def _cascade_disease_selection(self, disease_name, selected_state):
        """
        When a disease is selected/deselected, cascade to related metabolites/pathways/upstream.
        Only cascades if the 'Auto-connect to other categories' checkbox is enabled.
        """
        # Check if auto-connect is enabled for diseases
        if not self.disease_auto_connect.get():
            # Auto-connect is disabled, only select the disease itself (no cascading)
            logger.info(f"ℹ️ Disease auto-connect disabled: {disease_name} selected independently")
            return
        
        # Get all metabolites with this disease
        metabolites = self.disease_to_metabolites.get(disease_name, set())
        
        if selected_state:
            # SELECTING: cascade selection to related items
            pathways_selected = 0
            upstream_selected = 0
            
            # Track unique pathways to avoid double-counting
            unique_pathways = set()
            
            for met in metabolites:
                # Select all pathways containing this metabolite
                pathways = self.metabolite_to_pathways.get(met, set())
                for pathway in pathways:
                    if pathway in self.pathway_selections:
                        self.pathway_selections[pathway] = True
                        unique_pathways.add(pathway)
                
                # Select all upstream for this metabolite
                upstreams = self.metabolite_to_upstream.get(met, set())
                for upstream in upstreams:
                    # Match plain upstream name to composite selection keys
                    for sel_key in self.upstream_selections.keys():
                        if '_' in sel_key:
                            key_parts = sel_key.split('_', 1)
                            if len(key_parts) == 2:
                                key_type, key_name = key_parts
                                if key_name == upstream and key_type in ['Enzyme', 'Transporter']:
                                    self.upstream_selections[sel_key] = True
                                    upstream_selected += 1
                        elif sel_key == upstream:
                            self.upstream_selections[sel_key] = True
                            upstream_selected += 1
            
            pathways_selected = len(unique_pathways)
            
            logger.info(f"✨ Cascaded selection for disease: {disease_name}")
            logger.info(f"   → Selected {len(metabolites)} metabolites")
            logger.info(f"   → Selected {pathways_selected} unique pathways")
            logger.info(f"   → Selected {upstream_selected} upstream regulators")
        else:
            # DESELECTING: Just deselect the disease itself
            # DO NOT cascade deselection to pathways/upstream
            # Compounds will still appear in network if they have other connections
            # (to other selected pathways, diseases, or upstream)
            
            # Count how many metabolites are "orphaned" (lose their only connection)
            all_other_connected_metabolites = set()
            
            # Compounds from other selected diseases
            for other_disease, other_mets in self.disease_to_metabolites.items():
                if other_disease != disease_name and self.disease_selections.get(other_disease, False):
                    all_other_connected_metabolites.update(other_mets)
            
            # Metabolites from selected pathways
            for pathway_name_key, pathway_mets in self.pathway_to_metabolites.items():
                if self.pathway_selections.get(pathway_name_key, False):
                    all_other_connected_metabolites.update(pathway_mets)
            
            # Compounds from selected upstream regulators
            for upstream_key, upstream_mets in self.upstream_to_metabolites.items():
                is_selected = False
                for sel_key, sel_value in self.upstream_selections.items():
                    if sel_value and (sel_key == upstream_key or sel_key.endswith(f"_{upstream_key}")):
                        is_selected = True
                        break
                if is_selected:
                    all_other_connected_metabolites.update(upstream_mets)
            
            # Count orphaned metabolites (those that will lose network connection)
            orphaned_count = 0
            protected_count = 0
            for met in metabolites:
                if met not in all_other_connected_metabolites:
                    orphaned_count += 1
                else:
                    protected_count += 1
            
            logger.info(f"✨ Cascaded deselection for disease: {disease_name}")
            logger.info(f"   → {protected_count} metabolites still connected (via pathways/upstream/other diseases)")
            logger.info(f"   → {orphaned_count} metabolites will be hidden (no other connections)")
        
        # Refresh all trees
        self._refresh_tree_selections()
    
    def _cascade_upstream_selection(self, upstream_name, selected_state):
        """
        When an upstream regulator is selected/deselected, cascade to related metabolites/pathways/diseases.
        Note: upstream_name might be a composite key like "Enzyme_ABC123" or just "ABC123"
        Only cascades if the 'Auto-connect to other categories' checkbox is enabled.
        """
        # Check if auto-connect is enabled for upstream regulators
        if not self.upstream_auto_connect.get():
            # Auto-connect is disabled, only select the upstream itself (no cascading)
            logger.info(f"ℹ️ Upstream auto-connect disabled: {upstream_name} selected independently")
            return
        
        # Get all metabolites using this enzyme/transporter
        # The upstream_to_metabolites dict uses plain names, so extract the name part
        plain_name = upstream_name
        if '_' in upstream_name:
            # Could be "Enzyme_ABC" or "Transporter_ABC", extract the actual name
            parts = upstream_name.split('_', 1)
            if len(parts) == 2 and parts[0] in ['Enzyme', 'Transporter']:
                plain_name = parts[1]
        
        metabolites = self.upstream_to_metabolites.get(plain_name, set())
        if not metabolites:
            # Try with the original key
            metabolites = self.upstream_to_metabolites.get(upstream_name, set())
        
        if selected_state:
            # SELECTING: cascade selection to related items
            for met in metabolites:
                # Select all pathways containing this metabolite
                pathways = self.metabolite_to_pathways.get(met, set())
                for pathway in pathways:
                    if pathway in self.pathway_selections:
                        self.pathway_selections[pathway] = True
                
                # Select all diseases for this metabolite
                diseases = self.metabolite_to_diseases.get(met, set())
                for disease in diseases:
                    if disease in self.disease_selections:
                        self.disease_selections[disease] = True
        else:
            # DESELECTING: Just deselect the upstream itself
            # DO NOT cascade deselection to pathways/diseases
            # Compounds will still appear in network if they have other connections
            
            # Count how many metabolites are "orphaned" (lose their only connection)
            all_other_connected_metabolites = set()
            
            # Compounds from other selected upstream regulators
            for other_upstream, other_mets in self.upstream_to_metabolites.items():
                is_selected = False
                for sel_key in self.upstream_selections:
                    if self.upstream_selections[sel_key] and (sel_key == other_upstream or sel_key.endswith(f"_{other_upstream}")):
                        if sel_key != upstream_name:  # Don't count the one we're deselecting
                            is_selected = True
                            break
                if is_selected:
                    all_other_connected_metabolites.update(other_mets)
            
            # Compounds from selected pathways
            for pathway_name_key, pathway_mets in self.pathway_to_metabolites.items():
                if self.pathway_selections.get(pathway_name_key, False):
                    all_other_connected_metabolites.update(pathway_mets)
            
            # Compounds from selected diseases
            for disease_name_key, disease_mets in self.disease_to_metabolites.items():
                if self.disease_selections.get(disease_name_key, False):
                    all_other_connected_metabolites.update(disease_mets)
            
            # Count orphaned vs protected metabolites
            orphaned_count = 0
            protected_count = 0
            for met in metabolites:
                if met not in all_other_connected_metabolites:
                    orphaned_count += 1
                else:
                    protected_count += 1
            
            logger.info(f"✨ Cascaded deselection for upstream: {upstream_name}")
            logger.info(f"   → {protected_count} metabolites still connected (via pathways/diseases/other upstream)")
            logger.info(f"   → {orphaned_count} metabolites will be hidden (no other connections)")
        
        # Refresh all trees
        self._refresh_tree_selections()
    
    def _refresh_tree_selections(self):
        """Refresh all tree views to reflect current selection states"""
        # Refresh pathway tree
        for item in self.pathway_tree.get_children():
            values = list(self.pathway_tree.item(item, 'values'))
            if len(values) >= 2:
                pathway_name = values[1]
                is_selected = self.pathway_selections.get(pathway_name, False)
                values[0] = "☑️" if is_selected else "☐"
                self.pathway_tree.item(item, values=values)
        
        # Refresh upstream tree - use composite key (Type_Name)
        for item in self.upstream_tree.get_children():
            values = list(self.upstream_tree.item(item, 'values'))
            if len(values) >= 3:
                # values: [check, name, type, count]
                upstream_name = values[1]
                upstream_type = values[2]
                # Build composite key like we do in toggle/populate
                composite_key = f"{upstream_type}_{upstream_name}"
                is_selected = self.upstream_selections.get(composite_key, False)
                values[0] = "☑️" if is_selected else "☐"
                self.upstream_tree.item(item, values=values)
        
        # Refresh disease tree
        for item in self.disease_tree.get_children():
            values = list(self.disease_tree.item(item, 'values'))
            if len(values) >= 2:
                disease_name = values[1]
                is_selected = self.disease_selections.get(disease_name, False)
                values[0] = "☑️" if is_selected else "☐"
                self.disease_tree.item(item, values=values)
    
    def _collect_selection_snapshot(self):
        """Return the current checkbox selections for pathways/upstream/diseases."""
        snapshot = {'pathways': set(), 'upstream': set(), 'diseases': set()}

        def _is_checked(value):
            try:
                return str(value).strip() in self.CHECKED_SYMBOLS
            except Exception:
                return False

        tree = getattr(self, 'pathway_tree', None)
        if tree is not None:
            try:
                for iid in tree.get_children(''):
                    vals = tree.item(iid, 'values')
                    if len(vals) >= 2 and _is_checked(vals[0]):
                        name = str(vals[1]).strip()
                        if name:
                            snapshot['pathways'].add(name)
            except Exception:
                pass

        tree = getattr(self, 'upstream_tree', None)
        if tree is not None:
            try:
                for iid in tree.get_children(''):
                    vals = tree.item(iid, 'values')
                    if len(vals) >= 3 and _is_checked(vals[0]):
                        name = str(vals[1]).strip()
                        if name:
                            snapshot['upstream'].add(name)
            except Exception:
                pass

        tree = getattr(self, 'disease_tree', None)
        if tree is not None:
            try:
                for iid in tree.get_children(''):
                    vals = tree.item(iid, 'values')
                    if len(vals) >= 2 and _is_checked(vals[0]):
                        name = str(vals[1]).strip()
                        if name:
                            snapshot['diseases'].add(name)
            except Exception:
                pass

        return snapshot

    def _reset_all_selections(self):
        """Reset all selections to unselected state (called after network generation)"""
        # Reset pathway selections
        for pathway_name in self.pathway_selections:
            self.pathway_selections[pathway_name] = False
        
        # Reset upstream selections
        for upstream_name in self.upstream_selections:
            self.upstream_selections[upstream_name] = False
        
        # Reset disease selections
        for disease_name in self.disease_selections:
            self.disease_selections[disease_name] = False
        
        # Refresh all trees to show unselected state
        self._refresh_tree_selections()
        
        logger.info("✨ All selections reset - ready for new network")
    
    # ============ PATHWAY SELECTION ============
    
    def toggle_pathway_selection(self, event):
        """Toggle pathway selection on double-click and cascade to related items"""
        try:
            if not self.pathway_tree.selection():
                return
            item = self.pathway_tree.selection()[0]
            values = self.pathway_tree.item(item, 'values')
            if len(values) >= 2:
                pathway_name = values[1]
                current_state = self.pathway_selections.get(pathway_name, True)
                new_state = not current_state
                self.pathway_selections[pathway_name] = new_state
                new_symbol = "☑️" if new_state else "☐"
                new_values = list(values)
                new_values[0] = new_symbol
                self.pathway_tree.item(item, values=new_values)
                
                # CASCADE to related metabolites/diseases/upstream
                self._cascade_pathway_selection(pathway_name, new_state)
        except Exception as e:
            print(f"Error toggling pathway selection: {e}")
    
    def select_all_pathways(self):
        """Select all pathways and cascade to related upstream/diseases if auto-connect enabled"""
        selected_count = 0
        for item in self.pathway_tree.get_children():
            values = list(self.pathway_tree.item(item, 'values'))
            if len(values) >= 2:
                pathway_name = values[1]
                self.pathway_selections[pathway_name] = True
                values[0] = "☑️"
                self.pathway_tree.item(item, values=values)
                selected_count += 1
                
                # CRITICAL FIX: Trigger cascading selection for each pathway
                # This ensures upstream and diseases are also selected
                self._cascade_pathway_selection(pathway_name, selected_state=True)
        
        # Refresh UI to show cascaded selections
        self._refresh_tree_selections()
        
        self._log_network(f"Selected all {selected_count} pathways\n")
        if self.pathway_auto_connect.get():
            self._log_network(f"   (Auto-connect enabled: cascaded to related upstream & diseases)\n")
    
    def deselect_all_pathways(self):
        """Deselect all pathways and cascade to related upstream/diseases if auto-connect enabled"""
        deselected_count = 0
        for item in self.pathway_tree.get_children():
            values = list(self.pathway_tree.item(item, 'values'))
            if len(values) >= 2:
                pathway_name = values[1]
                self.pathway_selections[pathway_name] = False
                values[0] = "☐"
                self.pathway_tree.item(item, values=values)
                deselected_count += 1
                
                # CRITICAL FIX: Trigger cascading deselection for each pathway
                # This ensures orphaned upstream and diseases are also deselected
                self._cascade_pathway_selection(pathway_name, selected_state=False)
        
        # Refresh UI to show cascaded deselections
        self._refresh_tree_selections()
        
        self._log_network(f"Deselected all {deselected_count} pathways\n")
        if self.pathway_auto_connect.get():
            self._log_network(f"   (Auto-connect enabled: cascaded deselection to related items)\n")
    
    def remove_selected_pathways(self):
        """Remove selected pathways from the table"""
        if self.pathway_filtered_pathways_data is None:
            messagebox.showwarning("Warning", "No pathway data available")
            return
        
        selected_pathways = [name for name, selected in self.pathway_selections.items() if selected]
        
        if not selected_pathways:
            messagebox.showwarning("Warning", "No pathways selected for removal")
            return
        
        for pathway_name in selected_pathways:
            if pathway_name in self.pathway_filtered_pathways_data:
                del self.pathway_filtered_pathways_data[pathway_name]
        
        self.update_pathway_table()
        self._log_network(f"Removed {len(selected_pathways)} pathways\n")
    
    # ============ UPSTREAM SELECTION ============
    
    def toggle_upstream_selection(self, event):
        """Toggle upstream selection on double-click and cascade to related items"""
        try:
            if not self.upstream_tree.selection():
                return
            item = self.upstream_tree.selection()[0]
            values = list(self.upstream_tree.item(item, 'values'))
            if len(values) >= 3:
                name = values[1]
                reg_type = values[2]
                key = f"{reg_type}_{name}"
                new_state = not self.upstream_selections.get(key, True)
                self.upstream_selections[key] = new_state
                new_symbol = "☑️" if new_state else "☐"
                new_values = list(values)
                new_values[0] = new_symbol
                self.upstream_tree.item(item, values=new_values)
                
                # CASCADE to related metabolites/pathways/diseases
                self._cascade_upstream_selection(key, new_state)
        except Exception as e:
            print(f"Error toggling upstream selection: {e}")
    
    def select_all_upstream(self):
        """Select all upstream regulators"""
        selected_count = 0
        for item in self.upstream_tree.get_children():
            values = list(self.upstream_tree.item(item, 'values'))
            if len(values) >= 3:
                name = values[1]
                reg_type = values[2]
                key = f"{reg_type}_{name}"
                self.upstream_selections[key] = True
                values[0] = "☑️"
                self.upstream_tree.item(item, values=values)
                selected_count += 1
        self._log_network(f"Selected all {selected_count} upstream regulators\n")
    
    def deselect_all_upstream(self):
        """Deselect all upstream regulators"""
        deselected_count = 0
        for item in self.upstream_tree.get_children():
            values = list(self.upstream_tree.item(item, 'values'))
            if len(values) >= 3:
                name = values[1]
                reg_type = values[2]
                key = f"{reg_type}_{name}"
                self.upstream_selections[key] = False
                values[0] = "☐"
                self.upstream_tree.item(item, values=values)
                deselected_count += 1
        self._log_network(f"Deselected all {deselected_count} upstream regulators\n")
    
    def remove_selected_upstream(self):
        """Remove selected upstream regulators"""
        items_to_remove = []
        for item in self.upstream_tree.get_children():
            values = list(self.upstream_tree.item(item, 'values'))
            if len(values) >= 3:
                name = values[1]
                reg_type = values[2]
                key = f"{reg_type}_{name}"
                if self.upstream_selections.get(key, False):
                    items_to_remove.append(item)
        
        for item in items_to_remove:
            self.upstream_tree.delete(item)
        
        remaining = len(self.upstream_tree.get_children())
        self._log_network(f"Removed {len(items_to_remove)} upstream regulators, {remaining} remaining\n")
    
    # ============ DISEASE SELECTION ============
    
    def toggle_disease_selection(self, event):
        """Toggle disease selection on double-click and cascade to related items"""
        try:
            if not self.disease_tree.selection():
                return
            item = self.disease_tree.selection()[0]
            values = list(self.disease_tree.item(item, 'values'))
            if len(values) >= 2:
                disease_name = values[1]
                new_state = not self.disease_selections.get(disease_name, True)
                self.disease_selections[disease_name] = new_state
                new_symbol = "☑️" if new_state else "☐"
                new_values = list(values)
                new_values[0] = new_symbol
                self.disease_tree.item(item, values=new_values)
                
                # CASCADE to related metabolites/pathways/upstream
                self._cascade_disease_selection(disease_name, new_state)
        except Exception as e:
            print(f"Error toggling disease selection: {e}")
    
    def select_all_diseases(self):
        """Select all diseases"""
        selected_count = 0
        for item in self.disease_tree.get_children():
            values = list(self.disease_tree.item(item, 'values'))
            if len(values) >= 2:
                disease_name = values[1]
                self.disease_selections[disease_name] = True
                values[0] = "☑️"
                self.disease_tree.item(item, values=values)
                selected_count += 1
        self._log_network(f"Selected all {selected_count} diseases\n")
    
    def deselect_all_diseases(self):
        """Deselect all diseases"""
        deselected_count = 0
        for item in self.disease_tree.get_children():
            values = list(self.disease_tree.item(item, 'values'))
            if len(values) >= 2:
                disease_name = values[1]
                self.disease_selections[disease_name] = False
                values[0] = "☐"
                self.disease_tree.item(item, values=values)
                deselected_count += 1
        self._log_network(f"Deselected all {deselected_count} diseases\n")
    
    def remove_selected_diseases(self):
        """Remove selected diseases"""
        items_to_remove = []
        for item in self.disease_tree.get_children():
            values = list(self.disease_tree.item(item, 'values'))
            if len(values) >= 2:
                disease_name = values[1]
                if self.disease_selections.get(disease_name, False):
                    items_to_remove.append(item)
        
        for item in items_to_remove:
            self.disease_tree.delete(item)
        
        remaining = len(self.disease_tree.get_children())
        self._log_network(f"Removed {len(items_to_remove)} diseases, {remaining} remaining\n")
    
    # ============ FILTER OPERATIONS ============
    
    def _filter_metabolites_by_pvalue(self, metabolites_df, pvalue_threshold):
        """Filter metabolites by p-value threshold.
        
        This is a critical filter applied FIRST to ensure only significant metabolites
        are included in pathway calculations and network analysis.
        
        Parameters:
        -----------
        metabolites_df : pd.DataFrame
            DataFrame of metabolites with pvalue column
        pvalue_threshold : float
            P-value threshold for metabolite significance
            
        Returns:
        --------
        pd.DataFrame
            Filtered DataFrame containing only metabolites with pvalue < threshold
        """
        if metabolites_df is None or len(metabolites_df) == 0:
            return metabolites_df
        
        p_col = self.pvalue_col if self.pvalue_col in metabolites_df.columns else self._resolve_assigned_column(
            'P-Value', metabolites_df, ['pvalue', 'P-Value', 'p_value', 'p-value'])

        if not p_col or p_col not in metabolites_df.columns:
            logger.warning("⚠️  Verified p-value column not found in metabolites data")
            return metabolites_df
        
        try:
            component_label = self._component_label(plural=True)
            component_label_lower = component_label.lower()
            original_count = len(metabolites_df)
            numeric_p = pd.to_numeric(metabolites_df[p_col], errors='coerce')
            mask = numeric_p < pvalue_threshold
            filtered_df = metabolites_df.loc[mask].copy()
            filtered_count = len(filtered_df)
            
            if filtered_count < original_count:
                removed = original_count - filtered_count
                logger.info(f"📊 {component_label} p-value filter: {original_count} → {filtered_count} {component_label_lower} (removed {removed} with p ≥ {pvalue_threshold})")
                self._log_network(f"   Removed {removed} {component_label_lower} with p-value ≥ {pvalue_threshold}\n")
            else:
                logger.info(f"📊 {component_label} p-value filter: No change ({original_count} {component_label_lower} already pass p < {pvalue_threshold})")
                self._log_network(f"   All {original_count} {component_label_lower} already pass p < {pvalue_threshold} threshold\n")
            
            return filtered_df
        except Exception as e:
            logger.error(f"Error filtering metabolites by p-value: {str(e)}")
            return metabolites_df
    
    def apply_pathway_filter(self):
        """Apply pathway filtering based on parameters
        
        IMPORTANT: Metabolite p-value threshold is applied FIRST to ensure 
        only significant metabolites are included in all calculations.
        """
        try:
            self._log_network("🔍 Applying pathway filters...\n")

            component_label = self._component_label(plural=True)
            component_label_lower = component_label.lower()
            
            # Get filter parameters
            min_metabolites = int(self.network_min_metabolites.get())
            p_threshold = float(self.network_p_threshold.get())
            z_activation = float(self.network_zscore_activation.get())
            z_inhibition = float(self.network_zscore_inhibition.get())
            status_filter = self.network_status_filter.get()
            
            # FIRST: Apply metabolite p-value threshold filter
            # This is critical - only significant metabolites should be used in pathway analysis
            try:
                metabolite_pvalue_threshold = float(self.component_pvalue_threshold.get())
                self._log_network(f"📊 {component_label} p-value threshold: p < {metabolite_pvalue_threshold}\n")
            except (ValueError, AttributeError):
                metabolite_pvalue_threshold = 0.05
                logger.warning(f"Could not parse {component_label_lower} p-value threshold, using default 0.05")
            
            # Filter metabolites by p-value first
            # NOTE: If data already comes filtered from annotation step, this may not change anything
            filtered_metabolites = self._filter_metabolites_by_pvalue(
                self.pathway_original_metabolites_data,
                metabolite_pvalue_threshold
            )
            
            # Update filtered metabolites for downstream use
            self.pathway_filtered_metabolites_data = filtered_metabolites
            
            # Update metabolite table with filtered data
            self._populate_metabolite_tree()
            
            # Check if we have original data
            if not self.pathway_original_pathways_data:
                messagebox.showwarning("No Data", "No pathway data available to filter")
                return
            
            # Import disease detection function
            try:
                from main_script.metabolite_pathway_network import is_disease_pathway_global
            except ImportError:
                import re
                _disease_suffix_pattern = re.compile(r'(emia|uria|osis)(?:\s|$)', re.IGNORECASE)
                def is_disease_pathway_global(pathway_name):
                    pathway_lower = pathway_name.lower()
                    pathway_exclusions = [
                        "warburg effect", "metabolism", "pathway", "trna charging", "urea cycle",
                        "metabolic", "biosynthesis", "degradation", "oxidation", "reduction"
                    ]
                    for exclusion in pathway_exclusions:
                        if exclusion in pathway_lower:
                            return False
                    if _disease_suffix_pattern.search(pathway_name):
                        return True
                    disease_keywords = ["disease", "disorder", "syndrome", "deficiency", "newborn", "type", "transient"]
                    for keyword in disease_keywords:
                        if keyword in pathway_lower:
                            return True
                    return False
            
            # Clear existing tree
            for item in self.pathway_tree.get_children():
                self.pathway_tree.delete(item)
            
            # Preserve selections
            prev_selections = getattr(self, 'pathway_selections', {}).copy()
            self.pathway_selections = {}
            
            # Create a new filtered pathways dict with recalculated status
            self.pathway_filtered_pathways_data = {}
            
            # Apply filters
            displayed = 0
            filtered_out = 0
            
            for pathway_name, stats in self.pathway_original_pathways_data.items():
                if pathway_name.startswith('_'):
                    continue
                
                # TEMPORARILY DISABLED: Skip disease pathways (debugging 25->24 drop)
                # if is_disease_pathway_global(pathway_name):
                #     filtered_out += 1
                #     self._log_network(f"   [DEBUG] Filtered out disease pathway: {pathway_name}\n")
                #     continue
                
                if not isinstance(stats, dict):
                    continue
                
                p_value = stats.get('combined_pvalue', 1.0)
                z_score = stats.get('z_score', 0.0)
                n_metabolites = stats.get('n_metabolites', 0)
                
                # CRITICAL FIX: Re-determine pathway status based on CURRENT z-score thresholds
                # Don't use the status from annotation - recalculate it here!
                if z_score >= z_activation:
                    status = 'Activated'
                elif z_score <= z_inhibition:
                    status = 'Inhibited'
                else:
                    status = 'No Change'
                
                # Apply filters
                if n_metabolites < min_metabolites:
                    filtered_out += 1
                    continue
                
                if p_value > p_threshold:
                    filtered_out += 1
                    continue
                
                # Apply status filter
                if status_filter == "Activated" and status != "Activated":
                    filtered_out += 1
                    continue
                elif status_filter == "Inhibited" and status != "Inhibited":
                    filtered_out += 1
                    continue
                elif status_filter == "Activated + Inhibited" and status == "No Change":
                    filtered_out += 1
                    continue
                elif status_filter == "No Change" and status != "No Change":
                    filtered_out += 1
                    continue
                
                # Store in filtered pathways with updated status
                updated_stats = stats.copy()
                updated_stats['status'] = status  # Use recalculated status
                self.pathway_filtered_pathways_data[pathway_name] = updated_stats
                
            
            self._log_network(f"✅ Pathway filter applied: {len(self.pathway_filtered_pathways_data)} pathways pass filters, {filtered_out} filtered out\n")
            self._log_network(f"   Filters: Metabolite p < {metabolite_pvalue_threshold}, min={min_metabolites}, p≤{p_threshold}, status={status_filter}\n")
            
            # Repopulate pathway tree with ML scores (instead of manual insertion above)
            self._populate_pathway_tree()
            
            # Update upstream and disease tables with filtered metabolites
            self._populate_upstream_tree()
            self._populate_disease_tree()
            
            # CRITICAL: Update metabolite table to recalculate counts and filter zeros
            self._populate_metabolite_tree()
            
            self._update_stats_labels()
            

        except ValueError as e:
            messagebox.showerror("Error", f"Invalid filter value: {str(e)}")
        except Exception as e:
            logger.error(f"Failed to apply pathway filter: {str(e)}")
            import traceback
            traceback.print_exc()
            messagebox.showerror("Error", f"Failed to apply pathway filter: {str(e)}")
    
    def reset_pathway_filters(self):
        """Reset pathway filters to original data"""
        self._log_network("🔄 Resetting pathway filters...\n")
        if self.pathway_original_pathways_data is not None:
            self.pathway_filtered_pathways_data = self.pathway_original_pathways_data.copy()
            self._populate_pathway_tree()
            self._log_network("✅ Pathway filters reset to original data\n")
    
    def toggle_upstream_filter(self):
        """Toggle upstream auto filter on/off and enable/disable controls"""
        if self.upstream_apply_filter.get():
            # Auto Filter is ON - DISABLE controls and use default top 25 filter
            self.upstream_min_label.config(state='disabled')
            self.upstream_min_entry.config(state='disabled')
            self.upstream_ml_label.config(state='disabled')
            self.upstream_ml_entry.config(state='disabled')
            self.upstream_consensus_label.config(state='disabled')
            self.upstream_consensus_entry.config(state='disabled')
            self.upstream_apply_button.config(state='disabled')
            self.upstream_reset_button.config(state='disabled')
            self._log_network("✅ Upstream auto filter enabled (showing top 25, GUI controls disabled)\n")
        else:
            # Auto Filter is OFF - ENABLE controls to use GUI settings with no top limit
            self.upstream_min_label.config(state='normal')
            self.upstream_min_entry.config(state='normal')
            self.upstream_ml_label.config(state='normal')
            self.upstream_ml_entry.config(state='normal')
            self.upstream_consensus_label.config(state='normal')
            self.upstream_consensus_entry.config(state='normal')
            self.upstream_apply_button.config(state='normal')
            self.upstream_reset_button.config(state='normal')
            self._log_network("✅ Upstream auto filter disabled (using GUI settings, no top limit)\n")
        self._populate_upstream_tree()
        self._populate_metabolite_tree()
    
    def apply_upstream_filter(self):
        """Apply filters to upstream regulators by repopulating tree"""
        try:
            self._populate_upstream_tree()
            # Update metabolite table to recalculate counts
            self._populate_metabolite_tree()
            self._log_network(f"✅ Upstream filters applied\n")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to apply upstream filter: {str(e)}")
    
    def reset_upstream_filter(self):
        """Reset upstream filters to defaults"""
        self.upstream_min_metabolites.set("3")
        self.upstream_ml_score_min.set("0.0")
        self.upstream_consensus_rank_max.set("999999")
        self._populate_upstream_tree()
        # Update metabolite table to recalculate counts
        self._populate_metabolite_tree()
        self._log_network("✅ Upstream filters reset (Min=3, ML=0, Consensus=999999)\n")
    
    def toggle_disease_filter(self):
        """Toggle disease auto filter on/off and enable/disable controls"""
        if self.disease_apply_filter.get():
            # Auto Filter is ON - DISABLE controls and use default top 25 filter
            self.disease_min_label.config(state='disabled')
            self.disease_min_entry.config(state='disabled')
            self.disease_ml_label.config(state='disabled')
            self.disease_ml_entry.config(state='disabled')
            self.disease_consensus_label.config(state='disabled')
            self.disease_consensus_entry.config(state='disabled')
            self.disease_apply_button.config(state='disabled')
            self.disease_reset_button.config(state='disabled')
            self._log_network("✅ Disease auto filter enabled (showing top 25, GUI controls disabled)\n")
        else:
            # Auto Filter is OFF - ENABLE controls to use GUI settings with no top limit
            self.disease_min_label.config(state='normal')
            self.disease_min_entry.config(state='normal')
            self.disease_ml_label.config(state='normal')
            self.disease_ml_entry.config(state='normal')
            self.disease_consensus_label.config(state='normal')
            self.disease_consensus_entry.config(state='normal')
            self.disease_apply_button.config(state='normal')
            self.disease_reset_button.config(state='normal')
            self._log_network("✅ Disease auto filter disabled (using GUI settings, no top limit)\n")
        self._populate_disease_tree()
        self._populate_metabolite_tree()
    
    def apply_disease_filter(self):
        """Apply filters to diseases by repopulating tree"""
        try:
            self._populate_disease_tree()
            # Update metabolite table to recalculate counts
            self._populate_metabolite_tree()
            self._log_network(f"✅ Disease filters applied\n")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to apply disease filter: {str(e)}")
    
    def reset_disease_filter(self):
        """Reset disease filters to defaults"""
        self.disease_min_metabolites.set("1")
        self.disease_ml_score_min.set("0.0")
        self.disease_consensus_rank_max.set("999999")
        self._populate_disease_tree()
        # Update metabolite table to recalculate counts
        self._populate_metabolite_tree()
        self._log_network("✅ Disease filters reset (Min=1, ML=0, Consensus=999999)\n")
    
    def open_visualization_settings(self):
        """Open a window to configure network visualization parameters"""
        # Create new window
        settings_window = tk.Toplevel(self.root)
        settings_window.title("Network Visualization Settings")
        settings_window.geometry("500x650")
        settings_window.resizable(True, True)
        settings_window.grab_set()
        
        # Create main frame with scrollbar
        canvas = tk.Canvas(settings_window, bg='white', highlightthickness=0)
        scrollbar = ttk.Scrollbar(settings_window, orient='vertical', command=canvas.yview)
        scrollable_frame = tk.Frame(canvas, bg='white')
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Title
        title_label = tk.Label(scrollable_frame, text="Network Visualization Settings", 
                              font=('Arial', 12, 'bold'), bg='white', fg='#2c3e50')
        title_label.pack(pady=(10, 20), padx=10, fill='x')
        
        # Create parameter controls
        params = [
            ("Metabolite Font Size", self.network_metabolite_font_size, 6, 24),
            ("Pathway Font Size", self.network_pathway_font_size, 6, 24),
            ("Enzyme/Transporter Font Size", self.network_enzyme_font_size, 6, 24),
            ("Metabolite Max Characters/Line", self.network_metabolite_max_chars, 5, 30),
            ("Metabolite Max Lines", self.network_metabolite_max_lines, 1, 5),
            ("Pathway Max Characters/Line", self.network_pathway_max_chars, 8, 50),
            ("Pathway Max Lines", self.network_pathway_max_lines, 1, 8),
        ]
        
        controls = {}
        trace_bindings = []

        def make_update_func(current_lbl, variable):
            def update_display(*args):
                try:
                    if current_lbl.winfo_exists():
                        current_lbl.config(text=f"(Current: {variable.get()})")
                except tk.TclError:
                    # Widget was destroyed; callback will be removed on close.
                    pass
            return update_display

        def _bind_trace(variable, callback):
            trace_id = variable.trace_add('write', callback)
            trace_bindings.append((variable, trace_id))
        
        for label_text, var, min_val, max_val in params:
            # Create parameter frame
            param_frame = tk.Frame(scrollable_frame, bg='#ecf0f1', relief='flat', padx=10, pady=8)
            param_frame.pack(fill='x', padx=10, pady=(0, 5))
            
            # Label
            label = tk.Label(param_frame, text=label_text, font=('Arial', 9), bg='#ecf0f1', fg='#2c3e50')
            label.pack(anchor='w', pady=(0, 3))
            
            # Value frame (spinbox + current value display)
            value_frame = tk.Frame(param_frame, bg='#ecf0f1')
            value_frame.pack(fill='x')
            
            spinbox = tk.Spinbox(value_frame, from_=min_val, to=max_val, textvariable=var,
                                font=('Arial', 10), width=8, bg='white', fg='#2c3e50')
            spinbox.pack(side='left', padx=(0, 10))
            
            # Current value display
            current_label = tk.Label(value_frame, text=f"(Current: {var.get()})", 
                                    font=('Arial', 8), bg='#ecf0f1', fg='#7f8c8d')
            current_label.pack(side='left')
            
            _bind_trace(var, make_update_func(current_label, var))
            controls[label_text] = (spinbox, var)

        bargraph_frame = tk.LabelFrame(
            scrollable_frame,
            text="Bar Graph Label Settings (Top 20 & Selected)",
            font=('Arial', 9, 'bold'),
            bg='white',
            fg='#2c3e50'
        )
        bargraph_frame.pack(fill='x', padx=10, pady=(10, 5))

        for label_text, var, min_val, max_val in [
            ("Title Font Size", self.bargraph_title_font_size, 10, 40),
            ("Axis Title Font Size", self.bargraph_axis_title_font_size, 8, 36),
            ("X-Axis Tick Label Font Size", self.bargraph_tick_font_size, 6, 30),
            ("Pathway Label Font Size", self.bargraph_pathway_label_font_size, 6, 30),
            ("Plot Width (px)", self.bargraph_width_px, 900, 4000),
            ("Plot Height (px)", self.bargraph_height_px, 500, 4000),
            ("Legend Title Font Size", self.bargraph_legend_title_font_size, 8, 36),
            ("Legend Tick Font Size", self.bargraph_legend_tick_font_size, 6, 30),
            ("Legend Thickness", self.bargraph_legend_thickness, 6, 40),
            ("Legend Length (%)", self.bargraph_legend_length_percent, 30, 100),
            ("Pathway Max Characters/Line", self.bargraph_pathway_max_chars, 8, 50),
            ("Pathway Max Lines", self.bargraph_pathway_max_lines, 1, 8),
        ]:
            param_frame = tk.Frame(bargraph_frame, bg='#ecf0f1', relief='flat', padx=10, pady=8)
            param_frame.pack(fill='x', padx=10, pady=(0, 5))
            tk.Label(param_frame, text=label_text, font=('Arial', 9), bg='#ecf0f1', fg='#2c3e50').pack(anchor='w', pady=(0, 3))
            value_frame = tk.Frame(param_frame, bg='#ecf0f1')
            value_frame.pack(fill='x')
            spinbox = tk.Spinbox(value_frame, from_=min_val, to=max_val, textvariable=var,
                                font=('Arial', 10), width=8, bg='white', fg='#2c3e50')
            spinbox.pack(side='left', padx=(0, 10))
            current_label = tk.Label(value_frame, text=f"(Current: {var.get()})", font=('Arial', 8), bg='#ecf0f1', fg='#7f8c8d')
            current_label.pack(side='left')
            _bind_trace(var, make_update_func(current_label, var))

        chord_frame = tk.LabelFrame(
            scrollable_frame,
            text="Chord Diagram Label Settings",
            font=('Arial', 9, 'bold'),
            bg='white',
            fg='#2c3e50'
        )
        chord_frame.pack(fill='x', padx=10, pady=(10, 5))

        for label_text, var, min_val, max_val in [
            ("Metabolite Label Font Size", self.chord_metabolite_font_size, 6, 30),
            ("Pathway Label Font Size", self.chord_pathway_font_size, 6, 30),
            ("Metabolite Max Characters/Line", self.chord_metabolite_max_chars, 5, 30),
            ("Metabolite Max Lines", self.chord_metabolite_max_lines, 1, 5),
            ("Pathway Max Characters/Line", self.chord_pathway_max_chars, 8, 50),
            ("Pathway Max Lines", self.chord_pathway_max_lines, 1, 8),
            ("Chord Figure Width (inches)", self.chord_figure_width_inches, 8, 30),
            ("Chord Figure Height (inches)", self.chord_figure_height_inches, 8, 30),
        ]:
            param_frame = tk.Frame(chord_frame, bg='#ecf0f1', relief='flat', padx=10, pady=8)
            param_frame.pack(fill='x', padx=10, pady=(0, 5))
            tk.Label(param_frame, text=label_text, font=('Arial', 9), bg='#ecf0f1', fg='#2c3e50').pack(anchor='w', pady=(0, 3))
            value_frame = tk.Frame(param_frame, bg='#ecf0f1')
            value_frame.pack(fill='x')
            spinbox = tk.Spinbox(value_frame, from_=min_val, to=max_val, textvariable=var,
                                font=('Arial', 10), width=8, bg='white', fg='#2c3e50')
            spinbox.pack(side='left', padx=(0, 10))
            current_label = tk.Label(value_frame, text=f"(Current: {var.get()})", font=('Arial', 8), bg='#ecf0f1', fg='#7f8c8d')
            current_label.pack(side='left')
            _bind_trace(var, make_update_func(current_label, var))

        # Chord label casing option
        chord_case_frame = tk.Frame(chord_frame, bg='#ecf0f1', relief='flat', padx=10, pady=8)
        chord_case_frame.pack(fill='x', padx=10, pady=(0, 5))

        chord_case_title = tk.Label(
            chord_case_frame,
            text="Chord Diagram Label Casing",
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#2c3e50'
        )
        chord_case_title.pack(anchor='w', pady=(0, 3))

        chord_case_cb = tk.Checkbutton(
            chord_case_frame,
            text="🔤 Keep metabolite labels as-is",
            variable=self.chord_keep_metabolite_case,
            font=('Arial', 9),
            bg='#ecf0f1',
            fg='#2c3e50',
            activebackground='#ecf0f1',
            selectcolor='#ecf0f1'
        )
        chord_case_cb.pack(anchor='w')

        chord_case_note = tk.Label(
            chord_case_frame,
            text="Unchecked = convert metabolite labels to sentence case",
            font=('Arial', 8),
            bg='#ecf0f1',
            fg='#7f8c8d'
        )
        chord_case_note.pack(anchor='w', pady=(2, 0))
        
        # Add separator
        separator = ttk.Separator(scrollable_frame, orient='horizontal')
        separator.pack(fill='x', padx=10, pady=15)
        
        # Default values info box
        info_frame = tk.Frame(scrollable_frame, bg='#d5f4e6', relief='flat', padx=10, pady=8)
        info_frame.pack(fill='x', padx=10, pady=(0, 10))
        
        info_text = tk.Label(info_frame, text="Default Values:\n"
                            "• Metabolite Font: 10  •  Pathway Font: 9  •  Enzyme Font: 8\n"
                            "• Network labels: Metabolite Chars 10 / Lines 2, Pathway Chars 20 / Lines 3\n"
                            "• Bar graphs: Title 22, Axis Title 20, X Tick 15, Pathway Label 20, Size 1200x700, Legend 14/12/10/70%, Pathway Chars 40 / Lines 3\n"
                            "• Chord diagram: Metabolite Font 12, Pathway Font 12, Metabolite Chars 40 / Lines 3, Pathway Chars 40 / Lines 3, Size 12x12 in",
                            font=('Arial', 8), bg='#d5f4e6', fg='#27ae60', justify='left')
        info_text.pack(anchor='w')
        
        # Buttons frame
        button_frame = tk.Frame(scrollable_frame, bg='white')
        button_frame.pack(fill='x', padx=10, pady=(10, 10))
        
        def reset_to_defaults():
            """Reset all parameters to defaults"""
            self.network_metabolite_font_size.set(10)
            self.network_pathway_font_size.set(9)
            self.network_enzyme_font_size.set(8)
            self.network_metabolite_max_chars.set(10)
            self.network_metabolite_max_lines.set(2)
            self.network_pathway_max_chars.set(20)
            self.network_pathway_max_lines.set(3)
            self.bargraph_title_font_size.set(22)
            self.bargraph_axis_title_font_size.set(20)
            self.bargraph_tick_font_size.set(15)
            self.bargraph_pathway_label_font_size.set(20)
            self.bargraph_width_px.set(1200)
            self.bargraph_height_px.set(700)
            self.bargraph_legend_title_font_size.set(14)
            self.bargraph_legend_tick_font_size.set(12)
            self.bargraph_legend_thickness.set(10)
            self.bargraph_legend_length_percent.set(70)
            self.bargraph_pathway_max_chars.set(40)
            self.bargraph_pathway_max_lines.set(3)
            self.chord_metabolite_font_size.set(12)
            self.chord_pathway_font_size.set(12)
            self.chord_metabolite_max_chars.set(40)
            self.chord_metabolite_max_lines.set(3)
            self.chord_pathway_max_chars.set(40)
            self.chord_pathway_max_lines.set(3)
            self.chord_figure_width_inches.set(12)
            self.chord_figure_height_inches.set(12)
            self.chord_keep_metabolite_case.set(False)
            messagebox.showinfo("Reset", "All parameters reset to defaults!")
        
        def _close_settings_window():
            for variable, trace_id in trace_bindings:
                try:
                    variable.trace_remove('write', trace_id)
                except Exception:
                    pass
            try:
                if settings_window.winfo_exists():
                    settings_window.destroy()
            except Exception:
                pass

        settings_window.protocol("WM_DELETE_WINDOW", _close_settings_window)

        # Reset button
        reset_btn = tk.Button(button_frame, text="↺ Reset to Defaults", command=reset_to_defaults,
                             bg='#95a5a6', fg='white', font=('Arial', 9, 'bold'),
                             padx=10, pady=5, cursor='hand2')
        reset_btn.pack(side='left', padx=(0, 5))
        
        # Close button
        close_btn = tk.Button(button_frame, text="✓ Close", command=_close_settings_window,
                             bg='#27ae60', fg='white', font=('Arial', 9, 'bold'),
                             padx=10, pady=5, cursor='hand2')
        close_btn.pack(side='left')
        
        # Pack canvas and scrollbar
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')
    
    def _create_statistics_explanation_file(self, output_path):
        """Create a simple text file explaining all statistics and methods."""
        try:
            explanation_text = """METABOLITE PATHWAY ANALYSIS - STATISTICS GUIDE
===============================================

This file explains all statistical methods, scores, and rankings used in your analysis.
Written in simple terms without jargon.


1. PATHWAY ENRICHMENT METHODS
==============================

ORA (Over-Representation Analysis)
-----------------------------------
What it does:
- Tests if your changed metabolites appear in a pathway more than expected by chance
- Uses Fisher's exact test (a standard statistical test)

How to read:
- P-Value: Lower is better. Below 0.05 means statistically significant
- # Compounds: How many of your compounds are in this pathway

When to use:
- You have a list of significantly changed metabolites
- You want to know which pathways are "enriched" with your metabolites


IWPA (Integrated Weighted Pathway Analysis)
--------------------------------------------
What it does:
- Considers ALL metabolites (not just significant ones)
- Weights metabolites by their fold change AND p-value
- Looks at how consistently metabolites change in same direction

How to read:
- Z-Score: How strongly the pathway is affected
  • Positive Z-Score (e.g., +2.5) = Pathway ACTIVATED (metabolites increased)
  • Negative Z-Score (e.g., -2.5) = Pathway INHIBITED (metabolites decreased)
  • Higher absolute value = stronger effect
- P-Value: Lower is better. Below 0.05 means significant
- Status: Activated, Inhibited, or No Change (based on your Z-score thresholds)

When to use:
- You want to capture subtle but consistent changes
- You care about direction of change (up vs down)
- More sensitive than ORA for detecting pathway activity


2. MACHINE LEARNING SCORES
==========================

ML Score (Random Walk Diffusion)
---------------------------------
What it does:
- Simulates how metabolite signals "flow" through the network
- Starts at your changed metabolites, spreads to connected pathways
- Pathways receiving more signal have higher scores

How to read:
- Range: 0.0 to 1.0 (probability score)
- Good score: Above 0.1
- Strong score: Above 0.3
- Higher = pathway receives more signal from your metabolites

Why it matters:
- Finds pathways that are well-connected to your metabolites
- Can detect important pathways missed by traditional statistics
- Publication-quality network-based scoring


PageRank
--------
What it does:
- Measures how "central" or "hub-like" a pathway is in the network
- Same algorithm Google uses to rank web pages

How to read:
- Range: 0.0 to 1.0
- Good score: Above 0.01
- Strong score: Above 0.05
- Higher = pathway is more structurally important in the network

Why it matters:
- Identifies key regulatory pathways
- Hub pathways often control many metabolic processes


Betweenness Centrality
----------------------
What it does:
- Counts how often a pathway sits "between" other metabolites
- Measures if pathway acts as a bridge or connector

How to read:
- Absolute count (depends on network size)
- Good score: Above 100
- Strong score: Above 1000
- Higher = pathway connects more parts of the network

Why it matters:
- Bridge pathways are critical for network communication
- Disrupting these can have widespread effects


Consensus Rank
--------------
What it does:
- Combines traditional statistics with ML scores
- Average of two ranks:
  1. Rank by Z-Score (or p-value for diseases/upstream)
  2. Rank by ML Score

How to read:
- Integer rank: 1, 2, 3, etc.
- Good: Rank ≤ 10
- Strong: Rank ≤ 5
- LOWER IS BETTER (Rank 1 = best)

Why it matters:
- Highlights pathways strong in BOTH traditional stats AND network ML
- Best for prioritizing pathways for follow-up experiments
- Reduces false positives


3. RANKING MODES
================

Default Ranking
---------------
Pathways: Sorted by P-Value first, then |Z-Score|
Upstream/Diseases: Sorted by # Compounds, then ML Score

ML Ranking
----------
Sorted by ML Score (highest first), then by default metric
Use when: You want network-based prioritization

Consensus Ranking
-----------------
Sorted by Consensus Rank (lowest first)
Use when: You want balanced prioritization between stats and ML


4. FILTERING GUIDELINES
=======================

Statistical Filters
-------------------
- P-Value: Start with ≤ 0.05 (standard significance)
- Min Metabolites: At least 3-5 (avoids spurious results)
- Z-Score Activation: Typically ≥ 2.0
- Z-Score Inhibition: Typically ≤ -2.0

ML Filters
----------
- ML Score: Start with ≥ 0.1 (removes noise)
  • For very strict: ≥ 0.2
  • For very permissive: ≥ 0.05
- Consensus Rank: Top 10-20 pathways (≤ 10 or ≤ 20)


5. IMPORTANT METABOLITES
========================

How determined:
- Weight = sign(log2FC) × -log10(p-value)
  OR
- Weight = log2FC × -log10(p-value)

What makes a metabolite "important":
- Large fold change (high |log2FC|)
- High confidence (low p-value)
- Both together = highest weight

Example:
- log2FC = 2.0, p-value = 0.001
- Weight = 2.0 × -log10(0.001) = 2.0 × 3 = 6.0
- This metabolite strongly influences ML scores


6. NETWORK GENERATION
=====================

What's included:
- Metabolites (nodes colored by fold change)
- Pathways (nodes colored by status: red=inhibited, blue=activated)
- Upstream Regulators (enzymes/transporters, yellow nodes)
- Diseases (associated conditions, purple nodes)
- Edges (connections between metabolites and pathways/diseases/upstream)

Node colors:
- Red gradient: Downregulated/Inhibited
- Blue gradient: Upregulated/Activated
- Yellow: Upstream regulators
- Purple: Diseases

Node sizes:
- Larger nodes = more connections OR stronger effect


7. EXPORT FILES EXPLAINED
=========================

Excel Sheets:
-------------
1. Pathways: All pathway statistics + ML scores + metabolite list
2. Metabolites: Individual metabolite data with pathway associations
3. Upstream_Regulators: Enzymes/transporters + ML scores + metabolite list
4. Diseases: Associated diseases + ML scores + metabolite list

Network Files:
--------------
1. network.html: Interactive network (open in web browser)
2. network.graphml: For Cytoscape (advanced network visualization)
3. network_legend.png: Color/shape legend for the network
4. pathway_*.png: Statistical plots (if generated)


8. TIPS FOR PUBLICATION
=======================

1. Always report:
   - Statistical method used (ORA vs IWPA)
   - P-value threshold (e.g., p < 0.05)
   - Multiple testing correction (FDR/Bonferroni if applied)
   - Z-score thresholds for activation/inhibition

2. For ML scores:
   - Mention "Random Walk with Restart for diffusion scoring"
   - Cite network centrality metrics (PageRank, Betweenness)
   - Note: "Consensus ranking integrates statistical and network-based evidence"

3. Filtering transparency:
   - Report minimum metabolites threshold
   - State ML Score cutoff if used
   - Document any manual curation

4. Visualization:
   - Include network legend in figures
   - Color code by biological meaning (activation/inhibition)
   - Size nodes by importance (degree or effect size)


9. TROUBLESHOOTING
==================

No pathways showing up?
- Lower P-value threshold (try 0.1)
- Reduce minimum metabolites (try 2)
- Check if ML filters are too strict

Too many pathways?
- Stricter P-value (try 0.01)
- Higher ML Score filter (try 0.2)
- Use Consensus Rank ≤ 10

ML Scores all zero?
- Need at least 5-10 metabolites with good fold changes
- Check metabolite p-values are significant
- Ensure pathways have connections to metabolites

Consensus Rank seems wrong?
- It's an average of two ranks, not a direct score
- Lower is always better (rank 1 = best)
- Compare to both Z-Score and ML Score columns


10. GLOSSARY
============

log2FC (Log2 Fold Change):
- How much a metabolite increased or decreased
- Positive = increased, Negative = decreased
- Example: log2FC = 1.0 means 2x increase (doubled)

P-Value:
- Probability result happened by chance
- Lower = more confident it's real
- 0.05 = 5% chance it's random (standard cutoff)

FDR (False Discovery Rate):
- Corrects p-values when testing many pathways at once
- Controls false positives

Bipartite Network:
- Network with two types of nodes (metabolites and pathways)
- Metabolites connect to pathways, not to other metabolites

Random Walk:
- Simulation that "walks" through network connections
- Spreads signal from starting points (metabolites) to destinations (pathways)


=== END OF GUIDE ===

Generated by Metabolite Pathway Analysis Tool
For questions or issues, refer to the main documentation.
"""
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(explanation_text)
            self._log_network(f"📄 Statistics guide saved: {output_path}\n")
            return True
        except Exception as e:
            self._log_network(f"⚠️ Could not create statistics guide: {e}\n")
            return False
    
    def export_pathway_table(self):
        """Export all currently visible tables (pathways, metabolites, upstream, diseases) to Excel"""
        from tkinter import filedialog, messagebox
        import pandas as pd
        from datetime import datetime
        
        try:
            # Check if we have any data
            has_data = False
            
            # Check pathway tree
            pathway_items = []
            if hasattr(self, 'pathway_tree') and self.pathway_tree is not None:
                pathway_items = self.pathway_tree.get_children()
                if pathway_items:
                    has_data = True
            
            if not has_data:
                messagebox.showwarning("No Data", "No data to export. Apply filters or generate network first.")
                return
            
            # Ask for save location
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            default_filename = f"pathway_network_all_tables_{timestamp}.xlsx"
            
            filepath = filedialog.asksaveasfilename(
                defaultextension=".xlsx",
                filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")],
                initialfile=default_filename,
                title="Export All Network Tables"
            )
            
            if not filepath:
                return
            
            # Collect all data
            sheets_data = {}
            
            # 1. Export Pathways with metabolite list
            if pathway_items:
                pathway_data = []
                for item in pathway_items:
                    values = list(self.pathway_tree.item(item)['values'])
                    # Get pathway name (column index 1)
                    pathway_name = values[1] if len(values) > 1 else None
                    # Get metabolite list from pathway stats
                    metabolite_list = ''
                    if pathway_name and self.pathway_filtered_pathways_data:
                        stats = self.pathway_filtered_pathways_data.get(pathway_name, {})
                        if isinstance(stats, dict):
                            mets = stats.get('metabolites', [])
                            if isinstance(mets, (list, tuple, set)):
                                metabolite_list = ', '.join(sorted(set(str(m) for m in mets if m)))
                    values.append(metabolite_list)
                    pathway_data.append(values)
                if pathway_data:
                    pathway_columns = list(self.pathway_tree['columns']) + ['Metabolites']
                    sheets_data['Pathways'] = pd.DataFrame(pathway_data, columns=pathway_columns)
            
            # 2. Export Metabolites (no change needed)
            if hasattr(self, 'metabolite_tree') and self.metabolite_tree is not None:
                metabolite_items = self.metabolite_tree.get_children()
                if metabolite_items:
                    metabolite_data = []
                    for item in metabolite_items:
                        values = self.metabolite_tree.item(item)['values']
                        if len(values) == 7:
                            metabolite_data.append(values)
                    if metabolite_data:
                        metabolite_columns = ['Compounds', 'P-Value', 'log2FC', 'Pathways', 'Upstream', 'Diseases', 'Regulation']
                        sheets_data['Metabolites'] = pd.DataFrame(metabolite_data, columns=metabolite_columns)
            
            # 3. Export Upstream Regulators with metabolite list
            if hasattr(self, 'upstream_tree') and self.upstream_tree is not None:
                upstream_items = self.upstream_tree.get_children()
                if upstream_items:
                    upstream_data = []
                    for item in upstream_items:
                        values = list(self.upstream_tree.item(item)['values'])
                        # Get upstream name (column index 1)
                        upstream_name = values[1] if len(values) > 1 else None
                        # Get metabolite list from upstream_to_metabolites mapping
                        metabolite_list = ''
                        if upstream_name and hasattr(self, 'upstream_to_metabolites'):
                            mets = self.upstream_to_metabolites.get(upstream_name, set())
                            if mets:
                                metabolite_list = ', '.join(sorted(mets))
                        values.append(metabolite_list)
                        upstream_data.append(values)
                    upstream_columns = list(self.upstream_tree['columns']) + ['Metabolites']
                    sheets_data['Upstream_Regulators'] = pd.DataFrame(upstream_data, columns=upstream_columns)
            
            # 4. Export Diseases with metabolite list
            if hasattr(self, 'disease_tree') and self.disease_tree is not None:
                disease_items = self.disease_tree.get_children()
                if disease_items:
                    disease_data = []
                    for item in disease_items:
                        values = list(self.disease_tree.item(item)['values'])
                        # Get disease name (column index 1)
                        disease_name = values[1] if len(values) > 1 else None
                        # Get metabolite list from disease_to_metabolites mapping
                        metabolite_list = ''
                        if disease_name and hasattr(self, 'disease_to_metabolites'):
                            mets = self.disease_to_metabolites.get(disease_name, set())
                            if mets:
                                metabolite_list = ', '.join(sorted(mets))
                        values.append(metabolite_list)
                        disease_data.append(values)
                    disease_columns = list(self.disease_tree['columns']) + ['Metabolites']
                    sheets_data['Diseases'] = pd.DataFrame(disease_data, columns=disease_columns)
            
            # Convert numeric columns from text to actual numbers
            def _convert_numeric_columns(df):
                """Convert scientific notation strings to proper numeric values."""
                for col in df.columns:
                    # Skip checkbox, name, status, regulation, type columns
                    if col in ['Select', 'Pathway', 'Disease', 'Name', 'Status', 'Regulation', 'Type', 'Metabolites']:
                        continue
                    try:
                        # Attempt to convert to numeric, handling scientific notation
                        df[col] = pd.to_numeric(df[col], errors='coerce')  # Changed from 'ignore' for type safety
                    except Exception:
                        pass
                return df
            
            # Apply numeric conversion to all sheets
            for name in sheets_data:
                sheets_data[name] = _convert_numeric_columns(sheets_data[name])
            
            # Save to Excel with multiple sheets
            with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
                for sheet_name, df in sheets_data.items():
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
            
            # Create statistics explanation text file in same directory
            output_dir = os.path.dirname(filepath)
            stats_guide_path = os.path.join(output_dir, 'STATISTICS_GUIDE.txt')
            self._create_statistics_explanation_file(stats_guide_path)
            
            # Log results
            total_sheets = len(sheets_data)
            self._log_network(f"✅ All tables exported to: {filepath}\n")
            for sheet_name, df in sheets_data.items():
                self._log_network(f"   • {sheet_name}: {len(df)} rows\n")
            
            messagebox.showinfo("Export Complete", 
                              f"Successfully exported {total_sheets} tables to:\n{filepath}\n\n" +
                              "\n".join([f"• {name}: {len(df)} rows" for name, df in sheets_data.items()]))
        
        except Exception as e:
            error_msg = f"Error exporting pathway table: {str(e)}"
            self._log_network(f"❌ {error_msg}\n")
            messagebox.showerror("Export Error", error_msg)
    
    # ============ TABLE UPDATE OPERATIONS ============
    
    def update_pathway_table(self):
        """Update the pathway table display"""
        for item in self.pathway_tree.get_children():
            self.pathway_tree.delete(item)
        
        if self.pathway_filtered_pathways_data is None:
            self.pathway_stats_label.config(text="No pathway data loaded")
            return
        
        # Rebuild rows, including ML columns and consensus rank
        self.pathway_selections = {}

        # Collect items to compute ranks
        items = []
        for pathway_name, stats in self.pathway_filtered_pathways_data.items():
            if isinstance(pathway_name, str) and pathway_name.startswith('_'):
                continue
            if not isinstance(stats, dict):
                continue
            try:
                zscore = float(stats.get('z_score', 0.0) or 0.0)
            except Exception:
                zscore = 0.0
            try:
                pvalue = float(stats.get('combined_pvalue', 1.0) or 1.0)
            except Exception:
                pvalue = 1.0
            n_metabolites = stats.get('n_metabolites', 0) or 0
            status = stats.get('status', 'No Change')
            try:
                ml_score = float(stats.get('ml_rwr_score', 0.0) or 0.0)
            except Exception:
                ml_score = 0.0
            try:
                pr = float(stats.get('ml_pagerank', 0.0) or 0.0)
            except Exception:
                pr = 0.0
            try:
                btw = float(stats.get('ml_betweenness', 0.0) or 0.0)
            except Exception:
                btw = 0.0

            items.append({
                'name': pathway_name,
                'status': status,
                'z': zscore,
                'p': pvalue,
                'n': n_metabolites,
                'ml': ml_score,
                'pr': pr,
                'btw': btw,
            })

        z_sorted = sorted(items, key=lambda x: abs(x['z']), reverse=True)
        ml_sorted = sorted(items, key=lambda x: x['ml'], reverse=True)
        rank_z = {it['name']: idx + 1 for idx, it in enumerate(z_sorted)}
        rank_ml = {it['name']: idx + 1 for idx, it in enumerate(ml_sorted)}

        for it in items:
            name = it['name']
            consensus = int(round((rank_z.get(name, 0) + rank_ml.get(name, 0)) / 2.0)) if items else 0
            self.pathway_selections[name] = False
            select_symbol = "☐"
            self.pathway_tree.insert('', 'end', values=(
                select_symbol,
                name,
                it['status'],
                f"{it['z']:.3f}",
                f"{it['p']:.2e}",
                it['n'],
                f"{it['ml']:.4f}",
                f"{it['pr']:.4f}",
                f"{it['btw']:.4f}",
                consensus
            ))

        total_pathways = len([k for k in self.pathway_filtered_pathways_data.keys() if not (isinstance(k, str) and k.startswith('_'))])
        self.pathway_stats_label.config(text=f"Total: {total_pathways}")
    
    def generate_network_visualization(self):
        """Generate interactive network and exports (HTML/GraphML), using current selections.

        Adapted from old GUI implementation to work within Pathway Network tab without external dependencies.
        """
        try:
            # Basic validations
            if self.pathway_filtered_metabolites_data is None or self.pathway_filtered_pathways_data is None:
                self.frame.after(0, lambda: messagebox.showerror("Error", "No pathway data available. Load data from the annotation tab first."))
                return

            # Ensure output folder
            base_output_dir = None
            if hasattr(self, 'network_output_folder') and self.network_output_folder.get():
                base_output_dir = self.network_output_folder.get()
            if not base_output_dir:
                folder_path = filedialog.askdirectory(title="Select Output Folder for Network Results")
                if not folder_path:
                    return
                if not hasattr(self, 'network_output_folder'):
                    self.network_output_folder = tk.StringVar()
                self.network_output_folder.set(folder_path)
                base_output_dir = folder_path

            # Enforce STRICT selections based on the current table state (checkmarks only)
            selected_names = []
            try:
                for iid in self.pathway_tree.get_children(''):
                    vals = self.pathway_tree.item(iid, 'values')
                    if len(vals) >= 2 and str(vals[0]).strip() in ('☑️', '✓'):
                        selected_names.append(vals[1])
            except Exception:
                pass

            # Build selected pathways dict strictly from selected names
            selected_pathways = {}
            source_dict = self.pathway_original_pathways_data or {}
            for name in selected_names:
                if name in source_dict and isinstance(source_dict[name], dict):
                    selected_pathways[name] = source_dict[name]
            if not selected_pathways:
                self.frame.after(0, lambda: messagebox.showwarning("Warning", "No pathways selected. Please select at least one pathway."))
                return
            
            # CRITICAL FIX: Recalculate pathway status using CURRENT z-score thresholds
            # This ensures network colors always reflect the latest threshold settings
            try:
                current_z_activation = float(self.network_zscore_activation.get())
                current_z_inhibition = float(self.network_zscore_inhibition.get())
                self._log_network(f"Applying current z-score thresholds to network:\n")
                self._log_network(f"  Activation: z >= {current_z_activation}\n")
                self._log_network(f"  Inhibition: z <= {current_z_inhibition}\n")
                
                # Update status for each selected pathway based on CURRENT thresholds
                for pathway_name, stats in selected_pathways.items():
                    if isinstance(stats, dict):
                        z_score = stats.get('z_score', 0.0)
                        # Recalculate status using current UI thresholds
                        if z_score >= current_z_activation:
                            stats['status'] = 'Activated'
                        elif z_score <= current_z_inhibition:
                            stats['status'] = 'Inhibited'
                        else:
                            stats['status'] = 'No Change'
            except (ValueError, AttributeError) as e:
                self._log_network(f"⚠️ Could not apply z-score thresholds: {e}\n")
                # Continue with original status values

            # Log start
            self._log_network("\n=== NETWORK GENERATION STARTED ===\n")
            self._log_network(f"Using {len(selected_pathways)} selected pathways\n")

            # Start with empty dataframe - we'll ONLY include metabolites connected to SELECTED items
            df_combined = pd.DataFrame()
            
            # CRITICAL: Get list of SIGNIFICANT metabolites (those passing filters in the Lipids tab)
            # We should ONLY include metabolites from pathway stats that are ALSO in the filtered data
            significant_metabolites = set()
            
            self._log_network(f"\n{'='*70}\n")
            self._log_network(f"🔍 DEBUG: IDENTIFYING SIGNIFICANT METABOLITES\n")
            self._log_network(f"{'='*70}\n")
            
            if self.pathway_filtered_metabolites_data is not None and not self.pathway_filtered_metabolites_data.empty:
                self._log_network(f"✓ pathway_filtered_metabolites_data exists\n")
                self._log_network(f"  Shape: {self.pathway_filtered_metabolites_data.shape}\n")
                self._log_network(f"  Columns: {self.pathway_filtered_metabolites_data.columns.tolist()}\n")
                
                id_col = self.feature_col if hasattr(self, 'feature_col') and self.feature_col in self.pathway_filtered_metabolites_data.columns else 'Name'
                self._log_network(f"  ID column to use: '{id_col}'\n")
                
                if id_col in self.pathway_filtered_metabolites_data.columns:
                    significant_metabolites = set(self.pathway_filtered_metabolites_data[id_col].tolist())
                    self._log_network(f"\n✅ SIGNIFICANT METABOLITES EXTRACTED: {len(significant_metabolites)}\n")
                    self._log_network(f"   First 10 IDs: {list(significant_metabolites)[:10]}\n")
                    
                    # Show p-value range if available
                    if self.pvalue_col and self.pvalue_col in self.pathway_filtered_metabolites_data.columns:
                        pvals = self.pathway_filtered_metabolites_data[self.pvalue_col]
                        self._log_network(f"   P-value range: {pvals.min():.6f} to {pvals.max():.6f}\n")
                else:
                    self._log_network(f"❌ ERROR: ID column '{id_col}' not found in filtered data!\n")
            else:
                self._log_network(f"❌ ERROR: pathway_filtered_metabolites_data is None or empty!\n")
            
            self._log_network(f"\n{'='*70}\n")
            self._log_network(f"🔍 DEBUG: COLLECTING METABOLITES FROM SELECTED PATHWAYS\n")
            self._log_network(f"{'='*70}\n")
            
            # STEP 1: Collect metabolites from SELECTED pathways (intersected with significant metabolites)
            pathway_metabolites = set()
            for pathway_name, stats in selected_pathways.items():
                # CRITICAL: Use 'hits' (significant members) instead of 'metabolites' (all measured members)
                # 'hits' contains only metabolites that passed significance threshold in Fisher ORA
                sig_mets = stats.get('hits', stats.get('metabolites', []))  # Fallback to 'metabolites' for compatibility
                all_mets = stats.get('metabolites', [])
                
                self._log_network(f"\n📊 Pathway: '{pathway_name}'\n")
                self._log_network(f"   Total measured metabolites: {len(all_mets)}\n")
                self._log_network(f"   Significant metabolites (hits): {len(sig_mets)}\n")
                self._log_network(f"   First 5 significant IDs: {sig_mets[:5]}\n")
                
                pathway_metabolites.update(sig_mets)
                self._log_network(f"   ✓ Added {len(sig_mets)} SIGNIFICANT metabolites to pathway_metabolites set\n")
            
            self._log_network(f"\n{'='*70}\n")
            self._log_network(f"📊 SUMMARY: PATHWAY METABOLITES\n")
            self._log_network(f"{'='*70}\n")
            self._log_network(f"✅ Total significant metabolites from selected pathways: {len(pathway_metabolites)}\n")
            self._log_network(f"   First 10 IDs: {list(pathway_metabolites)[:10]}\n")
            
            # CRITICAL: Add metabolites from selected diseases AND upstream regulators
            # This ensures ALL connected metabolites appear in network regardless of pathway connections
            include_diseases = getattr(self, 'network_include_diseases', tk.BooleanVar(value=True)).get()
            include_upstream = getattr(self, 'network_include_upstream', tk.BooleanVar(value=True)).get()
            
            all_required_metabolites = set(pathway_metabolites)  # Start with pathway metabolites
            disease_metabolites = set()  # Initialize
            upstream_metabolites = set()  # Initialize
            
            # STEP 2: Collect metabolites from SELECTED diseases (intersected with significant metabolites)
            if include_diseases and hasattr(self, 'disease_to_metabolites') and self.disease_to_metabolites:
                selected_disease_count = 0
                for disease_name, disease_mets in self.disease_to_metabolites.items():
                    is_selected = self.disease_selections.get(disease_name, False)
                    if is_selected:
                        selected_disease_count += 1
                        # CRITICAL: Only include metabolites that are ALSO in the significant set
                        sig_disease_mets = [m for m in disease_mets if m in significant_metabolites] if significant_metabolites else disease_mets
                        disease_metabolites.update(sig_disease_mets)
                        self._log_network(f"   ✓ Including metabolites from selected disease: {disease_name} ({len(sig_disease_mets)}/{len(disease_mets)} significant)\n")
                
                all_required_metabolites.update(disease_metabolites)
                
                self._log_network(f"\n🔍 SELECTED DISEASE METABOLITES:\n")
                self._log_network(f"   Selected diseases: {selected_disease_count} / {len(self.disease_to_metabolites)}\n")
                self._log_network(f"   Additional significant metabolites from diseases: {len(disease_metabolites)}\n")
                if selected_disease_count == 0:
                    self._log_network(f"   ⚠️ NO DISEASES SELECTED - disease metabolites will not be added\n")
                    self._log_network(f"   💡 TIP: Double-click diseases in Disease table to select them before generating network\n\n")
            
            # STEP 3: Collect metabolites from SELECTED upstream regulators (intersected with significant metabolites)
            if include_upstream and hasattr(self, 'upstream_to_metabolites') and self.upstream_to_metabolites:
                selected_upstream_count = 0
                for upstream_name, upstream_mets in self.upstream_to_metabolites.items():
                    is_selected = self.upstream_selections.get(upstream_name, False)
                    if is_selected:
                        selected_upstream_count += 1
                        # CRITICAL: Only include metabolites that are ALSO in the significant set
                        sig_upstream_mets = [m for m in upstream_mets if m in significant_metabolites] if significant_metabolites else upstream_mets
                        upstream_metabolites.update(sig_upstream_mets)
                        self._log_network(f"   ✓ Including metabolites from selected upstream: {upstream_name} ({len(sig_upstream_mets)}/{len(upstream_mets)} significant)\n")
                
                all_required_metabolites.update(upstream_metabolites)
                
                self._log_network(f"\n🔍 SELECTED UPSTREAM METABOLITES:\n")
                self._log_network(f"   Selected upstream regulators: {selected_upstream_count} / {len(self.upstream_to_metabolites)}\n")
                self._log_network(f"   Additional significant metabolites from upstream: {len(upstream_metabolites)}\n")
            
            # STEP 4: Build df_combined with ONLY the required metabolites
            self._log_network(f"\n{'='*70}\n")
            self._log_network(f"📊 FINAL METABOLITE COUNT FOR NETWORK\n")
            self._log_network(f"{'='*70}\n")
            self._log_network(f"✅ TOTAL SIGNIFICANT METABOLITES TO INCLUDE: {len(all_required_metabolites)}\n")
            self._log_network(f"   From pathways: {len(pathway_metabolites)}\n")
            self._log_network(f"   Additional from diseases: {len(disease_metabolites) if include_diseases else 0}\n")
            self._log_network(f"   Additional from upstream: {len(upstream_metabolites) if include_upstream else 0}\n")
            self._log_network(f"   Combined unique metabolites: {len(all_required_metabolites)}\n")
            self._log_network(f"   First 10 IDs: {list(all_required_metabolites)[:10]}\n\n")
            
            # CRITICAL FIX: Use FILTERED data (significant metabolites only), NOT original data
            # This ensures network only shows metabolites that passed the p-value threshold
            source_df = self.pathway_filtered_metabolites_data
            
            self._log_network(f"{'='*70}\n")
            self._log_network(f"🔍 DEBUG: BUILDING df_combined FROM SOURCE DATA\n")
            self._log_network(f"{'='*70}\n")
            
            if source_df is not None and not source_df.empty:
                self._log_network(f"✓ Source DataFrame exists (pathway_filtered_metabolites_data)\n")
                self._log_network(f"  Shape: {source_df.shape}\n")
            if source_df is not None and not source_df.empty:
                self._log_network(f"✓ Source DataFrame exists (pathway_filtered_metabolites_data)\n")
                self._log_network(f"  Shape: {source_df.shape}\n")
                # Use cached feature column for filtering
                id_col = self.feature_col if hasattr(self, 'feature_col') and self.feature_col in source_df.columns else 'Name'
                
                self._log_network(f"\n🔍 DataFrame filtering debug:\n")
                self._log_network(f"   Source DF shape: {source_df.shape}\n")
                self._log_network(f"   Using ID column: '{id_col}'\n")
                self._log_network(f"   Available columns: {source_df.columns.tolist()}\n")
                self._log_network(f"   Sample IDs in source: {source_df[id_col].head(10).tolist() if id_col in source_df.columns else 'N/A'}\n")
                self._log_network(f"   Required metabolites count: {len(all_required_metabolites)}\n")
                self._log_network(f"   Sample required IDs: {list(all_required_metabolites)[:10]}\n")
                
                if id_col in source_df.columns:
                    # Filter to only required metabolites
                    df_combined = source_df[source_df[id_col].isin(all_required_metabolites)].copy()
                    self._log_network(f"\n✅ FILTERED df_combined:\n")
                    self._log_network(f"   Rows in df_combined: {len(df_combined)}\n")
                    self._log_network(f"   Expected: {len(all_required_metabolites)}\n")
                    self._log_network(f"   Match: {'✓ YES' if len(df_combined) == len(all_required_metabolites) else f'❌ NO (difference: {len(all_required_metabolites) - len(df_combined)})'}\n")
                    self._log_network(f"   IDs in df_combined: {df_combined[id_col].tolist()[:10]}\n")
                    
                    # Check for metabolites not found in source data
                    found_metabolites = set(df_combined[id_col].tolist())
                    missing_metabolites = all_required_metabolites - found_metabolites
                    
                    if missing_metabolites:
                        self._log_network(f"\n⚠️ WARNING: {len(missing_metabolites)} metabolites not found in source data:\n")
                        for met_name in list(missing_metabolites)[:10]:
                            self._log_network(f"   - {met_name}\n")
                        for met_name in missing_metabolites:
                            new_row = {id_col: met_name}
                            # Use cached column names for other fields
                            if hasattr(self, 'pvalue_col'):
                                new_row[self.pvalue_col] = 1.0
                            if hasattr(self, 'log2fc_col'):
                                new_row[self.log2fc_col] = 0.0
                            new_row['Regulation'] = 'Unknown'
                            
                            df_combined = pd.concat([df_combined, pd.DataFrame([new_row])], ignore_index=True)
                            self._log_network(f"   ⚠️ {met_name} (p=1.0, log2FC=0.0)\n")
                else:
                    self._log_network(f"❌ ERROR: Column '{id_col}' not found in source data. Available columns: {source_df.columns.tolist()}\n")
                    df_combined = pd.DataFrame()
            else:
                self._log_network(f"⚠️ WARNING: No source data available - creating df_combined with default values\n")
                rows = []
                for met_name in all_required_metabolites:
                    rows.append({
                        'Name': met_name,
                        'pvalue': 1.0,
                        'log2FC': 0.0,
                        'Regulation': 'Unknown'
                    })
                df_combined = pd.DataFrame(rows)
            
            self._log_network(f"\n{'='*70}\n")
            self._log_network(f"✅ FINAL df_combined READY FOR NETWORK GENERATION\n")
            self._log_network(f"{'='*70}\n")
            self._log_network(f"   Shape: {df_combined.shape}\n")
            self._log_network(f"   Columns: {df_combined.columns.tolist()}\n")
            if not df_combined.empty:
                self._log_network(f"   First 5 rows:\n")
                for idx, row in df_combined.head(5).iterrows():
                    lipid_id = row.get(id_col, 'N/A')
                    pval = row.get(self.pvalue_col, 'N/A') if hasattr(self, 'pvalue_col') else 'N/A'
                    log2fc = row.get(self.log2fc_col, 'N/A') if hasattr(self, 'log2fc_col') else 'N/A'
                    self._log_network(f"      {lipid_id}: p={pval}, log2FC={log2fc}\n")
            
            # ========================================================
            # CRITICAL: ENSURE STANDARDIZED COLUMNS FOR DOWNSTREAM
            # The network generation expects: 'Name', 'pvalue', 'log2FC'
            # ========================================================
            if not df_combined.empty:
                self._log_network(f"\n🔄 ENSURING STANDARDIZED COLUMNS:\n")
                
                # Standardize 'Name' column (copy from id_col if different)
                if 'Name' not in df_combined.columns:
                    if id_col and id_col in df_combined.columns:
                        df_combined['Name'] = df_combined[id_col]
                        self._log_network(f"   ✓ Created 'Name' from '{id_col}'\n")
                    else:
                        # Try other candidates
                        for candidate in ['LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite']:
                            if candidate in df_combined.columns:
                                df_combined['Name'] = df_combined[candidate]
                                self._log_network(f"   ✓ Created 'Name' from '{candidate}'\n")
                                break
                else:
                    self._log_network(f"   ✓ 'Name' column already exists\n")
                
                # Standardize 'pvalue' column
                if 'pvalue' not in df_combined.columns:
                    pval_col = self.pvalue_col if hasattr(self, 'pvalue_col') else None
                    if pval_col and pval_col in df_combined.columns:
                        df_combined['pvalue'] = df_combined[pval_col]
                        self._log_network(f"   ✓ Created 'pvalue' from '{pval_col}'\n")
                    else:
                        for candidate in ['p_value', 'P-Value', 'p-value', 'adj_p']:
                            if candidate in df_combined.columns:
                                df_combined['pvalue'] = df_combined[candidate]
                                self._log_network(f"   ✓ Created 'pvalue' from '{candidate}'\n")
                                break
                else:
                    self._log_network(f"   ✓ 'pvalue' column already exists\n")
                
                # Standardize 'log2FC' column
                if 'log2FC' not in df_combined.columns:
                    log2_col = self.log2fc_col if hasattr(self, 'log2fc_col') else None
                    if log2_col and log2_col in df_combined.columns:
                        df_combined['log2FC'] = df_combined[log2_col]
                        self._log_network(f"   ✓ Created 'log2FC' from '{log2_col}'\n")
                    else:
                        for candidate in ['log2_FC', 'Log2FC', 'log2FoldChange', 'FC']:
                            if candidate in df_combined.columns:
                                df_combined['log2FC'] = df_combined[candidate]
                                self._log_network(f"   ✓ Created 'log2FC' from '{candidate}'\n")
                                break
                else:
                    self._log_network(f"   ✓ 'log2FC' column already exists\n")
                
                # Verify standardization
                has_name = 'Name' in df_combined.columns
                has_pvalue = 'pvalue' in df_combined.columns
                has_log2fc = 'log2FC' in df_combined.columns
                self._log_network(f"   📋 Standardized columns present: Name={has_name}, pvalue={has_pvalue}, log2FC={has_log2fc}\n")
            
            self._log_network(f"{'='*70}\n\n")

            # Build upstream data from selections
            upstream_data = None
            selected_upstream_names = set()
            if include_upstream and hasattr(self, 'original_upstream_data') and self.original_upstream_data:
                # Collect strictly from checked rows in the upstream table
                selected_upstream = set()
                try:
                    for iid in self.upstream_tree.get_children(''):
                        vals = self.upstream_tree.item(iid, 'values')
                        if len(vals) >= 3 and str(vals[0]).strip() in self.CHECKED_SYMBOLS:
                            # values: [check, name, type, count]
                            selected_upstream.add((vals[1], vals[2]))
                except Exception:
                    pass
                
                # DEBUG: Log upstream selection status
                self._log_network(f"\n🔍 UPSTREAM SELECTION STATUS:\n")
                self._log_network(f"   Include upstream: {include_upstream}\n")
                self._log_network(f"   Total upstream available: {len(self.original_upstream_data)}\n")
                self._log_network(f"   Upstream with checkboxes: {len(self.upstream_tree.get_children(''))}\n")
                self._log_network(f"   Selected upstream (checked ☑️): {len(selected_upstream)}\n")
                if selected_upstream:
                    for (name, reg_type) in list(selected_upstream)[:10]:  # Show first 10
                        self._log_network(f"      ✓ {name} ({reg_type})\n")
                else:
                    self._log_network(f"   ⚠️ NO UPSTREAM SELECTED - Network will not include any upstream nodes!\n")
                    self._log_network(f"   💡 TIP: Click 'Select All' in the Upstream table or double-click individual items\n\n")
                
                selected_upstream_names = {name for (name, _) in selected_upstream}
                rows = []
                for (name, reg_type) in selected_upstream:
                    info = self.original_upstream_data.get(name) or self.original_upstream_data.get(f"{reg_type}_{name}")
                    # Our stored dict uses name as key; ensure we have metabolites
                    if isinstance(info, dict):
                        mets = sorted(set(info.get('metabolites', [])))
                        rows.append({
                            'Name': name,
                            'Type': reg_type,
                            'Associated_Metabolites': ' | '.join(mets),
                            '# Compounds': len(mets)
                        })
                if rows:
                    upstream_data = pd.DataFrame(rows)

            # Build diseases data from selections
            disease_data = None
            selected_disease_names = set()
            if include_diseases and hasattr(self, 'original_disease_data') and self.original_disease_data:
                # Collect strictly from checked rows in the disease table
                selected_diseases = []
                try:
                    for iid in self.disease_tree.get_children(''):
                        vals = self.disease_tree.item(iid, 'values')
                        if len(vals) >= 2 and str(vals[0]).strip() in self.CHECKED_SYMBOLS:
                            selected_diseases.append(vals[1])
                except Exception:
                    pass
                
                # DEBUG: Log disease selection status
                self._log_network(f"\n🔍 DISEASE SELECTION STATUS:\n")
                self._log_network(f"   Include diseases: {include_diseases}\n")
                self._log_network(f"   Total diseases available: {len(self.original_disease_data)}\n")
                self._log_network(f"   Diseases with checkboxes: {len(self.disease_tree.get_children(''))}\n")
                self._log_network(f"   Selected diseases (checked ☑️): {len(selected_diseases)}\n")
                if selected_diseases:
                    for disease_name in selected_diseases:
                        mets = self.original_disease_data.get(disease_name, [])
                        self._log_network(f"      ✓ {disease_name} ({len(mets)} metabolites)\n")
                else:
                    self._log_network(f"   ⚠️ NO DISEASES SELECTED - Network will not include any disease nodes!\n")
                    self._log_network(f"   💡 TIP: Double-click diseases in the Disease table to select them\n\n")
                
                rows = []
                for disease_name in selected_diseases:
                    mets = self.original_disease_data.get(disease_name, [])
                    smets = sorted(set(mets))
                    rows.append({
                        'Disease Name': disease_name,
                        'Associated_Metabolites': ' | '.join(smets),
                        '# Compounds': len(smets)
                    })
                if rows:
                    disease_data = pd.DataFrame(rows)
                selected_disease_names = set(selected_diseases)

            # Create figure using shared network module
            from main_script.metabolite_pathway_network import (
                create_interactive_metabolite_pathway_plot_from_stats,
                export_network_to_cytoscape
            )

            # Reasonable defaults for sizes
            network_width = 1400
            network_height = 900

            # Build outputs based on user selection (validation already done in start_generate_network)
            gen_network = getattr(self, 'network_gen_network', tk.BooleanVar(value=True)).get()
            gen_plots = getattr(self, 'network_gen_plots', tk.BooleanVar(value=True)).get()

            fig = None
            if gen_network:
                fig = create_interactive_metabolite_pathway_plot_from_stats(
                df_combined,
                selected_pathways,
                node_spacing=1.8,
                layout_iterations=200,
                metabolite_font_size=10,
                pathway_font_size=9,
                enzyme_font_size=8,
                metabolite_max_chars=10,
                metabolite_max_lines=2,
                pathway_max_chars=20,
                pathway_max_lines=3,
                network_width=network_width,
                network_height=network_height,
                include_upstream=include_upstream,
                include_diseases=include_diseases,
                upstream_data=upstream_data,
                disease_data=disease_data
                )

            # Prepare output folder with timestamp
            timestamp = time.strftime("%Y%m%d_%H%M%S")
            pathways_networks_dir = os.path.join(base_output_dir, f'Pathways_networks_{timestamp}')
            os.makedirs(pathways_networks_dir, exist_ok=True)
            self.pathways_networks_dir = pathways_networks_dir

            # File names
            # Save network outputs
            if gen_network and fig is not None:
                graphml_path = os.path.join(pathways_networks_dir, 'network.graphml')
                legend_png_path = os.path.join(pathways_networks_dir, 'network_legend.png')

                # Export GraphML / Cytoscape (only if checkbox enabled)
                if hasattr(self, 'export_cytoscape') and self.export_cytoscape.get():
                    try:
                        export_network_to_cytoscape(
                            df_combined,
                            selected_pathways,
                            output_file=graphml_path,
                            include_upstream=include_upstream,
                            include_diseases=include_diseases,
                            upstream_data=upstream_data,
                            disease_data=disease_data
                        )
                        self._log_network(f"Exported GraphML/Cytoscape: {graphml_path}\n")
                    except Exception as e:
                        self._log_network(f"GraphML export failed: {e}\n")
                else:
                    self._log_network(f"Cytoscape export skipped (checkbox unchecked)\n")

                # Always generate a PNG legend for the network (no HTML fallback)
                try:
                    from main_script.metabolite_pathway_network import create_network_legend_figure  # type: ignore
                    leg_fig = create_network_legend_figure(include_upstream=include_upstream, include_diseases=include_diseases)
                    leg_fig.write_image(legend_png_path, format='png', width=520, height=360)
                    self._log_network(f"Saved Legend PNG: {legend_png_path}\n")
                except Exception as e:
                    self._log_network(f"Legend PNG export failed: {e}\n")
                    self._log_network(f"Note: Install kaleido for PNG export: pip install kaleido\n")
                
                # Create statistics explanation text file
                stats_guide_path = os.path.join(pathways_networks_dir, 'STATISTICS_GUIDE.txt')
                self._create_statistics_explanation_file(stats_guide_path)

            # Generate statistical plots if requested
            if gen_plots:
                try:
                    from main_script.metabolite_pathway_network import create_pathway_summary_plots
                    import math

                    # Build metabolites dict from df_combined
                    metabolites = {}
                    for _, r in df_combined.iterrows():
                        met_name = self._get_metabolite_value(r, 'name', 'Unknown')
                        metabolites[met_name] = {
                            'log2FC': self._as_float(self._get_metabolite_value(r, 'log2fc', 0.0), 0.0),
                            'pvalue': self._as_float(self._get_metabolite_value(r, 'pvalue', 1.0), 1.0)
                        }

                    # Generate plots with SELECTED pathways (for states, significance, selected enrichment)
                    logger.info(f"[PLOT DEBUG] Generating plots with {len(selected_pathways)} SELECTED pathways")
                    figs_selected = create_pathway_summary_plots(
                        metabolites,
                        selected_pathways,
                        font_size_scale=self.pathway_font_size_scale.get(),
                        bargraph_pathway_max_chars=getattr(self, 'bargraph_pathway_max_chars', tk.IntVar(value=40)).get(),
                        bargraph_pathway_max_lines=getattr(self, 'bargraph_pathway_max_lines', tk.IntVar(value=3)).get(),
                        bargraph_title_font_size=getattr(self, 'bargraph_title_font_size', tk.IntVar(value=22)).get(),
                        bargraph_axis_title_font_size=getattr(self, 'bargraph_axis_title_font_size', tk.IntVar(value=20)).get(),
                        bargraph_tick_font_size=getattr(self, 'bargraph_tick_font_size', tk.IntVar(value=15)).get(),
                        bargraph_pathway_label_font_size=getattr(self, 'bargraph_pathway_label_font_size', tk.IntVar(value=20)).get(),
                        bargraph_width_px=getattr(self, 'bargraph_width_px', tk.IntVar(value=1200)).get(),
                        bargraph_height_px=getattr(self, 'bargraph_height_px', tk.IntVar(value=700)).get(),
                        bargraph_legend_title_font_size=getattr(self, 'bargraph_legend_title_font_size', tk.IntVar(value=14)).get(),
                        bargraph_legend_tick_font_size=getattr(self, 'bargraph_legend_tick_font_size', tk.IntVar(value=12)).get(),
                        bargraph_legend_thickness=getattr(self, 'bargraph_legend_thickness', tk.IntVar(value=10)).get(),
                        bargraph_legend_length_percent=getattr(self, 'bargraph_legend_length_percent', tk.IntVar(value=70)).get(),
                        chord_metabolite_max_chars=getattr(self, 'chord_metabolite_max_chars', tk.IntVar(value=40)).get(),
                        chord_metabolite_max_lines=getattr(self, 'chord_metabolite_max_lines', tk.IntVar(value=2)).get(),
                        chord_pathway_max_chars=getattr(self, 'chord_pathway_max_chars', tk.IntVar(value=40)).get(),
                        chord_pathway_max_lines=getattr(self, 'chord_pathway_max_lines', tk.IntVar(value=2)).get(),
                        chord_metabolite_font_size=getattr(self, 'chord_metabolite_font_size', tk.IntVar(value=12)).get(),
                        chord_pathway_font_size=getattr(self, 'chord_pathway_font_size', tk.IntVar(value=12)).get(),
                        chord_width_inches=getattr(self, 'chord_figure_width_inches', tk.IntVar(value=12)).get(),
                        chord_height_inches=getattr(self, 'chord_figure_height_inches', tk.IntVar(value=12)).get(),
                        chord_keep_metabolite_case=getattr(self, 'chord_keep_metabolite_case', tk.BooleanVar(value=False)).get(),
                    )
                    # figs_selected now includes chord legend as 8th element
                    
                    # Extract plots that should use SELECTED pathways
                    fig1 = figs_selected[0] if len(figs_selected) > 0 else None  # Pathway States
                    fig2 = figs_selected[1] if len(figs_selected) > 1 else None  # Pathway Significance (Vertical)
                    fig3 = figs_selected[2] if len(figs_selected) > 2 else None  # Pathway Significance (Horizontal)
                    fig4 = figs_selected[3] if len(figs_selected) > 3 else None  # Bubble plot (may be disabled)
                    fig6 = figs_selected[5] if len(figs_selected) > 5 else None  # Pathway Enrichment (Selected)
                    fig7 = figs_selected[6] if len(figs_selected) > 6 else None  # Chord Diagram (Selected)
                    fig_legend = figs_selected[7] if len(figs_selected) > 7 else None  # Chord Legend
                    
                    # Generate Top 20 enrichment plot with ALL pathways
                    all_pathways_for_top20 = self.pathway_original_pathways_data or selected_pathways
                    logger.info(f"[PLOT DEBUG] Generating Top 20 enrichment with {len(all_pathways_for_top20)} ALL pathways")
                    
                    fig5 = None  # Pathway Enrichment (Top 20)
                    if len(all_pathways_for_top20) > len(selected_pathways):
                        # Generate plots with ALL pathways only to extract Top 20 enrichment
                        figs_all = create_pathway_summary_plots(
                            metabolites,
                            all_pathways_for_top20,
                            font_size_scale=self.pathway_font_size_scale.get(),
                            bargraph_pathway_max_chars=getattr(self, 'bargraph_pathway_max_chars', tk.IntVar(value=40)).get(),
                            bargraph_pathway_max_lines=getattr(self, 'bargraph_pathway_max_lines', tk.IntVar(value=3)).get(),
                            bargraph_title_font_size=getattr(self, 'bargraph_title_font_size', tk.IntVar(value=22)).get(),
                            bargraph_axis_title_font_size=getattr(self, 'bargraph_axis_title_font_size', tk.IntVar(value=20)).get(),
                            bargraph_tick_font_size=getattr(self, 'bargraph_tick_font_size', tk.IntVar(value=15)).get(),
                            bargraph_pathway_label_font_size=getattr(self, 'bargraph_pathway_label_font_size', tk.IntVar(value=20)).get(),
                            bargraph_width_px=getattr(self, 'bargraph_width_px', tk.IntVar(value=1200)).get(),
                            bargraph_height_px=getattr(self, 'bargraph_height_px', tk.IntVar(value=700)).get(),
                            bargraph_legend_title_font_size=getattr(self, 'bargraph_legend_title_font_size', tk.IntVar(value=14)).get(),
                            bargraph_legend_tick_font_size=getattr(self, 'bargraph_legend_tick_font_size', tk.IntVar(value=12)).get(),
                            bargraph_legend_thickness=getattr(self, 'bargraph_legend_thickness', tk.IntVar(value=10)).get(),
                            bargraph_legend_length_percent=getattr(self, 'bargraph_legend_length_percent', tk.IntVar(value=70)).get(),
                            chord_metabolite_max_chars=getattr(self, 'chord_metabolite_max_chars', tk.IntVar(value=40)).get(),
                            chord_metabolite_max_lines=getattr(self, 'chord_metabolite_max_lines', tk.IntVar(value=2)).get(),
                            chord_pathway_max_chars=getattr(self, 'chord_pathway_max_chars', tk.IntVar(value=40)).get(),
                            chord_pathway_max_lines=getattr(self, 'chord_pathway_max_lines', tk.IntVar(value=2)).get(),
                            chord_metabolite_font_size=getattr(self, 'chord_metabolite_font_size', tk.IntVar(value=12)).get(),
                            chord_pathway_font_size=getattr(self, 'chord_pathway_font_size', tk.IntVar(value=12)).get(),
                            chord_width_inches=getattr(self, 'chord_figure_width_inches', tk.IntVar(value=12)).get(),
                            chord_height_inches=getattr(self, 'chord_figure_height_inches', tk.IntVar(value=12)).get(),
                            chord_keep_metabolite_case=getattr(self, 'chord_keep_metabolite_case', tk.BooleanVar(value=False)).get(),
                        )
                        # figs_all now includes chord legend as 8th element
                        fig5 = figs_all[4] if len(figs_all) > 4 else None
                        logger.info(f"[PLOT DEBUG] Top 20 enrichment uses {len(all_pathways_for_top20)} pathways")
                    else:
                        # If no additional pathways, use selected pathways for Top 20 as well
                        fig5 = figs_selected[4] if len(figs_selected) > 4 else None
                        logger.info(f"[PLOT DEBUG] Top 20 enrichment uses {len(selected_pathways)} selected pathways (no additional pathways available)")

                    def _color_name(c):
                        """Classify a color string as 'red' or 'blue' when using hex gradients.

                        - If literal names present, return them.
                        - If hex, parse RGB and compare red vs blue channels.
                        - Otherwise, return raw string for upstream handling.
                        """
                        s = str(c).strip().lower()
                        if 'red' in s:
                            return 'red'
                        if 'blue' in s:
                            return 'blue'
                        # Hex like #rrggbb
                        if s.startswith('#') and len(s) == 7:
                            try:
                                r = int(s[1:3], 16)
                                g = int(s[3:5], 16)
                                b = int(s[5:7], 16)
                                return 'red' if r >= b else 'blue'
                            except Exception:
                                return s
                        return s

                    def _postprocess_vertical_states(fig):
                        # Ensure blue (inhibited) bars face down, red (activated) face up and reds come first
                        if fig is None or not getattr(fig, 'data', None):
                            return fig
                        tr = fig.data[0]
                        xs = list(getattr(tr, 'x', []))
                        ys = list(getattr(tr, 'y', []))
                        cols = getattr(getattr(tr, 'marker', None), 'color', [])
                        if not isinstance(cols, (list, tuple)):
                            cols = [cols] * len(xs)
                        print("[DEBUG] Vertical States: initial colors=", cols[:10])
                        items = []
                        for i, name in enumerate(xs):
                            classified = _color_name(cols[i])
                            orig_col = cols[i]
                            val = ys[i] if i < len(ys) else 0
                            try:
                                v = float(val)
                            except Exception:
                                v = 0.0
                            if classified == 'blue':
                                v = -abs(v)
                            elif classified == 'red':
                                v = abs(v)
                            else:
                                # derive from sign
                                v = abs(v) if v >= 0 else -abs(v)
                                classified = 'red' if v >= 0 else 'blue'
                            items.append((0 if classified == 'red' else 1, name, v, orig_col))
                        # Reds first then blues
                        items.sort(key=lambda t: (t[0], -abs(t[2])))
                        new_x = [it[1] for it in items]
                        new_y = [it[2] for it in items]
                        new_c = [it[3] for it in items]
                        print("[DEBUG] Vertical States: reordered first colors=", new_c[:10])
                        tr.x = new_x
                        tr.y = new_y
                        tr.marker.color = new_c
                        tr.orientation = 'v'
                        fig.update_layout(
                            barmode='relative',
                            yaxis=dict(zeroline=True, zerolinewidth=1, zerolinecolor='#333'),
                            margin=dict(b=160, l=80, r=40, t=60)
                        )
                        return fig

                    def _postprocess_vertical_significance(fig):
                        # Sort reds first then blues based on colors on the first trace.
                        if fig is None or not getattr(fig, 'data', None):
                            return fig
                        tr = fig.data[0]
                        xs = list(getattr(tr, 'x', []))
                        ys = list(getattr(tr, 'y', []))
                        cols = getattr(getattr(tr, 'marker', None), 'color', [])
                        if not isinstance(cols, (list, tuple)):
                            cols = [cols] * len(xs)
                        print("[DEBUG] Vertical Sig: initial colors=", cols[:10])
                        items = []
                        for i, name in enumerate(xs):
                            classified = _color_name(cols[i])
                            orig_col = cols[i]
                            val = ys[i] if i < len(ys) else 0
                            try:
                                v = abs(float(val))
                            except Exception:
                                v = 0.0
                            items.append((0 if classified == 'red' else 1, name, v, orig_col))
                        items.sort(key=lambda t: (t[0], -t[2]))
                        tr.x = [it[1] for it in items]
                        tr.y = [it[2] for it in items]
                        tr.marker.color = [it[3] for it in items]
                        print("[DEBUG] Vertical Sig: reordered first colors=", tr.marker.color[:10])
                        tr.orientation = 'v'
                        fig.update_layout(margin=dict(b=160, l=80, r=40, t=60))
                        return fig

                    def _postprocess_horizontal_significance(fig):
                        # Ensure horizontal bars with correct -log10(p) to the right and reds first
                        if fig is None or not getattr(fig, 'data', None):
                            return fig
                        tr = fig.data[0]
                        ys = list(getattr(tr, 'y', []))
                        xs = list(getattr(tr, 'x', []))
                        cols = getattr(getattr(tr, 'marker', None), 'color', [])
                        if not isinstance(cols, (list, tuple)):
                            cols = [cols] * len(ys)
                        print("[DEBUG] Horizontal Sig: initial colors=", cols[:10])
                        items = []
                        for i, name in enumerate(ys):
                            classified = _color_name(cols[i])
                            orig_col = cols[i]
                            val = xs[i] if i < len(xs) else 0
                            try:
                                v = abs(float(val))
                            except Exception:
                                v = 0.0
                            items.append((0 if classified == 'red' else 1, name, v, orig_col))
                        items.sort(key=lambda t: (t[0], -t[2]))
                        tr.y = [it[1] for it in items]
                        tr.x = [it[2] for it in items]
                        tr.marker.color = [it[3] for it in items]
                        print("[DEBUG] Horizontal Sig: reordered first colors=", tr.marker.color[:10])
                        tr.orientation = 'h'
                        xmax = max([it[2] for it in items], default=1.0)
                        fig.update_layout(
                            margin=dict(l=260, r=40, t=60, b=40),
                            xaxis=dict(range=[0, xmax * 1.1], domain=[0.5, 1.0])  # >=50% area to bars
                        )
                        return fig

                    # Apply postprocessing tweaks
                    if fig1 is not None and self.network_plot_selections['pathway_states'].get():
                        fig1 = _postprocess_vertical_states(fig1)
                    if fig2 is not None and self.network_plot_selections['pathway_sig_v'].get():
                        fig2 = _postprocess_vertical_significance(fig2)
                    if fig3 is not None and self.network_plot_selections['pathway_sig_h'].get():
                        fig3 = _postprocess_horizontal_significance(fig3)
                    # fig5, fig6, and fig7 (Pathway Enrichment and Chord Diagram) don't need postprocessing

                    # Choose which plots to save
                    sel = self.network_plot_selections
                    plots = []
                    if sel['pathway_states'].get() and fig1 is not None:
                        plots.append((fig1, 'pathway_states.png', 'Pathway States'))
                    if sel['pathway_sig_v'].get() and fig2 is not None:
                        plots.append((fig2, 'pathway_significance_vertical.png', 'Pathway Significance (Vertical)'))
                    if sel['pathway_sig_h'].get() and fig3 is not None:
                        plots.append((fig3, 'pathway_significance_horizontal.png', 'Pathway Significance (Horizontal)'))
                    if sel['pathway_enrichment'].get() and fig5 is not None:
                        plots.append((fig5, 'pathway_enrichment_top20.png', 'Pathway Enrichment (Top 20)'))
                    if sel['pathway_enrichment_selected'].get() and fig6 is not None:
                        plots.append((fig6, 'pathway_enrichment_selected.png', 'Pathway Enrichment (Selected Pathways)'))
                    if sel['chord_diagram'].get() and fig7 is not None:
                        plots.append((fig7, 'pathway_chord_diagram.png', 'Pathway-Metabolite Chord Diagram'))
                    if sel['chord_diagram'].get() and fig_legend is not None:
                        plots.append((fig_legend, 'pathway_chord_legend.png', 'Chord Diagram Legend'))

                    saved = 0
                    for pfig, fname, desc in plots:
                        try:
                            outp = os.path.join(pathways_networks_dir, fname)
                            
                            # Check if it's a matplotlib figure (fig7 chord diagram)
                            if hasattr(pfig, 'savefig'):
                                # Matplotlib figure
                                pfig.savefig(outp, format='png', dpi=300, bbox_inches='tight')
                                self._log_network(f"Saved {desc}: {outp}\n")
                                saved += 1
                            else:
                                # Plotly figure
                                pfig.write_image(outp, format='png')
                                self._log_network(f"Saved {desc}: {outp}\n")
                                saved += 1
                        except Exception as e1:
                            # Fallback to HTML if image export not available
                            try:
                                outp = os.path.join(pathways_networks_dir, fname.replace('.png', '.html'))
                                if hasattr(pfig, 'write_html'):
                                    # Plotly figure
                                    pfig.write_html(outp)
                                    self._log_network(f"Saved {desc} as HTML: {outp}\n")
                                    saved += 1
                                else:
                                    self._log_network(f"Cannot save {desc} (unsupported format)\n")
                            except Exception as e2:
                                self._log_network(f"Failed to save {desc}: {e2}\n")
                except Exception as e:
                    self._log_network(f"Plot generation failed: {e}\n")

            # Always export a filtered Excel workbook capturing selections and stats
            self._log_network("Starting Excel export...\n")
            try:
                import numpy as np
                selected_norm = {self._normalize_pathway_key(k) for k in selected_pathways.keys()}

                def _row_in_selected(row):
                    names = []
                    if isinstance(row.get('All_Pathways_Display'), list):
                        names = row.get('All_Pathways_Display')
                    elif isinstance(row.get('All_Pathways_Display'), str):
                        names = self._split_pathways(row.get('All_Pathways_Display'))
                    elif isinstance(row.get('All_Pathways'), list):
                        names = row.get('All_Pathways')
                    elif isinstance(row.get('All_Pathways'), str):
                        names = self._split_pathways(row.get('All_Pathways'))
                    cleaned = {self._normalize_pathway_key(n) for n in names if not self._should_exclude_pathway_name(str(n))}
                    return bool(cleaned.intersection(selected_norm))

                df_selected_mets = df_combined[df_combined.apply(_row_in_selected, axis=1)].copy()
                # Add -log10 pvalue for convenience
                with np.errstate(divide='ignore'):
                    pval_col = self.pvalue_col if self.pvalue_col in df_selected_mets.columns else 'pvalue'
                    df_selected_mets['neg_log10_p'] = -np.log10(df_selected_mets[pval_col].fillna(1.0))

                # Pathways sheet
                path_rows = []
                for pname, st in selected_pathways.items():
                    row = {
                        'Pathway': pname,
                        'Status': st.get('status'),
                        'z_score': st.get('z_score'),
                        'combined_pvalue': st.get('combined_pvalue'),
                        'n_metabolites': st.get('n_metabolites'),
                        'mean_log2fc': st.get('mean_log2fc'),
                    }
                    path_rows.append(row)
                df_paths = pd.DataFrame(path_rows)

                excel_path = os.path.join(pathways_networks_dir, 'selected_network_export.xlsx')
                component_label = self._component_label(plural=True)
                sheet_name = component_label

                with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
                    df_paths.to_excel(writer, sheet_name='Selected_Pathways', index=False)
                    (upstream_data if upstream_data is not None else pd.DataFrame()).to_excel(writer, sheet_name='Upstream', index=False)
                    (disease_data if disease_data is not None else pd.DataFrame()).to_excel(writer, sheet_name='Diseases', index=False)
                    df_selected_mets.to_excel(writer, sheet_name=sheet_name, index=False)
                self._log_network(f"✓ Saved filtered Excel: {excel_path}\n")
                self._log_network(f"  - Selected Pathways: {len(df_paths)} rows\n")
                self._log_network(f"  - Upstream: {len(upstream_data) if upstream_data is not None else 0} rows\n")
                self._log_network(f"  - Diseases: {len(disease_data) if disease_data is not None else 0} rows\n")
                self._log_network(f"  - {component_label}: {len(df_selected_mets)} rows\n")
            except Exception as e:
                self._log_network(f"Excel export failed: {e}\n")

            self.frame.after(0, lambda: messagebox.showinfo("Success", f"Network generated.\nSaved to: {pathways_networks_dir}"))
            
            # Auto-launch interactive viewer if enabled
            if INTERACTIVE_VIEWER_AVAILABLE and hasattr(self, 'auto_launch_interactive') and self.auto_launch_interactive.get():
                self._log_network(f"\n🌐 Auto-launching interactive viewer...\n")
                self.frame.after(500, self.launch_interactive_viewer)
            
            # Snapshot selections for interactive viewer fallback
            try:
                self._last_viewer_selection = {
                    'pathways': set(selected_pathways.keys()),
                    'upstream': set(selected_upstream_names),
                    'diseases': set(selected_disease_names)
                }
            except Exception:
                self._last_viewer_selection = {'pathways': set(), 'upstream': set(), 'diseases': set()}

            # NOTE: Removed auto-deselection - users can manually deselect if needed
            # self._reset_all_selections()
        except Exception as e:
            import traceback
            self._log_network(f"\n❌ ERROR: {e}\n{traceback.format_exc()}\n")
            err_msg = str(e)
            self.frame.after(0, lambda msg=err_msg: messagebox.showerror("Error", f"Network generation failed: {msg}"))

    def launch_interactive_viewer(self):
        """Launch the interactive network viewer in a web browser."""
        if not INTERACTIVE_VIEWER_AVAILABLE:
            messagebox.showerror("Not Available", 
                                 "Interactive viewer not available.\n\n"
                                 "Install with: pip install dash dash-cytoscape")
            return
        
        # Launch viewer directly - mode can be switched in browser
        self._launch_viewer_with_mode()
    
    def _launch_viewer_with_mode(self):
        """Actually launch the viewer with the selected mode."""
        # Collect pathway data from the filtered results
        pathway_data = getattr(self, 'pathway_filtered_pathways_data', None)
        metabolite_data = getattr(self, 'pathway_filtered_metabolites_data', None)
        
        if not pathway_data:
            messagebox.showwarning("No Data", 
                                   "No pathway data available.\n\n"
                                   "Please run Pathway Annotation first.")
            return
        
        snapshot = self._collect_selection_snapshot()
        selected_pathways = snapshot['pathways']
        selected_upstream = snapshot['upstream']
        selected_diseases = snapshot['diseases']

        if selected_pathways:
            self._last_viewer_selection = snapshot
        elif self._last_viewer_selection.get('pathways'):
            selected_pathways = set(self._last_viewer_selection.get('pathways', set()))
            selected_upstream = set(self._last_viewer_selection.get('upstream', set()))
            selected_diseases = set(self._last_viewer_selection.get('diseases', set()))
            self._log_network("ℹ️ No current checkboxes selected; reusing last network selection for viewer.\n")
        # If still no pathways selected, ask user if they want to use all
        if not selected_pathways:
            use_all = messagebox.askyesno(
                "No Pathways Selected",
                "No pathways are selected with checkboxes.\n\n"
                "Do you want to visualize ALL pathways?\n\n"
                "(Click No to cancel and select specific pathways)"
            )
            if not use_all:
                return
            # Use all pathways
            selected_pathways = set(k for k in pathway_data.keys() if not k.startswith('_'))
            self._last_viewer_selection = {
                'pathways': set(selected_pathways),
                'upstream': set(selected_upstream),
                'diseases': set(selected_diseases)
            }
        
        try:
            # Convert pathway dict to DataFrame for the viewer (ONLY selected pathways)
            rows = []
            for pathway_name, pathway_info in pathway_data.items():
                if pathway_name.startswith('_'):
                    continue
                
                # Filter to selected pathways only
                if pathway_name not in selected_pathways:
                    continue
                    
                # Extract pathway info
                if isinstance(pathway_info, dict):
                    hits = pathway_info.get('k_hits', pathway_info.get('count', pathway_info.get('n_metabolites', 1)))
                    pvalue = pathway_info.get('combined_pvalue', pathway_info.get('pvalue', pathway_info.get('p_value', 1.0)))
                    fdr = pathway_info.get('fdr', pathway_info.get('FDR', pvalue))
                    z_score = pathway_info.get('z_score', 0.0)
                    status = pathway_info.get('status', 'No Change')
                    metabolites = pathway_info.get('metabolites', [])
                    if isinstance(metabolites, list):
                        metabolites = "; ".join(metabolites)
                else:
                    hits = 1
                    pvalue = 1.0
                    fdr = 1.0
                    z_score = 0.0
                    status = 'No Change'
                    metabolites = ""
                
                rows.append({
                    'Pathway': pathway_name,
                    'Hits': hits,
                    'P-Value': pvalue,
                    'FDR': fdr,
                    'Z-Score': z_score,
                    'Status': status,
                    'Metabolites': metabolites,
                })
            
            if not rows:
                messagebox.showwarning("No Data", "No pathways to display.")
                return
            
            df_pathways = pd.DataFrame(rows)
            
            # Log the action
            self._log_network(f"\n🌐 Launching Interactive Network Viewer...\n")
            self._log_network(f"   • {len(df_pathways)} selected pathways\n")
            if selected_upstream:
                self._log_network(f"   • {len(selected_upstream)} selected upstream regulators\n")
            if selected_diseases:
                self._log_network(f"   • {len(selected_diseases)} selected diseases\n")
            
            # Prepare metabolite data for viewer
            viewer_metabolite_data = None
            print(f"\n🔍 DEBUG: Preparing metabolite data for viewer")
            print(f"   metabolite_data is None: {metabolite_data is None}")
            if metabolite_data is not None:
                print(f"   metabolite_data type: {type(metabolite_data)}")
                print(f"   metabolite_data shape: {metabolite_data.shape if hasattr(metabolite_data, 'shape') else 'N/A'}")
                print(f"   metabolite_data empty: {metabolite_data.empty if hasattr(metabolite_data, 'empty') else 'N/A'}")
                if hasattr(metabolite_data, 'columns'):
                    print(f"   metabolite_data columns: {metabolite_data.columns.tolist()}")
            
            if metabolite_data is not None and not metabolite_data.empty:
                # Ensure metabolite data has required columns: Name, log2FC, pvalue
                viewer_metabolite_data = metabolite_data.copy()
                # Check and standardize column names
                if 'Name' not in viewer_metabolite_data.columns:
                    name_candidates = ['name', 'Metabolite', 'metabolite', 'LipidID', 'Lipid_ID']
                    for col in name_candidates:
                        if col in viewer_metabolite_data.columns:
                            viewer_metabolite_data['Name'] = viewer_metabolite_data[col]
                            break
                if 'log2FC' not in viewer_metabolite_data.columns:
                    fc_candidates = ['log2_FC', 'Log2FC', 'log2FoldChange', 'FC']
                    for col in fc_candidates:
                        if col in viewer_metabolite_data.columns:
                            viewer_metabolite_data['log2FC'] = viewer_metabolite_data[col]
                            break
                if 'pvalue' not in viewer_metabolite_data.columns:
                    pval_candidates = ['p_value', 'P-Value', 'p-value', 'adj_p']
                    for col in pval_candidates:
                        if col in viewer_metabolite_data.columns:
                            viewer_metabolite_data['pvalue'] = viewer_metabolite_data[col]
                            break
                
                print(f"   ✅ After standardization:")
                print(f"      viewer_metabolite_data shape: {viewer_metabolite_data.shape}")
                print(f"      viewer_metabolite_data columns: {viewer_metabolite_data.columns.tolist()}")
                print(f"      Has 'Name': {'Name' in viewer_metabolite_data.columns}")
                print(f"      Has 'log2FC': {'log2FC' in viewer_metabolite_data.columns}")
                print(f"      Has 'pvalue': {'pvalue' in viewer_metabolite_data.columns}")
                
                # CRITICAL: Filter to only metabolites connected to selected pathways/upstream/diseases
                # This matches the GraphML export behavior
                metabolites_to_include = set()
                
                # Get metabolites from selected pathways
                for pathway_name in selected_pathways:
                    pathway_info = pathway_data.get(pathway_name, {})
                    if isinstance(pathway_info, dict):
                        mets = pathway_info.get('metabolites', [])
                        metabolites_to_include.update(mets)
                
                print(f"\n   🔍 Filtering metabolites to only those connected to selections:")
                print(f"      Metabolites from {len(selected_pathways)} pathways: {len(metabolites_to_include)}")
                
                # Get metabolites from selected upstream (will be added in next section)
                upstream_metabolites = set()
                disease_metabolites = set()
                
                # Store these for later addition
                initial_pathway_metabolites = len(metabolites_to_include)
                
                self._log_network(f"   • {len(viewer_metabolite_data)} total metabolites in data\n")
            else:
                print(f"   ⚠️ No metabolite data available for viewer!")
                self._log_network(f"   ⚠️ No metabolite data available\n")
                metabolites_to_include = set()
                upstream_metabolites = set()
                disease_metabolites = set()
                initial_pathway_metabolites = 0
            
            # Prepare upstream and disease data
            upstream_data = None
            disease_data = None
            
            if True:  # Always prepare data for both modes
                upstream_data = None
                disease_data = None
                
                if selected_upstream:
                    # Use upstream_to_metabolites mapping which is populated during _populate_upstream_tree
                    upstream_list = []
                    for upstream_name in selected_upstream:
                        # Get metabolites from upstream_to_metabolites mapping
                        mets = self.upstream_to_metabolites.get(upstream_name, set())
                        if isinstance(mets, set):
                            mets = list(mets)
                        upstream_metabolites.update(mets)  # Add to filter set
                        
                        # Get type info from original_upstream_data if available
                        upstream_type = 'enzyme'  # default
                        if hasattr(self, 'original_upstream_data') and self.original_upstream_data:
                            info = self.original_upstream_data.get(upstream_name, {})
                            if isinstance(info, dict):
                                upstream_type = info.get('type', 'enzyme')
                        
                        upstream_list.append({
                            'name': upstream_name,
                            'type': upstream_type,
                            'metabolites': '; '.join(mets) if mets else '',
                        })
                    if upstream_list:
                        upstream_data = pd.DataFrame(upstream_list)
                    print(f"      Metabolites from {len(selected_upstream)} upstream: {len(upstream_metabolites)}")
                
                if selected_diseases:
                    # Use disease_to_metabolites mapping which is populated during _populate_disease_tree
                    disease_list = []
                    for disease_name in selected_diseases:
                        # Get metabolites from disease_to_metabolites mapping
                        mets = self.disease_to_metabolites.get(disease_name, set())
                        if isinstance(mets, set):
                            mets = list(mets)
                        disease_metabolites.update(mets)  # Add to filter set
                        
                        disease_list.append({
                            'name': disease_name,
                            'metabolites': '; '.join(mets) if mets else '',
                        })
                    if disease_list:
                        disease_data = pd.DataFrame(disease_list)
                    print(f"      Metabolites from {len(selected_diseases)} diseases: {len(disease_metabolites)}")
                
                # Combine all metabolite sources
                metabolites_to_include.update(upstream_metabolites)
                metabolites_to_include.update(disease_metabolites)
                
                print(f"      TOTAL metabolites to include: {len(metabolites_to_include)}")
                
                # Filter viewer_metabolite_data to only include connected metabolites
                if viewer_metabolite_data is not None and not viewer_metabolite_data.empty:
                    if 'Name' in viewer_metabolite_data.columns and metabolites_to_include:
                        before_count = len(viewer_metabolite_data)
                        viewer_metabolite_data = viewer_metabolite_data[
                            viewer_metabolite_data['Name'].isin(metabolites_to_include)
                        ].copy()
                        after_count = len(viewer_metabolite_data)
                        print(f"      ✅ Filtered metabolites: {before_count} → {after_count}")
                        self._log_network(f"   • {after_count} connected metabolites (filtered from {before_count})\n")
                    else:
                        print(f"      ⚠️ Could not filter metabolites (no Name column or no metabolites to include)")
                        self._log_network(f"   • ⚠️ Metabolite filtering skipped\n")
                
            # Create viewer with both modes available (user can switch in browser)
            # Pass the FULL pathway_data dict (not DataFrame) as pathway_stats for z-scores/status
            # Get focused mode setting from GUI
            is_focused_mode = getattr(self, 'focused_mode', tk.BooleanVar(value=True)).get()
            
            print(f"\n🚀 Creating PathwayNetworkViewer:")
            print(f"   df_pathways shape: {df_pathways.shape}")
            print(f"   viewer_metabolite_data: {viewer_metabolite_data.shape if viewer_metabolite_data is not None else 'None'}")
            print(f"   upstream_data: {upstream_data.shape if upstream_data is not None else 'None'}")
            print(f"   disease_data: {disease_data.shape if disease_data is not None else 'None'}")
            print(f"   Focused mode: {is_focused_mode}")
            print(f"   Starting mode: 'full_network' (default)\n")
            
            viewer = PathwayNetworkViewer(
                pathway_data=df_pathways,
                metabolite_data=viewer_metabolite_data if viewer_metabolite_data is not None else pd.DataFrame(),
                upstream_data=upstream_data if upstream_data is not None else pd.DataFrame(),
                disease_data=disease_data if disease_data is not None else pd.DataFrame(),
                pathway_stats=pathway_data,  # Pass full stats dict (not DataFrame) for z-scores
                title="Interactive Pathway Network",
                mode='full_network',  # Start in full network mode (default)
                show_edges=True,
                # Use the lowest edge threshold so ANY shared metabolite
                # between two pathways will create a similarity edge when
                # "Pathway Similarity" mode is selected in the viewer.
                min_edge_weight=0.0,
                focused_mode=is_focused_mode  # Pass focused mode setting from GUI
            )
            
            self._log_network(f"   • Opening in browser...\n\n")
            
            # Run threaded so it doesn't block the GUI
            port = viewer.run_threaded(open_browser=True)
            
            if port:
                self._log_network(f"✅ Interactive viewer running at: http://127.0.0.1:{port}\n")
                self._log_network(f"   Browser should open automatically.\n")
                self._log_network(f"   Toggle between Pathway Similarity and Full Network modes in the browser\n")
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            messagebox.showerror("Error", f"Failed to launch interactive viewer:\n{str(e)}")

    def run_pathway_enrichment(self):
        """
        Run pathway enrichment analysis (Fisher ORA or IWPA) on verified metabolites data.
        This calculates pathway statistics and populates pathway_filtered_pathways_data.
        Called after column verification when enrichment data doesn't already exist.
        """
        import threading
        
        def _enrichment_worker():
            try:
                self._log_network(f"{'='*70}\n")
                self._log_network(f"🧬 PATHWAY ENRICHMENT ANALYSIS\n")
                self._log_network(f"{'='*70}\n\n")
                
                # Get verified metabolites data
                self._log_network(f"🔍 Checking verified metabolites data...\n")
                # CRITICAL: Always use original uploaded raw data as source for enrichment
                # This ensures each enrichment run starts fresh with unprocessed data
                df_raw = getattr(self, 'original_uploaded_raw_data', None)
                if df_raw is None:
                    # Fallback to verified_network_dataframe_raw if original not stored (backwards compatibility)
                    df_raw = getattr(self, 'verified_network_dataframe_raw', None)
                df_display = getattr(self, 'verified_network_dataframe', None)
                
                # Verify we're using original raw data, not processed results
                if df_raw is not None:
                    self._log_network(f"✓ Using original_uploaded_raw_data ({len(df_raw)} rows, {len(df_raw.columns)} columns)\n")
                    if hasattr(self, 'pathway_filtered_metabolites_data') and self.pathway_filtered_metabolites_data is not None:
                        self._log_network(f"   (pathway_filtered_metabolites_data also exists with {len(self.pathway_filtered_metabolites_data)} rows - NOT used for enrichment)\n")
                
                if df_raw is None:
                    self._log_network(f"❌ ERROR: No metabolites dataframe available (raw input missing)\n")
                    self.root.after(0, lambda: messagebox.showerror(
                        "Error",
                        "No metabolites data available for enrichment analysis.\n"
                        "Please upload and verify data first."
                    ))
                    return
                if len(df_raw) == 0:
                    self._log_network(f"❌ ERROR: original raw data is empty\n")
                    self.root.after(0, lambda: messagebox.showerror(
                        "Error",
                        "No metabolites data available for enrichment analysis."
                    ))
                    return
                
                # Use original raw data for enrichment calculation
                df = df_raw.copy()
                # Use display dataframe for row count reporting only
                display_df = df_display.copy() if df_display is not None else df.copy()

                self._log_network(f"✓ Data shape (pre-compression): {df.shape}\n")
                if len(display_df) != len(df):
                    self._log_network(f"✓ Display nodes after class compression: {display_df.shape}\n")
                
                # Apply column renaming using stored rename map (maps original column names to standard names)
                if hasattr(self, 'verified_network_rename_map') and self.verified_network_rename_map:
                    # Only apply renames for columns that exist in current df
                    applicable_renames = {old: new for old, new in self.verified_network_rename_map.items() if old in df.columns}
                    if applicable_renames:
                        df.rename(columns=applicable_renames, inplace=True)
                        self._log_network(f"✓ Applied column renaming: {applicable_renames}\n")
                
                self._log_network(f"✓ Columns available: {list(df.columns)}\n\n")
                
                # Ensure required columns exist
                self._log_network(f"🔍 Checking required columns...\n")
                missing_cols = []
                if 'Name' not in df.columns:
                    missing_cols.append('Name')
                if 'pvalue' not in df.columns:
                    missing_cols.append('pvalue')
                if 'log2FC' not in df.columns:
                    missing_cols.append('log2FC')
                
                if missing_cols:
                    self._log_network(f"❌ ERROR: Missing columns: {missing_cols}\n")
                    self.root.after(0, lambda: messagebox.showerror(
                        "Error",
                        f"Required columns missing: {', '.join(missing_cols)}\n"
                        "Please verify columns first."
                    ))
                    return
                
                self._log_network(f"✓ Required columns present: Name, pvalue, log2FC\n\n")
                
                # Check if data has All_Pathways column (required for enrichment)
                self._log_network(f"🔍 Checking for pathway annotations...\n")
                if 'All_Pathways' not in df.columns:
                    # Check if All_Pathways_Display exists as alternative
                    if 'All_Pathways_Display' in df.columns:
                        self._log_network(f"⚠️ All_Pathways not found, but All_Pathways_Display exists\n")
                        self._log_network(f"   Converting separator from semicolons to pipes for enrichment...\n")
                        # Convert semicolons to pipes and DEDUPLICATE (enrichment function expects pipe-separated)
                        def deduplicate_pathways(pathway_str):
                            """Remove duplicate pathways (case-insensitive, comma-insensitive) from pipe/semicolon-separated string"""
                            if not pathway_str or pd.isna(pathway_str):
                                return ''
                            # Split by both ; and | to handle any format
                            pathways = [p.strip() for p in str(pathway_str).replace(';', '|').split('|') if p.strip()]
                            # Deduplicate using normalized comparison (remove commas, case-insensitive)
                            # Keep the FIRST occurrence with most characters (likely most complete name)
                            seen = {}
                            unique = []
                            for p in pathways:
                                # Normalize for comparison: remove commas and lowercase
                                p_normalized = p.lower().replace(',', '')
                                if p_normalized not in seen:
                                    seen[p_normalized] = p
                                    unique.append(p)
                                else:
                                    # If we already have this pathway, keep the longer version (more complete)
                                    if len(p) > len(seen[p_normalized]):
                                        # Replace with longer version
                                        idx = unique.index(seen[p_normalized])
                                        unique[idx] = p
                                        seen[p_normalized] = p
                            return '|'.join(unique)
                        
                        df['All_Pathways'] = df['All_Pathways_Display'].apply(deduplicate_pathways)
                        self._log_network(f"✓ All_Pathways column created with pipe separators (duplicates removed)\n\n")
                    else:
                        self._log_network(f"❌ ERROR: Neither All_Pathways nor All_Pathways_Display found\n")
                        self._log_network(f"   Available columns: {list(df.columns)}\n")
                        self.root.after(0, lambda: messagebox.showerror(
                            "Error",
                            "All_Pathways column not found.\n"
                            "Enrichment requires pathway-annotated data.\n\n"
                            "Please run Pathway Annotation first."
                        ))
                        return
                else:
                    # Ensure existing All_Pathways uses pipes and DEDUPLICATE
                    def deduplicate_pathways(pathway_str):
                        """Remove duplicate pathways (case-insensitive, comma-insensitive) from pipe/semicolon-separated string"""
                        if not pathway_str or pd.isna(pathway_str):
                            return ''
                        # Split by both ; and | to handle any format
                        pathways = [p.strip() for p in str(pathway_str).replace(';', '|').split('|') if p.strip()]
                        # Deduplicate using normalized comparison (remove commas, case-insensitive)
                        # Keep the FIRST occurrence with most characters (likely most complete name)
                        seen = {}
                        unique = []
                        for p in pathways:
                            # Normalize for comparison: remove commas and lowercase
                            p_normalized = p.lower().replace(',', '')
                            if p_normalized not in seen:
                                seen[p_normalized] = p
                                unique.append(p)
                            else:
                                # If we already have this pathway, keep the longer version (more complete)
                                if len(p) > len(seen[p_normalized]):
                                    # Replace with longer version
                                    idx = unique.index(seen[p_normalized])
                                    unique[idx] = p
                                    seen[p_normalized] = p
                        return '|'.join(unique)
                    
                    df['All_Pathways'] = df['All_Pathways'].apply(deduplicate_pathways)
                    self._log_network(f"✓ All_Pathways column found (duplicates removed, separators normalized to pipes)\n\n")
                
                # Get enrichment method and parameters
                self._log_network(f"\n{'='*60}\n")
                self._log_network(f"🔬 STARTING NEW PATHWAY ENRICHMENT CALCULATION\n")
                self._log_network(f"{'='*60}\n\n")
                
                method = self.pathway_analysis_method.get() if hasattr(self, 'pathway_analysis_method') else "Fisher ORA"
                apply_fdr = self.pathway_apply_fdr.get() if hasattr(self, 'pathway_apply_fdr') else True
                
                # Get parameters - prioritize pathway_min_metabolites from settings panel
                try:
                    min_metabolites = int(self.pathway_min_metabolites.get()) if hasattr(self, 'pathway_min_metabolites') else int(self.network_min_metabolites.get())
                    p_threshold = float(self.metabolite_pvalue_filter.get()) if hasattr(self, 'metabolite_pvalue_filter') else 0.05
                    log2fc_threshold = float(self.metabolite_log2fc_filter.get()) if hasattr(self, 'metabolite_log2fc_filter') else 0.0
                    pathway_p_threshold = float(self.pathway_pvalue_threshold.get()) if hasattr(self, 'pathway_pvalue_threshold') else 0.05
                    # CRITICAL FIX: Read from correct variable name (pathway_z_threshold, not pathway_zscore_threshold)
                    z_threshold = float(self.pathway_z_threshold.get()) if hasattr(self, 'pathway_z_threshold') else 1.0
                except (ValueError, AttributeError) as e:
                    self._log_network(f"⚠️ Error extracting parameters (using defaults): {e}\n")
                    min_metabolites = 2
                    p_threshold = 0.05
                    log2fc_threshold = 0.0
                    pathway_p_threshold = 0.05
                    z_threshold = 1.0
                
                # Get additional settings
                try:
                    species = self.pathway_organism_var.get() if hasattr(self, 'pathway_organism_var') else 'Homo sapiens'
                    weight_mode = self.iwpa_weight_mode.get() if hasattr(self, 'iwpa_weight_mode') else 'signed_p'
                    max_pathways = self.pathway_max_total.get() if hasattr(self, 'pathway_max_total') else 25
                    limit_enabled = self.pathway_limit_pathways.get() if hasattr(self, 'pathway_limit_pathways') else True
                except:
                    species = 'Homo sapiens'
                    weight_mode = 'signed_p'
                    max_pathways = 25
                    limit_enabled = True
                
                self._log_network(f"📋 ENRICHMENT METHOD: {method}\n")
                self._log_network(f"{'─'*60}\n")
                self._log_network(f"📊 ALL CURRENT SETTINGS:\n")
                self._log_network(f"\n   🔧 Analysis Parameters:\n")
                self._log_network(f"      • Method: {method}\n")
                self._log_network(f"      • Species/Organism: {species}\n")
                self._log_network(f"      • Apply FDR correction (BH): {apply_fdr}\n")
                if method == "IWPA":
                    self._log_network(f"      • IWPA Weight Mode: {weight_mode}\n")
                self._log_network(f"\n   🎯 Filter Thresholds:\n")
                self._log_network(f"      • Min metabolites per pathway: {min_metabolites}\n")
                if method == "Fisher ORA":
                    self._log_network(f"      • Metabolite p-value filter: ≤ {p_threshold} (defines significant metabolites)\n")
                    self._log_network(f"      • Metabolite |log2FC| filter: ≥ {log2fc_threshold} (absolute fold change threshold)\n")
                else:  # IWPA
                    self._log_network(f"      • Metabolite p-value filter: ≤ {p_threshold} (same as Fisher ORA for consistency)\n")
                    self._log_network(f"      • ✓ IWPA uses CONTINUOUS WEIGHTS within significant set (not binary)\n")
                    self._log_network(f"      • Metabolite |log2FC| filter: Not used (continuous weights in IWPA)\n")
                self._log_network(f"      • Pathway p-value threshold: ≤ {pathway_p_threshold} (for DISPLAY filtering only)\n")
                self._log_network(f"      • Z-score threshold (±): {z_threshold}\n")
                self._log_network(f"      • Max pathways displayed: {max_pathways if limit_enabled else 'Unlimited'}\n")
                self._log_network(f"\n   📈 Input Data:\n")
                self._log_network(f"      • Total metabolites loaded (pre-compression): {len(df)}\n")
                if len(display_df) != len(df):
                    self._log_network(f"      • Display nodes after class compression: {len(display_df)}\n")
                # Count significant metabolites (applying both p-value and log2FC filters)
                sig_count = 0
                pvals_numeric = None
                sig_mask = None
                if 'pvalue' in df.columns:
                    pvals_numeric = pd.to_numeric(df['pvalue'], errors='coerce')
                    pval_mask = pvals_numeric <= p_threshold
                    sig_mask = pval_mask
                    if log2fc_threshold > 0 and 'log2FC' in df.columns:
                        fc_series = pd.to_numeric(df['log2FC'], errors='coerce').abs()
                        fc_mask = fc_series >= log2fc_threshold
                        sig_mask = pval_mask & fc_mask
                        sig_count = int(sig_mask.fillna(False).sum())
                        self._log_network(f"      • Significant metabolites (p≤{p_threshold} AND |log2FC|≥{log2fc_threshold}): {sig_count}\n")
                    else:
                        sig_count = int(pval_mask.fillna(False).sum())
                        self._log_network(f"      • Significant metabolites (p≤{p_threshold}): {sig_count}\n")
                else:
                    self._log_network("      • ⚠️ 'pvalue' column not found; skipping significance count\n")

                significant_names = []
                if sig_count > 0 and sig_mask is not None and 'Name' in df.columns:
                    sig_mask_bool = sig_mask.fillna(False)
                    significant_names = (
                        df.loc[sig_mask_bool & df['Name'].notna(), 'Name']
                        .astype(str)
                        .str.strip()
                        .tolist()
                    )
                significant_display_nodes = set()
                if significant_names and self.class_compression_active:
                    for name in significant_names:
                        mapped = self.metabolite_class_lookup.get(name, name)
                        if mapped:
                            significant_display_nodes.add(mapped)
                if self.class_compression_active:
                    self._log_network(
                        f"      • Significant class nodes after compression: {len(significant_display_nodes)}\n"
                    )

                if sig_count == 0:
                    smallest_p = None
                    if pvals_numeric is not None:
                        try:
                            smallest_p = float(pvals_numeric.min(skipna=True))
                        except Exception:
                            smallest_p = None
                    if smallest_p is not None and not np.isnan(smallest_p):
                        self._log_network(
                            f"      ⚠️ No metabolites met the current significance threshold. "
                            f"Smallest observed p-value: {smallest_p:.2e}.\n"
                        )
                        self._log_network(
                            "      💡 Tip: Increase the 'Metabolite p-value filter' setting or disable class compression to include more metabolites.\n"
                        )
                    else:
                        self._log_network(
                            "      ⚠️ No metabolites met the current significance threshold. Adjust the 'Metabolite p-value filter' to continue.\n"
                        )
                self._log_network(f"\n{'─'*60}\n\n")
                
                # Run enrichment calculation based on selected method
                self._log_network(f"🔬 Calculating pathway statistics using {method}...\n")
                self._log_network(f"   This may take a few moments...\n\n")
                
                # Progress callback
                def progress_update(percent, message=""):
                    self.root.after(0, lambda p=percent, m=message: self._log_network(f"   [{p:3d}%] {m}\n"))
                
                self._log_network(f"🚀 Starting {method} calculation...\n")
                self._log_network(f"\n📋 Function Call Parameters:\n")
                self._log_network(f"   • min_metabolites={min_metabolites}\n")
                if method == "Fisher ORA":
                    self._log_network(f"   • metabolite_pvalue_threshold={p_threshold} (defines which metabolites are 'significant')\n")
                    self._log_network(f"   • metabolite_log2fc_threshold={log2fc_threshold} (|log2FC| >= threshold)\n")
                else:  # IWPA
                    self._log_network(f"   • metabolite_pvalue_threshold={p_threshold} (same as Fisher ORA)\n")
                    self._log_network(f"   • ✓ IWPA uses CONTINUOUS WEIGHTS (not binary hit/no-hit)\n")
                self._log_network(f"   • z_threshold={z_threshold} ⚠️ CRITICAL FOR STATUS\n")
                self._log_network(f"   • fdr_method={'fdr_bh' if apply_fdr else None}\n")
                self._log_network(f"\n   📊 Display Filtering (applied AFTER calculation):\n")
                self._log_network(f"   • pathway_p_threshold={pathway_p_threshold} (only show pathways with p≤{pathway_p_threshold})\n")
                self._log_network(f"   • max_pathways={max_pathways if limit_enabled else None}\n")
                if method == "Fisher ORA":
                    self._log_network(f"   • species={getattr(self, 'pathway_organism_var', tk.StringVar(value='Homo sapiens')).get()}\n")
                else:
                    self._log_network(f"   • weight_mode={weight_mode}\n")
                self._log_network(f"\n")
                
                try:
                    if method == "Fisher ORA":
                        # Import Fisher ORA function
                        from main_script.fisher_ora_pathway_analysis import calculate_pathway_statistics_fisher_ora
                        
                        # Get species selection
                        species = getattr(self, 'pathway_organism_var', tk.StringVar(value='Homo sapiens')).get()
                        
                        # CRITICAL: Fisher ORA pvalue_threshold parameter is for METABOLITE significance (0.05)
                        # pathway_p_threshold is only used for DISPLAY filtering AFTER calculation
                        # Use max_pathways from GUI if limit is enabled
                        max_pathways_param = max_pathways if limit_enabled else None
                        
                        self._log_network(f"🧪 Calling Fisher ORA with:\n")
                        self._log_network(f"   metabolite_p_threshold={p_threshold} (defines significant metabolites)\n")
                        self._log_network(f"   metabolite_log2fc_threshold={log2fc_threshold} (|log2FC| >= threshold)\n")
                        self._log_network(f"   z_threshold={z_threshold}\n")
                        self._log_network(f"   max_pathways={max_pathways_param}\n")
                        self._log_network(f"   pathway_p_threshold={pathway_p_threshold} will be used for display filtering AFTER calculation\n\n")
                        
                        pathways_data = calculate_pathway_statistics_fisher_ora(
                            df,
                            min_metabolites=min_metabolites,
                            pvalue_threshold=p_threshold,  # CRITICAL: Use metabolite filter (0.05), NOT pathway threshold
                            log2fc_threshold=log2fc_threshold,  # Absolute log2FC threshold for significance
                            z_threshold=z_threshold,
                            fdr_method='fdr_bh' if apply_fdr else None,
                            max_total_pathways=None,  # Don't limit in calculation, filter for display later
                            species=species,
                            auto_adjust_pvalue_threshold=False,
                            progress_callback=progress_update
                        )
                    elif method == "IWPA":
                        # Import IWPA function
                        from main_script.fisher_ora_pathway_analysis import calculate_pathway_statistics_iwpa
                        
                        # CRITICAL: Use max_pathways from GUI if limit is enabled
                        max_pathways_param = max_pathways if limit_enabled else None
                        
                        # CRITICAL: IWPA uses the same p-value threshold as Fisher ORA for consistency
                        # But applies continuous weights (not binary hit/no-hit like Fisher ORA)
                        self._log_network(f"🧪 Calling IWPA with:\n")
                        self._log_network(f"   metabolite_p_threshold={p_threshold} (same as Fisher ORA for consistency)\n")
                        self._log_network(f"   ✓ IWPA uses CONTINUOUS WEIGHTS (not binary like Fisher ORA)\n")
                        self._log_network(f"   z_threshold={z_threshold}, weight_mode={weight_mode}, max_pathways={max_pathways_param}\n\n")
                        
                        pathways_data = calculate_pathway_statistics_iwpa(
                            df,
                            min_metabolites=min_metabolites,
                            pvalue_threshold=p_threshold,  # CRITICAL: Pass GUI threshold
                            fdr_method='fdr_bh' if apply_fdr else None,
                            weight_mode=weight_mode,
                            z_threshold=z_threshold,
                            max_total_pathways=max_pathways_param,  # FIXED: Use GUI value
                            progress_callback=progress_update
                        )
                    else:
                        raise ValueError(f"Unknown enrichment method: {method}")
                    
                    self._log_network(f"✓ {method} calculation completed\n\n")
                except Exception as calc_error:
                    self._log_network(f"❌ ERROR in {method} calculation: {calc_error}\n")
                    import traceback
                    self._log_network(f"{traceback.format_exc()}\n")
                    raise
                
                # CRITICAL: Normalize field names for display compatibility
                # Fisher ORA returns: fisher_pvalue, fdr_qvalue, direction_zscore, direction
                # IWPA returns: pvalue, adjusted_pvalue, Z, status
                # Network tab display expects: combined_pvalue, z_score, status
                self._log_network(f"🔧 Normalizing field names for display...\n")
                if pathways_data:
                    for pathway_name, stats in pathways_data.items():
                        if isinstance(stats, dict):
                            # Map direction_zscore → z_score (Fisher ORA)
                            if 'direction_zscore' in stats and 'z_score' not in stats:
                                stats['z_score'] = stats['direction_zscore']
                            
                            # Map Z → z_score (IWPA)
                            if 'Z' in stats and 'z_score' not in stats:
                                stats['z_score'] = stats['Z']
                            
                            # Map direction → status (Fisher ORA)
                            if 'direction' in stats and 'status' not in stats:
                                stats['status'] = stats['direction']
                            
                            # Map fisher_pvalue/fdr_qvalue → combined_pvalue (Fisher ORA)
                            if apply_fdr and 'fdr_qvalue' in stats:
                                stats['combined_pvalue'] = stats['fdr_qvalue']
                                stats['adjusted_pvalue'] = stats['fdr_qvalue']
                            elif 'fisher_pvalue' in stats:
                                stats['combined_pvalue'] = stats['fisher_pvalue']
                                stats['adjusted_pvalue'] = stats['fisher_pvalue']
                            # Map pvalue/adjusted_pvalue → combined_pvalue (IWPA)
                            elif apply_fdr and 'adjusted_pvalue' in stats:
                                stats['combined_pvalue'] = stats['adjusted_pvalue']
                            elif 'pvalue' in stats:
                                stats['combined_pvalue'] = stats['pvalue']
                                if 'adjusted_pvalue' not in stats:
                                    stats['adjusted_pvalue'] = stats['pvalue']
                    
                    self._log_network(f"✓ Field normalization complete\n\n")

                if pathways_data and self.class_compression_active:
                    self._log_network("🔁 Mapping significant metabolites to class nodes for display...\n")
                    # DEBUG: Log counts before and after remapping
                    sample_pathway = list(pathways_data.keys())[0] if pathways_data else None
                    if sample_pathway and isinstance(pathways_data[sample_pathway], dict):
                        before_k = pathways_data[sample_pathway].get('k_hits', 0)
                        before_n = pathways_data[sample_pathway].get('n_metabolites', 0)
                        self._log_network(f"   [DEBUG BEFORE] {sample_pathway}: k_hits={before_k}, n_metabolites={before_n}\n")
                    
                    self._map_pathway_hits_to_display_nodes(pathways_data)
                    
                    if sample_pathway and isinstance(pathways_data[sample_pathway], dict):
                        after_k = pathways_data[sample_pathway].get('k_hits', 0)
                        after_n = pathways_data[sample_pathway].get('n_metabolites', 0)
                        self._log_network(f"   [DEBUG AFTER] {sample_pathway}: k_hits={after_k}, n_metabolites={after_n}\n")
                    self._log_network("\n")
                
                # CRITICAL: Consolidate amino acid-specific pathways into single pathways
                # All "tRNA Charging: <amino acid>" → "tRNA Charging"
                # All "Protein Synthesis: <amino acid>" → "Protein Synthesis"
                self._log_network(f"🔄 Consolidating amino acid-specific pathways...\n")
                consolidated_pathways = {}
                
                # Track tRNA Charging pathways
                trna_metabolites = set()
                trna_stats_list = []
                trna_count = 0
                
                # Track Protein Synthesis pathways
                protein_metabolites = set()
                protein_stats_list = []
                protein_count = 0
                
                # Track Sphingolipid Metabolism pathways
                sphingolipid_metabolites = set()
                sphingolipid_stats_list = []
                sphingolipid_count = 0
                
                for pathway_name, stats in pathways_data.items():
                    if not isinstance(stats, dict):
                        consolidated_pathways[pathway_name] = stats
                        continue
                    
                    pathway_lower = pathway_name.lower()
                    
                    # Check if this is a tRNA Charging pathway
                    if 'trna charging:' in pathway_lower:
                        trna_count += 1
                        # Collect metabolites from this pathway
                        mets = stats.get('metabolites', []) or stats.get('Metabolites', [])
                        if isinstance(mets, (list, set)):
                            trna_metabolites.update(mets)
                        # Store stats for aggregation
                        trna_stats_list.append(stats)
                    # Check if this is a Protein Synthesis pathway
                    elif 'protein synthesis:' in pathway_lower:
                        protein_count += 1
                        # Collect metabolites from this pathway
                        mets = stats.get('metabolites', []) or stats.get('Metabolites', [])
                        if isinstance(mets, (list, set)):
                            protein_metabolites.update(mets)
                        # Store stats for aggregation
                        protein_stats_list.append(stats)
                    # Check if this is a Sphingolipid Metabolism variant (keep only base version)
                    elif pathway_lower.startswith('sphingolipid metabolism'):
                        sphingolipid_count += 1
                        # Collect metabolites from this pathway
                        mets = stats.get('metabolites', []) or stats.get('Metabolites', [])
                        if isinstance(mets, (list, set)):
                            sphingolipid_metabolites.update(mets)
                        # Store stats for aggregation
                        sphingolipid_stats_list.append(stats)
                    else:
                        # Keep other pathways as-is
                        consolidated_pathways[pathway_name] = stats
                
                # Merge tRNA Charging pathways if found
                if trna_count > 0:
                    self._log_network(f"   • Found {trna_count} tRNA Charging pathways, consolidating...\n")
                    
                    # Create consolidated tRNA Charging pathway
                    consolidated_name = "tRNA Charging"
                    
                    # Aggregate statistics: use the one with lowest p-value as base
                    best_stats = min(trna_stats_list, key=lambda s: s.get('combined_pvalue', 1.0))
                    merged_stats = dict(best_stats)
                    
                    # Update metabolite list to include all unique metabolites
                    merged_stats['metabolites'] = list(trna_metabolites)
                    merged_stats['Metabolites'] = list(trna_metabolites)
                    
                    # Update counts
                    merged_stats['k_hits'] = len(trna_metabolites)
                    merged_stats['n_metabolites'] = len(trna_metabolites)
                    
                    # Average z-scores if multiple
                    if len(trna_stats_list) > 1:
                        z_scores = [s.get('z_score', 0.0) for s in trna_stats_list]
                        merged_stats['z_score'] = sum(z_scores) / len(z_scores)
                    
                    consolidated_pathways[consolidated_name] = merged_stats
                    self._log_network(f"   • Merged into '{consolidated_name}' with {len(trna_metabolites)} unique metabolites\n")
                
                # Merge Protein Synthesis pathways if found
                if protein_count > 0:
                    self._log_network(f"   • Found {protein_count} Protein Synthesis pathways, consolidating...\n")
                    
                    # Create consolidated Protein Synthesis pathway
                    consolidated_name = "Protein Synthesis"
                    
                    # Aggregate statistics: use the one with lowest p-value as base
                    best_stats = min(protein_stats_list, key=lambda s: s.get('combined_pvalue', 1.0))
                    merged_stats = dict(best_stats)
                    
                    # Update metabolite list to include all unique metabolites
                    merged_stats['metabolites'] = list(protein_metabolites)
                    merged_stats['Metabolites'] = list(protein_metabolites)
                    
                    # Update counts
                    merged_stats['k_hits'] = len(protein_metabolites)
                    merged_stats['n_metabolites'] = len(protein_metabolites)
                    
                    # Average z-scores if multiple
                    if len(protein_stats_list) > 1:
                        z_scores = [s.get('z_score', 0.0) for s in protein_stats_list]
                        merged_stats['z_score'] = sum(z_scores) / len(z_scores)
                    
                    consolidated_pathways[consolidated_name] = merged_stats
                    self._log_network(f"   • Merged into '{consolidated_name}' with {len(protein_metabolites)} unique metabolites\n")
                
                # Merge Sphingolipid Metabolism pathways if found
                if sphingolipid_count > 0:
                    self._log_network(f"   • Found {sphingolipid_count} Sphingolipid Metabolism pathways, consolidating...\n")
                    
                    # Create consolidated Sphingolipid Metabolism pathway
                    consolidated_name = "Sphingolipid Metabolism"
                    
                    # Aggregate statistics: use the one with lowest p-value as base
                    best_stats = min(sphingolipid_stats_list, key=lambda s: s.get('combined_pvalue', 1.0))
                    merged_stats = dict(best_stats)
                    
                    # Update metabolite list to include all unique metabolites
                    merged_stats['metabolites'] = list(sphingolipid_metabolites)
                    merged_stats['Metabolites'] = list(sphingolipid_metabolites)
                    
                    # Update counts
                    merged_stats['k_hits'] = len(sphingolipid_metabolites)
                    merged_stats['n_metabolites'] = len(sphingolipid_metabolites)
                    
                    # Average z-scores if multiple
                    if len(sphingolipid_stats_list) > 1:
                        z_scores = [s.get('z_score', 0.0) for s in sphingolipid_stats_list]
                        merged_stats['z_score'] = sum(z_scores) / len(z_scores)
                    
                    consolidated_pathways[consolidated_name] = merged_stats
                    self._log_network(f"   • Merged into '{consolidated_name}' with {len(sphingolipid_metabolites)} unique metabolites\n")
                
                # Use consolidated pathways for further processing
                pathways_data = consolidated_pathways
                self._log_network("\n")
                
                # CRITICAL: Deduplicate pathways by normalized name (case-insensitive, whitespace-normalized)
                # Keep the entry with the highest |z-score| for duplicates (prioritize biological effect over significance)
                self._log_network(f"🔍 Deduplicating pathways by normalized name...\n")
                normalized_pathways = {}
                duplicate_count = 0
                for pathway_name, stats in pathways_data.items():
                    if not isinstance(stats, dict):
                        continue
                    # Normalize: lowercase, strip whitespace, collapse multiple spaces
                    normalized_name = ' '.join(pathway_name.lower().strip().split())
                    
                    if normalized_name in normalized_pathways:
                        duplicate_count += 1
                        # Keep the one with higher |z-score| (stronger biological effect)
                        existing_z = abs(normalized_pathways[normalized_name][1].get('z_score', 0.0))
                        current_z = abs(stats.get('z_score', 0.0))
                        if current_z > existing_z:
                            # Current entry has stronger effect, replace
                            self._log_network(f"   • Duplicate found: '{pathway_name}' (kept |z|={current_z:.3f}, discarded |z|={existing_z:.3f})\n")
                            normalized_pathways[normalized_name] = (pathway_name, stats)
                        else:
                            # Existing entry has stronger effect, keep it
                            self._log_network(f"   • Duplicate found: '{pathway_name}' (kept existing |z|={existing_z:.3f}, discarded |z|={current_z:.3f})\n")
                    else:
                        normalized_pathways[normalized_name] = (pathway_name, stats)
                
                if duplicate_count > 0:
                    self._log_network(f"   • Removed {duplicate_count} duplicate pathway entries\n")
                
                # Reconstruct pathways_data from deduplicated entries
                pathways_data = {name: stats for name, stats in normalized_pathways.values()}
                self._log_network(f"   • Unique pathways after deduplication: {len(pathways_data)}\n\n")
                
                # CRITICAL: Apply pathway p-value threshold for DISPLAY filtering
                self._log_network(f"🔍 Applying pathway p-value threshold for display: ≤ {pathway_p_threshold}\n")
                original_count = len(pathways_data)
                filtered_pathways_data = {}
                for pathway_name, stats in pathways_data.items():
                    if isinstance(stats, dict):
                        pathway_pval = stats.get('combined_pvalue', 1.0)
                        if pathway_pval <= pathway_p_threshold:
                            filtered_pathways_data[pathway_name] = stats
                
                filtered_count = len(filtered_pathways_data)
                removed_count = original_count - filtered_count
                self._log_network(f"   • Total pathways calculated: {original_count}\n")
                self._log_network(f"   • Pathways passing threshold (p≤{pathway_p_threshold}): {filtered_count}\n")
                self._log_network(f"   • Pathways filtered out: {removed_count}\n\n")
                
                # CRITICAL: Filter out disease pathways and overly general pathways BEFORE applying max_pathways limit
                # This ensures unwanted pathways don't count against the top N limit
                try:
                    from main_script.metabolite_pathway_network import should_filter_pathway
                except ImportError:
                    def should_filter_pathway(pathway_name):
                        """Minimal fallback filtering"""
                        if not pathway_name:
                            return True
                        pathway_lower = pathway_name.lower()
                        disease_keywords = ['disease', 'syndrome', 'disorder', 'deficiency', 'cancer']
                        return any(kw in pathway_lower for kw in disease_keywords)
                
                filtered_out_pathways = {}
                non_disease_pathways = {}
                for pathway_name, stats in filtered_pathways_data.items():
                    if should_filter_pathway(pathway_name):
                        filtered_out_pathways[pathway_name] = stats
                    else:
                        non_disease_pathways[pathway_name] = stats
                
                if filtered_out_pathways:
                    self._log_network(f"🏥 Disease/General pathways excluded from top N: {len(filtered_out_pathways)}\n")
                    self._log_network(f"   • Examples: {list(filtered_out_pathways.keys())[:3]}\n")
                    self._log_network(f"   • Specific pathways remaining: {len(non_disease_pathways)}\n\n")
                
                # Now apply max pathways limit to non-disease pathways only
                # CRITICAL: Prioritize Activated/Inhibited pathways, then fill remaining slots with No Change
                if limit_enabled and len(non_disease_pathways) > max_pathways:
                    self._log_network(f"🔍 Applying max pathways limit: {max_pathways} (disease and general pathways already excluded)\n")
                    
                    # Separate by status
                    activated_inhibited = {}
                    no_change = {}
                    for p_name, p_stats in non_disease_pathways.items():
                        status = p_stats.get('status', 'No Change')
                        z_score = p_stats.get('z_score', 0.0)
                        if status in ('Activated', 'Inhibited') and abs(z_score) > 0:
                            activated_inhibited[p_name] = p_stats
                        else:
                            no_change[p_name] = p_stats
                    
                    # Sort each group by p-value
                    sorted_act_inh = sorted(activated_inhibited.items(), key=lambda x: x[1].get('combined_pvalue', 1.0))
                    sorted_no_change = sorted(no_change.items(), key=lambda x: x[1].get('combined_pvalue', 1.0))
                    
                    # Take all Activated/Inhibited first, then fill remaining slots with No Change
                    filtered_pathways_data = dict(sorted_act_inh[:max_pathways])
                    remaining_slots = max_pathways - len(filtered_pathways_data)
                    
                    if remaining_slots > 0 and sorted_no_change:
                        filtered_pathways_data.update(dict(sorted_no_change[:remaining_slots]))
                        self._log_network(f"   • Took {len(sorted_act_inh[:max_pathways])} Activated/Inhibited + {min(remaining_slots, len(sorted_no_change))} No Change pathways\n")
                    else:
                        self._log_network(f"   • Limited to {len(filtered_pathways_data)} Activated/Inhibited pathways only (max limit reached)\n")
                    self._log_network(f"   • Total: {len(filtered_pathways_data)}/{max_pathways} pathways\n\n")
                else:
                    # Use non-disease pathways without limit
                    filtered_pathways_data = non_disease_pathways
                
                if not filtered_pathways_data or len(filtered_pathways_data) == 0:
                    self.root.after(0, lambda: self._log_network(f"\n⚠️ No pathways passed the display threshold (p≤{pathway_p_threshold}).\n"))
                    self.root.after(0, lambda: messagebox.showwarning(
                        "No Pathways Found",
                        f"No pathways met the display criteria:\n\n"
                        f"• Min metabolites: {min_metabolites}\n"
                        f"• Metabolite p-value: ≤ {p_threshold}\n"
                        f"• Pathway p-value threshold: ≤ {pathway_p_threshold}\n"
                        f"• Z-score threshold: ±{z_threshold}\n\n"
                        f"Try increasing the pathway p-value threshold (currently {pathway_p_threshold})."
                    ))
                    return
                
                # Store enrichment results (use filtered data for display)
                self.pathway_original_pathways_data = pathways_data.copy()  # Keep all calculated pathways
                pathways_data = filtered_pathways_data  # Use filtered for display
                self.pathway_filtered_pathways_data = pathways_data.copy()
                
                # CRITICAL: Update only pathway_filtered_metabolites_data (significant metabolites for display)
                # NEVER overwrite pathway_original_metabolites_data or original_uploaded_raw_data here
                # These should only be set when data is first loaded, not after enrichment
                # Filter display_df to only include significant metabolites
                try:
                    # Apply the same p-value filter used for enrichment
                    significant_mask = display_df['pvalue'] < p_threshold
                    self.pathway_filtered_metabolites_data = display_df[significant_mask].copy()
                    sig_count = len(self.pathway_filtered_metabolites_data)
                    total_count = len(display_df)
                    self._log_network(f"📊 Metabolite data stored for display:\n")
                    self._log_network(f"   • Total metabolites (from display_df): {total_count}\n")
                    self._log_network(f"   • Significant metabolites (p<{p_threshold}): {sig_count}\n")
                    self._log_network(f"   • Non-significant filtered out: {total_count - sig_count}\n")
                    self._log_network(f"   • Original raw data preserved (not overwritten)\n\n")
                except Exception as filter_err:
                    logger.warning(f"Failed to filter metabolites, using full dataset: {filter_err}")
                    self.pathway_filtered_metabolites_data = display_df.copy()
                
                # Count pathway types (after p-value & max-pathways filters)
                total_significant = len(pathways_data)
                activated = sum(1 for p in pathways_data.values() if isinstance(p, dict) and p.get('status') == 'Activated')
                inhibited = sum(1 for p in pathways_data.values() if isinstance(p, dict) and p.get('status') == 'Inhibited')
                no_change = sum(1 for p in pathways_data.values() if isinstance(p, dict) and p.get('status') == 'No Change')

                # Reproducibility snapshot: log settings and exact pathway list
                self._log_network(f"\n{'='*70}\n")
                self._log_network("📎 ENRICHMENT RUN SNAPSHOT (for reproducibility)\n")
                self._log_network(f"   • Method: {method}\n")
                self._log_network(f"   • min_metabolites={min_metabolites}\n")
                self._log_network(f"   • metabolite_pvalue_threshold={p_threshold}\n")
                self._log_network(f"   • metabolite_log2fc_threshold={log2fc_threshold}\n")
                self._log_network(f"   • pathway_p_threshold={pathway_p_threshold}\n")
                self._log_network(f"   • z_threshold={z_threshold}\n")
                self._log_network(f"   • FDR={'ON' if apply_fdr else 'OFF'}\n")
                self._log_network(f"   • total_shown={total_significant}, activated={activated}, inhibited={inhibited}, no_change={no_change}\n")

                # Log per-pathway details in a deterministic order (sorted by adjusted/combined p-value, then name)
                try:
                    def _get_sort_key(item):
                        name, stats = item
                        if not isinstance(stats, dict):
                            return (1.0, str(name))
                        p_adj = stats.get('adjusted_pvalue')
                        if p_adj is None:
                            p_adj = stats.get('combined_pvalue', 1.0)
                        try:
                            p_adj = float(p_adj)
                        except Exception:
                            p_adj = 1.0
                        return (p_adj, str(name))

                    sorted_items = sorted(pathways_data.items(), key=_get_sort_key)
                    self._log_network("   • Pathways (sorted by adj/combined p-value):\n")
                    for name, stats in sorted_items:
                        if not isinstance(stats, dict):
                            self._log_network(f"      - {name} (non-dict stats, raw={stats})\n")
                            continue
                        status = stats.get('status', 'NA')
                        z_val = stats.get('z_score', stats.get('Z', 0.0))
                        try:
                            z_str = f"{float(z_val):.3f}"
                        except Exception:
                            z_str = str(z_val)
                        p_adj = stats.get('adjusted_pvalue', None)
                        p_comb = stats.get('combined_pvalue', stats.get('pvalue', None))
                        try:
                            p_adj_str = f"{float(p_adj):.2e}" if p_adj is not None else "NA"
                        except Exception:
                            p_adj_str = str(p_adj)
                        try:
                            p_comb_str = f"{float(p_comb):.2e}" if p_comb is not None else "NA"
                        except Exception:
                            p_comb_str = str(p_comb)
                        n_mets = stats.get('n_metabolites', stats.get('count', 'NA'))
                        self._log_network(
                            f"      - {name} | status={status}, z={z_str}, adj_p={p_adj_str}, p={p_comb_str}, n_metabolites={n_mets}\n"
                        )
                except Exception as snapshot_err:
                    # Never break enrichment if snapshot logging fails
                    self._log_network(f"   ⚠️ Snapshot logging error: {snapshot_err}\n")

                self._log_network(f"\n{'='*70}\n")
                self._log_network(f"✅ ENRICHMENT ANALYSIS COMPLETE\n")
                self._log_network(f"{'='*70}\n")
                # Note: Detailed summary that matches the table exactly is logged
                # after _populate_pathway_tree runs, using the displayed counts.
                
                # Build relationship mappings for cascading selections
                # This MUST complete BEFORE populating tables, so call it directly
                self._log_network(f"📋 Building relationship mappings...\n")
                self._build_relationship_mappings(display_df, pathways_data)
                
                # Populate ALL tables (pathway, metabolite, upstream, disease)
                # Stagger the calls to avoid UI freezing
                self._log_network(f"📋 Populating all data tables...\n")
                self.root.after(0, self._populate_pathway_tree)
                self.root.after(100, self._populate_metabolite_tree)
                self.root.after(200, self._populate_upstream_tree)
                self.root.after(300, self._populate_disease_tree)
                self.root.after(400, self._update_stats_labels)
                
                # Enable Generate button
                if hasattr(self, 'generate_network_button'):
                    self.root.after(0, lambda: self.generate_network_button.config(state='normal', bg='#3498db'))
                
                self._log_network(f"✅ Ready to generate network visualization!\n")
                self._log_network(f"   Click 'Generate Network' button to create interactive network.\n\n")
                
                # Show completion message using the counts that will appear in the table.
                # We compute these lazily when the callback runs, after the table is populated.
                def _show_enrichment_popup():
                    try:
                        displayed = getattr(self, '_displayed_pathway_counts', None) or {}
                        shown_total = int(displayed.get('total', 0))
                        shown_activated = int(displayed.get('activated', 0))
                        shown_inhibited = int(displayed.get('inhibited', 0))

                        messagebox.showinfo(
                            "Enrichment Complete",
                            f"Pathway enrichment analysis completed!\n\n"
                            f"• {shown_total} significant pathways shown with current settings\n"
                            f"• Activated: {shown_activated}\n"
                            f"• Inhibited: {shown_inhibited}\n\n"
                            f"Click 'Generate Network' to visualize results."
                        )
                    except Exception:
                        # Fallback: show a generic message
                        messagebox.showinfo(
                            "Enrichment Complete",
                            "Pathway enrichment analysis completed! Click 'Generate Network' to visualize results."
                        )

                self.root.after(0, _show_enrichment_popup)
                
            except Exception as e:
                error_msg = f"Enrichment calculation failed: {str(e)}"
                self._log_network(f"\n❌ ERROR: {error_msg}\n")
                import traceback
                self._log_network(f"{traceback.format_exc()}\n")
                self.root.after(0, lambda: messagebox.showerror("Enrichment Error", error_msg))
        
        # Run in background thread
        self._log_network(f"Starting enrichment calculation in background...\n")
        threading.Thread(target=_enrichment_worker, daemon=True).start()

    def start_generate_network(self):
        """Run network generation in a background thread to keep the GUI responsive."""
        
        # ═══════════════════════════════════════════════════════════════════
        # VALIDATION: At least one top-level output must be selected
        # ═══════════════════════════════════════════════════════════════════
        gen_network = getattr(self, 'network_gen_network', tk.BooleanVar(value=True)).get()
        gen_plots = getattr(self, 'network_gen_plots', tk.BooleanVar(value=True)).get()
        
        if not gen_network and not gen_plots:
            messagebox.showwarning(
                "Selection Required",
                "Please select at least one output option:\n\n"
                "• Generate Network Graph\n"
                "• Generate Statistics\n\n"
                "At least one must be checked to proceed."
            )
            return
        
        # ═══════════════════════════════════════════════════════════════════
        # CYTOSCAPE PROMPT: Only show if Cytoscape export is enabled
        # ═══════════════════════════════════════════════════════════════════
        export_cytoscape = getattr(self, 'export_cytoscape', tk.BooleanVar(value=False)).get()
        
        if export_cytoscape:
            try:
                proceed = messagebox.askyesno(
                    title="Open Cytoscape",
                    message=(
                        "Please ensure Cytoscape is open before generating the network.\n\n"
                        "- Open Cytoscape now\n"
                        "- Return here and click Proceed\n\n"
                        "Click Yes to Proceed or No to Cancel."
                    ),
                    icon=messagebox.WARNING
                )
                if not proceed:
                    return
            except Exception:
                # If messagebox fails for any reason, continue silently
                pass
        btn = getattr(self, 'generate_network_button', None)
        if btn:
            try:
                btn.config(state='disabled', text='⏳ Generating...')
            except Exception:
                pass
        self._log_network("Starting background network generation...\n")
        def _worker():
            try:
                # Silence third-party logs during generation to surface our debugs
                self._silence_third_party_logs(True)
                self.generate_network_visualization()
            finally:
                # Restore logging levels
                try:
                    self._silence_third_party_logs(False)
                except Exception:
                    pass
                if btn:
                    self.frame.after(0, lambda: btn.config(state='normal', text='🚀 Generate Network') if btn else None)  # type: ignore
        threading.Thread(target=_worker, daemon=True).start()
