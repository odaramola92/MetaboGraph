import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np
import re
import colorsys
from scipy import stats
from sklearn.cluster import AgglomerativeClustering
from scipy.stats import pearsonr
import re
import logging
from collections import defaultdict
# silence py4cytoscape verbose responses and related HTTP debug noise
for name in ('py4cytoscape', 'py4cytoscape.cyrest', 'requests', 'urllib3'):
    lg = logging.getLogger(name)
    lg.setLevel(logging.WARNING)
    lg.propagate = False
    
# Import Fisher ORA module (parallel implementation for comparison)
try:
    from main_script.fisher_ora_pathway_analysis import calculate_pathway_statistics_fisher_ora
    FISHER_ORA_AVAILABLE = True
except ImportError:
    FISHER_ORA_AVAILABLE = False
    print("⚠️  Fisher ORA module not available - using standard method only")

          
# Simple in-memory cache for pathway statistics to avoid repeated recalculation
# Keyed by (n_rows, tuple(columns), sum_log2fc, sum_pvalues, params_tuple)
_PATHWAY_STATS_CACHE = {}

# --- Disease pathway detection patterns (module-level) ---
# Define disease suffix patterns that match at end OR followed by type/variant info
_DISEASE_SUFFIX_KEYWORDS = ["emia", "uria", "osis"]
# Compile regex pattern: suffix must be followed by whitespace or end-of-string (case-insensitive)
# This catches "Tyrosinemia" and "Tyrosinemia Type I" but not "Metabolism" (osis in middle)
_DISEASE_SUFFIX_PATTERN = re.compile(r'(' + '|'.join(map(re.escape, _DISEASE_SUFFIX_KEYWORDS)) + r')(?:\s|$)', re.IGNORECASE)

# Hardcoded pathway name corrections/normalizations (maps variations to canonical disease name)
_DISEASE_PATHWAY_CORRECTIONS = {
    # Tyrosinemia variations - normalize to standard name
    "tyrosinemia, transient, of the newborn": "Transient Tyrosinemia of the Newborn",
    "transient;of the newborn": "Transient Tyrosinemia of the Newborn",
    "tyrosinemia,transient;of the newborn": "Transient Tyrosinemia of the Newborn",
}

# Pathway-context keywords that should override disease keywords when the phrase
# is clearly a biological pathway rather than a disease label.
_PATHWAY_OVERRIDE_KEYWORDS = [
    "metabolism",
    "metabolic",
    "biosynthesis",
    "synthesis",
    "degradation",
    "catabolism",
    "anabolism",
    "transport",
    "cycle",
    "pathway",
    "signaling",
    "oxidation",
    "reduction",
    "homeostasis",
    "processing",
    "uptake",
    "turnover",
]

# ============================================================================
# CENTRALIZED PATHWAY FILTERING ARCHITECTURE
# ============================================================================
# All pathway filtering (disease detection, overly general exclusions) is
# centralized in this module to ensure consistency across the entire pipeline:
#
#   1. Annotation Stage (metabolite_pathways_annotator.py)
#      → Filters pathways during database lookup/annotation
#
#   2. Network Stage (pathway_network_tab.py)
#      → Filters pathways before displaying in tree view
#
#   3. Enrichment Stage (pathway_network_tab.py enrichment thread)
#      → Filters pathways after Fisher ORA calculation
#
# USE should_filter_pathway() for all filtering - it combines:
#   • is_disease_pathway_global() - detects disease pathways
#   • is_overly_general_pathway() - excludes overly broad pathways
#
# This approach eliminates code duplication and ensures all exclusions
# (encephalopathy, x-linked, triosephosphate isomerase, etc.) are applied
# consistently at all stages.
# ============================================================================

def normalize_disease_pathway_name(pathway_name: str) -> str:
    """Normalize disease pathway names using hardcoded corrections"""
    pathway_lower = pathway_name.lower().strip()
    # Check if we have a hardcoded correction for this pathway
    if pathway_lower in _DISEASE_PATHWAY_CORRECTIONS:
        return _DISEASE_PATHWAY_CORRECTIONS[pathway_lower]
    return pathway_name

def is_pathway_overriding_disease(pathway_name: str) -> bool:
    """Return True when pathway-context keywords should override disease matching."""
    pathway_lower = pathway_name.lower()
    return any(keyword in pathway_lower for keyword in _PATHWAY_OVERRIDE_KEYWORDS)

def is_disease_pathway_global(pathway_name):
    """Check if pathway should be treated as a disease (suffix match or keyword match)"""
    pathway_lower = pathway_name.lower()

    # PRIORITY 0: If the name is clearly a pathway phrase, keep it as a pathway
    # even if it contains disease-related words like "cancer".
    if is_pathway_overriding_disease(pathway_name):
        return False
    
    # PRIORITY 1: Check for explicit disease pathways first (before any exclusions)
    # These are definitive disease pathways that should ALWAYS be classified as diseases
    explicit_disease_keywords = [
        "anemia", "homocystinuria", "megaloblastic", "mngie", "defect in cobala",
        "disease", "disorder", "syndrome", "carcinoma", "cancer", "tumor",
        "diabetes", "alzheimer", "parkinson", "huntington", "hypoacetylaspartia",
        "leukodystrophy", "cardiomyopathy", "leishmaniasis", "oncogenic action",
        "copy number variation", "15q25", "2-hydroxyglutarate", "encephalopathy"
    ]
    for keyword in explicit_disease_keywords:
        if keyword in pathway_lower:
            return True  # Definitive disease, return immediately
    
    # PRIORITY 2: Exclude legitimate pathways that contain disease-like keywords but are NOT diseases
    # These are normal metabolic pathways - only check if not already identified as disease above
    pathway_exclusions = [
        "warburg effect",              # Cancer metabolism pathway, not a disease
        "trna charging",               # Normal tRNA aminoacylation process
        "urea cycle",                  # Normal nitrogen disposal pathway
    ]
    
    # Check if pathway matches any exclusion (these are NOT diseases)
    for exclusion in pathway_exclusions:
        if exclusion in pathway_lower:
            return False  # This is a normal pathway, not a disease
    
    # PRIORITY 3: Check suffix patterns (must be at end)
    if _DISEASE_SUFFIX_PATTERN.search(pathway_name):
        return True
    
    # PRIORITY 4: Check remaining disease keywords
    disease_secondary_keywords = [
        "deficiency", "newborn", "type", "transient", "diabetic", "copy number variation"
    ]
    for keyword in disease_secondary_keywords:
        if keyword in pathway_lower:
            return True
    
    return False

def is_overly_general_pathway(pathway_name):
    """Check if pathway is too broad/general and should be excluded from results"""
    pathway_lower = pathway_name.lower().strip()
    
    # Exact matches for overly general pathways
    general_exact = [
        "biochemical pathways: part i",
        "biochemical pathways: part ii",
        "biochemical pathways: part iii",
        "metabolic pathways",
        "metabolism",
        "biosynthesis",
        "degradation",
        "signaling pathways",
    ]
    
    if pathway_lower in general_exact:
        return True
    
    # Specific pathways to exclude (user-requested removals)
    specific_exclusions = [
        "pgk1",
        "bitter",
        "umami taste",
        "and umami taste",
        "sensory perception of sweet",
        "fructose-1",
        "x-linked",
        "triosephosphate isomerase",
        "pleural mesothelioma"
    ]
    
    if pathway_lower in specific_exclusions:
        return True
    
    # Substring matches for specific unwanted pathways
    substring_exclusions = [
        "mesec into h2se",
        "metabolism of selenium into h2se",
        "acylcarnitine",  # Not real pathways, just compound names
    ]
    
    for exclusion in substring_exclusions:
        if exclusion in pathway_lower:
            return True
    
    return False


def should_filter_pathway(pathway_name):
    """
    Centralized pathway filtering function.
    Returns True if pathway should be EXCLUDED (filtered out).
    Checks both disease pathways and overly general pathways.
    
    Use this function in all stages (annotation, network, enrichment) to ensure
    consistent filtering across the entire pipeline.
    
    Args:
        pathway_name: Name of the pathway to check
        
    Returns:
        True if pathway should be excluded, False if it should be kept
    """
    if not pathway_name or pd.isna(pathway_name):
        return True
        
    # Check if it's a disease pathway (should be filtered)
    if is_disease_pathway_global(pathway_name):
        return True
    
    # Check if it's overly general pathway (should be filtered)
    if is_overly_general_pathway(pathway_name):
        return True
    
    return False

def normalize_pathway_name(pathway_name: str) -> str:
    """
    Normalize pathway names to handle case variations and formatting differences.
    This ensures 'Nicotine and Nicotinamide Metabolims' and 'Nicotine and Nicotinamide metabolims'
    are treated as the same pathway.
    
    Returns the normalized name with proper title casing.
    """
    if not pathway_name or pd.isna(pathway_name):
        return ""
    
    # Strip whitespace
    normalized = str(pathway_name).strip()
    
    # Title case the pathway name (capitalizes first letter of each word)
    # This handles: "Nicotine and Nicotinamide Metabolims" vs "Nicotine and Nicotinamide metabolims"
    normalized = normalized.title()
    
    # Handle common acronyms/abbreviations that should stay uppercase
    # Replace after title casing
    common_acronyms = {
        'Atp': 'ATP', 'Adp': 'ADP', 'Amp': 'AMP', 'Gtp': 'GTP', 'Gdp': 'GDP',
        'Nadh': 'NADH', 'Nad+': 'NAD+', 'Nadph': 'NADPH', 'Nadp+': 'NADP+',
        'Fadh2': 'FADH2', 'Fad': 'FAD', 'Dna': 'DNA', 'Rna': 'RNA',
        'Trna': 'tRNA', 'Mrna': 'mRNA', 'Rrna': 'rRNA',
        'Coq': 'CoQ', 'Coa': 'CoA', 'Smp': 'SMP',
        'Tca': 'TCA', 'Etc': 'ETC', 'Ppp': 'PPP',
        'Hiv': 'HIV', 'Hpv': 'HPV', 'Hsv': 'HSV',
        'Dna': 'DNA', 'Rna': 'RNA', 'Atp': 'ATP',
        'Ii': 'II', 'Iii': 'III', 'Iv': 'IV', 'Vi': 'VI', 'Vii': 'VII',
    }
    
    for wrong, correct in common_acronyms.items():
        # Use word boundaries to avoid replacing parts of words
        normalized = re.sub(r'\b' + re.escape(wrong) + r'\b', correct, normalized)
    
    return normalized

# --- Global pathway name normalization helper (used across functions) ---
def normalize_pathway_name_global(pathway_name):
    """Normalize pathway names to merge duplicates with slight variations.

    This mirrors the normalization used during pathway stats computation so we can
    reliably match pathway names found in the input rows with the keys present in
    the computed pathway_stats dict.
    """
    if not isinstance(pathway_name, str):
        return str(pathway_name)

    # Lowercase + strip
    normalized = pathway_name.lower().strip()

    # Normalize common patterns and spacing
    normalized = normalized.replace('beta-', 'beta ')
    normalized = normalized.replace('beta_', 'beta ')
    normalized = normalized.replace(' pathway', '')
    normalized = normalized.replace(' pathways', '')
    normalized = normalized.replace('-chain', ' chain')
    normalized = normalized.replace('branched-chain', 'branched chain')
    normalized = normalized.replace('branched_chain', 'branched chain')
    normalized = ' '.join(normalized.split())

    # Title case for consistency
    normalized = normalized.title()

    # Handle specific known duplicates
    if 'pyrimidine metabolism' in normalized.lower():
        normalized = 'Pyrimidine Metabolism'
    elif 'beta ureidopropionase deficiency' in normalized.lower():
        normalized = 'Beta Ureidopropionase Deficiency'
    elif 'oxidation of branched' in normalized.lower() and 'fatty acid' in normalized.lower():
        normalized = 'Oxidation Of Branched Chain Fatty Acids'

    return normalized


def build_class_level_network_dataframe(
    df: pd.DataFrame,
    class_col: str = 'Main_Class',
    compress_classes: bool = True,
    name_candidates=None,
    pathway_col_candidates=None,
):
    """Aggregate lipid rows to the class level when class-aware mode is enabled.

    Parameters
    ----------
    df : pandas.DataFrame
        Annotated metabolite dataframe already normalized by the GUI.
    class_col : str
        Column containing the lipid class display name (default 'Main_Class').
    compress_classes : bool
        Whether class compression is requested. When False, the dataframe is
        returned untouched and the mapping table is ``None``.
    name_candidates : list[str], optional
        Ordered list of candidate identifier columns. Defaults to the same
        candidates used elsewhere in the network module.
    pathway_col_candidates : list[str], optional
        Ordered list of candidate pathway columns to aggregate. Defaults to the
        standard ``All_Pathways``→``Metabolika_Pathways`` priority.

    Returns
    -------
    (processed_df, mapping_df)
        ``processed_df`` is safe to feed to the network code. ``mapping_df`` is
        a long-form table describing Pathway→Class→Member connections that can
        be displayed or exported for audit.
    """

    if df is None or not compress_classes or df.empty:
        return df, None

    if class_col not in df.columns:
        legacy_class_col = 'Main_Class'
        if legacy_class_col in df.columns:
            class_col = legacy_class_col
        else:
            return df, None

    if name_candidates is None:
        name_candidates = [
            'Name', 'LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID',
            'Metabolite', 'Metabolite_Name', 'Feature', 'Identifier'
        ]

    if pathway_col_candidates is None:
        pathway_col_candidates = [
            'All_Pathways', 'All_Pathways_Display', 'Pathways', 'HMDB_Pathways',
            'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways',
            'Metabolika_Pathways', 'Reactome_Pathways', 'KEGG_Pathways'
        ]

    name_col = next((c for c in name_candidates if c in df.columns), None)
    if not name_col:
        name_col = 'Name' if 'Name' in df.columns else df.columns[0]

    pathway_col = next((c for c in pathway_col_candidates if c in df.columns), None)
    if pathway_col is None:
        return df, None

    class_values = df[class_col].astype(str).str.strip()
    valid_mask = df[class_col].notna() & class_values.str.len().gt(0) & class_values.str.lower().ne('na')
    if not valid_mask.any():
        return df, None

    compressible_df = df.loc[valid_mask].copy()
    individual_df = df.loc[~valid_mask].copy()

    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    disease_cols = [col for col in ['Associated_Diseases', 'Reactome_Disease'] if col in df.columns]
    upstream_cols = [
        col for col in [
            'Enzyme_Gene_name', 'Transporter_Gene_name', 'Transporters',
            'Transporter_Uniprot_ID', 'Enzymes_Accession'
        ] if col in df.columns
    ]
    pathway_source_cols = [col for col in ['pathway_sources', 'Pathway_sources'] if col in df.columns]

    split_pattern = re.compile(r'[|;,]+')

    def _unique(sequence):
        seen = set()
        ordered = []
        for item in sequence:
            if item not in seen:
                ordered.append(item)
                seen.add(item)
        return ordered

    text_collapse_cols = _unique(
        disease_cols + upstream_cols + pathway_source_cols + ['All_Pathways_Display']
    )

    def _tokenize(value):
        if value is None or (isinstance(value, (float, np.floating)) and np.isnan(value)):
            return []
        if isinstance(value, (list, tuple, set)):
            iterable = value
        else:
            text = str(value).strip()
            lowered = text.lower()
            if not text or lowered in {'nan', '<na>'}:
                return []
            iterable = split_pattern.split(text)
        tokens = [str(item).strip() for item in iterable if str(item).strip()]
        return tokens

    def _tokenize_series(series):
        collected = []
        for value in series.dropna().tolist():
            collected.extend(_tokenize(value))
        return collected

    class_rows = []
    mapping_rows = []

    for class_name, block in compressible_df.groupby(class_col):
        block = block.copy()
        agg_row = block.iloc[0].copy()
        member_names = [
            str(val).strip()
            for val in block[name_col].astype(str).tolist()
            if str(val).strip() and str(val).lower() != 'nan'
        ]
        if not member_names:
            member_names = [class_name]

        agg_row['Name'] = class_name
        agg_row[class_col] = class_name
        agg_row['Class_members'] = '|'.join(_unique(member_names))
        agg_row['Class_member_count'] = len(member_names)

        if numeric_cols:
            numeric_block = block[numeric_cols].apply(pd.to_numeric, errors='coerce')
            agg_row.loc[numeric_cols] = numeric_block.mean(skipna=True)

        pathway_member_map = defaultdict(set)
        for _, row in block.iterrows():
            member = str(row.get(name_col, '')).strip()
            if member.lower() == 'nan':
                member = ''
            for pathway in _tokenize(row.get(pathway_col)):
                pathway_member_map[pathway].add(member)

        aggregated_pathways = _unique([p for p in pathway_member_map.keys() if p])
        agg_row[pathway_col] = '; '.join(aggregated_pathways)
        if 'All_Pathways_Display' in block.columns:
            agg_row['All_Pathways_Display'] = '; '.join(aggregated_pathways)

        column_token_cache = {}
        for col in text_collapse_cols:
            if col in block.columns:
                tokens = _unique(_tokenize_series(block[col]))
                column_token_cache[col] = tokens
                if tokens:
                    agg_row[col] = '; '.join(tokens)

        diseases = _unique([
            token for col in disease_cols for token in column_token_cache.get(col, [])
        ]) if disease_cols else []
        upstreams = _unique([
            token for col in upstream_cols for token in column_token_cache.get(col, [])
        ]) if upstream_cols else []
        source_tokens = _unique([
            token for col in pathway_source_cols for token in column_token_cache.get(col, [])
        ]) if pathway_source_cols else []

        class_rows.append(agg_row)

        for pathway_name, members in pathway_member_map.items():
            if not pathway_name:
                continue
            member_list = [m for m in members if m]
            mapping_rows.append({
                'Pathway': pathway_name,
                'Class': class_name,
                'Member_Lipids': '|'.join(_unique(member_list)) or agg_row['Class_members'],
                'Class_Member_Count': len(member_names),
                'Pathway_Member_Count': len(member_list),
                'Diseases': '; '.join(diseases),
                'Upstreams': '; '.join(upstreams),
                'Sources': '; '.join(source_tokens),
            })

    class_df = pd.DataFrame(class_rows)
    combined = pd.concat([class_df, individual_df], ignore_index=True, sort=False)

    ordered_cols = list(df.columns)
    for extra_col in ['Class_members', 'Class_member_count']:
        if extra_col in combined.columns and extra_col not in ordered_cols:
            ordered_cols.append(extra_col)
    combined = combined.loc[:, ordered_cols]

    mapping_df = pd.DataFrame(mapping_rows) if mapping_rows else None
    return combined, mapping_df


def wrap_pathway_name(name: str, max_chars: int = 45, max_lines: int = None, break_token: str = '<br>') -> str:
    """
    Wrap long pathway names by inserting line breaks at word boundaries.
    
    Parameters:
    -----------
    name : str
        The pathway name to wrap
    max_chars : int
        Maximum characters per line before wrapping (default: 45)
    
    Returns:
    --------
    str
        The pathway name with '<br>' inserted at appropriate word boundaries
    """
    if not name or len(name) <= max_chars:
        return name
    
    words = name.split()
    lines = []
    remaining_words = words[:]

    def _split_long_word(word: str, chunk_size: int):
        if chunk_size <= 0:
            return [word]
        return [word[i:i + chunk_size] for i in range(0, len(word), chunk_size)]

    while remaining_words:
        if max_lines is not None and len(lines) >= max_lines - 1:
            lines.append(' '.join(remaining_words))
            break

        current_line = []
        current_length = 0

        while remaining_words:
            word = remaining_words[0]
            if len(word) > max_chars:
                chunks = _split_long_word(word, max_chars)
                remaining_words.pop(0)
                for chunk in reversed(chunks[1:]):
                    remaining_words.insert(0, chunk)
                word = chunks[0]
            separator_length = 1 if current_line else 0
            next_length = current_length + separator_length + len(word)

            if current_line and next_length > max_chars:
                break

            current_line.append(word)
            current_length = next_length
            remaining_words.pop(0)

        if current_line:
            lines.append(' '.join(current_line))
        else:
            lines.append(remaining_words.pop(0))

        if max_lines is not None and len(lines) >= max_lines and remaining_words:
            overflow = ' '.join(remaining_words)
            if len(lines[-1]) + len(overflow) + 1 > max_chars:
                if max_chars > 3:
                    lines[-1] = lines[-1][:max_chars - 3].rstrip() + '...'
                else:
                    lines[-1] = lines[-1][:max_chars]
            else:
                lines[-1] = (lines[-1] + ' ' + overflow).strip()
            break
    
    return break_token.join(lines)


def import_data(file_path1, sheet_name1=None, 
                pvalue_threshold=0.05, 
                pos_set_fc=0, neg_set_fc=0,
                exclude_pathway_keywords=None,
                name_column=None):
    """
    Read metabolite data from Excel file.
    
    Parameters:
    -----------
    file_path1 : str
        Path to the Excel file
    sheet_name1 : str, optional
        Sheet name to read
    pvalue_threshold : float
        P-value threshold for filtering
    pos_set_fc : float
        Positive log2FC threshold
    neg_set_fc : float
        Negative log2FC threshold
    exclude_pathway_keywords : list, optional
        Keywords to exclude from pathways
    name_column : str, optional
        User-verified name column. If provided, will be used instead of auto-detection.
        Should typically be 'Name' (the standardized name after annotation).
    """
    try:
        # Set default for excluded keywords
        if exclude_pathway_keywords is None:
            exclude_pathway_keywords = []
        elif isinstance(exclude_pathway_keywords, str):
            exclude_pathway_keywords = [exclude_pathway_keywords]
        
        print(f"Import filtering criteria:")
        print(f" - p-value threshold: {pvalue_threshold}")
        print(f" - Positive log2FC threshold: {pos_set_fc}")
        print(f" - Negative log2FC threshold: {neg_set_fc}")
        #if exclude_pathway_keywords:
            #print(f" - Excluded pathway keywords: {exclude_pathway_keywords}")

        df1 = pd.read_excel(file_path1, sheet_name=sheet_name1) if sheet_name1 else pd.read_excel(file_path1)
        # Normalize column names to avoid hidden whitespace/unicode issues
        df1.columns = df1.columns.map(lambda c: str(c).strip())
        print(f"File 1 loaded. Shape: {df1.shape}")

        # Select relevant columns for df1 - accept Name, Metabolite or Molecule
        # Build desired column list and pick only those present to avoid KeyError
        desired_cols = ["Name", "pvalue", "log2FC", "HMDB_Pathways", "PathBank_Pathways", "SMPDB_Pathways",
                        'WikiPathways', 'Metabolika_Pathways', 'Enzymes_Accession', 'Enzyme_Gene_name', 'Reaction_Role',
                        'Transporters', 'Transporter_Gene_name', 'Transporter_Uniprot_ID']

        # Normalize primary identifier column
        # PRIORITY: Use user-verified name_column parameter if provided (removes hardcoding)
        if name_column and name_column in df1.columns:
            if name_column != 'Name':
                df1 = df1.rename(columns={name_column: 'Name'})
                print(f"✓ Using verified name column '{name_column}' → renamed to 'Name'")
        # FALLBACK: Auto-detect for backward compatibility (legacy behavior)
        elif 'Name' in df1.columns:
            pass
        elif 'Metabolite' in df1.columns:
            df1 = df1.rename(columns={'Metabolite': 'Name'})
        elif 'Molecule' in df1.columns:
            df1 = df1.rename(columns={'Molecule': 'Name'})
        elif 'LipidID' in df1.columns:
            # For lipid mode, use LipidID as Name
            df1 = df1.rename(columns={'LipidID': 'Name'})
        elif 'metabolite_id' in df1.columns:
            df1 = df1.rename(columns={'metabolite_id': 'Name'})
        else:
            raise ValueError("First file must contain a column named 'Name', 'Metabolite', 'Molecule', 'LipidID', or provide name_column parameter")

        available_cols = [c for c in desired_cols if c in df1.columns]
        # Ensure we always include the identifier and pvalue/log2FC if present
        for required in ["Name", "pvalue", "log2FC"]:
            if required in df1.columns and required not in available_cols:
                available_cols.insert(0, required)

        df1 = df1[available_cols]
        print(f"Columns selected successfully")

        df = df1.copy()
        # Drop rows with missing metabolite identifier to avoid NaN nodes later
        if 'Name' in df.columns:
            before = df.shape[0]
            df = df.dropna(subset=['Name'])
            after = df.shape[0]
            if before != after:
                print(f"Dropped {before - after} rows with missing 'Name' values to avoid NaN nodes")
        print(f"Data merged. Shape: {df.shape}")
        print(f"Columns in merged data: {df.columns.tolist()}")

        # Drop rows only if all available pathway columns are missing
        pathway_cols = ['HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways', 'Metabolika_Pathways']
        present_pathway_cols = [c for c in pathway_cols if c in df.columns]
        if present_pathway_cols:
            df = df.dropna(subset=present_pathway_cols, how='all')
        else:
            # No pathway columns present at all; keep df but warn
            print("⚠️ No pathway columns (HMDB_Pathways/PathBank_Pathways/SMPDB_Pathways) found in data.")
        print(f"After dropping rows with no pathways. Shape: {df.shape}")

        # Merge Reactome_Disease column into Associated_Diseases so diseases are
        # available to the network (even when annotation was loaded from file)
        try:
            if 'Reactome_Disease' in df.columns:
                # Normalize existing diseases to list
                def _merge_disease(a: str, b: str) -> str:
                    parts = []
                    for src in (a, b):
                        if pd.notna(src) and str(src).strip() and str(src).lower() != 'nan':
                            if '|' in str(src):
                                parts.extend([x.strip() for x in str(src).split('|') if x.strip()])
                            else:
                                parts.extend([x.strip() for x in str(src).split(';') if x.strip()])
                    # Deduplicate, preserve order
                    seen = set()
                    uniq = []
                    for p in parts:
                        pl = p.lower()
                        if pl not in seen:
                            seen.add(pl)
                            uniq.append(p)
                    return ' | '.join(uniq)

                existing = df['Associated_Diseases'] if 'Associated_Diseases' in df.columns else ''
                df['Associated_Diseases'] = [
                    _merge_disease(existing.iloc[i] if isinstance(existing, pd.Series) else '', rd)
                    for i, rd in enumerate(df['Reactome_Disease'])
                ]
        except Exception as _e:
            # Non-fatal; continue without blocking the pipeline
            print(f"WARN: Could not merge Reactome_Disease into Associated_Diseases: {_e}")
        
        # Define disease suffix patterns that must match at the end of pathway names
        disease_suffix_keywords = ["emia", "uria", "osis"]
        # Compile regex pattern for end-of-string matching (case-insensitive)
        disease_suffix_pattern = re.compile(r'(' + '|'.join(map(re.escape, disease_suffix_keywords)) + r')$', re.IGNORECASE)
        
        def is_disease_pathway(pathway_name):
            """Check if pathway should be treated as a disease (suffix match or keyword match)"""
            pathway_lower = pathway_name.lower()
            # Check suffix patterns (must be at end)
            if disease_suffix_pattern.search(pathway_name):
                return True
            # Check general disease keywords (substring match)
            disease_general_keywords = ["disease", "disorder", "syndrome", "deficiency"]
            for keyword in disease_general_keywords:
                if keyword in pathway_lower:
                    return True
            return False
        
        def should_exclude_pathway(pathway_name):
            """Check if pathway should be excluded based on keywords (non-mammalian, plant, etc.)"""
            if not exclude_pathway_keywords:
                return False
            pathway_lower = pathway_name.lower()
            for keyword in exclude_pathway_keywords:
                if keyword in pathway_lower:
                    return True
            return False
        
        # Split pathway columns into lists with keyword filtering
        def split_pathways(row):
            # First, check if All_Pathways already exists as a list (from new annotator)
            all_pathways = row.get('All_Pathways')
            if isinstance(all_pathways, list) and all_pathways:
                pathways = all_pathways.copy()
            else:
                # Fallback: build from individual pathway columns
                pathways = []
                # HMDB pathways use '|' separator (with spaces: ' | ')
                hmdb_val = row.get('HMDB_Pathways') if isinstance(row, dict) else row.get('HMDB_Pathways', None)
                if pd.notna(hmdb_val):
                    # Handle both '|' and ' | ' separators
                    hmdb_str = str(hmdb_val)
                    if ' | ' in hmdb_str:
                        pathways += [p.strip() for p in hmdb_str.split(' | ') if p.strip()]
                    else:
                        pathways += [p.strip() for p in hmdb_str.split('|') if p.strip()]

                # PathBank pathways are joined with '; ' in the annotator, so split by ';'
                pathbank_val = row.get('PathBank_Pathways') if isinstance(row, dict) else row.get('PathBank_Pathways', None)
                if pd.notna(pathbank_val):
                    pathways += [p.strip() for p in str(pathbank_val).split(';') if p.strip()]

                # SMPDB pathways use '|' separator
                smpdb_val = row.get('SMPDB_Pathways') if isinstance(row, dict) else row.get('SMPDB_Pathways', None)
                if pd.notna(smpdb_val):
                    pathways += [p.strip() for p in str(smpdb_val).split('|') if p.strip()]

            # WikiPathways are joined with '; ' in the annotator, so split by ';'
            wikipath_val = row.get('WikiPathways') if isinstance(row, dict) else row.get('WikiPathways', None)  
            if pd.notna(wikipath_val):
                pathways += [p.strip() for p in str(wikipath_val).split(';') if p.strip()]  

            metabolika_val = row.get('Metabolika_Pathways') if isinstance(row, dict) else row.get('Metabolika_Pathways', None)
            if pd.notna(metabolika_val):
                pathways += [p.strip() for p in str(metabolika_val).split(';') if p.strip()]

            # Normalize pathway names and remove case-insensitive duplicates
            # This handles: "Nicotine and Nicotinamide Metabolims" vs "Nicotine and Nicotinamide metabolims"
            normalized_pathways = []
            seen_normalized = {}  # Maps normalized lowercase name to canonical (first seen) name
            
            for pathway in pathways:
                if not pathway:
                    continue
                
                normalized = normalize_pathway_name(pathway)
                normalized_lower = normalized.lower()
                
                # Check if we've seen this pathway name (case-insensitive)
                if normalized_lower not in seen_normalized:
                    # First time seeing this pathway - use the normalized name
                    normalized_pathways.append(normalized)
                    seen_normalized[normalized_lower] = normalized
            
            # Replace pathways with normalized, deduplicated list
            pathways = normalized_pathways

            # Extract disease-related pathways and add to Associated_Diseases column
            # instead of filtering them out
            if exclude_pathway_keywords:
                filtered_pathways = []
                disease_pathways = []
                for pathway in pathways:
                    # Check if pathway is disease-related (using suffix pattern or disease keywords)
                    if is_disease_pathway(pathway):
                        disease_pathways.append(pathway)
                    # Check if pathway should be excluded (non-mammalian, plant, etc.)
                    elif should_exclude_pathway(pathway):
                        # Skip this pathway entirely (don't add to filtered or disease)
                        continue
                    else:
                        filtered_pathways.append(pathway)
                
                # Store disease pathways in a temporary column to merge later
                if disease_pathways:
                    # Get the metabolite name from the row data
                    metabolite_name = row.get('Name', '')
                    if metabolite_name and metabolite_name in df.index:
                        idx = metabolite_name
                    else:
                        # If not in index, skip storing to temp column (data structure issue)
                        idx = None
                    
                    if metabolite_name and idx and idx in df.index:
                        # Merge with existing Associated_Diseases
                        existing_diseases = df.at[idx, 'Associated_Diseases'] if 'Associated_Diseases' in df.columns else ''
                        existing_list = []
                        if pd.notna(existing_diseases) and str(existing_diseases).strip():
                            # Split by | or ; separator
                            if '|' in str(existing_diseases):
                                existing_list = [d.strip() for d in str(existing_diseases).split('|') if d.strip()]
                            else:
                                existing_list = [d.strip() for d in str(existing_diseases).split(';') if d.strip()]
                        
                        # Combine and deduplicate
                        all_diseases = existing_list + disease_pathways
                        unique_diseases = []
                        seen = set()
                        for d in all_diseases:
                            d_lower = d.lower()
                            if d_lower not in seen:
                                unique_diseases.append(d)
                                seen.add(d_lower)
                        
                        df.at[idx, 'Associated_Diseases'] = ' | '.join(unique_diseases)
                
                pathways = filtered_pathways
            
            return pathways

        df['All_Pathways'] = df.apply(split_pathways, axis=1)
        
        # Calculate pathway filtering statistics
        if exclude_pathway_keywords:
            total_unique_pathways_before = 0
            total_unique_pathways_after = 0
            excluded_pathways = set()
            
            # Collect all pathways to calculate filtering statistics
            for _, row in df.iterrows():
                # Get original pathways (without filtering)
                original_pathways = []
                for col in ['HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways', 'Metabolika_Pathways']:
                    if col in df.columns and pd.notna(row.get(col)):
                        if col == 'HMDB_Pathways' or col == 'SMPDB_Pathways':
                            original_pathways += [p.strip() for p in str(row[col]).split('|') if p.strip()]
                        else:
                            original_pathways += [p.strip() for p in str(row[col]).split(';') if p.strip()]
                
                # Check which pathways were excluded
                for pathway in original_pathways:
                    for keyword in exclude_pathway_keywords:
                        if keyword.lower() in pathway.lower():
                            excluded_pathways.add(pathway)
                            break
            
            total_unique_pathways_before = len(set([p for pathways in [
                [p.strip() for p in str(row.get(col, '')).split('|' if col in ['HMDB_Pathways', 'SMPDB_Pathways'] else ';') if p.strip()]
                for col in ['HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways', 'Metabolika_Pathways']
                if col in df.columns
            ] for _, row in df.iterrows() for pathways in [[p.strip() for p in str(row.get(col, '')).split('|' if col in ['HMDB_Pathways', 'SMPDB_Pathways'] else ';') if p.strip()] for col in ['HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways', 'Metabolika_Pathways'] if col in df.columns] for p in pathways if p]))
            
            total_unique_pathways_after = df['All_Pathways'].explode().nunique()
            
            print(f"Pathway filtering summary:")
            print(f" - Total unique pathways before filtering: {total_unique_pathways_before}")
            print(f" - Total unique pathways after filtering: {total_unique_pathways_after}")
            print(f" - Pathways excluded: {len(excluded_pathways)}")
            if excluded_pathways:
                print(f" - Sample excluded pathways: {list(excluded_pathways)[:5]}{'...' if len(excluded_pathways) > 5 else ''}")
        
        print("Data imported successfully with shape:", df.shape)
        print(f"Number of Metabolites with p-value < {pvalue_threshold}:", df[df['pvalue'] < pvalue_threshold].shape[0])

        # Filter metabolites with p-value < threshold
        pvalue_filtered_df = df[df['pvalue'] < pvalue_threshold].copy()

        # Count unique pathways among these metabolites
        num_pathways_pval = pvalue_filtered_df['All_Pathways'].explode().nunique()
        print(f"Number of pathways with p-value < {pvalue_threshold}:", num_pathways_pval)

        print("Data imported successfully after p-value filtering with shape:", pvalue_filtered_df.shape)
        
        # Count upregulated and downregulated metabolites from p-value filtered data
        upregulated_count = (pvalue_filtered_df['log2FC'] >= pos_set_fc).sum()
        downregulated_count = (pvalue_filtered_df['log2FC'] <= neg_set_fc).sum()

        print(f"Upregulated metabolites (log2FC >= {pos_set_fc}): {upregulated_count}")

        # Count upregulated metabolites that have at least one pathway
        upregulated_with_pathways = pvalue_filtered_df[
            (pvalue_filtered_df['log2FC'] >= pos_set_fc) & 
            (pvalue_filtered_df['All_Pathways'].apply(len) > 0)
        ].shape[0]

        print(f"Upregulated metabolites with at least 1 Pathway: {upregulated_with_pathways}")

        print(f"Downregulated metabolites (log2FC <= {neg_set_fc}): {downregulated_count}")

        # Count downregulated metabolites that have at least one pathway  
        downregulated_with_pathways = pvalue_filtered_df[
            (pvalue_filtered_df['log2FC'] <= neg_set_fc) & 
            (pvalue_filtered_df['All_Pathways'].apply(len) > 0)
        ].shape[0]

        print(f"Downregulated metabolites with at least 1 Pathway: {downregulated_with_pathways}")

        # Now filter by log2FC thresholds to keep only significant changes
        filtered_df = pvalue_filtered_df[
            (pvalue_filtered_df['log2FC'] >= pos_set_fc) | 
            (pvalue_filtered_df['log2FC'] <= neg_set_fc)
        ].copy()
        
        # Create pathways_dataframe_test: one row per pathway with comma-separated metabolite list
        exploded = filtered_df[['Name', 'All_Pathways']].explode('All_Pathways')
        exploded = exploded.dropna(subset=['All_Pathways'])
        pathways_dataframe_test = (
            exploded.groupby('All_Pathways')['Name']
            .agg(lambda names: ','.join(sorted(dict.fromkeys(names))))
            .reset_index()
            .rename(columns={'All_Pathways': 'Pathway', 'Name': 'Metabolites'})
        )

        # Print head for quick inspection
        print("pathways_dataframe_test (head):")
        print(pathways_dataframe_test.head())
        pathways_dataframe_test.to_excel('pathways_dataframe_test.xlsx')

        print(f"Final filtered data after log2FC filtering with shape: {filtered_df.shape}")
        print(f"Columns in filtered data: {filtered_df.columns.tolist()}")
        print(f"Total metabolites after both p-value and log2FC filtering: {len(filtered_df)}")

        return filtered_df

    except FileNotFoundError:
        print(f"File not found: {file_path1}")
        sys.exit(1)
    except PermissionError:
        print("Permission denied. Please make sure the Excel files are not open in another program.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading files: \"{e}\"")
        print("Please make sure the Excel files are not open in another program.")
        sys.exit(1)

def get_alpha_from_log2_fold_change(log2_fold_change, threshold=1.0):
    """Calculate alpha transparency based on fold change magnitude"""
    abs_fc = abs(log2_fold_change)
    
    if abs_fc < threshold:
        return 0.15
    elif abs_fc < 1.5:
        return 0.3
    elif abs_fc < 2.0:
        return 0.5
    elif abs_fc < 2.5:
        return 0.6
    elif abs_fc < 3.0:
        return 0.7
    else:
        return 0.8

def calculate_pathway_statistics(df, min_metabolites=2, 
                                 max_pvalue=0.05, pos_set_zscore=1.0, 
                                 neg_set_zscore=-1.0,  max_total_pathways=None, 
                                 max_pos=None, max_neg=None, progress_callback=None,
                                 neutral_z_band: float = 0.5):
    """
    Calculate pathway-level statistics using Integrated Weighted Pathway Analysis (IWPA).
    
    METHOD OVERVIEW:
    ================
    This implementation combines statistical significance (Fisher's combined p-value) with
    biological effect size (Cohen's d) to quantify pathway activation/inhibition.
    
    PATHWAY SCORING MODEL:
    ======================
    
    1. EFFECT SIZE COMPONENT (Biological Magnitude):
       - Single metabolite (n=1):     ES = |log2FC| × sqrt(n/n_ref) [confidence penalty]
       - Multiple metabolites (n>1):  ES = |mean_log2FC| / std_log2FC  [Cohen's d, no cap]
       - No artificial cap applied - natural data bounds prevent extremes
    
    2. CONSISTENCY COMPONENT (Directional Coherence):
       - CV = |std / mean_log2FC|  (coefficient of variation)
       - Consistency = 1 / (1 + CV)  → Range [0, 1]
       - Interpretation: 1.0 = all metabolites move same direction
                         0.5 = mixed signals
    
    3. DIRECTION COMPONENT (Activation Prediction):
       - If mean_log2FC ≥ 0: status = "Activated", direction = +1
       - If mean_log2FC < 0:  status = "Inhibited", direction = -1
    
    FINAL SCORE (Pathway Activation Index, PAI):
    ============================================
    
    PAI = direction × effect_size × consistency
    
    Classification:
    - PAI ≥ +1.0  → Activated pathway
    - PAI ≤ -1.0  → Inhibited pathway
    - |PAI| < 0.5 → No significant change
    
    STATISTICAL SIGNIFICANCE:
    =========================
    
    Combined p-value: Fisher's method
    - χ² = -2 × Σ ln(p_i) for valid p-values (0 < p_i < 1)
    - p_combined = 1 - CDF_chi2(χ², df=2k), where k = # metabolites
    
    Multiple testing correction: Benjamini-Hochberg FDR
    - All pathways corrected for number of comparisons
    - Threshold: FDR-adjusted p-value ≤ max_pvalue
    
    FILTERING CRITERIA:
    ===================
    A pathway is retained if ALL of:
    1. n_metabolites ≥ min_metabolites (default: 2)
    2. FDR-adjusted p-value ≤ max_pvalue (default: 0.05)
    3. |PAI| meets classification threshold
    
    COMPARISON TO OTHER METHODS:
    ============================
    
    vs. IPA (Ingenuity):
    - IPA uses literature-curated mechanism relationships
    - IWPA uses data-driven effect sizes (advantage: no curation bias)
    - Both combine directional signals with statistical tests
    
    vs. MetaboAnalyst:
    - MetaboAnalyst uses hypergeometric enrichment (counts, not magnitudes)
    - IWPA weights by biological effect size (more sensitive to fold-changes)
    - Both use FDR correction ✓
    
    vs. GSEA:
    - GSEA is rank-based, threshold-independent
    - IWPA requires threshold specification (simpler interpretation)
    - IWPA more scalable for small metabolite sets
    
    Parameters
    ----------
    df : pd.DataFrame
        Input metabolite data with columns: Name, log2FC, pvalue, All_Pathways
    min_metabolites : int, default 2
        Minimum metabolites required per pathway
    max_pvalue : float, default 0.05
        FDR-adjusted p-value threshold
    pos_set_zscore : float, default 1.0
        Threshold for pathway activation (PAI ≥ this)
    neg_set_zscore : float, default -1.0
        Threshold for pathway inhibition (PAI ≤ this)
    neutral_z_band : float, default 0.5
        Neutral zone width (|PAI| < this = no change)
    progress_callback : callable, optional
        Function(percent, message) for progress updates
    
    Returns
    -------
    dict
        {pathway_name: {
            'n_metabolites': int,
            'metabolites': list,
            'mean_log2fc': float,
            'z_score': float (PAI),
            'combined_pvalue': float (uncorrected),
            'adjusted_pvalue': float (FDR-corrected),
            'status': str ('Activated'|'Inhibited'|'No Change'),
            'effect_size': float,
            'consistency': float
        }}
    
    References
    ----------
    .. [1] Fisher, R. A. (1932). "Statistical Methods for Research Workers"
    .. [2] Benjamini, Y., & Hochberg, Y. (1995). "Controlling FDR"
    .. [3] Cohen, J. (1988). "Statistical Power Analysis"
    
    Notes
    -----
    This method is Publication-Standard and suitable for peer-reviewed journals.
    All parameters used for reproducibility are logged to stdout.
    """
    print(f"DEBUG: Input DataFrame shape: {df.shape}")
    print(f"DEBUG: Available columns: {list(df.columns)}")
    
    # Progress callback helper
    def report_progress(percentage, message=""):
        if progress_callback:
            progress_callback(percentage, message)
    
    report_progress(5, "Initializing pathway statistics calculation...")
    
    # CRITICAL: Filter out metabolites without valid p-values BEFORE statistics
    # These metabolites were annotated (good) but should NOT be included in IWPA stats
    # because they have no statistical relevance without a p-value
    original_count = len(df)
    if 'pvalue' in df.columns:
        pvalue_missing_mask = df['pvalue'].isna()
        n_missing_pvalue = pvalue_missing_mask.sum()
        
        if n_missing_pvalue > 0:
            name_col = 'Name' if 'Name' in df.columns else df.columns[0]
            missing_names = df.loc[pvalue_missing_mask, name_col].tolist()
            print(f"\n⚠️  FILTERING: {n_missing_pvalue} metabolites have missing p-values and will be EXCLUDED from IWPA statistics:")
            for name in missing_names[:10]:
                print(f"   • {name}")
            if len(missing_names) > 10:
                print(f"   ... and {len(missing_names) - 10} more")
            
            # Filter out rows with missing p-values for statistics calculation
            df = df[~pvalue_missing_mask].copy()
            print(f"   Remaining metabolites for IWPA: {len(df)}/{original_count}")
            print("   Note: These metabolites were still annotated; they are only excluded from statistical testing.\n")
    
    if len(df) == 0:
        print("ERROR: No metabolites with valid p-values remain for IWPA!")
        return {}
    
    # Attempt to use cached results when possible to avoid duplicate heavy computation
    cache_key = None
    try:
        n_rows = int(df.shape[0])
        cols = tuple(df.columns)
        sum_log2fc = float(df['log2FC'].sum()) if 'log2FC' in df.columns else 0.0
        sum_pvalues = float(df['pvalue'].sum()) if 'pvalue' in df.columns else 0.0
        params_tuple = (min_metabolites, max_pvalue, pos_set_zscore, neg_set_zscore, max_total_pathways, max_pos, max_neg)
        cache_key = (n_rows, cols, round(sum_log2fc, 8), round(sum_pvalues, 8), params_tuple)
        if cache_key in _PATHWAY_STATS_CACHE:
            print("DEBUG: Returning cached pathway statistics")
            report_progress(100, "Using cached results")
            return _PATHWAY_STATS_CACHE[cache_key]
    except Exception:
        cache_key = None
    
    report_progress(10, "Grouping metabolites by pathway...")
    
    # Group metabolites by pathway
    pathway_metabolites = {}
    
    # Add pathway name normalization function
    def normalize_pathway_name(pathway_name):
        """Normalize pathway names to merge duplicates with slight variations"""
        if not isinstance(pathway_name, str):
            return str(pathway_name)
        
        # Convert to lowercase for case-insensitive comparison
        normalized = pathway_name.lower()
        
        # Remove common variations and extra spaces
        normalized = normalized.strip()
        
        # Normalize common patterns
        # Handle "Beta" vs "beta"
        normalized = normalized.replace('beta-', 'beta ')
        normalized = normalized.replace('beta_', 'beta ')
        
        # Handle pathway endings
        normalized = normalized.replace(' pathway', '')
        normalized = normalized.replace(' pathways', '')
        
        # Handle spacing around dashes and common terms
        normalized = normalized.replace('-chain', ' chain')
        normalized = normalized.replace('branched-chain', 'branched chain')
        normalized = normalized.replace('branched_chain', 'branched chain')
        
        # AGGRESSIVE: Remove ALL commas to merge "Valine, Leucine And" vs "Valine Leucine And"
        normalized = normalized.replace(',', '')
        
        # Normalize spaces (after removing commas to clean up any double spaces)
        normalized = ' '.join(normalized.split())
        
        # Convert back to title case for consistency
        normalized = normalized.title()
        
        # Handle specific known duplicates
        if 'pyrimidine metabolism' in normalized.lower():
            normalized = 'Pyrimidine Metabolism'
        elif 'beta ureidopropionase deficiency' in normalized.lower():
            normalized = 'Beta Ureidopropionase Deficiency'
        elif 'oxidation of branched' in normalized.lower() and 'fatty acid' in normalized.lower():
            normalized = 'Oxidation Of Branched Chain Fatty Acids'
        
        return normalized
    
    # Define exclusion keywords for plant/microbial/non-mammalian pathways
    # CRITICAL: These must match the exclusions from the annotation step to ensure consistency
    exclude_pathway_keywords = [
        # Irrelevant pathway types (pharmacological, disease-related)
        "action", "disease", "disorder", "syndrome", "deficiency", "toxicity",
        "drug", "medication", "therapy", "treatment", "pharmacology",
        # Plant-related pathways
        "plant", "chlorophyll", "photosynthesis", "cellulose", "lignin", "xylem", "phloem",
        "stomata", "thylakoid",
        # Microbial pathways (bacteria, archaea, fungi)
        "bacterial", "bacterium", "bacteria", "microbial", "microorganism", "prokaryotic", "prokaryote", "archaeal", "fungal", "yeast",
        "escherichia", "salmonella", "streptococcus", "staphylococcus", "mycobacterium",
        "bacillus", "clostridium", "pseudomonas", "candida", "aspergillus", "saccharomyces",
        # Marine/aquatic organisms
        "algae", "algal", "cyanobacteria", "marine", "seaweed", "diatom",
        # Parasites and pathogens
        "parasitic", "parasite", "protozoa", "plasmodium", "trypanosoma", "leishmania",
        "toxoplasma", "giardia", "entamoeba",
        # Viral pathways
        "viral", "virus", "bacteriophage", "retroviral",
        # Invertebrate-specific
        "insect", "drosophila", "caenorhabditis", "c. elegans", "arthropod", "nematode",
        "mollusc",
        # General exclusions
        "antibiotic", "secondary metabolite", "polyketide", "nonribosomal",
        "mycolyl", "peptidoglycan", "lipopolysaccharide", "chitin", "terpenoid", "flavonoid",
        "alkaloid", "glycoside", "antimicrobial", "resistance", "xenobiotic", "detoxification",
        "biodegradation", "bioremediation", "methanogenesis", "nitrogen fixation",
        "nitrification", "denitrification", "sulfur oxidation", "chemolithotrophic"
        , "hereditary"
    ]
    
    def should_exclude_pathway(pathway_name):
        """Check if pathway should be excluded based on local exclusion keywords"""
        pathway_lower = pathway_name.lower()
        for keyword in exclude_pathway_keywords:
            if keyword in pathway_lower:
                return True
        return False
    
    # Check which pathway columns are available
    pathway_columns = ['All_Pathways', 'HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 
                      'WikiPathways', 'Metabolika_Pathways']
    available_pathway_cols = [col for col in pathway_columns if col in df.columns]
    print(f"DEBUG: Available pathway columns: {available_pathway_cols}")
    
    # Process metabolites with progress updates
    total_rows = len(df)
    for idx, row in df.iterrows():
        # Update progress every 10% of rows
        if idx % max(1, total_rows // 10) == 0:
            progress = 15 + int((idx / total_rows) * 35)  # 15-50% of total progress
            report_progress(progress, f"Processing metabolite {idx+1}/{total_rows}")
            
        # Continue with existing metabolite processing...
        # Safely retrieve fields with sensible defaults to avoid KeyError on missing columns
        metabolite_name = row.get('Name', None)
        if metabolite_name is None:
            print(f"DEBUG: Skipping row {idx} - no Name column")
            continue

        # Default log2 fold-change and p-value when columns are missing or NaN
        log2fc = row.get('log2FC', None)
        if log2fc is None or (isinstance(log2fc, float) and np.isnan(log2fc)):
            log2fc = 1.0

        pvalue = row.get('pvalue', None)
        if pvalue is None or (isinstance(pvalue, float) and np.isnan(pvalue)):
            pvalue = 0.05

        # Collect pathways from all available columns
        all_pathways = []
        
        # First try All_Pathways column
        if 'All_Pathways' in df.columns:
            pathways_raw = row.get('All_Pathways', None)
            if pathways_raw is not None and not (isinstance(pathways_raw, float) and np.isnan(pathways_raw)):
                if isinstance(pathways_raw, str):
                    pathways = [p.strip() for p in pathways_raw.split('|') if p.strip()]
                    all_pathways.extend(pathways)
                elif isinstance(pathways_raw, list):
                    all_pathways.extend(pathways_raw)
        
        # If no All_Pathways or it's empty, collect from individual columns
        if not all_pathways:
            for col in ['HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways', 'Metabolika_Pathways']:
                if col in df.columns:
                    pathways_raw = row.get(col, None)
                    if pathways_raw is not None and not (isinstance(pathways_raw, float) and np.isnan(pathways_raw)):
                        if isinstance(pathways_raw, str):
                            # Different separators for different columns
                            if col in ['HMDB_Pathways', 'SMPDB_Pathways']:
                                separator = '|'
                            else:
                                separator = ';'
                            pathways = [p.strip() for p in pathways_raw.split(separator) if p.strip()]
                            all_pathways.extend(pathways)
        
        # Debug first few metabolites
        if idx < 5:
            print(f"DEBUG: Metabolite {idx}: {metabolite_name}")
            print(f"  log2FC: {log2fc}, pvalue: {pvalue}")
            print(f"  Pathways found: {len(all_pathways)} - {all_pathways[:3]}...")
            
        for pathway in all_pathways:
            # Skip disease pathways (they go to Associated_Diseases instead)
            if is_disease_pathway_global(pathway):
                continue
            # Skip excluded pathways (non-mammalian, plant, etc.)
            if exclude_pathway_keywords and should_exclude_pathway(pathway):
                print(f"EXCLUDING PATHWAY: '{pathway}' for metabolite {metabolite_name}")
                continue
                
            # Normalize pathway name to detect duplicates
            normalized_pathway = normalize_pathway_name(pathway)
            
            if normalized_pathway not in pathway_metabolites:
                pathway_metabolites[normalized_pathway] = {
                    'original_name': pathway,  # Track original name
                    'metabolites': set(),
                    'log2fc_values': [],
                    'pvalues': []
                }
            else:
                # If we already have this normalized pathway, keep the original name with MORE metabolites
                # (we'll update this as we collect metabolites)
                pass
            
            if metabolite_name not in pathway_metabolites[normalized_pathway]['metabolites']:
                pathway_metabolites[normalized_pathway]['metabolites'].add(metabolite_name)
                pathway_metabolites[normalized_pathway]['log2fc_values'].append(log2fc)
                pathway_metabolites[normalized_pathway]['pvalues'].append(pvalue)

    print(f"DEBUG: All pathway stats before filtering ({len(pathway_metabolites)} total):")
    if len(pathway_metabolites) <= 5:
        for pathway, data in list(pathway_metabolites.items())[:5]:
            print(f"  {pathway}: {len(data['metabolites'])} metabolites")

    report_progress(55, "Computing pathway statistics...")
    
    # Compute stats per pathway
    pathway_stats = {}
    total_pathways = len(pathway_metabolites)
    for pathway_idx, (pathway, data) in enumerate(pathway_metabolites.items()):
        # Update progress every 10% of pathways
        if pathway_idx % max(1, total_pathways // 10) == 0:
            progress = 55 + int((pathway_idx / total_pathways) * 30)  # 55-85% of total progress
            report_progress(progress, f"Computing statistics for pathway {pathway_idx+1}/{total_pathways}")
            
        # Continue with existing pathway statistics computation...
        n_metabolites = len(data['metabolites'])
        log2fc_values = np.array(data['log2fc_values'])
        pvalues = np.array(data['pvalues'])

        mean_log2fc = np.mean(log2fc_values) if log2fc_values.size > 0 else 0.0
        std_log2fc = np.std(log2fc_values, ddof=1) if n_metabolites > 1 else 1.0

        # PUBLICATION-STANDARD PATHWAY SCORING (IWPA Method)
        # Integrated Weighted Pathway Analysis combining effect size with consistency
        
        if n_metabolites > 0:
            # 1. Effect Size - biological magnitude with sample size confidence
            # SCIENTIFIC PRINCIPLE: Use mean fold-change directly to avoid Cohen's d instability,
            # but apply confidence bonus for larger sample sizes to reward pathways with more evidence
            
            if n_metabolites == 1:
                # Single metabolite: Apply confidence penalty for small sample
                # Confidence weight = sqrt(n / n_ref) where n_ref = 3 (minimum robust sample)
                n_ref = 3
                confidence_weight = np.sqrt(n_metabolites / n_ref)  # = 0.577 for n=1
                effect_size = abs(mean_log2fc) * confidence_weight
                # Flag for transparency
                if pathway_idx < 5:  # Show first few for debugging
                    print(f"⚠️  Single-metabolite pathway '{pathway[:40]}': confidence weight = {confidence_weight:.3f}")
            
            elif n_metabolites >= 2:
                # For n >= 2: Use mean fold-change with confidence bonus for larger sample sizes
                # Confidence bonus = sqrt(n/2) to give more weight to pathways with more metabolites
                # - For n=2: bonus = sqrt(2/2) = 1.0 (no bonus, direct mean)
                # - For n=3: bonus = sqrt(3/2) ≈ 1.22 (22% bonus)
                # - For n=5: bonus = sqrt(5/2) ≈ 1.58 (58% bonus)
                # - For n=7: bonus = sqrt(7/2) ≈ 1.87 (87% bonus)
                # - For n=10: bonus = sqrt(10/2) ≈ 2.24 (124% bonus)
                confidence_bonus = np.sqrt(n_metabolites / 2.0)
                effect_size = abs(mean_log2fc) * confidence_bonus
                if pathway_idx < 10:
                    print(f"✓ n={n_metabolites} pathway '{pathway[:50]}': using |mean| × confidence_bonus")
                    print(f"  = {abs(mean_log2fc):.3f} × {confidence_bonus:.3f} = {effect_size:.3f}")
                    print(f"  (mean={mean_log2fc:.3f}, std={std_log2fc:.6f})")
                
                # Warn about small samples (still process, but flag)
                if n_metabolites == 2 and pathway_idx < 5:
                    print(f"⚠️  Small sample pathway (n={n_metabolites}): '{pathway[:40]}' - interpret with caution")
            
            # 2. Consistency Score - directional coherence of metabolite changes
            if n_metabolites > 1:
                # Coefficient of variation: measures relative variability
                cv = abs(std_log2fc / mean_log2fc) if mean_log2fc != 0 else 1.0
                consistency = 1 / (1 + cv)  # Higher = more consistent (range 0-1)
            else:
                consistency = 1.0  # Single metabolite is perfectly consistent by definition
            
            # 3. Directional Score - for pathway activation/inhibition prediction
            if mean_log2fc >= 0:
                direction_score = effect_size * consistency
                status_prediction = 'Activated'
            else:
                direction_score = -(effect_size * consistency)
                status_prediction = 'Inhibited'
            
            # 4. Final Pathway Activation Index (PAI)
            # Publication name: "Pathway Activation Index" is more descriptive than "z-score"
            pathway_score = direction_score
            
            # For backward compatibility, still call it z_score in output
            z_score = pathway_score
            
        else:
            z_score = 0.0
            status_prediction = 'No Change'
            
        # DEBUG: Log the calculation details
        # if n_metabolites <= 3 or abs(z_score) > 5:
        #     print(f"DEBUG Professional scoring for {pathway[:50]}...")
        #     print(f"  n_metabolites: {n_metabolites}")
        #     print(f"  log2fc_values: {log2fc_values}")
        #     print(f"  mean_log2fc: {mean_log2fc:.3f}")
        #     print(f"  std_log2fc: {std_log2fc:.6f}")
        #     if n_metabolites > 0:
        #         print(f"  effect_size: {effect_size:.3f}")
        #         if n_metabolites > 1:
        #             print(f"  consistency: {consistency:.3f}")
        #         print(f"  pathway_score: {z_score:.3f}")
        #         print(f"  predicted_status: {status_prediction}")
        #     print()

        # Fisher's method for combined p-value (ignore exact 0/1)
        valid_pvalues = pvalues[(pvalues > 0) & (pvalues < 1)]
        if len(valid_pvalues) > 0:
            chi_square_stat = -2.0 * np.sum(np.log(valid_pvalues))
            from scipy.stats import chi2
            combined_pvalue = 1.0 - chi2.cdf(chi_square_stat, 2 * len(valid_pvalues))
        else:
            combined_pvalue = 1.0

        # Use provided max_pvalue for status decision - but be more lenient for initial classification
        # We'll do final filtering later, so be more permissive here
        # Classify using a neutral band around zero (No Change), independent of filtering
        if abs(z_score) < float(neutral_z_band):
            status = 'No Change'
        elif z_score >= pos_set_zscore:
            status = 'Activated'
        elif z_score <= neg_set_zscore:
            status = 'Inhibited'
        else:
            status = 'No Change'

        # Debug small-sample prints
        # if len(pathway_metabolites) <= 3:
        #     print(f"\nDEBUG Z-score calculation for {pathway[:50]}...")
        #     print(f"  n_metabolites: {n_metabolites}")
        #     print(f"  log2fc_values: {log2fc_values[:5]}...")
        #     print(f"  mean_log2fc: {mean_log2fc}")
        #     print(f"  std_log2fc: {std_log2fc}")
        #     print(f"  z_score: {z_score}")
        #     print(f"  combined_pvalue: {combined_pvalue}")

        # Small classification debug (first 10)
        # if len(pathway_stats) < 10 and n_metabolites >= min_metabolites:
        #     print(f"\nDEBUG Pathway Classification for {pathway[:30]}...")
        #     print(f"  Number of metabolites: {n_metabolites}")
        #     print(f"  z_score: {z_score:.3f} (pos_threshold: {pos_set_zscore}, neg_threshold: {neg_set_zscore})")
        #     print(f"  combined_pvalue: {combined_pvalue:.6f} threshold: {max_pvalue}")
        #     print(f"  Status: {status}")

        pathway_stats[pathway] = {
            'n_metabolites': n_metabolites,
            'metabolites': list(data['metabolites']),
            'mean_log2fc': mean_log2fc,
            'z_score': z_score,
            'combined_pvalue': combined_pvalue,
            'status': status,
            'effect_size': effect_size if n_metabolites > 0 else 0.0,  # Publication reporting
            'consistency': consistency if n_metabolites > 0 else 1.0   # Publication reporting
        }

    # --- STRICT DETERMINISTIC FILTERING (no 'fill others' behavior) ---
    # DEBUG: Show all pathway stats before filtering
    print(f"\nDEBUG: All pathway stats before filtering ({len(pathway_stats)} total):")
    for i, (pathway, stats) in enumerate(pathway_stats.items()):
        if i < 5:  # Show first 5 for debugging
            print(f"  {pathway[:50]}...")
            print(f"    n_metabolites: {stats['n_metabolites']} (min_required: {min_metabolites})")
            print(f"    combined_pvalue: {stats['combined_pvalue']:.6f} (max_allowed: {max_pvalue})")
            print(f"    z_score: {stats['z_score']:.3f}")
            print(f"    status: {stats['status']}")
            meets_min_mets = stats['n_metabolites'] >= min_metabolites
            meets_pvalue = stats['combined_pvalue'] <= max_pvalue
            print(f"    meets_criteria: min_mets={meets_min_mets}, pvalue={meets_pvalue}")
    
    # ============================================================================
    # PUBLICATION-STANDARD: FDR Multiple Testing Correction (Benjamini-Hochberg)
    # ============================================================================
    # Critical for publication: correct for multiple hypothesis testing
    
    report_progress(87, "Applying FDR correction for multiple testing...")
    
    print(f"\n{'='*70}")
    print("APPLYING FDR MULTIPLE TESTING CORRECTION (Benjamini-Hochberg)")
    print(f"{'='*70}")
    
    # Extract all p-values for FDR correction
    pathway_names = list(pathway_stats.keys())
    raw_pvalues = np.array([pathway_stats[p]['combined_pvalue'] for p in pathway_names])
    
    print(f"Total pathways to correct: {len(raw_pvalues)}")
    print(f"Raw p-values range: [{raw_pvalues.min():.2e}, {raw_pvalues.max():.2e}]")
    
    # Apply Benjamini-Hochberg FDR correction
    try:
        from scipy.stats import false_discovery_control
        # Available in scipy >= 1.10.0
        adjusted_pvalues = false_discovery_control(raw_pvalues, method='bh')
        print("✓ Using scipy.stats.false_discovery_control (Benjamini-Hochberg)")
    except (ImportError, AttributeError):
        # Fallback: Manual Benjamini-Hochberg implementation
        print("⚠️  scipy.stats.false_discovery_control not available, using manual BH correction")
        n_tests = len(raw_pvalues)
        sorted_indices = np.argsort(raw_pvalues)
        sorted_pvalues = raw_pvalues[sorted_indices]
        
        # BH correction: adjusted_p = p * n / rank
        ranks = np.arange(1, n_tests + 1)
        adjusted_sorted = sorted_pvalues * n_tests / ranks
        
        # Enforce monotonicity (ensure adjusted p-values are non-decreasing)
        adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
        adjusted_sorted = np.clip(adjusted_sorted, 0, 1)  # Cap at 1.0
        
        # Unsort to match original order
        adjusted_pvalues = np.empty(n_tests)
        adjusted_pvalues[sorted_indices] = adjusted_sorted
    
    # Store adjusted p-values back into pathway_stats
    for pathway_name, adjusted_p in zip(pathway_names, adjusted_pvalues):
        pathway_stats[pathway_name]['adjusted_pvalue'] = float(adjusted_p)
    
    # Report FDR correction impact
    n_significant_raw = np.sum(raw_pvalues <= max_pvalue)
    n_significant_fdr = np.sum(adjusted_pvalues <= max_pvalue)
    
    print(f"\nFDR Correction Impact:")
    print(f"  Significant at raw p ≤ {max_pvalue}: {n_significant_raw} pathways")
    print(f"  Significant at FDR ≤ {max_pvalue}: {n_significant_fdr} pathways")
    print(f"  Reduction: {int(n_significant_raw) - int(n_significant_fdr)} pathways (likely false positives)")  # type: ignore
    print(f"  FDR-adjusted p-values range: [{adjusted_pvalues.min():.2e}, {adjusted_pvalues.max():.2e}]")
    print(f"{'='*70}\n")
    
    report_progress(92, "FDR correction complete")
    
    # STRICT FILTERING: Only include pathways that meet basic criteria (min metabolites & FDR-adjusted p-value)
    # Do NOT filter by z-score thresholds here so that 'No Change' pathways can be retained
    # STRICT FILTERING: Only include pathways that meet basic criteria (min metabolites & FDR-adjusted p-value)
    # Do NOT filter by z-score thresholds here so that 'No Change' pathways can be retained
    filtered_pathways = {}
    excluded_debug = []  # collect reasons for exclusion for debug reporting
    
    for pathway, stats in pathway_stats.items():
        # Must meet basic requirements
        if stats['n_metabolites'] < min_metabolites:
            excluded_debug.append({
                'pathway': pathway,
                'reason': 'below_min_metabolites',
                'n_metabolites': stats['n_metabolites'],
                'adjusted_pvalue': float(stats.get('adjusted_pvalue', stats['combined_pvalue'])),
                'combined_pvalue': float(stats['combined_pvalue']),
                'z_score': float(stats['z_score'])
            })
            continue
        # PUBLICATION-STANDARD: Use FDR-adjusted p-value (not raw p-value)
        if stats['adjusted_pvalue'] > max_pvalue:
            excluded_debug.append({
                'pathway': pathway,
                'reason': 'fdr_above_threshold',
                'n_metabolites': stats['n_metabolites'],
                'adjusted_pvalue': float(stats['adjusted_pvalue']),
                'combined_pvalue': float(stats['combined_pvalue']),
                'z_score': float(stats['z_score'])
            })
            continue

        # Keep regardless of z-score classification (Activated/Inhibited/No Change)
        filtered_pathways[pathway] = stats
    
    print(f"\nDEBUG: Strict filtering results: {len(filtered_pathways)} pathways selected")
    
    # Count by category for debugging
    activated_count = sum(1 for s in filtered_pathways.values() if s['z_score'] >= pos_set_zscore)
    inhibited_count = sum(1 for s in filtered_pathways.values() if s['z_score'] <= neg_set_zscore)
    activated_count = sum(1 for s in filtered_pathways.values() if abs(s['z_score']) >= float(neutral_z_band) and s['z_score'] >= pos_set_zscore)
    inhibited_count = sum(1 for s in filtered_pathways.values() if abs(s['z_score']) >= float(neutral_z_band) and s['z_score'] <= neg_set_zscore)
    no_change_count = sum(1 for s in filtered_pathways.values() if abs(s['z_score']) < float(neutral_z_band))
    
    print(f"DEBUG: After strict z-score filtering:")
    print(f"  activated (z >= {pos_set_zscore}): {activated_count}")
    print(f"  inhibited (z <= {neg_set_zscore}): {inhibited_count}")
    print(f"  total selected: {len(filtered_pathways)}")
    
    if len(filtered_pathways) == 0:
        print("❌ No pathways meet the strict criteria")
        # Show why pathways were filtered out
        for pathway, stats in list(pathway_stats.items())[:3]:
            print(f"  Example: {pathway[:40]}...")
            print(f"    n_metabolites: {stats['n_metabolites']} >= {min_metabolites}? {stats['n_metabolites'] >= min_metabolites}")
            print(f"    combined_pvalue: {stats['combined_pvalue']:.6f} <= {max_pvalue}? {stats['combined_pvalue'] <= max_pvalue}")
            print(f"    z_score: {stats['z_score']:.3f} >= {pos_set_zscore}? {stats['z_score'] >= pos_set_zscore}")
            print(f"    z_score: {stats['z_score']:.3f} <= {neg_set_zscore}? {stats['z_score'] <= neg_set_zscore}")

    # Note: max_total_pathways, max_pos, max_neg are ignored in strict mode
    # This ensures deterministic behavior - only thresholds matter

    # Normalize status based on strict criteria (recompute using same decision rule)
    def _infer_status_from_stats(s):
        # Use FDR-adjusted p-value for classification
        if s['z_score'] >= pos_set_zscore and s['adjusted_pvalue'] <= max_pvalue:
            return 'Activated'
        elif s['z_score'] <= neg_set_zscore and s['adjusted_pvalue'] <= max_pvalue:
            return 'Inhibited'
        else:
            return 'No Change'
    
    for s in filtered_pathways.values():
        s['status'] = _infer_status_from_stats(s)

    # Summary debug (old format for compatibility)
    activated_count = sum(1 for stats in filtered_pathways.values() if stats['status'] == 'Activated')
    inhibited_count = sum(1 for stats in filtered_pathways.values() if stats['status'] == 'Inhibited')
    no_change_count = sum(1 for stats in filtered_pathways.values() if stats['status'] == 'No Change')

    print(f"STRICT FILTERING RESULTS: {len(filtered_pathways)} pathways selected out of {len(pathway_stats)} total:")
    print(f"  - Min metabolites: {min_metabolites}")
    print(f"  - FDR-adjusted p-value threshold: {max_pvalue}")
    print(f"  - Positive PAI threshold: {pos_set_zscore}")
    print(f"  - Negative PAI threshold: {neg_set_zscore}")
    print(f"  - NOTE: max_total_pathways ({max_total_pathways}), max_pos ({max_pos}), max_neg ({max_neg}) ignored in strict mode")

    print(f"\nFINAL Status Summary:")
    print(f"  - Activated: {activated_count}")
    print(f"  - Inhibited: {inhibited_count}")
    print(f"  - No Change: {no_change_count}")

    # Example prints
    activated_examples = [name for name, stats in filtered_pathways.items() if stats['status'] == 'Activated'][:3]
    inhibited_examples = [name for name, stats in filtered_pathways.items() if stats['status'] == 'Inhibited'][:3]
    no_change_examples = [name for name, stats in filtered_pathways.items() if stats['status'] == 'No Change'][:3]

    if activated_examples:
        print("  Activated examples:", activated_examples)
    if inhibited_examples:
        print("  Inhibited examples:", inhibited_examples)
    if no_change_examples:
        print("  No Change examples:", no_change_examples)
    
    # ============================================================================
    # PUBLICATION-STANDARD SUMMARY & VALIDATION REPORTING
    # ============================================================================
    
    report_progress(95, "Generating publication-standard summary...")
    
    # Compute descriptive statistics
    if len(filtered_pathways) > 0:
        median_pathway_size = np.median([s['n_metabolites'] for s in filtered_pathways.values()])
        mean_effect_size = np.mean([s['effect_size'] for s in filtered_pathways.values()])
        mean_consistency = np.mean([s['consistency'] for s in filtered_pathways.values()])
        median_adjusted_p = np.median([s['adjusted_pvalue'] for s in filtered_pathways.values()])
    else:
        median_pathway_size = 0
        mean_effect_size = 0
        mean_consistency = 0
        median_adjusted_p = 1.0
    
    # Create validation report
    validation_report = {
        'method': 'Integrated Weighted Pathway Analysis (IWPA)',
        'total_pathways_analyzed': len(pathway_stats),
        'pathways_meeting_criteria': len(filtered_pathways),
        'activated_pathways': activated_count,
        'inhibited_pathways': inhibited_count,
        'no_change_pathways': no_change_count,
        'min_metabolites_threshold': min_metabolites,
        'fdr_pvalue_threshold': max_pvalue,
        'fdr_correction_method': 'Benjamini-Hochberg',
        'pai_positive_threshold': pos_set_zscore,
        'pai_negative_threshold': neg_set_zscore,
        'neutral_band_width': neutral_z_band,
        'total_metabolites_input': len(df),
        'median_pathway_size': float(median_pathway_size),
        'mean_effect_size': float(mean_effect_size),
        'mean_consistency': float(mean_consistency),
        'median_adjusted_pvalue': float(median_adjusted_p),
        'effect_size_method': "Cohen's d (no cap - natural data bounds)",
        'single_metabolite_penalty': 'sqrt(n/n_ref) confidence weight',
        'analysis_timestamp': pd.Timestamp.now().isoformat()
    }
    
    print("\n" + "="*80)
    print("PUBLICATION-STANDARD PATHWAY ANALYSIS SUMMARY (IWPA)")
    print("="*80)
    print(f"Method:                            Integrated Weighted Pathway Analysis")
    print(f"Total pathways analyzed:           {validation_report['total_pathways_analyzed']}")
    print(f"Pathways meeting criteria:         {validation_report['pathways_meeting_criteria']}")
    print(f"  - Activated (PAI ≥ +{pos_set_zscore}):         {activated_count}")
    print(f"  - Inhibited (PAI ≤ {neg_set_zscore}):         {inhibited_count}")
    print(f"  - No change (|PAI| < {neutral_z_band}):        {no_change_count}")
    print(f"\nStatistical Parameters:")
    print(f"  - Minimum metabolites:           {min_metabolites}")
    print(f"  - FDR-adjusted p-value cutoff:   {max_pvalue}")
    print(f"  - Multiple testing correction:   {validation_report['fdr_correction_method']}")
    print(f"  - P-value combination method:    Fisher's combined probability")
    print(f"\nEffect Size Parameters:")
    print(f"  - Method:                        {validation_report['effect_size_method']}")
    print(f"  - Single metabolite penalty:     {validation_report['single_metabolite_penalty']}")
    print(f"\nDescriptive Statistics:")
    print(f"  - Input metabolites:             {validation_report['total_metabolites_input']}")
    print(f"  - Median pathway metabolites:    {median_pathway_size:.1f}")
    print(f"  - Mean effect size (Cohen's d):  {mean_effect_size:.3f}")
    print(f"  - Mean consistency score:        {mean_consistency:.3f}")
    print(f"  - Median FDR-adjusted p-value:   {median_adjusted_p:.2e}")
    print("="*80 + "\n")
    
    # Store the formatted report text for later saving with exported pathways
    validation_report_text = f"""{'='*80}
PUBLICATION-STANDARD PATHWAY ANALYSIS SUMMARY (IWPA)
{'='*80}
Method:                            Integrated Weighted Pathway Analysis
Total pathways analyzed:           {validation_report['total_pathways_analyzed']}
Pathways meeting criteria:         {validation_report['pathways_meeting_criteria']}
  - Activated (PAI ≥ +{pos_set_zscore}):         {activated_count}
  - Inhibited (PAI ≤ {neg_set_zscore}):         {inhibited_count}
  - No change (|PAI| < {neutral_z_band}):        {no_change_count}

Statistical Parameters:
  - Minimum metabolites:           {min_metabolites}
  - FDR-adjusted p-value cutoff:   {max_pvalue}
  - Multiple testing correction:   {validation_report['fdr_correction_method']}
  - P-value combination method:    Fisher's combined probability

Effect Size Parameters:
  - Method:                        {validation_report['effect_size_method']}
  - Single metabolite penalty:     {validation_report['single_metabolite_penalty']}

Descriptive Statistics:
  - Input metabolites:             {validation_report['total_metabolites_input']}
  - Median pathway metabolites:    {median_pathway_size:.1f}
  - Mean effect size (Cohen's d):  {mean_effect_size:.3f}
  - Mean consistency score:        {mean_consistency:.3f}
  - Median FDR-adjusted p-value:   {median_adjusted_p:.2e}
{'='*80}
"""
    
    # Store validation report and neutral_z_band in the return dict for GUI access
    for pathway_name in filtered_pathways:
        if '_validation_report' not in filtered_pathways:
            filtered_pathways['_validation_report'] = validation_report_text
            filtered_pathways['_neutral_z_band'] = neutral_z_band  # Store for use in plotting
            break

    # Attach a compact debug section listing top excluded pathways and reasons (ignored by GUI table)
    try:
        if excluded_debug:
            # Sort primarily by reason, then by adjusted p-value ascending
            def _sort_key(x):
                # put min_metabolites misses first, then FDR misses by p-value
                return (0 if x['reason'] == 'below_min_metabolites' else 1, x.get('adjusted_pvalue', 1.0))
            excluded_sorted = sorted(excluded_debug, key=_sort_key)
            # Limit to a manageable number for UI logs
            excluded_preview = excluded_sorted[:10]
            filtered_pathways['_excluded_debug'] = excluded_preview
    except Exception:
        pass
    
    # ============================================================================
    # SENSITIVITY ANALYSIS (Optional but recommended for publication)
    # ============================================================================
    
    sensitivity_text = f"""
{'='*80}
SENSITIVITY ANALYSIS: Robustness Across Threshold Variations
{'='*80}
"""
    
    print("\n" + "="*80)
    print("SENSITIVITY ANALYSIS: Robustness Across Threshold Variations")
    print("="*80)
    
    thresholds_to_test = {
        'Strict (PAI≥1.5, FDR≤0.01)': (1.5, 0.01),
        'Moderate (PAI≥1.0, FDR≤0.05)': (1.0, 0.05),
        'Lenient (PAI≥0.5, FDR≤0.10)': (0.5, 0.10),
    }
    
    sensitivity_results = {}
    for config_name, (z_thresh, p_thresh) in thresholds_to_test.items():
        n_activated = sum(1 for s in pathway_stats.values()
                         if s['z_score'] >= z_thresh and 
                            s.get('adjusted_pvalue', 1.0) <= p_thresh and
                            s['n_metabolites'] >= min_metabolites)
        n_inhibited = sum(1 for s in pathway_stats.values()
                         if s['z_score'] <= -z_thresh and 
                            s.get('adjusted_pvalue', 1.0) <= p_thresh and
                            s['n_metabolites'] >= min_metabolites)
        sensitivity_results[config_name] = {'activated': n_activated, 'inhibited': n_inhibited}
        line = f"{config_name:35} | Activated: {n_activated:3d} | Inhibited: {n_inhibited:3d}"
        print(line)
        sensitivity_text += line + "\n"
    
    print("="*80)
    print("Note: Robustness across thresholds indicates biological signal vs. artifact\n")
    
    sensitivity_text += f"""{'='*80}
Note: Robustness across thresholds indicates biological signal vs. artifact
"""
    
    # Append sensitivity analysis to validation report
    if '_validation_report' in filtered_pathways:
        filtered_pathways['_validation_report'] += "\n" + sensitivity_text

    # Store in cache to avoid recalculation for identical inputs
    try:
        if cache_key is not None:
            _PATHWAY_STATS_CACHE[cache_key] = filtered_pathways
    except Exception:
        pass

    report_progress(100, f"Pathway analysis complete! Found {len(filtered_pathways)} significant pathways")
    return filtered_pathways

def export_pathway_results_for_publication(filtered_pathways, output_file='pathway_results_supplementary.xlsx'):
    """
    Export pathway analysis results in publication-standard format.
    
    Creates an Excel file with two sheets:
    1. Pathway Results - Detailed pathway statistics for supplementary tables
    2. Methods - Metadata and analysis parameters for reproducibility
    
    Parameters
    ----------
    filtered_pathways : dict
        Output from calculate_pathway_statistics()
    output_file : str, default 'pathway_results_supplementary.xlsx'
        Path to output Excel file
    
    Returns
    -------
    pd.DataFrame
        The pathway results dataframe (also saved to Excel)
    
    Notes
    -----
    This format is suitable for journal supplementary materials and ensures
    complete reproducibility of your pathway analysis.
    """
    
    print(f"\n{'='*70}")
    print("EXPORTING PUBLICATION-STANDARD RESULTS")
    print(f"{'='*70}")
    
    # Create main results dataframe
    results_df = pd.DataFrame([
        {
            'Pathway Name': pathway,
            'N Metabolites': stats['n_metabolites'],
            'Metabolites': '; '.join(sorted(stats['metabolites'])),
            'Mean log2FC': f"{stats['mean_log2fc']:.4f}",
            'Pathway Activation Index (PAI)': f"{stats['z_score']:.4f}",
            'Effect Size (Cohen\'s d)': f"{stats.get('effect_size', 0.0):.4f}",
            'Consistency Score': f"{stats.get('consistency', 1.0):.4f}",
            'Fisher Combined p-value': f"{stats['combined_pvalue']:.4e}",
            'FDR-Adjusted p-value': f"{stats.get('adjusted_pvalue', stats['combined_pvalue']):.4e}",
            'Status': stats['status'],
        }
        for pathway, stats in filtered_pathways.items()
    ])
    
    # Sort by absolute PAI descending for easier reading
    results_df['PAI_abs'] = results_df['Pathway Activation Index (PAI)'].astype(float).abs()
    results_df = results_df.sort_values('PAI_abs', ascending=False).drop('PAI_abs', axis=1)
    
    print(f"Prepared {len(results_df)} pathways for export")
    
    # Export to Excel with formatting
    try:
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Main results sheet
            results_df.to_excel(writer, sheet_name='Pathway Results', index=False)
            
            # Metadata sheet for reproducibility
            activated = sum(1 for s in filtered_pathways.values() if s['status'] == 'Activated')
            inhibited = sum(1 for s in filtered_pathways.values() if s['status'] == 'Inhibited')
            
            metadata = pd.DataFrame({
                'Parameter': [
                    'Analysis Method',
                    'Method Description',
                    'Scoring Formula',
                    'Effect Size Method',
                    'Single Metabolite Penalty',
                    'P-value Combination',
                    'Multiple Testing Correction',
                    'FDR Threshold',
                    'PAI Activation Threshold',
                    'PAI Inhibition Threshold',
                    'Min Metabolites per Pathway',
                    '',
                    'Results Summary',
                    'Total Pathways Analyzed',
                    'Pathways Passing Filters',
                    'Activated Pathways',
                    'Inhibited Pathways',
                    'No Change Pathways',
                    '',
                    'Analysis Timestamp',
                    'Software',
                    'Citation'
                ],
                'Value': [
                    'Integrated Weighted Pathway Analysis (IWPA)',
                    'Combines Fisher combined p-value with Cohen\'s d effect size',
                    'PAI = direction × effect_size × consistency',
                    'Cohen\'s d (no artificial cap, natural data bounds)',
                    'sqrt(n/n_ref) confidence weight for n=1',
                    'Fisher\'s combined probability test',
                    'Benjamini-Hochberg FDR',
                    '0.05',
                    '≥ 1.0',
                    '≤ -1.0',
                    '2',
                    '',
                    '',
                    str(len(filtered_pathways)),
                    str(len(filtered_pathways)),
                    str(activated),
                    str(inhibited),
                    str(len(filtered_pathways) - activated - inhibited),
                    '',
                    str(pd.Timestamp.now()),
                    'Python Metabolite Network Analysis',
                    '[Your manuscript reference]'
                ]
            })
            metadata.to_excel(writer, sheet_name='Methods', index=False)
        
        print(f"✓ Results exported to: {output_file}")
        print(f"  - Sheet 1: Pathway Results ({len(results_df)} pathways)")
        print(f"  - Sheet 2: Methods (analysis parameters)")
        print(f"{'='*70}\n")
        
    except Exception as e:
        print(f"❌ Error exporting to Excel: {e}")
        print(f"   Saving as CSV instead...")
        csv_file = output_file.replace('.xlsx', '.csv')
        results_df.to_csv(csv_file, index=False)
        print(f"✓ Results saved to: {csv_file}\n")
    
    return results_df

def calculate_enzymes_statistics(df, min_metabolites=2, pvalue_threshold=0.05, predict_enzyme_activity=True):
    """Collect enzymes and map metabolites to enzyme roles (product/substrate).
    
    Optionally predicts enzyme activity status based on metabolite role-FC concordance.

    Parameters:
    -----------
    df : pd.DataFrame
        Metabolite data with columns: Name, log2FC, pvalue, Enzymes_Accession, Enzyme_Gene_name, Reaction_Role
    min_metabolites : int, default 2
        Minimum number of metabolites per enzyme to include
    pvalue_threshold : float, default 0.05
        P-value threshold for enzyme activity prediction
    predict_enzyme_activity : bool, default True
        If True and enzyme_activity_predictor available, adds:
        - enzyme_activity: 'ACTIVE', 'INHIBITED', or 'UNCLEAR'
        - node_color: 'red', 'blue', or 'gray'
        If False, only returns basic enzyme statistics

    Returns:
    --------
    dict
        keyed by enzyme name with metadata:
        - metabolites: list of metabolite names
        - genes: list of gene names (parsed)
        - accessions: list of accessions
        - roles: mapping metabolite -> role
        - positive_metabolites: metabolites contributing positively to enzyme
        - negative_metabolites: metabolites contributing negatively to enzyme
        - concat: single string 'Enzyme:Gene:Accession:Role' for excel
        - enzyme_activity: (optional) 'ACTIVE'/'INHIBITED'/'UNCLEAR'
        - node_color: (optional) 'red'/'blue'/'gray'
    """
    enzyme_stats = {}

    # If columns do not exist, return empty
    if 'Enzymes_Accession' not in df.columns and 'Enzyme_Gene_name' not in df.columns and 'Reaction_Role' not in df.columns:
        return {}
    
    # Find the name column (support multiple column name variations)
    name_candidates = ['Name', 'LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite', 'name', 'Metabolite_Name']
    name_col = next((c for c in name_candidates if c in df.columns), 'Name')
    
    for _, row in df.iterrows():
        met = row.get(name_col, row.get('Name', ''))
        log2fc = row.get('log2FC', 0)
        access_col = row.get('Enzymes_Accession', '')
        gene_col = row.get('Enzyme_Gene_name', '')
        role_col = row.get('Reaction_Role', '')

        # Parse possible multiple enzymes separated by commas
        access_items = [s.strip() for s in str(access_col).split(',') if s.strip()]
        gene_items = [s.strip() for s in str(gene_col).replace('|', ',').split(',') if s.strip()]
        role_items = [s.strip() for s in str(role_col).split(',') if s.strip()]

        # If access_items look like 'Name: Accession', normalize
        parsed = []
        if access_items:
            for item in access_items:
                # If item contains ':', assume format 'Enzyme: Accession'
                if ':' in item:
                    parts = [p.strip() for p in item.split(':')]
                    enzyme_name = parts[0]
                    accession = parts[1] if len(parts) > 1 else ''
                else:
                    # try to infer enzyme by index
                    enzyme_name = item
                    accession = ''
                parsed.append((enzyme_name, accession))
        elif gene_items:
            # Fallback: use gene names as enzyme identifiers
            for g in gene_items:
                parsed.append((g, ''))
        else:
            continue

        # Align roles per enzyme if given; otherwise mark unknown
        for i, (enzyme_name, accession) in enumerate(parsed):
            role = role_items[i] if i < len(role_items) else 'unknown'
            
            # Use gene name as key if available to prevent duplicates of same gene
            # Otherwise fall back to enzyme name
            gene_key = gene_items[0] if gene_items else enzyme_name
            key = gene_key
            
            if key not in enzyme_stats:
                enzyme_stats[key] = {
                    'metabolites': [],
                    'genes': [],
                    'accessions': [],
                    'roles': {},
                    'positive_metabolites': [],
                    'negative_metabolites': [],
                    'enzyme_names': [],  # Track all enzyme names for this gene
                    'concat': '',
                }

            enzyme_stats[key]['metabolites'].append(met)
            enzyme_stats[key]['accessions'].append(accession)
            # Track enzyme names for this gene
            if enzyme_name not in enzyme_stats[key]['enzyme_names']:
                enzyme_stats[key]['enzyme_names'].append(enzyme_name)
            # Try to attach a gene if available (take first gene)
            if gene_items:
                enzyme_stats[key]['genes'].extend([g for g in gene_items if g not in enzyme_stats[key]['genes']])
            enzyme_stats[key]['roles'][met] = role

            # Determine contribution sign based on role and metabolite direction
            if role.lower() == 'product':
                # product: upregulated -> positive, downregulated -> negative
                if log2fc > 0:
                    enzyme_stats[key]['positive_metabolites'].append(met)
                elif log2fc < 0:
                    enzyme_stats[key]['negative_metabolites'].append(met)
            elif role.lower() == 'substrate':
                # substrate: reversed
                if log2fc > 0:
                    enzyme_stats[key]['negative_metabolites'].append(met)
                elif log2fc < 0:
                    enzyme_stats[key]['positive_metabolites'].append(met)
            else:
                # unknown: don't attribute direction, but include in metabolites
                pass

    # Post-process and filter by min_metabolites
    filtered = {}
    for enzyme, data in enzyme_stats.items():
        unique_mets = list(dict.fromkeys(data['metabolites']))
        n = len(unique_mets)
        if n >= min_metabolites:
            concat_items = []
            # create a single concatenated string using first gene/accession/role as exemplar
            genes = ','.join(data['genes']) if data['genes'] else ''
            accessions = ','.join([a for a in data['accessions'] if a])
            # For role summary, list unique roles
            roles = ','.join(sorted(set(data['roles'].values())))
            # Include all enzyme names for this gene
            enzyme_names = ','.join(data['enzyme_names']) if data['enzyme_names'] else enzyme
            concat = f"{enzyme_names}: {genes}: {accessions}: {roles}"

            filtered[enzyme] = {
                'n_metabolites': n,
                'metabolites': unique_mets,
                'genes': data['genes'],
                'accessions': data['accessions'],
                'roles': data['roles'],
                'positive_metabolites': data['positive_metabolites'],
                'negative_metabolites': data['negative_metabolites'],
                'concat': concat
            }


    return filtered

def calculate_transporters_statistics(df, min_metabolites=2):
    """Collect transporters and their associated metabolites.

    Returns dict keyed by transporter name with genes, accessions and metabolites.
    """
    transporter_stats = {}
    if 'Transporters' not in df.columns and 'Transporter_Gene_name' not in df.columns and 'Transporter_Uniprot_ID' not in df.columns:
        return {}
    
    # Find the name column (support multiple column name variations)
    name_candidates = ['Name', 'LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite', 'name', 'Metabolite_Name']
    name_col = next((c for c in name_candidates if c in df.columns), 'Name')
    
    for _, row in df.iterrows():
        met = row.get(name_col, row.get('Name', ''))
        log2fc = row.get('log2FC', 0)
        trans_col = row.get('Transporters', '')
        gene_col = row.get('Transporter_Gene_name', '')
        acc_col = row.get('Transporter_Uniprot_ID', '')

        trans_items = [s.strip() for s in str(trans_col).split('|') if s.strip()]
        gene_items = [s.strip() for s in str(gene_col).replace('|', ',').split(',') if s.strip()]
        acc_items = [s.strip() for s in str(acc_col).replace('|', ',').split(',') if s.strip()]

        # Try to align by index; otherwise map all metabolites to all transporters
        if trans_items:
            for i, t in enumerate(trans_items):
                # Use gene name as key if available to prevent duplicates of same gene
                # Otherwise fall back to transporter name
                gene_key = gene_items[0] if gene_items else t
                key = gene_key
                
                if key not in transporter_stats:
                    transporter_stats[key] = {
                        'metabolites': [],
                        'genes': [],
                        'accessions': [],
                        'transporter_names': [],  # Track all transporter names for this gene
                    }
                transporter_stats[key]['metabolites'].append(met)
                # Track transporter names for this gene
                if t not in transporter_stats[key]['transporter_names']:
                    transporter_stats[key]['transporter_names'].append(t)
                if gene_items:
                    transporter_stats[key]['genes'].extend([g for g in gene_items if g not in transporter_stats[key]['genes']])
                if acc_items:
                    transporter_stats[key]['accessions'].extend([a for a in acc_items if a not in transporter_stats[key]['accessions']])

    # Filter
    filtered = {}
    for t, data in transporter_stats.items():
        unique_mets = list(dict.fromkeys(data['metabolites']))
        if len(unique_mets) >= min_metabolites:
            # Include all transporter names for this gene
            transporter_names = ','.join(data['transporter_names']) if data['transporter_names'] else t
            concat = f"{transporter_names}: {','.join(data['genes']) if data['genes'] else ''}: {','.join(data['accessions']) if data['accessions'] else ''}"
            filtered[t] = {
                'n_metabolites': len(unique_mets),
                'metabolites': unique_mets,
                'genes': data['genes'],
                'accessions': data['accessions'],
                'concat': concat,
                'transporter_names': data['transporter_names']
            }

    return filtered

def calculate_associated_diseases_statistics(df, min_metabolites=2):
    """Collect associated diseases mapping to metabolites.

    Returns dict keyed by disease name with metabolites list and count.
    """
    disease_stats = {}
    if 'Associated_Diseases' not in df.columns:
        return {}

    # Find the name column (support multiple column name variations)
    name_candidates = ['Name', 'LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite', 'name', 'Metabolite_Name']
    name_col = next((c for c in name_candidates if c in df.columns), 'Name')

    for _, row in df.iterrows():
        met = row.get(name_col, row.get('Name', ''))
        disease_col = row.get('Associated_Diseases', '')
        if pd.isna(disease_col) or str(disease_col).strip() == '':
            continue
        items = [s.strip() for s in str(disease_col).split('|') if s.strip()]
        for d in items:
            if d not in disease_stats:
                disease_stats[d] = {'metabolites': []}
            disease_stats[d]['metabolites'].append(met)

    # Filter by min_metabolites
    filtered = {}
    for d, data in disease_stats.items():
        unique_mets = list(dict.fromkeys(data['metabolites']))
        if len(unique_mets) >= min_metabolites:
            filtered[d] = {
                'n_metabolites': len(unique_mets),
                'metabolites': unique_mets,
                'concat': d
            }

    return filtered

def create_metabolite_pathway_network(df, pos_set_zscore=1, neg_set_zscore=-1, 
                                       max_pvalue=1.0, min_metabolites=2,
                                       max_total_pathways=None, max_pos=None, max_neg=None,
                                       include_enzymes=True, include_transporters=True, include_diseases=True):
    """Create a network connecting metabolites to pathways based on fold changes
    
    Parameters:
    -----------
    include_enzymes : bool, default=True
        Whether to include enzyme nodes and their connections in the network
    include_transporters : bool, default=True
        Whether to include transporter nodes and their connections in the network
    include_diseases : bool, default=True
        Whether to include disease nodes and their connections in the network
    """
    
    # Calculate pathway statistics first (respect the limits passed in)
    pathway_stats = calculate_pathway_statistics(
        df,
        min_metabolites=min_metabolites,
        max_pvalue=max_pvalue,
        pos_set_zscore=pos_set_zscore,
        neg_set_zscore=neg_set_zscore,
        max_total_pathways=max_total_pathways,
        max_pos=max_pos,
        max_neg=max_neg
    )

    # Calculate enzymes, transporters and diseases using helper functions (conditionally)
    enzyme_stats = calculate_enzymes_statistics(df, min_metabolites=min_metabolites) if include_enzymes else {}
    transporter_stats = calculate_transporters_statistics(df, min_metabolites=min_metabolites) if include_transporters else {}
    disease_stats = calculate_associated_diseases_statistics(df, min_metabolites=min_metabolites) if include_diseases else {}

    # Determine which column contains pathway information (support multiple possible column names)
    pathway_col_candidates = ['All_Pathways', 'Pathways', 'HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways', 'Metabolika_Pathways']
    pathway_col = None
    for col in pathway_col_candidates:
        if col in df.columns:
            pathway_col = col
            break

    # Find the name column (support multiple column name variations)
    name_candidates = ['Name', 'LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite', 'name', 'Metabolite_Name']
    name_col = next((c for c in name_candidates if c in df.columns), 'Name')

    # Filter metabolites to only include those connected to filtered pathways
    connected_metabolites = set()
    for _, row in df.iterrows():
        metabolite = row.get(name_col, row.get('Name', None))
        if metabolite is None:
            continue

        # Safely retrieve pathways list/string
        pathways_raw = row.get(pathway_col, None) if pathway_col is not None else None
        if pathways_raw is None or (isinstance(pathways_raw, float) and np.isnan(pathways_raw)):
            continue

        # Normalize to list
        if isinstance(pathways_raw, str):
            pathways_list = [p.strip() for p in pathways_raw.split('|') if p.strip()]
        elif isinstance(pathways_raw, list):
            pathways_list = pathways_raw
        else:
            # Unexpected type - skip
            continue

        for pathway in pathways_list:
            # Normalize pathway name for matching
            norm_path = normalize_pathway_name_global(pathway)
            if norm_path in pathway_stats:
                connected_metabolites.add(metabolite)
                break
    
    # Create metabolites dict only for connected metabolites
    metabolites = {}
    for metabolite in connected_metabolites:
        metabolite_data = df[df[name_col] == metabolite].iloc[0]
        metabolites[metabolite] = {
            'log2FC': metabolite_data['log2FC'],
            'pvalue': metabolite_data['pvalue']
        }
    
    print(f"DEBUG: Connected metabolites: {len(metabolites)} out of {len(df)} total")
    print(f"DEBUG: Filtered pathways: {len(pathway_stats)}")

    # Create graph
    G = nx.Graph()

    # Helper functions for shading and alpha
    def _alpha_from_mag(m, max_mag=2.0, min_a=0.35, max_a=1.0):
        try:
            m = abs(float(m))
        except Exception:
            m = 0.0
        m = max(0.0, min(max_mag, m))
        a = min_a + (max_a - min_a) * (m / max_mag)
        return round(a, 3)

    def _get_metabolite_color_solid(fc):
        """Return solid color for metabolites based on regulation"""
        try:
            v = float(fc)
        except Exception:
            v = 0.0
        if v > 0:
            return '#FF8C00'  # Solid orange for upregulated
        elif v < 0:
            return '#228B22'  # Solid forest green for downregulated
        else:
            return '#CCCCCC'  # Gray for no change

    def _get_pathway_color_solid(z):
        """Return solid color for pathways based on z-score"""
        try:
            zf = float(z)
        except Exception:
            zf = 0.0
        if zf > 0:
            return '#FF0000'  # Solid red for activated
        elif zf < 0:
            return '#3498DB'  # Solid blue for inhibited
        return '#BDC3C7'  # Gray for no change

    # Add metabolites as nodes (circles)
    for metabolite, data in metabolites.items():
        log2FC = data['log2FC']
        pvalue = data['pvalue']

        if log2FC > 0:  # Upregulated
            regulation = 'upregulated'
        elif log2FC < 0:
            regulation = 'downregulated'
        else:
            regulation = 'no_change'

        color = _get_metabolite_color_solid(log2FC)
        alpha = _alpha_from_mag(log2FC)

        G.add_node(metabolite,
                   node_type='metabolite',
                   color=color,
                   log2_fold_change=log2FC,
                   regulation=regulation,
                   shape='circle',
                   alpha=alpha)

    # Add pathways as nodes (diamonds) - only if they still have at least one connected metabolite
    for pathway, data in pathway_stats.items():
        # Ensure this pathway has at least one metabolite present in the filtered set
        connected_mets_for_path = [m for m in data.get('metabolites', []) if m in metabolites]
        if len(connected_mets_for_path) == 0:
            continue
        # Use solid color based on z-score direction
        color = _get_pathway_color_solid(data['z_score'])
        alpha = _alpha_from_mag(data['z_score'])
        avg_log2_fold_change = data['mean_log2fc']
        
        G.add_node(pathway, 
                  node_type='pathway', 
                  color=color, 
                  status=data['status'],
                  shape='diamond',
                  z_score=data['z_score'],
                  combined_pvalue=data['combined_pvalue'],
                  alpha=alpha,
                  n_metabolites=data['n_metabolites'],
                  avg_log2_fold_change=avg_log2_fold_change)

    # Add enzymes as trapeziums (grey) - only if enabled
    if include_enzymes:
        for enzyme, data in enzyme_stats.items():
            # Limit to SIGNIFICANT metabolites that are actually present in the filtered metabolite set
            connected_mets = [m for m in data.get('metabolites', []) 
                              if m in metabolites and metabolites[m].get('regulation', 'no_change') != 'no_change']
            if len(connected_mets) == 0:
                # Skip orphan enzyme nodes
                continue
            color = '#888888'
            gene_label = data['genes'][0] if data['genes'] else enzyme
            G.add_node(enzyme,
                       node_type='enzyme',
                       color=color,
                       shape='trapezium',
                       genes=data['genes'],
                       accessions=data['accessions'],
                       roles=data['roles'],
                       concat=data['concat'],
                       textlabel=gene_label,
                       n_metabolites=len(connected_mets))

    # Add transporters as rectangles (ash) - only if enabled
    if include_transporters:
        for transporter, data in transporter_stats.items():
            # Limit to SIGNIFICANT metabolites only
            connected_mets = [m for m in data.get('metabolites', []) 
                              if m in metabolites and metabolites[m].get('regulation', 'no_change') != 'no_change']
            if len(connected_mets) == 0:
                continue
            color = '#B2B2B2'  # ash
            gene_label = data['genes'][0] if data['genes'] else transporter
            G.add_node(transporter,
                       node_type='transporter',
                       color=color,
                       shape='rectangle',
                       genes=data['genes'],
                       accessions=data['accessions'],
                       concat=data['concat'],
                       textlabel=gene_label,
                       n_metabolites=len(connected_mets))

    # Add diseases as purple hexagons - only if enabled
    if include_diseases:
        for disease, data in disease_stats.items():
            # Limit to SIGNIFICANT metabolites only
            connected_mets = [m for m in data.get('metabolites', []) 
                              if m in metabolites and metabolites[m].get('regulation', 'no_change') != 'no_change']
            if len(connected_mets) == 0:
                continue
            # Standardize disease color to purple
            color = '#8e44ad'
            G.add_node(disease,
                       node_type='disease',
                       color=color,
                       shape='hexagon',
                       concat=data['concat'],
                       n_metabolites=len(connected_mets),
                       metabolites=connected_mets)

    # Add edges (connect metabolites to pathways)
    for _, row in df.iterrows():
        metabolite = row.get('Name', None)
        if metabolite is None:
            continue

        # Safe log2FC and pvalue
        log2FC = row.get('log2FC', None)
        if log2FC is None or (isinstance(log2FC, float) and np.isnan(log2FC)):
            log2FC = 1.0

        # Determine pathways column (reuse earlier detection if available)
        # Try common field names
        pathways_raw = None
        for col in ['All_Pathways', 'HMDB_Pathways', 'Pathways', 'HMDB_Pathway']:
            if col in row.index:
                pathways_raw = row.get(col)
                break

        if pathways_raw is None or (isinstance(pathways_raw, float) and np.isnan(pathways_raw)):
            continue

        if isinstance(pathways_raw, str):
            pathways_list = [p.strip() for p in pathways_raw.split('|') if p.strip()]
        elif isinstance(pathways_raw, list):
            pathways_list = pathways_raw
        else:
            continue

        for pathway in pathways_list:
            norm_path = normalize_pathway_name_global(pathway)
            if norm_path in pathway_stats and metabolite in metabolites:
                if log2FC > 0:  # Positive fold change - Activation (red)
                    edge_color = 'red'
                    edge_type = 'activating'
                elif log2FC < 0:  # Negative fold change - Inhibition (blue)
                    edge_color = 'blue'
                    edge_type = 'inhibiting'
                else:  # No change (log2FC exactly 0) - Neutral (light gray)
                    edge_color = '#cccccc'
                    edge_type = 'neutral'

                G.add_edge(metabolite, norm_path, 
                          color=edge_color, 
                          edge_type=edge_type,
                          weight=abs(log2FC),
                          log2_fold_change=log2FC)

    # Add edges for enzymes: metabolites -> enzyme with sign depending on role
    # Only connect SIGNIFICANT metabolites (not gray/no_change)
    if include_enzymes:
        for enzyme, data in enzyme_stats.items():
            for met in [m for m in data.get('metabolites', []) 
                        if m in metabolites and metabolites[m].get('regulation', 'no_change') != 'no_change']:
                # determine sign
                role = data['roles'].get(met, 'unknown')
                # find log2fc for metabolite
                log2fc = metabolites[met]['log2FC']

                if role.lower() == 'product':
                    # Map to plotting edge groups so they render
                    if log2fc > 0:
                        edge_color = 'red'
                        edge_type = 'activating'
                    elif log2fc < 0:
                        edge_color = 'blue'
                        edge_type = 'inhibiting'
                    else:
                        edge_color = '#cccccc'
                        edge_type = 'neutral'
                elif role.lower() == 'substrate':
                    # Reverse interpretation for substrates
                    if log2fc > 0:
                        edge_color = 'blue'
                        edge_type = 'inhibiting'
                    elif log2fc < 0:
                        edge_color = 'red'
                        edge_type = 'activating'
                    else:
                        edge_color = '#cccccc'
                        edge_type = 'neutral'
                else:
                    edge_color = '#cccccc'
                    edge_type = 'neutral'

                G.add_edge(met, enzyme,
                           color=edge_color,
                           edge_type=edge_type,
                           weight=abs(log2fc),
                           role=role,
                           log2_fold_change=log2fc)

    # Add edges for transporters: link similarly to pathways (up -> red, down -> blue)
    # Only connect SIGNIFICANT metabolites (not gray/no_change)
    if include_transporters:
        for transporter, data in transporter_stats.items():
            for met in [m for m in data.get('metabolites', []) 
                        if m in metabolites and metabolites[m].get('regulation', 'no_change') != 'no_change']:
                log2fc = metabolites[met]['log2FC']
                if log2fc > 0:
                    edge_color = 'red'
                    edge_type = 'activating'
                elif log2fc < 0:
                    edge_color = 'blue'
                    edge_type = 'inhibiting'
                else:
                    edge_color = '#cccccc'
                    edge_type = 'neutral'
                G.add_edge(met, transporter,
                           color=edge_color,
                           edge_type=edge_type,
                           weight=abs(log2fc),
                           log2_fold_change=log2fc)

    # Add edges for diseases: link like transporters (no gene labels)
    # Only connect SIGNIFICANT metabolites (not gray/no_change)
    if include_diseases:
        for disease, data in disease_stats.items():
            for met in [m for m in data.get('metabolites', []) 
                        if m in metabolites and metabolites[m].get('regulation', 'no_change') != 'no_change']:
                log2fc = metabolites[met]['log2FC']
                if log2fc > 0:
                    edge_color = 'red'
                    edge_type = 'activating'
                elif log2fc < 0:
                    edge_color = 'blue'
                    edge_type = 'inhibiting'
                else:
                    edge_color = '#cccccc'
                    edge_type = 'neutral'
                G.add_edge(met, disease,
                           color=edge_color,
                           edge_type=edge_type,
                           weight=abs(log2fc),
                           log2_fold_change=log2fc)
    
    # Remove orphan gray metabolite nodes (no_change metabolites with no connections)
    orphan_gray_nodes = []
    for node, data in G.nodes(data=True):
        if data.get('node_type') == 'metabolite' and data.get('regulation') == 'no_change':
            if G.degree(node) == 0:  # type: ignore
                orphan_gray_nodes.append(node)
    
    if orphan_gray_nodes:
        print(f"🧹 Removing {len(orphan_gray_nodes)} orphan gray (non-significant) metabolite nodes")
        G.remove_nodes_from(orphan_gray_nodes)
        # Also remove from metabolites dict
        for node in orphan_gray_nodes:
            if node in metabolites:
                del metabolites[node]
    
    # Debug information about connections
    print(f"🔗 Network Connection Summary:")
    metabolite_nodes = [n for n, d in G.nodes(data=True) if d['node_type'] == 'metabolite']
    pathway_nodes = [n for n, d in G.nodes(data=True) if d['node_type'] == 'pathway']
    enzyme_nodes = [n for n, d in G.nodes(data=True) if d['node_type'] == 'enzyme']
    transporter_nodes = [n for n, d in G.nodes(data=True) if d['node_type'] == 'transporter']
    disease_nodes = [n for n, d in G.nodes(data=True) if d['node_type'] == 'disease']
    
    print(f"   Metabolites: {len(metabolite_nodes)}")
    print(f"   Pathways: {len(pathway_nodes)}")
    print(f"   Enzymes: {len(enzyme_nodes)}")
    print(f"   Transporters: {len(transporter_nodes)}")
    print(f"   Diseases: {len(disease_nodes)}")
    print(f"   Total edges: {len(G.edges())}")
    
    # Count edge types
    edge_types = {}
    for u, v, d in G.edges(data=True):
        u_type = G.nodes[u]['node_type']
        v_type = G.nodes[v]['node_type']
        connection_type = f"{u_type}-{v_type}"
        edge_types[connection_type] = edge_types.get(connection_type, 0) + 1
    
    print(f"🔗 Edge types:")
    for edge_type, count in edge_types.items():
        print(f"   {edge_type}: {count}")
    
    return G, metabolites, pathway_stats

def create_metabolite_pathway_network_from_stats(df, pathway_stats, 
                                                  pos_set_fc=0, neg_set_fc=0,
                                                  include_upstream=True, include_diseases=True,
                                                  upstream_data=None, disease_data=None):
    """Create a network from pre-calculated pathway statistics
    
    This function creates the network directly from already filtered pathway statistics
    to ensure the visualization shows exactly the pathways that met the user's criteria.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        The original metabolomics data
    pathway_stats : dict
        Already filtered pathway statistics from calculate_pathway_statistics
    pos_set_fc : float, default=0
        Positive fold change threshold for determining upregulation
    neg_set_fc : float, default=0
        Negative fold change threshold for determining downregulation
    include_upstream : bool, default=True
        Whether to include upstream regulator nodes (enzymes/transporters) and their connections
    include_diseases : bool, default=True
        Whether to include disease nodes and their connections in the network
    upstream_data : pandas.DataFrame, optional
        DataFrame containing upstream regulator data with columns: Name, Gene, Type, Associated_Metabolites, # Metabolites
    disease_data : pandas.DataFrame, optional
        DataFrame containing disease data with columns: Disease, Associated_Metabolites, # Metabolites
    """
    
    G = nx.Graph()
    
    # Helper to safely extract metabolites list from a stats dict (supports alternate keys and string lists)
    import re as _re
    def _extract_stats_metabolites(stats_obj):
        if not isinstance(stats_obj, dict):
            return []
        for key in ['metabolites', 'Metabolites', 'metabolites_list', 'Metabolites_List', 'Metabolites_in_Pathway', 'Metabolite_Names']:
            if key in stats_obj and stats_obj[key] is not None:
                val = stats_obj[key]
                if isinstance(val, (list, tuple, set)):
                    return [str(x).strip() for x in val if str(x).strip()]
                if isinstance(val, str):
                    parts = _re.split(r"[|,;]", val)
                    return [p.strip() for p in parts if p and p.strip()]
        return []
    
    # Get all metabolites that appear in the filtered pathways
    all_metabolites_in_pathways = set()
    for _, stats in pathway_stats.items():
        all_metabolites_in_pathways.update(_extract_stats_metabolites(stats))
    
    # If DataFrame already has canonical columns (provided by GUI verification),
    # prefer those and avoid noisy debug prints.
    # Try multiple name variations for flexibility
    name_candidates = ['Name', 'LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite', 'name', 'Metabolite_Name']
    log2fc_candidates = ['log2FC', 'log2_FC', 'Log2FC', 'log2FoldChange', 'FC']
    pvalue_candidates = ['pvalue', 'P-Value', 'p-value', 'p_value', 'P_value', 'pValue']
    
    name_col = next((c for c in name_candidates if c in df.columns), None)
    log2fc_col = next((c for c in log2fc_candidates if c in df.columns), None)
    pvalue_col = next((c for c in pvalue_candidates if c in df.columns), None)
    
    if name_col and log2fc_col and pvalue_col:
        # Minimal info-level logging instead of verbose prints
        try:
            import logging as _logging
            _logging.getLogger(__name__).info(f"Using columns: Name='{name_col}', log2FC='{log2fc_col}', pvalue='{pvalue_col}'")
        except Exception:
            pass
    else:
        print(f"📊 DEBUG: Pathway stats contain {len(all_metabolites_in_pathways)} unique metabolite names")
        print(f"   - Metabolites in pathways: {sorted(list(all_metabolites_in_pathways))[:10]}...")  # Show first 10

        # Detect column names (handle variations)
        print(f"\n📋 DEBUG: DataFrame columns: {df.columns.tolist()}")
        name_col = None
        log2fc_col = None
        pvalue_col = None

        # Try to find Name column (comprehensive list including lipid columns)
        for col in ['Name', 'LipidID', 'Lipid_ID', 'name', 'Metabolite', 'metabolite', 'Metabolite_Name', 'Compound_Name', 'Feature_ID']:
            if col in df.columns:
                name_col = col
                break

        # Try to find log2FC column  
        for col in ['log2FC', 'log2_FC', 'Log2FC', 'Log2_FC', 'log2FoldChange', 'FC']:
            if col in df.columns:
                log2fc_col = col
                break

        # Try to find pvalue column
        for col in ['pvalue', 'p-value', 'p_value', 'P-value', 'P_value', 'pValue']:
            if col in df.columns:
                pvalue_col = col
                break

        if not name_col or not log2fc_col or not pvalue_col:
            print(f"⚠️ WARNING: Could not find required columns!")
            print(f"   - Name column: {name_col}")
            print(f"   - log2FC column: {log2fc_col}")
            print(f"   - pvalue column: {pvalue_col}")
            # Use defaults as fallback
            name_col = name_col or 'Name'
            log2fc_col = log2fc_col or 'log2FC'
            pvalue_col = pvalue_col or 'pvalue'

        print(f"✅ Using columns: Name='{name_col}', log2FC='{log2fc_col}', pvalue='{pvalue_col}'")
    
    # Create metabolites dictionary - include ALL metabolites in df
    # This ensures disease/upstream metabolites appear even without pathway connections
    metabolites = {}
    skipped_count = 0
    metabolites_with_pathways = 0
    metabolites_without_pathways = 0
    
    for index, row in df.iterrows():
        try:
            metabolite_name = row[name_col]
        except KeyError:
            print(f"⚠️ Row {index}: Cannot find name column '{name_col}'")
            skipped_count += 1
            continue
        
        # Include ALL metabolites in df (pathways, diseases, upstream)
        try:
            log2fc = float(row[log2fc_col])
            pvalue = float(row[pvalue_col])
            
            if log2fc >= pos_set_fc:  # Use pos_set_fc threshold for upregulation
                regulation = 'upregulated'
                color = '#FF8C00'  # orange for upregulated metabolite nodes
            elif log2fc <= neg_set_fc:  # Use neg_set_fc threshold for downregulation
                regulation = 'downregulated'
                color = 'green'
            else:
                regulation = 'no_change'
                color = 'lightgray'
            
            alpha = max(0.3, min(1.0, -np.log10(pvalue) / 10))
            metabolites[metabolite_name] = {
                'log2FC': log2fc,
                'p_value': pvalue,
                'regulation': regulation,
                'color': color,
                'alpha': alpha
            }
            
            # Track if metabolite has pathway connections
            if metabolite_name in all_metabolites_in_pathways:
                metabolites_with_pathways += 1
            else:
                metabolites_without_pathways += 1
                
        except (ValueError, TypeError, KeyError) as e:
            # Skip metabolites with invalid data
            print(f"⚠️ Skipping {metabolite_name}: {e}")
            skipped_count += 1
            continue
    
    if skipped_count > 0:
        print(f"⚠️ Skipped {skipped_count} metabolites due to missing/invalid data")
    
    print(f"📊 Creating network from filtered stats:")
    print(f"   - {len(metabolites)} total metabolites included in network")
    print(f"   - {metabolites_with_pathways} metabolites with pathway connections")
    print(f"   - {metabolites_without_pathways} metabolites without pathways (disease/upstream only)")
    print(f"   - {len(pathway_stats)} pathways")
    print(f"   - {len(all_metabolites_in_pathways)} metabolites expected from pathway stats")
    print(f"   - {len(metabolites)} total metabolites included in network")
    print(f"   - {metabolites_with_pathways} metabolites with pathway connections")
    print(f"   - {metabolites_without_pathways} metabolites without pathways (disease/upstream only)")
    print(f"   - {len(pathway_stats)} pathways")
    print(f"   - {len(all_metabolites_in_pathways)} metabolites expected from pathway stats")
    
    # DEBUG: Show metabolites without pathway connections
    if metabolites_without_pathways > 0:
        no_pathway_mets = [m for m in metabolites.keys() if m not in all_metabolites_in_pathways]
        print(f"\n🔍 DEBUG: Metabolites WITHOUT pathway connections (should connect to disease/upstream):")
        for met in no_pathway_mets[:10]:  # Show first 10
            print(f"   - {met}: log2FC={metabolites[met]['log2FC']:.3f}, pvalue={metabolites[met]['p_value']:.6f}")
    
    # Show which metabolites from pathways are missing
    missing_metabolites = all_metabolites_in_pathways - set(metabolites.keys())
    if missing_metabolites:
        print(f"⚠️ WARNING: {len(missing_metabolites)} metabolites in pathways but not in network:")
        for met in list(missing_metabolites)[:5]:  # Show first 5
            print(f"   - {met}")
    
    # DEBUG: Show metabolite log2FC values and regulation status
    print(f"\n🔍 DEBUG: Metabolite regulation status:")
    for met_name, met_data in metabolites.items():
        print(f"  📍 {met_name}:")
        print(f"    ↳ log2FC: {met_data['log2FC']:.3f}, |log2FC|: {abs(met_data['log2FC']):.3f}")
        print(f"    ↳ regulation: {met_data['regulation']}, color: {met_data['color']}")
        print(f"    ↳ p-value: {met_data['p_value']:.6f}")
    
    # DEBUG: Show which metabolites are expected to be in each pathway
    print(f"\n🔍 DEBUG: Pathways and their metabolites from pathway_stats:")
    for pathway_name, stats in pathway_stats.items():
        mets_list = _extract_stats_metabolites(stats)
        st = stats.get('status', 'No Change')
        z = stats.get('z_score', 0.0)
        print(f"  📍 {str(pathway_name)[:60]}...")
        print(f"    ↳ Status: {st}, Z-score: {float(z):.3f}")
        print(f"    ↳ Metabolites ({len(mets_list)}): {mets_list}")
        # Check which of these metabolites are actually in our metabolites dict
        in_dict = [m for m in mets_list if m in metabolites]
        not_in_dict = [m for m in mets_list if m not in metabolites]
        if in_dict:
            print(f"    ✅ In metabolites dict: {in_dict}")
        if not_in_dict:
            print(f"    ❌ NOT in metabolites dict: {not_in_dict}")
        print()
    
    # Add metabolite nodes
    for metabolite, data in metabolites.items():
        G.add_node(metabolite, 
                  node_type='metabolite',
                  log2_fold_change=data['log2FC'],
                  regulation=data['regulation'],
                  color=data['color'],
                  alpha=data['alpha'])
    
    # Add pathway nodes and connections
    for pathway_name, stats in pathway_stats.items():
        pathway_metabolites = _extract_stats_metabolites(stats)
        
        # Find which metabolites/classes in our network connect to this pathway
        # SMART LOGIC: Class nodes connect if ≥50% of members are in this pathway
        valid_mets = []
        for met_name in metabolites.keys():
            met_data = metabolites[met_name]
            is_connected = False
            
            if met_data.get('is_class', False):
                # CLASS NODE: Connect if this pathway is in pathways_meeting_threshold (≥50% of class members in pathway)
                pathways_meeting_threshold = met_data.get('pathways_meeting_threshold', {})
                if pathway_name in pathways_meeting_threshold:
                    is_connected = True
            else:
                # INDIVIDUAL NODE: Connect if this individual is in the pathway metabolites
                # (Individual nodes are only present if they're NOT part of any class)
                if met_name in pathway_metabolites:
                    is_connected = True
            
            if is_connected:
                valid_mets.append(met_name)
        
        if not valid_mets:
            continue
            
        # Determine pathway status and color
        if stats.get('status') == 'Activated':
            color = 'red'   # pathways: red for activated
        elif stats.get('status') == 'Inhibited':
            color = 'blue'
        else:
            color = '#cccccc'
        
        # Add pathway node
        G.add_node(pathway_name,
                  node_type='pathway',
                  status=stats['status'],
                  combined_pvalue=stats['combined_pvalue'],
                  z_score=stats['z_score'],
                  color=color,
                  alpha=0.8)
        
        # Connect pathway to its metabolites/classes based ONLY on metabolite regulation
        # Upregulated => activating (orange); Downregulated => inhibiting (blue)
        for metabolite in valid_mets:
            if metabolite in metabolites:
                met_regulation = metabolites[metabolite]['regulation']
                if met_regulation == 'upregulated':
                    edge_type = 'activating'
                    edge_color = '#FF8C00'  # Orange for upregulated (dotted)
                    G.add_edge(metabolite, pathway_name, edge_type=edge_type, color=edge_color, weight=1)
                elif met_regulation == 'downregulated':
                    edge_type = 'inhibiting'
                    edge_color = '#3498DB'  # Blue for downregulated
                    G.add_edge(metabolite, pathway_name, edge_type=edge_type, color=edge_color, weight=1)
                # Note: neutral metabolites not connected (by design)
    
    # DEBUG: Print detailed connection information for each metabolite
    print(f"\n🔍 DEBUG: Detailed metabolite-pathway connections:")
    for metabolite in metabolites.keys():
        connected_pathways = []
        for neighbor in G.neighbors(metabolite):
            if G.nodes[neighbor]['node_type'] == 'pathway':
                edge_data = G.get_edge_data(metabolite, neighbor)
                # Show full pathway name instead of truncating
                connected_pathways.append(f"{neighbor} (edge_type: {edge_data.get('edge_type', 'unknown')})")
        
        if connected_pathways:
            print(f"  📍 {metabolite}: Connected to {len(connected_pathways)} pathways")
            for pathway in connected_pathways:
                print(f"    ↳ {pathway}")
        else:
            print(f"  ⚠️  {metabolite}: NO PATHWAY CONNECTIONS!")
    
    # Add upstream regulator and disease nodes if requested
    if include_upstream or include_diseases:
        print(f"🔍 DEBUG: Network inclusion statistics:")
        print(f"   - Include upstream regulators: {include_upstream}")
        print(f"   - Include diseases: {include_diseases}")
        
        # Add upstream regulator nodes (purple) - only if enabled and data provided
        if include_upstream and upstream_data is not None and len(upstream_data) > 0:
            print(f"   - Processing {len(upstream_data)} upstream regulators")
            print(f"   - Upstream data columns: {upstream_data.columns.tolist()}")
            
            # Dictionary to consolidate genes: {gene_name: {metabolites: set(), full_names: set(), type: str}}
            consolidated_upstream = {}
            
            for idx, row in upstream_data.iterrows():
                upstream_name = row.get('Name', '')
                # Try to get Gene column, fallback to Name if not present
                gene_name = row.get('Gene', upstream_name) 
                upstream_type = row.get('Type', 'Unknown')
                associated_metabolites_str = row.get('Associated_Metabolites', '')
                
                # Skip if name is NaN or empty
                if pd.isna(upstream_name) or str(upstream_name).strip() == '' or str(upstream_name).lower() == 'nan':
                    print(f"   ⚠️ Skipping empty/NaN upstream name at row {idx}")
                    continue
                
                # Parse associated metabolites (support '|', ',', ';' delimiters)
                if pd.isna(associated_metabolites_str) or str(associated_metabolites_str).strip() == '':
                    associated_mets = []
                else:
                    parts = re.split(r"[|,;]", str(associated_metabolites_str))
                    associated_mets = [m.strip() for m in parts if m.strip()]
                
                # Only include upstream regulators with at least one SIGNIFICANT metabolite that is ALSO connected to pathways
                # Filter out 'no_change' (gray) metabolites AND metabolites not connected to pathways
                # This ensures FOCUSED MODE: upstream only connect to pathway-linked significant metabolites
                connected_mets = [m for m in associated_mets 
                                  if m in metabolites 
                                  and metabolites[m]['regulation'] != 'no_change'
                                  and m in all_metabolites_in_pathways]  # CRITICAL: Only pathway-connected
                if len(connected_mets) == 0:
                    print(f"   ⚠️ Skipping {upstream_name}: no SIGNIFICANT metabolites in filtered network")
                    continue
                
                # Use gene name as primary key; if Gene column exists use it, otherwise use Name
                primary_key = gene_name if (gene_name and not pd.isna(gene_name) and str(gene_name).strip() != '' and str(gene_name).lower() != 'nan') else upstream_name
                
                # Consolidate by gene name
                if primary_key not in consolidated_upstream:
                    consolidated_upstream[primary_key] = {
                        'metabolites': set(),
                        'full_names': set(),
                        'type': upstream_type,
                        'gene': gene_name if (gene_name and not pd.isna(gene_name) and str(gene_name).strip() != '') else upstream_name
                    }
                
                # Add metabolites and full names to this gene's entry
                consolidated_upstream[primary_key]['metabolites'].update(connected_mets)
                consolidated_upstream[primary_key]['full_names'].add(upstream_name)
                # Keep the type (Enzyme or Transporter)
                if upstream_type in ['Enzyme', 'Transporter']:
                    consolidated_upstream[primary_key]['type'] = upstream_type
            
            # Now add consolidated nodes to the graph
            print(f"   - Consolidated to {len(consolidated_upstream)} unique genes/regulators")
            if len(consolidated_upstream) == 0:
                print(f"   ⚠️ WARNING: No upstream regulators added to network!")
            
            for gene_key, data in consolidated_upstream.items():
                # Use purple color for all upstream regulators (#8e44ad)
                color = '#8e44ad'
                
                connected_mets = list(data['metabolites'])
                full_names_list = list(data['full_names'])
                
                # Create hover text with all full names
                full_names_str = '; '.join(full_names_list)
                
                print(f"   ✅ Adding upstream node: {gene_key} ({data['type']}) with {len(connected_mets)} metabolites")
                
                G.add_node(gene_key,
                           node_type='upstream',
                           upstream_type=data['type'],  # 'Enzyme' or 'Transporter'
                           color=color,
                           gene=data['gene'],
                           textlabel=gene_key,  # Display the gene name
                           full_names=full_names_str,  # Store all full names for hover text
                           metabolites=connected_mets,
                           n_metabolites=len(connected_mets))
                
                # Connect upstream regulators to SIGNIFICANT metabolites only
                for met in connected_mets:
                    if met in metabolites and metabolites[met]['regulation'] != 'no_change':
                        edge_color = '#888888'  # Gray edge color for associations
                        edge_type = 'associated'
                        G.add_edge(met, gene_key,
                                   color=edge_color,
                                   edge_type=edge_type,
                                   weight=1)
        
        # Add disease nodes (purple) - only if enabled and data provided
        if include_diseases and disease_data is not None and len(disease_data) > 0:
            print(f"   - Processing {len(disease_data)} diseases")
            print(f"   - Disease data columns: {disease_data.columns.tolist()}")
            
            for idx, row in disease_data.iterrows():
                # Try both 'Disease' and 'Disease Name' columns
                disease_name = row.get('Disease Name', row.get('Disease', ''))
                associated_metabolites_str = row.get('Associated_Metabolites', '')
                
                # Skip if disease name is NaN or empty
                if pd.isna(disease_name) or str(disease_name).strip() == '' or str(disease_name).lower() == 'nan':
                    print(f"   ⚠️ Skipping empty/NaN disease name at row {idx}")
                    continue
                
                # Parse associated metabolites (support '|', ',', ';' delimiters)
                if pd.isna(associated_metabolites_str) or str(associated_metabolites_str).strip() == '':
                    associated_mets = []
                else:
                    parts = re.split(r"[|,;]", str(associated_metabolites_str))
                    associated_mets = [m.strip() for m in parts if m.strip()]
                
                # Only include diseases with at least one SIGNIFICANT metabolite that is ALSO connected to pathways
                # Filter out 'no_change' (gray) metabolites AND metabolites not connected to pathways
                # This ensures FOCUSED MODE: diseases only connect to pathway-linked significant metabolites
                connected_mets = [m for m in associated_mets 
                                  if m in metabolites 
                                  and metabolites[m]['regulation'] != 'no_change'
                                  and m in all_metabolites_in_pathways]  # CRITICAL: Only pathway-connected
                
                # DEBUG: Show which metabolites were found vs not found
                if len(associated_mets) != len(connected_mets):
                    missing_mets = [m for m in associated_mets if m not in metabolites]
                    nonsig_mets = [m for m in associated_mets if m in metabolites and metabolites[m]['regulation'] == 'no_change']
                    print(f"   ⚠️ {disease_name}: {len(associated_mets)} total metabolites, {len(connected_mets)} significant in network")
                    if missing_mets:
                        print(f"      Missing from network: {missing_mets}")
                    if nonsig_mets:
                        print(f"      Non-significant (excluded): {nonsig_mets}")
                
                if len(connected_mets) == 0:
                    print(f"   ⚠️ Skipping {disease_name}: no SIGNIFICANT metabolites in filtered network")
                    continue
                
                # Use purple color for diseases
                color = '#8e44ad'
                
                print(f"   ✅ Adding disease node: {disease_name} with {len(connected_mets)} metabolites")
                
                G.add_node(disease_name,
                           node_type='disease',
                           color=color,
                           metabolites=connected_mets,
                           n_metabolites=len(connected_mets))
                
                # Connect diseases to SIGNIFICANT metabolites only
                for met in connected_mets:
                    if met in metabolites and metabolites[met]['regulation'] != 'no_change':
                        edge_color = '#888888'  # Gray edge color for associations
                        edge_type = 'associated'
                        G.add_edge(met, disease_name,
                                   color=edge_color,
                                   edge_type=edge_type,
                                   weight=1)
    
    # Remove orphan gray metabolite nodes (no_change metabolites with no connections)
    orphan_gray_nodes = []
    for node, data in G.nodes(data=True):
        if data.get('node_type') == 'metabolite' and data.get('regulation') == 'no_change':
            if G.degree(node) == 0:  # type: ignore
                orphan_gray_nodes.append(node)
    
    if orphan_gray_nodes:
        print(f"🧹 Removing {len(orphan_gray_nodes)} orphan gray (non-significant) metabolite nodes")
        G.remove_nodes_from(orphan_gray_nodes)
        # Also remove from metabolites dict
        for node in orphan_gray_nodes:
            if node in metabolites:
                del metabolites[node]
    
    # Print network statistics
    print(f"🔗 Network created with:")
    print(f"   - {len([n for n, d in G.nodes(data=True) if d['node_type'] == 'metabolite'])} metabolites")
    print(f"   - {len([n for n, d in G.nodes(data=True) if d['node_type'] == 'pathway'])} pathways")
    print(f"   - {len([n for n, d in G.nodes(data=True) if d['node_type'] == 'upstream'])} upstream regulators")
    print(f"   - {len([n for n, d in G.nodes(data=True) if d['node_type'] == 'disease'])} diseases")
    print(f"   - {G.number_of_edges()} total edges")
    
    # Count edge types
    edge_types = {}
    for edge in G.edges(data=True):
        edge_type = edge[2].get('edge_type', 'unknown')
        edge_types[edge_type] = edge_types.get(edge_type, 0) + 1
    
    print(f"🔗 Edge types:")
    for edge_type, count in edge_types.items():
        print(f"   {edge_type}: {count}")
    
    return G, metabolites, pathway_stats

def get_layout_positions(G, layout_type='spring', k_spacing=1.5, iterations=100):
    """Get node positions for different layout algorithms
    
    Parameters:
    -----------
    G : networkx.Graph
        The network graph
    layout_type : str
        Type of layout ('spring', 'circular', 'shell', 'kamada_kawai', 'bipartite')
    k_spacing : float, default=1.5
        Optimal distance between nodes (lower = more compact)
    iterations : int, default=100
        Number of iterations for spring layout (higher = better positioning)
    """
    if layout_type == 'spring':
        return nx.spring_layout(G, k=k_spacing, iterations=iterations)
    elif layout_type == 'circular':
        return nx.circular_layout(G, scale=float(k_spacing))
    elif layout_type == 'shell':
        # Group nodes by type for shell layout
        metabolites = [n for n, d in G.nodes(data=True) if d['node_type'] == 'metabolite']
        pathways = [n for n, d in G.nodes(data=True) if d['node_type'] == 'pathway']
        upstream = [n for n, d in G.nodes(data=True) if d['node_type'] == 'upstream']
        diseases = [n for n, d in G.nodes(data=True) if d['node_type'] == 'disease']
        
        # Create shells: inner=metabolites, middle=pathways, outer=upstream+diseases
        shells = []
        if metabolites:
            shells.append(metabolites)
        if pathways:
            shells.append(pathways)
        others = upstream + diseases
        if others:
            shells.append(others)
            
        return nx.shell_layout(G, nlist=shells, scale=float(k_spacing))
    elif layout_type == 'kamada_kawai':
        return nx.kamada_kawai_layout(G, scale=float(k_spacing))
    elif layout_type == 'bipartite':
        # Separate metabolites and non-metabolites
        metabolites = [n for n, d in G.nodes(data=True) if d['node_type'] == 'metabolite']
        pathways = [n for n, d in G.nodes(data=True) if d['node_type'] == 'pathway']
        enzymes = [n for n, d in G.nodes(data=True) if d['node_type'] == 'enzyme']
        transporters = [n for n, d in G.nodes(data=True) if d['node_type'] == 'transporter']
        diseases = [n for n, d in G.nodes(data=True) if d['node_type'] == 'disease']
        
        pos = {}
        # Adjust spacing based on k_spacing parameter
        vertical_spacing = k_spacing * 0.5  # Vertical spacing between nodes
        horizontal_spacing = k_spacing      # Horizontal spacing between columns
        
        # Place metabolites on the left
        for i, met in enumerate(metabolites):
            pos[met] = (0, (i - len(metabolites)/2) * vertical_spacing)
        
        # Place pathways in the middle
        for i, path in enumerate(pathways):
            pos[path] = (horizontal_spacing, (i - len(pathways)/2) * vertical_spacing)
            
        # Place enzymes, transporters, diseases on the right
        others = enzymes + transporters + diseases
        for i, other in enumerate(others):
            pos[other] = (2 * horizontal_spacing, (i - len(others)/2) * vertical_spacing)
            
        return pos
    else:
        return nx.spring_layout(G, k=k_spacing, iterations=iterations)

def adjust_positions_for_text_length(G, pos, metabolite_max_chars=30, pathway_max_chars=30, 
                                    metabolite_max_lines=1, pathway_max_lines=3):
    """Adjust node positions to prevent text cutoff by moving nodes away from boundaries
    
    This function scales and shifts node positions to create a safe zone where text
    won't be cut off, regardless of text length.
    """
    if not pos:
        return pos
    
    # Get all positions
    x_coords = [coord[0] for coord in pos.values()]
    y_coords = [coord[1] for coord in pos.values()]
    
    if not x_coords or not y_coords:
        return pos
    
    # Find current bounds
    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)
    
    # Calculate text extent factors - be more generous with space
    max_text_length = max(metabolite_max_chars, pathway_max_chars)
    max_lines = max(metabolite_max_lines, pathway_max_lines)
    
    # More conservative scaling to provide ample text space
    # Shrink the network more to leave space for text
    base_scale = 0.6  # Start with 60% of original size
    text_penalty = max_text_length / 100.0  # Reduce scale based on text length
    line_penalty = max_lines / 10.0         # Reduce scale based on number of lines
    
    # Final scale - never go below 30% but be more aggressive about shrinking
    final_scale = max(0.3, base_scale - text_penalty - line_penalty)
    
    print(f"📐 Position adjustment: scale={final_scale:.2f}, max_text={max_text_length}, max_lines={max_lines}")
    
    # Calculate ranges
    x_range = x_max - x_min if x_max != x_min else 1
    y_range = y_max - y_min if y_max != y_min else 1
    
    # Calculate centers
    x_center = (x_min + x_max) / 2
    y_center = (y_min + y_max) / 2
    
    # Create adjusted positions with more spacing
    adjusted_pos = {}
    import random
    for node, (x, y) in pos.items():
        # Scale relative to center, then shift to ensure no boundary touching
        scaled_x = x_center + (x - x_center) * final_scale
        scaled_y = y_center + (y - y_center) * final_scale

        # Add some randomness to prevent perfect overlap
        random.seed(hash(node) % 1000)  # Consistent randomness per node
        jitter_x = random.uniform(-0.02, 0.02)
        jitter_y = random.uniform(-0.02, 0.02)

        adjusted_pos[node] = (scaled_x + jitter_x, scaled_y + jitter_y)

    print(f"🔧 Position adjustment for text safety:")
    print(f"   - Max text length: {max_text_length} chars, Max lines: {max_lines}")
    print(f"   - Scale factor applied: {final_scale:.2f}")
    print(f"   - Original bounds: x=[{x_min:.2f}, {x_max:.2f}], y=[{y_min:.2f}, {y_max:.2f}]")

    return adjusted_pos

def create_diamond_trace(x, y, size: float = 10.0, color: str = 'blue', alpha: float = 0.5, text: str = '', node_id: str = ''):
    """Create a diamond shape using scatter plot"""
    # Diamond coordinates (relative to center)
    diamond_x = np.array([0, size/20, 0, -size/20, 0]) + x
    diamond_y = np.array([size/20, 0, -size/20, 0, size/20]) + y
    
    # Convert color to rgba with alpha (accept hex)
    if isinstance(color, str) and color.startswith('#') and len(color) == 7:
        try:
            r = int(color[1:3], 16)
            g = int(color[3:5], 16)
            b = int(color[5:7], 16)
            rgba_color = f'rgba({r}, {g}, {b}, {alpha})'
        except Exception:
            rgba_color = f'rgba(128, 128, 128, {alpha})'
    elif color == 'red':
        rgba_color = f'rgba(255, 0, 0, {alpha})'
    elif color == 'blue':
        rgba_color = f'rgba(0, 0, 255, {alpha})'
    elif color == 'purple':
        rgba_color = f'rgba(128, 0, 128, {alpha})'
    elif color == 'green':
        rgba_color = f'rgba(0, 128, 0, {alpha})'
    else:
        rgba_color = f'rgba(128, 128, 128, {alpha})'
    
    return go.Scatter(
        x=diamond_x, y=diamond_y,
        fill='toself',
        fillcolor=rgba_color,
        line=dict(color=rgba_color, width=0),  # No border
        mode='lines',
        hoverinfo='text',
        hovertext=text,
        showlegend=False,
        name=node_id
    )

def create_triangle_trace(x, y, size=10.0, color='#8e44ad', alpha=1.0, text='', node_id=''):
    """Create an upward-pointing triangle shape"""
    w = size / 40.0
    h = size / 40.0
    # Triangle points: top center, bottom left, bottom right
    xs = np.array([x, x - w, x + w, x])
    ys = np.array([y + h, y - h, y - h, y + h])
    # Convert color to rgba
    if color.startswith('#'):
        # Parse hex color
        r = int(color[1:3], 16)
        g = int(color[3:5], 16)
        b = int(color[5:7], 16)
        rgba = f'rgba({r},{g},{b},{alpha})'
    else:
        rgba = color
    return go.Scatter(
        x=xs, y=ys,
        fill='toself',
        fillcolor=rgba,
        line=dict(color=rgba, width=0),
        mode='lines',
        hoverinfo='text',
        hovertext=text,
        showlegend=False,
        name=node_id
    )

def create_trapezium_trace(x, y, size=10, color='#888888', alpha=1.0, text='', node_id=''):
    # Simple trapezium centered at x,y
    w = size / 40.0
    h = size / 40.0
    top_w = w * 0.6
    xs = np.array([x - top_w, x + top_w, x + w, x - w, x - top_w])
    ys = np.array([y + h, y + h, y - h, y - h, y + h])
    rgba = f'rgba(136,136,136,{alpha})'
    return go.Scatter(x=xs, y=ys, fill='toself', fillcolor=rgba, line=dict(color=rgba, width=0), mode='lines', hoverinfo='text', hovertext=text, showlegend=False, name=node_id)

def create_rectangle_trace(x, y, size=10, color='#B2B2B2', alpha=1.0, text='', node_id=''):
    w = size / 40.0
    h = size / 40.0
    xs = np.array([x - w, x + w, x + w, x - w, x - w])
    ys = np.array([y + h, y + h, y - h, y - h, y + h])
    rgba = f'rgba(178,178,178,{alpha})'
    return go.Scatter(x=xs, y=ys, fill='toself', fillcolor=rgba, line=dict(color=rgba, width=0), mode='lines', hoverinfo='text', hovertext=text, showlegend=False, name=node_id)

def create_hexagon_trace(x, y, size: float = 10.0, color: str = '#8e44ad', alpha: float = 1.0, text: str = '', node_id: str = ''):
    # regular hexagon
    r = size / 40.0
    angles = np.linspace(0, 2 * np.pi, 7)[:-1] + np.pi/6
    xs = x + r * np.cos(angles)
    ys = y + r * np.sin(angles)
    # Convert color to rgba
    if color.startswith('#'):
        # Parse hex color
        r_val = int(color[1:3], 16)
        g_val = int(color[3:5], 16)
        b_val = int(color[5:7], 16)
        rgba = f'rgba({r_val},{g_val},{b_val},{alpha})'
    elif color == 'lightblue':
        rgba = f'rgba(173,216,230,{alpha})'
    else:
        rgba = color
    return go.Scatter(x=xs, y=ys, fill='toself', fillcolor=rgba, line=dict(color=rgba, width=0), mode='lines', hoverinfo='text', hovertext=text, showlegend=False, name=node_id)

def create_interactive_metabolite_pathway_plot(df, pos_set_zscore=1, neg_set_zscore=-1, 
                                              max_pvalue=1.0, min_metabolites=2,
                                              max_total_pathways=None, max_pos=None, max_neg=None,
                                              node_spacing=1.5, layout_iterations=100,
                                              include_enzymes=True, include_transporters=False, include_diseases=False,
                                              pathway_stats=None):
    """Create interactive plot for metabolite-pathway network

    This function now accepts pathway-limiting parameters and forwards them to
    `create_metabolite_pathway_network` so plotted networks respect the same
    max_total_pathways / max_pos / max_neg selection used elsewhere.
    
    Parameters:
    -----------
    include_enzymes : bool, default=True
        Whether to include enzyme nodes and their connections in the network visualization
    include_transporters : bool, default=True
        Whether to include transporter nodes and their connections in the network visualization
    include_diseases : bool, default=True
        Whether to include disease nodes and their connections in the network visualization
    """
    # If the caller has already computed and filtered pathway statistics, delegate
    # to the _from_stats variant to avoid recalculation and ensure consistency.
    if pathway_stats is not None:
        # Map old parameter names to new ones: include_enzymes/transporters -> include_upstream
        return create_interactive_metabolite_pathway_plot_from_stats(
            df,
            pathway_stats,
            pos_set_fc=0,
            neg_set_fc=0,
            node_spacing=node_spacing,
            layout_iterations=layout_iterations,
            include_upstream=(include_enzymes or include_transporters),
            include_diseases=include_diseases
        )

    G, metabolites, pathways = create_metabolite_pathway_network(
        df,
        pos_set_zscore=pos_set_zscore,
        neg_set_zscore=neg_set_zscore,
        max_pvalue=max_pvalue,
        min_metabolites=min_metabolites,
        max_total_pathways=max_total_pathways,
        max_pos=max_pos,
        max_neg=max_neg,
        include_enzymes=include_enzymes,
        include_transporters=include_transporters,
        include_diseases=include_diseases
    )

    # Available layouts
    layouts = ['shell', 'spring', 'circular', 'kamada_kawai', 'bipartite']
    
    fig = go.Figure()
    
    # Create traces for each layout
    for i, layout_name in enumerate(layouts):
        pos = get_layout_positions(G, layout_name, k_spacing=node_spacing, iterations=layout_iterations)
        
        # Add edges grouped by type and color (no neutral)
        edge_types = ['activating', 'inhibiting', 'associated']
        edge_colors = {
            'activating': '#FF0000',  # Red dotted
            'inhibiting': '#3498DB',  # Blue dotted
            # 'neutral': '#808080',   # commented out everywhere per request
            'associated': '#888888'   # Gray for associations
        }
        
        for edge_type in edge_types:
            edge_x = []
            edge_y = []
            edge_weights = []
            # Set thinner lines
            line_width = 2.5
            for edge in G.edges(data=True):
                # skip edges whose nodes do not have layout positions (can happen if nodes are invalid/NaN)
                try:
                    et = edge[2].get('edge_type') if isinstance(edge[2], dict) else edge[2]['edge_type']
                except Exception:
                    continue
                if et == edge_type:
                    if edge[0] not in pos or edge[1] not in pos:
                        # unexpected/missing node position (e.g. NaN); skip this edge
                        continue
                    x0, y0 = pos[edge[0]]
                    x1, y1 = pos[edge[1]]
                    edge_x.extend([x0, x1, None])
                    edge_y.extend([y0, y1, None])
                    edge_weights.append(edge[2].get('weight', 1))
            if edge_x:  # Only add trace if there are edges of this type
                fig.add_trace(go.Scatter(
                    x=edge_x, y=edge_y,
                    line=dict(
                        width=line_width,
                        color=edge_colors[edge_type],
                        dash='8,8'
                    ),
                    hoverinfo='none',
                    mode='lines',
                    visible=(i == 0),
                    name=f'edges_{edge_type}_{layout_name}',
                    showlegend=False
                ))
        
        # Add metabolite nodes (circles)
        metabolite_nodes = [(n, d) for n, d in G.nodes(data=True) if d['node_type'] == 'metabolite' and n in pos]
        
        if metabolite_nodes:
            met_x = [pos[n][0] for n, d in metabolite_nodes]
            met_y = [pos[n][1] for n, d in metabolite_nodes]
            met_colors = []
            for n, d in metabolite_nodes:
                color = d["color"]
                alpha = d["alpha"]
                # Upregulated metabolite should be orange
                if color in ("#FF8C00", "orange", "Orange", "#FF0000", "red", "Red"):
                    rgba_color = f'rgba(255, 140, 0, {alpha})'  # orange
                elif color == "green":
                    rgba_color = f'rgba(0, 128, 0, {alpha})'
                elif color == "red":
                    rgba_color = f'rgba(255, 0, 0, {alpha})'
                elif color == "lightgray":
                    rgba_color = f'rgba(211, 211, 211, {alpha})'
                else:
                    # Default to grey for any other color
                    rgba_color = f'rgba(128, 128, 128, {alpha})'
                met_colors.append(rgba_color)
            met_texts = [wrap_text_smart(n, max_chars_per_line=12, max_lines=2) for n, d in metabolite_nodes]
            met_hover = [f'{n} (FC: {d["log2_fold_change"]:.1f}, {d["regulation"]})' 
                        for n, d in metabolite_nodes]
            
            fig.add_trace(go.Scatter(
                x=met_x, y=met_y,
                mode='markers+text',
                marker=dict(
                    size=15,
                    color=met_colors,
                    line=dict(width=0),  # No border
                    symbol='circle'
                ),
                text=met_texts,
                textposition="middle center",
                textfont=dict(color='black', size=11, family='Arial Black'),
                hoverinfo='text',
                hovertext=met_hover,
                visible=(i == 0),
                name=f'metabolites_{layout_name}',
                showlegend=False
            ))
        
        # Add pathway nodes (diamonds)
        pathway_nodes = [(n, d) for n, d in G.nodes(data=True) if d['node_type'] == 'pathway' and n in pos]
            
        for node, data in pathway_nodes:
                x, y = pos[node]
                diamond_trace = create_diamond_trace(
                    x, y, 
                    size=0.2,  # Small, data-units-based diamond size
                    color=data['color'],
                    alpha=data['alpha'],
                    text=f"{node} ({data['status']}, p={data['combined_pvalue']:.3f})",
                    node_id=node
                )
                diamond_trace.visible = (i == 0)
                diamond_trace.name = f'pathway_{node}_{layout_name}'
                fig.add_trace(diamond_trace)
                
                # Add text labels for pathways
                fig.add_trace(go.Scatter(
                    x=[x], y=[y],
                    mode='text',
                    text=[node.replace(' ', '<br>')],  # Break long pathway names
                    textfont=dict(color='black', size=11, family='Arial Black'),
                    textposition='middle center',
                    hoverinfo='text',
                    hovertext=f'{node} ({data["status"]})',
                    visible=(i == 0),
                    name=f'pathway_text_{node}_{layout_name}',
                    showlegend=False
                ))

        # Add upstream regulator nodes (diamond shape) and show gene label - only if enabled
        # In this function, upstream refers to enzymes/transporters flags
        if (include_enzymes or include_transporters):
            upstream_nodes = [(n, d) for n, d in G.nodes(data=True)
                              if d.get('node_type') in ('upstream', 'enzyme', 'transporter') and n in pos]
            for node, data in upstream_nodes:
                    x, y = pos[node]
                    node_type = data.get('node_type', 'Unknown')
                    gene_label = data.get('textlabel', node)  # Gene symbol for display
                    genes = data.get('genes', [])
                    accessions = data.get('accessions', [])
                    n_metabolites = len(data.get('metabolites', []))
                    
                    # Create detailed hover text with all relevant information
                    hover_text = f"Type: {node_type.capitalize()}<br>"
                    hover_text += f"Gene: {gene_label}<br>"
                    
                    # Add full names for transporters (from transporter_names field)
                    if node_type == 'transporter':
                        transporter_names = data.get('transporter_names', [])
                        if transporter_names:
                            hover_text += f"Names: {', '.join(transporter_names)}<br>"
                    elif node_type == 'enzyme':
                        # For enzymes, use the concat field which has all names
                        concat_info = data.get('concat', '')
                        if concat_info and ':' in concat_info:
                            # Extract enzyme name(s) before first colon
                            enzyme_names = concat_info.split(':')[0].strip()
                            if enzyme_names:
                                hover_text += f"Names: {enzyme_names}<br>"
                    
                    # Add genes if available
                    if genes:
                        hover_text += f"Genes: {', '.join(genes)}<br>"
                    
                    # Add accessions if available
                    if accessions:
                        # Filter out empty accessions
                        valid_accessions = [a for a in accessions if a and str(a).strip()]
                        if valid_accessions:
                            hover_text += f"Accessions: {', '.join(valid_accessions)}<br>"
                    
                    hover_text += f"Metabolites: {n_metabolites}"
                    
                    # Use purple color for all upstream regulators (#8e44ad)
                    trace = create_triangle_trace(x, y, size=0.25, color='#8e44ad', alpha=1.0,
                                                 text=hover_text,
                                                 node_id=node)
                    trace.visible = (i == 0)
                    trace.name = f'upstream_{node}_{layout_name}'
                    fig.add_trace(trace)

                    # Label with gene name
                    label_text = gene_label
                    fig.add_trace(go.Scatter(
                        x=[x], y=[y],
                        mode='text',
                        text=[label_text],
                        textfont=dict(color='black', size=10, family='Arial Black'),
                        hoverinfo='text',
                        hovertext=hover_text,
                        visible=(i == 0),
                        name=f'upstream_text_{node}_{layout_name}',
                        showlegend=False
                    ))

        # Add disease nodes (hexagon) but DO NOT label on network - only if enabled
        if include_diseases:
            disease_nodes = [(n, d) for n, d in G.nodes(data=True) if d['node_type'] == 'disease' and n in pos]
            for node, data in disease_nodes:
                    x, y = pos[node]
                    trace = create_hexagon_trace(x, y, size=0.25, color=data.get('color', '#8e44ad'), alpha=1.0,
                                                 text=f"{node} | Metabolites: {len(data.get('metabolites',[]))}",
                                                 node_id=node)
                    trace.visible = (i == 0)
                    trace.name = f'disease_{node}_{layout_name}'
                    fig.add_trace(trace)
                    # No text label added for diseases per request
    
    # Create dropdown menu for layout selection
    dropdown_buttons = []
    for i, layout_name in enumerate(layouts):
        total_traces = len(list(fig.data))
        traces_per_layout = total_traces // len(layouts)
        start_idx = i * traces_per_layout
        end_idx = min(start_idx + traces_per_layout, total_traces)
        
        visibility = [False] * total_traces
        for j in range(start_idx, end_idx):
            visibility[j] = True
        
        dropdown_buttons.append(
            dict(
                args=[{"visible": visibility}],
                label=layout_name.title(),
                method="restyle"
            )
        )
    
    # Add dummy traces for legend (Dashed Lines and diamonds)
    # NEW: Edge legend based on metabolite regulation
    # Red line for Upregulated metabolite (remove word 'Edge')
    fig.add_trace(go.Scatter(
        x=[None, None], y=[None, None],
        mode='lines',
        line=dict(color='orange', width=2.5, dash='8,8'),
        name='Predicted Activation',
        showlegend=True
    ))
    # Green line for Downregulated metabolite (remove word 'Edge')
    fig.add_trace(go.Scatter(
        x=[None, None], y=[None, None],
        mode='lines',
        line=dict(color='green', width=2.5, dash='8,8'),
        name='Predicted Inhibition',
        showlegend=True
    ))
    # Commented out: gray line for Neutral (not used with current metabolite filtering)
    # fig.add_trace(go.Scatter(
    #     x=[None, None], y=[None, None],
    #     mode='lines',
    #     line=dict(color='#808080', width=2.5, dash='8,8'),
    #     name='Neutral',
    #     showlegend=True
    # ))
    # Red diamond for Activated pathway
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='diamond', size=12, color='orange', line=dict(width=0)),
        name='Activated Pathway',
        showlegend=True
    ))
    # Blue diamond for Inhibited pathway
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='diamond', size=12, color='blue', line=dict(width=0)),
        name='Inhibited Pathway',
        showlegend=True
    ))
    # Orange circle for Upregulated metabolite
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='circle', size=12, color='#FF8C00', line=dict(width=0)),
        name='Upregulated Metabolites',
        showlegend=True
    ))
    # Green circle for Downregulated
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='circle', size=12, color='green', line=dict(width=0)),
        name='Downregulated Metabolites',
        showlegend=True
    ))
    # Light gray circle for No Change
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='diamond', size=12, color='#cccccc', line=dict(width=0)),
        name='No Change',
        showlegend=True
    ))

    # Legend entries for upstream regulators and diseases - only if enabled
    # Upstream refers to enzymes/transporters in this function's API
    if (include_enzymes or include_transporters):
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(symbol='triangle-up', size=12, color='#8e44ad', line=dict(width=0)),
            name='Upstream Regulator',
            showlegend=True
        ))
    
    if include_diseases:
        # hexagon symbol usually available as 'hexagon' or 'hexagon2'
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(symbol='hexagon', size=12, color='#8e44ad', line=dict(width=0)),
            name='Associated Disease',
            showlegend=True
        ))

    # Legend for associations - only if upstream or diseases enabled
    if (include_enzymes or include_transporters) or include_diseases:
        fig.add_trace(go.Scatter(
            x=[None, None], y=[None, None],
            mode='lines',
            line=dict(color='#888888', width=2.5, dash='8,8'),  # Gray for associations
            name='Association',
            showlegend=True
        ))

    # Update layout
    fig.update_layout(
        title={
            'text': "Metabolite-Pathway Interaction Network<br><sub>Circles: Metabolites | Diamonds: Pathways | Dashed Lines: Interactions</sub>",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 16}
        },
        updatemenus=[
            dict(
                buttons=dropdown_buttons,
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.1,
                xanchor="left",
                y=1.15,
                yanchor="top",
            ),
        ],
        annotations=[
            dict(
                text="Layout:",
                x=0.02, y=1.12,
                xref="paper", yref="paper",
                align="left",
                showarrow=False,
                font=dict(size=12)
            )
        ],
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
            font=dict(size=14)
        ),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, constrain='domain'),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, constrain='domain'),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=1600,
        height=1000,
        margin=dict(l=150, r=400, t=150, b=150),
        autosize=True
    )
    
    return fig

def wrap_text_smart(text, max_chars_per_line=15, max_lines=3):
    """Smart text wrapping that tries to avoid breaking unless absolutely necessary"""
    # Don't wrap short text at all
    if len(text) <= max_chars_per_line:
        return text
    
    # Be more generous with longer text - only wrap if significantly exceeding limits
    # This helps prevent unnecessary text breaking that makes names harder to read
    generous_limit = max_chars_per_line * 1.8  # Allow 80% overflow before breaking
    
    if len(text) <= generous_limit:
        return text  # Keep as single line even if a bit long
    
    # Only break very long pathway names
    words = text.split(' ')
    if len(words) <= 1:
        # Single word (metabolite names, etc.) - avoid breaking unless extremely long
        if len(text) <= max_chars_per_line * 2.5:  # Very generous for single words
            return text
        else:
            # Only break extremely long single words
            mid_point = len(text) // 2
            return text[:mid_point] + '<br>' + text[mid_point:]
    
    # Multiple words - break intelligently only for very long names
    lines = []
    current_line = ""
    
    for word in words:
        test_line = current_line + ' ' + word if current_line else word
        # Use generous limit for breaking decisions too
        if len(test_line) <= generous_limit:
            current_line = test_line
        else:
            if current_line:
                lines.append(current_line)
                current_line = word
            else:
                lines.append(word)
                current_line = ""
            
            if len(lines) >= max_lines:
                break
    
    if current_line and len(lines) < max_lines:
        lines.append(current_line)
    
    # If we have multiple lines, return them joined with <br>
    if len(lines) > 1:
        return '<br>'.join(lines[:max_lines])
    elif len(lines) == 1:
        return lines[0]
    else:
        return text  # Fallback to original text
    
def create_interactive_metabolite_pathway_plot_from_stats(df, pathway_stats, 
                                                          pos_set_fc=0, neg_set_fc=0,
                                                          node_spacing=1.5, layout_iterations=100,
                                                          metabolite_font_size=14, pathway_font_size=12, enzyme_font_size=10,
                                                          metabolite_max_chars=12, metabolite_max_lines=2,
                                                          pathway_max_chars=15, pathway_max_lines=3,
                                                          network_width=1400, network_height=800,
                                                          include_upstream=True, include_diseases=False,
                                                          upstream_data=None, disease_data=None,
                                                          compress_by_class=True):
    """Create interactive plot for metabolite-pathway network using pre-filtered pathway stats
    
    This function creates the network from already filtered pathway statistics to ensure
    the visualization shows exactly the pathways that met the user's criteria.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        The original metabolomics data 
    pathway_stats : dict
        Already filtered pathway statistics from calculate_pathway_statistics
    pos_set_fc : float, default=0
        Positive fold change threshold for determining upregulation
    neg_set_fc : float, default=0
        Negative fold change threshold for determining downregulation
    node_spacing : float, default=1.5
        Controls spacing between nodes (lower = more compact, higher = more spread out)
        Recommended values: 0.5-3.0
    layout_iterations : int, default=100
        Number of iterations for spring layout (higher = better positioning)
    metabolite_font_size : int, default=10
        Font size for metabolite labels (recommended: 8-16)
    pathway_font_size : int, default=9
        Font size for pathway labels (recommended: 8-14)
    enzyme_font_size : int, default=8
        Font size for enzyme and transporter labels (recommended: 6-12)
    metabolite_max_chars : int, default=12
        Maximum characters per line for metabolite text wrapping (recommended: 8-20)
    metabolite_max_lines : int, default=2
        Maximum lines for metabolite text wrapping (recommended: 1-3)
    pathway_max_chars : int, default=15
        Maximum characters per line for pathway text wrapping (recommended: 10-25)
    pathway_max_lines : int, default=3
        Maximum lines for pathway text wrapping (recommended: 2-4)
    include_upstream : bool, default=True
        Whether to include upstream regulator nodes (enzymes/transporters) and their connections
    include_diseases : bool, default=False
        Whether to include disease nodes and their connections in the network visualization
    upstream_data : pandas.DataFrame, optional
        DataFrame containing upstream regulator data with columns: Name, Gene, Type, Associated_Metabolites, # Metabolites
    disease_data : pandas.DataFrame, optional
        DataFrame containing disease data with columns: Disease, Associated_Metabolites, # Metabolites
    compress_by_class : bool, default=True
        If True and Main_Class column exists, compress lipids by class (show one node per class)
        If False, show all individual features
    """
    
    # Create network directly from filtered pathway stats instead of re-filtering
    G, metabolites, pathways = create_metabolite_pathway_network_from_stats(
        df, pathway_stats,
        pos_set_fc=pos_set_fc,
        neg_set_fc=neg_set_fc,
        include_upstream=include_upstream,
        include_diseases=include_diseases,
        upstream_data=upstream_data,
        disease_data=disease_data
    )

    # Calculate dynamic margins based on maximum text lengths
    # Get all metabolite and pathway names from the network
    metabolite_names = [n for n, d in G.nodes(data=True) if d['node_type'] == 'metabolite']
    pathway_names = [n for n, d in G.nodes(data=True) if d['node_type'] == 'pathway']
    
    # Calculate maximum text lengths after wrapping
    max_metabolite_length = 0
    max_pathway_length = 0
    
    if metabolite_names:
        for name in metabolite_names:
            wrapped_text = wrap_text_smart(name, max_chars_per_line=metabolite_max_chars, max_lines=metabolite_max_lines)
            # Count longest line in wrapped text
            lines = wrapped_text.split('<br>')
            max_line_length = max(len(line) for line in lines)
            max_metabolite_length = max(max_metabolite_length, max_line_length)
    
    if pathway_names:
        for name in pathway_names:
            wrapped_text = wrap_text_smart(name, max_chars_per_line=pathway_max_chars, max_lines=pathway_max_lines)
            # Count longest line in wrapped text
            lines = wrapped_text.split('<br>')
            max_line_length = max(len(line) for line in lines)
            max_pathway_length = max(max_pathway_length, max_line_length)
    
    # Calculate generous margins to prevent any text cutoff
    # Reduced base margins to use more space for actual content
    base_margin_left = 80      # Reduced from 250 to 80
    base_margin_right = 150    # Reduced from 500 to 150 (still space for legend)
    base_margin_top = 60       # Reduced from 180 to 60
    base_margin_bottom = 60    # Reduced from 180 to 60
    
    # Add moderate margin based on text length
    char_width_estimate = 8   # More generous character width estimate
    max_text_length = max(max_metabolite_length, max_pathway_length)
    text_margin = min(100, max_text_length * char_width_estimate)  # Reduced cap from 300px to 100px
    
    dynamic_left = base_margin_left + text_margin // 2      # Reduced text margin influence
    dynamic_right = base_margin_right + text_margin // 4    # Even less margin on right
    dynamic_top = base_margin_top + (metabolite_max_lines + pathway_max_lines) * 10  # Reduced line spacing from 20 to 10
    dynamic_bottom = base_margin_bottom + text_margin // 4  # Less margin on bottom
    
    print(f"📐 Dynamic margin calculation:")
    print(f"   - Max metabolite text length: {max_metabolite_length} chars")
    print(f"   - Max pathway text length: {max_pathway_length} chars")
    print(f"   - Calculated margins: l={dynamic_left}, r={dynamic_right}, t={dynamic_top}, b={dynamic_bottom}")

    # Available layouts
    layouts = ['shell', 'spring', 'circular', 'kamada_kawai', 'bipartite']
    
    fig = go.Figure()
    
    # Create traces for each layout
    for i, layout_name in enumerate(layouts):
        pos = get_layout_positions(G, layout_name, k_spacing=node_spacing, iterations=layout_iterations)
        
        # Adjust positions to prevent text cutoff
        pos = adjust_positions_for_text_length(G, pos, 
                                              metabolite_max_chars=metabolite_max_chars,
                                              pathway_max_chars=pathway_max_chars,
                                              metabolite_max_lines=metabolite_max_lines,
                                              pathway_max_lines=pathway_max_lines)
        
        # Add edges grouped by type and color (include 'associated' edges for upstream/disease links)
        # Edge types: 'activating', 'inhibiting', 'associated' (must match edge creation in network builder)
        edge_types = ['activating', 'inhibiting', 'associated']
        edge_colors = {
            'activating': '#FF0000',   # Red for activating edges (upregulated metabolites)
            'inhibiting': '#3498DB',   # Blue for inhibiting edges (downregulated metabolites)
            'neutral': "#808080",      # Neutral (not currently used)
            'associated': '#888888'    # Gray for upstream/disease associations
        }
        
        for edge_type in edge_types:
            edge_x = []
            edge_y = []
            edge_weights = []
            # Set thinner lines
            line_width = 2.5
            for edge in G.edges(data=True):
                # skip edges whose nodes do not have layout positions (can happen if nodes are invalid/NaN)
                try:
                    et = edge[2].get('edge_type') if isinstance(edge[2], dict) else edge[2]['edge_type']
                except Exception:
                    continue
                if et == edge_type:
                    if edge[0] not in pos or edge[1] not in pos:
                        # unexpected/missing node position (e.g. NaN); skip this edge
                        continue
                    x0, y0 = pos[edge[0]]
                    x1, y1 = pos[edge[1]]
                    edge_x.extend([x0, x1, None])
                    edge_y.extend([y0, y1, None])
                    edge_weights.append(edge[2].get('weight', 1))
            if edge_x:  # Only add trace if there are edges of this type
                fig.add_trace(go.Scatter(
                    x=edge_x, y=edge_y,
                    line=dict(
                        width=line_width,
                        color=edge_colors[edge_type],
                        dash='8,8'
                    ),
                    hoverinfo='none',
                    mode='lines',
                    visible=(i == 0),
                    name=f'edges_{edge_type}_{layout_name}',
                    showlegend=False
                ))
        
        # Add metabolite nodes (circles)
        metabolite_nodes = [(n, d) for n, d in G.nodes(data=True) if d['node_type'] == 'metabolite' and n in pos]
        
        if metabolite_nodes:
            met_x = [pos[n][0] for n, d in metabolite_nodes]
            met_y = [pos[n][1] for n, d in metabolite_nodes]
            met_colors = []
            for n, d in metabolite_nodes:
                color = d["color"]
                alpha = d["alpha"]
                # Upregulated metabolite should be orange
                if color in ("#FF8C00", "orange", "Orange", "#FF0000", "red", "Red"):
                    rgba_color = f'rgba(255, 140, 0, {alpha})'  # orange
                elif color == "green":
                    rgba_color = f'rgba(0, 128, 0, {alpha})'
                elif color == "red":
                    rgba_color = f'rgba(255, 0, 0, {alpha})'
                elif color == "lightgray":
                    rgba_color = f'rgba(211, 211, 211, {alpha})'
                else:
                    # Default to grey for any other color
                    rgba_color = f'rgba(128, 128, 128, {alpha})'
                met_colors.append(rgba_color)
            met_texts = [wrap_text_smart(n, max_chars_per_line=metabolite_max_chars, max_lines=metabolite_max_lines) for n, d in metabolite_nodes]
            met_hover = [f'{n} (FC: {d["log2_fold_change"]:.1f}, {d["regulation"]})' 
                        for n, d in metabolite_nodes]
            
            fig.add_trace(go.Scatter(
                x=met_x, y=met_y,
                mode='markers+text',
                marker=dict(
                    size=15,
                    color=met_colors,
                    line=dict(width=0),  # No border
                    symbol='circle'
                ),
                text=met_texts,
                textposition="middle center",
                textfont=dict(color='black', size=metabolite_font_size, family='Arial Black'),
                hoverinfo='text',
                hovertext=met_hover,
                visible=(i == 0),
                name=f'metabolites_{layout_name}',
                showlegend=False
            ))
        
        # Add pathway nodes (diamonds)
        pathway_nodes = [(n, d) for n, d in G.nodes(data=True) if d['node_type'] == 'pathway' and n in pos]
            
        for node, data in pathway_nodes:
                x, y = pos[node]
                diamond_trace = create_diamond_trace(
                    x, y, 
                    size=0.2,  # Small, data-units-based diamond size
                    color=data['color'],
                    alpha=data['alpha'],
                    text=f"{node} ({data['status']}, p={data['combined_pvalue']:.3f})",
                    node_id=node
                )
                diamond_trace.visible = (i == 0)
                diamond_trace.name = f'pathway_{node}_{layout_name}'
                fig.add_trace(diamond_trace)
                
                # Add text labels for pathways
                fig.add_trace(go.Scatter(
                    x=[x], y=[y],
                    mode='text',
                    text=[wrap_text_smart(node, max_chars_per_line=pathway_max_chars, max_lines=pathway_max_lines)],  # Break long pathway names
                    textfont=dict(color='black', size=pathway_font_size, family='Arial Black'),
                    hoverinfo='text',
                    hovertext=f'{node} ({data["status"]})',
                    visible=(i == 0),
                    name=f'pathway_text_{node}_{layout_name}',
                    showlegend=False
                ))

        # Add upstream regulator nodes (purple diamonds) with labels - only if enabled
        if include_upstream:
            upstream_nodes = [(n, d) for n, d in G.nodes(data=True) if d.get('node_type') == 'upstream' and n in pos]
            for node, data in upstream_nodes:
                x, y = pos[node]
                color = data.get('color', '#8e44ad')
                node_type = data.get('upstream_type', 'Upstream')
                genes = data.get('genes', [])
                accessions = data.get('accessions', [])
                n_metabolites = len(data.get('metabolites', []))
                
                # Create detailed hover text
                hover_text = f"{node} ({node_type})"
                if genes:
                    hover_text += f"<br>Genes: {', '.join(genes)}"
                if accessions:
                    valid_accessions = [a for a in accessions if a and str(a).strip()]
                    if valid_accessions:
                        hover_text += f"<br>Accessions: {', '.join(valid_accessions)}"
                hover_text += f"<br>Metabolites: {n_metabolites}"
                
                # Draw upstream node as triangle
                u_trace = create_triangle_trace(
                    x, y,
                    size=0.18,
                    color=color,
                    alpha=1.0,
                    text=hover_text,
                    node_id=node
                )
                u_trace.visible = (i == 0)
                u_trace.name = f'upstream_{node}_{layout_name}'
                fig.add_trace(u_trace)

                # Label upstream node by gene/name
                fig.add_trace(go.Scatter(
                    x=[x], y=[y],
                    mode='text',
                    text=[wrap_text_smart(str(data.get('gene', node)), max_chars_per_line=metabolite_max_chars, max_lines=metabolite_max_lines)],
                    textfont=dict(color='black', size=enzyme_font_size, family='Arial Black'),
                    hoverinfo='text',
                    hovertext=f"{node} ({data.get('upstream_type','Upstream')})",
                    visible=(i == 0),
                    name=f'upstream_text_{node}_{layout_name}',
                    showlegend=False
                ))

        # Add disease nodes (hexagon) and label by name - only if enabled
        if include_diseases:
            disease_nodes = [(n, d) for n, d in G.nodes(data=True) if d['node_type'] == 'disease' and n in pos]
            for node, data in disease_nodes:
                    x, y = pos[node]
                    trace = create_hexagon_trace(x, y, size=0.25, color=data.get('color', '#8e44ad'), alpha=1.0,
                                                 text=f"{node} | Metabolites: {len(data.get('metabolites',[]))}",
                                                 node_id=node)
                    trace.visible = (i == 0)
                    trace.name = f'disease_{node}_{layout_name}'
                    fig.add_trace(trace)
                    # Add disease name label
                    fig.add_trace(go.Scatter(
                        x=[x], y=[y],
                        mode='text',
                        text=[wrap_text_smart(node, max_chars_per_line=pathway_max_chars, max_lines=pathway_max_lines)],
                        textfont=dict(color='black', size=pathway_font_size, family='Arial Black'),
                        hoverinfo='text',
                        hovertext=node,
                        visible=(i == 0),
                        name=f'disease_text_{node}_{layout_name}',
                        showlegend=False
                    ))
    
    # Create dropdown menu for layout selection
    dropdown_buttons = []
    for i, layout_name in enumerate(layouts):
        total_traces = len(list(fig.data))
        traces_per_layout = total_traces // len(layouts)
        start_idx = i * traces_per_layout
        end_idx = min(start_idx + traces_per_layout, total_traces)
        
        visibility = [False] * total_traces
        for j in range(start_idx, end_idx):
            visibility[j] = True
        
        dropdown_buttons.append(
            dict(
                args=[{"visible": visibility}],
                label=layout_name.title(),
                method="restyle"
            )
        )
    
    # Add dummy traces for legend (Dashed Lines and diamonds)
    # NEW: Edge legend based on metabolite regulation
    # Red line for Upregulated metabolite (remove word 'Edge')
    fig.add_trace(go.Scatter(
        x=[None, None], y=[None, None],
        mode='lines',
        line=dict(color='orange', width=2.5, dash='8,8'),
        name='Predicted Activation',
        showlegend=True
    ))
    # Green line for Downregulated metabolite (remove word 'Edge')
    fig.add_trace(go.Scatter(
        x=[None, None], y=[None, None],
        mode='lines',
        line=dict(color='green', width=2.5, dash='8,8'),
        name='Predicted Inhibition',
        showlegend=True
    ))
    # Commented out: gray line for Neutral (not used with current metabolite filtering)
    # fig.add_trace(go.Scatter(
    #     x=[None, None], y=[None, None],
    #     mode='lines',
    #     line=dict(color='#808080', width=2.5, dash='8,8'),
    #     name='Neutral',
    #     showlegend=True
    # ))
    # Red diamond for Activated pathway
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='diamond', size=12, color='orange', line=dict(width=0)),
        name='Activated Pathway',
        showlegend=True
    ))
    # Blue diamond for Inhibited pathway
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='diamond', size=12, color='blue', line=dict(width=0)),
        name='Inhibited Pathway',
        showlegend=True
    ))
    # Orange circle for Upregulated metabolite
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='circle', size=12, color='#FF8C00', line=dict(width=0)),
        name='Upregulated Metabolites',
        showlegend=True
    ))
    # Green circle for Downregulated
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='circle', size=12, color='green', line=dict(width=0)),
        name='Downregulated Metabolites',
        showlegend=True
    ))
    # Light gray circle for No Change
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode='markers',
        marker=dict(symbol='diamond', size=12, color='#cccccc', line=dict(width=0)),
        name='No Change',
        showlegend=True
    ))

    # Legend entries for upstream regulators and diseases - only if enabled
    if include_upstream:
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(symbol='triangle-up', size=12, color='#8e44ad', line=dict(width=0)),
            name='Upstream Regulator',
            showlegend=True
        ))
    
    if include_diseases:
        # hexagon symbol usually available as 'hexagon' or 'hexagon2'
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(symbol='hexagon', size=12, color='#8e44ad', line=dict(width=0)),
            name='Associated Disease',
            showlegend=True
        ))

    # Legend for associations - only if upstream or diseases enabled
    if include_upstream or include_diseases:
        fig.add_trace(go.Scatter(
            x=[None, None], y=[None, None],
            mode='lines',
            line=dict(color='#888888', width=2.5, dash='8,8'),  # Gray for associations
            name='Association',
            showlegend=True
        ))

    # Update layout
    fig.update_layout(
        title={
            'text': "Metabolite-Pathway Interaction Network<br><sub>Circles: Metabolites | Diamonds: Pathways | Dashed Lines: Interactions</sub>",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 16}
        },
        updatemenus=[
            dict(
                buttons=dropdown_buttons,
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.1,
                xanchor="left",
                y=1.15,
                yanchor="top",
            ),
        ],
        annotations=[
            dict(
                text="Layout:",
                x=0.02, y=1.12,
                xref="paper", yref="paper",
                align="left",
                showarrow=False,
                font=dict(size=12)
            )
        ],
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
            font=dict(size=14)
        ),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=network_width,
        height=network_height,
        margin=dict(l=dynamic_left, r=dynamic_right, t=dynamic_top, b=dynamic_bottom)
    )
    
    return fig

def create_pathway_summary_plots(
    metabolites,
    pathways,
    font_size_scale=1.0,
    bargraph_pathway_max_chars=40,
    bargraph_pathway_max_lines=3,
    bargraph_title_font_size=22,
    bargraph_axis_title_font_size=20,
    bargraph_tick_font_size=15,
    bargraph_pathway_label_font_size=20,
    bargraph_width_px=1200,
    bargraph_height_px=700,
    bargraph_legend_title_font_size=14,
    bargraph_legend_tick_font_size=12,
    bargraph_legend_thickness=10,
    bargraph_legend_length_percent=70,
    chord_metabolite_max_chars=40,
    chord_metabolite_max_lines=2,
    chord_pathway_max_chars=40,
    chord_pathway_max_lines=2,
    chord_metabolite_font_size=12,
    chord_pathway_font_size=12,
    chord_width_inches=12,
    chord_height_inches=12,
    chord_keep_metabolite_case=False,
):
    """Create individual bar plots and bubble plots for pathway analysis
    
    Parameters:
    -----------
    metabolites : dict
        Metabolite data
    pathways : dict
        Pathway data
    font_size_scale : float, optional (default=1.0)
        Multiplier for all font sizes. Use 1.5 for 50% larger, 2.0 for 2x larger, etc.
    bargraph_pathway_max_chars : int, optional (default=40)
        Maximum characters per line for pathway labels in the bar graphs.
    bargraph_pathway_max_lines : int, optional (default=3)
        Maximum wrapped lines for pathway labels in the bar graphs.
    bargraph_title_font_size : int, optional (default=22)
        Title font size for enrichment bar graphs.
    bargraph_axis_title_font_size : int, optional (default=20)
        Axis-title font size for enrichment bar graphs.
    bargraph_tick_font_size : int, optional (default=15)
        Tick-label font size for x-axis values in enrichment bar graphs.
    bargraph_pathway_label_font_size : int, optional (default=20)
        Tick-label font size for pathway labels (y-axis) in enrichment bar graphs.
    bargraph_width_px : int, optional (default=1200)
        Minimum width for enrichment bar graph figures in pixels.
    bargraph_height_px : int, optional (default=700)
        Minimum height for enrichment bar graph figures in pixels.
    bargraph_legend_title_font_size : int, optional (default=14)
        Colorbar title font size for enrichment bar graphs.
    bargraph_legend_tick_font_size : int, optional (default=12)
        Colorbar tick font size for enrichment bar graphs.
    bargraph_legend_thickness : int, optional (default=10)
        Colorbar thickness for enrichment bar graphs.
    bargraph_legend_length_percent : int, optional (default=70)
        Colorbar length as a percent of axis height for enrichment bar graphs.
    chord_metabolite_max_chars : int, optional (default=40)
        Maximum characters per line for metabolite labels in the chord diagram.
    chord_metabolite_max_lines : int, optional (default=2)
        Maximum wrapped lines for metabolite labels in the chord diagram.
    chord_pathway_max_chars : int, optional (default=40)
        Maximum characters per line for pathway labels in the chord diagram.
    chord_pathway_max_lines : int, optional (default=2)
        Maximum wrapped lines for pathway labels in the chord diagram.
    chord_metabolite_font_size : int, optional (default=12)
        Font size for metabolite labels in the chord diagram.
    chord_pathway_font_size : int, optional (default=12)
        Font size for pathway labels in the chord diagram.
    chord_width_inches : float, optional (default=12)
        Chord diagram figure width in inches.
    chord_height_inches : float, optional (default=12)
        Chord diagram figure height in inches.
    chord_keep_metabolite_case : bool, optional (default=False)
        If True, preserve original metabolite/compound label casing in the chord diagram.
        If False, convert metabolite/compound labels to sentence case.
    """
    
    # Extract neutral_z_band from metadata (use 0.5 as default if not found)
    neutral_z_band = pathways.get('_neutral_z_band', 0.5)
    
    # Convert pathways dict to list for easier processing, and compute z_score
    pathway_list = []
    
    # DEBUG: Check for 'hits' field in input data
    print("\n" + "="*70)
    print("📊 DEBUG: CHECKING 'hits' FIELD IN PATHWAY DATA FOR PLOTS")
    print("="*70)
    
    for idx, (name, data) in enumerate(pathways.items(), 1):
        # Skip metadata keys (start with underscore)
        if name.startswith('_'):
            continue
            
        # DEBUG: Check what fields are available
        if idx <= 3:  # Show first 3 pathways
            has_hits = 'hits' in data if isinstance(data, dict) else False
            has_metabolites = 'metabolites' in data if isinstance(data, dict) else False
            hits_count = len(data.get('hits', [])) if has_hits else 0
            mets_count = len(data.get('metabolites', [])) if has_metabolites else 0
            print(f"  [{idx}] {name[:50]}")
            print(f"      has 'hits': {has_hits} (count: {hits_count})")
            print(f"      has 'metabolites': {has_metabolites} (count: {mets_count})")
        
        # Calculate average fold change for metabolites in this pathway
        # CRITICAL: Use 'hits' (significant metabolites) if available, otherwise 'metabolites'
        hits_list = data.get('hits', data.get('metabolites', []))
        log2_fold_changes = [metabolites[met]['log2FC'] for met in hits_list if met in metabolites]
        avg_fc = np.mean(log2_fold_changes) if log2_fold_changes else 0
        # Use z_score from pathway data (this is what determines Activated/Inhibited status)
        z_score = data.get('z_score', 0)
        pathway_list.append({
            'name': name,
            'status': data['status'],
            'combined_pvalue': data['combined_pvalue'],
            'neg_log_pvalue': -np.log10(data['combined_pvalue']) if data['combined_pvalue'] > 0 else 0,
            'metabolites': data.get('metabolites', []),  # Keep original for reference
            'hits': hits_list,  # CRITICAL: Include hits for bar charts
            'avg_log2_fold_change': avg_fc,
            'z_score': z_score,
            'bubble_number': idx
        })
    
    print("="*70 + "\n")
    
    # Calculate dynamic dimensions based on number of pathways
    num_pathways = len(pathway_list)
    
    # For vertical bar plots: width scales with pathway count
    # Base: 600px, add 80px per pathway, min 600, max 2000
    vertical_width = max(600, min(2000, 600 + num_pathways * 80))
    vertical_height = 600  # Fixed height for vertical plots
    
    # For horizontal bar plots: height scales with pathway count
    # Base: 400px, add 60px per pathway, min 400, max 1500
    horizontal_width = 900  # Fixed width for horizontal plots
    horizontal_height = max(400, min(1500, 400 + num_pathways * 60))
    
    # For bubble plot: adjust based on pathway count
    bubble_width = max(900, min(1600, 900 + num_pathways * 40))
    bubble_height = max(600, min(1000, 600 + num_pathways * 20))

    # Figure 1: Vertical Bar plot for Pathway States (using z-score)
    # STATISTICS PLOTS COLOR POLICY: Only two colors, no fading
    # - Red for Activated (z > 0)
    # - Blue for Inhibited (z < 0)
    # - Optional gray for exactly zero
    def get_color_by_zscore(z_score, status, neutral_z_band):
        if z_score > 0:
            return '#FF4444'  # Red
        elif z_score < 0:
            return '#3A7CFF'  # Blue
        else:
            return '#BDC3C7'  # Gray (rare)
    
    # Generate colors based on z-score intensity (use same neutral_z_band from classification)
    pathway_colors = [get_color_by_zscore(p['z_score'], p['status'], neutral_z_band) for p in pathway_list]
    
    # Debug: Print z-scores and colors for verification
    print("\n=== PATHWAY COLOR DEBUG ===")
    print(f"neutral_z_band: {neutral_z_band}")
    for i, p in enumerate(pathway_list[:5]):  # Show first 5
        print(f"{p['name'][:50]}: z={p['z_score']:.3f}, status='{p['status']}', color={pathway_colors[i]}")
    print("===========================\n")
    
    fig1 = go.Figure()
    
    fig1.add_trace(
        go.Bar(
            x=[wrap_pathway_name(p['name'], max_chars=bargraph_pathway_max_chars, max_lines=bargraph_pathway_max_lines) for p in pathway_list],
            y=[p['z_score'] for p in pathway_list],
            orientation='v',
            marker=dict(color=pathway_colors),
            showlegend=False
        )
    )
    
    fig1.update_layout(
        title={
            'text': "Pathway States by Significance",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': int(16 * font_size_scale)}
        },
        yaxis=dict(
            title=dict(text="z-score", font=dict(weight='bold', size=int(14 * font_size_scale))),
            showline=True,
            linewidth=2,
            linecolor='black',
            tickfont=dict(weight='bold', size=int(12 * font_size_scale)),
            tickformat='.1f',
            zeroline=True,
            zerolinecolor='black',
            zerolinewidth=2,
            range=[min(-3, min([p['z_score'] for p in pathway_list], default=-3) * 1.1),
                   max(3, max([p['z_score'] for p in pathway_list], default=3) * 1.1)]
        ),
        xaxis=dict(
            title=dict(text="Pathways", font=dict(weight='bold', size=int(14 * font_size_scale))),
            tickangle=-45,
            showline=False,
            linewidth=2,
            linecolor='black',
            tickfont=dict(weight='bold', size=int(12 * font_size_scale)),
            mirror=False
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=vertical_width,
        height=vertical_height,
        margin=dict(l=80, r=50, t=80, b=150),
        bargap=0.3  # Add gap between bars (0.0 = no gap, 1.0 = maximum gap)
    )
    # Legend: only two entries (Activated/Red, Inhibited/Blue)
    fig1.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                              marker=dict(color='#FF4444', size=10),
                              name='Activated (Red)', showlegend=True))
    fig1.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                              marker=dict(color='#3A7CFF', size=10),
                              name='Inhibited (Blue)', showlegend=True))
    fig1.update_layout(legend=dict(orientation='v', yanchor='top', y=1, xanchor='left', x=1.02, font=dict(size=int(12 * font_size_scale))))
    
    # Figure 2: Pathway Significance using -log(p-value) - Vertical Bar Plot
    fig2 = go.Figure()

    # Ensure colors are mapped to the exact bar order
    fig2_names = [wrap_pathway_name(p['name'], max_chars=bargraph_pathway_max_chars, max_lines=bargraph_pathway_max_lines) for p in pathway_list]
    fig2_values = [p['neg_log_pvalue'] for p in pathway_list]
    fig2_colors = [get_color_by_zscore(p['z_score'], p['status'], neutral_z_band) for p in pathway_list]

    # Debug mapping for fig2
    try:
        print("[DEBUG] Fig2 names count:", len(fig2_names), "values count:", len(fig2_values), "colors count:", len(fig2_colors))
        print("[DEBUG] Fig2 first 5:")
        for i in range(min(5, len(fig2_names))):
            print(f"   {i+1}. {fig2_names[i]} | -log10p={fig2_values[i]:.3f} | color={fig2_colors[i]}")
    except Exception as _e:
        print("[DEBUG] Fig2 print failed:", _e)

    fig2.add_trace(
        go.Bar(
            x=fig2_names,
            y=fig2_values,
            marker=dict(color=fig2_colors),
            showlegend=False
        )
    )
    
    fig2.update_layout(
        title={
            'text': "Pathway Significance (-log p-value) - Vertical",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': int(16 * font_size_scale)}
        },
        xaxis=dict(
            title=dict(text="Pathways", font=dict(weight='bold', size=int(14 * font_size_scale))),
            tickangle=45,
            showline=True,
            linewidth=2,
            linecolor='black',
            categoryorder='array',
            categoryarray=fig2_names,
            tickfont=dict(weight='bold', size=int(12 * font_size_scale))
        ),
        yaxis=dict(
            title=dict(text="-log10(p-value)", font=dict(weight='bold', size=int(14 * font_size_scale))),
            showline=True,
            linewidth=2,
            linecolor='black',
            tickfont=dict(weight='bold', size=int(12 * font_size_scale)),
            tickformat='.1f',
            range=[0, max(5, max(fig2_values, default=5) * 1.1)]
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=vertical_width,
        height=vertical_height,
        margin=dict(l=100, r=50, t=80, b=150),
        bargap=0.3  # Add gap between bars
    )
    # Legend: only two entries (Activated/Red, Inhibited/Blue)
    fig2.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                              marker=dict(color='#FF4444', size=10),
                              name='Activated (Red)', showlegend=True))
    fig2.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                              marker=dict(color='#3A7CFF', size=10),
                              name='Inhibited (Blue)', showlegend=True))
    fig2.update_layout(showlegend=True,
                       legend=dict(orientation='v', yanchor='top', y=1, xanchor='left', x=1.02, font=dict(size=int(12 * font_size_scale))))
    
    # Figure 3: Pathway Significance using -log(p-value) - Horizontal Bar Plot
    fig3 = go.Figure()

    # Ensure colors are mapped to the exact bar order
    fig3_names = [wrap_pathway_name(p['name'], max_chars=bargraph_pathway_max_chars, max_lines=bargraph_pathway_max_lines) for p in pathway_list]
    fig3_values = [p['neg_log_pvalue'] for p in pathway_list]
    fig3_colors = [get_color_by_zscore(p['z_score'], p['status'], neutral_z_band) for p in pathway_list]

    # Debug mapping for fig3
    try:
        print("[DEBUG] Fig3 names count:", len(fig3_names), "values count:", len(fig3_values), "colors count:", len(fig3_colors))
        print("[DEBUG] Fig3 first 5:")
        for i in range(min(5, len(fig3_names))):
            print(f"   {i+1}. {fig3_names[i]} | -log10p={fig3_values[i]:.3f} | color={fig3_colors[i]}")
    except Exception as _e:
        print("[DEBUG] Fig3 print failed:", _e)

    fig3.add_trace(
        go.Bar(
            y=fig3_names,
            x=fig3_values,
            orientation='h',
            marker=dict(color=fig3_colors),
            showlegend=False
        )
    )
    
    fig3.update_layout(
        title={
            'text': "Pathway Significance (-log p-value) - Horizontal",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': int(16 * font_size_scale)}
        },
        xaxis=dict(
            title=dict(text="-log(p-value)", font=dict(weight='bold', size=int(14 * font_size_scale))),
            showline=True,
            linewidth=2,
            linecolor='black',
            tickfont=dict(weight='bold', size=int(12 * font_size_scale)),
            tickformat='.1f',
            range=[0, max(5, max(fig3_values, default=5) * 1.1)]
        ),
        yaxis=dict(
            title=dict(text="Pathways", font=dict(weight='bold', size=int(14 * font_size_scale))),
            showline=True,
            linewidth=2,
            linecolor='black',
            tickfont=dict(weight='bold', size=int(12 * font_size_scale)),
            categoryorder='array',
            categoryarray=fig3_names
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=horizontal_width,
        height=horizontal_height,
        margin=dict(l=300, r=50, t=80, b=80),
        bargap=0.3  # Add gap between bars
    )
    # Legend: only two entries (Activated/Red, Inhibited/Blue)
    fig3.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                              marker=dict(color='#FF4444', size=10),
                              name='Activated (Red)', showlegend=True))
    fig3.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                              marker=dict(color='#3A7CFF', size=10),
                              name='Inhibited (Blue)', showlegend=True))
    fig3.update_layout(showlegend=True,
                       legend=dict(orientation='v', yanchor='top', y=1, xanchor='left', x=1.02, font=dict(size=int(12 * font_size_scale))))
    
    
    # Figure 4: Bubble plot - Pathway impact (x: avg fold change, size: number of metabolites)
    fig4 = go.Figure()

    # Adjustable knobs so designers can tweak bubble prominence/spacing later without hunting
    bubble_size_scale = 1.5          # Increase for larger dots
    bubble_vertical_padding_pct = 0.2  # Increase for more vertical breathing room
    
    # Calculate bubble sizes based on levels instead of direct numbers
    # CRITICAL: Use 'hits' (significant only) not 'metabolites' (all measured)
    metabolite_counts = [len(p.get('hits', p.get('metabolites', []))) for p in pathway_list]
    if metabolite_counts:
        # Create size levels (quantiles) for better visualization
        unique_counts = sorted(set(metabolite_counts))
        
        def get_size_level(count):
            # Map counts to 5 levels with reasonable sizes
            if len(unique_counts) == 1:
                return 30  # Default size if all pathways have same metabolite count
            
            # Find position in sorted unique counts
            position = unique_counts.index(count)
            level = int((position / (len(unique_counts) - 1)) * 4) + 1  # 1-5 levels
            
            # Map levels to sizes: level 1 = 20, level 5 = 60
            return 20 + (level - 1) * 10
        
        # CRITICAL: Use 'hits' (significant only) not 'metabolites' (all measured)
        bubble_sizes = [get_size_level(len(p.get('hits', p.get('metabolites', [])))) for p in pathway_list]
    else:
        bubble_sizes = [30] * len(pathway_list)

    bubble_sizes = [size * bubble_size_scale for size in bubble_sizes]
    
    bubble_y_values = [-np.log10(p['combined_pvalue']) if p['combined_pvalue'] > 0 else 0 for p in pathway_list]

    fig4.add_trace(
        go.Scatter(
            x=[p['avg_log2_fold_change'] for p in pathway_list],
            y=bubble_y_values,
            mode='markers+text',
            marker=dict(
                size=bubble_sizes,
                color=pathway_colors,  # Use z-score based colors
                line=dict(width=1, color='white'),  # White border for better visibility
                opacity=0.8
            ),
            text=[str(p['bubble_number']) for p in pathway_list],
            textposition="middle center",
            textfont=dict(size=int(12 * font_size_scale), family='Arial Black', color='white'),
            name='Pathway Impact',
            hovertemplate='<b>Pathway %{text}</b><br>' +
                         'Avg Fold Change: %{x:.2f}<br>' +
                         'Num Metabolites: %{customdata}<br>' +
                         '-log(p-value): %{y:.2f}<br>' +
                         '<extra></extra>',
            # CRITICAL: Use 'hits' (significant only) not 'metabolites' (all measured)
            customdata=[len(p.get('hits', p.get('metabolites', []))) for p in pathway_list]
        )
    )
    
    # Create legend text that fits better
    legend_items = []
    for p in pathway_list:
        short_name = p['name'][:40] + '...' if len(p['name']) > 40 else p['name']
        legend_items.append(f"{p['bubble_number']}: {short_name}")
    
    # Split legend into columns if too many items
    legend_text = '<b>Bubble Legend:</b><br>'
    if len(legend_items) > 15:  # If more than 15 items, split into 2 columns
        mid = len(legend_items) // 2
        col1 = legend_items[:mid]
        col2 = legend_items[mid:]
        
        legend_text += '<br>'.join(col1)
        legend_text += f'<br><br><b>Continued:</b><br>'
        legend_text += '<br>'.join(col2)
    else:
        legend_text += '<br>'.join(legend_items)

    if bubble_y_values:
        y_min = min(bubble_y_values)
        y_max = max(bubble_y_values)
        y_padding = max(0.5, (y_max - y_min) * bubble_vertical_padding_pct)
    else:
        y_min, y_max, y_padding = 0, 1, 0.5

    fig4.update_layout(
        title={
            'text': "Metabolites Pathway Bubble Plot Analysis<br><sub>Bubble size represents metabolite count level</sub>",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': int(16 * font_size_scale)}
        },
        xaxis=dict(
            title=dict(text="Average log2(Fold Change)", font=dict(weight='bold', size=int(14 * font_size_scale))),
            zeroline=True,
            zerolinecolor='black',
            zerolinewidth=1,
            showline=True,
            linewidth=1,
            linecolor='black',
            tickfont=dict(weight='bold', size=int(12 * font_size_scale)),
            tickformat='.1f'
        ),
        yaxis=dict(
            title=dict(text="-log10(p-value)", font=dict(weight='bold', size=int(14 * font_size_scale))),
            showline=True,
            linewidth=1,
            linecolor='black',
            tickfont=dict(weight='bold', size=int(12 * font_size_scale)),
            tickformat='.1f',
            range=[max(0, y_min - y_padding), y_max + y_padding]
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=bubble_width,
        height=bubble_height,
        margin=dict(l=80, r=450, t=120, b=80),
        showlegend=False  # Remove default legend since we're using custom annotation
    )
    
    # Add custom legend as annotation outside plot area
    fig4.add_annotation(
        text=legend_text,
        x=1.05, y=1,
        xref="paper", yref="paper",
        align="left",
        showarrow=False,
        font=dict(size=int(11 * font_size_scale), family='Arial'),
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="rgba(0,0,0,0.2)",
        borderwidth=1,
        valign="top",
        xanchor="left",
        yanchor="top"
    )
    
    # Add vertical line at x=0
    fig4.add_vline(x=0, line_dash="dash", line_color="gray", opacity=0.5)
    
    # Figure 5: Pathway Enrichment Plot - Top 20 pathways by count and -log p-value
    # Layout knobs keep multi-line pathway labels readable (adjust here later if needed)
    enrichment_row_base_height = 36
    enrichment_extra_per_line = 18
    enrichment_bar_gap = 0.32
    enrichment_bar_width = 0.82

    def _count_wrapped_lines(label):
        """Count display lines in labels that may use either <br> or \n wrapping."""
        normalized = str(label).replace('<br>', '\n')
        return max(1, len([ln for ln in normalized.split('\n') if ln != '']))

    def _compute_enrichment_height(label_list):
        total = 260
        for label in label_list:
            lines = _count_wrapped_lines(label)
            total += enrichment_row_base_height + enrichment_extra_per_line * (lines - 1)
        return max(650, min(1700, total))

    # Sort pathways by -log(p-value) descending to get most enriched
    sorted_pathways = sorted(pathway_list, key=lambda p: p['neg_log_pvalue'], reverse=True)
    
    # Take top 20 (or fewer if less than 20 pathways exist)
    top_pathways = sorted_pathways[:min(20, len(sorted_pathways))]
    
    print(f"\n[ENRICHMENT PLOT - TOP 20]")
    print(f"  Total pathways in pathway_list: {len(pathway_list)}")
    print(f"  Requesting top 20, actually showing: {len(top_pathways)} pathways")
    if len(pathway_list) < 20:
        print(f"  NOTE: Less than 20 pathways available, showing all {len(pathway_list)}")
    
    # Get counts (number of SIGNIFICANT metabolites in each pathway)
    # CRITICAL: Use 'hits' (significant only) not 'metabolites' (all measured)
    counts = [len(p.get('hits', p.get('metabolites', []))) for p in top_pathways]
    names = [
        wrap_pathway_name(
            p['name'],
            max_chars=bargraph_pathway_max_chars,
            max_lines=bargraph_pathway_max_lines,
            break_token='<br>',
        )
        for p in top_pathways
    ]
    neg_log_pvals = [p['neg_log_pvalue'] for p in top_pathways]
    
    # Create color gradient from deep orange to gray based on -log p-values
    # Normalize -log p-values to 0-1 range for color mapping
    if neg_log_pvals:
        min_pval = min(neg_log_pvals)
        max_pval = max(neg_log_pvals)
        if max_pval > min_pval:
            normalized = [(val - min_pval) / (max_pval - min_pval) for val in neg_log_pvals]
        else:
            normalized = [1.0] * len(neg_log_pvals)
    else:
        normalized = []
    
    # Generate colors: deep orange (#FF6B00) for high significance, gray (#808080) for low
    def interpolate_color(value):
        """Interpolate between deep orange and gray based on value (0-1)"""
        # Deep orange RGB: (255, 107, 0)
        # Gray RGB: (128, 128, 128)
        r = int(255 - value * (255 - 128))
        g = int(107 + value * (128 - 107))
        b = int(0 + value * (128 - 0))
        return f'rgb({r},{g},{b})'
    
    enrichment_colors = [interpolate_color(1 - norm) for norm in normalized]  # Invert so highest = deepest orange
    
    # Reverse order so most significant is at the top
    names = names[::-1]
    counts = counts[::-1]
    enrichment_colors = enrichment_colors[::-1]
    neg_log_pvals = neg_log_pvals[::-1]
    
    # Calculate dynamic height and left margin based on longest pathway name
    enrichment_height = _compute_enrichment_height(names)
    
    # Calculate left margin dynamically based on longest display line.
    max_line_length = 20
    for name in names:
        lines = str(name).replace('<br>', '\n').split('\n')
        for line in lines:
            max_line_length = max(max_line_length, len(line))
    # Average character width in pixels at font size 11 is approximately 7 pixels
    # Add buffer for safety
    left_margin = min(800, max(350, max_line_length * 7 + 50))
    font_size_y = max(8, min(12, 11)) if max_line_length < 100 else max(7, 11 - (max_line_length - 100) // 50)
    
    fig5 = go.Figure()
    
    fig5.add_trace(
        go.Bar(
            y=names,
            x=counts,
            orientation='h',
            marker=dict(
                color=enrichment_colors,
                line=dict(width=1, color='white')
            ),
            showlegend=False,
            hovertemplate='<b>%{y}</b><br>' +
                         'Count: %{x}<br>' +
                         '-log10(p-value): %{customdata:.2f}<br>' +
                         '<extra></extra>',
            customdata=neg_log_pvals,
            width=enrichment_bar_width
        )
    )
    
    # Add a dummy scatter trace for the color bar
    max_neg_log_pval = max(neg_log_pvals) if neg_log_pvals else 1.0
    min_neg_log_pval = 0  # Always start from 0
    bargraph_legend_len = max(0.3, min(1.0, float(bargraph_legend_length_percent) / 100.0))
    
    fig5.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(
                colorscale=[
                    [0, 'rgb(128,128,128)'],      # Gray at 0
                    [1, 'rgb(255,107,0)']          # Deep orange at max
                ],
                cmin=min_neg_log_pval,
                cmax=max_neg_log_pval,
                colorbar=dict(
                    title=dict(
                        text="-log10(p-value)",
                        font=dict(size=int(bargraph_legend_title_font_size * font_size_scale), weight='bold')
                    ),
                    thickness=max(6, int(bargraph_legend_thickness)),
                    len=bargraph_legend_len,
                    x=1.02,
                    xanchor='left',
                    tickfont=dict(size=int(bargraph_legend_tick_font_size * font_size_scale), weight='bold'),
                    tickformat='.1f'
                ),
                showscale=True
            ),
            showlegend=False,
            hoverinfo='skip'
        )
    )
    
    fig5.update_layout(
        title={
            'text': "Pathway Enrichment Analysis (Top 20)",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': int(bargraph_title_font_size * font_size_scale), 'weight': 'bold'}
        },
        xaxis=dict(
            title=dict(text="Metabolite Count", font=dict(weight='bold', size=int(bargraph_axis_title_font_size * font_size_scale))),
            showline=True,
            linewidth=2,
            linecolor='black',
            tickfont=dict(weight='bold', size=int(bargraph_tick_font_size * font_size_scale)),
            tickformat='d',
            range=[0, max(counts, default=1) * 1.1] if counts else [0, 1]
        ),
        yaxis=dict(
            title=dict(text="Pathways", font=dict(weight='bold', size=int(bargraph_axis_title_font_size * font_size_scale))),
            showline=True,
            linewidth=2,
            linecolor='black',
            tickfont=dict(weight='bold', size=int(max(6, bargraph_pathway_label_font_size) * font_size_scale)),
            automargin=True,
            categoryorder='array',
            categoryarray=names
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=max(int(bargraph_width_px), left_margin + 600),
        height=max(int(bargraph_height_px), enrichment_height),
        margin=dict(l=left_margin, r=150, t=80, b=60),
        bargap=enrichment_bar_gap
    )

    fig5.update_traces(selector=dict(type='bar'), width=enrichment_bar_width)
    
    # Figure 6: Pathway Enrichment Plot - Selected Pathways (uses all pathways in pathway_list)
    # Sort pathways by -log(p-value) descending
    selected_sorted = sorted(pathway_list, key=lambda p: p['neg_log_pvalue'], reverse=True)
    
    # Get counts (number of SIGNIFICANT metabolites in each pathway)
    # CRITICAL: Use 'hits' (significant only) not 'metabolites' (all measured)
    selected_counts = [len(p.get('hits', p.get('metabolites', []))) for p in selected_sorted]
    # Keep selected pathways wrapped the same way as Top 20 for consistent splitting.
    selected_names = [
        wrap_pathway_name(
            p['name'],
            max_chars=bargraph_pathway_max_chars,
            max_lines=bargraph_pathway_max_lines,
            break_token='<br>',
        )
        for p in selected_sorted
    ]
    selected_neg_log_pvals = [p['neg_log_pvalue'] for p in selected_sorted]
    
    # Create color gradient from deep orange to gray based on -log p-values
    if selected_neg_log_pvals:
        sel_min_pval = min(selected_neg_log_pvals)
        sel_max_pval = max(selected_neg_log_pvals)
        if sel_max_pval > sel_min_pval:
            sel_normalized = [(val - sel_min_pval) / (sel_max_pval - sel_min_pval) for val in selected_neg_log_pvals]
        else:
            sel_normalized = [1.0] * len(selected_neg_log_pvals)
    else:
        sel_normalized = []
    
    selected_enrichment_colors = [interpolate_color(1 - norm) for norm in sel_normalized]
    
    # Reverse order so most significant is at the top
    selected_names = selected_names[::-1]
    selected_counts = selected_counts[::-1]
    selected_enrichment_colors = selected_enrichment_colors[::-1]
    selected_neg_log_pvals = selected_neg_log_pvals[::-1]
    
    # Calculate dynamic height and left margin based on longest pathway name
    selected_enrichment_height = _compute_enrichment_height(selected_names)
    
    # Calculate left margin dynamically based on longest display line
    # Handle both <br> and \n line break styles.
    selected_max_line_length = 20
    for name in selected_names:
        lines = str(name).replace('<br>', '\n').split('\n')
        for line in lines:
            selected_max_line_length = max(selected_max_line_length, len(line))
    # Average character width in pixels at font size 11 is approximately 7 pixels
    # Add buffer for safety
    selected_left_margin = min(800, max(350, selected_max_line_length * 7 + 50))
    selected_font_size_y = max(8, min(12, 11)) if selected_max_line_length < 100 else max(7, 11 - (selected_max_line_length - 100) // 50)
    
    fig6 = go.Figure()
    
    fig6.add_trace(
        go.Bar(
            y=selected_names,
            x=selected_counts,
            orientation='h',
            marker=dict(
                color=selected_enrichment_colors,
                line=dict(width=1, color='white')
            ),
            showlegend=False,
            hovertemplate='<b>%{y}</b><br>' +
                         'Count: %{x}<br>' +
                         '-log10(p-value): %{customdata:.2f}<br>' +
                         '<extra></extra>',
            customdata=selected_neg_log_pvals,
            width=enrichment_bar_width
        )
    )
    
    # Add a dummy scatter trace for the color bar
    sel_max_neg_log_pval = max(selected_neg_log_pvals) if selected_neg_log_pvals else 1.0
    sel_min_neg_log_pval = 0
    
    fig6.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(
                colorscale=[
                    [0, 'rgb(128,128,128)'],
                    [1, 'rgb(255,107,0)']
                ],
                cmin=sel_min_neg_log_pval,
                cmax=sel_max_neg_log_pval,
                colorbar=dict(
                    title=dict(
                        text="-log10(p-value)",
                        font=dict(size=int(bargraph_legend_title_font_size * font_size_scale), weight='bold')
                    ),
                    thickness=max(6, int(bargraph_legend_thickness)),
                    len=bargraph_legend_len,
                    x=1.02,
                    xanchor='left',
                    tickfont=dict(size=int(bargraph_legend_tick_font_size * font_size_scale), weight='bold'),
                    tickformat='.1f'
                ),
                showscale=True
            ),
            showlegend=False,
            hoverinfo='skip'
        )
    )
    
    fig6.update_layout(
        title={
            'text': "Pathway Enrichment Analysis (Selected Pathways)",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': int(bargraph_title_font_size * font_size_scale), 'weight': 'bold'}
        },
        xaxis=dict(
            title=dict(text="Metabolite Count", font=dict(weight='bold', size=int(bargraph_axis_title_font_size * font_size_scale))),
            showline=True,
            linewidth=2,
            linecolor='black',
            tickfont=dict(weight='bold', size=int(bargraph_tick_font_size * font_size_scale)),
            tickformat='d',
            range=[0, max(selected_counts, default=1) * 1.1] if selected_counts else [0, 1]
        ),
        yaxis=dict(
            title=dict(text="Pathways", font=dict(weight='bold', size=int(bargraph_axis_title_font_size * font_size_scale))),
            showline=True,
            linewidth=2,
            linecolor='black',
            tickfont=dict(weight='bold', size=int(max(6, bargraph_pathway_label_font_size) * font_size_scale)),
            automargin=True,
            categoryorder='array',
            categoryarray=selected_names
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=max(int(bargraph_width_px), selected_left_margin + 600),
        height=max(int(bargraph_height_px), selected_enrichment_height),
        margin=dict(l=selected_left_margin, r=150, t=80, b=60),
        bargap=enrichment_bar_gap
    )

    fig6.update_traces(selector=dict(type='bar'), width=enrichment_bar_width)
    
    # Figure 7: Chord Diagram - Pathway to Metabolite Connections
    # Use top pathways (same as fig5 for consistency) and their 'hits' metabolites
    fig_legend = None  # Initialize chord diagram legend figure
    print(f"\n[CHORD DIAGRAM - TOP PATHWAYS]")
    print(f"  Creating chord diagram for {len(top_pathways)} pathways")
    
    try:
        import pycirclize
        
        # Build chord data: collect all pathway->metabolite edges with metabolite log2FC
        chord_edges = []
        metabolite_labels = []
        metabolite_seen = set()
        metabolite_log2fc_map = {}  # Map metabolite name to its log2FC
        
        # Assign distinct colors to pathways
        num_pathways_chord = len(top_pathways)
        pathway_colors_chord = {}
        
        for pidx, pathway in enumerate(top_pathways):
            hue = (pidx / max(1, num_pathways_chord - 1)) if num_pathways_chord > 1 else 0.5
            saturation = 0.7
            lightness = 0.5
            rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
            hex_color = '#{:02x}{:02x}{:02x}'.format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
            pathway_colors_chord[pathway['name']] = hex_color
            
            # Collect edges for this pathway (using 'hits' for significant metabolites)
            for met_name in pathway['hits']:
                if met_name in metabolites:
                    met_data = metabolites[met_name]
                    log2fc = met_data.get('log2FC', 0)
                    
                    chord_edges.append({
                        'pathway': pathway['name'],
                        'metabolite': met_name,
                        'log2fc': log2fc
                        # NOTE: Chord color will be pathway color, not metabolite color
                    })
                    
                    # Track unique metabolites and their log2FC
                    if met_name not in metabolite_seen:
                        metabolite_labels.append(met_name)
                        metabolite_seen.add(met_name)
                        metabolite_log2fc_map[met_name] = log2fc
                    else:
                        # Update with latest log2FC value (in case metabolite appears in multiple pathways)
                        metabolite_log2fc_map[met_name] = log2fc
        
        print(f"  Total chord edges: {len(chord_edges)}")
        print(f"  Pathway colors assigned: {len(pathway_colors_chord)}")
        print(f"  Unique metabolites: {len(metabolite_seen)}")
        
        if chord_edges and len(metabolite_labels) > 0:
            # Calculate actual log2FC range from metabolites for proper color scaling
            if metabolite_log2fc_map:
                log2fc_values = list(metabolite_log2fc_map.values())
                min_log2fc = min(log2fc_values)
                max_log2fc = max(log2fc_values)
                
                # Ensure symmetric range for better color representation
                max_abs_log2fc = max(abs(min_log2fc), abs(max_log2fc))
                # Use at least 0.5 to avoid too narrow range
                max_abs_log2fc = max(0.5, max_abs_log2fc)
                
                log2fc_range_min = -max_abs_log2fc
                log2fc_range_max = max_abs_log2fc
                
                print(f"  Log2FC range: {min_log2fc:.3f} to {max_log2fc:.3f}")
                print(f"  Using symmetric range for colors: {log2fc_range_min:.3f} to {log2fc_range_max:.3f}")
            else:
                # Fallback defaults
                log2fc_range_min = -2.0
                log2fc_range_max = 2.0
            
            # Create sectors: pathways and metabolites
            sectors = {}
            for pathway in top_pathways:
                sectors[pathway['name']] = len([e for e in chord_edges if e['pathway'] == pathway['name']])
            
            for met in metabolite_labels:
                sectors[met] = len([e for e in chord_edges if e['metabolite'] == met])
            
            # Initialize circos plot with sectors
            # pycirclize treats `space` as the gap between every sector.
            # With many sectors, a fixed gap can exceed the 360-degree circle and fail.
            num_sectors = max(1, len(sectors))
            safe_sector_space = min(5, max(0, 360 // num_sectors))
            circos = pycirclize.Circos(sectors, space=safe_sector_space)
            
            # Add colored arc tracks and labels for all sectors
            for sector in circos.sectors:
                # First: Add colored arc track for this sector
                try:
                    # Determine the appropriate color for this sector
                    if sector.name in metabolite_log2fc_map:
                        # Metabolite: color based on log2FC gradient
                        log2fc = metabolite_log2fc_map[sector.name]
                        
                        # Map log2FC to color gradient using actual data range
                        # Normalize to 0-1 range using calculated min/max
                        norm_value = (log2fc - log2fc_range_min) / (log2fc_range_max - log2fc_range_min)
                        norm_value = max(0, min(1, norm_value))  # Clamp to 0-1
                        
                        # Color gradient from green (#22BB22) through gray (#CCCCCC) to red (#FF4444)
                        if norm_value < 0.5:
                            # Green to Gray: interpolate 0->0.5
                            t = norm_value * 2
                            r = int(0x22 + (t * 0xAA))  # 34 to 204
                            g = int(0xBB + (t * 0x11))  # 187 to 204
                            b = int(0x22 + (t * 0xAA))  # 34 to 204
                        else:
                            # Gray to Red: interpolate 0.5->1
                            t = (norm_value - 0.5) * 2
                            r = int(0xCC + (t * 0x33))  # 204 to 255
                            g = int(0xCC - (t * 0xAA))  # 204 to 34
                            b = int(0xCC - (t * 0xAA))  # 204 to 34
                        
                        sector_color = f'#{r:02x}{g:02x}{b:02x}'
                    
                    elif sector.name in pathway_colors_chord:
                        # Pathway: use same color as its chords
                        sector_color = pathway_colors_chord[sector.name]
                    
                    else:
                        # Default color if not found
                        sector_color = '#CCCCCC'
                    
                    # Add colored arc track at the perimeter (radius 95-100)
                    track = sector.add_track((95, 100))
                    track.rect(sector.start, sector.end, fc=sector_color, ec='black', lw=0.5)
                    
                except Exception as track_err:
                    print(f"  Warning: Could not add color track for {sector.name}: {track_err}")
                
                # Second: Add text label with proper rotation (matching chord_diagram.py exactly)
                label_text = sector.name
                
                # Apply formatting: optional sentence case for metabolites, uppercase for pathways
                if sector.name in metabolite_log2fc_map:
                    # Metabolite: preserve source case when requested, otherwise sentence case.
                    if not chord_keep_metabolite_case:
                        label_text = label_text[0].upper() + label_text[1:].lower() if len(label_text) > 0 else label_text
                else:
                    # Pathway: keep in all caps
                    label_text = label_text.upper()

                label_text = wrap_pathway_name(
                    label_text,
                    max_chars=chord_metabolite_max_chars if sector.name in metabolite_log2fc_map else chord_pathway_max_chars,
                    max_lines=chord_metabolite_max_lines if sector.name in metabolite_log2fc_map else chord_pathway_max_lines,
                    break_token="\n",
                )
                
                # Calculate proper text rotation and alignment based on sector position
                try:
                    # Get sector degree position
                    deg_start, deg_end = sector.deg_lim
                    mid_deg = (deg_start + deg_end) / 2.0
                    
                    # Compute rotation/ha from side only: all RIGHT sectors share the same
                    # baseline direction, and all LEFT sectors share the same baseline.
                    # This keeps orientation consistent and does not rely on names.
                    if 90 < mid_deg < 270:
                        # LEFT side: radial outward + 180 so baseline reads L→R consistently
                        rotation = 90 - mid_deg + 180
                        ha = 'right'
                    else:
                        # RIGHT side: radial outward, keep left alignment
                        rotation = 90 - mid_deg
                        ha = 'left'

                    # normalize rotation to [-180, 180] then to readable range [-90, 90]
                    rotation = ((rotation + 180) % 360) - 180
                    if rotation < -90:
                        rotation += 180
                        ha = 'right' if ha == 'left' else 'left'
                    elif rotation > 90:
                        rotation -= 180
                        ha = 'right' if ha == 'left' else 'left'

                    # Use consistent label radius for all positions
                    r_used = 112
                    label_font_size = (
                        chord_metabolite_font_size if sector.name in metabolite_log2fc_map else chord_pathway_font_size
                    )
                    
                    # Add label at sector position with proper rotation
                    sector.text(
                        label_text,
                        r=r_used,
                        size=int(max(6, label_font_size) * font_size_scale),
                        weight='bold',
                        rotation=rotation,
                        ha=ha,
                        va='center',
                        rotation_mode='anchor',
                        adjust_rotation=False  # respect provided rotation/ha instead of auto-adjust
                    )
                except Exception as label_err:
                    print(f"  Warning: Could not add label for {sector.name}: {label_err}")
            
            # Draw chords for each edge using PATHWAY colors
            for edge in chord_edges:
                try:
                    pathway_sector = circos.get_sector(edge['pathway'])
                    met_sector = circos.get_sector(edge['metabolite'])
                    
                    # Draw chord with PATHWAY color (not metabolite color)
                    pathway_color = pathway_colors_chord[edge['pathway']]
                    linewidth = max(0.5, min(2.0, abs(edge['log2fc']) * 0.4 + 0.5))
                    circos.link(
                        (pathway_sector.name, pathway_sector.start, pathway_sector.end),
                        (met_sector.name, met_sector.start, met_sector.end),
                        color=pathway_color,
                        alpha=0.5,
                        linewidth=linewidth
                    )
                except Exception as link_err:
                    print(f"  Warning: Could not draw link {edge['pathway']} -> {edge['metabolite']}: {link_err}")
            
            # Render the circos diagram
            fig7 = circos.plotfig()
            
            # If plotfig() doesn't return a figure, get current figure
            if fig7 is None:
                fig7 = plt.gcf()
            
            # Set background color
            fig7.patch.set_facecolor('white')
            
            # Set dynamic size based on number of sectors for better label accommodation
            # More sectors = need more space for labels
            num_sectors = len(sectors)
            # Base size 8 inches, add 0.1 inch per sector (up to max 12 inches)
            dynamic_size = min(12, max(8, 8 + (num_sectors - 20) * 0.1))
            target_w = max(float(chord_width_inches), dynamic_size) if chord_width_inches else dynamic_size
            target_h = max(float(chord_height_inches), dynamic_size) if chord_height_inches else dynamic_size
            fig7.set_size_inches(target_w, target_h)
            
            print("[CHORD DIAGRAM] Successfully created\n")
            
            # Generate separate legend figure for the chord diagram
            try:
                fig_legend = create_chord_diagram_legend(pathway_colors_chord, log2fc_range_min, log2fc_range_max)
                print("[CHORD DIAGRAM LEGEND] Successfully created\n")
            except Exception as legend_err:
                print(f"[⚠️ CHORD DIAGRAM LEGEND] Error creating legend: {legend_err}")
                fig_legend = None
            
        else:
            # If no edges, create a placeholder figure
            fig7 = go.Figure()
            fig7.add_annotation(
                text="No edges found for chord diagram",
                xref="paper", yref="paper",
                x=0.5, y=0.5,
                showarrow=False
            )
            fig7.update_layout(
                title="Pathway-Metabolite Chord Diagram",
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                plot_bgcolor='white',
                paper_bgcolor='white',
                width=int(max(600, float(chord_width_inches) * 100)),
                height=int(max(450, float(chord_height_inches) * 100))
            )
            print("[CHORD DIAGRAM] No edges - placeholder created\n")
            print("[CHORD DIAGRAM] No edges - placeholder created\n")
    
    except ImportError:
        # pycirclize not installed - create placeholder with warning
        print("[⚠️ CHORD DIAGRAM] pycirclize not installed. Install with: pip install pycirclize")
        fig7 = go.Figure()
        fig7.add_annotation(
            text="⚠️ pycirclize library not installed\nInstall with: pip install pycirclize",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=14, color='red')
        )
        fig7.update_layout(
            title="Pathway-Metabolite Chord Diagram (Requires pycirclize)",
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            plot_bgcolor='white',
            paper_bgcolor='white',
            width=int(max(600, float(chord_width_inches) * 100)),
            height=int(max(450, float(chord_height_inches) * 100))
        )
    except Exception as e:
        # Catch any other errors and show warning
        print(f"[⚠️ CHORD DIAGRAM] Error creating chord diagram: {e}")
        fig7 = go.Figure()
        fig7.add_annotation(
            text=f"❌ Error creating chord diagram:\n{str(e)}",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=12, color='red')
        )
        fig7.update_layout(
            title="Pathway-Metabolite Chord Diagram (Error)",
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            plot_bgcolor='white',
            paper_bgcolor='white',
            width=int(max(600, float(chord_width_inches) * 100)),
            height=int(max(450, float(chord_height_inches) * 100))
        )
    
    return fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig_legend

def create_chord_diagram_legend(pathway_colors_dict, log2fc_min, log2fc_max):
    """Create a standalone legend figure for the chord diagram.
    
    Args:
        pathway_colors_dict: Dictionary mapping pathway names to their hex colors
        log2fc_min: Minimum log2FC value for metabolite color scale
        log2fc_max: Maximum log2FC value for metabolite color scale
    
    Returns:
        matplotlib Figure object with the legend
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.colors as mcolors
    
    # Create figure
    fig, (ax_pathways, ax_metabolites) = plt.subplots(1, 2, figsize=(10, 6), 
                                                        gridspec_kw={'width_ratios': [2, 1]})
    
    # --- Left panel: Pathway Colors (Chord Colors) ---
    ax_pathways.axis('off')
    ax_pathways.set_xlim(0, 1)
    ax_pathways.set_ylim(0, 1)
    
    # Title
    ax_pathways.text(0.5, 0.95, 'Chord Colors = Pathways', 
                     ha='center', va='top', fontsize=14, fontweight='bold')
    
    # Pathway legend entries
    num_pathways = len(pathway_colors_dict)
    if num_pathways > 0:
        y_start = 0.88
        y_step = min(0.06, 0.75 / max(1, num_pathways))  # Adjust spacing based on number
        
        for idx, (pathway_name, color) in enumerate(sorted(pathway_colors_dict.items())):
            y_pos = y_start - idx * y_step
            
            # Color box
            rect = mpatches.Rectangle((0.05, y_pos - 0.015), 0.08, 0.03, 
                                       facecolor=color, edgecolor='black', linewidth=1)
            ax_pathways.add_patch(rect)
            
            # Pathway name (truncate if too long)
            display_name = pathway_name[:45] + '...' if len(pathway_name) > 45 else pathway_name
            ax_pathways.text(0.16, y_pos, display_name, 
                            va='center', ha='left', fontsize=10, fontweight='bold')
    
    # --- Right panel: Metabolite Log2FC Color Scale ---
    ax_metabolites.axis('off')
    ax_metabolites.set_xlim(0, 1)
    ax_metabolites.set_ylim(0, 1)
    
    # Title
    ax_metabolites.text(0.5, 0.95, 'Metabolite Colors', 
                        ha='center', va='top', fontsize=14, fontweight='bold')
    
    # Create color gradient bar
    gradient_height = 0.5
    gradient_bottom = 0.25
    n_colors = 100
    
    for i in range(n_colors):
        norm_value = i / (n_colors - 1)
        y_pos = gradient_bottom + (i / n_colors) * gradient_height
        
        # Calculate color (matching chord diagram gradient: green -> gray -> red)
        if norm_value < 0.5:
            t = norm_value * 2
            r = int(0x22 + (t * 0xAA))
            g = int(0xBB + (t * 0x11))
            b = int(0x22 + (t * 0xAA))
        else:
            t = (norm_value - 0.5) * 2
            r = int(0xCC + (t * 0x33))
            g = int(0xCC - (t * 0xAA))
            b = int(0xCC - (t * 0xAA))
        
        color = f'#{r:02x}{g:02x}{b:02x}'
        rect = mpatches.Rectangle((0.3, y_pos), 0.4, gradient_height / n_colors,
                                   facecolor=color, edgecolor='none')
        ax_metabolites.add_patch(rect)
    
    # Border around gradient
    border = mpatches.Rectangle((0.3, gradient_bottom), 0.4, gradient_height,
                                 facecolor='none', edgecolor='black', linewidth=2)
    ax_metabolites.add_patch(border)
    
    # Labels
    ax_metabolites.text(0.15, gradient_bottom + gradient_height, f'{log2fc_max:.2f}',
                        va='center', ha='right', fontsize=11, fontweight='bold')
    ax_metabolites.text(0.15, gradient_bottom + gradient_height/2, '0.00',
                        va='center', ha='right', fontsize=11, fontweight='bold')
    ax_metabolites.text(0.15, gradient_bottom, f'{log2fc_min:.2f}',
                        va='center', ha='right', fontsize=11, fontweight='bold')
    
    # Label text
    ax_metabolites.text(0.5, 0.10, 'Log2 Fold Change',
                        ha='center', va='top', fontsize=11, fontweight='bold')
    ax_metabolites.text(0.5, 0.83, '(Upregulated)',
                        ha='center', va='top', fontsize=9, color='#E53935')
    ax_metabolites.text(0.5, gradient_bottom - 0.02, '(Downregulated)',
                        ha='center', va='top', fontsize=9, color='#18A558')
    
    fig.patch.set_facecolor('white')
    fig.tight_layout()
    
    return fig

def _safe_div(numerator, denominator, default=0.0):
    """Safely divide two numbers; returns default if denominator is zero or NaN."""
    try:
        if denominator is None:
            return default
        if isinstance(denominator, (int, float)) and denominator == 0:
            return default
        val = numerator / denominator
        if isinstance(val, float) and (np.isnan(val) or np.isinf(val)):
            return default
        return val
    except Exception:
        return default

def compute_pathway_enrichment(metabolites, pathway_stats, df, pvalue_threshold=0.05, 
                               use_hypergeometric=True, fdr_method='fdr_bh'):
    """Compute pathway enrichment with proper statistical testing.

    ENRICHMENT CALCULATION:
    -----------------------
    Enrichment ratio = (k/K) / (M/N)
      - N: total measured metabolites (unique in df['Name'])
      - K: total significant metabolites in dataset (pvalue <= pvalue_threshold)
      - M: pathway size (metabolites present in dataset for this pathway)
      - k: significant metabolites within this pathway

    STATISTICAL TESTING:
    -------------------
    Uses hypergeometric test (standard in pathway enrichment analysis):
      - Null hypothesis: metabolites are randomly distributed across pathways
      - Alternative: pathway is enriched for significant metabolites
      - P-value calculated from hypergeometric distribution
      - FDR correction applied across all pathways (Benjamini-Hochberg by default)

    This approach aligns with standard methods used by:
      - DAVID, Metaboanalyst, Enrichr, KEGG pathway analysis
      - Gene Ontology (GO) enrichment tools
      - Most published metabolomics pathway analysis studies

    Parameters:
    -----------
    metabolites : dict
        Dictionary of metabolite data
    pathway_stats : dict
        Dictionary of pathway statistics
    df : pandas.DataFrame
        Metabolomics dataset with 'Name' and 'pvalue' columns
    pvalue_threshold : float, default=0.05
        Threshold for considering metabolites significant
    use_hypergeometric : bool, default=True
        Whether to calculate hypergeometric p-values (recommended)
    fdr_method : str, default='fdr_bh'
        FDR correction method: 'fdr_bh' (Benjamini-Hochberg), 'bonferroni', or None

    Returns:
    --------
    list of dict
        Each dict contains: Pathway_Name, k, M, K, N, Enrichment_Ratio, 
        Hypergeometric_Pvalue (if enabled), FDR_Qvalue (if correction applied), Status
    """
    from scipy.stats import hypergeom
    from statsmodels.stats.multitest import multipletests
    
    # Validate inputs
    df = df.copy()
    
    # Find the name column (support multiple column name variations)
    name_candidates = ['Name', 'LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite', 'name', 'Metabolite_Name']
    name_col = next((c for c in name_candidates if c in df.columns), None)
    if name_col is None:
        # Can't compute without names
        return []
    
    # Ensure pvalue column exists; if not, treat as all non-significant
    if 'pvalue' not in df.columns:
        df['pvalue'] = np.nan

    # Unique measured metabolites
    measured = df[name_col].dropna().astype(str).unique().tolist()
    measured_set = set(measured)
    N = len(measured_set)

    # Significant metabolites in dataset
    try:
        sig_series = df['pvalue'] <= float(pvalue_threshold)
    except Exception:
        sig_series = pd.Series([False] * len(df))
    sig_mets = set(df.loc[sig_series.fillna(False), name_col].dropna().astype(str).tolist())
    K = len(sig_mets)

    enrichment_rows = []
    for pname, stats in pathway_stats.items():
        mets_in_pathway = [m for m in stats.get('metabolites', []) if m in measured_set]
        M = len(mets_in_pathway)
        if M == 0:
            ratio = 0.0
            k = 0
            hypergeo_pval = 1.0
        else:
            k = sum(1 for m in mets_in_pathway if m in sig_mets)
            # Enrichment ratio = (k/K) / (M/N)
            ratio = _safe_div(_safe_div(k, K, default=0.0), _safe_div(M, N, default=np.nan), default=0.0)
            
            # Calculate hypergeometric p-value (one-tailed test for enrichment)
            # P(X >= k) where X follows hypergeometric distribution
            # This is the standard method used in pathway enrichment analysis
            if use_hypergeometric and K > 0:
                try:
                    # hypergeom.sf(k-1, N, M, K) = P(X >= k)
                    # N = population size (all measured metabolites)
                    # M = number of success states in population (pathway size)
                    # K = number of draws (significant metabolites)
                    # k = number of observed successes (sig mets in pathway)
                    hypergeo_pval = hypergeom.sf(k-1, N, M, K)
                except Exception as e:
                    print(f"Warning: Could not calculate hypergeometric p-value for {pname}: {e}")
                    hypergeo_pval = 1.0
            else:
                hypergeo_pval = 1.0

        row_data = {
            'Pathway_Name': pname,
            'k_hits_in_pathway': k,
            'M_pathway_size': M,
            'K_total_significant': K,
            'N_total_measured': N,
            'Enrichment_Ratio': ratio,
            'Status': stats.get('status', ''),
        }
        
        if use_hypergeometric:
            row_data['Hypergeometric_Pvalue'] = hypergeo_pval
            
        enrichment_rows.append(row_data)

    # Apply FDR correction if requested
    if use_hypergeometric and fdr_method and len(enrichment_rows) > 0:
        pvalues = [r['Hypergeometric_Pvalue'] for r in enrichment_rows]
        try:
            if fdr_method in ['fdr_bh', 'bonferroni']:
                reject, qvalues, _, _ = multipletests(pvalues, method=fdr_method)
                for i, row in enumerate(enrichment_rows):
                    row['FDR_Qvalue'] = qvalues[i]
                    row['Significant_FDR'] = reject[i]
            else:
                print(f"Warning: Unknown FDR method '{fdr_method}', no correction applied")
        except Exception as e:
            print(f"Warning: Could not apply FDR correction: {e}")

    # Sort by enrichment ratio descending, then by p-value if available
    if use_hypergeometric:
        enrichment_rows.sort(key=lambda r: (-r['Enrichment_Ratio'], r.get('Hypergeometric_Pvalue', 1.0)))
    else:
        enrichment_rows.sort(key=lambda r: (r['Enrichment_Ratio'], r['k_hits_in_pathway']), reverse=True)
    
    return enrichment_rows

def calculate_pathway_significance(df, pos_set_zscore=1, neg_set_zscore=-1, 
                                  max_pvalue=0.05, min_metabolites=2,
                                  max_total_pathways=None, max_pos=None, max_neg=None):
    """Calculate pathway significance statistics (wrapper for calculate_pathway_statistics)
    
    This function provides the same functionality as calculate_pathway_statistics
    but with a different name for backward compatibility.
    """
    return calculate_pathway_statistics(df, pos_set_zscore=pos_set_zscore, 
                                      neg_set_zscore=neg_set_zscore,
                                      max_pvalue=max_pvalue, 
                                      min_metabolites=min_metabolites,
                                      max_total_pathways=max_total_pathways,
                                      max_pos=max_pos, max_neg=max_neg)

def run_visualizations(excel_file="metabolite_pathways_annotated.xlsx",
                      sheet_name="Metabolites_with_Pathways",
                      pvalue_threshold=0.05,
                      pos_set_fc=0,
                      neg_set_fc=0,
                      pos_set_zscore=0,
                      neg_set_zscore=0,
                      min_metabolites=1,
                      max_pvalue=0.05,
                      max_total_pathways=20,
                      max_pos=15,
                      max_neg=5,
                      node_spacing=1.0,
                      layout_iterations=100,
                      network_width=1800,
                      network_height=1000,
                      metabolite_font_size=16,
                      pathway_font_size=14,
                      enzyme_font_size=11,
                      metabolite_max_chars=30,
                      metabolite_max_lines=1,
                      pathway_max_chars=30,
                      pathway_max_lines=3,
                      exclude_pathway_keywords=None):
    """Main function to create all visualizations
    
    Parameters:
    -----------
    excel_file : str
        Path to the Excel file containing metabolite pathway data
    sheet_name : str
        Name of the sheet to read from the Excel file
    pvalue_threshold : float
        P-value threshold for data filtering
    pos_set_fc : float
        Positive fold change threshold
    neg_set_fc : float
        Negative fold change threshold
    pos_set_zscore : float
        Positive z-score threshold for pathway grouping
    neg_set_zscore : float
        Negative z-score threshold for pathway grouping
    min_metabolites : int
        Minimum number of metabolites required per pathway
    max_pvalue : float
        Maximum p-value for pathway significance
    max_total_pathways : int
        Maximum total pathways to include
    max_pos : int
        Maximum positive pathways to include
    max_neg : int
        Maximum negative pathways to include
    node_spacing : float
        Spacing between nodes in network visualization
    layout_iterations : int
        Number of iterations for network layout algorithm
    network_width : int
        Width of network plot in pixels
    network_height : int
        Height of network plot in pixels
    metabolite_font_size : int
        Font size for metabolite labels
    pathway_font_size : int
        Font size for pathway labels
    enzyme_font_size : int
        Font size for enzyme/transporter labels
    metabolite_max_chars : int
        Maximum characters per line for metabolite text wrapping
    metabolite_max_lines : int
        Maximum lines for metabolite text wrapping
    pathway_max_chars : int
        Maximum characters per line for pathway text wrapping
    pathway_max_lines : int
        Maximum lines for pathway text wrapping
    exclude_pathway_keywords : list
        List of keywords to exclude pathways containing these terms
    """
        # Pathway filtering - Add keywords to exclude pathways containing these terms
    # All filtering is case-insensitive (will match upper, lower, or mixed case)
    exclude_pathway_keywords = [
        # Plant-related pathways
        "plant",
        "chlorophyll", 
        "photosynthesis",
        "cellulose",
        "lignin",
        "xylem",
        "phloem",
        "stomata",
        "thylakoid",
        
        # Microbial pathways (bacteria, archaea, fungi)
        "bacterial",
        "bacterium",
        "prokaryotic",
        "archaeal",
        "fungal",
        "yeast",
        "escherichia",
        "salmonella",
        "streptococcus",
        "staphylococcus",
        "mycobacterium",
        "bacillus",
        "clostridium",
        "pseudomonas",
        "candida",
        "aspergillus",
        "saccharomyces",
        
        # Marine/aquatic organisms
        "algae",
        "algal",
        "cyanobacteria",
        "marine",
        "seaweed",
        "diatom",
        
        # Parasites and pathogens
        "parasitic",
        "parasite",
        "protozoa",
        "plasmodium",
        "trypanosoma",
        "leishmania",
        "toxoplasma",
        "giardia",
        "entamoeba",
        
        # Viral pathways
        "viral",
        "virus",
        "bacteriophage",
        "retroviral",
        
        # Invertebrate-specific
        "insect",
        "drosophila",
        "caenorhabditis",
        "c. elegans",
        "arthropod",
        "nematode",
        "mollusc",
        
        # General exclusions for complex/superpathways
        "degradation",
        "antibiotic",       # Antibiotic pathways (neomycin, kanamycin, etc.)
        "secondary metabolite",  # Plant or microbial secondary metabolism
        "polyketide",       # Polyketide biosynthesis (microbial)
        "nonribosomal",     # NRPS peptide biosynthesis (microbial)
        "mycolyl",          # Mycolic acid (bacterial cell wall)
        "peptidoglycan",    # Bacterial/fungal cell wall
        "lipopolysaccharide", # Gram-negative bacterial outer membrane
        "chitin",           # Invertebrate/fungal chitin
        "terpenoid",        # Plant/microbial terpenoids
        "flavonoid",        # Plant metabolites
        "alkaloid",         # Plant alkaloids
        "glycoside",        # Plant glycosides
        "antimicrobial",    # General microbial metabolite pathways
        
        # Antibiotic and drug resistance (often microbial)
        "antibiotic",
        "resistance",
        "antimicrobial",
        
        # Environmental/industrial processes
        "xenobiotic",
        "detoxification",
        "biodegradation",
        "bioremediation",
        
        # Specific metabolic processes not relevant to mammalian studies
        "methanogenesis",
        "nitrogen fixation",
        "nitrification",
        "denitrification",
        "sulfur oxidation",
        "chemolithotrophic",
        
        # Add more keywords as needed for your specific study
    ]
    
    
    # Set default exclude keywords if none provided
    if exclude_pathway_keywords is None:
        exclude_pathway_keywords = []

    try:
        df = import_data(excel_file, 
                         sheet_name1=sheet_name,
                         pvalue_threshold=pvalue_threshold,
                         pos_set_fc=pos_set_fc,
                         neg_set_fc=neg_set_fc,
                         exclude_pathway_keywords=exclude_pathway_keywords
                         )
    except Exception as e:
        print(f"Error reading files: {e}")
        print("Please make sure the Excel files are not open in another program.")
        return

    print(f"  - Max p-value: {max_pvalue}")
    print(f"  - Pos z-score threshold for grouping: {pos_set_zscore}")
    print(f"  - Neg z-score threshold for grouping: {neg_set_zscore}")
  
    
    # Get pathway statistics first (no top pathway filtering here)
    metabolites = df.set_index('Name')[['log2FC', 'pvalue']].to_dict(orient='index')
    pathway_stats = calculate_pathway_statistics(df, 
                                               min_metabolites=min_metabolites,
                                               max_pvalue=max_pvalue,
                                               pos_set_zscore=pos_set_zscore,
                                               neg_set_zscore=neg_set_zscore,
                                               max_total_pathways=max_total_pathways,
                                               max_pos=max_pos, max_neg=max_neg)

    # Debug pathway statistics
    print(f"\nDEBUG: Final pathway statistics:")
    activated_count = 0
    inhibited_count = 0
    no_change_count = 0
    
    for pathway, stats in pathway_stats.items():
        print(f"Pathway: {pathway[:50]}...")
        print(f"  Status: {stats['status']}")
        print(f"  Z-score: {stats['z_score']:.3f}")
        print(f"  P-value: {stats['combined_pvalue']:.6f}")
        print(f"  Metabolites: {stats['n_metabolites']}")
        print(f"  Metabolite names: {stats['metabolites'][:3]}...")  # Show first 3 metabolites
        print("---")
        
        if stats['status'] == 'Activated':
            activated_count += 1
        elif stats['status'] == 'Inhibited':
            inhibited_count += 1
        else:
            no_change_count += 1
    
    print(f"\nPathway Status Summary:")
    print(f"- Activated: {activated_count}")
    print(f"- Inhibited: {inhibited_count}")
    print(f"- No Change: {no_change_count}")
    
    # Debug: Check unique metabolites across all pathways
    all_metabolites_in_pathways = set()
    for pathway, stats in pathway_stats.items():
        all_metabolites_in_pathways.update(stats['metabolites'])
    
    print(f"\nDEBUG: Total unique metabolites across all {len(pathway_stats)} pathways: {len(all_metabolites_in_pathways)}")
    print(f"First 10 unique metabolites: {list(all_metabolites_in_pathways)[:10]}")
    
    # Always keep pathways regardless of z-score band (include 'No Change')
    print("\nKeeping pathways regardless of z-score band (including 'No Change').")


        # Create individual summary plots for this group
    fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig_legend = create_pathway_summary_plots(metabolites, pathway_stats)

    # # Enrichment analysis: horizontal bar chart with enrichment ratios
    # fig5, enrichment_df = create_enrichment_bar_chart(
    #     metabolites, pathway_stats, df, pvalue_threshold=pvalue_threshold, top_n=30
    # )

    # Create the interactive network plot for this group using already filtered pathway_stats
    network_fig = create_interactive_metabolite_pathway_plot_from_stats(
        df,
        pathway_stats,  # Use the already filtered pathway stats
        pos_set_fc=pos_set_fc,
        neg_set_fc=neg_set_fc,
        node_spacing=node_spacing,
        layout_iterations=layout_iterations,
        metabolite_font_size=metabolite_font_size,
        pathway_font_size=pathway_font_size,
        enzyme_font_size=enzyme_font_size,
        metabolite_max_chars=metabolite_max_chars,
        metabolite_max_lines=metabolite_max_lines,
        pathway_max_chars=pathway_max_chars,
        pathway_max_lines=pathway_max_lines,
        network_width=network_width,
        network_height=network_height,
        # IMPORTANT: include upstream (enzymes/transporters) and diseases when requested
        # create_interactive_metabolite_pathway_plot_from_stats expects include_upstream/include_diseases
        include_upstream=True,
        include_diseases=True
    )

        # Update the network plot title to include group information
    network_fig.update_layout(
            title={
                'text': f"Metabolite-Pathway Network - <br><sub>Circles: Metabolites | Diamonds: Pathways | Dashed Lines: Interactions</sub>",
                'x': 0.5,
                'xanchor': 'center',
                'font': {'size': 16}
            },
            width=1600,
            height=1000,
            margin=dict(l=150, r=400, t=150, b=150),
            autosize=True,
            xaxis=dict(constrain='domain'),
            yaxis=dict(constrain='domain')
        )
        
        # Save plots with group identifier

        
    # Save each plot individually with error handling - SVG format only
    plots_to_save = [
        (fig1, "pathway_states_horizontal.svg", "Pathway States Horizontal"),
        (fig2, "pathway_significance_vertical.svg", "Pathway Significance Vertical"),
        (fig3, "pathway_significance_horizontal.svg", "Pathway Significance Horizontal"),
       #(fig4, "pathway_impact_bubble.svg", "Pathway Impact Bubble"),
        #(fig5, "pathway_enrichment_ratio.svg", "Pathway Enrichment Ratio"),
        (network_fig, "pathway_network_interactive.svg", "Pathway Network Interactive")
    ]
    
    print("📊 Saving plots as SVG files...")
    saved_count = 0
    for fig, filename, description in plots_to_save:
        try:
            if "bubble" in filename or "network" in filename:
                fig.write_image(filename, format="svg", width=1200, height=600)
            else:
                fig.write_image(filename, format="svg")
            print(f"✅ Saved {description}: {filename}")
            saved_count += 1
        except Exception as e:
            # Common cause: plotly image export requires 'kaleido'. Fallback to HTML if not available.
            print(f"⚠️ Failed to save {description} as SVG: {e}")
            try:
                html_name = filename.replace('.svg', '.html')
                fig.write_html(html_name)
                print(f"➡️  Fallback saved {description} as HTML: {html_name}")
            except Exception as e_html:
                print(f"❌ Also failed to save {description} as HTML: {e_html}")
    
    # Save HTML version only for network pathway (interactive features)
    try:
        network_fig.write_html("pathway_network_interactive.html")
        print("✅ Saved interactive network HTML: pathway_network_interactive.html")
    except Exception as e:
        print(f"⚠️ Failed to save network HTML: {e}")
    
    print(f"📊 Successfully saved {saved_count} out of {len(plots_to_save)} SVG plots")
    
    # Create and save Excel file with pathway analysis results
    try:
        # Convert pathway_stats to DataFrame for Excel export
        pathway_data = []
        for pathway_name, stats in pathway_stats.items():
            pathway_data.append({
                'Pathway_Name': pathway_name,
                'Status': stats['status'],
                'Z_Score': stats['z_score'],
                'Combined_P_Value': stats['combined_pvalue'],
                'Num_Metabolites': stats['n_metabolites'],
                'Mean_Log2FC': stats['mean_log2fc'],
                'Metabolites': '; '.join(stats['metabolites'])
            })
        
        pathway_df = pd.DataFrame(pathway_data)
        
        # Convert metabolites to DataFrame for Excel export
        metabolite_data = []
        for metabolite_name, data in metabolites.items():
            metabolite_data.append({
                'Metabolite_Name': metabolite_name,
                'Log2FC': data['log2FC'],
                'P_Value': data['pvalue']
            })
        
        metabolite_df = pd.DataFrame(metabolite_data)
        
        # Prepare enzyme/transporter/disease data if available from network creation
        # Recompute the network to get enzyme/transporter/disease stats (we already built them earlier but not preserved here)
        G_all, metabolites_all, pathways_all = create_metabolite_pathway_network(df, pos_set_zscore, 
                                                                                 neg_set_zscore, max_pvalue, 
                                                                                 min_metabolites, 
                                                                                 max_total_pathways=max_total_pathways, 
                                                                                 max_pos=max_pos, max_neg=max_neg)

        # Extract enzyme, transporter, disease info from graph nodes
        enzyme_rows = []
        transporter_rows = []
        disease_rows = []
        for n, d in G_all.nodes(data=True):
            if d.get('node_type') == 'enzyme':
                enzyme_rows.append({
                    'Enzyme': n,
                    'Genes': ','.join(d.get('genes', [])) if d.get('genes') else '',
                    'Accessions': ','.join([a for a in d.get('accessions', []) if a]) if d.get('accessions') else '',
                    'Roles_by_metabolite': '; '.join([f"{m}:{r}" for m, r in d.get('roles', {}).items()]) if d.get('roles') else '',
                    'Num_Metabolites': d.get('n_metabolites', 0),
                    'Concat': d.get('concat', '')
                })
            elif d.get('node_type') == 'transporter':
                transporter_rows.append({
                    'Transporter': n,
                    'Genes': ','.join(d.get('genes', [])) if d.get('genes') else '',
                    'Accessions': ','.join(d.get('accessions', [])) if d.get('accessions') else '',
                    'Num_Metabolites': d.get('n_metabolites', 0),
                    'Concat': d.get('concat', '')
                })
            elif d.get('node_type') == 'disease':
                disease_rows.append({
                    'Disease': n,
                    'Num_Metabolites': d.get('n_metabolites', 0),
                    'Metabolites': '; '.join(d.get('metabolites', [])) if d.get('metabolites') else '',
                    'Concat': d.get('concat', '')
                })

        enzyme_df = pd.DataFrame(enzyme_rows)
        transporter_df = pd.DataFrame(transporter_rows)
        disease_df = pd.DataFrame(disease_rows)

        # Save to Excel with multiple sheets
        with pd.ExcelWriter('pathway_groups_analysis.xlsx', engine='openpyxl') as writer:
            pathway_df.to_excel(writer, sheet_name='Pathway_Analysis', index=False)
            metabolite_df.to_excel(writer, sheet_name='Metabolite_Data', index=False)
            if not enzyme_df.empty:
                enzyme_df.to_excel(writer, sheet_name='Enzymes', index=False)
            if not transporter_df.empty:
                transporter_df.to_excel(writer, sheet_name='Transporters', index=False)
            if not disease_df.empty:
                disease_df.to_excel(writer, sheet_name='Associated_Diseases', index=False)
            # # Add enrichment results
            # if enrichment_df is not None and not enrichment_df.empty:
            #     enrichment_df.to_excel(writer, sheet_name='Pathway_Enrichment', index=False)
        
        print(f"Excel file successfully saved: pathway_groups_analysis.xlsx")
        
    except Exception as e:
        print(f"Error saving Excel file: {e}")
    
    print(f"\nFinal Summary:")
    print(f"- Total unique metabolites across all groups: {len(metabolites)}")
    print(f"- Total unique pathways across all groups: {len(pathway_stats)}")
    print(f"- Excel file saved: pathway_groups_analysis.xlsx")

def export_network_to_graphml(df, pathway_stats, output_file="network.graphml", 
                               include_upstream=True, include_diseases=True,
                               upstream_data=None, disease_data=None,
                               pos_set_fc=0, neg_set_fc=0):
    """
    Export metabolite-pathway network as GraphML for use in Cytoscape, Gephi, or yEd.
    
    USES EXACT SAME COLOR SCHEME AS HTML NETWORK VISUALIZATION
    
    This creates a fully interactive network file that can be opened in professional
    network visualization tools for:
    - Dragging and repositioning nodes
    - Changing colors, fonts, sizes interactively
    - Applying different layout algorithms
    - Publication-quality exports (SVG, PNG, PDF)
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Metabolite data with columns including Name, log2FC, pvalue, pathways
    pathway_stats : dict
        Dictionary of pathway statistics (FILTERED - only selected pathways)
    output_file : str
        Output GraphML filename
    include_upstream : bool
        Include enzyme/transporter nodes
    include_diseases : bool
        Include disease nodes
    upstream_data : pandas.DataFrame, optional
        DataFrame with upstream regulator data (filtered by user selection)
    disease_data : pandas.DataFrame, optional
        DataFrame with disease data (filtered by user selection)
    pos_set_fc : float
        Positive fold change threshold
    neg_set_fc : float
        Negative fold change threshold
    
    Returns:
    --------
    str : Path to saved GraphML file
    
    Usage in Cytoscape:
    -------------------
    1. File → Import → Network from File
    2. Select the generated .graphml file
    3. Layout → Apply layout (e.g., Force-Directed)
    4. Edit nodes/edges interactively
    5. File → Export as Image (SVG/PNG/PDF)
    """
    
    # STANDARDIZED COLOR SCHEME - MATCHES HTML NETWORK EXACTLY
    PATHWAY_COLORS = {
        'Activated': 'red',         # Pathways: red for activated
        'Inhibited': '#3498DB',     # Blue
        'No Change': '#BDC3C7'      # Gray
    }
    
    METABOLITE_COLORS = {
        'Upregulated': '#FF8C00',   # Orange metabolite nodes
        'Downregulated': 'green',   # Green
        'No Change': '#cccccc'      # Light gray
    }
    
    UPSTREAM_COLOR = '#8e44ad'   # Purple (matches network)
    DISEASE_COLOR = '#8e44ad'    # Purple (matches upstream)
    
    EDGE_COLORS = {
        'activating': '#FF0000',   # Red edges for activating (dotted)
        'inhibiting': '#3498DB',   # Blue
        # 'neutral': '#808080',    # Gray (commented out everywhere)
        'associated': '#888888'    # Gray for upstream/disease associations
    }
    
    G = nx.Graph()
    
    # Track what we're adding
    metabolite_count = 0
    pathway_count = 0
    upstream_count = 0
    disease_count = 0
    
    print(f"\n🔄 Creating GraphML network for interactive editing...")
    print(f"   • Filtering: Only selected pathways, metabolites, upstream, and diseases")
    
    # CRITICAL: Detect the name column (may be 'LipidID' for lipid mode)
    name_candidates = ['Name', 'LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite', 'name', 'Metabolite_Name']
    name_col = next((c for c in name_candidates if c in df.columns), 'Name')
    print(f"   • Using column '{name_col}' for metabolite names")
    
    # Get list of metabolites that are in the selected pathways
    metabolites_in_pathways = set()
    for pathway_name, pathway_data in pathway_stats.items():
        metabolites_in_pathways.update(pathway_data.get('metabolites', []))
    
    print(f"   • {len(metabolites_in_pathways)} metabolites connected to selected pathways")
    
    # FOCUSED MODE: Do NOT include additional metabolites from upstream/diseases
    # Only metabolites that are connected to pathways are included
    # This ensures diseases/upstream only connect to pathway-linked metabolites
    print(f"   • FOCUSED MODE: Only including pathway-connected metabolites")
    
    # Add metabolite nodes (ONLY metabolites in selected pathways)
    for _, row in df.iterrows():
        # CRITICAL FIX: Use detected name_col instead of hardcoded 'Name'
        met_name = row.get(name_col, row.get('Name', 'Unknown'))
        
        # CRITICAL: Only include metabolites that are in selected pathways
        if met_name not in metabolites_in_pathways:
            continue
            
        log2fc = row.get('log2FC', 0)
        pvalue = row.get('pvalue', 1)
        
        # Determine metabolite status using SAME logic as network visualization
        if pd.notna(log2fc):
            if log2fc > pos_set_fc:
                status = 'Upregulated'
                color = METABOLITE_COLORS['Upregulated']
            elif log2fc < neg_set_fc:
                status = 'Downregulated'
                color = METABOLITE_COLORS['Downregulated']
            else:
                status = 'No Change'
                color = METABOLITE_COLORS['No Change']
        else:
            status = 'No Change'
            color = METABOLITE_COLORS['No Change']
        
        G.add_node(met_name,
                   node_type='metabolite',
                   shape='ellipse',
                   log2FC=float(log2fc) if pd.notna(log2fc) else 0.0,
                   pvalue=float(pvalue) if pd.notna(pvalue) else 1.0,
                   status=status,
                   color=color,
                   size=30)
        metabolite_count += 1
    
    print(f"   • Added {metabolite_count} metabolite nodes")
    
    # Add pathway nodes and edges to metabolites (ONLY selected pathways)
    for pathway_name, pathway_data in pathway_stats.items():
        # Determine pathway status using EXACT same colors as network
        status = pathway_data.get('status', 'No Change')
        color = PATHWAY_COLORS.get(status, PATHWAY_COLORS['No Change'])
        
        G.add_node(pathway_name,
                   node_type='pathway',
                   shape='diamond',
                   status=status,
                   combined_pvalue=float(pathway_data.get('combined_pvalue', 1)),
                   z_score=float(pathway_data.get('z_score', 0)),
                   n_metabolites=int(pathway_data.get('n_metabolites', 0)),
                   color=color,
                   size=40)
        pathway_count += 1
        
        # Add edges from pathway to its metabolites
    # EDGE COLOR BASED ON METABOLITE REGULATION STATUS ONLY
    # Upregulated metabolite = ORANGE edge (dotted line)
    # Downregulated metabolite = BLUE edge (dotted line)
        # Neutral metabolites are NOT connected (no gray lines)
        for met in pathway_data.get('metabolites', []):
            if G.has_node(met):
                # Get metabolite status from node attributes
                met_status = G.nodes[met].get('status', 'No Change')
                
                # Only create edges for upregulated or downregulated metabolites
                if met_status == 'Upregulated':
                    edge_color = '#FF8C00'  # Upregulated metabolite = orange edge (dotted)
                    edge_type = 'upregulated'
                    line_type = 'dotted'
                    
                    G.add_edge(pathway_name, met,
                               edge_type=edge_type,
                               interaction=edge_type,
                               color=edge_color,
                               line_type=line_type)
                               
                elif met_status == 'Downregulated':
                    edge_color = '#3498DB'  # Downregulated metabolite = blue edge
                    edge_type = 'downregulated'
                    line_type = 'dotted'
                    
                    G.add_edge(pathway_name, met,
                               edge_type=edge_type,
                               interaction=edge_type,
                               color=edge_color,
                               line_type=line_type)
                # else:
                #     # Neutral/No Change - COMMENTED OUT (no neutral edges shown)
                #     pass
    
    print(f"   • Added {pathway_count} pathway nodes")
    
    # Add upstream regulators (ONLY those in upstream_data if provided)
    if include_upstream and upstream_data is not None and not upstream_data.empty:
        for _, row in upstream_data.iterrows():
            regulator_name = row.get('Name', '')
            regulator_type = row.get('Type', 'Unknown')
            associated_mets = row.get('Associated_Metabolites', '')
            
            if not regulator_name:
                continue
            
            # Parse associated metabolites (support '|', ',', ';' delimiters)
            if pd.notna(associated_mets) and str(associated_mets).strip():
                met_list = [m.strip() for m in re.split(r"[|,;]", str(associated_mets)) if m.strip()]
                
                # Only add if at least one metabolite is in our network
                connected_mets = [m for m in met_list if G.has_node(m)]
                if not connected_mets:
                    continue
                
                # Add upstream node if not already present
                if not G.has_node(regulator_name):
                    G.add_node(regulator_name,
                              node_type='upstream',
                              regulator_type=regulator_type,
                              shape='rectangle',
                              color=UPSTREAM_COLOR,
                              size=25)
                    upstream_count += 1
                
                # Add edges to connected metabolites
                for met in connected_mets:
                    G.add_edge(regulator_name, met,
                               edge_type='association',
                               interaction='regulates',
                               color=EDGE_COLORS['associated'],
                               line_type='dotted')
    
    if upstream_count > 0:
        print(f"   • Added {upstream_count} upstream regulator nodes")
    
    # Add diseases (ONLY those in disease_data if provided)
    if include_diseases and disease_data is not None and not disease_data.empty:
        for _, row in disease_data.iterrows():
            disease_name = row.get('Disease', '') or row.get('Disease Name', '')
            associated_mets = row.get('Associated_Metabolites', '')
            
            if not disease_name:
                continue
            
            # Parse associated metabolites (support '|', ',', ';' delimiters)
            if pd.notna(associated_mets) and str(associated_mets).strip():
                met_list = [m.strip() for m in re.split(r"[|,;]", str(associated_mets)) if m.strip()]
                
                # Only add if at least one metabolite is in our network
                connected_mets = [m for m in met_list if G.has_node(m)]
                if not connected_mets:
                    continue
                
                # Add disease node if not already present
                if not G.has_node(disease_name):
                    G.add_node(disease_name,
                              node_type='disease',
                              shape='hexagon',
                              color=DISEASE_COLOR,
                              size=30)
                    disease_count += 1
                
                # Add edges to connected metabolites
                for met in connected_mets:
                    G.add_edge(disease_name, met,
                               edge_type='association',
                               interaction='associated_with',
                               color=EDGE_COLORS['associated'],
                               line_type='dotted')
    
    # Add diseases (ONLY those in disease_data if provided)
    if disease_count > 0:
        print(f"   • Added {disease_count} disease nodes")
    
    # -----------------------------------------------------------------
    # ✅ Add Cytoscape-compatible visual attributes
    # -----------------------------------------------------------------
    for n, d in G.nodes(data=True):
        # Colors
        if "color" in d:
            d["fillColor"] = d["color"]
            d["node_color"] = d["color"]
        
        # Node size (explicit)
        if "size" in d:
            d["node_size"] = d["size"]
        
        # Node shape mapping for Cytoscape
        shape = d.get("shape", "ellipse")
        if shape == "diamond":
            cy_shape = "diamond"
        elif shape == "rectangle":
            cy_shape = "rectangle"
        elif shape == "hexagon":
            cy_shape = "hexagon"
        else:
            cy_shape = "ellipse"
        d["cy_shape"] = cy_shape

    # Edges - add Cytoscape-compatible attributes
    for u, v, d in G.edges(data=True):
        if "color" in d:
            d["lineColor"] = d["color"]
            d["edge_color"] = d["color"]
            d["EDGE_PAINT"] = d["color"]  # Additional Cytoscape attribute
        d["edge_width"] = 2.0  # uniform edge thickness for now
        
        # Map line_type to Cytoscape line style
        line_type = d.get("line_type", "solid")
        if line_type == "dotted":
            d["EDGE_LINE_TYPE"] = "DOT"
            d["cy_line_style"] = "DOT"
        else:
            d["EDGE_LINE_TYPE"] = "SOLID"
            d["cy_line_style"] = "SOLID"
    # -----------------------------------------------------------------

    # Save as GraphML
    nx.write_graphml(G, output_file)

    
    print(f"✅ GraphML network exported: {output_file}")
    print(f"\n📊 Network Statistics:")
    print(f"   • Metabolites: {metabolite_count}")
    print(f"   • Pathways: {pathway_count}")
    if upstream_count > 0:
        print(f"   • Upstream Regulators: {upstream_count}")
    if disease_count > 0:
        print(f"   • Diseases: {disease_count}")
    print(f"   • Total nodes: {G.number_of_nodes()}")
    print(f"   • Total edges: {G.number_of_edges()}")
    
    print(f"\n🎨 Color Scheme (matches HTML network exactly):")
    print(f"   PATHWAY NODES:")
    print(f"   • Activated Pathways: Orange ({PATHWAY_COLORS['Activated']})")
    print(f"   • Inhibited Pathways: Blue ({PATHWAY_COLORS['Inhibited']})")
    print(f"   • No Change Pathways: Gray ({PATHWAY_COLORS['No Change']})")
    print(f"\n   METABOLITE NODES:")
    print(f"   • Upregulated Metabolites: Orange")
    print(f"   • Downregulated Metabolites: Green")
    print(f"\n   EDGES (based on metabolite regulation):")
    print(f"   • Upregulated Metabolite → Pathway: Orange edge")
    print(f"   • Downregulated Metabolite → Pathway: Green edge")
    # print(f"   • Neutral Metabolite → Pathway: Gray edge (currently not used)")
    print(f"\n   OTHER NODES:")
    print(f"   • Upstream Regulators: Purple ({UPSTREAM_COLOR})")
    print(f"   • Diseases: Purple ({DISEASE_COLOR})")
    
    print(f"\n🎨 To edit interactively in Cytoscape:")
    print(f"   1. Download Cytoscape: https://cytoscape.org/download.html")
    print(f"   2. Open Cytoscape → File → Import → Network from File")
    print(f"   3. Select: {output_file}")
    print(f"   4. Apply edge styles:")
    print(f"      • Style panel → Edge → Stroke Color → Passthrough Mapping from 'color'")
    print(f"      • Style panel → Edge → Line Type → Discrete Mapping from 'line_type':")
    print(f"        - dotted → DOT")
    print(f"        - solid → SOLID")
    print(f"   5. Layout → Apply Layout → Force Directed")
    print(f"   6. Export: File → Export as Image (SVG/PNG/PDF)")
    
    return output_file

def export_network_to_cytoscape(df, pathway_stats, output_file="network.graphml", 
                               include_upstream=True, include_diseases=True,
                               upstream_data=None, disease_data=None,
                               pos_set_fc=0, neg_set_fc=0,
                               compress_by_class=True):
    """
    Create and export the metabolite–pathway network directly to Cytoscape using its API.
    Preserves colors, node types, and edges exactly as in the HTML network visualization.
    Falls back to GraphML if Cytoscape (or py4cytoscape) is not available.
    """

    import pandas as pd
    import networkx as nx
    try:
        import py4cytoscape as p4c
    except Exception as e:
        print(f"⚠️ py4cytoscape not available ({e}); exporting GraphML instead.")
        return export_network_to_graphml(df, pathway_stats, output_file=output_file,
                         include_upstream=include_upstream,
                         include_diseases=include_diseases,
                         upstream_data=upstream_data,
                         disease_data=disease_data,
                         pos_set_fc=pos_set_fc, neg_set_fc=neg_set_fc)

    # Standardized colors to match HTML network
    PATHWAY_COLORS = {
        'Activated': 'red',      # Pathways: red for activated
        'Inhibited': '#3498DB',  # Blue
        'No Change': '#BDC3C7'
    }
    METABOLITE_COLORS = {
        'Upregulated': '#FF8C00',   # Orange metabolite nodes
        'Downregulated': 'green',
        'No Change': '#cccccc'
    }
    UPSTREAM_COLOR = '#8e44ad'
    DISEASE_COLOR = '#8e44ad'  # Purple (matches upstream)
    EDGE_COLORS = {
        'activating': '#FF0000',   # Red edges for activating (upregulated)
        'inhibiting': '#228B22',   # Green for inhibiting (downregulated)
        'neutral': '#808080',
        'associated': '#888888'    # Gray for upstream/disease associations
    }

    # Build graph nodes/edges similar to GraphML export
    G = nx.Graph()

    # CRITICAL: Use ONLY metabolites present in the filtered df
    # The df has already been filtered by the GUI to include only selected metabolites
    # Support multiple column name variations
    name_candidates = ['Name', 'LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite', 'name', 'Metabolite_Name']
    name_col = next((c for c in name_candidates if c in df.columns), None)
    available_metabolites = set(df[name_col].dropna().tolist()) if name_col else set()
    
    print(f"   • {len(available_metabolites)} metabolites available in filtered data (using column '{name_col}')")

    # Simplified solid colors for Cytoscape (no shading variations)
    # Colors are based on status only, using z-score/fold-change for size instead
    def get_metabolite_color(fc, pos_threshold, neg_threshold):
        """Return solid color based on metabolite regulation status"""
        try:
            v = float(fc)
        except Exception:
            v = 0.0
        if v > pos_threshold:
            return '#FF8C00'  # Solid orange for upregulated
        elif v < neg_threshold:
            return '#228B22'  # Solid forest green for downregulated
        else:
            return '#CCCCCC'  # Gray for no change
    
    def get_pathway_color(status):
        """Return solid color based on pathway status"""
        if status == 'Activated':
            return '#FF0000'  # Solid red for activated
        elif status == 'Inhibited':
            return '#3498DB'  # Solid blue for inhibited
        else:
            return '#BDC3C7'  # Gray for no change

    # ═══════════════════════════════════════════════════════════════════════════
    # FOCUSED MODE: First, get all metabolites that are in selected pathways
    # Only these metabolites should be included in the network
    # ═══════════════════════════════════════════════════════════════════════════
    metabolites_in_pathways = set()
    for pname, pdata in pathway_stats.items():
        metabolites_in_pathways.update(pdata.get('metabolites', []))
    
    print(f"   • {len(metabolites_in_pathways)} metabolites connected to selected pathways (FOCUSED MODE)")

    # Metabolite nodes - ONLY include metabolites that are in pathways
    for _, row in df.iterrows():
        # CRITICAL FIX: Use detected name_col instead of hardcoded 'Name'
        met = row.get(name_col, row.get('Name', 'Unknown')) if name_col else row.get('Name', 'Unknown')
        if met not in available_metabolites:
            continue
        
        # FOCUSED MODE: Only include metabolites that are connected to pathways
        if met not in metabolites_in_pathways:
            continue
            
        log2fc = row.get('log2FC', 0)
        pval = row.get('pvalue', 1)
        if pd.notna(log2fc):
            if log2fc > pos_set_fc:
                status = 'Upregulated'
            elif log2fc < neg_set_fc:
                status = 'Downregulated'
            else:
                status = 'No Change'
        else:
            status = 'No Change'
            log2fc = 0.0
        
        # Use solid colors based on status
        color = get_metabolite_color(log2fc, pos_set_fc, neg_set_fc)
        G.add_node(met, node_type='metabolite', color=color, status=status,
                   log2FC=float(log2fc) if pd.notna(log2fc) else 0.0,
                   pvalue=float(pval) if pd.notna(pval) else 1.0)

    # Pathway nodes and edges
    print("\n🎨 Pathway Color Assignment (Solid Colors):")
    for pname, pdata in pathway_stats.items():
        status = pdata.get('status', 'No Change')
        # Use solid color based on status only
        color = get_pathway_color(status)
        print(f"  {pname[:60]}: status='{status}' → color='{color}'")
        G.add_node(pname, node_type='pathway', color=color, status=status,
                   combined_pvalue=float(pdata.get('combined_pvalue', 1)),
                   z_score=float(pdata.get('z_score', 0)),
                   n_metabolites=int(pdata.get('n_metabolites', 0)))
        for m in pdata.get('metabolites', []):
            if G.has_node(m):
                # Edge type depends ONLY on metabolite regulation status
                met_status = G.nodes[m].get('status', 'No Change')
                if met_status == 'Upregulated':
                    etype, ecolor = 'activating', '#FF0000'  # Red for activating edges
                    G.add_edge(pname, m, edge_type=etype, color=ecolor)
                elif met_status == 'Downregulated':
                    etype, ecolor = 'inhibiting', '#228B22'  # Green for inhibiting edges
                    G.add_edge(pname, m, edge_type=etype, color=ecolor)
                # else: neutral commented out

    # FOCUSED MODE: Identify metabolites that have pathway connections
    # This is used to filter upstream/disease connections
    metabolites_with_pathway = set()
    for pname, pdata in pathway_stats.items():
        for m in pdata.get('metabolites', []):
            if G.has_node(m):
                metabolites_with_pathway.add(m)
    
    print(f"\n🎯 FOCUSED MODE: {len(metabolites_with_pathway)} metabolites have pathway connections")

    # Upstream nodes and edges - FOCUSED: only connect to metabolites with pathway connections
    if include_upstream and upstream_data is not None and not upstream_data.empty:
        for _, r in upstream_data.iterrows():
            reg = r.get('Name', '')
            associated_mets = r.get('Associated_Metabolites', '')
            mets = []
            if pd.notna(associated_mets) and str(associated_mets).strip():
                mets = [m.strip() for m in re.split(r"[|,;]", str(associated_mets)) if m.strip()]
            # FOCUSED MODE: Only include metabolites that exist in graph AND have pathway connections
            mets = [m for m in mets if G.has_node(m) and m in metabolites_with_pathway]
            if reg and mets:
                G.add_node(reg, node_type='upstream', color=UPSTREAM_COLOR)
                for m in mets:
                    G.add_edge(reg, m, edge_type='association', color=EDGE_COLORS['associated'])

    # Disease nodes and edges - FOCUSED: only connect to metabolites with pathway connections
    if include_diseases and disease_data is not None and not disease_data.empty:
        for _, r in disease_data.iterrows():
            dis = r.get('Disease', '') or r.get('Disease Name', '')
            associated_mets = r.get('Associated_Metabolites', '')
            mets = []
            if pd.notna(associated_mets) and str(associated_mets).strip():
                mets = [m.strip() for m in re.split(r"[|,;]", str(associated_mets)) if m.strip()]
            # FOCUSED MODE: Only include metabolites that exist in graph AND have pathway connections
            mets = [m for m in mets if G.has_node(m) and m in metabolites_with_pathway]
            if dis and mets:
                G.add_node(dis, node_type='disease', color=DISEASE_COLOR)
                for m in mets:
                    G.add_edge(dis, m, edge_type='association', color=EDGE_COLORS['associated'])

    # Convert to data frames
    nodes_df = pd.DataFrame([{'id': n, **G.nodes[n]} for n in G.nodes])
    edges_df = pd.DataFrame([{'source': u, 'target': v, **G.edges[u, v]} for u, v in G.edges])

    # Try Cytoscape
    try:
        p4c.cytoscape_ping()
        res = p4c.create_network_from_data_frames(nodes_df, edges_df, title='Metabolite-Pathway Network')
        net_suid = res.get('networkSUID') if isinstance(res, dict) else None

        # Create a visual style to map color and shape
        style_name = 'MetabolitePathwayStyle'
        try:
            # Delete if exists
            p4c.delete_visual_style(style_name)
        except Exception:
            pass

        # Create style with defaults
        try:
            p4c.create_visual_style(style_name)
        except Exception as e:
            print(f"⚠️ Could not create style: {e}")
            # Fall back to GraphML
            return export_network_to_graphml(df, pathway_stats, output_file=output_file,
                                           include_upstream=include_upstream,
                                           include_diseases=include_diseases,
                                           upstream_data=upstream_data,
                                           disease_data=disease_data,
                                           pos_set_fc=pos_set_fc, neg_set_fc=neg_set_fc)
        
        # Set the new style as current
        p4c.set_visual_style(style_name)
        
        # Set default values
        p4c.set_node_shape_default('ELLIPSE')
        p4c.set_node_size_default(30)
        p4c.set_node_color_default('#cccccc')
        p4c.set_edge_color_default('#888888')
        p4c.set_edge_line_width_default(2.0)
        
        # Map node shapes by node_type using discrete mapping
        try:
            p4c.set_node_shape_mapping('node_type', 
                                      ['metabolite', 'pathway', 'upstream', 'disease'],
                                      ['ELLIPSE', 'DIAMOND', 'TRIANGLE', 'HEXAGON'])
            print("✅ Node shape mapping applied")
        except Exception as e:
            print(f"⚠️ Node shape mapping failed: {e}")
        
        # Map node colors using passthrough mapping (column value = color)
        try:
            p4c.set_node_color_mapping('color', mapping_type='p')
            print("✅ Node color mapping applied")
        except Exception as e:
            print(f"⚠️ Node color mapping failed: {e}")
        
        # Map edge colors using passthrough mapping
        try:
            p4c.set_edge_color_mapping('color', mapping_type='p')
            print("✅ Edge color mapping applied")
        except Exception as e:
            print(f"⚠️ Edge color mapping failed: {e}")
        
        # Map edge line style by edge_type (without mapping_type parameter)
        try:
            p4c.set_edge_line_style_mapping('edge_type',
                                           ['activating', 'inhibiting', 'association'],
                                           ['DOT', 'DOT', 'DOT'])
            print("✅ Edge line style mapping applied")
        except Exception as e:
            print(f"⚠️ Edge line style mapping failed: {e}")

        print("✅ Network successfully created and styled in Cytoscape.")
        print(f"📊 Network contains:")
        print(f"   • {len([n for n in G.nodes if G.nodes[n].get('node_type') == 'metabolite'])} metabolites")
        print(f"   • {len([n for n in G.nodes if G.nodes[n].get('node_type') == 'pathway'])} pathways")
        print(f"   • {len(G.edges)} edges")
        print(f"\n⚠️ IMPORTANT: Visual style 'MetabolitePathwayStyle' has been applied")
        print(f"   Colors are mapped from the 'color' column in node/edge tables")
        print("💾 To save: File → Save Session As → .cys file")
        return output_file
    except Exception as e:
        print(f"⚠️ Cytoscape connection failed ({e}); exporting GraphML instead.")
        return export_network_to_graphml(df, pathway_stats, output_file=output_file,
                                         include_upstream=include_upstream,
                                         include_diseases=include_diseases,
                                         upstream_data=upstream_data,
                                         disease_data=disease_data,
                                         pos_set_fc=pos_set_fc, neg_set_fc=neg_set_fc)

def main():
    """Main function with configuration parameters"""
    
    # ===================== CONFIGURATION SECTION =====================
    # Data import parameters
    excel_file = "metabolite_pathways_annotated.xlsx"
    sheet_name = "Metabolites_with_Pathways"
    pvalue_threshold = 0.05
    pos_set_fc = 0
    neg_set_fc = 0
    
    # Pathway analysis parameters
    pos_set_zscore = 0
    neg_set_zscore = 0
    min_metabolites = 1
    max_pvalue = 0.05
    max_total_pathways = 50
    max_pos = 25
    max_neg = 25
    
    # Network visualization parameters
    node_spacing = 1.0          # Lower = more compact (0.5-3.0 recommended)
    layout_iterations = 100      # Higher = better positioning (50-300 recommended)
    
    # Network plot dimensions
    network_width = 1800        # Width of network plot in pixels (increased to prevent cutoff)
    network_height = 1000       # Height of network plot in pixels (increased to prevent cutoff)
    
    # Font size controls
    metabolite_font_size = 12   # Font size for metabolite labels (8-16 recommended)
    pathway_font_size = 12       # Font size for pathway labels (8-14 recommended)  
    enzyme_font_size = 8        # Font size for enzyme/transporter labels (6-12 recommended)
    
    # Text wrapping controls - Conservative settings to avoid unnecessary breaking
    metabolite_max_chars = 30   # Max characters per line for metabolite names (allow longer single-line names)
    metabolite_max_lines = 1    # Max lines for metabolite text wrapping (1 line to preserve chemical names)
    pathway_max_chars = 30      # Max characters per line for pathway names (allow longer before breaking)
    pathway_max_lines = 3       # Max lines for pathway text wrapping (2-4 recommended)
    

    # ===================== END CONFIGURATION =====================
    
    # Call run_visualizations with all parameters
    run_visualizations(
        excel_file=excel_file,
        sheet_name=sheet_name,
        pvalue_threshold=pvalue_threshold,
        pos_set_fc=pos_set_fc,
        neg_set_fc=neg_set_fc,
        pos_set_zscore=pos_set_zscore,
        neg_set_zscore=neg_set_zscore,
        min_metabolites=min_metabolites,
        max_pvalue=max_pvalue,
        max_total_pathways=max_total_pathways,
        max_pos=max_pos,
        max_neg=max_neg,
        node_spacing=node_spacing,
        layout_iterations=layout_iterations,
        network_width=network_width,
        network_height=network_height,
        metabolite_font_size=metabolite_font_size,
        pathway_font_size=pathway_font_size,
        enzyme_font_size=enzyme_font_size,
        metabolite_max_chars=metabolite_max_chars,
        metabolite_max_lines=metabolite_max_lines,
        pathway_max_chars=pathway_max_chars,
        pathway_max_lines=pathway_max_lines
    )

def create_network_legend_figure(include_upstream=False, include_diseases=False):
    """
    Create a Plotly figure showing the legend for pathway network visualization.
    
    Args:
        include_upstream: Whether upstream pathways are included
        include_diseases: Whether disease pathways are included
    
    Returns:
        plotly.graph_objects.Figure: A figure with legend information
    """
    import plotly.graph_objects as go
    
    # Create legend text
    legend_text = "Network Legend<br><br>"
    
    legend_text += "<b>Node Colors:</b><br>"
    legend_text += "🟢 Green: Upregulated (positive log2FC)<br>"
    legend_text += "🔴 Red: Downregulated (negative log2FC)<br>"
    legend_text += "⚪ Gray: No data/neutral<br><br>"
    
    legend_text += "<b>Node Size:</b><br>"
    legend_text += "Larger nodes = higher statistical significance<br><br>"
    
    legend_text += "<b>Connections:</b><br>"
    legend_text += "Edges show pathway-metabolite relationships<br>"
    
    if include_upstream:
        legend_text += "🟦 Blue nodes: Upstream pathways<br>"
    
    if include_diseases:
        legend_text += "🟥 Pink nodes: Disease-associated pathways<br>"
    
    # Create a simple Plotly figure with legend text
    fig = go.Figure()
    
    # Add annotation with legend text
    fig.add_annotation(
        text=legend_text,
        xref="paper", yref="paper",
        x=0.5, y=0.5,
        showarrow=False,
        font=dict(size=12),
        align="left",
        bgcolor="rgba(240, 240, 240, 0.8)",
        bordercolor="black",
        borderwidth=1,
        borderpad=10
    )
    
    # Configure layout
    fig.update_layout(
        title="Pathway Network Legend",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        plot_bgcolor="white",
        paper_bgcolor="white",
        width=520,
        height=360,
        margin=dict(l=20, r=20, t=40, b=20)
    )
    
    return fig


if __name__ == "__main__":
    main()











