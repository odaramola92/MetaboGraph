#!/usr/bin/env python3
"""
Metabolite ID Annotator - Phase 1
Focused on ID annotation using API-first approach (PubChem/KEGG) then offline databases
"""

import pandas as pd
import numpy as np
import requests
import time
import pickle
import os
import logging
import re
import threading
from typing import Dict, List, Optional, Any
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def fetch_pubchem_cross_references(pubchem_cid: int, max_retries: int = 3, delay: float = 0.25) -> Dict[str, str]:
    """
    Fetch cross-references (KEGG, HMDB, ChEBI, CAS) for a PubChem CID.
    
    Args:
        pubchem_cid: The PubChem Compound ID
        max_retries: Maximum number of retry attempts
        delay: Delay between requests in seconds
        
    Returns:
        Dictionary with cross-reference IDs: {'KEGG': 'CLP104', 'HMDB': 'HMDB0007970', ...}
    """
    xrefs = {'KEGG': '', 'HMDB': '', 'ChEBI': '', 'CAS': '', 'LipidMaps': ''}
    
    # Use RegistryID endpoint to get all cross-references
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_cid}/xrefs/RegistryID/JSON"
    
    for attempt in range(max_retries):
        try:
            time.sleep(delay)  # Rate limiting
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                if 'InformationList' in data and 'Information' in data['InformationList']:
                    info = data['InformationList']['Information']
                    if info and len(info) > 0 and 'RegistryID' in info[0]:
                        registry_ids = info[0]['RegistryID']
                        
                        # Parse the registry IDs to find cross-references
                        for reg_id in registry_ids:
                            reg_id_str = str(reg_id).strip()
                            
                            # KEGG IDs (format: Cxxxxx or CLPxxx)
                            # First clean the ID to remove any cpd: or cpd+ prefixes
                            cleaned_reg_id = clean_kegg_id(reg_id_str)
                            
                            # Standard KEGG compound IDs: C followed by exactly 5 digits (C00001-C99999)
                            if cleaned_reg_id.startswith('C') and len(cleaned_reg_id) == 6:
                                # Check if it's a KEGG compound ID (C followed by 5 digits)
                                if cleaned_reg_id[1:].isdigit():
                                    xrefs['KEGG'] = cleaned_reg_id
                                    continue
                            
                            # KEGG lipid IDs: CLP followed by digits (CLP001-CLP999...)
                            if cleaned_reg_id.startswith('CLP') and len(cleaned_reg_id) >= 6:
                                if cleaned_reg_id[3:].isdigit():
                                    xrefs['KEGG'] = cleaned_reg_id
                                    continue
                            
                            # HMDB IDs
                            if reg_id_str.startswith('HMDB'):
                                xrefs['HMDB'] = reg_id_str
                                continue
                            
                            # ChEBI IDs
                            if reg_id_str.startswith('CHEBI:'):
                                xrefs['ChEBI'] = reg_id_str
                                continue
                            
                            # LipidMaps IDs
                            if reg_id_str.startswith('LM'):
                                xrefs['LipidMaps'] = reg_id_str
                                continue
                            
                            # CAS numbers (format: xxx-xx-x)
                            if '-' in reg_id_str and reg_id_str.replace('-', '').replace(' ', '').isdigit():
                                # Basic CAS validation: should have 2-3 dashes and be numeric when dashes removed
                                parts = reg_id_str.split('-')
                                if len(parts) >= 3 and all(part.isdigit() for part in parts):
                                    xrefs['CAS'] = reg_id_str
                                    continue
                        
                break
            elif response.status_code == 404:
                # No cross-references found
                break
            else:
                logger.warning(f"PubChem API error for CID {pubchem_cid}: {response.status_code}")
                
        except requests.RequestException as e:
            logger.warning(f"Request failed for CID {pubchem_cid} (attempt {attempt + 1}): {e}")
            if attempt < max_retries - 1:
                time.sleep(delay * (attempt + 1))  # Exponential backoff
                
        except Exception as e:
            logger.error(f"Unexpected error for CID {pubchem_cid}: {e}")
            break
    
    return xrefs

def clean_kegg_id(kegg_id: str) -> str:
    """
    Clean KEGG ID by removing 'cpd:' and 'cpd+' prefixes.
    
    Args:
        kegg_id: The KEGG ID to clean (e.g., 'cpd:C00435', 'cpd+C00435', 'C00435')
        
    Returns:
        Cleaned KEGG ID without prefixes (e.g., 'C00435')
    """
    if not kegg_id or not isinstance(kegg_id, str):
        return kegg_id
    
    kegg_id = kegg_id.strip()
    
    # Remove 'cpd:' prefix (case-insensitive)
    if kegg_id.lower().startswith('cpd:'):
        kegg_id = kegg_id[4:]
    
    # Remove 'cpd+' prefix (case-insensitive)  
    elif kegg_id.lower().startswith('cpd+'):
        kegg_id = kegg_id[4:]
    
    return kegg_id.strip().upper()


def _is_non_empty_value(value: Any) -> bool:
    """Return True when a value is present and not an empty/null-like placeholder."""
    if value is None:
        return False
    try:
        if pd.isna(value):
            return False
    except Exception:
        pass
    text = str(value).strip()
    if not text:
        return False
    return text.lower() not in ('nan', 'none', 'null')


def _normalize_yes_no(value: Any) -> str:
    """Normalize endogenous-style values to strict Yes/No text."""
    if value is None:
        return 'No'
    text = str(value).strip()
    if not text or text.lower() in ('nan', 'none', 'null'):
        return 'No'
    if text.lower() in ('yes', 'y', 'true', '1'):
        return 'Yes'
    if text.lower() in ('no', 'n', 'false', '0'):
        return 'No'
    return 'No'


def _enforce_endogenous_rules(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure Endogenous/Endogenous_Source are strict Yes/No and never blank.

    Rule: if LipidMaps_ID is present for a row, Endogenous must be Yes.
    """
    if df is None or not isinstance(df, pd.DataFrame):
        return df

    out = df.copy()

    if 'Endogenous' not in out.columns:
        if 'Endogenous_Source' in out.columns:
            out['Endogenous'] = out['Endogenous_Source']
        else:
            out['Endogenous'] = 'No'

    out['Endogenous'] = out['Endogenous'].apply(_normalize_yes_no)

    if 'Endogenous_Source' in out.columns:
        out['Endogenous_Source'] = out['Endogenous_Source'].apply(_normalize_yes_no)
    else:
        out['Endogenous_Source'] = out['Endogenous']

    if 'LipidMaps_ID' in out.columns:
        has_lipidmaps = out['LipidMaps_ID'].astype(str).str.strip().ne('') & out['LipidMaps_ID'].notna()
        out.loc[has_lipidmaps, 'Endogenous'] = 'Yes'
        out.loc[has_lipidmaps, 'Endogenous_Source'] = 'Yes'

    return out


def _compute_lipid_msi_level(lipid_id_value: Any) -> str:
    """Classify lipid annotation specificity for lipid-mode outputs.

    Rules requested:
    - Sum composition only (e.g., AcCa(26:1), BisMePA(36:3)) -> Level 3
    - Chain-resolved/molecular species (e.g., BisMePA(19:3_22:6), Cer(d12:2_14:2)) -> Level 2
    """
    if lipid_id_value is None:
        return ''

    text = str(lipid_id_value).strip()
    if not text or text.lower() in ('nan', 'none', 'null'):
        return ''

    match = re.search(r'\(([^)]*)\)', text)
    if not match:
        return ''

    inside = match.group(1).strip()
    if not inside:
        return ''

    # Molecular species / chain-resolved notation commonly uses '_' (or '/').
    if '_' in inside or '/' in inside:
        return 'Level 2'

    # Sum-composition: one total like 36:3 (optionally prefixed like d36:3).
    if re.match(r'^\s*[A-Za-z]*\d+:\d+\s*$', inside):
        return 'Level 3'

    # Conservative default for lipid formula-like parenthetical patterns.
    if re.search(r'\d+:\d+', inside):
        return 'Level 3'

    # Complete single identifier (no chain separation, no ratio notation) -> Level 2
    # Matches: Q10, h26, ggh3444 (simple alphanumeric tokens)
    # Rejects: n_hj_j_11_788_3 (underscores/slashes indicate multiple components)
    if re.match(r'^[A-Za-z0-9\-\.]+$', inside):
        return 'Level 2'

    return ''


def _normalize_header_name(name: Any) -> str:
    """Normalize a column header for tolerant matching."""
    try:
        text = str(name).strip().lower()
    except Exception:
        return ''
    return re.sub(r'[^a-z0-9]+', '', text)


def _resolve_column_by_aliases(columns: List[Any], aliases: List[str]) -> Optional[str]:
    """Find first matching column using exact then normalized alias matching."""
    if not columns:
        return None

    # Exact (case-insensitive, trimmed) pass first.
    col_map_exact = {str(c).strip().lower(): str(c) for c in columns}
    for alias in aliases:
        key = str(alias).strip().lower()
        if key in col_map_exact:
            return col_map_exact[key]

    # Normalized fallback pass.
    col_map_norm = {_normalize_header_name(c): str(c) for c in columns}
    for alias in aliases:
        key = _normalize_header_name(alias)
        if key in col_map_norm:
            return col_map_norm[key]

    return None


def _resolve_confidence_columns(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    """Resolve required columns for MSI/confidence annotation via aliases."""
    cols = list(df.columns)
    return {
        'name': _resolve_column_by_aliases(cols, ['Name', 'Metabolite', 'Molecule', 'LipidID']),
        'ms2': _resolve_column_by_aliases(cols, [
            'MS2', 'MS/MS', 'MS2_Type', 'MS2 Type', 'MS2_Status', 'MS2 Status'
        ]),
        'ms2_purity': _resolve_column_by_aliases(cols, [
            'MS2_Purity', 'MS2 Purity [%]', 'MS2 Purity', 'MS2 Purity %', 'MS2Purity'
        ]),
        'hmdb': _resolve_column_by_aliases(cols, ['HMDB_ID', 'HMDB ID', 'HMDB']),
        'kegg': _resolve_column_by_aliases(cols, ['KEGG_ID', 'KEGG ID', 'KEGG']),
        'lipidmaps': _resolve_column_by_aliases(cols, [
            'LipidMaps_ID', 'LipidMaps ID', 'LipidMaps', 'LMID'
        ]),
    }


def _normalize_ms2_text(ms2_value: Any) -> str:
    try:
        return str(ms2_value).strip().lower()
    except Exception:
        return ''


def _parse_ms2_purity(row: pd.Series, schema_cols: Dict[str, Optional[str]]) -> float:
    """Parse MS2 purity from resolved schema column names."""
    purity_col = schema_cols.get('ms2_purity') if isinstance(schema_cols, dict) else None
    raw = row.get(purity_col, np.nan) if purity_col else np.nan
    if raw is None:
        return np.nan
    if isinstance(raw, str):
        raw = raw.replace('%', '').replace(',', '.').strip()
    try:
        return float(pd.to_numeric(raw, errors='coerce'))
    except Exception:
        return np.nan


def _compute_msi_level_and_confidence(row: pd.Series, schema_cols: Dict[str, Optional[str]]) -> tuple[str, str]:
    """Classify Level 2 evidence and assign MetaboGraph confidence."""
    name_col = schema_cols.get('name') if isinstance(schema_cols, dict) else None
    hmdb_col = schema_cols.get('hmdb') if isinstance(schema_cols, dict) else None
    kegg_col = schema_cols.get('kegg') if isinstance(schema_cols, dict) else None
    lipidmaps_col = schema_cols.get('lipidmaps') if isinstance(schema_cols, dict) else None
    ms2_col = schema_cols.get('ms2') if isinstance(schema_cols, dict) else None

    name_present = _is_non_empty_value(row.get(name_col)) if name_col else False
    curated_db_present = any(
        _is_non_empty_value(row.get(col))
        for col in (hmdb_col, kegg_col, lipidmaps_col)
        if col
    )

    if not (name_present or curated_db_present):
        return 'Level 4', 'low'

    ms2_value = row.get(ms2_col, '') if ms2_col else ''
    ms2_text = _normalize_ms2_text(ms2_value)
    ms2_purity = _parse_ms2_purity(row, schema_cols)

    is_preferred = 'dda for preferred ion' in ms2_text
    is_other = ('dda for other ion' in ms2_text) or ('dda for non-preferred ion' in ms2_text)
    has_ms2_signal = _is_non_empty_value(ms2_value) and ('no ms2' not in ms2_text)
    useful_ms2 = is_preferred or is_other or has_ms2_signal

    # Confidence is only elevated above low when both curated ID support and useful MS2 exist.
    if curated_db_present and useful_ms2:
        # Level 2 (main): preferred ion + good purity + curated DB support.
        if is_preferred and not np.isnan(ms2_purity) and ms2_purity >= 80:
            return 'Level 2', 'high'

        # Level 2 (weaker): curated ID support plus less certain MS2 evidence.
        if is_preferred and (np.isnan(ms2_purity) or ms2_purity < 80):
            return 'Level 2', 'medium'
        if is_other and (np.isnan(ms2_purity) or ms2_purity < 80):
            return 'Level 2', 'medium'
        if is_other and not np.isnan(ms2_purity) and ms2_purity >= 80:
            return 'Level 2', 'medium'
        if has_ms2_signal:
            return 'Level 2', 'medium'

    # Name-only, curated-only, or MS2-only rows remain Level 2 but low confidence.
    return 'Level 2', 'low'


def add_msi_confidence_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add MSI_levels and Metabograph_Confidence columns to a DataFrame."""
    if df is None:
        return df
    if not isinstance(df, pd.DataFrame):
        return df
    if df.empty:
        out = df.copy()
        out['MSI_levels'] = ''
        out['Metabograph_Confidence'] = ''
        return out

    out = df.copy()
    schema_cols = _resolve_confidence_columns(out)
    msi_conf = out.apply(lambda row: _compute_msi_level_and_confidence(row, schema_cols), axis=1)
    out['MSI_levels'] = msi_conf.apply(lambda x: x[0] if isinstance(x, tuple) and len(x) == 2 else '')
    out['Metabograph_Confidence'] = msi_conf.apply(lambda x: x[1] if isinstance(x, tuple) and len(x) == 2 else '')
    return out


def apply_confidence_filter(df: pd.DataFrame, confidence_filter_mode: str = 'exclude_low') -> pd.DataFrame:
    """Apply confidence-based filtering.

    Modes:
    - 'exclude_low': remove rows where Metabograph_Confidence == 'low'
    - 'none': keep all rows
    """
    if df is None or not isinstance(df, pd.DataFrame):
        return df
    if df.empty:
        return df.copy()

    mode = confidence_filter_mode if confidence_filter_mode in ('exclude_low', 'none') else 'exclude_low'
    out = df.copy()
    if mode == 'none':
        return out

    # Safety: if no MS2 evidence column exists, confidence cannot be meaningfully stratified.
    # In this case, keep rows unchanged instead of dropping everything as "low".
    schema_cols = _resolve_confidence_columns(out)
    if not schema_cols.get('ms2'):
        logger.warning(
            "Confidence filter skipped: MS2 column not found in current dataframe; "
            "keeping all rows for this sheet/output."
        )
        return out.reset_index(drop=True)

    if 'Metabograph_Confidence' not in out.columns:
        out = add_msi_confidence_columns(out)

    conf = out['Metabograph_Confidence'].astype(str).str.strip().str.lower()
    out = out[(conf != 'low') & (conf != '') & (conf != 'nan')].copy()
    return out.reset_index(drop=True)

class MetaboliteIDAnnotator:
    """
    Phase 1: Metabolite ID Annotation Only
    Priority: PubChem/KEGG APIs first, then offline HMDB databases
    Output: Metabolites with all available IDs (no pathway information)
    """

    def __init__(self, input_file: str = "metabolomics.xlsx", output_file: str = "metabolite_ids_annotated.xlsx", 
                 progress_callback=None,
                 stop_check_callback=None,
                 lipid_mode: bool = False,
                 skip_pubchem: bool = False,
                 skip_kegg: bool = False,
                 skip_hmdb: bool = False,
                 pubchem_enrichment: bool = True,
                 cleaned_metabolites_df: Optional[pd.DataFrame] = None,
                 input_df_override: Optional[pd.DataFrame] = None,
                 metabolite_ids_df: Optional[pd.DataFrame] = None,
                 force_file_input: bool = False,
                 # New: polarity DataFrames from memory + filtering options
                 pos_df: Optional[pd.DataFrame] = None,
                 neg_df: Optional[pd.DataFrame] = None,
                 id_filter_mode: str = 'none',
                 selected_id_columns: Optional[List[str]] = None,
                 skip_id_filtering: bool = False,
                 ms2_filter_mode: str = 'none',
                 confidence_filter_mode: str = 'exclude_low',
                 require_endogenous_yes: bool = False,
                 # New: configurable RT window (minutes) for Formula+ID dedup during ID annotation
                 dedup_rt_window_minutes: float = 2.0):
        self.input_file = input_file
        self.output_file = output_file
        self.progress_callback = progress_callback
        self.stop_check_callback = stop_check_callback if stop_check_callback else lambda: False
        # Lipid offline mode (uses lipid_search.py)
        self.lipid_mode = bool(lipid_mode)
        # hardcode PPM; not configurable in GUI/CLI
        self._lipid_ppm = 10.0
        
        # Database search options
        self.skip_pubchem = bool(skip_pubchem)
        self.skip_kegg = bool(skip_kegg)
        self.skip_hmdb = bool(skip_hmdb)
        
        # PubChem enrichment for cross-references
        # Force PubChem enrichment ON when in lipid_mode (always include pubchem crossrefs)
        if self.lipid_mode:
            self.pubchem_enrichment = True
        else:
            self.pubchem_enrichment = bool(pubchem_enrichment)
        
        # Store cleaned metabolites DataFrame from data cleaning process
        self.cleaned_metabolites_df = cleaned_metabolites_df

        # Optional direct dataframe override from the upload-time column assignment dialog
        self.input_df_override = input_df_override
        
        # Store metabolite IDs DataFrame for merging (from MzCloud/Metabolika files)
        self.metabolite_ids_df = metabolite_ids_df
        
        # Option to force file input over memory data (for custom user files)
        self.force_file_input = bool(force_file_input)

        # Polarity DataFrames (optional) for post-merge ID assignment
        self.pos_df = pos_df
        self.neg_df = neg_df

        # Filtering options for polarity DataFrames
        # id_filter_mode: 'none' | 'any_selected'
        self.id_filter_mode = id_filter_mode if id_filter_mode in ('none', 'any_selected') else 'none'
        self.selected_id_columns = selected_id_columns or []
        self.skip_id_filtering = bool(skip_id_filtering)
        self.ms2_filter_mode = ms2_filter_mode if ms2_filter_mode in ('none', 'preferred_only', 'preferred_or_other') else 'none'
        self.confidence_filter_mode = confidence_filter_mode if confidence_filter_mode in ('exclude_low', 'none') else 'exclude_low'
        self.require_endogenous_yes = bool(require_endogenous_yes)
        # RT window for Formula+ID summing during ID annotation stage
        try:
            self.id_dedup_rt_window_minutes = float(dedup_rt_window_minutes)
            if not (0 <= self.id_dedup_rt_window_minutes <= 20_000):
                self.id_dedup_rt_window_minutes = 2.0
        except Exception:
            self.id_dedup_rt_window_minutes = 2.0

        # Load offline databases for ID lookup
        logger.info("Loading offline databases for ID annotation...")
        
        # Import database finder utility
        from gui.shared.utils import find_database
        
        # Find HMDB database
        hmdb_path = find_database("hmdb_database.feather")
        
        if hmdb_path is None:
            error_msg = (
                "❌ Required database file 'hmdb_database.feather' not found!\n\n"
                "📖 Please refer to the 'Help' tab in the GUI for instructions on:\n"
                "   1. Downloading the HMDB database (hmdb_metabolites.xml)\n"
                "   2. Processing it into feather format\n"
                "   3. Placing the files in the correct directory\n\n"
                "💡 Quick Fix: Use the GUI's Help tab to automatically process your database files.\n"
                "📚 Or see DATABASE_SETUP_GUIDE.md for detailed instructions."
            )
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        try:
            self.hmdb_df = pd.read_feather(hmdb_path)
        except Exception as e:
            logger.error(f"Error loading HMDB database: {e}")
            raise
        
        # Load synonyms database for faster synonym lookups
        synonyms_path = find_database("All_metabolites_synonyms_hmdb.feather")
        
        if synonyms_path:
            try:
                self.hmdb_synonyms_df = pd.read_feather(synonyms_path)
                # Create lowercase synonym column for faster matching
                self.hmdb_synonyms_df['Synonyms_lower'] = self.hmdb_synonyms_df['Synonyms'].str.lower()
                logger.info(f"Loaded HMDB synonyms: {len(self.hmdb_synonyms_df)} entries")
            except Exception as e:
                logger.warning(f"⚠️ Failed to load synonyms file: {e}. Will fall back to main database for synonym matching.")
                self.hmdb_synonyms_df = None
        else:
            logger.warning("⚠️ Synonyms file 'All_metabolites_synonyms_hmdb.feather' not found. Will fall back to main database for synonym matching.")
            self.hmdb_synonyms_df = None
        
        # Load caches
        self.pubchem_cache = self._load_cache("pubchem_cache.pkl")
        self.kegg_cache = self._load_cache("kegg_cache.pkl")
        self.hmdb_cache = self._load_cache("hmdb_cache.pkl")
        self.lipidmaps_cache = self._load_cache("lipidmaps_cache.pkl")
        
        # Load MZCloud cache from data cleaner
        self.mzcloud_cache = self._load_mzcloud_cache()

        # Sanitize caches
        self._sanitize_caches()

        logger.info(f"Loaded HMDB: {len(self.hmdb_df)} entries")

        # Rate-limiting controls (thread-safe)
        self._pubchem_lock = threading.Lock()
        self._pubchem_min_interval = 1.0 / 5.0  # 5 requests per second
        self._pubchem_last_call = 0.0

        self._kegg_lock = threading.Lock()
        self._kegg_min_interval = 1.0 / 3.0  # 3 requests per second
        self._kegg_last_call = 0.0

    # -----------------------------
    # Lipid-mode (offline) helpers
    # -----------------------------
    def _build_lipidmaps_cache_key(self, data: Dict[str, Any]) -> str:
        """Build a comprehensive cache key for LipidMaps results using available IDs"""
        # Use multiple identifiers to create a robust cache key
        # CRITICAL: Include LipidID for lipid mode to prevent cache collisions when Name is empty
        identifiers = []
        for key in ['LipidID', 'LipidMaps_ID', 'InChIKey', 'SMILES', 'ChEBI_ID', 'PubChem_CID', 'Systematic_Name', 'Name']:
            val = data.get(key, '').strip()
            if val:
                identifiers.append(f"{key}:{val}")
        
        # If no identifiers found, use name as fallback
        if not identifiers:
            identifiers.append(f"Name:{data.get('Name', 'unknown')}")
        
        return f"lipidmaps_{'|'.join(identifiers)}"
    
    def _load_lipid_feather(self) -> pd.DataFrame:
        """Auto-detect an LMSD feather file in the working directory and load it."""
        from gui.shared.utils import find_database_paths
        
        candidates = [
            "lipidmap.feather", "lipidmaps.feather", "LMSD.feather", "lipidmaps_lmsd.feather", "lipidmaps_database.feather"
        ]
        
        # Check in all standard database locations
        for base_name in candidates:
            for path in find_database_paths():
                full_path = os.path.join(path, base_name)
                if os.path.exists(full_path):
                    try:
                        df = pd.read_feather(full_path)
                        logger.info(f"[LIPID MODE] Using lipid feather: {full_path}")
                        return df
                    except Exception as e:
                        logger.warning(f"[LIPID MODE] Failed reading {full_path}: {e}")
        
        # fallback: scan any .feather for LMSD-like columns in all search directories
        for search_dir in find_database_paths():
            if not os.path.exists(search_dir):
                continue
            try:
                for fn in os.listdir(search_dir):
                    if fn.lower().endswith('.feather'):
                        full_path = os.path.join(search_dir, fn)
                        try:
                            df = pd.read_feather(full_path)
                            cols = set(c.upper() for c in df.columns)
                            signature = {"LM_ID", "LMID", "LIPIDMAPS_ID", "CATEGORY", "MAIN_CLASS", "SUB_CLASS", "FORMULA", "EXACT_MASS"}
                            if cols & signature:
                                logger.info(f"[LIPID MODE] Auto-detected lipid feather: {full_path}")
                                return df
                        except Exception:
                            continue
            except Exception:
                pass
        
        error_msg = (
            "❌ Required database file 'lipidmap.feather' (or similar) not found!\n\n"
            "📖 Please refer to the 'Help' tab in the GUI for instructions on:\n"
            "   1. Downloading the LipidMaps database (structures.sdf)\n"
            "   2. Processing it into feather format\n"
            "   3. Placing the file in the correct directory\n\n"
            "💡 Quick Fix: Use the GUI's Help tab to automatically process your database files.\n"
            "📚 Or see DATABASE_SETUP_GUIDE.md for detailed instructions."
        )
        logger.error(error_msg)
        raise FileNotFoundError(error_msg)

    def _pick_first(self, d: Any, keys: List[str]) -> str:
        getter = getattr(d, 'get', None)
        for k in keys:
            try:
                v = getter(k) if getter else (d[k] if k in d else None)  # supports dict-like and pandas Series
            except Exception:
                v = None
            if v is not None and str(v).strip() != "":
                if pd.notna(v) and str(v).strip().lower() not in ("nan", "none", ""):
                    # Convert float to int for ID fields that should be integers
                    if k.upper() in ('PUBCHEM_CID', 'PUBCHEM CID', 'CHEBI_ID', 'CHEBI ID'):
                        try:
                            # If it's a float like 15511094.0, convert to int
                            float_val = float(v)
                            if float_val.is_integer():
                                return str(int(float_val))
                        except (ValueError, TypeError):
                            pass
                    return str(v).strip()
        return ""

    def _extract_lipid_id_fields(self, row: Any) -> Dict[str, Any]:
        """Map lipid_search output row to our standard annotation fields."""
        # Common variants across different LMSD exports
        lmid = self._pick_first(row, [
            'LM_ID', 'LMID', 'LM_ID ', 'LIPIDMAPS_ID', 'LIPIDMAPS ID', 'LM ID'
        ])
        # Prefer a stable LipidID originating from the input, then synthesized by lipid_search
        lipidid = self._pick_first(row, [
            'input_LipidID', 'LipidID',  # lipid_search inserts LipidID as first col
            'ABBREVIATION', 'input_ABBREVIATION',  # reasonable fallbacks
            'NAME', 'input_NAME'
        ])
    # Capture how lipid was matched in lipid_search (e.g., NAME_exact, EXACT_MASS_confirmed_by_NAME, etc.)
        match_type = self._pick_first(row, ['match_type', 'MATCH_TYPE'])
        # DO NOT extract KEGG_ID from LipidMaps - will be populated exclusively by KEGG API search
        # kegg = self._pick_first(row, ['KEGG_ID', 'KEGG ID', 'KEGG'])
        hmdb = self._pick_first(row, ['HMDBID', 'HMDB_ID', 'HMDB ID', 'HMDB'])
        chebi = self._pick_first(row, ['CHEBI_ID', 'ChEBI_ID', 'ChEBI ID', 'CHEBI'])
        pubchem = self._pick_first(row, ['PUBCHEM_CID', 'PubChem_CID', 'PUBCHEM CID', 'PubChem CID'])
        cas = self._pick_first(row, ['CAS', 'CAS_NUMBER', 'CAS Number'])
        inchikey = self._pick_first(row, ['INCHI_KEY', 'INCHI KEY', 'InChiKey', 'InChIKey'])
        inchi = self._pick_first(row, ['INCHI', 'InChI'])
        smiles = self._pick_first(row, ['SMILES'])
        formula = self._pick_first(row, ['FORMULA', 'Formula'])
        weight = self._pick_first(row, ['EXACT_MASS', 'Average Mass', 'AVERAGE_MASS'])
        category = self._pick_first(row, ['CATEGORY'])
        main_class = self._pick_first(row, ['MAIN_CLASS'])
        sub_class = self._pick_first(row, ['SUB_CLASS'])
        # Preferred names
        sys_name = self._pick_first(row, ['SYSTEMATIC_NAME', 'Systematic Name'])
        pref_name = self._pick_first(row, ['NAME', 'Preferred Name', 'PREFERRED_NAME'])
        abbr = self._pick_first(row, ['ABBREVIATION'])

        return {
            'LipidID': lipidid,
            'LipidMaps_ID': lmid,
            'LipidMaps_ID_Match_Type': match_type,
            'KEGG_ID': '',  # Will be populated exclusively by KEGG API search using Main_Class
            'PubChem_CID': pubchem,
            'HMDB_ID': hmdb,
            'ChEBI_ID': chebi,
            'CAS': cas,
            'InChIKey': inchikey,
            'InChI': inchi,
            'SMILES': smiles,
            'Molecular_Formula': formula,
            'Molecular_Weight': weight,
            'Super_Class': category,
            'Class': main_class,
            'Sub_Class': sub_class,
            'Systematic_Name': sys_name,
            'Preferred_Name': pref_name,
            'Abbreviation': abbr,
        }

    def run_lipid_id_annotation(self, max_workers: int = 4):
        """Run lipid ID annotation using offline lipid_search and enrich via PubChem, KEGG, and HMDB."""
        from main_script.lipid_search import run_search as lipid_run_search  # local import to avoid hard dep in metabolite mode

        # PRIORITY 1: Use explicit dataframe override when provided
        # PRIORITY 2: Use cleaned_metabolites_df from memory (unless force_file_input=True)
        # PRIORITY 3: Load from file
        # PRIORITY 4: Use pos_df only if nothing else available
        
        input_df = None
        polarity_sheets = {}
        
        if self.input_df_override is not None and isinstance(self.input_df_override, pd.DataFrame) and not self.input_df_override.empty:
            logger.info(f"[LIPID MODE] Using user-assigned input dataframe override")
            input_df = self.input_df_override.copy()
            logger.info(f"[LIPID MODE] Override dataframe: {len(input_df)} rows, columns: {list(input_df.columns)}")

        elif self.cleaned_metabolites_df is not None and isinstance(self.cleaned_metabolites_df, pd.DataFrame) and not self.cleaned_metabolites_df.empty and not self.force_file_input:
            logger.info(f"[LIPID MODE] Using pre-loaded COMBINED dataframe from memory (lipid_combined_df)")
            input_df = self.cleaned_metabolites_df.copy()
            logger.info(f"[LIPID MODE] COMBINED sheet: {len(input_df)} rows, columns: {list(input_df.columns)}")
            
            # Store polarity sheets from memory for later mapping
            if self.pos_df is not None and isinstance(self.pos_df, pd.DataFrame) and not self.pos_df.empty:
                polarity_sheets['Positive_Lipids'] = self.pos_df.copy()
                logger.info(f"[LIPID MODE] Positive sheet from memory: {len(self.pos_df)} rows")
            
            if self.neg_df is not None and isinstance(self.neg_df, pd.DataFrame) and not self.neg_df.empty:
                polarity_sheets['Negative_Lipids'] = self.neg_df.copy()
                logger.info(f"[LIPID MODE] Negative sheet from memory: {len(self.neg_df)} rows")
            
            if polarity_sheets:
                logger.info(f"[LIPID MODE] Will map annotated IDs back to {len(polarity_sheets)} polarity sheet(s)")
        
        elif self.input_file and os.path.exists(self.input_file):
            if self.force_file_input:
                logger.info(f"[LIPID MODE] Force file input enabled - bypassing pre-loaded memory dataframe")
            logger.info(f"[LIPID MODE] Loading input file: {self.input_file}")
            try:
                # Load Excel file and check sheets
                excel_file = pd.ExcelFile(self.input_file)
                sheet_names = excel_file.sheet_names
                logger.info(f"[LIPID MODE] Excel file contains {len(sheet_names)} sheet(s): {', '.join(sheet_names)}")
                
                # Try to find the best sheet for annotation (prefer sheets with LipidID, Class, Main_Class)
                annotation_sheet = None
                
                # Check for known combined/annotation sheet names
                preferred_sheets = ['Combined', 'Lipids', 'Annotation', 'All_Lipids']
                for pref in preferred_sheets:
                    if pref in sheet_names:
                        annotation_sheet = pref
                        logger.info(f"[LIPID MODE] Found preferred annotation sheet: '{annotation_sheet}'")
                        break
                
                # If no preferred sheet, use first sheet
                if annotation_sheet is None:
                    annotation_sheet = sheet_names[0]
                    logger.info(f"[LIPID MODE] Using first sheet for annotation: '{annotation_sheet}'")
                
                # Load annotation sheet
                input_df = pd.read_excel(self.input_file, sheet_name=annotation_sheet)
                logger.info(f"[LIPID MODE] Loaded sheet '{annotation_sheet}': {len(input_df)} rows, columns: {list(input_df.columns)}")
                
                # Check if this sheet has the minimal required columns (LipidID, Class)
                required_cols = {'LipidID', 'Class'}
                available_cols = set(input_df.columns)
                
                if not required_cols.issubset(available_cols):
                    logger.warning(f"[LIPID MODE] Sheet '{annotation_sheet}' missing required columns. Looking for: {required_cols}, Found: {available_cols}")
                    # Try to find a sheet with required columns
                    for sheet in sheet_names:
                        test_df = pd.read_excel(self.input_file, sheet_name=sheet, nrows=0)
                        if required_cols.issubset(set(test_df.columns)):
                            logger.info(f"[LIPID MODE] Found sheet with required columns: '{sheet}'")
                            input_df = pd.read_excel(self.input_file, sheet_name=sheet)
                            annotation_sheet = sheet
                            break
                
                # Store polarity sheets if they exist (for mapping results back)
                polarity_sheet_names = {'Positive_Lipids', 'Negative_Lipids', 'Pos_Lipids', 'Neg_Lipids'}
                class_sheets = {}  # Store class sheets separately (they don't need ID annotation)
                
                for sheet in sheet_names:
                    if sheet in polarity_sheet_names:
                        try:
                            polarity_sheets[sheet] = pd.read_excel(self.input_file, sheet_name=sheet)
                            logger.info(f"[LIPID MODE] Found polarity sheet: '{sheet}' with {len(polarity_sheets[sheet])} rows")
                        except Exception as e:
                            logger.warning(f"[LIPID MODE] Could not load polarity sheet '{sheet}': {e}")
                    # Detect class sheets (they start with "CLASS_" or contain "class" in name)
                    elif sheet.startswith('CLASS_') or '_class' in sheet.lower() or sheet.lower().startswith('class'):
                        try:
                            class_sheets[sheet] = pd.read_excel(self.input_file, sheet_name=sheet)
                            logger.info(f"[LIPID MODE] Found class sheet: '{sheet}' with {len(class_sheets[sheet])} rows (will be preserved as-is)")
                        except Exception as e:
                            logger.warning(f"[LIPID MODE] Could not load class sheet '{sheet}': {e}")
                
                if polarity_sheets:
                    logger.info(f"[LIPID MODE] Will map annotated IDs back to {len(polarity_sheets)} polarity sheet(s)")
                else:
                    logger.info(f"[LIPID MODE] No polarity sheets found - will output only annotated results")
                
                if class_sheets:
                    logger.info(f"[LIPID MODE] Found {len(class_sheets)} class sheet(s) - will preserve without ID annotation")
                    
            except Exception as e:
                logger.error(f"[LIPID MODE] Error loading file: {e}")
                raise
        
        elif self.pos_df is not None and isinstance(self.pos_df, pd.DataFrame) and not self.pos_df.empty:
            logger.info(f"[LIPID MODE] Fallback: Using pos_df as input (no combined sheet or file available)")
            input_df = self.pos_df.copy()
            # No polarity sheets in this fallback mode
        
        else:
            raise ValueError("[LIPID MODE] No input data available: cleaned_metabolites_df, input_file, and pos_df are all unavailable")
        
        # Validate input dataframe has required columns
        if input_df is None or input_df.empty:
            raise ValueError("[LIPID MODE] Input dataframe is empty")
        
        required_cols = ['LipidID', 'Class']
        available_cols = set(input_df.columns)
        if not set(required_cols).issubset(available_cols):
            raise ValueError(f"[LIPID MODE] Input dataframe missing required columns. Need: {required_cols}, Have: {list(available_cols)}")
        
        required_cols = ['LipidID', 'Class']
        available_cols = set(input_df.columns)
        if not set(required_cols).issubset(available_cols):
            raise ValueError(f"[LIPID MODE] Input dataframe missing required columns. Need: {required_cols}, Have: {list(available_cols)}")
        
        # Report database search options
        logger.info(f"[LIPID MODE] Enabled databases: LipidMaps + PubChem cross-references")
        # Ensure PubChem cross-reference enrichment is always enabled in lipid mode
        if not self.pubchem_enrichment:
            logger.info(f"[LIPID MODE] Forcing PubChem cross-reference enrichment: ENABLED (overriding GUI setting)")
        else:
            logger.info(f"[LIPID MODE] PubChem cross-reference enrichment: ENABLED")
        self.pubchem_enrichment = True

        # Identify name columns available (case-insensitive matching)
        name_columns = []
        expected_name_cols = ['Name', 'Metabolite', 'Molecule', 'ABBREVIATION', 'SYSTEMATIC_NAME', 'Abbreviation']
        input_cols_lower = {col.lower(): col for col in input_df.columns}
        
        for expected in expected_name_cols:
            expected_lower = expected.lower()
            if expected_lower in input_cols_lower:
                actual_col = input_cols_lower[expected_lower]
                if actual_col not in name_columns:  # avoid duplicates
                    name_columns.append(actual_col)
        logger.info(f"[LIPID MODE] Name columns available for search: {name_columns}")

        # Count non-empty metabolites in name columns (for reporting only)
        for col in name_columns:
            if col in input_df.columns:
                non_empty = input_df[col].notna().sum()
                logger.info(f"[LIPID MODE] Column '{col}': {non_empty} non-empty entries")

        # Use actual row count as total for progress tracking (not max non-empty values)
        valid_count = len(input_df)
        logger.info(f"[LIPID MODE] {valid_count} total lipids (rows) to process")

        # Notify GUI about the analysis: report as (current=0, total=valid_count)
        if self.progress_callback:
            try:
                total_for_progress = max(valid_count, 1)
                self.progress_callback(0, total_for_progress, f"📊 Found {valid_count} metabolites in columns: {', '.join(name_columns)}")
            except Exception:
                pass

        # Load lipid feather
        df_feather = self._load_lipid_feather()

        # Run matching
        logger.info(f"[LIPID MODE] Running offline lipid search")

        # Build enrichment message based on enabled databases
        enrichment_dbs = []
        if not self.skip_pubchem and not self.lipid_mode:
            enrichment_dbs.append("PubChem")
        if not self.skip_kegg and not self.lipid_mode:
            enrichment_dbs.append("KEGG")
        if not self.skip_hmdb and not self.lipid_mode:
            enrichment_dbs.append("HMDB")
        
        # Add PubChem cross-reference info
        cross_ref_msg = ""
        if self.pubchem_enrichment:
            cross_ref_msg = " + PubChem cross-refs (KEGG/HMDB/ChEBI)"
        
        if self.lipid_mode:
            enrichment_msg = f"Starting LipidMaps search{cross_ref_msg}..."
            logger.info(f"[LIPID MODE] {enrichment_msg}")
        elif enrichment_dbs:
            enrichment_msg = f"Starting {', '.join(enrichment_dbs)} enrichment{cross_ref_msg} for matched lipids..."
            logger.info(f"[LIPID MODE] {enrichment_msg}")
        else:
            enrichment_msg = f"LipidMaps-only mode (all enrichment databases skipped){cross_ref_msg}"
            logger.info(f"[LIPID MODE] {enrichment_msg}")

        # Notify GUI about enrichment phase: keep using valid_count as the total so percent is meaningful
        if self.progress_callback:
            try:
                total_for_progress = max(valid_count, 1)
                if self.lipid_mode:
                    progress_msg = f"🔍 Phase 1/3: Offline lipid matching{cross_ref_msg}..."
                    self.progress_callback(0, total_for_progress, progress_msg)
                elif enrichment_dbs or self.pubchem_enrichment:
                    progress_msg = f"🔍 Phase 1/3: Offline lipid matching and {'/'.join(enrichment_dbs)} enrichment"
                    if self.pubchem_enrichment:
                        progress_msg += " + PubChem cross-refs"
                    progress_msg += "..."
                    self.progress_callback(0, total_for_progress, progress_msg)
                else:
                    self.progress_callback(0, total_for_progress, "🔍 Phase 1/3: Offline lipid matching (enrichment databases disabled)...")
            except Exception:
                pass

        # Create wrapper callback for Phase 1: Map 0-100% of lipid_run_search to 0-40% overall
        def phase1_progress_callback(current, total, msg):
            if self.progress_callback and total > 0:
                # Map Phase 1 progress to overall 0-40% range
                # Calculate the actual position within 0-40% based on current/total
                phase1_progress = current / total  # 0.0 to 1.0
                overall_current = int(phase1_progress * 40)  # 0 to 40
                overall_total = 100  # Total is always 100 for unified progress
                
                phase_msg = f"Phase 1/3 (Lipid Search): {msg}"
                
                # Call GUI with current count mapped to 0-40 range, total=100
                self.progress_callback(overall_current, overall_total, phase_msg)
        
        # Pass through wrapped progress_callback for lipid matching
        df_all_results = lipid_run_search(input_df, df_feather, ppm=self._lipid_ppm, progress_callback=phase1_progress_callback, max_workers=max_workers)

        # Separate matched and unmatched rows
        df_matches = pd.DataFrame()
        df_unmatched = pd.DataFrame()
        
        try:
            if 'matched' in df_all_results.columns:
                df_matches = df_all_results[df_all_results['matched'] == True].reset_index(drop=True)
                df_unmatched = df_all_results[df_all_results['matched'] == False].reset_index(drop=True)
                
                lipidmaps_matched = len(df_matches)
                lipidmaps_unmatched = len(df_unmatched)
                
                logger.info(f"[LIPID MODE] LipidMaps matching complete:")
                logger.info(f"  ✅ Matched: {lipidmaps_matched}")
                logger.info(f"  ❌ Unmatched: {lipidmaps_unmatched}")
            else:
                df_matches = df_all_results
                logger.warning("[LIPID MODE] 'matched' column not found; treating all rows as matched")
        except Exception as e:
            logger.warning(f"[LIPID MODE] Error separating matched/unmatched rows: {e}")
            df_matches = df_all_results

        # KEGG search for ALL lipids using Main_Class to populate KEGG_ID column
        # (Do NOT use LipidMaps or PubChem for KEGG_ID - only KEGG API by Main_Class)
        class_col = 'Main_Class' if 'Main_Class' in input_df.columns else None
        if class_col:
            logger.info(f"[LIPID MODE] Starting KEGG search for ALL {len(input_df)} lipids using {class_col} to populate KEGG_ID...")
            
            if self.progress_callback:
                try:
                    self.progress_callback(40, 100, f"🔍 Phase 2/3 (KEGG Search): searching ALL lipids by {class_col}...")
                except Exception:
                    pass
            
            # Build a mapping of LipidID -> Main_Class from input_df for faster lookup
            lipid_class_map = {}
            if 'LipidID' in input_df.columns and class_col in input_df.columns:
                for _, input_row in input_df.iterrows():
                    lipid_id = input_row.get('LipidID')
                    class_name = input_row.get(class_col)
                    if pd.notna(lipid_id) and pd.notna(class_name):
                        lipid_class_map[str(lipid_id).strip()] = str(class_name).strip()
                
                logger.info(f"[KEGG SEARCH] ✅ Built {class_col} lookup map with {len(lipid_class_map)} entries")
                
                # Log sample entries for debugging
                if lipid_class_map:
                    sample_entries = list(lipid_class_map.items())[:5]
                    logger.info(f"[KEGG SEARCH] 📋 Sample class map (first 5): {sample_entries}")
                    
                    # Show unique class values
                    unique_classes = set(lipid_class_map.values())
                    logger.info(f"[KEGG SEARCH] 🏷️  Unique Main_Class values found: {sorted(unique_classes)}")
            else:
                logger.error(f"[KEGG SEARCH] ❌ Cannot build class map!")
                logger.error(f"  - LipidID in columns: {'LipidID' in input_df.columns}")
                logger.error(f"  - {class_col} in columns: {class_col in input_df.columns if class_col else 'N/A'}")
                logger.error(f"  - Available columns: {list(input_df.columns)}")
                lipid_class_map = {}
            
            # Build KEGG results for ALL input lipids (not just unmatched)
            # Store KEGG metadata so we can also report how the ID was matched later
            kegg_id_map = {}  # LipidID -> {'KEGG_ID': str, 'match_type': 'Class'}
            processed_count = 0
            kegg_class_cache = {}  # Cache results by class taxonomy to avoid redundant API calls
            
            # Iterate through unique LipidIDs from input to search KEGG
            unique_lipids = input_df['LipidID'].dropna().unique() if 'LipidID' in input_df.columns else []
            
            for lipid_id in unique_lipids:
                try:
                    lipid_id_str = str(lipid_id).strip()
                    # Get class taxonomy from the lookup map
                    class_name = lipid_class_map.get(lipid_id_str)
                    
                    if class_name:
                        # Normalize class_name: remove trailing 's' for KEGG search (singular form)
                        # e.g., "Ceramides" -> "Ceramide", "Ceramide" -> "Ceramide"
                        normalized_class_name = class_name.rstrip('s') if class_name.endswith('s') else class_name
                        
                        # Check cache first (use normalized name for cache key)
                        cache_key = normalized_class_name
                        if cache_key in kegg_class_cache:
                            kegg_result = kegg_class_cache[cache_key]
                            logger.debug(f"[KEGG CACHE HIT] LipidID: {lipid_id_str}, Class: '{class_name}' -> Normalized: '{normalized_class_name}'")
                        else:
                            logger.info(f"[KEGG SEARCH] 🔍 Searching for LipidID: {lipid_id_str} | Class: '{class_name}' -> Normalized: '{normalized_class_name}'")
                            kegg_result = self.search_kegg_api(normalized_class_name)
                            logger.debug(f"[KEGG RESULT] {normalized_class_name}: {kegg_result}")
                            kegg_class_cache[cache_key] = kegg_result
                        
                        if kegg_result and kegg_result.get('found') and kegg_result.get('KEGG_ID'):
                            kegg_id_map[lipid_id_str] = {
                                'KEGG_ID': kegg_result['KEGG_ID'],
                                'match_type': class_col or 'Main_Class'
                            }
                            logger.info(f"[KEGG SEARCH] ✅ Found KEGG_ID {kegg_result['KEGG_ID']} for {class_col}: '{normalized_class_name}' (original: '{class_name}', LipidID: {lipid_id_str})")
                        else:
                            logger.debug(f"[KEGG SEARCH] ❌ No KEGG match for {class_col}: '{normalized_class_name}' (original: '{class_name}')")
                    else:
                        logger.debug(f"[KEGG SEARCH] ⚠️ No {class_col} found for LipidID: {lipid_id_str}")
                    
                    processed_count += 1
                    # Progress update every 10 items or at the end
                    if self.progress_callback and (processed_count % 10 == 0 or processed_count == len(unique_lipids)):
                        try:
                            # Map Phase 2 progress (0-100%) to overall 40-70% (30% range)
                            phase2_percent = (processed_count / len(unique_lipids)) * 30  # Phase 2 is 30% of total
                            overall_percent = int(40 + phase2_percent)  # Start at 40%
                            self.progress_callback(overall_percent, 100, f"Phase 2/3 (KEGG): {processed_count}/{len(unique_lipids)}")
                        except Exception:
                            pass
                except Exception as e:
                    logger.warning(f"[KEGG SEARCH] Error processing LipidID {lipid_id}: {e}")
                    continue
            
            # Report unique classes searched
            unique_classes_searched = set(kegg_class_cache.keys())
            classes_with_hits = {k: v for k, v in kegg_class_cache.items() if v.get('found')}
            
            logger.info(f"[KEGG SEARCH] Summary:")
            logger.info(f"  Total unique lipid IDs: {len(unique_lipids)}")
            logger.info(f"  Unique class names searched: {len(unique_classes_searched)}")
            logger.info(f"  Classes with KEGG hits: {len(classes_with_hits)}")
            logger.info(f"  Lipids assigned KEGG IDs: {len(kegg_id_map)}")
            
            if classes_with_hits:
                logger.info(f"  Classes with KEGG matches: {', '.join(list(classes_with_hits.keys())[:10])}")
                if len(classes_with_hits) > 10:
                    logger.info(f"  ... and {len(classes_with_hits) - 10} more classes")
            else:
                logger.warning(f"[KEGG SEARCH] ⚠️ NO CLASSES RETURNED KEGG HITS! Searched classes: {', '.join(list(unique_classes_searched)[:10])}")
            
            # Store kegg_id_map for later population (after enrichment to avoid being overwritten)
            # We'll populate KEGG_ID in annotated_df after the enrichment loop completes
        else:
            kegg_id_map = {}  # Empty map if no class taxonomy column
            logger.info("[LIPID MODE] Main_Class column not found; skipping KEGG search")

        # NOTE: Backup file creation removed per user request. We proceed
        # directly to enrichment without writing intermediate backups.

        # Build standard annotation output in parallel
        total_rows = len(df_matches)
        out_rows: List[Optional[Dict[str, Any]]] = [None] * total_rows

        def enrich_row(row: Any) -> Dict[str, Any]:
            base: Dict[str, Any] = {}
            # Preserve the KEGG fallback marker if present
            is_kegg_fallback = False
            if isinstance(row, (dict,)):
                marker_val = row.get('_kegg_fallback_match', False)
                is_kegg_fallback = marker_val if marker_val is not None and marker_val is not False else False
            elif hasattr(row, '_kegg_fallback_match'):
                marker_val = getattr(row, '_kegg_fallback_match', False)
                # Convert pandas NaN or any non-True value to False
                import math
                if marker_val is True or (isinstance(marker_val, (bool,)) and marker_val):
                    is_kegg_fallback = True
                elif isinstance(marker_val, float) and math.isnan(marker_val):
                    is_kegg_fallback = False
                else:
                    is_kegg_fallback = bool(marker_val) if marker_val is not None else False
            
            # Choose a display name from input columns
            name_val = self._pick_first(row, [
                'input_ABBREVIATION', 'input_NAME', 'input_SYSTEMATIC_NAME', 'input_Name', 'input_Metabolite', 'input_Molecule'
            ])
            base['Name'] = name_val or ''

            mapped = self._extract_lipid_id_fields(row)
            base.update(mapped)
            base.setdefault('KEGG_Match_Type', '')
            
            # Restore KEGG fallback marker
            if is_kegg_fallback:
                base['_kegg_fallback_match'] = True

            # Determine which input column produced the match (from match_type)
            match_type = (row.get('match_type') or '') if isinstance(row, (dict,)) else (getattr(row, 'match_type', '') or '')
            match_source = ''
            try:
                if match_type:
                    # match_type patterns like 'ABBREVIATION_exact' or 'NAME_exact_confirmed_by_EXACT_MASS'
                    primary = str(match_type).split('_')[0]
                    if 'fallback_from_ABBREVIATION' in str(match_type):
                        match_source = f"{primary} (fallback from ABBREVIATION)"
                    else:
                        match_source = primary
            except Exception:
                match_source = ''
            base['match_source'] = match_source
            if match_source:
                logger.debug(f"[LIPID MODE] Match source for '{base.get('Name','')}' -> {match_source}")

            # Check LipidMaps cache first (comprehensive ID-based caching)
            lipidmaps_cache_key = self._build_lipidmaps_cache_key(base)
            if lipidmaps_cache_key in self.lipidmaps_cache:
                cached_result = self.lipidmaps_cache[lipidmaps_cache_key]
                logger.debug(f"[CACHE] LipidMaps cache hit for: {base['Name']}")
                base.update(cached_result)
                base['annotation_sources'] = ['LipidMaps_Cache']
                # Ensure marker is preserved correctly (don't let cache overwrite it)
                if is_kegg_fallback:
                    base['_kegg_fallback_match'] = True
                elif '_kegg_fallback_match' in base:
                    # This row is NOT a KEGG fallback, so remove the marker if it came from cache
                    del base['_kegg_fallback_match']
                return base

            # PubChem enrichment - For lipid mode, we skip individual PubChem searches and only do cross-references
            pubchem_sources = []
            if self.lipid_mode:
                # In lipid mode, we don't do individual PubChem searches, only cross-references later
                logger.debug(f"[LIPID MODE] Skipping individual PubChem search for: {base['Name']}")
            elif not self.skip_pubchem and not base.get('PubChem_CID'):
                # Use high-quality identifiers first: SYSTEMATIC_NAME, InChIKey, InChI, SMILES, ChEBI_ID
                query_identifiers = [
                    ('systematic_name', base.get('Systematic_Name')),
                    ('inchikey', base.get('InChIKey')),
                    ('inchi', base.get('InChI')),
                    ('smiles', base.get('SMILES')),
                    ('chebi_id', base.get('ChEBI_ID'))
                ]
                
                for id_type, identifier in query_identifiers:
                    if identifier and identifier.strip():
                        expected_formula = base.get('Molecular_Formula')
                        if id_type == 'inchikey':
                            pubchem_res = self.search_pubchem_by_inchikey(identifier.strip())
                        elif id_type == 'inchi':
                            pubchem_res = self.search_pubchem_by_inchi(identifier.strip())
                        elif id_type == 'smiles':
                            pubchem_res = self.search_pubchem_by_smiles(identifier.strip())
                        elif id_type == 'chebi_id':
                            pubchem_res = self.search_pubchem_by_chebi(identifier.strip())
                        else:  # systematic_name
                            pubchem_res = self.search_pubchem_api(identifier.strip(), expected_formula)
                        
                        if pubchem_res.get('found'):
                            base['PubChem_CID'] = pubchem_res.get('PubChem_CID', '')
                            base['SMILES'] = base.get('SMILES') or pubchem_res.get('SMILES', '')
                            base['InChI'] = base.get('InChI') or pubchem_res.get('InChI', '')
                            base['InChIKey'] = base.get('InChIKey') or pubchem_res.get('InChIKey', '')
                            base['IUPAC_Name'] = base.get('IUPAC_Name') or pubchem_res.get('IUPAC_Name', '')
                            base['Molecular_Formula'] = base.get('Molecular_Formula') or pubchem_res.get('Molecular_Formula', '')
                            base['Molecular_Weight'] = base.get('Molecular_Weight') or pubchem_res.get('Molecular_Weight', '')
                            # Cross-refs
                            base['ChEBI_ID'] = base.get('ChEBI_ID') or pubchem_res.get('ChEBI_ID', '')
                            base['KEGG_ID'] = base.get('KEGG_ID') or pubchem_res.get('KEGG_ID', '')
                            base['HMDB_ID'] = base.get('HMDB_ID') or pubchem_res.get('HMDB_ID', '')
                            base['CAS'] = base.get('CAS') or pubchem_res.get('CAS', '')
                            pubchem_sources.append('PubChem_API')
                            logger.info(f"[LIPID MODE] PubChem match found via {id_type} for: {base['Name']}")
                            break
            elif self.skip_pubchem:
                logger.debug(f"[SKIP] PubChem search skipped for: {base['Name']}")

            # KEGG enrichment (only if still missing after PubChem AND if skip_kegg is not enabled)
            if not self.skip_kegg and not base.get('KEGG_ID') and not self.lipid_mode:
                # Use systematic name first, then other identifiers
                kegg_identifiers = [
                    base.get('Systematic_Name'),
                    base.get('InChIKey'),
                    base.get('ChEBI_ID'),
                    base.get('Preferred_Name'),
                    base['Name']
                ]
                for query_name in kegg_identifiers:
                    if query_name and query_name.strip():
                        kegg_res = self.search_kegg_api(query_name.strip())
                        if kegg_res.get('found'):
                            base['KEGG_ID'] = kegg_res.get('KEGG_ID', base.get('KEGG_ID', ''))
                            base['KEGG_Match_Type'] = 'Name'
                            if not base.get('Molecular_Formula') and kegg_res.get('Molecular_Formula'):
                                base['Molecular_Formula'] = kegg_res['Molecular_Formula']
                            break
            elif self.skip_kegg or self.lipid_mode:
                logger.debug(f"[SKIP] KEGG search skipped for: {base['Name']}")

            # HMDB enrichment using enhanced ID-based context (only if skip_hmdb is not enabled)
            if not self.skip_hmdb and not self.lipid_mode:
                existing_ids_ctx = {
                    'HMDB_ID': base.get('HMDB_ID', ''),
                    'InChIKey': base.get('InChIKey', ''),
                    'Molecular_Formula': base.get('Molecular_Formula', ''),
                    'PubChem_CID': base.get('PubChem_CID', ''),
                    'KEGG_ID': base.get('KEGG_ID', ''),
                    'ChEBI_ID': base.get('ChEBI_ID', ''),
                    'SMILES': base.get('SMILES', ''),
                    'InChI': base.get('InChI', ''),
                }
                query_name_hmdb = base.get('Systematic_Name') or base.get('Preferred_Name') or base['Name']
                hmdb_res = self.search_hmdb_offline(query_name_hmdb or '', existing_ids_ctx)
                if hmdb_res.get('found'):
                    base['HMDB_ID'] = base.get('HMDB_ID') or hmdb_res.get('HMDB_ID', '')
                    base['InChIKey'] = base.get('InChIKey') or hmdb_res.get('InChIKey', '')
                    base['Molecular_Formula'] = base.get('Molecular_Formula') or hmdb_res.get('Molecular_Formula', '')
                    base['Endogenous'] = _normalize_yes_no(hmdb_res.get('Endogenous', 'No'))
                    base['Super_Class'] = base.get('Super_Class') or hmdb_res.get('Super_Class', '')
                    base['Class'] = base.get('Class') or hmdb_res.get('Class', '')
                    base['Sub_Class'] = base.get('Sub_Class') or hmdb_res.get('Sub_Class', '')
            elif self.skip_hmdb or self.lipid_mode:
                logger.debug(f"[SKIP] HMDB search skipped for: {base['Name']}")

            # PubChem cross-reference enrichment for missing IDs (only if enabled)
            if self.pubchem_enrichment and base.get('PubChem_CID'):
                try:
                    pubchem_cid = int(float(str(base['PubChem_CID']).strip()))
                    
                    # Only fetch cross-refs for missing IDs to avoid overwriting existing data
                    # In lipid mode, do NOT populate KEGG_ID from PubChem (only from KEGG API by Main_Class)
                    needs_enrichment = (
                        not base.get('HMDB_ID') or 
                        not base.get('ChEBI_ID') or 
                        not base.get('CAS') or
                        not base.get('LipidMaps_ID')
                    )
                    
                    if needs_enrichment:
                        logger.debug(f"[PUBCHEM_XREF] Fetching cross-references for CID {pubchem_cid}: {base['Name']}")
                        xrefs = fetch_pubchem_cross_references(pubchem_cid)
                        
                        # DO NOT populate KEGG_ID from PubChem in lipid mode - only from KEGG API search
                        # if (not base.get('KEGG_ID') or str(base.get('KEGG_ID')).strip() == '') and xrefs.get('KEGG'):
                        #     base['KEGG_ID'] = xrefs['KEGG']
                        if not base.get('HMDB_ID') and xrefs.get('HMDB'):
                            base['HMDB_ID'] = xrefs['HMDB']
                        if not base.get('ChEBI_ID') and xrefs.get('ChEBI'):
                            base['ChEBI_ID'] = xrefs['ChEBI']
                        if not base.get('CAS') and xrefs.get('CAS'):
                            base['CAS'] = xrefs['CAS']
                        if not base.get('LipidMaps_ID') and xrefs.get('LipidMaps'):
                            base['LipidMaps_ID'] = xrefs['LipidMaps']
                            
                        # Log successful enrichments
                        enriched_fields = []
                        for field, value in [('HMDB_ID', xrefs.get('HMDB')), 
                                           ('ChEBI_ID', xrefs.get('ChEBI')), ('CAS', xrefs.get('CAS')), 
                                           ('LipidMaps_ID', xrefs.get('LipidMaps'))]:
                            if value:
                                enriched_fields.append(f"{field}={value}")
                        if enriched_fields:
                            logger.info(f"[PUBCHEM_XREF] Enriched {base['Name']}: {', '.join(enriched_fields)}")
                            
                except (ValueError, TypeError) as e:
                    logger.warning(f"[PUBCHEM_XREF] Invalid PubChem_CID '{base.get('PubChem_CID')}' for {base['Name']}: {e}")
                except Exception as e:
                    logger.error(f"[PUBCHEM_XREF] Error fetching cross-references for {base['Name']}: {e}")

            if not base.get('IUPAC_Name'):
                base['IUPAC_Name'] = base.get('Systematic_Name', '')

            # Check if row already has annotation_sources from KEGG fallback
            existing_sources = []
            if isinstance(row, (dict,)) and 'annotation_sources' in row:
                existing_sources_val = row.get('annotation_sources')
                if isinstance(existing_sources_val, str):
                    existing_sources = [existing_sources_val]
                elif isinstance(existing_sources_val, list):
                    existing_sources = existing_sources_val
            elif hasattr(row, 'annotation_sources'):
                existing_sources_val = getattr(row, 'annotation_sources', None)
                if isinstance(existing_sources_val, str):
                    existing_sources = [existing_sources_val]
                elif isinstance(existing_sources_val, list):
                    existing_sources = existing_sources_val
            
            # Build sources list
            if 'KEGG_API_ClassNameFallback' in existing_sources:
                # This row was matched via KEGG fallback, not LipidMaps
                sources = ['KEGG_API_ClassNameFallback']
            else:
                sources = ['LipidMaps_Feather']
            
            if pubchem_sources:
                sources.extend(pubchem_sources)

            # Only report KEGG_API if we actually obtained a KEGG ID and KEGG wasn't skipped
            if not self.skip_kegg and base.get('KEGG_ID') and not self.lipid_mode:
                if 'KEGG_API_ClassNameFallback' not in sources:  # Don't duplicate
                    sources.append('KEGG_API')

            # Only report HMDB_Offline if we actually obtained an HMDB ID and HMDB wasn't skipped
            if not self.skip_hmdb and base.get('HMDB_ID') and not self.lipid_mode:
                sources.append('HMDB_Offline')

            # Track PubChem cross-reference enrichment only when cross-refs were actually retrieved
            if self.pubchem_enrichment and base.get('PubChem_CID'):
                # Determine whether cross-ref fields were filled by PubChem fetch
                xref_fields = ('KEGG_ID', 'HMDB_ID', 'ChEBI_ID', 'CAS', 'LipidMaps_ID')
                has_xref = any(base.get(f) for f in xref_fields)
                if pubchem_sources or has_xref:
                    sources.append('PubChem_CrossRef')

            # Do not add skip markers for lipid mode; only report actual data sources.
                
            base['annotation_sources'] = sources

            # Cache the enriched result (EXCLUDE _kegg_fallback_match marker to prevent cache pollution)
            cache_data = base.copy()
            cache_data.pop('_kegg_fallback_match', None)  # Remove marker from cache
            self.lipidmaps_cache[lipidmaps_cache_key] = cache_data
            return base

        # Submit tasks
        tasks = []
        rows_list = list(df_matches.iterrows())
        for idx, (_, row) in enumerate(rows_list):
            tasks.append((idx, row))

        processed = 0
        last_logged = 0
        lock = threading.Lock()
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {executor.submit(enrich_row, row): idx for idx, row in tasks}
            for future in as_completed(future_map):
                # Check if user requested stop
                if self.stop_check_callback():
                    logger.info("[LIPID MODE] Stop requested by user - cancelling remaining tasks")
                    # Cancel pending futures
                    for fut in future_map.keys():
                        fut.cancel()
                    break
                
                idx = future_map[future]
                try:
                    base = future.result()
                    out_rows[idx] = base
                except Exception as e:
                    logger.error(f"[LIPID MODE] Worker error at row {idx}: {e}")
                    out_rows[idx] = {'Name': '', 'annotation_sources': ['ERROR']}
                # Progress update
                with lock:
                    processed += 1
                    if self.progress_callback:
                        try:
                            name_disp = out_rows[idx].get('Name', '') if out_rows[idx] else ''
                            # Map Phase 3 progress (0-100%) to overall 70-100% (30% range)
                            phase3_percent = (processed / total_rows) * 30  # Phase 3 is 30% of total
                            overall_percent = int(70 + phase3_percent)  # Start at 70%
                            self.progress_callback(overall_percent, 100, f"Phase 3/3 (Enrichment): {name_disp}")
                        except Exception:
                            pass
                    # Log every 5% or every 100 rows
                    percent = (processed / total_rows * 100) if total_rows else 0
                    if percent - last_logged >= 5 or processed % 100 == 0 or processed == total_rows:
                        logger.info(f"[LIPID MODE] Progress: {processed}/{total_rows} ({percent:.1f}%)")
                        last_logged = percent

        # Convert results to DataFrame (preserve order)
        annotated_df = pd.DataFrame([r if r is not None else {} for r in out_rows])
        
        # CRITICAL: Dedup annotated_df BEFORE merging to ensure best match per LipidID survives
        # This prevents merge from losing KEGG fallback IDs when there are multiple matches per LipidID
        if not annotated_df.empty and 'LipidID' in annotated_df.columns:
            # Apply same dedup logic as later: prefer rows with more IDs, preserve KEGG fallback marker
            dedup_basis = 'LipidID'
            
            def _norm_key_inline(x):
                import re as _re
                if pd.isna(x):
                    return ''
                s = str(x).strip()
                s = _re.sub(r"\+.*$", "", s)  # remove adducts
                s = s.replace("_", "/")
                if '(' not in s and ')' not in s:
                    m = _re.match(r"^([A-Za-z0-9+-]+)\s+(.+)$", s)
                    if m:
                        rest = m.group(2).strip()
                        if ':' in rest or '/' in rest:
                            s = f"{m.group(1)}({rest})"
                s = _re.sub(r"\s+", "", s)
                return s
            
            annotated_df['__dedup_key'] = annotated_df[dedup_basis].apply(_norm_key_inline)
            
            # Separate KEGG fallback from non-KEGG
            if '_kegg_fallback_match' in annotated_df.columns:
                is_kegg = annotated_df['_kegg_fallback_match'].fillna(False).astype(bool)
            else:
                is_kegg = pd.Series([False] * len(annotated_df), index=annotated_df.index)
            
            kegg_rows = annotated_df[is_kegg].copy() if is_kegg.any() else pd.DataFrame(columns=annotated_df.columns)
            non_kegg = annotated_df[~is_kegg].copy()
            
            # CRITICAL FIX: Merge KEGG fallback IDs into LipidMaps rows when both exist for same LipidID
            # This prevents losing KEGG_IDs when pandas merge picks only one row
            if not kegg_rows.empty and not non_kegg.empty:
                # Build a lookup: LipidID -> KEGG_ID from fallback rows
                kegg_lookup = {}
                for _, kr in kegg_rows.iterrows():
                    lipid_id_key = _norm_key_inline(kr.get('LipidID'))
                    kegg_id_val = kr.get('KEGG_ID')
                    if lipid_id_key and kegg_id_val and str(kegg_id_val).strip() != '':
                        kegg_lookup[lipid_id_key] = str(kegg_id_val).strip()
                
                # Apply KEGG_ID from fallback to non-KEGG rows with matching LipidID (only if KEGG_ID is empty)
                for idx, row in non_kegg.iterrows():
                    lipid_id_key = row.get('__dedup_key')
                    if lipid_id_key in kegg_lookup:
                        existing_kegg = row.get('KEGG_ID')
                        if pd.isna(existing_kegg) or str(existing_kegg).strip() == '':
                            non_kegg.at[idx, 'KEGG_ID'] = kegg_lookup[lipid_id_key]
                            logger.debug(f"[PRE-MERGE DEDUP] Merged KEGG_ID {kegg_lookup[lipid_id_key]} from fallback into {row.get('LipidID')}")
                
                # Now keep only KEGG rows that have NO corresponding LipidMaps row (orphan KEGG fallbacks)
                kegg_keys = set(kegg_rows['__dedup_key'].tolist())
                non_kegg_keys = set(non_kegg['__dedup_key'].tolist())
                orphan_kegg_keys = kegg_keys - non_kegg_keys
                
                kegg_rows = kegg_rows[kegg_rows['__dedup_key'].isin(orphan_kegg_keys)].copy()
                logger.info(f"[PRE-MERGE DEDUP] Merged KEGG IDs into LipidMaps rows; {len(kegg_rows)} orphan KEGG fallbacks remain")
            
            # Add ID presence scoring to non-KEGG rows (now includes merged KEGG_IDs)
            id_cols_score = ['KEGG_ID', 'HMDB_ID', 'PubChem_CID', 'LipidMaps_ID']
            id_weights_score = {'KEGG_ID': 4, 'HMDB_ID': 3, 'PubChem_CID': 2, 'LipidMaps_ID': 1}
            
            for col in id_cols_score:
                if col in non_kegg.columns:
                    non_kegg[f'__{col}_present'] = non_kegg[col].apply(lambda v: 1 if pd.notna(v) and str(v).strip() != '' else 0)
            
            non_kegg['__id_score'] = 0
            for col, weight in id_weights_score.items():
                if f'__{col}_present' in non_kegg.columns:
                    non_kegg['__id_score'] += non_kegg[f'__{col}_present'] * weight
            
            # Sort and dedup non-KEGG by key, keeping best
            sort_cols_inline = ['__dedup_key', '__id_score']
            sort_asc_inline = [True, False]
            non_kegg = non_kegg.sort_values(by=sort_cols_inline, ascending=sort_asc_inline, kind='mergesort')
            non_kegg = non_kegg.drop_duplicates(subset=['__dedup_key'], keep='first')
            
            # Dedup orphan KEGG rows by LipidID (should already be unique, but safeguard)
            if not kegg_rows.empty and 'LipidID' in kegg_rows.columns:
                kegg_rows = kegg_rows.drop_duplicates(subset=['LipidID'], keep='first')
            
            # Recombine (now includes merged IDs in non_kegg + orphan KEGG rows)
            annotated_df = pd.concat([non_kegg, kegg_rows], ignore_index=True)
            
            # Drop helper columns AND _kegg_fallback_match marker
            helper_cols_inline = [col for col in annotated_df.columns if col.startswith('__') or col == '_kegg_fallback_match']
            annotated_df = annotated_df.drop(columns=helper_cols_inline, errors='ignore')
            
            logger.info(f"[LIPID MODE] Pre-merge dedup: {len(out_rows)} enriched rows -> {len(annotated_df)} unique (orphan KEGG: {len(kegg_rows)})")

        # Merge with original data, preferring LipidID for lipids; fallback to Name if needed
        # CRITICAL FIX: When annotated_df is empty (no LipidMaps matches), still use input_df as base
        # so KEGG IDs from kegg_id_map can be populated into final_df
        if annotated_df.empty:
            # No LipidMaps matches, but we still want to preserve input rows for KEGG population
            logger.info("[LIPID MODE] No LipidMaps matches; using input_df as base for KEGG population")
            final_df = input_df.copy()
        else:
            can_merge_on_lipidid = ('LipidID' in input_df.columns) and ('LipidID' in annotated_df.columns)
            if can_merge_on_lipidid:
                final_df = input_df.merge(annotated_df, on='LipidID', how='left', suffixes=('', '_annotated'))
            else:
                if 'Name' in input_df.columns and 'Name' in annotated_df.columns:
                    final_df = input_df.merge(annotated_df, on='Name', how='left', suffixes=('', '_annotated'))
                elif 'Name' in input_df.columns:
                    final_df = input_df.copy()
                else:
                    # If neither LipidID nor Name merge is possible, return annotations as-is
                    final_df = annotated_df

        # Update columns from enriched data
        enriched_columns = ['LipidID', 'LipidMaps_ID', 'PubChem_CID', 'KEGG_ID', 'KEGG_Match_Type', 'HMDB_ID', 'ChEBI_ID', 'CAS',
                   'SMILES', 'InChI', 'InChIKey', 'IUPAC_Name', 'Molecular_Formula', 'Molecular_Weight',
                   'Endogenous', 'Super_Class', 'Class', 'Sub_Class', 'annotation_sources', 'LipidMaps_ID_Match_Type']

        for col in enriched_columns:
            annotated_col = f"{col}_annotated"
            if annotated_col in final_df.columns:
                # Use enriched value if it's not empty, otherwise keep original
                mask = final_df[annotated_col].notna() & (final_df[annotated_col].astype(str).str.strip() != '')
                final_df.loc[mask, col] = final_df.loc[mask, annotated_col]
                # Drop the annotated column
                final_df.drop(columns=[annotated_col], inplace=True)

        # Ensure expected columns exist
        for col in enriched_columns:
            if col not in final_df.columns:
                final_df[col] = ''

        # CRITICAL: Populate KEGG_ID from kegg_id_map AFTER all merging/deduping completes
        # This ensures we populate the final 1:1 rows with input, not the intermediate annotated_df
        if kegg_id_map and 'LipidID' in final_df.columns:
            if 'KEGG_Match_Type' not in final_df.columns:
                final_df['KEGG_Match_Type'] = ''

            lipidid_series = final_df['LipidID'].astype(str).str.strip()
            assigned_rows = 0

            for lipid_id, info in kegg_id_map.items():
                mask = lipidid_series == lipid_id
                if not mask.any():
                    continue
                final_df.loc[mask, 'KEGG_ID'] = info.get('KEGG_ID', '')
                final_df.loc[mask, 'KEGG_Match_Type'] = info.get('match_type', 'Class')
                assigned_rows += int(mask.sum())

            populated_count = final_df['KEGG_ID'].astype(str).str.strip().ne('').sum()
            logger.info(
                f"[KEGG SEARCH] ✅ Populated KEGG_ID column in final output: {populated_count}/{len(final_df)} rows have KEGG IDs "
                f"({assigned_rows} set via Main_Class matches)"
            )

        # Reorder columns to place PubChem_CID beside KEGG_ID for clarity
        try:
            cols = list(final_df.columns)
            if 'PubChem_CID' in cols and 'KEGG_ID' in cols:
                cols.remove('PubChem_CID')
                kegg_idx = cols.index('KEGG_ID')
                cols.insert(kegg_idx + 1, 'PubChem_CID')
                final_df = final_df[cols]
        except Exception:
            pass

        # Clean KEGG IDs before saving (remove cpd: and cpd+ prefixes)
        if 'KEGG_ID' in final_df.columns:
            try:
                final_df['KEGG_ID'] = final_df['KEGG_ID'].apply(lambda x: clean_kegg_id(x) if pd.notna(x) else x)
                logger.debug("Cleaned KEGG IDs in final DataFrame")
            except Exception as e:
                logger.warning(f"Could not clean KEGG IDs: {e}")

        # Reorder columns to place LipidID first when present
        try:
            cols = list(final_df.columns)
            if 'LipidID' in cols:
                cols.remove('LipidID')
                cols.insert(0, 'LipidID')
                final_df = final_df[cols]
        except Exception:
            pass

        # Lipid-only MSI classification (Level 2 vs Level 3) based on LipidID specificity.
        try:
            if 'LipidID' in final_df.columns:
                final_df['MSI_levels'] = final_df['LipidID'].apply(_compute_lipid_msi_level)
            else:
                final_df['MSI_levels'] = ''
        except Exception as _msi_err:
            logger.warning(f"[LIPID MODE] Could not compute MSI_levels: {_msi_err}")
            if 'MSI_levels' not in final_df.columns:
                final_df['MSI_levels'] = ''

        # Persist caches so PubChem/KEGG/HMDB lookups are reused next runs
        try:
            self._save_cache(self.pubchem_cache, "pubchem_cache.pkl")
            self._save_cache(self.kegg_cache, "kegg_cache.pkl")
            self._save_cache(self.hmdb_cache, "hmdb_cache.pkl")
            self._save_cache(self.lipidmaps_cache, "lipidmaps_cache.pkl")
        except Exception as e:
            logger.warning(f"[LIPID MODE] Could not save caches: {e}")

        # NOTE: Dedup already done before merge (see pre-merge dedup above)
        # This section is kept for final safety check in case merge produced any duplicates
        try:
            before_final_dedup = len(final_df)
            if 'LipidID' in final_df.columns:
                final_df = final_df.drop_duplicates(subset=['LipidID'], keep='first')
                after_final_dedup = len(final_df)
                if before_final_dedup != after_final_dedup:
                    logger.warning(f"[LIPID MODE] Post-merge dedup removed {before_final_dedup - after_final_dedup} unexpected duplicates")
                else:
                    logger.info(f"[LIPID MODE] Post-merge: {len(final_df)} unique rows (no duplicates from merge)")
            else:
                logger.info(f"[LIPID MODE] Post-merge: {len(final_df)} rows")
        except Exception as _dedup_e:
            logger.warning(f"[LIPID MODE] Post-merge dedup check failed: {_dedup_e}")

        # FILTER: Keep only rows that have at least one of the primary IDs
        # (LipidMaps_ID, KEGG_ID, PubChem_CID, HMDB_ID)
        primary_id_cols = ['LipidMaps_ID', 'KEGG_ID', 'PubChem_CID', 'HMDB_ID']
        def _has_primary_id(row):
            for col in primary_id_cols:
                val = row.get(col, '')
                if pd.notna(val) and str(val).strip() != '':
                    return True
            return False
        
        total_before_filter = len(final_df)
        # Apply filter to keep only rows with at least one primary ID unless user disables filtering
        if getattr(self, 'skip_id_filtering', False):
            logger.info("[LIPID MODE] Skip ID filtering enabled - keeping all rows regardless of ID presence")
        else:
            mask_with_ids = final_df.apply(_has_primary_id, axis=1)
            if mask_with_ids.sum() == 0:
                logger.warning("[LIPID MODE] No IDs found after annotation - keeping all rows to avoid empty output")
            else:
                final_df = final_df[mask_with_ids].reset_index(drop=True)
                total_after_filter = len(final_df)
                filtered_out = total_before_filter - total_after_filter
                logger.info(f"[LIPID MODE] Filtering: {total_before_filter} rows → {total_after_filter} rows with IDs (removed {filtered_out} without IDs)")
        
        # Report on all ID columns for summary
        id_cols = ['LipidMaps_ID', 'KEGG_ID', 'PubChem_CID', 'HMDB_ID', 'ChEBI_ID', 'CAS']
        def _has_any_id(row):
            for col in id_cols:
                val = row.get(col, '')
                if pd.notna(val) and str(val).strip() != '':
                    return True
            return False
        total_post_dedup = len(final_df)
        with_ids = int(final_df.apply(_has_any_id, axis=1).sum()) if total_post_dedup else 0
        logger.info(f"[LIPID MODE] Post-filter rows: {total_post_dedup}; rows with at least one ID: {with_ids}")

        try:
            # Remove empty Name column if it exists and is completely empty
            if 'Name' in final_df.columns:
                name_values = final_df['Name'].fillna('').astype(str).str.strip()
                if (name_values == '').all():
                    final_df = final_df.drop(columns=['Name'])
                    logger.info("[LIPID MODE] Removed empty Name column from output")
            
            logger.info(f"[LIPID MODE] Saving deduplicated annotated results to: {self.output_file}")
            
            # If polarity sheets exist, map the annotated IDs back to them
            if 'polarity_sheets' in locals() and polarity_sheets:
                logger.info(f"[LIPID MODE] Mapping annotated IDs to polarity sheets...")
                with pd.ExcelWriter(self.output_file, engine='openpyxl') as writer:
                    # Write annotated first sheet
                    sheet_name = 'Annotated_Lipids'
                    final_df.to_excel(writer, sheet_name=sheet_name, index=False)
                    logger.info(f"[LIPID MODE]   ✓ Wrote '{sheet_name}' with {len(final_df)} rows")
                    
                    # Map annotated IDs to each polarity sheet
                    # Create ID mapping dictionary: LipidID -> annotated row data
                    if 'LipidID' in final_df.columns:
                        id_columns = [col for col in final_df.columns if col not in ['LipidID']]
                        
                        for sheet_name, sheet_df in polarity_sheets.items():
                            if 'LipidID' in sheet_df.columns:
                                # Merge annotated data with polarity sheet on LipidID
                                merged_sheet = sheet_df.merge(
                                    final_df[['LipidID'] + id_columns], 
                                    on='LipidID', 
                                    how='left', 
                                    suffixes=('', '_annotated')
                                )
                                
                                # Update original columns with annotated values where available
                                for col in id_columns:
                                    if col in sheet_df.columns:
                                        annotated_col = f"{col}_annotated"
                                        if annotated_col in merged_sheet.columns:
                                            # Use annotated value if not empty, otherwise keep original
                                            mask = merged_sheet[annotated_col].notna() & (merged_sheet[annotated_col].astype(str).str.strip() != '')
                                            merged_sheet.loc[mask, col] = merged_sheet.loc[mask, annotated_col]
                                            merged_sheet.drop(columns=[annotated_col], inplace=True)
                                    elif f"{col}_annotated" in merged_sheet.columns:
                                        # New column from annotation
                                        merged_sheet[col] = merged_sheet[f"{col}_annotated"]
                                        merged_sheet.drop(columns=[f"{col}_annotated"], inplace=True)
                                
                                merged_sheet.to_excel(writer, sheet_name=sheet_name, index=False)
                                logger.info(f"[LIPID MODE]   ✓ Wrote '{sheet_name}' with {len(merged_sheet)} rows")
                            else:
                                # No LipidID column to merge on, write original sheet
                                sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)
                                logger.info(f"[LIPID MODE]   ⚠️  '{sheet_name}' has no LipidID column, wrote original data")
                    else:
                        # No LipidID in annotated data, just write polarity sheets as-is
                        for sheet_name, sheet_df in polarity_sheets.items():
                            sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)
                            logger.info(f"[LIPID MODE]   ⚠️  No LipidID mapping possible, wrote original '{sheet_name}'")
                    
                    # Write class sheets AS-IS (no ID annotation needed)
                    if 'class_sheets' in locals() and class_sheets:
                        logger.info(f"[LIPID MODE] Writing {len(class_sheets)} class sheet(s) without ID annotation...")
                        for sheet_name, sheet_df in class_sheets.items():
                            sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)
                            logger.info(f"[LIPID MODE]   ✓ Wrote class sheet '{sheet_name}' with {len(sheet_df)} rows (preserved as-is)")
                
                logger.info(f"[LIPID MODE] Saved annotated results with {len(polarity_sheets)} polarity sheet(s) and {len(class_sheets) if 'class_sheets' in locals() else 0} class sheet(s)")
            else:
                # No polarity sheets, just write the annotated data as single sheet
                final_df.to_excel(self.output_file, index=False)
                logger.info(f"[LIPID MODE] Deduplicated rows saved: {len(final_df)}")
        except Exception as e:
            logger.error(f"[LIPID MODE] Failed to save annotated results: {e}")
            import traceback
            logger.error(f"[LIPID MODE] Traceback: {traceback.format_exc()}")

        # Summary
        # Summaries: report original unique input LipidIDs and unique matched by source
        input_total = len(input_df) if 'input_df' in locals() else None
        unique_input_lipids = None
        if input_df is not None:
            if 'LipidID' in input_df.columns:
                unique_input_lipids = int(input_df['LipidID'].dropna().astype(str).str.strip().nunique())
            else:
                unique_input_lipids = input_total
        # Report post-dedup row count for matched rows
        matched_count = len(final_df)
        total = matched_count

        def _count_nonempty_col(df, col):
            try:
                return int((df[col].notna() & (df[col].astype(str).str.strip() != '')).sum())
            except Exception:
                return 0

        lipid_ids = _count_nonempty_col(final_df, 'LipidMaps_ID') if 'LipidMaps_ID' in final_df.columns else 0
        kegg_hits = _count_nonempty_col(final_df, 'KEGG_ID')
        pubchem_hits = _count_nonempty_col(final_df, 'PubChem_CID')
        hmdb_hits = _count_nonempty_col(final_df, 'HMDB_ID')
        logger.info("="*60)
        logger.info("LIPID ID ANNOTATION SUMMARY")
        logger.info("="*60)
        if input_total is not None:
            logger.info(f"Original input rows: {input_total}")
        if unique_input_lipids is not None:
            logger.info(f"Original unique LipidIDs: {unique_input_lipids}")
        logger.info(f"Matched rows (primary output): {matched_count}")
        # Unique match counts by source (LipidMaps vs KEGG fallback)
        try:
            if ('_kegg_fallback_match' in final_df.columns):
                # Build masks on final_df safely
                kegg_mask = final_df['_kegg_fallback_match'].fillna(False).astype(bool)
                lm_mask = ~kegg_mask
                lm_unique = int(final_df.loc[lm_mask, 'LipidID'].dropna().astype(str).str.strip().nunique()) if 'LipidID' in final_df.columns else 0
                kegg_unique = int(final_df.loc[kegg_mask, 'LipidID'].dropna().astype(str).str.strip().nunique()) if 'LipidID' in final_df.columns else 0
                logger.info(f"Unique Lipids matched - LipidMaps: {lm_unique}, KEGG fallback: {kegg_unique}, Combined unique: {lm_unique + kegg_unique}")
        except Exception:
            pass
        logger.info(f"LipidMaps ID found: {lipid_ids} ({(lipid_ids/(input_total if input_total else 1)*100 if input_total else 0):.1f}% of input) ")
        if self.lipid_mode:
            logger.info(f"PubChem cross-refs found: {kegg_hits + hmdb_hits} (KEGG: {kegg_hits}, HMDB: {hmdb_hits})")
            try:
                # Report how many rows have at least one ID (post-dedup)
                id_cols = ['LipidMaps_ID', 'KEGG_ID', 'PubChem_CID', 'HMDB_ID', 'ChEBI_ID', 'CAS']
                has_any_id = final_df[id_cols].apply(lambda s: s.astype(str).str.strip() != '' , axis=0).any(axis=1)
                logger.info(f"Rows with ≥1 ID: {int(has_any_id.sum())} / {len(final_df)}")
            except Exception:
                pass
        else:
            logger.info(f"KEGG ID found: {kegg_hits} ({(kegg_hits/total*100 if total else 0):.1f}%)")
            logger.info(f"PubChem CID found: {pubchem_hits} ({(pubchem_hits/total*100 if total else 0):.1f}%)")
            logger.info(f"HMDB ID found: {hmdb_hits} ({(hmdb_hits/total*100 if total else 0):.1f}%)")
        logger.info("="*60)

        return final_df

    def _load_cache(self, filename: str) -> Dict:
        """Load cache from pickle file"""
        if os.path.exists(filename):
            try:
                with open(filename, 'rb') as f:
                    return pickle.load(f)
            except Exception as e:
                logger.warning(f"Could not load cache {filename}: {e}")
        return {}
    
    def _save_cache(self, cache: Dict, filename: str):
        """Save cache to pickle file"""
        try:
            with open(filename, 'wb') as f:
                pickle.dump(cache, f)
        except Exception as e:
            logger.warning(f"Could not save cache {filename}: {e}")
    
    def _load_mzcloud_cache(self) -> pd.DataFrame:
        """Load MZCloud cache from data cleaner"""
        cache_file = "mzcloud_cache.pkl"
        try:
            if os.path.exists(cache_file):
                with open(cache_file, 'rb') as f:
                    cache_data = pickle.load(f)
                    if isinstance(cache_data, pd.DataFrame) and not cache_data.empty:
                        logger.info(f"Loaded {len(cache_data)} entries from MZCloud cache")
                        return cache_data
                    else:
                        logger.info(f"MZCloud cache file exists but is empty or invalid")
                        return pd.DataFrame()
            else:
                logger.info(f"MZCloud cache file not found")
                return pd.DataFrame()
        except Exception as e:
            logger.warning(f"Error loading MZCloud cache: {e}")
            return pd.DataFrame()

    def _sanitize_caches(self):
        """Ensure all caches are in proper dictionary format"""
        cache_names = ['pubchem_cache', 'kegg_cache', 'hmdb_cache', 'lipidmaps_cache']
            
        for cache_name in cache_names:
            cache = getattr(self, cache_name)
            if not isinstance(cache, dict):
                logger.warning(f"{cache_name} is not a dict, resetting to empty dict")
                setattr(self, cache_name, {})

    def _wait_rate_limit(self, service: str):
        """Thread-safe rate limiting for API calls"""
        if service == 'pubchem':
            lock = self._pubchem_lock
            min_interval = self._pubchem_min_interval
            last_attr = '_pubchem_last_call'
        elif service == 'kegg':
            lock = self._kegg_lock
            min_interval = self._kegg_min_interval
            last_attr = '_kegg_last_call'
        else:
            return

        with lock:
            now = time.time()
            last = getattr(self, last_attr, 0.0) or 0.0
            elapsed = now - last
            if elapsed < min_interval:
                to_sleep = min_interval - elapsed
                logger.debug(f"[RATE] Sleeping {to_sleep:.3f}s for {service}")
                time.sleep(to_sleep)
            setattr(self, last_attr, time.time())

    def _normalize_name(self, s: Optional[str]) -> str:
        """Case-insensitive, whitespace-collapsed string for exact name checks."""
        if not s:
            return ''
        return ' '.join(str(s).strip().lower().split())

    def _normalize_formula(self, s: Optional[str]) -> str:
        """Uppercased, whitespace-removed formula for exact formula checks."""
        if not s:
            return ''
        return str(s).replace(' ', '').upper()

    def _strip_stereochemistry_prefix(self, name: str) -> Optional[str]:
        """
        Strip stereochemistry prefixes from metabolite names for fallback search.
        
        Handles:
        - Symbol-based prefixes: "(-)-", "(+)-", "(±)-", "(R)-", "(S)-", "(E)-", "(Z)-"
        - Short letter prefixes (≤2 chars): "L-", "D-", "DL-", "N-", "O-", "S-"
        
        Returns stripped name if a prefix was found, None otherwise.
        The function is strict about letter prefixes (max 2 letters) to avoid
        accidentally removing meaningful parts of compound names.
        """
        import re
        
        if not name or not isinstance(name, str):
            return None
        
        original_name = name.strip()
        
        # Symbol-based stereochemistry prefixes (parenthesized)
        # These can have more characters inside parentheses
        symbol_prefixes = [
            r'^\(\-\)\-',      # (-)-
            r'^\(\+\)\-',      # (+)-
            r'^\(±\)\-',       # (±)-
            r'^\(R\)\-',       # (R)-
            r'^\(S\)\-',       # (S)-
            r'^\(E\)\-',       # (E)-
            r'^\(Z\)\-',       # (Z)-
            r'^\(R,R\)\-',     # (R,R)-
            r'^\(S,S\)\-',     # (S,S)-
            r'^\(R,S\)\-',     # (R,S)-
            r'^\(S,R\)\-',     # (S,R)-
        ]
        
        # Short letter prefixes (max 2 letters, case insensitive)
        # Only strip if it's exactly 1-2 letters followed by a hyphen
        letter_prefix_pattern = r'^([A-Za-z]{1,2})\-'
        
        # Try symbol prefixes first
        for prefix_pattern in symbol_prefixes:
            if re.match(prefix_pattern, original_name, re.IGNORECASE):
                stripped = re.sub(prefix_pattern, '', original_name, count=1, flags=re.IGNORECASE)
                if stripped and stripped != original_name:
                    return stripped.strip()
        
        # Try short letter prefix (strict: max 2 letters)
        match = re.match(letter_prefix_pattern, original_name)
        if match:
            prefix_letters = match.group(1)
            # Only strip common stereochemistry/locant prefixes
            common_prefixes = ['l', 'd', 'dl', 'ld', 'n', 'o', 's', 'r', 'p', 'm']
            if prefix_letters.lower() in common_prefixes:
                stripped = original_name[len(prefix_letters) + 1:]  # +1 for the hyphen
                if stripped and stripped != original_name:
                    return stripped.strip()
        
        return None

    def search_pubchem_api(self, name: str, expected_formula: Optional[str] = None) -> Dict[str, Any]:
        """
        Search PubChem API for metabolite IDs and cross-references
        Priority 1: API-first approach
        """
        name_key = self._normalize_name(name)
        formula_key = self._normalize_formula(expected_formula)
        cache_key = f"pubchem_api_{name_key}|formula={formula_key}"
        
        if cache_key in self.pubchem_cache:
            logger.debug(f"[CACHE] PubChem API cache hit for: {name}")
            cached = self.pubchem_cache[cache_key]
            # Normalize older cache entries that may lack newer keys
            if isinstance(cached, dict):
                cached.setdefault('LipidMaps_ID', '')
                cached.setdefault('found', False)
            return cached

        self._wait_rate_limit('pubchem')
        
        result = {
            'PubChem_CID': '',
            'SMILES': '',
            'InChI': '',
            'InChIKey': '',
            'IUPAC_Name': '',
            'Molecular_Formula': '',
            'Molecular_Weight': '',
            'ChEBI_ID': '',
            'KEGG_ID': '',
            'HMDB_ID': '',
            'CAS': '',
            'LipidMaps_ID': '',
            'found': False
        }

        logger.info(f"[PUBCHEM] Searching for: {name}")

        try:
            # Search by name
            base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name"
            url = f"{base_url}/{name}/cids/JSON"
            
            response = requests.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                if 'IdentifierList' in data and data['IdentifierList']['CID']:
                    cids = data['IdentifierList']['CID']
                    
                    # Use first CID (most relevant match)
                    cid = cids[0]
                    
                    # Get detailed information for this CID
                    details_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey,IUPACName/JSON"
                    details_response = requests.get(details_url, timeout=10)
                    
                    accepted_reason = "name_match"
                    # If formula is provided, check if it matches
                    if expected_formula and details_response.status_code == 200:
                        details_data = details_response.json()
                        if 'PropertyTable' in details_data and details_data['PropertyTable']['Properties']:
                            prop = details_data['PropertyTable']['Properties'][0]
                            pubchem_formula = self._normalize_formula(prop.get('MolecularFormula', ''))
                            
                            if formula_key and pubchem_formula != formula_key:
                                logger.info(f"[PUBCHEM] Formula mismatch for {name}: expected {formula_key}, got {pubchem_formula}")
                                result['found'] = False
                                self.pubchem_cache[cache_key] = result
                                return result
                            elif formula_key and pubchem_formula == formula_key:
                                accepted_reason = "name_and_formula_match"
                    
                    if details_response.status_code == 200:
                        details_data = details_response.json()
                        
                        if 'PropertyTable' in details_data and details_data['PropertyTable']['Properties']:
                            prop = details_data['PropertyTable']['Properties'][0]
                            
                            result.update({
                                'PubChem_CID': str(cid),
                                'SMILES': prop.get('CanonicalSMILES', ''),
                                'InChI': prop.get('InChI', ''),
                                'InChIKey': prop.get('InChIKey', ''),
                                'IUPAC_Name': prop.get('IUPACName', ''),
                                'Molecular_Formula': prop.get('MolecularFormula', ''),
                                'Molecular_Weight': str(prop.get('MolecularWeight', ''))
                            })

                    # Get cross-references to other databases using helper (captures LipidMaps IDs)
                    try:
                        xrefs = fetch_pubchem_cross_references(int(cid))
                        # Only update empty fields
                        if not result.get('KEGG_ID') and xrefs.get('KEGG'):
                            result['KEGG_ID'] = xrefs.get('KEGG')
                        if not result.get('HMDB_ID') and xrefs.get('HMDB'):
                            result['HMDB_ID'] = xrefs.get('HMDB')
                        if not result.get('ChEBI_ID') and xrefs.get('ChEBI'):
                            result['ChEBI_ID'] = xrefs.get('ChEBI')
                        if not result.get('CAS') and xrefs.get('CAS'):
                            result['CAS'] = xrefs.get('CAS')
                        if not result.get('LipidMaps_ID') and xrefs.get('LipidMaps'):
                            result['LipidMaps_ID'] = xrefs.get('LipidMaps')
                    except Exception:
                        # Best-effort enrichment; fall back to dedicated endpoints below
                        pass

                    # 2) Dedicated xref endpoints for higher fidelity
                    try:
                        if not result['KEGG_ID']:
                            kegg_xref = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/KEGG/JSON"
                            rx = requests.get(kegg_xref, timeout=10)
                            if rx.status_code == 200:
                                jd = rx.json()
                                info = jd.get('InformationList', {}).get('Information', [])
                                if info and 'KEGG' in info[0] and info[0]['KEGG']:
                                    result['KEGG_ID'] = info[0]['KEGG'][0]
                        if not result['HMDB_ID']:
                            hmdb_xref = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/HMDB/JSON"
                            rx = requests.get(hmdb_xref, timeout=10)
                            if rx.status_code == 200:
                                jd = rx.json()
                                info = jd.get('InformationList', {}).get('Information', [])
                                if info and 'HMDB' in info[0] and info[0]['HMDB']:
                                    result['HMDB_ID'] = info[0]['HMDB'][0]
                        if not result['ChEBI_ID']:
                            chebi_xref = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/ChEBI/JSON"
                            rx = requests.get(chebi_xref, timeout=10)
                            if rx.status_code == 200:
                                jd = rx.json()
                                info = jd.get('InformationList', {}).get('Information', [])
                                if info and 'ChEBI' in info[0] and info[0]['ChEBI']:
                                    result['ChEBI_ID'] = info[0]['ChEBI'][0]
                        if not result['CAS']:
                            cas_xref = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RN/JSON"
                            rx = requests.get(cas_xref, timeout=10)
                            if rx.status_code == 200:
                                jd = rx.json()
                                info = jd.get('InformationList', {}).get('Information', [])
                                if info and 'RN' in info[0] and info[0]['RN']:
                                    result['CAS'] = info[0]['RN'][0]
                    except Exception:
                        # Xref enrichment is best-effort
                        pass

                    result['found'] = True
                    logger.info(f"[OK] PubChem accepted CID {cid} for: {name} (by {accepted_reason})")
                else:
                    if formula_key:
                        logger.info(f"[MISS] PubChem: exact name+formula required, none matched for '{name}'")
                    else:
                        logger.info(f"[MISS] PubChem: exact name required, none matched for '{name}'")
            else:
                logger.info(f"[MISS] PubChem name search failed for: {name}")
                
                # Fallback: Try searching by IUPAC name if the input name looks like an IUPAC name
                # IUPAC names typically contain characteristic patterns like parentheses with stereochemistry
                if any(pattern in name.lower() for pattern in ['(', ')', '-', 'yl', 'ic acid', 'ol', 'one', 'ene', 'ane']):
                    logger.info(f"[FALLBACK] Trying IUPAC name search for: {name}")
                    try:
                        # Use compound/name endpoint but specify it might be an IUPAC name
                        iupac_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
                        iupac_response = requests.get(iupac_url, timeout=10)
                        
                        if iupac_response.status_code == 200:
                            iupac_data = iupac_response.json()
                            if 'IdentifierList' in iupac_data and iupac_data['IdentifierList']['CID']:
                                cids = iupac_data['IdentifierList']['CID']
                                cid = cids[0]
                                
                                # Get detailed information for this CID
                                details_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey,IUPACName/JSON"
                                details_response = requests.get(details_url, timeout=10)
                                
                                accepted_reason = "iupac_name_match"
                                # If formula is provided, check if it matches
                                if expected_formula and details_response.status_code == 200:
                                    details_data = details_response.json()
                                    if 'PropertyTable' in details_data and details_data['PropertyTable']['Properties']:
                                        prop = details_data['PropertyTable']['Properties'][0]
                                        pubchem_formula = self._normalize_formula(prop.get('MolecularFormula', ''))
                                        
                                        if formula_key and pubchem_formula != formula_key:
                                            logger.info(f"[IUPAC FALLBACK] Formula mismatch for {name}: expected {formula_key}, got {pubchem_formula}")
                                        else:
                                            if formula_key and pubchem_formula == formula_key:
                                                accepted_reason = "iupac_name_and_formula_match"
                                            
                                            # Process the successful IUPAC match
                                            if details_response.status_code == 200:
                                                details_data = details_response.json()
                                                
                                                if 'PropertyTable' in details_data and details_data['PropertyTable']['Properties']:
                                                    prop = details_data['PropertyTable']['Properties'][0]
                                                    
                                                    result.update({
                                                        'PubChem_CID': str(cid),
                                                        'SMILES': prop.get('CanonicalSMILES', ''),
                                                        'InChI': prop.get('InChI', ''),
                                                        'InChIKey': prop.get('InChIKey', ''),
                                                        'IUPAC_Name': prop.get('IUPACName', ''),
                                                        'Molecular_Formula': prop.get('MolecularFormula', ''),
                                                        'Molecular_Weight': str(prop.get('MolecularWeight', ''))
                                                    })

                                            # Get cross-references (same as primary search)
                                            try:
                                                xrefs = fetch_pubchem_cross_references(int(cid))
                                                if not result.get('KEGG_ID') and xrefs.get('KEGG'):
                                                    result['KEGG_ID'] = xrefs.get('KEGG')
                                                if not result.get('HMDB_ID') and xrefs.get('HMDB'):
                                                    result['HMDB_ID'] = xrefs.get('HMDB')
                                                if not result.get('ChEBI_ID') and xrefs.get('ChEBI'):
                                                    result['ChEBI_ID'] = xrefs.get('ChEBI')
                                                if not result.get('CAS') and xrefs.get('CAS'):
                                                    result['CAS'] = xrefs.get('CAS')
                                                if not result.get('LipidMaps_ID') and xrefs.get('LipidMaps'):
                                                    result['LipidMaps_ID'] = xrefs.get('LipidMaps')
                                            except Exception:
                                                pass

                                            result['found'] = True
                                            logger.info(f"[OK] PubChem IUPAC fallback success - CID {cid} for: {name} (by {accepted_reason})")
                                            
                    except Exception as e:
                        logger.debug(f"[IUPAC FALLBACK] Exception during IUPAC search for {name}: {e}")
                        
                if not result['found']:
                    logger.info(f"[MISS] PubChem no results (name + IUPAC fallback) for: {name}")
        except Exception as e:
            logger.error(f"[ERROR] PubChem API exception for {name}: {e}")

        self.pubchem_cache[cache_key] = result
        return result

    def search_pubchem_by_inchikey(self, inchikey: str) -> Dict[str, Any]:
        """Search PubChem by InChIKey"""
        cache_key = f"pubchem_inchikey_{inchikey.strip()}"
        if cache_key in self.pubchem_cache:
            return self.pubchem_cache[cache_key]
        
        self._wait_rate_limit('pubchem')
        result = {'found': False, 'PubChem_CID': '', 'SMILES': '', 'InChI': '', 'InChIKey': '', 'IUPAC_Name': '', 'Molecular_Formula': '', 'Molecular_Weight': '', 'ChEBI_ID': '', 'KEGG_ID': '', 'HMDB_ID': '', 'CAS': '', 'LipidMaps_ID': ''}
        
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if 'IdentifierList' in data and data['IdentifierList']['CID']:
                    cid = data['IdentifierList']['CID'][0]
                    result.update(self._get_pubchem_details(cid))
                    result['found'] = True
                    logger.info(f"[OK] PubChem found CID {cid} by InChIKey: {inchikey}")
        except Exception as e:
            logger.error(f"[ERROR] PubChem InChIKey search failed for {inchikey}: {e}")
        
        self.pubchem_cache[cache_key] = result
        return result

    def search_pubchem_by_inchi(self, inchi: str) -> Dict[str, Any]:
        """Search PubChem by InChI"""
        cache_key = f"pubchem_inchi_{self._normalize_name(inchi)}"
        if cache_key in self.pubchem_cache:
            return self.pubchem_cache[cache_key]
        
        self._wait_rate_limit('pubchem')
        result = {'found': False, 'PubChem_CID': '', 'SMILES': '', 'InChI': '', 'InChIKey': '', 'IUPAC_Name': '', 'Molecular_Formula': '', 'Molecular_Weight': '', 'ChEBI_ID': '', 'KEGG_ID': '', 'HMDB_ID': '', 'CAS': '', 'LipidMaps_ID': ''}
        
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/{inchi}/cids/JSON"
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if 'IdentifierList' in data and data['IdentifierList']['CID']:
                    cid = data['IdentifierList']['CID'][0]
                    result.update(self._get_pubchem_details(cid))
                    result['found'] = True
                    logger.info(f"[OK] PubChem found CID {cid} by InChI")
        except Exception as e:
            logger.error(f"[ERROR] PubChem InChI search failed: {e}")
        
        self.pubchem_cache[cache_key] = result
        return result

    def search_pubchem_by_smiles(self, smiles: str) -> Dict[str, Any]:
        """Search PubChem by SMILES"""
        cache_key = f"pubchem_smiles_{self._normalize_name(smiles)}"
        if cache_key in self.pubchem_cache:
            return self.pubchem_cache[cache_key]
        
        self._wait_rate_limit('pubchem')
        result = {'found': False, 'PubChem_CID': '', 'SMILES': '', 'InChI': '', 'InChIKey': '', 'IUPAC_Name': '', 'Molecular_Formula': '', 'Molecular_Weight': '', 'ChEBI_ID': '', 'KEGG_ID': '', 'HMDB_ID': '', 'CAS': '', 'LipidMaps_ID': ''}
        
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if 'IdentifierList' in data and data['IdentifierList']['CID']:
                    cid = data['IdentifierList']['CID'][0]
                    result.update(self._get_pubchem_details(cid))
                    result['found'] = True
                    logger.info(f"[OK] PubChem found CID {cid} by SMILES")
        except Exception as e:
            logger.error(f"[ERROR] PubChem SMILES search failed: {e}")
        
        self.pubchem_cache[cache_key] = result
        return result

    def search_pubchem_by_chebi(self, chebi_id: str) -> Dict[str, Any]:
        """Search PubChem by ChEBI ID"""
        cache_key = f"pubchem_chebi_{chebi_id.strip()}"
        if cache_key in self.pubchem_cache:
            return self.pubchem_cache[cache_key]
        
        self._wait_rate_limit('pubchem')
        result = {'found': False, 'PubChem_CID': '', 'SMILES': '', 'InChI': '', 'InChIKey': '', 'IUPAC_Name': '', 'Molecular_Formula': '', 'Molecular_Weight': '', 'ChEBI_ID': '', 'KEGG_ID': '', 'HMDB_ID': '', 'CAS': '', 'LipidMaps_ID': ''}
        
        try:
            # Try ChEBI cross-reference lookup
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/ChEBI/{chebi_id}/cids/JSON"
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if 'IdentifierList' in data and data['IdentifierList']['CID']:
                    cid = data['IdentifierList']['CID'][0]
                    result.update(self._get_pubchem_details(cid))
                    result['found'] = True
                    logger.info(f"[OK] PubChem found CID {cid} by ChEBI: {chebi_id}")
        except Exception as e:
            logger.error(f"[ERROR] PubChem ChEBI search failed for {chebi_id}: {e}")
        
        self.pubchem_cache[cache_key] = result
        return result

    def _get_pubchem_details(self, cid: str) -> Dict[str, Any]:
        """Get detailed properties for a PubChem CID"""
        details = {'PubChem_CID': str(cid), 'SMILES': '', 'InChI': '', 'InChIKey': '', 'IUPAC_Name': '', 'Molecular_Formula': '', 'Molecular_Weight': '', 'ChEBI_ID': '', 'KEGG_ID': '', 'HMDB_ID': '', 'CAS': '', 'LipidMaps_ID': ''}
        
        try:
            # Get properties
            details_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey,IUPACName/JSON"
            details_response = requests.get(details_url, timeout=10)
            
            if details_response.status_code == 200:
                details_data = details_response.json()
                if 'PropertyTable' in details_data and details_data['PropertyTable']['Properties']:
                    prop = details_data['PropertyTable']['Properties'][0]
                    details.update({
                        'SMILES': prop.get('CanonicalSMILES', ''),
                        'InChI': prop.get('InChI', ''),
                        'InChIKey': prop.get('InChIKey', ''),
                        'IUPAC_Name': prop.get('IUPACName', ''),
                        'Molecular_Formula': prop.get('MolecularFormula', ''),
                        'Molecular_Weight': str(prop.get('MolecularWeight', ''))
                    })

            # Get cross-references
            xrefs_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON"
            xrefs_response = requests.get(xrefs_url, timeout=10)
            
            if xrefs_response.status_code == 200:
                xrefs_data = xrefs_response.json()
                if 'InformationList' in xrefs_data and xrefs_data['InformationList']['Information']:
                    xrefs = xrefs_data['InformationList']['Information'][0]
                    if 'RegistryID' in xrefs and xrefs['RegistryID']:
                        for reg_id in xrefs['RegistryID']:
                            if reg_id.startswith('CHEBI') and not details['ChEBI_ID']:
                                details['ChEBI_ID'] = reg_id
                            elif reg_id.startswith('HMDB') and not details['HMDB_ID']:
                                details['HMDB_ID'] = reg_id
                            elif reg_id.startswith('LM') and not details['LipidMaps_ID']:
                                details['LipidMaps_ID'] = reg_id
                            elif reg_id.startswith('C') and len(reg_id) == 6 and not details['KEGG_ID']:
                                details['KEGG_ID'] = reg_id
                    if 'RN' in xrefs and xrefs['RN'] and not details['CAS']:
                        details['CAS'] = xrefs['RN'][0]

            # Try dedicated endpoints for better cross-reference coverage
            for endpoint, key in [('KEGG', 'KEGG_ID'), ('HMDB', 'HMDB_ID'), ('ChEBI', 'ChEBI_ID'), ('RN', 'CAS')]:
                if not details[key]:
                    try:
                        endpoint_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/{endpoint}/JSON"
                        endpoint_response = requests.get(endpoint_url, timeout=10)
                        if endpoint_response.status_code == 200:
                            endpoint_data = endpoint_response.json()
                            info = endpoint_data.get('InformationList', {}).get('Information', [])
                            if info and endpoint in info[0] and info[0][endpoint]:
                                details[key] = info[0][endpoint][0]
                    except Exception:
                        pass

        except Exception as e:
            logger.error(f"[ERROR] Failed to get PubChem details for CID {cid}: {e}")
        
        return details

    def search_kegg_api(self, name: str) -> Dict[str, Any]:
        """
        Search KEGG API for metabolite IDs
        Priority 2: Complement PubChem data
        
        Uses EXACT NAME MATCHING to ensure correct mapping.
        For class names like "ceramide", KEGG returns C00195 with "ceramide" as primary name.
        """
        cache_key = f"kegg_api_{self._normalize_name(name)}"
        
        if cache_key in self.kegg_cache:
            logger.debug(f"[CACHE] KEGG API cache hit for: {name}")
            return self.kegg_cache[cache_key]

        self._wait_rate_limit('kegg')
        
        result = {
            'KEGG_ID': '',
            'Molecular_Formula': '',
            'found': False
        }

        logger.info(f"[KEGG] Searching for: {name}")

        try:
            # Search in KEGG compound database
            search_url = f"https://rest.kegg.jp/find/compound/{name}"
            response = requests.get(search_url, timeout=10)
            
            if response.status_code == 200 and response.text.strip():
                lines = response.text.strip().split('\n')
                
                # Look for EXACT name matches only (no partial matches)
                best_match = None
                normalized_search = self._normalize_name(name)
                
                # KEGG returns results ranked by relevance
                # Check each result for exact match of the query name
                for line in lines:
                    if '\t' in line:
                        compound_id, compound_names = line.split('\t', 1)
                        # Split compound names by semicolon and normalize each
                        names_list = [self._normalize_name(n.strip()) for n in compound_names.split(';')]
                        
                        # Check if our search name exactly matches ANY of the compound names
                        # This handles cases like "ceramide" matching "ceramide" from KEGG entry
                        if normalized_search in names_list:
                            best_match = compound_id
                            logger.info(f"[KEGG] ✅ Exact match found: {compound_id} for '{name}'")
                            break
                
                # If no exact match found, reject all results (don't use fallback)
                if not best_match:
                    logger.debug(f"[KEGG] ❌ No exact name match for: {name} (found {len(lines)} results, but none matched exactly)")
                
                if best_match:
                    # Get detailed information for this compound
                    details_url = f"https://rest.kegg.jp/get/{best_match}"
                    details_response = requests.get(details_url, timeout=10)
                    
                    if details_response.status_code == 200:
                        details_text = details_response.text
                        
                        # Extract molecular formula
                        formula = ''
                        for line in details_text.split('\n'):
                            if line.startswith('FORMULA'):
                                formula = line.replace('FORMULA', '').strip()
                                break
                        
                        result.update({
                            'KEGG_ID': clean_kegg_id(best_match),
                            'Molecular_Formula': formula,
                            'found': True
                        })
                        
                        logger.info(f"[KEGG] ✅ Found {best_match} (formula: {formula}) for: {name}")
                    else:
                        logger.warning(f"[KEGG] Could not retrieve details for {best_match}")
                else:
                    logger.debug(f"[KEGG] No exact match results for: {name}")
            else:
                logger.debug(f"[KEGG] Empty response or error for: {name}")
                
        except Exception as e:
            logger.error(f"[ERROR] KEGG API exception for {name}: {e}")

        self.kegg_cache[cache_key] = result
        return result

    def search_hmdb_offline(self, name: str, existing_ids: Dict[str, Any]) -> Dict[str, Any]:
        """
        Search HMDB offline database for metabolite information
        Priority 3: Fill gaps and provide biological classification
        Enhanced with ID-based cache key to avoid stale negatives
        """
        # Build enhanced cache key using available IDs to avoid negative cache collisions
        cache_identifiers = []
        for key in ['InChIKey', 'HMDB_ID', 'PubChem_CID', 'KEGG_ID', 'ChEBI_ID', 'SMILES', 'InChI']:
            val = existing_ids.get(key, '').strip()
            if val:
                cache_identifiers.append(f"{key}:{val}")
        
        # Include normalized name as fallback
        normalized_name = self._normalize_name(name)
        cache_identifiers.append(f"name:{normalized_name}")
        
        cache_key = f"hmdb_offline_{'|'.join(cache_identifiers)}"
        
        if cache_key in self.hmdb_cache:
            logger.debug(f"[CACHE] HMDB offline cache hit for: {name}")
            return self.hmdb_cache[cache_key]
        
        result = {
            'HMDB_ID': '',
            'InChIKey': '',
            'Molecular_Formula': '',
            'Endogenous': 'No',
            'Super_Class': '',
            'Class': '',
            'Sub_Class': '',
            'found': False
        }

        logger.info(f"[HMDB] Searching offline for: {name}")

        try:
            # Strategy 1: Match by InChIKey if available
            if existing_ids.get('InChIKey'):
                inchikey = existing_ids['InChIKey']
                inchikey_matches = self.hmdb_df[
                    self.hmdb_df['InChIKey'].astype(str).str.strip().str.upper() == str(inchikey).strip().upper()
                ]
                if not inchikey_matches.empty:
                    match = inchikey_matches.iloc[0]
                    result.update({
                        'HMDB_ID': str(match.get('HMDB', '')),
                        'InChIKey': str(match.get('InChIKey', '')),
                        'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                        'Endogenous': _normalize_yes_no(match.get('Is_Endogenous', 'No')),
                        'Super_Class': str(match.get('Super_Class', '')),
                        'Class': str(match.get('Class', '')),
                        'Sub_Class': str(match.get('Sub_Class', '')),
                        'found': True
                    })
                    logger.info(f"[OK] HMDB matched by InChIKey for: {name}")
                    self.hmdb_cache[cache_key] = result
                    return result

            # Strategy 2: Match by HMDB ID if available
            if existing_ids.get('HMDB_ID'):
                hmdb_id = existing_ids['HMDB_ID']
                hmdb_matches = self.hmdb_df[
                    self.hmdb_df['HMDB'].astype(str).str.strip().str.upper() == str(hmdb_id).strip().upper()
                ]
                if not hmdb_matches.empty:
                    match = hmdb_matches.iloc[0]
                    result.update({
                        'HMDB_ID': str(match.get('HMDB', '')),
                        'InChIKey': str(match.get('InChIKey', '')),
                        'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                        'Endogenous': _normalize_yes_no(match.get('Is_Endogenous', 'No')),
                        'Super_Class': str(match.get('Super_Class', '')),
                        'Class': str(match.get('Class', '')),
                        'Sub_Class': str(match.get('Sub_Class', '')),
                        'found': True
                    })
                    logger.info(f"[OK] HMDB matched by HMDB_ID for: {name}")
                    self.hmdb_cache[cache_key] = result
                    return result

            # Strategy 2b: Match by PubChem CID (if HMDB export includes it)
            if existing_ids.get('PubChem_CID'):
                cid = existing_ids['PubChem_CID']
                for col in ['PubChem_CID', 'PubChem CID', 'PubChem']:
                    if col in self.hmdb_df.columns:
                        cid_matches = self.hmdb_df[self.hmdb_df[col].astype(str).str.strip() == str(cid).strip()]
                        if not cid_matches.empty:
                            match = cid_matches.iloc[0]
                            result.update({
                                'HMDB_ID': str(match.get('HMDB', '')),
                                'InChIKey': str(match.get('InChIKey', '')),
                                'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                                'Endogenous': _normalize_yes_no(match.get('Is_Endogenous', 'No')),
                                'Super_Class': str(match.get('Super_Class', '')),
                                'Class': str(match.get('Class', '')),
                                'Sub_Class': str(match.get('Sub_Class', '')),
                                'found': True
                            })
                            logger.info(f"[OK] HMDB matched by PubChem CID for: {name}")
                            self.hmdb_cache[cache_key] = result
                            return result

            # Strategy 2c: Match by KEGG ID (if present in HMDB export)
            if existing_ids.get('KEGG_ID'):
                kegg_id = existing_ids['KEGG_ID']
                for col in ['KEGG_ID', 'KEGG', 'KEGG ID']:
                    if col in self.hmdb_df.columns:
                        kegg_matches = self.hmdb_df[self.hmdb_df[col].astype(str).str.strip() == str(kegg_id).strip()]
                        if not kegg_matches.empty:
                            match = kegg_matches.iloc[0]
                            result.update({
                                'HMDB_ID': str(match.get('HMDB', '')),
                                'InChIKey': str(match.get('InChIKey', '')),
                                'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                                'Endogenous': _normalize_yes_no(match.get('Is_Endogenous', 'No')),
                                'Super_Class': str(match.get('Super_Class', '')),
                                'Class': str(match.get('Class', '')),
                                'Sub_Class': str(match.get('Sub_Class', '')),
                                'found': True
                            })
                            logger.info(f"[OK] HMDB matched by KEGG ID for: {name}")
                            self.hmdb_cache[cache_key] = result
                            return result

            # Strategy 2d: Match by ChEBI ID (if present in HMDB export)
            if existing_ids.get('ChEBI_ID'):
                chebi_id = existing_ids['ChEBI_ID']
                for col in ['ChEBI_ID', 'ChEBI', 'CHEBI_ID', 'CHEBI']:
                    if col in self.hmdb_df.columns:
                        chebi_matches = self.hmdb_df[self.hmdb_df[col].astype(str).str.strip() == str(chebi_id).strip()]
                        if not chebi_matches.empty:
                            match = chebi_matches.iloc[0]
                            result.update({
                                'HMDB_ID': str(match.get('HMDB', '')),
                                'InChIKey': str(match.get('InChIKey', '')),
                                'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                                'Endogenous': _normalize_yes_no(match.get('Is_Endogenous', 'No')),
                                'Super_Class': str(match.get('Super_Class', '')),
                                'Class': str(match.get('Class', '')),
                                'Sub_Class': str(match.get('Sub_Class', '')),
                                'found': True
                            })
                            logger.info(f"[OK] HMDB matched by ChEBI ID for: {name}")
                            self.hmdb_cache[cache_key] = result
                            return result

            # Strategy 2e: Match by SMILES (if present in HMDB export)
            if existing_ids.get('SMILES'):
                smiles = existing_ids['SMILES']
                for col in ['SMILES', 'Smiles']:
                    if col in self.hmdb_df.columns:
                        smiles_matches = self.hmdb_df[self.hmdb_df[col].astype(str).str.strip() == str(smiles).strip()]
                        if not smiles_matches.empty:
                            match = smiles_matches.iloc[0]
                            result.update({
                                'HMDB_ID': str(match.get('HMDB', '')),
                                'InChIKey': str(match.get('InChIKey', '')),
                                'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                                'Endogenous': _normalize_yes_no(match.get('Is_Endogenous', 'No')),
                                'Super_Class': str(match.get('Super_Class', '')),
                                'Class': str(match.get('Class', '')),
                                'Sub_Class': str(match.get('Sub_Class', '')),
                                'found': True
                            })
                            logger.info(f"[OK] HMDB matched by SMILES for: {name}")
                            self.hmdb_cache[cache_key] = result
                            return result

            # Strategy 2f: Match by InChI (if present in HMDB export)
            if existing_ids.get('InChI'):
                inchi = existing_ids['InChI']
                for col in ['InChI', 'INCHI']:
                    if col in self.hmdb_df.columns:
                        inchi_matches = self.hmdb_df[self.hmdb_df[col].astype(str).str.strip() == str(inchi).strip()]
                        if not inchi_matches.empty:
                            match = inchi_matches.iloc[0]
                            result.update({
                                'HMDB_ID': str(match.get('HMDB', '')),
                                'InChIKey': str(match.get('InChIKey', '')),
                                'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                                'Endogenous': str(match.get('Is_Endogenous', '')),
                                'Super_Class': str(match.get('Super_Class', '')),
                                'Class': str(match.get('Class', '')),
                                'Sub_Class': str(match.get('Sub_Class', '')),
                                'found': True
                            })
                            logger.info(f"[OK] HMDB matched by InChI for: {name}")
                            self.hmdb_cache[cache_key] = result
                            return result

            # Strategy 3: Match by name
            name_normalized = name.lower().strip()
            name_matches = self.hmdb_df[
                self.hmdb_df['Name'].astype(str).str.lower().str.strip() == name_normalized
            ]
            
            if not name_matches.empty:
                match = name_matches.iloc[0]
                result.update({
                    'HMDB_ID': str(match.get('HMDB', '')),
                    'InChIKey': str(match.get('InChIKey', '')),
                    'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                    'Endogenous': str(match.get('Is_Endogenous', '')),
                    'Super_Class': str(match.get('Super_Class', '')),
                    'Class': str(match.get('Class', '')),
                    'Sub_Class': str(match.get('Sub_Class', '')),
                    'found': True
                })
                logger.info(f"[OK] HMDB matched by name for: {name}")
            else:
                # Strategy 3b: Try synonyms - use dedicated synonyms file if available (EXACT MATCH ONLY)
                if self.hmdb_synonyms_df is not None:
                    # Use the fast synonyms lookup table - exact match only to avoid false positives
                    synonym_matches = self.hmdb_synonyms_df[
                        self.hmdb_synonyms_df['Synonyms_lower'] == name_normalized
                    ]
                    
                    if not synonym_matches.empty:
                        # Get the HMDB ID and look up full record in main database
                        hmdb_id = synonym_matches.iloc[0]['HMDB']
                        main_match = self.hmdb_df[self.hmdb_df['HMDB'] == hmdb_id]
                        
                        if not main_match.empty:
                            match = main_match.iloc[0]
                            result.update({
                                'HMDB_ID': str(match.get('HMDB', '')),
                                'InChIKey': str(match.get('InChIKey', '')),
                                'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                                'Endogenous': str(match.get('Is_Endogenous', '')),
                                'Super_Class': str(match.get('Super_Class', '')),
                                'Class': str(match.get('Class', '')),
                                'Sub_Class': str(match.get('Sub_Class', '')),
                                'found': True
                            })
                            logger.info(f"[OK] HMDB matched by synonym (fast lookup) for: {name}")
                        else:
                            logger.info(f"[MISS] HMDB synonym found but no main record for: {name}")
                    else:
                        logger.info(f"[MISS] HMDB no synonym matches for: {name}")
                
                elif 'Synonyms' in self.hmdb_df.columns:
                    # Fallback to old method if synonyms file not available (exact match in pipe-separated list)
                    # Split synonyms and check for exact match to avoid false positives
                    def has_exact_synonym(syn_str):
                        if pd.isna(syn_str):
                            return False
                        synonyms = [s.strip().lower() for s in str(syn_str).split('|')]
                        return name_normalized in synonyms
                    
                    syn_mask = self.hmdb_df['Synonyms'].apply(has_exact_synonym)
                    synonym_matches = self.hmdb_df[syn_mask]
                    if not synonym_matches.empty:
                        match = synonym_matches.iloc[0]
                        result.update({
                            'HMDB_ID': str(match.get('HMDB', '')),
                            'InChIKey': str(match.get('InChIKey', '')),
                            'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                            'Endogenous': str(match.get('Is_Endogenous', '')),
                            'Super_Class': str(match.get('Super_Class', '')),
                            'Class': str(match.get('Class', '')),
                            'Sub_Class': str(match.get('Sub_Class', '')),
                            'found': True
                        })
                        logger.info(f"[OK] HMDB matched by synonym (fallback) for: {name}")
                    else:
                        logger.info(f"[MISS] HMDB no matches for: {name}")
                else:
                    logger.info(f"[MISS] HMDB no matches for: {name}")
            
            # Strategy 4: Try stripping stereochemistry prefixes if no match found yet
            # Only attempt if original search failed
            if not result['found']:
                stripped_name = self._strip_stereochemistry_prefix(name)
                if stripped_name and stripped_name.lower() != name_normalized:
                    logger.info(f"[RETRY] Trying stripped name: '{stripped_name}' (original: '{name}')")
                    
                    # Try name match with stripped name
                    stripped_normalized = stripped_name.lower().strip()
                    stripped_matches = self.hmdb_df[
                        self.hmdb_df['Name'].astype(str).str.lower().str.strip() == stripped_normalized
                    ]
                    
                    if not stripped_matches.empty:
                        match = stripped_matches.iloc[0]
                        result.update({
                            'HMDB_ID': str(match.get('HMDB', '')),
                            'InChIKey': str(match.get('InChIKey', '')),
                            'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                            'Endogenous': str(match.get('Is_Endogenous', '')),
                            'Super_Class': str(match.get('Super_Class', '')),
                            'Class': str(match.get('Class', '')),
                            'Sub_Class': str(match.get('Sub_Class', '')),
                            'found': True
                        })
                        logger.info(f"[OK] HMDB matched by stripped name for: {name} -> {stripped_name}")
                    
                    # Try synonym match with stripped name if name match failed
                    elif self.hmdb_synonyms_df is not None:
                        stripped_syn_matches = self.hmdb_synonyms_df[
                            self.hmdb_synonyms_df['Synonyms_lower'] == stripped_normalized
                        ]
                        
                        if not stripped_syn_matches.empty:
                            hmdb_id = stripped_syn_matches.iloc[0]['HMDB']
                            main_match = self.hmdb_df[self.hmdb_df['HMDB'] == hmdb_id]
                            
                            if not main_match.empty:
                                match = main_match.iloc[0]
                                result.update({
                                    'HMDB_ID': str(match.get('HMDB', '')),
                                    'InChIKey': str(match.get('InChIKey', '')),
                                    'Molecular_Formula': str(match.get('Molecular_Formula', '')),
                                    'Endogenous': str(match.get('Is_Endogenous', '')),
                                    'Super_Class': str(match.get('Super_Class', '')),
                                    'Class': str(match.get('Class', '')),
                                    'Sub_Class': str(match.get('Sub_Class', '')),
                                    'found': True
                                })
                                logger.info(f"[OK] HMDB matched by stripped synonym for: {name} -> {stripped_name}")
                        else:
                            logger.info(f"[MISS] HMDB no matches for stripped name: {stripped_name}")
                    else:
                        logger.info(f"[MISS] HMDB no matches for stripped name: {stripped_name}")

        except Exception as e:
            logger.error(f"[ERROR] HMDB offline search exception for {name}: {e}")

        self.hmdb_cache[cache_key] = result
        return result

    def search_mzcloud_cache(self, name: str, existing_ids: Dict[str, Any]) -> Dict[str, Any]:
        """
        Search MZCloud cache for additional metabolite information
        Priority 4: Endogenous status and additional cross-references
        """
        result = {
            'Endogenous': '',
            'additional_ids': {},
            'found': False
        }

        # Check if MZCloud cache is available
        if self.mzcloud_cache.empty:
            return result

    # logger.info(f"[CACHE] Searching MZCloud cache for: {name}")  # noisy per-metabolite log (commented out)

        try:
            matches = None
            
            # Strategy 1: Exact name match
            name_normalized = name.lower().strip()
            if 'Name' in self.mzcloud_cache.columns:
                exact_matches = self.mzcloud_cache[
                    self.mzcloud_cache['Name'].str.lower().str.strip() == name_normalized
                ]
                if not exact_matches.empty:
                    matches = exact_matches
                    logger.info(f"[OK] MZCloud exact name match for: {name}")
            
            # Strategy 2: InChIKey match (if available)
            if matches is None and existing_ids.get('InChIKey') and 'InChIKey' in self.mzcloud_cache.columns:
                inchikey = str(existing_ids['InChIKey']).strip()
                inchikey_matches = self.mzcloud_cache[self.mzcloud_cache['InChIKey'].astype(str).str.strip() == inchikey]
                if not inchikey_matches.empty:
                    matches = inchikey_matches
                    logger.info(f"[OK] MZCloud InChIKey match for: {name}")
            
            # Strategy 3: HMDB ID match (if available)
            if matches is None and existing_ids.get('HMDB_ID') and 'HMDB ID' in self.mzcloud_cache.columns:
                hmdb_id = str(existing_ids['HMDB_ID']).strip()
                hmdb_matches = self.mzcloud_cache[self.mzcloud_cache['HMDB ID'].astype(str).str.strip() == hmdb_id]
                if not hmdb_matches.empty:
                    matches = hmdb_matches
                    logger.info(f"[OK] MZCloud HMDB ID match for: {name}")
            
            # Strategy 4: PubChem CID match (if available)
            if matches is None and existing_ids.get('PubChem_CID') and 'PubChem CID' in self.mzcloud_cache.columns:
                pubchem_cid = str(existing_ids['PubChem_CID']).strip()
                cid_matches = self.mzcloud_cache[self.mzcloud_cache['PubChem CID'].astype(str).str.strip() == pubchem_cid]
                if not cid_matches.empty:
                    matches = cid_matches
                    logger.info(f"[OK] MZCloud PubChem CID match for: {name}")
            
            if matches is not None and not matches.empty:
                match_row = matches.iloc[0]  # Take first match
                result['found'] = True
                
                # Set Endogenous to 'Yes' if found in MZCloud (since it's an endogenous metabolite database)
                result['Endogenous'] = 'Yes'
                
                # Add any missing IDs from MZCloud
                additional_ids = {}
                id_columns = ['KEGG ID', 'CAS Number', 'SMILES', 'InChI', 'InChIKey', 'HMDB ID', 'PubChem CID', 'Formula']
                
                for col in id_columns:
                    if col in self.mzcloud_cache.columns:
                        value = str(match_row[col]).strip() if pd.notna(match_row[col]) else ''
                        if value and value != 'nan':
                            # Map column names to expected result keys
                            key_mapping = {
                                'KEGG ID': 'KEGG_ID',
                                'CAS Number': 'CAS',
                                'HMDB ID': 'HMDB_ID',
                                'PubChem CID': 'PubChem_CID',
                                'Formula': 'Molecular_Formula'
                            }
                            key = key_mapping.get(col, col)
                            additional_ids[key] = value
                
                result['additional_ids'] = additional_ids
                logger.info(f"[OK] MZCloud found endogenous metabolite: {name} with {len(additional_ids)} additional IDs")
            else:
                logger.info(f"[MISS] MZCloud no matches for: {name}")

        except Exception as e:
            logger.error(f"[ERROR] MZCloud cache search exception for {name}: {e}")

        return result

    def annotate_metabolite_ids(self, name: str, expected_formula: Optional[str] = None) -> Dict[str, Any]:
        """
        Complete ID annotation for a single metabolite
        Priority: 
        - KEGG API (highest priority for KEGG IDs - most trusted source)
        - PubChem API (for structure, other IDs, and KEGG fallback)
        - HMDB offline (for biological classification)
        - MZCloud cache (for endogenous status)
        """
    # logger.info(f"[START] Annotating metabolite: {name}")  # noisy per-metabolite log (commented out)
        
        # Initialize result with metabolite name
        result = {
            'Name': name,
            'LipidMaps_ID': '',
            'PubChem_CID': '',
            'KEGG_ID': '',
            'KEGG_Match_Type': '',
            'HMDB_ID': '',
            'ChEBI_ID': '',
            'CAS': '',
            'SMILES': '',
            'InChI': '',
            'InChIKey': '',
            'IUPAC_Name': '',
            'Molecular_Formula': '',
            'Molecular_Weight': '',
            'Endogenous': '',
            'Super_Class': '',
            'Class': '',
            'Sub_Class': '',
            'annotation_sources': []
        }
        
        # Step 1: KEGG API search (HIGHEST PRIORITY for KEGG IDs - most trusted source)
        kegg_result = self.search_kegg_api(name)
        if kegg_result['found']:
            result['KEGG_ID'] = kegg_result['KEGG_ID']
            result['KEGG_Match_Type'] = 'Name'
            # Use KEGG formula if available
            if kegg_result['Molecular_Formula']:
                result['Molecular_Formula'] = kegg_result['Molecular_Formula']
            result['annotation_sources'].append('KEGG_API')
        
        # Step 2: PubChem API search (for structure and other cross-references)
        pubchem_result = self.search_pubchem_api(name, expected_formula)
        if pubchem_result.get('found'):
            result.update({
                'PubChem_CID': pubchem_result.get('PubChem_CID', ''),
                'SMILES': pubchem_result.get('SMILES', ''),
                'InChI': pubchem_result.get('InChI', ''),
                'InChIKey': pubchem_result.get('InChIKey', ''),
                'IUPAC_Name': pubchem_result.get('IUPAC_Name', ''),
                'Molecular_Weight': pubchem_result.get('Molecular_Weight', ''),
                'ChEBI_ID': pubchem_result.get('ChEBI_ID', ''),
                'CAS': pubchem_result.get('CAS', ''),
                'LipidMaps_ID': pubchem_result.get('LipidMaps_ID', ''),
            })
            # Use PubChem formula only if not already set by KEGG
            if not result['Molecular_Formula'] and pubchem_result.get('Molecular_Formula'):
                result['Molecular_Formula'] = pubchem_result['Molecular_Formula']
            
            # Use PubChem's KEGG cross-reference ONLY if KEGG API didn't find it
            if not result['KEGG_ID'] and pubchem_result.get('KEGG_ID'):
                result['KEGG_ID'] = pubchem_result['KEGG_ID']
            
            # Use PubChem's other cross-references if not already set
            if pubchem_result['HMDB_ID'] and not result['HMDB_ID']:
                result['HMDB_ID'] = pubchem_result['HMDB_ID']
            if pubchem_result['LipidMaps_ID'] and not result['LipidMaps_ID']:
                result['LipidMaps_ID'] = pubchem_result['LipidMaps_ID']
            result['annotation_sources'].append('PubChem_API')

            # If PubChem enrichment is enabled, fetch additional cross-references
            if self.pubchem_enrichment and result.get('PubChem_CID'):
                try:
                    cid_int = int(float(str(result.get('PubChem_CID'))))
                    xrefs = fetch_pubchem_cross_references(cid_int)
                    
                    # Add LipidMaps if missing
                    if not result.get('LipidMaps_ID') and xrefs.get('LipidMaps'):
                        result['LipidMaps_ID'] = xrefs.get('LipidMaps')
                        logger.info(f"[PUBCHEM_XREF] Added LipidMaps_ID from PubChem for {name}: {result['LipidMaps_ID']}")
                    
                    # Add KEGG ONLY if KEGG API didn't find it (fallback only)
                    if not result.get('KEGG_ID') and xrefs.get('KEGG'):
                        result['KEGG_ID'] = xrefs.get('KEGG')
                        logger.info(f"[PUBCHEM_XREF] Added KEGG_ID from PubChem xrefs (fallback) for {name}: {result['KEGG_ID']}")
                    
                    # Add other IDs if missing
                    if not result.get('HMDB_ID') and xrefs.get('HMDB'):
                        result['HMDB_ID'] = xrefs.get('HMDB')
                    if not result.get('ChEBI_ID') and xrefs.get('ChEBI'):
                        result['ChEBI_ID'] = xrefs.get('ChEBI')
                    if not result.get('CAS') and xrefs.get('CAS'):
                        result['CAS'] = xrefs.get('CAS')
                except Exception as e:
                    logger.debug(f"[PUBCHEM_XREF] Failed to enrich PubChem cross-refs for {name}: {e}")

        # Step 3: HMDB offline database search (if HMDB ID not found via APIs)
        if not result['HMDB_ID']:
            hmdb_result = self.search_hmdb_offline(name, result)
            if hmdb_result['found']:
                result['HMDB_ID'] = hmdb_result['HMDB_ID']
                # Fill in missing fields from HMDB
                if not result['Molecular_Formula'] and hmdb_result['Molecular_Formula']:
                    result['Molecular_Formula'] = hmdb_result['Molecular_Formula']
                if not result['InChIKey'] and hmdb_result['InChIKey']:
                    result['InChIKey'] = hmdb_result['InChIKey']
                
                # Add biological information from HMDB
                result['Endogenous'] = hmdb_result.get('Endogenous', '')
                result['Super_Class'] = hmdb_result.get('Super_Class', '')
                result['Class'] = hmdb_result.get('Class', '')
                result['Sub_Class'] = hmdb_result.get('Sub_Class', '')
                    
                result['annotation_sources'].append('HMDB_Offline')
        
        # Alternative: Get biological info from HMDB even if we already have HMDB_ID from PubChem
        elif result['HMDB_ID'] and not result.get('Endogenous'):
            # Try to enhance existing HMDB_ID with biological information
            hmdb_result = self.search_hmdb_offline(name, result)
            if hmdb_result['found']:
                endogenous_value = hmdb_result.get('Endogenous', '')
                result['Endogenous'] = endogenous_value
                result['Super_Class'] = hmdb_result.get('Super_Class', '')
                result['Class'] = hmdb_result.get('Class', '')
                result['Sub_Class'] = hmdb_result.get('Sub_Class', '')
                # logger.info(f"[ENHANCED] Added biological data for {name} from HMDB - Endogenous: '{endogenous_value}'")  # noisy (commented out)

        # Step 4: MZCloud cache search (for endogenous status and additional IDs)
    # logger.info(f"[CACHE] Searching MZCloud cache for: {name} (Current Endogenous: '{result.get('Endogenous', '')}')")  # noisy (commented out)
        mzcloud_result = self.search_mzcloud_cache(name, result)
        if mzcloud_result['found']:
            # Update endogenous status if found in MZCloud
            if mzcloud_result.get('Endogenous'):
                result['Endogenous'] = mzcloud_result['Endogenous']
                # logger.info(f"[CACHE SUCCESS] MZCloud: Set {name} as Endogenous: {mzcloud_result['Endogenous']}")  # noisy (commented out)
            
            # Add any missing IDs from MZCloud
            added_ids = []
            for key, value in mzcloud_result.get('additional_ids', {}).items():
                if not result.get(key) and value:  # Only update if current value is empty
                    result[key] = value
                    added_ids.append(f"{key}={value}")
            
            #if added_ids:
                # logger.info(f"[CACHE SUCCESS] MZCloud: Added IDs for {name}: {', '.join(added_ids)}")  # noisy (commented out)
            
            result['annotation_sources'].append('MZCloud_Cache')
            # logger.info(f"[ENHANCED] Added endogenous status and IDs for {name} from MZCloud")  # noisy (commented out)
        #else:
            # logger.info(f"[CACHE MISS] MZCloud: No match found for {name}")  # noisy (commented out)

        # Summary
        sources = ', '.join(result['annotation_sources']) if result['annotation_sources'] else 'None'
        ids_found = sum(1 for key in ['PubChem_CID', 'KEGG_ID', 'HMDB_ID', 'ChEBI_ID'] if result[key])
    # logger.info(f"[COMPLETE] {name}: {ids_found} IDs found from sources: {sources} - Final Endogenous: '{result.get('Endogenous', '')}'")  # noisy (commented out)

        return result

    def process_metabolites(self, input_df: pd.DataFrame, max_workers: int = 6) -> pd.DataFrame:
        """
        Process all metabolites for ID annotation
        """
        # Determine optional formula column
        formula_col = None
        for cand in ['Molecular_Formula', 'Molecular Formula', 'Formula']:
            if cand in input_df.columns:
                formula_col = cand
                break

        records = input_df[['Name'] + ([formula_col] if formula_col else [])].to_dict('records')
        logger.info(f"Starting ID annotation for {len(records)} metabolites using {max_workers} workers")

        results = []
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_metabolite = {}
            for rec in records:
                name = rec.get('Name')
                if not isinstance(name, str) or not name.strip():
                    continue
                expected_formula = rec.get(formula_col) if formula_col else None
                fut = executor.submit(self.annotate_metabolite_ids, str(name), expected_formula)
                future_to_metabolite[fut] = str(name)
            
            # Process completed tasks
            last_logged_progress = 0
            for i, future in enumerate(as_completed(future_to_metabolite), 1):
                # Check if user requested stop
                if self.stop_check_callback():
                    logger.info("Stop requested by user - cancelling remaining tasks")
                    # Cancel pending futures
                    for fut in future_to_metabolite.keys():
                        fut.cancel()
                    break
                
                metabolite = future_to_metabolite[future]
                try:
                    result = future.result()
                    results.append(result)
                    
                    # Progress reporting
                    progress_percent = (i / len(records)) * 100
                    
                    # Only log progress to console every 20% to reduce clutter
                    if progress_percent - last_logged_progress >= 20 or i == len(records):
                        progress_msg = f"Progress: {i}/{len(records)} ({progress_percent:.1f}%)"
                        logger.info(progress_msg)
                        last_logged_progress = int(progress_percent // 20) * 20
                except Exception as e:
                    logger.error(f"Error processing metabolite {metabolite}: {e}")
        
        # Convert results to DataFrame
        final_df = pd.DataFrame(results) if results else pd.DataFrame()
        
        # Debug: Check Endogenous values after merge and fill
        if 'Endogenous' in final_df.columns:
            endogenous_after = final_df['Endogenous'].tolist()
            #logger.info(f"[DEBUG] Endogenous values after merge: {endogenous_after}")
            
            # Fix any NaN values in Endogenous column
            import numpy as np
            final_df['Endogenous'] = final_df['Endogenous'].replace(np.nan, '')
            final_df['Endogenous'] = final_df['Endogenous'].astype(str).replace('nan', '')
            
            endogenous_fixed = final_df['Endogenous'].tolist()
            #logger.info(f"[DEBUG] Endogenous values after NaN fix: {endogenous_fixed}")
            
            # Post-merge fallback: if Endogenous still empty but HMDB_ID exists,
            # map HMDB_ID -> Is_Endogenous from loaded HMDB dataframe.
            try:
                hmdb_map = {}
                if hasattr(self, 'hmdb_df') and 'HMDB' in self.hmdb_df.columns and 'Is_Endogenous' in self.hmdb_df.columns:
                    # build a fast lookup dict
                    temp = self.hmdb_df[['HMDB', 'Is_Endogenous']].dropna()
                    for idx, r in temp.iterrows():
                        key = str(r['HMDB']).strip().upper()
                        if key:
                            hmdb_map[key] = str(r['Is_Endogenous']).strip()

                if hmdb_map:
                    filled = 0
                    for i, row in final_df.iterrows():
                        cur = str(row.get('Endogenous', '')).strip()
                        if not cur:
                            hmdb_id = str(row.get('HMDB_ID', '')).strip().upper()
                            if hmdb_id and hmdb_id in hmdb_map and hmdb_map[hmdb_id]:
                                final_df.at[i, 'Endogenous'] = hmdb_map[hmdb_id]
                                filled += 1
                    #logger.info(f"[FALLBACK] Filled Endogenous from HMDB map for {filled} metabolites")
            except Exception as e:
                logger.error(f"[FALLBACK] Error filling Endogenous from HMDB map: {e}")

        logger.info(f"ID annotation completed for {len(final_df)} metabolites")
        return final_df

    def run_id_annotation(self, max_workers: int = 6):
        """
        Main method to run complete ID annotation process
        """
        try:
            # If lipid mode: route to offline lipid workflow (includes PubChem first, then KEGG/HMDB)
            if getattr(self, 'lipid_mode', False):
                logger.info("[MODE] Running Lipid ID mode (offline LMSD + PubChem first, then KEGG/HMDB)")
                return self.run_lipid_id_annotation(max_workers=max_workers)

            logger.info("\n" + "="*80)
            logger.info("🔍 DEBUG: DETERMINING INPUT SOURCE FOR ID ANNOTATION")
            logger.info("="*80)

            if self.input_df_override is not None and isinstance(self.input_df_override, pd.DataFrame) and not self.input_df_override.empty:
                logger.info("✅ USING USER-ASSIGNED INPUT DATAFRAME FROM COLUMN ASSIGNMENT DIALOG")
                input_df = self.input_df_override.copy()
                logger.info(f"User-assigned input rows: {len(input_df)}")
                logger.info(f"Columns in user-assigned input: {list(input_df.columns)}")
            else:
                # Priority 1: Use cleaned_metabolites_df (Combined sheet) if available
                cleaned_data_in_memory = self._get_cleaned_data_from_memory()
                
                logger.info(f"Checking cleaned_metabolites_df: {cleaned_data_in_memory is not None}")
                if cleaned_data_in_memory is not None:
                    logger.info(f"  - Shape: {cleaned_data_in_memory.shape}")
                    logger.info(f"  - Columns: {list(cleaned_data_in_memory.columns)[:5]}...")
                
                logger.info(f"Checking pos_df: {self.pos_df is not None}")
                if self.pos_df is not None and isinstance(self.pos_df, pd.DataFrame):
                    logger.info(f"  - Shape: {self.pos_df.shape}")
                    logger.info(f"  - Columns: {list(self.pos_df.columns)[:5]}...")
                
                logger.info(f"force_file_input: {self.force_file_input}")
                logger.info("="*80 + "\n")
                
                if (cleaned_data_in_memory is not None and not cleaned_data_in_memory.empty and 
                    not self.force_file_input):
                    # Use cleaned data from memory as primary source (Combined sheet with ALL metabolites)
                    logger.info("="*60)
                    logger.info("✅ USING COMBINED SHEET FROM MEMORY AS PRIMARY INPUT")
                    logger.info("="*60)
                    logger.info(f"Using Combined sheet from memory: {len(cleaned_data_in_memory)} rows")
                    input_df = cleaned_data_in_memory.copy()
                    logger.info(f"Columns in Combined sheet: {list(input_df.columns)}")
                # Priority 2: Use pos_df if available (custom mode with pre-processed/renamed columns)
                elif self.pos_df is not None and isinstance(self.pos_df, pd.DataFrame) and not self.pos_df.empty:
                    logger.info("⚠️ [METABOLITE MODE] Using pre-loaded dataframe (Positive sheet - fallback)")
                    input_df = self.pos_df.copy()
                    logger.info(f"[METABOLITE MODE] Columns in pre-loaded dataframe: {list(input_df.columns)}")
                else:
                    # Priority 3: Use file input when no memory data available or when forced
                    if self.force_file_input:
                        logger.info(f"Force file input enabled - loading input file: {self.input_file}")
                    else:
                        logger.info(f"No cleaned data in memory - loading input file: {self.input_file}")
                    input_df = pd.read_excel(self.input_file)
            
            # Validate required columns
            if 'Name' not in input_df.columns:
                if 'Metabolite' in input_df.columns:
                    input_df = input_df.rename(columns={'Metabolite': 'Name'})
                elif 'Molecule' in input_df.columns:
                    input_df = input_df.rename(columns={'Molecule': 'Name'})
                elif 'LipidID' in input_df.columns:
                    input_df = input_df.rename(columns={'LipidID': 'Name'})
                elif 'metabolite_id' in input_df.columns:
                    input_df = input_df.rename(columns={'metabolite_id': 'Name'})
                else:
                    raise ValueError("Input file must contain a column named 'Name', 'Metabolite', 'Molecule', or 'LipidID'")
            
            # Remove rows with missing names
            input_df = input_df.dropna(subset=['Name'])
            logger.info(f"Processing {len(input_df)} metabolites for ID annotation")

            # Process metabolites
            annotated_df = self.process_metabolites(input_df, max_workers)

            # Drop metabolites with no IDs at all (core ID columns all empty)
            # SKIP this filter if skip_id_filtering is True (custom mode)
            if not self.skip_id_filtering:
                try:
                    id_core_cols = ['LipidMaps_ID','PubChem_CID','KEGG_ID','HMDB_ID','ChEBI_ID','CAS','SMILES','InChI','InChIKey']
                    existing_id_cols = [c for c in id_core_cols if c in annotated_df.columns]
                    if existing_id_cols:
                        def _has_any_id(row):
                            for c in existing_id_cols:
                                val = row.get(c, None)
                                if pd.notna(val) and str(val).strip() != '' and str(val).strip().lower() not in ('none','nan'):
                                    return True
                            return False
                        before_count = len(annotated_df)
                        annotated_df = annotated_df[annotated_df.apply(_has_any_id, axis=1)].reset_index(drop=True)
                        removed = before_count - len(annotated_df)
                        logger.info(f"Filtered out {removed} metabolites with no IDs across {len(existing_id_cols)} ID columns")
                except Exception as e:
                    logger.warning(f"Failed to filter empty-ID metabolites: {e}")
            else:
                logger.info("Skip ID filtering enabled - keeping all metabolites regardless of ID status")
            
            # Clean up annotation_sources format for better readability
            if 'annotation_sources' in annotated_df.columns:
                def format_sources(sources_list, row):
                    """Format sources based on what IDs were actually found"""
                    if isinstance(sources_list, list):
                        clean_sources = []
                        for source in sources_list:
                            # Only include source if we actually got meaningful data from it
                            if 'LipidMaps' in source and row.get('LipidMaps_ID', ''):
                                clean_sources.append('LipidMaps')
                            elif 'PubChem' in source and row.get('PubChem_CID', ''):
                                clean_sources.append('PubChem')
                            elif 'KEGG' in source and row.get('KEGG_ID', ''):
                                clean_sources.append('KEGG')
                            elif 'HMDB' in source and row.get('HMDB_ID', ''):
                                clean_sources.append('HMDB')
                            elif 'MZCloud' in source and (row.get('Endogenous', '') == 'Yes' or 
                                                        row.get('PubChem_CID', '') or 
                                                        row.get('HMDB_ID', '') or 
                                                        row.get('KEGG_ID', '') or
                                                        row.get('InChIKey', '')):
                                clean_sources.append('MZCloud')
                            elif 'ChEBI' in source and row.get('ChEBI_ID', ''):
                                # ChEBI IDs usually come from PubChem cross-references
                                if 'PubChem' not in clean_sources:
                                    clean_sources.append('PubChem')
                        return ', '.join(sorted(set(clean_sources))) if clean_sources else ''
                    return str(sources_list) if sources_list else ''
                
                annotated_df['annotation_sources'] = annotated_df.apply(
                    lambda row: format_sources(row['annotation_sources'], row), axis=1)

            # Try to merge with metabolite IDs data (from MzCloud/Metabolika files)
            try:
                if self.metabolite_ids_df is not None and not self.metabolite_ids_df.empty:
                    logger.info("="*60)
                    logger.info("MERGING WITH METABOLITE IDS DATA FROM MEMORY")
                    logger.info("="*60)
                    annotated_df = self._merge_with_cleaned_data(annotated_df, self.metabolite_ids_df)
                    # After merging, again drop rows with no IDs in case merge added empty rows
                    # SKIP this filter if skip_id_filtering is True (custom mode)
                    if not self.skip_id_filtering:
                        try:
                            id_core_cols2 = ['LipidMaps_ID','PubChem_CID','KEGG_ID','HMDB_ID','ChEBI_ID','CAS','SMILES','InChI','InChIKey']
                            existing2 = [c for c in id_core_cols2 if c in annotated_df.columns]
                            if existing2:
                                def _has_any_id2(row):
                                    for c in existing2:
                                        val = row.get(c, None)
                                        if pd.notna(val) and str(val).strip() != '' and str(val).strip().lower() not in ('none','nan'):
                                            return True
                                    return False
                                before2 = len(annotated_df)
                                annotated_df = annotated_df[annotated_df.apply(_has_any_id2, axis=1)].reset_index(drop=True)
                                rem2 = before2 - len(annotated_df)
                                if rem2:
                                    logger.info(f"Post-merge filtering removed {rem2} no-ID metabolites")
                        except Exception as e:
                            logger.warning(f"Post-merge empty-ID filter failed: {e}")
                else:
                    logger.info("No metabolite IDs data available - skipping merge step")
            except Exception as e:
                logger.warning(f"Failed to merge with cleaned data: {e}")
                logger.info("Continuing with original annotated data")

            # If polarity DataFrames are available, assign IDs and save multi-sheet Excel
            posneg_results = None
            try:
                posneg_results = self._assign_ids_to_polarity_frames(annotated_df)
            except Exception as e:
                logger.warning(f"Failed to assign IDs to pos/neg DataFrames: {e}")

            # Add MSI-level and confidence columns to merged output.
            annotated_df = add_msi_confidence_columns(annotated_df)
            before_conf = len(annotated_df)
            annotated_df = apply_confidence_filter(annotated_df, self.confidence_filter_mode)
            if self.confidence_filter_mode != 'none':
                logger.info(f"Confidence filter ({self.confidence_filter_mode}): kept {len(annotated_df)}/{before_conf} rows")

            # Save results (multi-sheet when applicable)
            logger.info(f"Saving annotated results to: {self.output_file}")
            final_merged_df = annotated_df
            try:
                # SMILES is used internally for some lookup paths, but it should not be
                # written to the final exported sheets.
                export_drop_cols = ['SMILES']
                annotated_df = annotated_df.drop(columns=[c for c in export_drop_cols if c in annotated_df.columns], errors='ignore')

                # Prepare polarity outputs with robust fallbacks so the sheets always exist
                pos_df_out: pd.DataFrame = pd.DataFrame()
                neg_df_out: pd.DataFrame = pd.DataFrame()

                if posneg_results and isinstance(posneg_results, dict):
                    maybe_pos = posneg_results.get('pos_id_df')
                    maybe_neg = posneg_results.get('neg_id_df')
                    if isinstance(maybe_pos, pd.DataFrame):
                        pos_df_out = maybe_pos.copy()
                    if isinstance(maybe_neg, pd.DataFrame):
                        neg_df_out = maybe_neg.copy()

                # Fallback: split annotated_df by Polarity when needed
                if (pos_df_out is None or pos_df_out.empty) and 'Polarity' in annotated_df.columns:
                    try:
                        pol = annotated_df['Polarity'].astype(str).str.strip().str.lower()
                        pos_df_out = annotated_df[pol.str.startswith('pos')].copy()
                    except Exception:
                        pos_df_out = pd.DataFrame()
                if (neg_df_out is None or neg_df_out.empty) and 'Polarity' in annotated_df.columns:
                    try:
                        pol = annotated_df['Polarity'].astype(str).str.strip().str.lower()
                        neg_df_out = annotated_df[pol.str.startswith('neg')].copy()
                    except Exception:
                        neg_df_out = pd.DataFrame()

                if isinstance(pos_df_out, pd.DataFrame) and 'SMILES' in pos_df_out.columns:
                    pos_df_out = pos_df_out.drop(columns=['SMILES'])
                if isinstance(neg_df_out, pd.DataFrame) and 'SMILES' in neg_df_out.columns:
                    neg_df_out = neg_df_out.drop(columns=['SMILES'])

                # Rebuild the merged sheet from Pos/Neg outputs when they are available.
                # This preserves the existing Pos/Neg priority workflow and avoids using
                # the Combined annotation dataframe as the final merged export.
                merged_df_out = pd.DataFrame()
                merged_sources = []
                if isinstance(pos_df_out, pd.DataFrame) and not pos_df_out.empty:
                    merged_sources.append(pos_df_out.copy())
                if isinstance(neg_df_out, pd.DataFrame) and not neg_df_out.empty:
                    merged_sources.append(neg_df_out.copy())

                if merged_sources:
                    merged_df_out = pd.concat(merged_sources, ignore_index=True)
                    if 'LipidID' in merged_df_out.columns:
                        merged_df_out = merged_df_out.drop_duplicates(subset=['LipidID'], keep='first')
                    elif 'Name' in merged_df_out.columns:
                        merged_df_out = merged_df_out.drop_duplicates(subset=['Name'], keep='first')

                    merged_df_out = add_msi_confidence_columns(merged_df_out)
                    merged_df_out = apply_confidence_filter(merged_df_out, self.confidence_filter_mode)
                else:
                    merged_df_out = add_msi_confidence_columns(annotated_df.copy())
                    merged_df_out = apply_confidence_filter(merged_df_out, self.confidence_filter_mode)

                merged_df_out = _enforce_endogenous_rules(merged_df_out)

                final_merged_df = merged_df_out

                # Ensure MSI columns exist in all output sheets.
                pos_df_out = add_msi_confidence_columns(pos_df_out)
                neg_df_out = add_msi_confidence_columns(neg_df_out)
                pos_df_out = apply_confidence_filter(pos_df_out, self.confidence_filter_mode)
                neg_df_out = apply_confidence_filter(neg_df_out, self.confidence_filter_mode)
                pos_df_out = _enforce_endogenous_rules(pos_df_out)
                neg_df_out = _enforce_endogenous_rules(neg_df_out)

                # Expose for GUI/memory consumers
                try:
                    if isinstance(pos_df_out, pd.DataFrame):
                        self.pos_id_df = pos_df_out
                    if isinstance(neg_df_out, pd.DataFrame):
                        self.neg_id_df = neg_df_out
                except Exception:
                    pass

                with pd.ExcelWriter(self.output_file, engine='openpyxl') as writer:
                    merged_df_out.to_excel(writer, index=False, sheet_name='Merged_IDs')
                    # Always create Pos_id and Neg_id sheets (even if empty) so downstream steps find them
                    (pos_df_out if isinstance(pos_df_out, pd.DataFrame) else pd.DataFrame()).to_excel(
                        writer, index=False, sheet_name='Pos_id'
                    )
                    (neg_df_out if isinstance(neg_df_out, pd.DataFrame) else pd.DataFrame()).to_excel(
                        writer, index=False, sheet_name='Neg_id'
                    )

                # QA log: Endogenous Yes/No distribution per output sheet.
                def _log_endogenous_counts(df_obj: Optional[pd.DataFrame], sheet_name: str) -> None:
                    if not isinstance(df_obj, pd.DataFrame):
                        logger.info(f"[QA] {sheet_name} Endogenous: sheet unavailable")
                        return
                    if df_obj.empty:
                        logger.info(f"[QA] {sheet_name} Endogenous: empty sheet")
                        return
                    if 'Endogenous' not in df_obj.columns:
                        logger.info(f"[QA] {sheet_name} Endogenous: column missing")
                        return

                    normalized = df_obj['Endogenous'].apply(_normalize_yes_no)
                    yes_count = int((normalized == 'Yes').sum())
                    no_count = int((normalized == 'No').sum())
                    total_count = int(len(normalized))
                    logger.info(f"[QA] {sheet_name} Endogenous -> Yes: {yes_count}, No: {no_count}, Total: {total_count}")

                _log_endogenous_counts(merged_df_out, 'Merged_IDs')
                _log_endogenous_counts(pos_df_out, 'Pos_id')
                _log_endogenous_counts(neg_df_out, 'Neg_id')
            except Exception as e:
                # Fallback to single-sheet save if multi-sheet failed
                logger.warning(f"Multi-sheet save failed ({e}); falling back to single sheet.")
                annotated_df.to_excel(self.output_file, index=False)

            # Save caches
            self._save_cache(self.pubchem_cache, "pubchem_cache.pkl")
            self._save_cache(self.kegg_cache, "kegg_cache.pkl")
            self._save_cache(self.hmdb_cache, "hmdb_cache.pkl")
            self._save_cache(self.lipidmaps_cache, "lipidmaps_cache.pkl")

            # Print summary statistics
            total_metabolites = len(final_merged_df) if isinstance(final_merged_df, pd.DataFrame) else len(annotated_df)
            
            # Count database hits
            lipidmaps_specific = 0
            lipidmaps_total = 0
            stats_df = final_merged_df if isinstance(final_merged_df, pd.DataFrame) else annotated_df

            pubchem_hits = stats_df['PubChem_CID'].astype(str).str.strip().ne('').sum() if 'PubChem_CID' in stats_df.columns else 0
            kegg_hits = stats_df['KEGG_ID'].astype(str).str.strip().ne('').sum() if 'KEGG_ID' in stats_df.columns else 0
            hmdb_hits = stats_df['HMDB_ID'].astype(str).str.strip().ne('').sum() if 'HMDB_ID' in stats_df.columns else 0
            chebi_hits = stats_df['ChEBI_ID'].astype(str).str.strip().ne('').sum() if 'ChEBI_ID' in stats_df.columns else 0

            logger.info("="*60)
            logger.info("ID ANNOTATION SUMMARY")
            logger.info("="*60)
            logger.info(f"Total metabolites processed: {total_metabolites}")
            logger.info(f"PubChem CID found: {pubchem_hits} ({(pubchem_hits/total_metabolites*100 if total_metabolites else 0):.1f}%)")
            logger.info(f"KEGG ID found: {kegg_hits} ({(kegg_hits/total_metabolites*100 if total_metabolites else 0):.1f}%)")
            logger.info(f"HMDB ID found: {hmdb_hits} ({(hmdb_hits/total_metabolites*100 if total_metabolites else 0):.1f}%)")
            logger.info(f"ChEBI ID found: {chebi_hits} ({(chebi_hits/total_metabolites*100 if total_metabolites else 0):.1f}%)")
            logger.info("="*60)

            # Expose polarity ID DataFrames (if produced) for external callers (GUI memory_store)
            try:
                if posneg_results and isinstance(posneg_results, dict):
                    self.posneg_results = posneg_results
                    pos_df_out = posneg_results.get('pos_id_df') if isinstance(posneg_results, dict) else None
                    neg_df_out = posneg_results.get('neg_id_df') if isinstance(posneg_results, dict) else None
                    if isinstance(pos_df_out, pd.DataFrame):
                        self.pos_id_df = pos_df_out
                    if isinstance(neg_df_out, pd.DataFrame):
                        self.neg_id_df = neg_df_out
            except Exception as _expose_err:
                logger.warning(f"Could not expose polarity ID DataFrames: {_expose_err}")

            return annotated_df

        except Exception as e:
            logger.error(f"Error in ID annotation process: {e}")
            raise

    # -----------------------------------------
    # Post-merge: assign IDs to pos/neg frames
    # -----------------------------------------
    def _assign_ids_to_polarity_frames(self, merged_ids_df: pd.DataFrame) -> Dict[str, Optional[pd.DataFrame]]:
        """Assign ID columns from merged IDs to pos/neg DataFrames in memory by Name.

        - Creates Name_Key (lowercased Name) on both sides
        - Merges selected ID columns immediately after feature columns
        - Applies optional filters according to settings
        
        Returns dict with keys 'pos_id_df' and 'neg_id_df' (if created).
        """
        try:
            results: Dict[str, Optional[pd.DataFrame]] = {'pos_id_df': None, 'neg_id_df': None}

            if merged_ids_df is None or merged_ids_df.empty:
                return results

            # Build Name_Key on merged ids
            ids_df = merged_ids_df.copy()
            if 'Name' not in ids_df.columns:
                return results
            ids_df['Name_Key'] = ids_df['Name'].astype(str).str.strip().str.lower()

            # Define ID columns and feature columns
            id_columns_all = [
                'LipidMaps_ID', 'PubChem_CID', 'KEGG_ID', 'HMDB_ID', 'ChEBI_ID', 'CAS', 'Endogenous',
                'SMILES', 'InChI', 'InChIKey', 'IUPAC_Name', 'Super_Class', 'Class', 'Sub_Class', 'Endogenous_Source'
            ]
            # Only keep those present in merged IDs
            id_cols_present = [c for c in id_columns_all if c in ids_df.columns]

            feature_cols_expected = [
                'Name', 'Formula', 'Annot. DeltaMass [ppm]', 'Calc. MW', 'm/z',
                'Reference Ion', 'RT [min]', 'Area (Max.)', 'Metabolika Pathways',
                'BioCyc Pathways', 'Polarity', 'MS2', 'MS2 Purity [%]'
            ]
            # Also consider renamed stage-2 names
            alt_feature_names = {'RT [min]': 'RT', 'Annot. DeltaMass [ppm]': 'ppm', 'Calc. MW': 'MW', 'MS2 Purity [%]': 'MS2_Purity'}

            # --- Helpers for duplicate detection based on Formula + IDs ---
            def _normalize_formula(val: Any) -> str:
                try:
                    s = str(val)
                    # remove spaces and standardize casing
                    return s.replace(" ", "").strip()
                except Exception:
                    return str(val)

            def _get_dup_key(row: pd.Series) -> Optional[str]:
                """Build a grouping key using Formula + any available ID in priority order.
                Priority: HMDB_ID > KEGG_ID > LipidMaps_ID > PubChem_CID.
                Returns None if no IDs are available for this row.
                """
                id_cols_priority = ['HMDB_ID', 'KEGG_ID', 'LipidMaps_ID', 'PubChem_CID']
                formula = _normalize_formula(row.get('Formula', ''))
                if not formula:
                    return None
                for col in id_cols_priority:
                    if col in row.index:
                        v = row.get(col)
                        if pd.notna(v) and str(v).strip() not in ('', 'nan', 'None'):
                            return f"{formula}|{col}|{str(v).strip()}"
                return None

            def _ion_rank(ion: Any, polarity: str) -> int:
                """Rank Reference Ion using priority lists from data cleaning (lower is better)."""
                pos_order = ['[M+H]+1', '[M+2H]+2', '[M+H-H2O]+1', '[M+H+MeOH]+1', '[M+FA+H]+1']
                neg_order = ['[M-H]-1', '[M-2H]-2', '[M-H-H2O]-1', '[M-H-MeOH]-1', '[M+FA-H]-1']
                order = pos_order if polarity == 'positive' else neg_order
                try:
                    ion_norm = str(ion).strip()
                    return order.index(ion_norm) if ion_norm in order else 999
                except Exception:
                    return 999

            def _ms2_rank(ms2: Any) -> int:
                """Rank MS2 quality: preferred < other < anything else (including No MS2)."""
                s = str(ms2).strip().lower()
                if 'dda for preferred ion' in s:
                    return 0
                if 'dda for other ion' in s or 'dda for non-preferred ion' in s:
                    return 1
                if 'no ms2' in s:
                    return 5
                return 9

            def _dedup_by_formula_ids_sum(df_in: pd.DataFrame, polarity: str) -> pd.DataFrame:
                """Collapse duplicates where Formula + any ID matches.
                - Only rows within ±{self.id_dedup_rt_window_minutes} minutes of group mean RT are combined (configurable via GUI).
                - Numeric columns (including Area (Max.)) are summed; RT averaged.
                - Representative non-numeric metadata chosen by priority:
                  Reference Ion (ion order) → MS2 rank → Area (Max.) desc.
                Rows outside the RT window are kept as-is.
                """
                if df_in is None or df_in.empty:
                    return df_in

                df = df_in.copy()

                # Ensure needed columns exist gracefully
                if 'Formula' not in df.columns:
                    return df

                # Build duplicate key
                df['_dup_key'] = df.apply(_get_dup_key, axis=1)
                if df['_dup_key'].isna().all():
                    # No usable dup keys
                    df = df.drop(columns=['_dup_key'])
                    return df

                # Identify numeric columns to sum (exclude feature columns)
                feature_like = {'Name', 'Name_Key', 'Formula', 'MW', 'ppm', 'Reference Ion', 'MS2', 'm/z', 'RT',
                                'Polarity', 'MS2_Purity', 'Metabolika Pathways', 'BioCyc Pathways',
                                'LipidMaps_ID', 'PubChem_CID', 'KEGG_ID', 'HMDB_ID', 'ChEBI_ID', 'CAS',
                                'SMILES', 'InChI', 'InChIKey', 'IUPAC_Name', 'Super_Class', 'Class', 'Sub_Class',
                                'Endogenous_Source'}
                numeric_cols: list[str] = []
                for col in df.columns:
                    if col == '_dup_key':
                        continue
                    if col in feature_like:
                        continue
                    if pd.api.types.is_numeric_dtype(df[col]):
                        numeric_cols.append(col)
                    else:
                        # Try coercion to detect numeric
                        coerced = pd.to_numeric(df[col], errors='coerce')
                        if coerced.notna().any():
                            df[col] = coerced
                            numeric_cols.append(col)

                # Ensure Area (Max.) counts as numeric to sum
                if 'Area (Max.)' in df.columns and 'Area (Max.)' not in numeric_cols:
                    df['Area (Max.)'] = pd.to_numeric(df['Area (Max.)'], errors='coerce')
                    numeric_cols.append('Area (Max.)')

                out_rows: list[pd.Series] = []

                # Process groups
                for key, grp in df.groupby('_dup_key'):
                    if key is None or (isinstance(key, float) and pd.isna(key)):
                        # No key, keep originals
                        for _, r in grp.drop(columns=['_dup_key']).iterrows():
                            out_rows.append(r)
                        continue

                    if len(grp) == 1:
                        out_rows.append(grp.drop(columns=['_dup_key']).iloc[0])
                        continue

                    # RT window filtering (±configured minutes around mean)
                    if 'RT' in grp.columns:
                        rt_mean = pd.to_numeric(grp['RT'], errors='coerce').mean()
                        in_win = grp[pd.to_numeric(grp['RT'], errors='coerce').sub(rt_mean).abs() <= self.id_dedup_rt_window_minutes].copy()
                        out_win = grp.drop(in_win.index)
                    else:
                        in_win = grp.copy()
                        out_win = grp.iloc[0:0]

                    # Keep out-of-window rows unchanged
                    for _, r in out_win.drop(columns=['_dup_key']).iterrows():
                        out_rows.append(r)

                    if in_win.empty:
                        continue

                    # Choose representative by ion rank, MS2, Area
                    in_win['_ion_rank'] = in_win['Reference Ion'].apply(lambda x: _ion_rank(x, polarity)) if 'Reference Ion' in in_win.columns else 999
                    in_win['_ms2_rank'] = in_win['MS2'].apply(_ms2_rank) if 'MS2' in in_win.columns else 9
                    if 'Area (Max.)' in in_win.columns:
                        in_win['Area (Max.)'] = pd.to_numeric(in_win['Area (Max.)'], errors='coerce')
                    in_win = in_win.sort_values(['_ion_rank', '_ms2_rank', 'Area (Max.)'], ascending=[True, True, False])
                    rep = in_win.iloc[0].copy()

                    # Aggregate numeric columns
                    agg_row = rep.copy()
                    for col in numeric_cols:
                        if col in in_win.columns:
                            agg_row[col] = pd.to_numeric(in_win[col], errors='coerce').sum(min_count=1)
                    if 'RT' in in_win.columns:
                        agg_row['RT'] = pd.to_numeric(in_win['RT'], errors='coerce').mean()

                    # Drop helper cols
                    for c in ['_dup_key', '_ion_rank', '_ms2_rank']:
                        if c in agg_row.index:
                            agg_row = agg_row.drop(labels=c)

                    out_rows.append(agg_row)

                result = pd.DataFrame(out_rows)
                # Keep original column order when possible
                result = result[[c for c in df_in.columns if c in result.columns]]
                return result.reset_index(drop=True)

            def _merge_ids(one_df: Optional[pd.DataFrame], polarity: str) -> Optional[pd.DataFrame]:
                if one_df is None or one_df is pd.NA:
                    return None
                if not isinstance(one_df, pd.DataFrame) or one_df.empty:
                    return None

                dfp = one_df.copy()
                if 'Name' not in dfp.columns:
                    return None
                dfp['Name_Key'] = dfp['Name'].astype(str).str.strip().str.lower()

                # Merge IDs
                merged = dfp.merge(ids_df[['Name_Key'] + id_cols_present], on='Name_Key', how='left')

                # Compute feature columns present, accounting for alternate names
                features_present = []
                for col in feature_cols_expected:
                    if col in merged.columns:
                        features_present.append(col)
                    elif col in alt_feature_names and alt_feature_names[col] in merged.columns:
                        features_present.append(alt_feature_names[col])

                # Determine insertion index: after last feature column
                if features_present:
                    last_feature_idx = max(merged.columns.get_loc(c) for c in features_present)
                else:
                    # If no feature columns recognized, put IDs after 'Name' if present
                    last_feature_idx = merged.columns.get_loc('Name') if 'Name' in merged.columns else -1

                # Reorder columns: keep existing columns but move ID columns to right after features
                # Build final order
                cols = list(merged.columns)
                # Ensure ID columns are unique and present
                id_cols_in_merged = [c for c in id_cols_present if c in merged.columns]
                # Remove ID columns from current order
                cols_wo_ids = [c for c in cols if c not in id_cols_in_merged]
                # Insert ID columns after last_feature_idx position
                if last_feature_idx >= 0:
                    before = cols_wo_ids[:last_feature_idx + 1]
                    after = cols_wo_ids[last_feature_idx + 1:]
                    ordered = before + id_cols_in_merged + after
                else:
                    ordered = id_cols_in_merged + cols_wo_ids
                merged = merged[ordered]

                # Apply Formula+ID duplicate consolidation within ±2 min RT window (sum numeric)
                try:
                    merged = _dedup_by_formula_ids_sum(merged, polarity)
                except Exception as _dedup_err:
                    logger.warning(f"Formula+ID dedup (sum) failed for {polarity} data: {_dedup_err}")

                # Apply optional filtering
                if self.id_filter_mode != 'none':
                    selected = self.selected_id_columns or ['HMDB_ID', 'PubChem_CID', 'KEGG_ID', 'LipidMaps_ID', 'ChEBI_ID', 'CAS', 'SMILES', 'InChI', 'InChIKey']
                    selected = [c for c in selected if c in merged.columns]
                    if selected:
                        mask_any = False
                        for c in selected:
                            colmask = merged[c].astype(str).str.strip().ne('') & merged[c].notna()
                            mask_any = colmask if isinstance(mask_any, (pd.Series)) is False else (mask_any | colmask)
                        if isinstance(mask_any, pd.Series):
                            merged = merged[mask_any]

                if self.ms2_filter_mode != 'none' and 'MS2' in merged.columns:
                    if self.ms2_filter_mode == 'preferred_only':
                        merged = merged[merged['MS2'].astype(str).str.strip() == 'DDA for preferred ion']
                    elif self.ms2_filter_mode == 'preferred_or_other':
                        merged = merged[
                            merged['MS2'].astype(str).str.strip().isin(['DDA for preferred ion', 'DDA for other ion'])
                        ]

                merged = _enforce_endogenous_rules(merged)

                if self.require_endogenous_yes and 'Endogenous_Source' in merged.columns:
                    merged = merged[merged['Endogenous_Source'].astype(str).str.strip().str.lower() == 'yes']

                # Keep endogenous values explicit and non-empty in the saved outputs.
                for col in ['Endogenous', 'Endogenous_Source']:
                    if col in merged.columns:
                        merged[col] = merged[col].apply(_normalize_yes_no)
                    elif col == 'Endogenous' and 'HMDB_ID' in merged.columns:
                        merged[col] = 'No'

                merged = add_msi_confidence_columns(merged)
                merged = apply_confidence_filter(merged, self.confidence_filter_mode)
                return merged.reset_index(drop=True)

            results['pos_id_df'] = _merge_ids(self.pos_df, 'positive')
            results['neg_id_df'] = _merge_ids(self.neg_df, 'negative')

            return results
        except Exception as e:
            logger.error(f"Error assigning IDs to polarity frames: {e}")
            return {'pos_id_df': None, 'neg_id_df': None}

    def _get_cleaned_data_from_memory(self) -> Optional[pd.DataFrame]:
        """
        Get cleaned metabolite data DataFrame from memory if available
        
        Returns:
            Cleaned metabolite DataFrame or None if not available
        """
        # Check if cleaned data was passed during initialization
        logger.info(f"🔍 Checking for cleaned_metabolites_df attribute: {hasattr(self, 'cleaned_metabolites_df')}")
        if hasattr(self, 'cleaned_metabolites_df'):
            logger.info(f"🔍 cleaned_metabolites_df is not None: {self.cleaned_metabolites_df is not None}")
            if self.cleaned_metabolites_df is not None:
                logger.info(f"🔍 cleaned_metabolites_df type: {type(self.cleaned_metabolites_df)}")
                logger.info(f"🔍 cleaned_metabolites_df shape: {self.cleaned_metabolites_df.shape}")
                logger.info(f"🔍 cleaned_metabolites_df columns: {list(self.cleaned_metabolites_df.columns)}")
                # Normalize to two-column subset (Name, Formula) when possible
                df_mem = self.cleaned_metabolites_df.copy()
                if 'Name' in df_mem.columns and 'Formula' in df_mem.columns:
                    df_mem = df_mem[['Name', 'Formula']]
                elif 'Name' in df_mem.columns:
                    df_mem = df_mem[['Name']]
                logger.info(f"Found cleaned data in memory: {len(df_mem)} rows, columns: {list(df_mem.columns)}")
                return df_mem
        
        # Check if it's stored in a global variable or class variable
        # This would be set by the data cleaning module
        try:
            import __main__
            if hasattr(__main__, 'cleaned_metabolites_df'):
                df = getattr(__main__, 'cleaned_metabolites_df')
                if df is not None and not df.empty:
                    # Normalize to two-column subset if available
                    df2 = df.copy()
                    if 'Name' in df2.columns and 'Formula' in df2.columns:
                        df2 = df2[['Name', 'Formula']]
                    elif 'Name' in df2.columns:
                        df2 = df2[['Name']]
                    logger.info(f"Found cleaned data in global memory: {len(df2)} rows, columns: {list(df2.columns)}")
                    return df2
        except:
            pass
        
        logger.info("No cleaned metabolite data found in memory")
        return None

    def _merge_with_cleaned_data(self, annotated_df: pd.DataFrame, metabolite_ids_df: pd.DataFrame) -> pd.DataFrame:
        """
        Merge annotated DataFrame with cleaned metabolite data
        
        Args:
            annotated_df: DataFrame with ID annotations
            metabolite_ids_df: DataFrame with cleaned metabolite data from memory
            
        Returns:
            Enhanced merged DataFrame
        """
        try:
            logger.info(f"Using cleaned data from memory: {len(metabolite_ids_df)} rows, {len(metabolite_ids_df.columns)} columns")
            
            # Create Name_Key columns for matching
            annotated_df = annotated_df.copy()
            metabolite_ids_df = metabolite_ids_df.copy()
            
            annotated_df['Name_Key'] = annotated_df['Name'].astype(str).str.strip().str.lower()
            metabolite_ids_df['Name_Key'] = metabolite_ids_df['Name'].astype(str).str.strip().str.lower()
            
            # Find common keys
            annotated_keys = set(annotated_df['Name_Key'].dropna())
            metabolite_keys = set(metabolite_ids_df['Name_Key'].dropna())
            common_keys = annotated_keys & metabolite_keys
            
            logger.info(f"Found {len(common_keys)} common metabolites for merging")
            
            # Step 1: Populate empty IDs in annotated_df from metabolite_ids_df
            annotated_df = self._populate_empty_ids(annotated_df, metabolite_ids_df)
            
            # Step 2: Add new metabolites from metabolite_ids_df
            final_df = self._add_new_metabolites(annotated_df, metabolite_ids_df)
            
            # Step 3: Create Endogenous_Source column
            final_df = self._create_endogenous_source_column(final_df, metabolite_ids_df)
            
            # Step 4: Clean up unnecessary columns before saving
            columns_to_remove = ['Name_Key', 'Molecular_Weight', 'Endogenous', 'annotation_sources', '_from_metabolite_ids']
            # Also remove any columns ending with '_annotated' suffix
            annotated_columns = [col for col in final_df.columns if col.endswith('_annotated')]
            columns_to_remove.extend(annotated_columns)
            
            columns_to_drop = [col for col in columns_to_remove if col in final_df.columns]
            
            if columns_to_drop:
                logger.info(f"Removing columns: {columns_to_drop}")
                final_df = final_df.drop(columns=columns_to_drop)
            
            logger.info(f"Merge completed: {len(final_df)} rows (added {len(final_df) - len(annotated_df)} new metabolites)")
            
            return final_df
            
        except Exception as e:
            logger.error(f"Error in merge process: {e}")
            raise

    def _populate_empty_ids(self, annotated_df: pd.DataFrame, metabolite_ids_df: pd.DataFrame) -> pd.DataFrame:
        """Populate empty IDs in annotated_df from metabolite_ids_df where Name_Key matches"""
        logger.info("Populating empty IDs from cleaned data...")
        
        # ID columns to populate
        id_columns = ['HMDB_ID', 'KEGG_ID', 'PubChem_CID', 'InChIKey', 'CAS', 'SMILES', 'InChI', 'IUPAC_Name']
        
        # Create a mapping from Name_Key to metabolite_ids_df data
        metabolite_mapping = metabolite_ids_df.set_index('Name_Key').to_dict('index')
        
        population_stats = {}
        
        for id_col in id_columns:
            if id_col in annotated_df.columns and id_col in metabolite_ids_df.columns:
                original_count = (annotated_df[id_col].notna() & (annotated_df[id_col].astype(str).str.strip() != '')).sum()
                
                # Populate empty values
                for idx, row in annotated_df.iterrows():
                    name_key = row['Name_Key']
                    current_value = row[id_col]
                    
                    # If current value is empty and we have data in metabolite_ids_df
                    if (pd.isna(current_value) or str(current_value).strip() == '') and name_key in metabolite_mapping:
                        new_value = metabolite_mapping[name_key].get(id_col, '')
                        if pd.notna(new_value) and str(new_value).strip() != '':
                            annotated_df.loc[idx, id_col] = new_value
                
                final_count = (annotated_df[id_col].notna() & (annotated_df[id_col].astype(str).str.strip() != '')).sum()
                populated = final_count - original_count
                
                if populated > 0:
                    population_stats[id_col] = populated
                    logger.info(f"  {id_col}: {original_count} -> {final_count} (+{populated} populated)")
        
        if population_stats:
            logger.info(f"ID Population Results: {population_stats}")
        else:
            logger.info("No IDs were populated from cleaned data")
        
        return annotated_df

    def _add_new_metabolites(self, annotated_df: pd.DataFrame, metabolite_ids_df: pd.DataFrame) -> pd.DataFrame:
        """Add new metabolites from metabolite_ids_df that are not in annotated_df"""
        logger.info("Adding new metabolites from cleaned data...")
        
        # Add tracking column to original annotated data
        annotated_df = annotated_df.copy()
        annotated_df['_from_metabolite_ids'] = False
        
        # Find metabolites in metabolite_ids_df that are not in annotated_df
        annotated_keys = set(annotated_df['Name_Key'].dropna())
        metabolite_keys = set(metabolite_ids_df['Name_Key'].dropna())
        
        new_keys = metabolite_keys - annotated_keys
        logger.info(f"New metabolites to add: {len(new_keys)}")
        
        if len(new_keys) == 0:
            return annotated_df
        
        # Get new metabolites
        new_metabolites = metabolite_ids_df[metabolite_ids_df['Name_Key'].isin(new_keys)].copy()
        
        # Prepare new metabolites to match annotated_df structure
        new_rows = []
        for _, row in new_metabolites.iterrows():
            new_row = {}
            
            # Map columns from metabolite_ids_df to annotated_df structure
            column_mapping = {
                'Name': 'Name',
                'Formula': 'Formula', 
                'IUPAC_Name': 'IUPAC_Name',
                'KEGG_ID': 'KEGG_ID',
                'CAS': 'CAS',
                'SMILES': 'SMILES',
                'InChI': 'InChI',
                'InChIKey': 'InChIKey',
                'HMDB_ID': 'HMDB_ID',
                'PubChem_CID': 'PubChem_CID'
                # Note: Skipping 'Compound Class' to avoid conflicts
            }
            
            for source_col, target_col in column_mapping.items():
                if source_col in row.index:
                    new_row[target_col] = row[source_col]
            
            # Fill other columns with empty values
            for col in annotated_df.columns:
                if col not in new_row:
                    new_row[col] = ''
            
            # Add source tracking for endogenous classification
            new_row['_from_metabolite_ids'] = True
            
            new_rows.append(new_row)
        
        # Convert to DataFrame and append
        if new_rows:
            new_metabolites_df = pd.DataFrame(new_rows)
            new_metabolites_df = new_metabolites_df.reindex(columns=annotated_df.columns, fill_value='')
            final_df = pd.concat([annotated_df, new_metabolites_df], ignore_index=True)
            logger.info(f"Added {len(new_metabolites_df)} new metabolites")
            return final_df
        
        return annotated_df

    def _create_endogenous_source_column(self, df: pd.DataFrame, metabolite_ids_df: pd.DataFrame) -> pd.DataFrame:
        """Create Endogenous_Source column based on metabolite presence and Endogenous column"""
        logger.info("Creating Endogenous_Source column...")
        
        df = df.copy()
        df['Endogenous_Source'] = 'No'
        
        # Step 1: Set "Yes" for ANY metabolite that came from metabolite_ids_df (using tracking column)
        if '_from_metabolite_ids' in df.columns:
            metabolite_ids_mask = df['_from_metabolite_ids'] == True
            metabolite_ids_count = metabolite_ids_mask.sum()
            
            if metabolite_ids_count > 0:
                df.loc[metabolite_ids_mask, 'Endogenous_Source'] = 'Yes'
                logger.info(f"Set Endogenous_Source=Yes for {metabolite_ids_count} rows from metabolite IDs data")
        else:
            logger.warning("Tracking column '_from_metabolite_ids' not found - using fallback name matching")
            # Fallback to old method if tracking column is missing
            if 'Name_Key' not in df.columns:
                df['Name_Key'] = df['Name'].astype(str).str.strip().str.lower()
            
            metabolite_name_keys = set(metabolite_ids_df['Name_Key'].dropna()) if 'Name_Key' in metabolite_ids_df.columns else set()
            
            if len(metabolite_name_keys) > 0:
                metabolite_ids_mask = df['Name_Key'].isin(metabolite_name_keys)
                metabolite_ids_count = metabolite_ids_mask.sum()
                
                if metabolite_ids_count > 0:
                    df.loc[metabolite_ids_mask, 'Endogenous_Source'] = 'Yes'
                    logger.info(f"Set Endogenous_Source=Yes for {metabolite_ids_count} rows present in cleaned data")
        
        # Step 2: Set "Yes" for rows with "Yes" in Endogenous column (don't overwrite existing "Yes")
        if 'Endogenous' in df.columns:
            endogenous_mask = (df['Endogenous'].astype(str).str.strip().str.lower() == 'yes') & (df['Endogenous_Source'] != 'Yes')
            endogenous_count = endogenous_mask.sum()
            if endogenous_count > 0:
                df.loc[endogenous_mask, 'Endogenous_Source'] = 'Yes'
                logger.info(f"Set Endogenous_Source=Yes for {endogenous_count} additional rows with Endogenous=Yes")
        
        # Report final statistics
        yes_count = (df['Endogenous_Source'] == 'Yes').sum()
        no_count = (df['Endogenous_Source'] == 'No').sum()
        
        logger.info(f"Final Endogenous_Source distribution: Yes={yes_count} ({yes_count/len(df)*100:.1f}%), No={no_count} ({no_count/len(df)*100:.1f}%)")
        
        return df

def main():
    """Command line interface for metabolite ID annotation"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Metabolite/Lipid ID Annotator')
    parser.add_argument('--input', '-i', required=True, help='Input Excel file')
    parser.add_argument('--output', '-o', default='metabolite_ids_annotated.xlsx', help='Output Excel file')
    parser.add_argument('--workers', '-w', type=int, default=6, help='Number of parallel workers')
    # Modes
    parser.add_argument('--lipid-mode', action='store_true', help='Run Lipid ID mode (offline lipid feather + PubChem first, then KEGG/HMDB)')
    # Database skip options
    parser.add_argument('--skip-pubchem', action='store_true', help='Skip PubChem API searches (use only LipidMaps data)')
    parser.add_argument('--skip-kegg', action='store_true', help='Skip KEGG API searches')
    parser.add_argument('--skip-hmdb', action='store_true', help='Skip HMDB offline searches')
    
    args = parser.parse_args()
    
    annotator = MetaboliteIDAnnotator(
        args.input,
        args.output,
        lipid_mode=args.lipid_mode,
        skip_pubchem=args.skip_pubchem,
        skip_kegg=args.skip_kegg,
        skip_hmdb=args.skip_hmdb,
    )
    
    if args.lipid_mode:
        annotator.run_lipid_id_annotation(args.workers)
    else:
        annotator.run_id_annotation(args.workers)

def refilter_annotated_excel(
    input_file: str,
    output_file: str,
    selected_id_columns: Optional[List[str]] = None,
    skip_id_filtering: bool = False,
    ms2_filter_mode: str = 'none',
    confidence_filter_mode: str = 'exclude_low',
    require_endogenous_yes: bool = False,
    progress_callback: Optional[callable] = None
) -> Dict[str, Any]:
    """
    Re-filter an existing annotated Excel file without re-running annotation.
    
    Args:
        input_file: Path to the annotated Excel file (with Merged_IDs, Pos_id, Neg_id sheets)
        output_file: Path to save the re-filtered results
        selected_id_columns: List of ID column names to use for filtering (e.g., ['HMDB_ID', 'KEGG_ID'])
                            If empty or None, no ID-based filtering is applied
        skip_id_filtering: If True, skip all ID-based filtering (save all metabolites)
        ms2_filter_mode: MS2 filter mode ('none', 'preferred_only', 'preferred_or_other')
        confidence_filter_mode: Confidence filter mode ('exclude_low' or 'none')
        require_endogenous_yes: If True, require Endogenous_Source == 'Yes'
        progress_callback: Optional callback function to report progress
    
    Returns:
        Dictionary with keys: 'success', 'error', 'summary'
    """
    def log(msg):
        logger.info(msg)
        if progress_callback:
            progress_callback(msg)
    
    try:
        log(f"Loading annotated Excel file: {input_file}")
        
        # Read all sheets
        xl_file = pd.ExcelFile(input_file)
        sheet_names = xl_file.sheet_names
        
        log(f"Found sheets: {', '.join(sheet_names)}")
        
    # Check for required sheets (at least one of Pos_id or Neg_id must exist)
        has_pos = 'Pos_id' in sheet_names
        has_neg = 'Neg_id' in sheet_names
        has_merged = 'Merged_IDs' in sheet_names
        
        if not has_pos and not has_neg:
            return {
                'success': False,
                'error': 'Excel file must contain at least one of: Pos_id or Neg_id sheets'
            }
        
        log(f"Valid sheets found - Pos_id: {has_pos}, Neg_id: {has_neg}, Merged_IDs: {has_merged}")
        
        # Determine if this annotated file is from lipid mode
        lipid_mode_detected = False
        try:
            probe_df = None
            if has_pos:
                probe_df = pd.read_excel(input_file, sheet_name='Pos_id', nrows=1)
            elif has_neg:
                probe_df = pd.read_excel(input_file, sheet_name='Neg_id', nrows=1)
            elif has_merged:
                probe_df = pd.read_excel(input_file, sheet_name='Merged_IDs', nrows=1)
            if probe_df is not None and not probe_df.empty:
                cols_norm = {str(c).strip() for c in probe_df.columns}
                # Heuristics: Lipid outputs contain a LipidID column and often ABBREVIATION/MAIN_CLASS
                if ('LipidID' in cols_norm) or ({'ABBREVIATION', 'MAIN_CLASS'} & cols_norm):
                    lipid_mode_detected = True
        except Exception:
            # Non-fatal; default to metabolite mode behavior
            pass

        if lipid_mode_detected:
            log("Lipid mode detected from annotated Excel schema → bypassing ID/MS2/Endogenous filters.")
            # Force-disable all filters for lipid-mode re-filtering
            selected_id_columns = []
            ms2_filter_mode = 'none'
            confidence_filter_mode = 'none'
            require_endogenous_yes = False

        # Helper function to apply filters to a dataframe
        def apply_filters(df: pd.DataFrame, sheet_name: str) -> pd.DataFrame:
            if df is None or df.empty:
                return df
            
            original_count = len(df)
            log(f"  {sheet_name}: Starting with {original_count} rows")
            
            # Apply ID column filter (skipped if skip_id_filtering is True, or for lipid mode)
            if (not lipid_mode_detected) and (not skip_id_filtering) and selected_id_columns and len(selected_id_columns) > 0:
                selected = [c for c in selected_id_columns if c in df.columns]
                if selected:
                    mask_any = pd.Series([False] * len(df), index=df.index)
                    for c in selected:
                        colmask = df[c].astype(str).str.strip().ne('') & df[c].notna()
                        mask_any = mask_any | colmask
                    df = df[mask_any].copy()
                    log(f"  {sheet_name}: After ID filter ({', '.join(selected)}): {len(df)} rows")
            elif skip_id_filtering and not lipid_mode_detected:
                log(f"  {sheet_name}: ID filtering skipped (all {original_count} rows retained)")
            
            # Apply MS2 filter (for Pos_id and Neg_id only) — skipped for lipid mode
            if (not lipid_mode_detected) and sheet_name in ['Pos_id', 'Neg_id'] and ms2_filter_mode != 'none':
                if 'MS2' in df.columns:
                    before_ms2 = len(df)
                    if ms2_filter_mode == 'preferred_only':
                        df = df[df['MS2'].astype(str).str.strip() == 'DDA for preferred ion'].copy()
                    elif ms2_filter_mode == 'preferred_or_other':
                        df = df[df['MS2'].astype(str).str.strip().isin(['DDA for preferred ion', 'DDA for other ion'])].copy()
                    log(f"  {sheet_name}: After MS2 filter ({ms2_filter_mode}): {len(df)} rows (removed {before_ms2 - len(df)})")
            
            # Apply Endogenous_Source filter — skipped for lipid mode
            if (not lipid_mode_detected) and require_endogenous_yes and 'Endogenous_Source' in df.columns:
                before_endo = len(df)
                df = df[df['Endogenous_Source'].astype(str).str.strip().str.lower() == 'yes'].copy()
                log(f"  {sheet_name}: After Endogenous filter: {len(df)} rows (removed {before_endo - len(df)})")

            # Add MSI/confidence annotation columns on filtered results.
            df = add_msi_confidence_columns(df)
            before_conf = len(df)
            df = apply_confidence_filter(df, confidence_filter_mode)
            if confidence_filter_mode != 'none':
                log(f"  {sheet_name}: After confidence filter ({confidence_filter_mode}): {len(df)} rows (removed {before_conf - len(df)})")
            
            removed = original_count - len(df)
            log(f"  {sheet_name}: Final count: {len(df)} rows (removed {removed} total)")
            
            return df.reset_index(drop=True)
        
        # Filter Pos and Neg sheets
        pos_filtered = None
        neg_filtered = None
        
        if has_pos:
            log("Filtering Pos_id sheet...")
            pos_df = pd.read_excel(input_file, sheet_name='Pos_id')
            pos_filtered = apply_filters(pos_df, 'Pos_id')
        
        if has_neg:
            log("Filtering Neg_id sheet...")
            neg_df = pd.read_excel(input_file, sheet_name='Neg_id')
            neg_filtered = apply_filters(neg_df, 'Neg_id')
        
        # Re-create Merged_IDs from filtered Pos and Neg
        log("Re-creating Merged_IDs sheet...")
        
        merged_dfs = []
        if pos_filtered is not None and not pos_filtered.empty:
            merged_dfs.append(pos_filtered)
        if neg_filtered is not None and not neg_filtered.empty:
            merged_dfs.append(neg_filtered)
        
        if not merged_dfs:
            return {
                'success': False,
                'error': 'No data remaining after applying filters'
            }
        
        # Concatenate and remove duplicates
        merged_filtered = pd.concat(merged_dfs, ignore_index=True)
        
        # Deduplicate appropriately: prefer 'LipidID' when present; otherwise fall back to 'Name'
        if 'LipidID' in merged_filtered.columns:
            before_dedup = len(merged_filtered)
            merged_filtered = merged_filtered.drop_duplicates(subset=['LipidID'], keep='first')
            log(f"Merged_IDs: Removed {before_dedup - len(merged_filtered)} duplicates based on 'LipidID'")
        elif 'Name' in merged_filtered.columns:
            before_dedup = len(merged_filtered)
            merged_filtered = merged_filtered.drop_duplicates(subset=['Name'], keep='first')
            log(f"Merged_IDs: Removed {before_dedup - len(merged_filtered)} duplicates based on 'Name'")

        merged_filtered = add_msi_confidence_columns(merged_filtered)
        merged_filtered = apply_confidence_filter(merged_filtered, confidence_filter_mode)
        
        log(f"Merged_IDs: Final count: {len(merged_filtered)} rows")
        
        # Save to Excel with all three sheets
        log(f"Saving re-filtered results to: {output_file}")
        
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            merged_filtered.to_excel(writer, index=False, sheet_name='Merged_IDs')
            if pos_filtered is not None and not pos_filtered.empty:
                pos_filtered.to_excel(writer, index=False, sheet_name='Pos_id')
            if neg_filtered is not None and not neg_filtered.empty:
                neg_filtered.to_excel(writer, index=False, sheet_name='Neg_id')
        
        # Create summary
        summary_lines = []
        summary_lines.append(f"Merged_IDs: {len(merged_filtered)} rows")
        if pos_filtered is not None:
            summary_lines.append(f"Pos_id: {len(pos_filtered)} rows")
        if neg_filtered is not None:
            summary_lines.append(f"Neg_id: {len(neg_filtered)} rows")
        
        summary = "\n".join(summary_lines)
        log("Re-filtering completed successfully!")
        
        return {
            'success': True,
            'summary': summary,
            'output_file': output_file
        }
        
    except Exception as e:
        error_msg = f"Re-filtering failed: {str(e)}"
        logger.error(error_msg)
        return {
            'success': False,
            'error': error_msg
        }


if __name__ == "__main__":
    main()
