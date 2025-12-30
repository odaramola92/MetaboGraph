#!/usr/bin/env python3
"""
Metabolite Pathway Mapper - Phase 2
Maps metabolites with pre-annotated IDs to pathway databases
"""

import pandas as pd
import numpy as np
import logging
import os
import pickle
from typing import Dict, List, Optional, Any, Tuple
import re
import time
import socket
import urllib.request
import urllib.error
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# Import centralized pathway filtering function
from main_script.metabolite_pathway_network import (
    should_filter_pathway, 
    is_disease_pathway_global,
    normalize_disease_pathway_name as normalize_disease_pathway_name_global
)
from threading import Lock

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
    "transient": "Transient Tyrosinemia of the Newborn",
    "of the newborn": "Transient Tyrosinemia of the Newborn",
    "tyrosinemia,transient;of the newborn": "Transient Tyrosinemia of the Newborn",
}

def normalize_disease_pathway_name(pathway_name: str) -> str:
    """Normalize disease pathway names using hardcoded corrections"""
    # Use the global function from metabolite_pathway_network
    return normalize_disease_pathway_name_global(pathway_name)

def is_disease_pathway(pathway_name: str) -> bool:
    """
    Check if pathway should be treated as a disease.
    Uses centralized logic from metabolite_pathway_network.
    """
    return is_disease_pathway_global(pathway_name)

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

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants for Reactome and KEGG APIs
REACTOME_BASE = "https://reactome.org/ContentService/data"
KEGG_BASE = "https://rest.kegg.jp"
SPECIES_TAXID = {
    "Homo sapiens": "9606",
    "Mus musculus": "10090",
    "Rattus norvegicus": "10116",
}
KEGG_ORGANISMS = {
    "Homo sapiens": "hsa",
    "Mus musculus": "mmu",
    "Rattus norvegicus": "rno",
}

# --- Centralized pathway filtering constants ---
# Default exclusion keywords for plant/microbial/non-mammalian and non-specific pathways
DEFAULT_EXCLUDE_PATHWAY_KEYWORDS: List[str] = [
    # Plant-related
    "plant", "chlorophyll", "photosynthesis", "cellulose", "lignin",
    # Microbial and pathogens
    "bacterial", "bacterium", "bacteria", "microbial", "microorganism", "prokaryotic", "prokaryote", "archaeal", "fungal", "yeast",
    "escherichia", "salmonella", "streptococcus", "staphylococcus",
    # Marine/aquatic
    "algae", "algal", "cyanobacteria", "marine", "seaweed",
    # Parasites/viruses
    "parasitic", "parasite", "protozoa", "plasmodium", "viral", "virus",
    # Invertebrates
    "insect", "drosophila", "caenorhabditis", "arthropod", "nematode",
    # Drug/xenobiotic/biodegradation
    "antibiotic", "xenobiotic", "biodegradation", "methanogenesis",
    # Drug action pathways (exclude from metabolic set)
    "action pathway", "drug action", "pharmacodynamics", "pharmacokinetics",
    # Superpathways and secondary metabolism (often non-specific or non-mammalian)
    "secondary metabolite", "polyketide", "nonribosomal", "terpenoid",
    # Generic catch-alls
    "other",
]

# Non-specific labels to drop (prefix/contains checks, lower-cased)
NON_SPECIFIC_PATHWAY_PATTERNS: List[str] = [
    "biochemical pathways: part",  # e.g., "Biochemical Pathways: Part I/II"
    "hereditary",                  # explicitly exclude generic 'Hereditary' labels
]

def apply_annotation_filters(pathways: List[str], exclude_keywords: List[str] = None,
                             non_specific_patterns: List[str] = None) -> Tuple[List[str], List[str]]:
    """Pure helper to apply pathway filtering and disease relocation.

    Returns (retained_paths, moved_diseases).

    Filtering logic:
      1. Apply centralized pathway filtering (diseases + overly general pathways)
      2. Exclude plant/microbial/drug/superpathway/non-specific labels.
      3. Exclude very short (<=2 chars) or numeric-only tokens.
      4. Move disease-like pathways to disease list.
    """
    if exclude_keywords is None:
        exclude_keywords = [kw.lower() for kw in DEFAULT_EXCLUDE_PATHWAY_KEYWORDS]
    else:
        exclude_keywords = [kw.lower() for kw in exclude_keywords]
    if non_specific_patterns is None:
        non_specific_patterns = [p.lower() for p in NON_SPECIFIC_PATHWAY_PATTERNS]
    else:
        non_specific_patterns = [p.lower() for p in non_specific_patterns]

    def _is_excluded(name_lower: str) -> bool:
        for kw in exclude_keywords:
            if kw in name_lower:
                return True
        for pat in non_specific_patterns:
            if pat in name_lower:
                return True
        if len(name_lower.strip()) <= 2:
            return True
        if name_lower.replace('.', '').replace('-', '').isdigit():
            return True
        return False

    # First pass: apply centralized filtering to catch diseases and overly general pathways
    pre_filtered: List[str] = []
    for p in pathways:
        if not should_filter_pathway(p):
            pre_filtered.append(p)
    
    # Second pass: apply annotation-specific exclusions
    filtered: List[str] = []
    for p in pre_filtered:
        pl = p.lower()
        if not _is_excluded(pl):
            filtered.append(p)

    # Separate diseases from regular pathways
    moved_diseases: List[str] = []
    retained: List[str] = []
    for p in filtered:
        if is_disease_pathway(p):
            moved_diseases.append(normalize_disease_pathway_name(p))
        else:
            retained.append(p)

    # Deduplicate disease list preserving order
    seen = set()
    dedup_diseases: List[str] = []
    for d in moved_diseases:
        dl = d.lower()
        if dl not in seen:
            seen.add(dl)
            dedup_diseases.append(d)

    return retained, dedup_diseases

def fetch_json(url: str, retries: int = 1, timeout: int = 15) -> Dict:
    """Fetch JSON from URL with retries (default single attempt to avoid long hangs)."""
    backoff = 1.0
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "metabolite-mapper/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                if resp.status != 200:
                    raise urllib.error.HTTPError(url, resp.status, resp.reason, resp.headers, None)
                raw = resp.read()
                return json.loads(raw.decode("utf-8"))
        except Exception as e:
            if attempt == retries:
                raise
            time.sleep(backoff * (2 ** (attempt - 1)))

class ReactomeKEGGCache:
    """Simple file-based cache for Reactome and KEGG API responses."""
    
    def __init__(self, cache_dir: str = ".reactome_kegg_cache"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
    
    def _make_key(self, identifier: str, source: str, organism: str) -> str:
        return f"{source}_{identifier}_{organism}.json"
    
    def get(self, identifier: str, source: str, organism: str) -> Optional[List]:
        """Retrieve cached pathways."""
        key_file = self.cache_dir / self._make_key(identifier, source, organism)
        if key_file.exists():
            try:
                with open(key_file) as f:
                    data = json.load(f)
                return data.get("pathways", [])
            except Exception:
                pass
        return None
    
    def set(self, identifier: str, source: str, organism: str, pathways: List) -> None:
        """Cache pathways."""
        key_file = self.cache_dir / self._make_key(identifier, source, organism)
        try:
            with open(key_file, 'w') as f:
                json.dump({"pathways": pathways}, f)
        except Exception:
            pass


def _safe_get_scalar(row, key, default=''):
    """
    Safely extract a scalar value from a pandas Series row.
    
    Handles the case where duplicate column names cause row[key] to return 
    a Series instead of a scalar value.
    
    Parameters:
    -----------
    row : pd.Series
        A row from a DataFrame
    key : str
        Column name to extract
    default : any
        Default value if key is missing or extraction fails
        
    Returns:
    --------
    Scalar value (string, number, etc.)
    """
    try:
        value = row.get(key, default)
        
        # If we got a Series (duplicate columns), take the first value
        if isinstance(value, pd.Series):
            if len(value) > 0:
                value = value.iloc[0]
            else:
                value = default
        
        # Handle NaN
        if pd.isna(value):
            return default
            
        return value
    except Exception:
        return default


class MetabolitePathwayMapper:
    """
    Phase 2: Pathway Mapping using Pre-annotated Metabolite IDs
    Input: Excel file with metabolite names and IDs from Phase 1
    Output: Metabolites with pathway information ready for network analysis
    """

    def __init__(self, input_file: str = "metabolite_ids_annotated.xlsx", 
                 output_file: str = "metabolite_pathways_mapped.xlsx",
                 sheet_name: Optional[str] = None,
                 progress_callback=None,
                 exclude_pathway_keywords: Optional[List[str]] = None,
                 organism: str = 'Homo sapiens',
                 max_workers: int = 6,
                 max_metabolite_workers: int = 6,
                 cancel_flag=None):
        """
        Initialize the pathway mapper with organism-specific database files
        
        Args:
            input_file: Path to input Excel file with metabolites
            output_file: Path to output Excel file for results
            sheet_name: Sheet name in Excel file
            progress_callback: Optional callback function(current, total, message)
            exclude_pathway_keywords: Keywords to exclude pathways (e.g., plant/bacteria pathways)
            organism: Organism for pathway filtering ('Homo sapiens', 'Rattus norvegicus', 'Mus musculus')
            max_workers: Number of parallel workers for database queries per metabolite (default: 6)
            max_metabolite_workers: Number of metabolites to process in parallel (default: 4)
            cancel_flag: Optional callable that returns True if processing should be cancelled
        """
        self.input_file = input_file
        self.output_file = output_file
        self.sheet_name = sheet_name
        self.progress_callback = progress_callback
        self.organism = organism
        self.max_workers = max_workers
        self.max_metabolite_workers = max_metabolite_workers
        self.progress_lock = Lock()  # Thread-safe progress tracking
        self.cancel_flag = cancel_flag  # Callable to check for cancellation
        
        # Set default exclusion keywords using centralized constants
        if exclude_pathway_keywords is None:
            exclude_pathway_keywords = DEFAULT_EXCLUDE_PATHWAY_KEYWORDS
        self.exclude_pathway_keywords = [kw.lower() for kw in exclude_pathway_keywords]

        # Load offline pathway databases
        logger.info(f"Loading offline pathway databases for {organism}...")
        
        # Helper function to find database file - import from shared utils
        from gui.shared.utils import find_database
        
        # Load HMDB database
        hmdb_path = find_database("hmdb_database.feather")
        if hmdb_path is None:
            error_msg = f"❌ Required database file 'hmdb_database.feather' not found! Please ensure it exists in a 'Databases' folder in the working directory."
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        try:
            self.hmdb_df = pd.read_feather(hmdb_path)
        except Exception as e:
            logger.error(f"Error loading HMDB database from {hmdb_path}: {e}")
            raise
        
        # Load PathBank database
        pathbank_path = find_database("pathbank_selected.feather")
        if pathbank_path is None:
            error_msg = f"❌ Required database file 'pathbank_selected.feather' not found! Please ensure it exists in a 'Databases' folder in the working directory."
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        try:
            self.pathbank_df = pd.read_feather(pathbank_path)
        except Exception as e:
            logger.error(f"Error loading PathBank database from {pathbank_path}: {e}")
            raise
        
        # Load SMPDB database
        smpdb_path = find_database("merged_SMP_metabolites.feather")
        if smpdb_path is None:
            error_msg = f"❌ Required database file 'merged_SMP_metabolites.feather' not found! Please ensure it exists in a 'Databases' folder in the working directory."
            logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        try:
            self.SMPDB_df = pd.read_feather(smpdb_path)
        except Exception as e:
            logger.error(f"Error loading SMPDB database from {smpdb_path}: {e}")
            raise
        
        # Load organism-specific WikiPathways database
        organism_map = {
            'Homo sapiens': 'homo_sapiens',
            'Rattus norvegicus': 'rattus_norvegicus',
            'Mus musculus': 'mus_musculus'
        }
        
        organism_file = organism_map.get(organism, 'homo_sapiens')
        wikipathways_filename = f"wikipathways_{organism_file}.feather"
        wikipathways_path = find_database(wikipathways_filename)
        
        if wikipathways_path:
            try:
                self.wikipathways_df = pd.read_feather(wikipathways_path)
                logger.info(f"Loaded WikiPathways ({organism}): {len(self.wikipathways_df)} entries from {wikipathways_path}")
            except Exception as e:
                logger.warning(f"Error loading WikiPathways from {wikipathways_path}: {e}, trying fallback")
                wikipathways_path = None
        
        if not wikipathways_path:
            # Try fallback combined database
            fallback_path = find_database("wikipathways_df.feather")
            if fallback_path:
                logger.warning(f"WikiPathways database not found: {wikipathways_filename}, using combined database at {fallback_path}")
                self.wikipathways_df = pd.read_feather(fallback_path)
            else:
                error_msg = f"❌ WikiPathways database files not found! Tried {wikipathways_filename} and wikipathways_df.feather"
                logger.error(error_msg)
                raise FileNotFoundError(error_msg)

        logger.info(f"Loaded HMDB: {len(self.hmdb_df)} entries")
        
        # Filter PathBank by species
        if 'Species' in self.pathbank_df.columns:
            original_count = len(self.pathbank_df)
            self.pathbank_df = self.pathbank_df[self.pathbank_df['Species'] == organism].copy()
            logger.info(f"Loaded PathBank ({organism}): {len(self.pathbank_df)} entries (filtered from {original_count})")
        else:
            logger.warning("PathBank database does not have Species column, using all species")
            logger.info(f"Loaded PathBank: {len(self.pathbank_df)} entries") 
            
        logger.info(f"Loaded SMPDB: {len(self.SMPDB_df)} entries")
        
        # Load Metabolika pathways cache from data cleaner
        self.metabolika_cache = self._load_metabolika_cache()

        # Standardize column names across databases
        self._standardize_column_names()

    def _standardize_column_names(self):
        """Standardize column names across different databases"""
        # Mapping of possible column name variations to standard names
        column_mappings = {
            'hmdb_df': {
                'HMDB': 'HMDB_ID', 'HMDB ID': 'HMDB_ID',
                'Pathways': 'HMDB_Pathways', 'Pathway': 'HMDB_Pathways',
                'Enzymes': 'Enzymes_Accession', 'Enzyme': 'Enzymes_Accession'
            },
            'pathbank_df': {
                'PathBank_ID': 'PathBank_ID', 'PathBank': 'PathBank_ID',
                'Pathway Name': 'PathBank_Pathways', 'Pathways': 'PathBank_Pathways'
            },
            'SMPDB_df': {
                'SMPDB_ID': 'SMPDB_ID', 'SMPDB': 'SMPDB_ID', 
                'Pathway': 'SMPDB_Pathways', 'Pathways': 'SMPDB_Pathways'
            },
            'wikipathways_df': {
                'WikiPathways_ID': 'WikiPathways_ID',
                'Pathway Name': 'WikiPathways', 'Pathways': 'WikiPathways'
            }
        }

        for df_name, mappings in column_mappings.items():
            df = getattr(self, df_name)
            for old_name, new_name in mappings.items():
                if old_name in df.columns and new_name not in df.columns:
                    df.rename(columns={old_name: new_name}, inplace=True)

    def _load_metabolika_cache(self) -> pd.DataFrame:
        """Load Metabolika pathways cache from data cleaner"""
        cache_file = "metabolika_cache.pkl"
        try:
            if os.path.exists(cache_file):
                with open(cache_file, 'rb') as f:
                    cache_data = pickle.load(f)
                    if isinstance(cache_data, pd.DataFrame) and not cache_data.empty:
                        logger.info(f"Loaded {len(cache_data)} entries from Metabolika cache")
                        return cache_data
                    else:
                        logger.info(f"Metabolika cache file exists but is empty or invalid")
                        return pd.DataFrame()
            else:
                logger.info(f"Metabolika cache file not found")
                return pd.DataFrame()
        except Exception as e:
            logger.warning(f"Error loading Metabolika cache: {e}")
            return pd.DataFrame()

    def _save_metabolika_cache(self) -> None:
        """Persist the in-memory Metabolika cache to disk"""
        cache_file = "metabolika_cache.pkl"
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(self.metabolika_cache, f)
            logger.info(f"Saved Metabolika cache with {len(self.metabolika_cache)} entries")
        except Exception as e:
            logger.warning(f"Error saving Metabolika cache: {e}")

    def _clean_pathway_string(self, pathway_str: str, preserve_separator: Optional[str] = None) -> str:
        """Clean and normalize pathway strings"""
        if pd.isna(pathway_str) or pathway_str == "":
            return ""
        
        # Remove database annotations in brackets
        cleaned = re.sub(r'\[.*?\]', '', str(pathway_str))
        cleaned = re.sub(r'\(.*?\)', '', cleaned)
        
        # Clean up malformed parentheses patterns like "PC)" or "PE)"
        cleaned = re.sub(r'\bPC\)', 'PC', cleaned)
        cleaned = re.sub(r'\bPE\)', 'PE', cleaned)
        
        # If we need to preserve a specific separator, don't normalize it
        if preserve_separator:
            # Only clean whitespace around the preserved separator, but don't change the separator itself
            # Clean up extra whitespace
            cleaned = re.sub(r'\s+', ' ', cleaned)
            # Clean whitespace around the preserved separator
            pattern = r'\s*' + re.escape(preserve_separator) + r'\s*'
            cleaned = re.sub(pattern, preserve_separator, cleaned)
            # Remove leading/trailing separators and whitespace
            cleaned = cleaned.strip(f' {preserve_separator}')
        else:
            # Normalize separators and whitespace (original behavior)
            cleaned = re.sub(r'[|,;]+', ';', cleaned)
            cleaned = re.sub(r'\s*;\s*', '; ', cleaned)
            cleaned = re.sub(r'\s+', ' ', cleaned)
            cleaned = cleaned.strip(' ;,|')
        
        return cleaned

    @staticmethod
    def group_semicolon_fragments(pathway_str: str) -> str:
        """Join semicolon-fragmented pathway names when they appear as sequences of short tokens

        Heuristic: If multiple short tokens (likely metabolite names) are followed by a token
        containing a pathway keyword (degradation, metabolism, synthesis, pathway, cycle, transport, etc.),
        join the run of tokens up through that keyword token into a single phrase by replacing
        the internal semicolons with spaces. Returns a new pathway string with grouped fragments
        separated by semicolons.
        """
        try:
            if pd.isna(pathway_str) or not pathway_str:
                return pathway_str

            # Split on semicolons but keep empty tokens out
            parts = [p.strip() for p in re.split(r';', str(pathway_str)) if p and p.strip()]
            if len(parts) <= 1:
                return pathway_str

            # Keywords that indicate a pathway phrase
            pathway_keywords = {
                'degradation', 'metabolism', 'biosynthesis', 'synthesis', 'pathway',
                'transport', 'cycle', 'signaling', 'catabolism', 'oxidation', 'metabolic'
            }

            grouped = []
            i = 0
            n = len(parts)
            while i < n:
                # Look ahead up to 3 tokens to find a keyword token
                joined = None
                max_look = min(i + 4, n)
                for k in range(i + 1, max_look):
                    token_k = parts[k].lower()
                    if any(kw in token_k for kw in pathway_keywords):
                        # Require at least two short tokens before the keyword OR presence of ' and ' in the keyword token
                        if (k - i) >= 2 or ' and ' in token_k or token_k.startswith('and '):
                            group = parts[i:k + 1]
                            # Join using spaces to form a single pathway name
                            joined = ' '.join(group)
                            i = k + 1
                        break

                if joined:
                    grouped.append(joined)
                else:
                    grouped.append(parts[i])
                    i += 1

            # Recreate string with semicolon separators between grouped entries
            return ';'.join(grouped)
        except Exception:
            # Fail safe: on any error, return original string
            return pathway_str

    def _split_pathways(self, pathway_str: str, separator: str = ';') -> List[str]:
        """Split pathway string into individual pathways and filter out excluded keywords"""
        if pd.isna(pathway_str) or pathway_str == "":
            return []
        
        # For HMDB data, preserve the | separator
        if separator == '|':
            cleaned = self._clean_pathway_string(pathway_str, preserve_separator='|')
        else:
            # Attempt to join semicolon-fragmented pathway names before cleaning
            try:
                grouped = self.group_semicolon_fragments(pathway_str)
            except Exception:
                grouped = pathway_str
            cleaned = self._clean_pathway_string(grouped)
            
        if not cleaned:
            return []
            
        pathways = [p.strip() for p in cleaned.split(separator) if p.strip()]
        
        # Filter out pathways containing excluded keywords
        filtered_pathways = []
        excluded_count = 0
        
        for pathway in pathways:
            # Check if pathway contains any excluded keywords (case-insensitive)
            pathway_lower = pathway.lower()
            should_exclude = False
            matched_keyword = None
            
            for keyword in self.exclude_pathway_keywords:
                if keyword in pathway_lower:
                    should_exclude = True
                    matched_keyword = keyword
                    break
            
            if should_exclude:
                excluded_count += 1
                #print(f"EXCLUDED PATHWAY: '{pathway}' (matched keyword: '{matched_keyword}')")
            else:
                # Drop obvious junk like a lone number or very short tokens
                if pathway.isdigit() or len(pathway) <= 2:
                    excluded_count += 1
                    #print(f"EXCLUDED PATHWAY: '{pathway}' (junk token)")
                    continue
                filtered_pathways.append(pathway)
        
        #if excluded_count > 0:
            #print(f"Filtered out {excluded_count} pathways from {len(pathways)} total")
        
        return filtered_pathways

    def _extract_biological_info(self, match: pd.Series) -> Dict[str, str]:
        """Extract biological locations and diseases from HMDB match"""
        # Extract biological locations (combine tissue and biofluid info)
        tissues = str(match.get('Tissues', '')).strip()
        
        # Check individual location columns
        locations = []
        location_columns = {
            'In_Saliva': 'Saliva',
            'In_Blood': 'Blood', 
            'In_Urine': 'Urine',
            'In_CSF': 'Cerebrospinal fluid',
            'In_Breast_Milk': 'Breast milk',
            'In_Feces': 'Feces',
            'In_Sweat': 'Sweat'
        }
        
        for col, location in location_columns.items():
            if str(match.get(col, '')).lower() in ['yes', 'true', '1']:
                locations.append(location)
        
        # Combine tissues and biofluid locations
        all_locations = []
        if tissues and tissues.lower() != 'nan':
            all_locations.append(f"Tissues: {tissues}")
        if locations:
            all_locations.append(f"Biofluids: {', '.join(locations)}")
        
        biological_locations = ' | '.join(all_locations)
        
        # Extract diseases
        diseases = str(match.get('Associated_Diseases', '')).strip()
        
        # Clean up empty or 'nan' values
        def clean_field(value):
            return '' if str(value).lower() in ['nan', 'none', ''] else str(value)
        
        return {
            'biological_locations': clean_field(biological_locations),
            'diseases': clean_field(diseases)
        }
    
    def map_hmdb_pathways(self, metabolite_row: pd.Series) -> Dict[str, Any]:
        """Map metabolite to HMDB pathways using available IDs"""
        result = {
            'HMDB_Pathways': '',
            'Enzymes_Accession': '',
            'Enzyme_Gene_name': '',
            'Reaction_Role': '',
            'Biological_Locations': '',
            'Associated_Diseases': '',
            'Transporters': '',
            'Transporter_Gene_name': '',
            'Transporter_Uniprot_ID': '',
            'hmdb_match_method': ''
        }

        hmdb_id = _safe_get_scalar(metabolite_row, 'HMDB_ID', '')
        inchikey = _safe_get_scalar(metabolite_row, 'InChIKey', '')
        name = _safe_get_scalar(metabolite_row, 'Name', '')

        matches = None
        match_method = ''

        try:
            # Strategy 1: Direct HMDB ID match
            # NOTE: Column gets renamed from 'HMDB' to 'HMDB_ID' after loading
            if hmdb_id and pd.notna(hmdb_id) and str(hmdb_id).strip() and 'HMDB_ID' in self.hmdb_df.columns:
                hmdb_matches = self.hmdb_df[self.hmdb_df['HMDB_ID'].str.strip() == str(hmdb_id).strip()]
                if not hmdb_matches.empty:
                    matches = hmdb_matches
                    match_method = 'HMDB_ID'

            # Strategy 2: InChIKey match
            elif inchikey and pd.notna(inchikey) and str(inchikey).strip() and 'InChIKey' in self.hmdb_df.columns:
                inchikey_matches = self.hmdb_df[self.hmdb_df['InChIKey'].str.strip() == str(inchikey).strip()]
                if not inchikey_matches.empty:
                    matches = inchikey_matches
                    match_method = 'InChIKey'

            # Strategy 3: Name match
            elif name and pd.notna(name) and str(name).strip():
                name_matches = self.hmdb_df[self.hmdb_df['Name'].str.lower().str.strip() == str(name).lower().strip()]
                if not name_matches.empty:
                    matches = name_matches
                    match_method = 'Name'

            if matches is not None and not matches.empty:
                match = matches.iloc[0]
                
                # Extract biological information
                bio_info = self._extract_biological_info(match)
                result['Biological_Locations'] = bio_info['biological_locations']
                result['Associated_Diseases'] = bio_info['diseases']
                
                # Extract pathway information - column is renamed to HMDB_Pathways after loading
                pathways = match.get('HMDB_Pathways', '')
                if pathways and str(pathways).strip():  # Check for non-empty strings
                    # Split by | separator and clean individual pathways
                    pathway_list = self._split_pathways(pathways, separator='|')
                    if pathway_list:
                        # Join back with | to preserve HMDB format
                        result['HMDB_Pathways'] = ' | '.join(pathway_list)

                # Extract enzyme information
                enzymes = match.get('Enzymes_Accession', '') or match.get('Enzymes', '')
                if enzymes:
                    result['Enzymes_Accession'] = str(enzymes)

                enzyme_genes = match.get('Enzyme_Gene_name', '') or match.get('Gene_Names', '')
                if enzyme_genes:
                    result['Enzyme_Gene_name'] = str(enzyme_genes)

                reaction_role = match.get('Reaction_Role', '')
                if reaction_role:
                    result['Reaction_Role'] = str(reaction_role)

                # Extract transporter information
                transporters = match.get('Transporters', '')
                if transporters:
                    result['Transporters'] = str(transporters)

                transporter_genes = match.get('Transporter_Gene_name', '')
                if transporter_genes:
                    result['Transporter_Gene_name'] = str(transporter_genes)

                transporter_uniprot = match.get('Transporter_Uniprot_ID', '')
                if transporter_uniprot:
                    result['Transporter_Uniprot_ID'] = str(transporter_uniprot)

                result['hmdb_match_method'] = match_method
                logger.debug(f"HMDB match for {name}: {match_method}")

        except Exception as e:
            logger.error(f"Error mapping HMDB pathways for {name}: {e}")

        return result

    def map_pathbank_pathways(self, metabolite_row: pd.Series) -> Dict[str, Any]:
        """Map metabolite to PathBank pathways"""
        result = {
            'PathBank_Pathways': '',
            'pathbank_match_method': ''
        }

        hmdb_id = _safe_get_scalar(metabolite_row, 'HMDB_ID', '')
        name = _safe_get_scalar(metabolite_row, 'Name', '')
        inchikey = _safe_get_scalar(metabolite_row, 'InChIKey', '')
        kegg_id = _safe_get_scalar(metabolite_row, 'KEGG_ID', '')

        matches = None
        match_method = ''

        try:
            # Strategy: Try all available IDs until a match is found
            # Priority: HMDB ID > InChI Key > KEGG ID > Name
            
            # Try HMDB ID first
            if hmdb_id and pd.notna(hmdb_id) and str(hmdb_id).strip() and 'HMDB ID' in self.pathbank_df.columns:
                hmdb_matches = self.pathbank_df[self.pathbank_df['HMDB ID'].str.strip() == str(hmdb_id).strip()]
                if not hmdb_matches.empty:
                    matches = hmdb_matches
                    match_method = 'HMDB ID'

            # If HMDB ID didn't match, try InChI Key
            if matches is None and inchikey and pd.notna(inchikey) and str(inchikey).strip() and 'InChI Key' in self.pathbank_df.columns:
                inchikey_matches = self.pathbank_df[self.pathbank_df['InChI Key'].str.strip() == str(inchikey).strip()]
                if not inchikey_matches.empty:
                    matches = inchikey_matches
                    match_method = 'InChI Key'

            # If still no match, try KEGG ID
            if matches is None and kegg_id and pd.notna(kegg_id) and str(kegg_id).strip() and 'KEGG ID' in self.pathbank_df.columns:
                kegg_matches = self.pathbank_df[self.pathbank_df['KEGG ID'].str.strip() == str(kegg_id).strip()]
                if not kegg_matches.empty:
                    matches = kegg_matches
                    match_method = 'KEGG ID'

            # Last resort: Name match
            if matches is None and name and pd.notna(name) and str(name).strip() and 'Metabolite Name' in self.pathbank_df.columns:
                name_matches = self.pathbank_df[self.pathbank_df['Metabolite Name'].str.lower().str.strip() == str(name).lower().strip()]
                if not name_matches.empty:
                    matches = name_matches
                    match_method = 'Metabolite Name'

            if matches is not None and not matches.empty:
                # Collect all pathways from multiple matches
                all_pathways = []
                for _, match in matches.iterrows():
                    pathways = match.get('PathBank_Pathways', '') 
                    if pathways:
                        pathway_list = self._split_pathways(pathways)
                        all_pathways.extend(pathway_list)

                if all_pathways:
                    # Remove duplicates while preserving order
                    unique_pathways = []
                    seen = set()
                    for pathway in all_pathways:
                        if pathway not in seen:
                            unique_pathways.append(pathway)
                            seen.add(pathway)
                    
                    result['PathBank_Pathways'] = '; '.join(unique_pathways)
                    result['pathbank_match_method'] = match_method
                    logger.debug(f"PathBank match for {name}: {match_method}")

        except Exception as e:
            logger.error(f"Error mapping PathBank pathways for {name}: {e}")

        return result

    def map_SMPDB_pathways(self, metabolite_row: pd.Series) -> Dict[str, Any]:
        """Map metabolite to SMPDB pathways"""
        result = {
            'SMPDB_Pathways': '',
            'SMPDB_match_method': ''
        }

        hmdb_id = _safe_get_scalar(metabolite_row, 'HMDB_ID', '')
        kegg_id = _safe_get_scalar(metabolite_row, 'KEGG_ID', '')
        name = _safe_get_scalar(metabolite_row, 'Name', '')

        matches = None
        match_method = ''

        try:
            # Strategy: Try all available IDs until a match is found
            # Priority: HMDB ID > KEGG ID > Name
            
            # Try HMDB ID first
            # NOTE: Database uses 'HMDB' column, not 'HMDB_ID'
            if hmdb_id and pd.notna(hmdb_id) and str(hmdb_id).strip() and 'HMDB' in self.SMPDB_df.columns:
                hmdb_matches = self.SMPDB_df[self.SMPDB_df['HMDB'].str.strip() == str(hmdb_id).strip()]
                if not hmdb_matches.empty:
                    matches = hmdb_matches
                    match_method = 'HMDB_ID'

            # If HMDB ID didn't match, try KEGG ID
            # NOTE: Database uses 'KEGG' column, not 'KEGG_ID'
            if matches is None and kegg_id and pd.notna(kegg_id) and str(kegg_id).strip() and 'KEGG' in self.SMPDB_df.columns:
                kegg_matches = self.SMPDB_df[self.SMPDB_df['KEGG'].str.strip() == str(kegg_id).strip()]
                if not kegg_matches.empty:
                    matches = kegg_matches
                    match_method = 'KEGG_ID'

            # Last resort: Name match
            if matches is None and name and pd.notna(name) and str(name).strip() and 'Name' in self.SMPDB_df.columns:
                name_matches = self.SMPDB_df[self.SMPDB_df['Name'].str.lower().str.strip() == str(name).lower().strip()]
                if not name_matches.empty:
                    matches = name_matches
                    match_method = 'Name'

            if matches is not None and not matches.empty:
                # Collect all pathways from multiple matches
                all_pathways = []
                for _, match in matches.iterrows():
                    # NOTE: Database uses 'SMP_Pathways' column, not 'SMPDB_Pathways'
                    pathways = match.get('SMP_Pathways', '') or match.get('Pathway', '')
                    if pathways:
                        pathway_list = self._split_pathways(pathways, separator='|')  # SMPDB uses | separator
                        all_pathways.extend(pathway_list)

                if all_pathways:
                    # Remove duplicates while preserving order
                    unique_pathways = []
                    seen = set()
                    for pathway in all_pathways:
                        if pathway not in seen:
                            unique_pathways.append(pathway)
                            seen.add(pathway)
                    
                    result['SMPDB_Pathways'] = '|'.join(unique_pathways)  # Keep SMPDB's original separator
                    result['SMPDB_match_method'] = match_method
                    logger.debug(f"SMPDB match for {name}: {match_method}")

        except Exception as e:
            logger.error(f"Error mapping SMPDB pathways for {name}: {e}")

        return result

    def map_wikipathways(self, metabolite_row: pd.Series) -> Dict[str, Any]:
        """Map metabolite to WikiPathways"""
        
        result = {
            'WikiPathways': '',
            'wikipathways_match_method': ''
        }

        kegg_id = _safe_get_scalar(metabolite_row, 'KEGG_ID', '')
        hmdb_id = _safe_get_scalar(metabolite_row, 'HMDB_ID', '')
        inchikey = _safe_get_scalar(metabolite_row, 'InChIKey', '')
        name = _safe_get_scalar(metabolite_row, 'Name', '')

        matches = None
        match_method = ''

        try:
            # Strategy 1: KEGG ID match
            if kegg_id and pd.notna(kegg_id) and str(kegg_id).strip() and 'KEGG' in self.wikipathways_df.columns:
                kegg_matches = self.wikipathways_df[self.wikipathways_df['KEGG'].str.strip() == str(kegg_id).strip()]
                if not kegg_matches.empty:
                    matches = kegg_matches
                    match_method = 'KEGG'

            # Strategy 2: HMDB ID match (only if not already matched)
            if matches is None and hmdb_id and pd.notna(hmdb_id) and str(hmdb_id).strip() and 'HMDB_ID' in self.wikipathways_df.columns:
                hmdb_matches = self.wikipathways_df[self.wikipathways_df['HMDB_ID'].str.strip() == str(hmdb_id).strip()]
                if not hmdb_matches.empty:
                    matches = hmdb_matches
                    match_method = 'HMDB_ID'

            # Strategy 3: InChIKey match (only if not already matched)
            if matches is None and inchikey and pd.notna(inchikey) and str(inchikey).strip() and 'InChIKey' in self.wikipathways_df.columns:
                inchikey_matches = self.wikipathways_df[self.wikipathways_df['InChIKey'].str.strip() == str(inchikey).strip()]
                if not inchikey_matches.empty:
                    matches = inchikey_matches
                    match_method = 'InChIKey'

            # Strategy 4: Name match (only if not already matched)
            if matches is None and name and pd.notna(name) and str(name).strip() and 'name' in self.wikipathways_df.columns:
                name_matches = self.wikipathways_df[self.wikipathways_df['name'].str.lower().str.strip() == str(name).lower().strip()]
                if not name_matches.empty:
                    matches = name_matches
                    match_method = 'name'

            if matches is not None and not matches.empty:
                # Collect all pathways from multiple matches
                all_pathways = []
                for _, match in matches.iterrows():
                    # NOTE: Database uses 'Pathways' column, not 'WikiPathways'
                    pathways = match.get('Pathways', '') or match.get('WikiPathways', '')
                    if pathways:
                        pathway_list = self._split_pathways(pathways)
                        all_pathways.extend(pathway_list)

                if all_pathways:
                    # Remove duplicates while preserving order
                    unique_pathways = []
                    seen = set()
                    for pathway in all_pathways:
                        if pathway not in seen:
                            unique_pathways.append(pathway)
                            seen.add(pathway)
                    
                    result['WikiPathways'] = '; '.join(unique_pathways)
                    result['wikipathways_match_method'] = match_method
                    logger.debug(f"WikiPathways match for {name}: {match_method}")

        except Exception as e:
            logger.error(f"Error mapping WikiPathways for {name}: {e}")

        return result

    def map_metabolika_pathways(self, metabolite_row: pd.Series) -> Dict[str, Any]:
        """
        Map metabolite to Metabolika pathways from data cleaner cache
        """
        name = str(_safe_get_scalar(metabolite_row, 'Name', '')).strip()
        result = {
            'Metabolika_Pathways': '',
            'metabolika_match_method': ''
        }

        try:
            if self.metabolika_cache.empty:
                logger.debug(f"Metabolika cache is empty for {name}")
                return result

            logger.debug(f"Searching Metabolika cache for: {name}")

            # Search by exact name match
            matches = self.metabolika_cache[
                self.metabolika_cache['Name'].str.lower().str.strip() == name.lower().strip()
            ]

            if not matches.empty:
                match_row = matches.iloc[0]
                pathways_str = str(match_row.get('Pathways', '')).strip()
                
                if pathways_str and pathways_str != 'nan':
                    result['Metabolika_Pathways'] = pathways_str
                    result['metabolika_match_method'] = 'exact_name'
                    logger.debug(f"Metabolika exact name match for {name}")
                else:
                    logger.debug(f"Metabolika match found but no pathways for {name}")
            else:
                logger.debug(f"No Metabolika match for {name}")

        except Exception as e:
            logger.error(f"Error mapping Metabolika pathways for {name}: {e}")

        return result

    def map_reactome_pathways(self, metabolite_row: pd.Series, cache: ReactomeKEGGCache) -> Dict[str, Any]:
        """Map metabolite to Reactome pathways using ChEBI ID"""
        result = {
            'Reactome_Pathways': '',
            'Reactome_Disease': '',
            'reactome_match_method': ''
        }

        chebi_id = _safe_get_scalar(metabolite_row, 'ChEBI_ID', '')
        name = _safe_get_scalar(metabolite_row, 'Name', 'Unknown')

        if not chebi_id or pd.isna(chebi_id) or str(chebi_id).strip() == '':
            return result

        try:
            # Normalize ChEBI ID (strip "CHEBI:" prefix)
            chebi_id = str(chebi_id).replace("CHEBI:", "").strip()
            
            # Get organism taxID
            taxid = SPECIES_TAXID.get(self.organism, "9606")
            
            # Check cache
            cached = cache.get(chebi_id, 'reactome', self.organism)
            if cached is not None:
                pathways = [p for p in cached if not p.get('isInDisease', False)]
                disease_pathways = [p for p in cached if p.get('isInDisease', False)]
            else:
                # Query Reactome API
                url = f"{REACTOME_BASE}/mapping/ChEBI/{chebi_id}/pathways?speciesId={taxid}"
                data = fetch_json(url)
                
                if isinstance(data, list):
                    cache.set(chebi_id, 'reactome', self.organism, data)
                    pathways = [p for p in data if not p.get('isInDisease', False)]
                    disease_pathways = [p for p in data if p.get('isInDisease', False)]
                else:
                    pathways = []
                    disease_pathways = []

            # Format pathways
            if pathways:
                pathway_names = [p.get('displayName', '') for p in pathways if p.get('displayName')]
                result['Reactome_Pathways'] = '; '.join(pathway_names)
                result['reactome_match_method'] = 'ChEBI_ID'

            if disease_pathways:
                disease_names = [p.get('displayName', '') for p in disease_pathways if p.get('displayName')]
                result['Reactome_Disease'] = '; '.join(disease_names)

        except urllib.error.HTTPError as e:
            if e.code != 404:  # 404 is expected for metabolites not in Reactome
                logger.warning(f"Reactome HTTP error for {name} (ChEBI {chebi_id}): {e.code}")
        except Exception as e:
            logger.error(f"Error mapping Reactome pathways for {name}: {e}")

        return result

    def map_kegg_pathways(self, metabolite_row: pd.Series, cache: ReactomeKEGGCache) -> Dict[str, Any]:
        """Map metabolite to KEGG pathways using KEGG Compound ID"""
        result = {
            'KEGG_Pathways': '',
            'kegg_match_method': ''
        }

        kegg_id = _safe_get_scalar(metabolite_row, 'KEGG_ID', '')
        name = _safe_get_scalar(metabolite_row, 'Name', 'Unknown')

        if not kegg_id or pd.isna(kegg_id) or str(kegg_id).strip() == '':
            return result

        try:
            kegg_id = str(kegg_id).strip()
            org_code = KEGG_ORGANISMS.get(self.organism, 'hsa')
            
            # Check cache
            cached = cache.get(kegg_id, 'kegg', self.organism)
            if cached is not None:
                organism_pathways = cached
            else:
                # Step 1: Get generic pathways from compound entry
                url = f"{KEGG_BASE}/get/{kegg_id}"
                req = urllib.request.Request(url, headers={"User-Agent": "metabolite-mapper/1.0"})
                with urllib.request.urlopen(req, timeout=15) as resp:
                    result_text = resp.read().decode('utf-8')
                
                # Extract PATHWAY lines. KEGG uses a 'PATHWAY' label on the first
                # line and indents subsequent pathway lines. Collect both formats:
                # - Lines that start with 'PATHWAY' where the second token is the map id
                # - Indented lines that start directly with 'mapXXXXX'
                generic_pathways = []
                for line in result_text.split('\n'):
                    if line.startswith('PATHWAY'):
                        parts = line.split()
                        if len(parts) >= 3:
                            pathway_id = parts[1]
                            pathway_name = ' '.join(parts[2:])
                            generic_pathways.append((pathway_id, pathway_name))
                        continue

                    # Also capture indented continuation lines like '    map01100   Metabolic pathways'
                    m = re.match(r"^\s*(map\d+)\s+(.*)$", line)
                    if m:
                        pathway_id = m.group(1)
                        pathway_name = m.group(2).strip()
                        generic_pathways.append((pathway_id, pathway_name))
                
                # Step 2: Convert to organism-specific pathway IDs
                organism_pathways = []
                for generic_pathway_id, generic_name in generic_pathways:
                    org_pathway_id = generic_pathway_id.replace('map', org_code)
                    
                    try:
                        url = f"{KEGG_BASE}/get/{org_pathway_id}"
                        req = urllib.request.Request(url, headers={"User-Agent": "metabolite-mapper/1.0"})
                        with urllib.request.urlopen(req, timeout=15) as resp:
                            org_result = resp.read().decode('utf-8')
                        
                        # Extract NAME field
                        org_name = None
                        for line in org_result.split('\n'):
                            if line.startswith('NAME'):
                                org_name = line.replace('NAME', '').strip()
                                # Remove species suffix (e.g., " - Homo sapiens (human)")
                                org_name = re.sub(r'\s*-\s*(Homo sapiens|Mus musculus|Rattus norvegicus)\s*\([^)]+\)\s*$', '', org_name)
                                break
                        
                        if org_name:
                            organism_pathways.append((org_pathway_id, org_name))
                        else:
                            organism_pathways.append((org_pathway_id, generic_name))
                    except urllib.error.HTTPError:
                        pass  # Pathway not available for this organism
                    
                    time.sleep(0.1)  # Rate limiting
                
                # Cache the result
                cache.set(kegg_id, 'kegg', self.organism, organism_pathways)

            # Format pathways
            if organism_pathways:
                pathway_names = [name for pathway_id, name in organism_pathways]
                result['KEGG_Pathways'] = '; '.join(pathway_names)
                result['kegg_match_method'] = 'KEGG_ID'

        except urllib.error.HTTPError as e:
            if e.code != 404 and e.code != 400:  # 404/400 expected for some compounds
                logger.warning(f"KEGG HTTP error for {name} (KEGG {kegg_id}): {e.code}")
        except Exception as e:
            logger.error(f"Error mapping KEGG pathways for {name}: {e}")

        return result

    def map_metabolite_pathways(self, metabolite_row: pd.Series) -> Dict[str, Any]:
        """Complete pathway mapping for a single metabolite"""
        name = _safe_get_scalar(metabolite_row, 'Name', 'Unknown')
        logger.debug(f"Mapping pathways for: {name}")

        # Initialize result with original metabolite data
        result = metabolite_row.to_dict()

        # Initialize cache for Reactome/KEGG
        cache = ReactomeKEGGCache()

        # Map pathways from each database (run Reactome & KEGG concurrently)
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                'hmdb': executor.submit(self.map_hmdb_pathways, metabolite_row),
                'pathbank': executor.submit(self.map_pathbank_pathways, metabolite_row),
                'SMPDB': executor.submit(self.map_SMPDB_pathways, metabolite_row),
                'wikipathways': executor.submit(self.map_wikipathways, metabolite_row),
                'metabolika': executor.submit(self.map_metabolika_pathways, metabolite_row),
                'reactome': executor.submit(self.map_reactome_pathways, metabolite_row, cache),
                'kegg': executor.submit(self.map_kegg_pathways, metabolite_row, cache)
            }
            
            # Collect results
            hmdb_result = futures['hmdb'].result()
            pathbank_result = futures['pathbank'].result()
            SMPDB_result = futures['SMPDB'].result()
            wikipathways_result = futures['wikipathways'].result()
            metabolika_result = futures['metabolika'].result()
            reactome_result = futures['reactome'].result()
            kegg_result = futures['kegg'].result()

        # if metabolika_result.get('Metabolika_Pathways'):
        #     logger.info(f"[CACHE SUCCESS] Metabolika: Found pathways for {name}: {metabolika_result['Metabolika_Pathways'][:100]}...")
        # else:
        #     logger.info(f"[CACHE MISS] Metabolika: No pathways found for {name}")

        # Update result with pathway information
        result.update(hmdb_result)
        result.update(pathbank_result)
        result.update(SMPDB_result)
        result.update(wikipathways_result)
        result.update(metabolika_result)
        result.update(reactome_result)
        result.update(kegg_result)

        # If Metabolika provided pathways, ensure cache is updated so future runs reuse it
        try:
            metabolika_pw = metabolika_result.get('Metabolika_Pathways')
            if metabolika_pw:
                # Normalize name for existence check
                name_norm = str(name).strip()
                exists = False
                if not self.metabolika_cache.empty and 'Name' in self.metabolika_cache.columns:
                    exists = (self.metabolika_cache['Name'].str.lower().str.strip() == name_norm.lower()).any()

                if not exists:
                    # Append new cache row and persist
                    new_row = {'Name': name_norm, 'Pathways': metabolika_pw}
                    try:
                        if self.metabolika_cache.empty:
                            self.metabolika_cache = pd.DataFrame([new_row])
                        else:
                            self.metabolika_cache = pd.concat([self.metabolika_cache, pd.DataFrame([new_row])], ignore_index=True)
                        self._save_metabolika_cache()
                        logger.info(f"[CACHE UPDATE] Added Metabolika cache entry for: {name_norm}")
                    except Exception as e:
                        logger.warning(f"Failed to append to Metabolika cache for {name_norm}: {e}")
        except Exception:
            # Non-critical: don't break mapping on cache persistence errors
            pass

        # Create combined pathways list for network analysis compatibility
        all_pathways = []
        
        # Collect pathways from all sources (including Reactome and KEGG)
        for source_col in ['HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways', 
                           'Metabolika_Pathways', 'Reactome_Pathways', 'KEGG_Pathways']:
            pathways_str = result.get(source_col, '')
            if pathways_str:
                if source_col == 'SMPDB_Pathways':
                    # SMPDB uses | separator
                    pathway_list = self._split_pathways(pathways_str, separator='|')
                else:
                    # For combining sources, directly split on semicolon WITHOUT fragment grouping
                    # Fragment grouping was causing separate pathways to be concatenated
                    # (e.g., "Purine Metabolism" + "Purine Metabolism And Related Disorders" 
                    #  would become "Purine Metabolism Purine Metabolism And Related Disorders")
                    cleaned = self._clean_pathway_string(pathways_str)
                    pathway_list = [p.strip() for p in cleaned.split(';') if p.strip()] if cleaned else []
                all_pathways.extend(pathway_list)

        # Normalize pathway names and remove duplicates (case-insensitive)
        # This handles: "Nicotine and Nicotinamide Metabolims" vs "Nicotine and Nicotinamide metabolims"
        unique_pathways = []
        seen_normalized = {}  # Maps normalized name to canonical (first seen) name
        
        for pathway in all_pathways:
            normalized = normalize_pathway_name(pathway)
            normalized_lower = normalized.lower()
            
            # Check if we've seen this pathway name (case-insensitive)
            if normalized_lower not in seen_normalized:
                # First time seeing this pathway - use the normalized name
                unique_pathways.append(normalized)
                seen_normalized[normalized_lower] = normalized

        # Step 1: Filter out non-specific and excluded pathways (plant/microbial/superpathway)
        def _is_excluded_non_specific(name_lower: str) -> bool:
            # Keyword exclusions
            for kw in self.exclude_pathway_keywords:
                if kw in name_lower:
                    return True
            # Non-specific patterns
            for pat in NON_SPECIFIC_PATHWAY_PATTERNS:
                if pat in name_lower:
                    return True
            # Very short or numeric-like
            if len(name_lower.strip()) <= 2:
                return True
            if name_lower.replace('.', '').replace('-', '').isdigit():
                return True
            return False

        filtered_paths: List[str] = []
        for p in unique_pathways:
            if not _is_excluded_non_specific(p.lower()):
                filtered_paths.append(p)

        # Step 2: Move disease-labeled pathways from pathways list to Associated_Diseases
        moved_diseases: List[str] = []
        retained_paths: List[str] = []
        for p in filtered_paths:
            if is_disease_pathway(p):
                moved_diseases.append(normalize_disease_pathway_name(p))
            else:
                retained_paths.append(p)
        
        # Step 3: Deduplicate pathways (case-insensitive) before storing
        unique_paths: List[str] = []
        seen_paths_lower = set()
        for p in retained_paths:
            p_lower = p.lower()
            if p_lower not in seen_paths_lower:
                seen_paths_lower.add(p_lower)
                unique_paths.append(p)

        # Store as list for programmatic use (network analysis expects this)
        result['All_Pathways'] = unique_paths if unique_paths else []
        
        # Create semicolon-separated string for Excel display in separate column
        result['All_Pathways_Display'] = '; '.join(unique_paths) if unique_paths else ''

        # Merge Reactome_Disease and moved disease pathways into Associated_Diseases
        associated_diseases = []
        
        # Get diseases from HMDB (if already in Associated_Diseases)
        hmdb_diseases = result.get('Associated_Diseases', '')
        if hmdb_diseases and str(hmdb_diseases).strip() and str(hmdb_diseases).lower() != 'nan':
            # Split by | or ; separator
            if '|' in str(hmdb_diseases):
                associated_diseases.extend([d.strip() for d in str(hmdb_diseases).split('|') if d.strip()])
            else:
                associated_diseases.extend([d.strip() for d in str(hmdb_diseases).split(';') if d.strip()])
        
        # Get diseases from Reactome_Disease column
        reactome_diseases = result.get('Reactome_Disease', '')
        if reactome_diseases and str(reactome_diseases).strip() and str(reactome_diseases).lower() != 'nan':
            # Split by ; separator (Reactome uses ; )
            associated_diseases.extend([d.strip() for d in str(reactome_diseases).split(';') if d.strip()])
        
        # Include diseases moved from pathway filtering
        if moved_diseases:
            associated_diseases.extend(moved_diseases)

        # Remove duplicates while preserving order
        unique_diseases = []
        seen_diseases = set()
        for disease in associated_diseases:
            disease_lower = disease.lower()
            if disease_lower not in seen_diseases:
                unique_diseases.append(disease)
                seen_diseases.add(disease_lower)
        
        # Update Associated_Diseases with merged data (use | separator for consistency with HMDB)
        result['Associated_Diseases'] = ' | '.join(unique_diseases) if unique_diseases else ''

        # Summary statistics
        pathway_sources = []
        for source, col in [('HMDB', 'HMDB_Pathways'), ('PathBank', 'PathBank_Pathways'), 
                           ('SMPDB', 'SMPDB_Pathways'), ('WikiPathways', 'WikiPathways'),
                           ('Metabolika', 'Metabolika_Pathways'), ('Reactome', 'Reactome_Pathways'),
                           ('KEGG', 'KEGG_Pathways')]:
            if result.get(col):
                pathway_sources.append(source)

        result['pathway_sources'] = ', '.join(pathway_sources) if pathway_sources else 'None'
        result['total_pathways'] = len(unique_pathways)

        logger.debug(f"Pathway mapping complete for {name}: {len(unique_pathways)} pathways from {len(pathway_sources)} sources")

        return result

    def _is_cancelled(self):
        """Check if processing should be cancelled"""
        if self.cancel_flag is None:
            return False
        try:
            # cancel_flag can be a callable or a boolean-like attribute
            if callable(self.cancel_flag):
                return self.cancel_flag()
            return bool(self.cancel_flag)
        except Exception:
            return False

    def process_metabolites(self, input_df: pd.DataFrame) -> pd.DataFrame:
        """Process all metabolites for pathway mapping with parallel execution"""
        total_metabolites = len(input_df)
        logger.info(f"Starting parallel pathway mapping for {total_metabolites} metabolites (workers: {self.max_metabolite_workers})")

        # Check for cancellation before starting
        if self._is_cancelled():
            logger.info("Processing cancelled before start")
            return pd.DataFrame()

        # Convert DataFrame rows to list for easier parallel processing
        metabolite_rows = [(idx, row) for idx, row in input_df.iterrows()]
        
        # Thread-safe counters for progress tracking
        completed_count = 0
        last_logged_progress = 0
        cancelled = False
        
        def process_single_metabolite(idx_row_tuple):
            """Process a single metabolite (for parallel execution)"""
            nonlocal completed_count, last_logged_progress, cancelled
            
            # Check for cancellation
            if self._is_cancelled():
                cancelled = True
                return None
            
            idx, row = idx_row_tuple
            try:
                result = self.map_metabolite_pathways(row)
                
                # Thread-safe progress reporting
                with self.progress_lock:
                    # Check cancellation again with lock held
                    if self._is_cancelled():
                        cancelled = True
                        return result
                    
                    completed_count += 1
                    progress_percent = (completed_count / total_metabolites) * 100
                    current_progress_threshold = int(progress_percent // 10) * 10
                    
                    # Log to console every 10%
                    if current_progress_threshold > last_logged_progress or completed_count == total_metabolites:
                        progress_msg = f"Progress: {completed_count}/{total_metabolites} ({progress_percent:.1f}%)"
                        logger.info(progress_msg)
                        last_logged_progress = current_progress_threshold
                    
                    # Always call progress callback for GUI (keeps progress bar smooth)
                    if self.progress_callback:
                        detailed_msg = f"🛣️ Processing {completed_count}/{total_metabolites} ({progress_percent:.1f}%)"
                        self.progress_callback(completed_count, total_metabolites, detailed_msg)
                
                return result
                    
            except Exception as e:
                logger.error(f"Error processing metabolite {row.get('Name', 'Unknown')}: {e}")
                # Add error result
                error_result = row.to_dict()
                error_result.update({
                    'HMDB_Pathways': '', 'PathBank_Pathways': '', 'SMPDB_Pathways': '', 'WikiPathways': '',
                    'All_Pathways': [], 'All_Pathways_Display': '', 'pathway_sources': 'ERROR', 'total_pathways': 0
                })
                return error_result
        
        # Process metabolites in parallel using ThreadPoolExecutor
        results = []
        with ThreadPoolExecutor(max_workers=self.max_metabolite_workers) as executor:
            # Submit all metabolites for processing
            futures = [executor.submit(process_single_metabolite, item) for item in metabolite_rows]
            
            # Collect results as they complete, checking for cancellation
            for future in futures:
                if self._is_cancelled():
                    logger.info("Processing cancelled - stopping executor")
                    executor.shutdown(wait=False, cancel_futures=True)
                    break
                result = future.result()
                if result is not None:
                    results.append(result)

        if self._is_cancelled():
            logger.info(f"Pathway mapping cancelled after {len(results)} metabolites")
            # Return partial results or empty DataFrame
            if results:
                return pd.DataFrame(results)
            return pd.DataFrame()

        # Convert results to DataFrame
        results_df = pd.DataFrame(results)

        logger.info(f"Pathway mapping completed for {len(results_df)} metabolites")
        return results_df

    def create_pathway_summary(self, mapped_df: pd.DataFrame) -> pd.DataFrame:
        """Create pathway summary compatible with existing network analysis"""
        logger.info("Creating pathway summary for network analysis...")

        # Filter metabolites with pathways
        metabolites_with_pathways = mapped_df[mapped_df['total_pathways'] > 0].copy()
        
        if metabolites_with_pathways.empty:
            logger.warning("No metabolites with pathways found")
            return pd.DataFrame(columns=['Pathway', 'Metabolites'])

        # Create pathway-to-metabolites mapping
        pathway_metabolites = {}
        
        for _, row in metabolites_with_pathways.iterrows():
            metabolite_name = row['Name']
            
            # Handle both list (All_Pathways) and string (All_Pathways_Display) formats
            pathways = row.get('All_Pathways')
            if pathways is None or not isinstance(pathways, list):
                # Fallback to string format and split it
                pathways_str = row.get('All_Pathways_Display', '') or row.get('All_Pathways', '')
                if isinstance(pathways_str, str) and pathways_str:
                    # Try multiple separators for robust parsing
                    if '; ' in pathways_str:
                        pathways = self._split_pathways(pathways_str, separator=';')
                    elif ' | ' in pathways_str:
                        pathways = self._split_pathways(pathways_str, separator='|')
                    elif '|' in pathways_str:
                        pathways = self._split_pathways(pathways_str, separator='|')
                    else:
                        pathways = self._split_pathways(pathways_str)
                else:
                    pathways = []
            
            # Additional check for individual database pathway columns if All_Pathways_List failed
            if not pathways:
                individual_pathways = []
                
                # Check HMDB pathways
                hmdb_pathways = row.get('HMDB_Pathways', '')
                if hmdb_pathways:
                    individual_pathways.extend(self._split_pathways(hmdb_pathways, separator='|'))
                
                # Check PathBank pathways  
                pathbank_pathways = row.get('PathBank_Pathways', '')
                if pathbank_pathways:
                    individual_pathways.extend(self._split_pathways(pathbank_pathways, separator=';'))
                
                # Check SMPDB pathways
                smpdb_pathways = row.get('SMPDB_Pathways', '')
                if smpdb_pathways:
                    individual_pathways.extend(self._split_pathways(smpdb_pathways, separator='|'))
                
                # Check WikiPathways
                wiki_pathways = row.get('WikiPathways', '')
                if wiki_pathways:
                    individual_pathways.extend(self._split_pathways(wiki_pathways, separator=';'))
                
                pathways = individual_pathways
            
            if isinstance(pathways, list):
                for pathway in pathways:
                    if pathway.strip():
                        if pathway not in pathway_metabolites:
                            pathway_metabolites[pathway] = []
                        pathway_metabolites[pathway].append(metabolite_name)

        # Convert to DataFrame
        pathway_summary = []
        for pathway, metabolites in pathway_metabolites.items():
            pathway_summary.append({
                'Pathway': pathway,
                'Metabolites': ','.join(sorted(set(metabolites))),  # Remove duplicates
                'Metabolite_Count': len(set(metabolites))
            })

        pathway_df = pd.DataFrame(pathway_summary)
        pathway_df = pathway_df.sort_values('Metabolite_Count', ascending=False)

        logger.info(f"Created pathway summary with {len(pathway_df)} pathways")
        return pathway_df

    def filter_metabolites_by_pvalue(self, input_df: pd.DataFrame, pvalue_threshold: float = 0.05) -> pd.DataFrame:
        """
        Filter metabolites by p-value threshold BEFORE pathway annotation begins.
        This ensures only significant metabolites are included in pathway calculations.
        
        Parameters
        ----------
        input_df : pd.DataFrame
            Input dataframe with metabolites
        pvalue_threshold : float
            P-value threshold for filtering (default: 0.05)
            Metabolites with p-value > threshold will be excluded
        
        Returns
        -------
        pd.DataFrame
            Filtered dataframe containing only significant metabolites
        """
        initial_count = len(input_df)
        
        # Check if pvalue column exists - fail fast per new requirement
        if 'pvalue' not in input_df.columns:
            raise ValueError("P-value filtering requested but 'pvalue' column not found in input DataFrame.")
        
        # Filter metabolites: keep only those with pvalue <= threshold
        # Also filter out NaN p-values
        filtered_df = input_df[
            (input_df['pvalue'].notna()) & 
            (input_df['pvalue'] <= pvalue_threshold)
        ].copy()
        
        filtered_count = len(filtered_df)
        removed_count = initial_count - filtered_count
        
        logger.info(f"="*60)
        logger.info(f"P-VALUE FILTERING (Threshold: {pvalue_threshold})")
        logger.info(f"="*60)
        logger.info(f"Initial metabolites: {initial_count}")
        logger.info(f"Significant metabolites (p ≤ {pvalue_threshold}): {filtered_count} ({filtered_count/initial_count*100:.1f}%)")
        logger.info(f"Non-significant metabolites removed: {removed_count} ({removed_count/initial_count*100:.1f}%)")
        logger.info(f"="*60)
        
        if filtered_count == 0:
            raise ValueError(
                f"No metabolites passed p-value threshold of {pvalue_threshold}. "
                f"All {initial_count} metabolites were filtered out. "
                f"Consider using a less stringent threshold or checking your data."
            )
        
        return filtered_df
    
    def validate_and_process_input_data(self, input_df: pd.DataFrame) -> pd.DataFrame:
        """
        Validate and process input data to ensure required columns for network analysis
        Required: Name, pvalue, log2FC (or FC to calculate log2FC)
        """
        logger.info("Validating input data for network analysis requirements...")

        # If the DataFrame already contains canonical columns provided by the GUI
        # (i.e. the user verified columns and they were renamed to 'pvalue' and 'log2FC'),
        # skip auto-detection and noisy informational prints. This prevents duplicate
        # terminal logs when the GUI has already normalized the input.
        if 'pvalue' in input_df.columns and 'log2FC' in input_df.columns:
            logger.info("Input already has canonical 'pvalue' and 'log2FC' columns - skipping auto-detection")
            return input_df
        
        # Strict requirement: 'Name' must be present and is expected to be
        # provided by the GUI Verify Columns dialog. Do NOT auto-detect or
        # rename other columns here — fail fast so the user can correct input.
        if 'Name' not in input_df.columns:
            raise ValueError(
                "Input file must contain a 'Name' column containing the feature identifier. "
                "Please use the GUI 'Verify Columns' dialog to assign the Feature ID before annotation."
            )
        
        # Check for required statistical columns - validate BOTH before raising errors
        missing_columns = []
        
        # Check for pvalue column with PRIORITY: adjusted p-value > regular p-value
        # Build case-insensitive column lookup map
        col_lower_map = {str(c).lower(): c for c in input_df.columns}
        
        # Priority 1: Adjusted p-values (preferred) - exact match
        adj_pvalue_cols = ['adj_p', 'adjp', 'adj.p', 'p_adj', 'padj', 'p.adj', 
                          'adjusted_p', 'adjusted_pvalue', 'p_adjusted', 'fdr', 
                          'q_value', 'qvalue', 'control_vs_tbi_adj_p']
        # Priority 2: Regular p-values (fallback) - exact match
        reg_pvalue_cols = ['pvalue', 'p_value', 'p-value', 'p.value', 'pval', 
                          'p_val', 'p', 'control_vs_tbi_p']
        
        # Suffixes for detecting prefixed columns like PD_vs_Control_adj_p
        adj_p_suffixes = ['_adj_p', '_adjp', '_p_adj', '_padj', '_fdr', '_qvalue', '_q_value']
        reg_p_suffixes = ['_pvalue', '_p_value', '_pval', '_p_val']
        
        pvalue_col = None
        pvalue_source = None
        
        # First try adjusted p-values (highest priority) - exact match
        for col in adj_pvalue_cols:
            if col.lower() in col_lower_map:
                pvalue_col = col_lower_map[col.lower()]
                pvalue_source = 'adjusted'
                break
        
        # Try suffix match for adjusted p-values (e.g., PD_vs_Control_adj_p)
        if pvalue_col is None:
            for col_lower, col_orig in col_lower_map.items():
                for suffix in adj_p_suffixes:
                    if col_lower.endswith(suffix):
                        pvalue_col = col_orig
                        pvalue_source = 'adjusted'
                        break
                if pvalue_col:
                    break
        
        # If no adjusted p-value found, try regular p-values
        if pvalue_col is None:
            for col in reg_pvalue_cols:
                if col.lower() in col_lower_map:
                    pvalue_col = col_lower_map[col.lower()]
                    pvalue_source = 'regular'
                    break
        
        # Try suffix match for regular p-values
        if pvalue_col is None:
            for col_lower, col_orig in col_lower_map.items():
                for suffix in reg_p_suffixes:
                    if col_lower.endswith(suffix):
                        pvalue_col = col_orig
                        pvalue_source = 'regular'
                        break
                if pvalue_col:
                    break
        
        # Check for log2FC or FC column with case-insensitive matching
        log2fc_cols = ['log2fc', 'log2_fc', 'logfc', 'log_fc', 
                      'control_vs_tbi_log2fc', 'log2foldchange']
        fc_cols = ['fc', 'fold_change', 'foldchange', 'control_vs_tbi_fc']
        
        # Suffixes for detecting prefixed fold change columns like PD_vs_Control_log2FC
        log2fc_suffixes = ['_log2fc', '_log2_fc', '_logfc', '_log_fc', '_log2foldchange']
        fc_suffixes = ['_fc', '_fold_change', '_foldchange']
        
        log2fc_col = None
        fc_col = None
        
        # Check for existing log2FC column (case-insensitive) - exact match first
        for col in log2fc_cols:
            if col.lower() in col_lower_map:
                log2fc_col = col_lower_map[col.lower()]
                break
        
        # Try suffix match for log2FC columns (e.g., PD_vs_Control_log2FC)
        if log2fc_col is None:
            for col_lower, col_orig in col_lower_map.items():
                for suffix in log2fc_suffixes:
                    if col_lower.endswith(suffix):
                        log2fc_col = col_orig
                        break
                if log2fc_col:
                    break
        
        # Check for FC column if no log2FC found (case-insensitive) - exact match first
        if log2fc_col is None:
            for col in fc_cols:
                if col.lower() in col_lower_map:
                    fc_col = col_lower_map[col.lower()]
                    break
        
        # Try suffix match for FC columns (e.g., PD_vs_Control_FC)
        if log2fc_col is None and fc_col is None:
            for col_lower, col_orig in col_lower_map.items():
                for suffix in fc_suffixes:
                    if col_lower.endswith(suffix):
                        fc_col = col_orig
                        break
                if fc_col:
                    break
        
        # Never use negative log10 p-value columns
        neg_log_pvalue_cols = ['neg_log10_adj_p', 'neg_log10_p', '-log10_p', 
                               '_neg_log10_adj_p']
        for col in input_df.columns:
            col_lower = str(col).lower()
            for neg_col in neg_log_pvalue_cols:
                if neg_col.lower() in col_lower:
                    print(f"⚠️ Warning: Ignoring negative log p-value column '{col}' - not suitable for pathway analysis")
        
        # Collect missing column information
        if pvalue_col is None:
            missing_columns.append({
                'type': 'P-value',
                'names': 'Adjusted p-value (preferred): adj_p, p_adj, padj, fdr, q_value\n  Regular p-value (fallback): pvalue, p_value, p-value',
                'description': 'statistical significance values for metabolite changes (adjusted p-values are preferred over raw p-values)'
            })
        
        if log2fc_col is None and fc_col is None:
            missing_columns.append({
                'type': 'Fold change',
                'names': 'log2FC, log2_FC, Log2FC, Log2_FC, logFC, log_FC, FC, Fold_Change, FoldChange, fold_change',
                'description': 'magnitude of metabolite changes (fold change or log2 fold change values)'
            })
        
        # Raise comprehensive error if any columns are missing
        if missing_columns:
            error_message = "❌ MISSING REQUIRED COLUMNS: Pathway analysis requires both p-value and fold change columns.\n\n"
            
            for i, col_info in enumerate(missing_columns, 1):
                error_message += f"MISSING {i}: {col_info['type']} column\n"
                error_message += f"REQUIRED: One of these column names:\n"
                error_message += f"  {col_info['names']}\n"
                error_message += f"DESCRIPTION: This column should contain {col_info['description']}.\n\n"
            
            error_message += "SOLUTION:\n"
            error_message += "1. Add the missing column(s) to your Excel file\n"
            error_message += "2. Ensure column names match exactly (case-insensitive matching is supported)\n"
            error_message += "3. Use 'Verify Columns' button before analysis to check your columns\n"
            error_message += "4. Reload the file and try again\n\n"
            error_message += f"Your current file columns: {', '.join(input_df.columns.tolist())}"
            
            raise ValueError(error_message)
        
        # Log which p-value column was used
        if pvalue_source == 'adjusted':
            print(f"✅ Using adjusted p-value column: '{pvalue_col}'")
        elif pvalue_col == 'pvalue':
            # Don't warn if the column is already named 'pvalue' - it might have been renamed from adj_p upstream
            print(f"✓ Using p-value column: '{pvalue_col}'")
        else:
            print(f"⚠️ Using regular p-value column: '{pvalue_col}' (no adjusted p-value found)")
        
        # Process valid columns - rename to standard names
        if pvalue_col != 'pvalue':
            input_df = input_df.rename(columns={pvalue_col: 'pvalue'})
        logger.info(f"✓ Found p-value column: {pvalue_col}")
        
        if log2fc_col is not None:
            if log2fc_col != 'log2FC':
                input_df = input_df.rename(columns={log2fc_col: 'log2FC'})
            logger.info(f"✓ Found log2FC column: {log2fc_col}")
        elif fc_col is not None:
            logger.info(f"✓ Found FC column: {fc_col}. Converting to log2FC...")
            # Convert FC to log2FC: log2FC = log2(FC)
            import numpy as np
            input_df['log2FC'] = np.log2(input_df[fc_col].replace(0, np.nan))  # Avoid log(0)
            logger.info(f"✓ Converted {fc_col} to log2FC")
        
        # Log summary of validation (AFTER renaming to standard names)
        logger.info(f"✓ Data validation complete:")
        try:
            pval_min = float(input_df['pvalue'].min())
            pval_max = float(input_df['pvalue'].max())
            logger.info(f"  - pvalue column: ✓ (range: {pval_min:.4f} to {pval_max:.4f})")
        except Exception:
            logger.info(f"  - pvalue column: ✓")
        try:
            fc_min = float(input_df['log2FC'].min())
            fc_max = float(input_df['log2FC'].max())
            logger.info(f"  - log2FC column: ✓ (range: {fc_min:.2f} to {fc_max:.2f})")
        except Exception:
            logger.info(f"  - log2FC column: ✓")
        
        return input_df

    def run_pathway_mapping_with_dataframe(self, input_df: pd.DataFrame, pvalue_threshold: float = None):
        """
        Run pathway mapping with a pre-loaded DataFrame
        This allows the GUI to pass filtered data directly
        
        Parameters
        ----------
        input_df : pd.DataFrame
            Input dataframe with metabolites
        pvalue_threshold : float, optional
            NOT USED for annotation. Only relevant for Network Analysis tab.
            Annotation does NOT require or use p-values.
        """
        try:
            # Log that we're using pre-loaded data
            logger.info(f"Using pre-loaded DataFrame with {len(input_df)} metabolites")
            
            # Check if data already has pathway information
            pathway_cols = ['HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways', 'Metabolika_Pathways']
            has_pathway_data = any(col in input_df.columns for col in pathway_cols)
            
            if has_pathway_data:
                logger.info("Data already contains pathway information - creating All_Pathways list column")
                mapped_df = self._create_all_pathways_column(input_df.copy())
            else:
                # ANNOTATION MODE: Skip validation - just need Name column
                # Do NOT validate p-value or log2FC - not required for annotation
                if 'Name' not in input_df.columns:
                    raise ValueError(
                        "Input must contain a 'Name' column. "
                        "Use 'Verify Columns' to assign Feature ID."
                    )
                
                logger.info("✓ Name column found - ready for pathway annotation")
                logger.info(f"✓ Annotating all {len(input_df)} metabolites (no filtering)")
                
                # Log data validation results - Name column only
                logger.info("  - Name column: ✓")
                
                logger.info(f"Processing {len(input_df)} metabolites for pathway mapping")

                # Process metabolites
                mapped_df = self.process_metabolites(input_df)

            # *** CRITICAL: Reclassify disease pathways BEFORE creating summary and saving ***
            # This ensures disease pathways are moved to Associated_Diseases at the source
            mapped_df = self.reclassify_disease_pathways(mapped_df)
            
            # Create pathway summary (after reclassification to ensure clean data)
            # Note: Summary will not be saved as a sheet; kept in-memory only
            pathway_summary_df = self.create_pathway_summary(mapped_df)

            # Save results
            logger.info(f"Saving pathway-mapped results to: {self.output_file}")
            with pd.ExcelWriter(self.output_file, engine='openpyxl') as writer:
                # Drop the All_Pathways column (contains list objects) before saving to Excel
                # Keep All_Pathways_Display which contains the same info in readable string format
                columns_to_drop = ['All_Pathways'] if 'All_Pathways' in mapped_df.columns else []
                excel_df = mapped_df.drop(columns=columns_to_drop)
                excel_df.to_excel(writer, sheet_name='Metabolites_with_Pathways', index=False)
            # Do not save 'Pathway_Summary' sheet or separate summary file (per user request)

            # Print summary statistics
            total_metabolites = len(mapped_df)
            metabolites_with_pathways = mapped_df['total_pathways'].gt(0).sum()
            avg_pathways = mapped_df[mapped_df['total_pathways'] > 0]['total_pathways'].mean()
            
            hmdb_coverage = mapped_df['HMDB_Pathways'].str.strip().ne('').sum()
            pathbank_coverage = mapped_df['PathBank_Pathways'].str.strip().ne('').sum()
            SMPDB_coverage = mapped_df['SMPDB_Pathways'].str.strip().ne('').sum()
            wikipathways_coverage = mapped_df['WikiPathways'].str.strip().ne('').sum()

            logger.info("="*60)
            logger.info("PATHWAY MAPPING SUMMARY")
            logger.info("="*60)
            logger.info(f"Total metabolites processed: {total_metabolites}")
            logger.info(f"Metabolites with pathways: {metabolites_with_pathways} ({metabolites_with_pathways/total_metabolites*100:.1f}%)")
            logger.info(f"Average pathways per metabolite: {avg_pathways:.1f}")
            logger.info(f"Total unique pathways (not written to Excel): {len(pathway_summary_df)}")
            logger.info("")
            logger.info("DATABASE COVERAGE:")
            logger.info(f"HMDB pathways: {hmdb_coverage} ({hmdb_coverage/total_metabolites*100:.1f}%)")
            logger.info(f"PathBank pathways: {pathbank_coverage} ({pathbank_coverage/total_metabolites*100:.1f}%)")
            logger.info(f"SMPDB pathways: {SMPDB_coverage} ({SMPDB_coverage/total_metabolites*100:.1f}%)")
            logger.info(f"WikiPathways: {wikipathways_coverage} ({wikipathways_coverage/total_metabolites*100:.1f}%)")
            logger.info("="*60)

            return mapped_df, pathway_summary_df

        except Exception as e:
            logger.error(f"Error in pathway mapping process: {e}")
            raise
    
    def reclassify_disease_pathways(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Reclassify disease pathways from pathway columns to Associated_Diseases column.
        This ensures disease pathways are correctly categorized at the data source level.
        
        Args:
            df: DataFrame with pathway annotations
            
        Returns:
            DataFrame with disease pathways moved to Associated_Diseases
        """
        logger.info("🔄 Reclassifying disease pathways to Associated_Diseases column...")

        # Debug trace for a commonly reported pathway
        def _count_occurrences(frame: pd.DataFrame, term: str) -> Tuple[int, int]:
            term_l = term.lower()
            in_pathways = 0
            in_diseases = 0
            pathway_cols = [
                'HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways',
                'Metabolika_Pathways', 'Reactome_Pathways', 'KEGG_Pathways', 'All_Pathways_Display'
            ]
            for col in pathway_cols:
                if col in frame.columns:
                    series = frame[col].fillna('').astype(str)
                    in_pathways += series.str.lower().str.contains(term_l).sum()
            if 'Associated_Diseases' in frame.columns:
                series = frame['Associated_Diseases'].fillna('').astype(str)
                in_diseases = series.str.lower().str.contains(term_l).sum()
            return in_pathways, in_diseases

        try:
            wb_before_p, wb_before_d = _count_occurrences(df, 'Warburg Effect')
            logger.info(f"[TRACE] 'Warburg Effect' before reclassify → pathways: {wb_before_p}, diseases: {wb_before_d}")
        except Exception:
            pass
        
        # Pathway columns to check (in order of priority)
        pathway_columns = [
            'HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways',
            'Metabolika_Pathways', 'Reactome_Pathways', 'KEGG_Pathways'
        ]
        
        total_reclassified = 0
        metabolites_affected = 0
        
        for idx, row in df.iterrows():
            disease_pathways_found = []
            modified_pathway_columns = {}
            
            # Check each pathway column
            for col in pathway_columns:
                if col not in df.columns:
                    continue
                    
                pathway_str = row.get(col, '')
                if not pathway_str or pd.isna(pathway_str) or str(pathway_str).strip() == '':
                    continue
                
                # Determine separator based on column
                if col in ['HMDB_Pathways', 'Metabolika_Pathways']:
                    separator = '|'
                else:
                    separator = ';'
                
                # Split pathways
                pathways = [p.strip() for p in str(pathway_str).split(separator) if p.strip()]
                
                # Separate disease pathways from regular pathways
                valid_pathways = []
                disease_pathways = []
                
                for pathway in pathways:
                    if is_disease_pathway(pathway):
                        # Normalize disease pathway name using hardcoded corrections
                        normalized_name = normalize_disease_pathway_name(pathway)
                        disease_pathways.append(normalized_name)
                        disease_pathways_found.append(normalized_name)
                        total_reclassified += 1
                    else:
                        valid_pathways.append(pathway)
                
                # Update pathway column (remove disease pathways)
                if disease_pathways:
                    modified_pathway_columns[col] = separator.join(valid_pathways) if valid_pathways else ''
            
            # If we found disease pathways, update the row
            if disease_pathways_found:
                metabolites_affected += 1
                
                # Update pathway columns (remove disease pathways)
                for col, new_value in modified_pathway_columns.items():
                    old_value = df.at[idx, col]
                    df.at[idx, col] = new_value
                    # Log the first few changes for debugging
                    if total_reclassified <= 10:
                        logger.info(f"  [DISEASE REMOVAL] {row.get('Name', 'Unknown')} | {col}:")
                        logger.info(f"    BEFORE: {old_value}")
                        logger.info(f"    AFTER:  {new_value}")
                        logger.info(f"    Diseases moved: {disease_pathways}")
                
                # Add disease pathways to Associated_Diseases
                existing_diseases = row.get('Associated_Diseases', '')
                if existing_diseases and str(existing_diseases).strip() and str(existing_diseases).lower() != 'nan':
                    # Split existing diseases
                    if '|' in str(existing_diseases):
                        existing_disease_list = [d.strip() for d in str(existing_diseases).split('|') if d.strip()]
                    else:
                        existing_disease_list = [d.strip() for d in str(existing_diseases).split(';') if d.strip()]
                else:
                    existing_disease_list = []
                
                # Combine with newly found disease pathways
                all_diseases = existing_disease_list + disease_pathways_found
                
                # Remove duplicates while preserving order
                unique_diseases = []
                seen = set()
                for disease in all_diseases:
                    disease_lower = disease.lower()
                    if disease_lower not in seen:
                        unique_diseases.append(disease)
                        seen.add(disease_lower)
                
                # Update Associated_Diseases with | separator
                df.at[idx, 'Associated_Diseases'] = ' | '.join(unique_diseases)
        
        logger.info(f"✓ Reclassification complete:")
        logger.info(f"  - {total_reclassified} disease pathways moved to Associated_Diseases")
        logger.info(f"  - {metabolites_affected} metabolites affected")
        
        # Verify individual pathway columns after reclassification
        logger.info("[VERIFY] Checking individual pathway columns for remaining diseases...")
        pathway_columns = [
            'HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways',
            'Metabolika_Pathways', 'Reactome_Pathways', 'KEGG_Pathways'
        ]
        for col in pathway_columns:
            if col in df.columns:
                # Count rows with disease keywords
                disease_count = df[col].fillna('').astype(str).str.lower().str.contains('disease|disorder|deficiency|emia|uria|osis|acidemia|aciduria').sum()
                if disease_count > 0:
                    logger.warning(f"[VERIFY] {col} still contains {disease_count} rows with potential disease pathways!")
                    # Show first few examples
                    samples = df[df[col].fillna('').astype(str).str.lower().str.contains('disease|disorder|deficiency|emia|uria|osis|acidemia|aciduria')][col].head(3).tolist()
                    for sample in samples:
                        logger.warning(f"  Example: {sample}")
        
        # Always rebuild All_Pathways and All_Pathways_Display columns after reclassification
        # This ensures diseases are removed from All_Pathways even if columns don't exist yet
        logger.info("🔄 Rebuilding All_Pathways columns after disease reclassification...")
        df = self._create_all_pathways_column(df)

        # Post-reclassification trace
        try:
            wb_after_p, wb_after_d = _count_occurrences(df, 'Warburg Effect')
            logger.info(f"[TRACE] 'Warburg Effect' after reclassify → pathways: {wb_after_p}, diseases: {wb_after_d}")
        except Exception:
            pass
        
        # Final verification: check All_Pathways for diseases after rebuild
        if 'All_Pathways' in df.columns:
            disease_in_all_pathways = 0
            for idx, row in df.iterrows():
                all_pw = row.get('All_Pathways', [])
                if isinstance(all_pw, str):
                    all_pw = [p.strip() for p in all_pw.split(';') if p.strip()]
                for pw in all_pw:
                    if is_disease_pathway(pw):
                        disease_in_all_pathways += 1
                        if disease_in_all_pathways <= 3:  # Log first 3 only
                            logger.error(f"[BUG!] Disease pathway found in All_Pathways after reclassification: '{pw}' in metabolite {row.get('Name', 'Unknown')}")
            if disease_in_all_pathways > 0:
                logger.error(f"[BUG!] Found {disease_in_all_pathways} disease pathways still in All_Pathways after reclassification!")
            else:
                logger.info("✓ No diseases found in All_Pathways after reclassification")
        
        return df
    
    def _create_all_pathways_column(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Create the All_Pathways list column from existing pathway columns
        This is needed for pathway network analysis
        """
        logger.info("Creating All_Pathways list column from existing pathway data")
        
        # Ensure we have the required columns for network analysis
        if 'total_pathways' not in df.columns:
            df['total_pathways'] = 0
        if 'pathway_sources' not in df.columns:
            df['pathway_sources'] = 'Unknown'
        
        all_pathways_list = []
        
        for idx, row in df.iterrows():
            pathways = []
            
            # Collect pathways from all sources
            pathway_sources = []
            
            # HMDB pathways (separator: |)
            hmdb_pathways = row.get('HMDB_Pathways', '')
            if hmdb_pathways and str(hmdb_pathways).strip() and str(hmdb_pathways).lower() != 'nan':
                hmdb_list = self._split_pathways(str(hmdb_pathways), separator='|')
                pathways.extend(hmdb_list)
                if hmdb_list:
                    pathway_sources.append('HMDB')
            
            # PathBank pathways (separator: ;)
            pathbank_pathways = row.get('PathBank_Pathways', '')
            if pathbank_pathways and str(pathbank_pathways).strip() and str(pathbank_pathways).lower() != 'nan':
                pathbank_list = self._split_pathways(str(pathbank_pathways), separator=';')
                pathways.extend(pathbank_list)
                if pathbank_list:
                    pathway_sources.append('PathBank')
            
            # SMPDB pathways (separator: |)
            smpdb_pathways = row.get('SMPDB_Pathways', '')
            if smpdb_pathways and str(smpdb_pathways).strip() and str(smpdb_pathways).lower() != 'nan':
                smpdb_list = self._split_pathways(str(smpdb_pathways), separator='|')
                pathways.extend(smpdb_list)
                if smpdb_list:
                    pathway_sources.append('SMPDB')
            
            # WikiPathways (separator: ;)
            wiki_pathways = row.get('WikiPathways', '')
            if wiki_pathways and str(wiki_pathways).strip() and str(wiki_pathways).lower() != 'nan':
                wiki_list = self._split_pathways(str(wiki_pathways), separator=';')
                pathways.extend(wiki_list)
                if wiki_list:
                    pathway_sources.append('WikiPathways')
            
            # Metabolika pathways (separator: ;)
            metabolika_pathways = row.get('Metabolika_Pathways', '')
            if metabolika_pathways and str(metabolika_pathways).strip() and str(metabolika_pathways).lower() != 'nan':
                metabolika_list = self._split_pathways(str(metabolika_pathways), separator=';')
                pathways.extend(metabolika_list)
                if metabolika_list:
                    pathway_sources.append('Metabolika')
            
            # Reactome pathways (separator: ;) – include in sources and All_Pathways
            reactome_pathways = row.get('Reactome_Pathways', '')
            if reactome_pathways and str(reactome_pathways).strip() and str(reactome_pathways).lower() != 'nan':
                reactome_list = self._split_pathways(str(reactome_pathways), separator=';')
                pathways.extend(reactome_list)
                if reactome_list:
                    pathway_sources.append('Reactome')

            # KEGG pathways (separator: ;) – include in sources and All_Pathways
            kegg_pathways = row.get('KEGG_Pathways', '')
            if kegg_pathways and str(kegg_pathways).strip() and str(kegg_pathways).lower() != 'nan':
                kegg_list = self._split_pathways(str(kegg_pathways), separator=';')
                pathways.extend(kegg_list)
                if kegg_list:
                    pathway_sources.append('KEGG')
            
            # Normalize pathway names and remove duplicates (case-insensitive)
            # This handles: "Nicotine and Nicotinamide Metabolims" vs "Nicotine and Nicotinamide metabolims"
            unique_pathways = []
            seen_normalized = {}  # Maps normalized lowercase name to canonical (first seen) name
            
            for pathway in pathways:
                if not pathway:
                    continue
                
                normalized = normalize_pathway_name(pathway)
                normalized_lower = normalized.lower()
                
                # Check if we've seen this pathway name (case-insensitive)
                if normalized_lower not in seen_normalized:
                    # First time seeing this pathway - use the normalized name
                    unique_pathways.append(normalized)
                    seen_normalized[normalized_lower] = normalized
            
            # Apply centralized annotation filters to remove non-specific/plant/superpathway
            # and move diseases to Associated_Diseases when rebuilding All_Pathways
            retained_paths, moved_diseases = apply_annotation_filters(
                unique_pathways,
                exclude_keywords=self.exclude_pathway_keywords,
                non_specific_patterns=NON_SPECIFIC_PATHWAY_PATTERNS
            )

            # Merge any newly moved diseases into Associated_Diseases
            if moved_diseases:
                existing = str(row.get('Associated_Diseases', '') or '').strip()
                ex_list = []
                if existing:
                    sep = '|' if '|' in existing else ';'
                    ex_list = [d.strip() for d in existing.split(sep) if d.strip()]
                # dedupe preserve order
                seen = set(d.lower() for d in ex_list)
                for d in moved_diseases:
                    dl = d.lower()
                    if dl not in seen:
                        ex_list.append(d)
                        seen.add(dl)
                df.at[idx, 'Associated_Diseases'] = ' | '.join(ex_list)
            
            # Deduplicate pathways (case-insensitive) before storing
            unique_paths = []
            seen_paths_lower = set()
            for p in retained_paths:
                p_lower = p.lower()
                if p_lower not in seen_paths_lower:
                    seen_paths_lower.add(p_lower)
                    unique_paths.append(p)

            # Store as list for programmatic use (network analysis expects this)
            all_pathways_list.append(unique_paths)
            
            # Update summary columns
            df.at[idx, 'total_pathways'] = len(unique_paths)
            df.at[idx, 'pathway_sources'] = ', '.join(pathway_sources) if pathway_sources else 'None'
            
            # Update or create All_Pathways_Display column
            df.at[idx, 'All_Pathways_Display'] = '; '.join(unique_paths) if unique_paths else ''
        
        # Add the All_Pathways column as a list
        df['All_Pathways'] = all_pathways_list
        
        logger.info(f"Created All_Pathways column for {len(df)} metabolites")
        metabolites_with_pathways = sum(1 for pathways in all_pathways_list if pathways)
        logger.info(f"Metabolites with pathways: {metabolites_with_pathways}/{len(df)} ({metabolites_with_pathways/len(df)*100:.1f}%)")
        
        return df

    def run_pathway_mapping(self):
        """Main method to run complete pathway mapping process"""
        try:
            # Load input data
            logger.info(f"Loading input file: {self.input_file}")
            if self.sheet_name:
                logger.info(f"Using sheet: {self.sheet_name}")
                input_df = pd.read_excel(self.input_file, sheet_name=self.sheet_name)
            else:
                input_df = pd.read_excel(self.input_file)
            
            # Validate and process input data for network analysis requirements
            input_df = self.validate_and_process_input_data(input_df)
            
            # Log data validation results
            logger.info("  - Name column: ✓")
            if 'pvalue' in input_df.columns:
                try:
                    min_pvalue = float(input_df['pvalue'].min())
                    max_pvalue = float(input_df['pvalue'].max())
                    logger.info(f"  - pvalue column: ✓ (range: {min_pvalue:.4f} to {max_pvalue:.4f})")
                except Exception:
                    logger.info("  - pvalue column: ✓")
            if 'log2FC' in input_df.columns:
                try:
                    min_fc = float(input_df['log2FC'].min())
                    max_fc = float(input_df['log2FC'].max())
                    logger.info(f"  - log2FC column: ✓ (range: {min_fc:.2f} to {max_fc:.2f})")
                except Exception:
                    logger.info("  - log2FC column: ✓")
            
            logger.info(f"Processing {len(input_df)} metabolites for pathway mapping")

            # Process metabolites
            mapped_df = self.process_metabolites(input_df)

            # *** CRITICAL: Reclassify disease pathways BEFORE creating summary and saving ***
            # This ensures disease pathways are moved to Associated_Diseases at the source
            mapped_df = self.reclassify_disease_pathways(mapped_df)
            
            # Create pathway summary (after reclassification to ensure clean data)
            # Note: Summary will not be saved as a sheet; kept in-memory only
            pathway_summary_df = self.create_pathway_summary(mapped_df)

            # Save results
            logger.info(f"Saving pathway-mapped results to: {self.output_file}")
            with pd.ExcelWriter(self.output_file, engine='openpyxl') as writer:
                # Drop the All_Pathways column (contains list objects) before saving to Excel
                # Keep All_Pathways_Display which contains the same info in readable string format
                columns_to_drop = ['All_Pathways'] if 'All_Pathways' in mapped_df.columns else []
                excel_df = mapped_df.drop(columns=columns_to_drop)
                excel_df.to_excel(writer, sheet_name='Metabolites_with_Pathways', index=False)
            # Do not save 'Pathway_Summary' sheet or separate summary file (per user request)

            # Print summary statistics
            total_metabolites = len(mapped_df)
            metabolites_with_pathways = mapped_df['total_pathways'].gt(0).sum()
            avg_pathways = mapped_df[mapped_df['total_pathways'] > 0]['total_pathways'].mean()
            
            hmdb_coverage = mapped_df['HMDB_Pathways'].str.strip().ne('').sum()
            pathbank_coverage = mapped_df['PathBank_Pathways'].str.strip().ne('').sum()
            SMPDB_coverage = mapped_df['SMPDB_Pathways'].str.strip().ne('').sum()
            wikipathways_coverage = mapped_df['WikiPathways'].str.strip().ne('').sum()

            logger.info("="*60)
            logger.info("PATHWAY MAPPING SUMMARY")
            logger.info("="*60)
            logger.info(f"Total metabolites processed: {total_metabolites}")
            logger.info(f"Metabolites with pathways: {metabolites_with_pathways} ({metabolites_with_pathways/total_metabolites*100:.1f}%)")
            logger.info(f"Average pathways per metabolite: {avg_pathways:.1f}")
            logger.info(f"Total unique pathways: {len(pathway_summary_df)}")
            logger.info("")
            logger.info("DATABASE COVERAGE:")
            logger.info(f"HMDB pathways: {hmdb_coverage} ({hmdb_coverage/total_metabolites*100:.1f}%)")
            logger.info(f"PathBank pathways: {pathbank_coverage} ({pathbank_coverage/total_metabolites*100:.1f}%)")
            logger.info(f"SMPDB pathways: {SMPDB_coverage} ({SMPDB_coverage/total_metabolites*100:.1f}%)")
            logger.info(f"WikiPathways: {wikipathways_coverage} ({wikipathways_coverage/total_metabolites*100:.1f}%)")
            logger.info("="*60)

            return mapped_df, pathway_summary_df

        except Exception as e:
            logger.error(f"Error in pathway mapping process: {e}")
            raise

def main():
    """Command line interface for pathway mapping"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Metabolite Pathway Mapper - Phase 2')
    parser.add_argument('--input', '-i', required=True, help='Input Excel file with annotated metabolite IDs')
    parser.add_argument('--output', '-o', default='metabolite_pathways_mapped.xlsx', help='Output Excel file')
    
    args = parser.parse_args()
    
    mapper = MetabolitePathwayMapper(args.input, args.output)
    mapper.run_pathway_mapping()

if __name__ == "__main__":
    main()
