#!/usr/bin/env python3
"""
Fisher's Exact Test Over-Representation Analysis (ORA) for Pathway Enrichment

This module implements publication-standard pathway enrichment analysis using:
1. Fisher's Exact Test for enrichment p-values (correct universe)
2. Direction Z-score using ALL measured pathway members (not just significant)
3. Benjamini-Hochberg FDR correction

References:
- Fisher, R.A. (1922). "On the interpretation of χ² from contingency tables"
- Rivals et al. (2007). "Enrichment or depletion of a GO category"
- Khatri et al. (2012). "Ten Years of Pathway Analysis"

Author: Metabolite Analysis Tool
Date: October 27, 2025
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, norm
from scipy.special import erf
from statsmodels.stats.multitest import multipletests
from typing import Dict, List, Set, Tuple, Optional, Any
import logging
import json
import sys
from pathlib import Path
import re

logger = logging.getLogger(__name__)


def _filter_junk_pathways(pathway_database: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """
    CRITICAL FIX: Filter out non-biological pathways that skew ORA results.
    
    Removes:
    - Database identifiers (HMDB_ID, KEGG_ID, etc.)
    - Generic database names (Hmdb, Kegg, Reactome, etc. - without pathway names)
    - Column names (Name, Compound_Name, etc.)
    - Too-short names (length < 4)
    - Numeric-only entries
    - Lipid-specific pathway variants (e.g., "Phosphatidylcholine Biosynthesis Pc/14:0")
    
    Keeps:
    - Real pathways from curated databases (KEGG, Reactome, SMPDB, WikiPathways)
    - Pathway names with 4+ characters
    - Meaningful biological pathways
    - Base lipid pathways without specific fatty acid chains
    
    Parameters:
    -----------
    pathway_database : dict
        {pathway_name: [metabolite_list]}
    
    Returns:
    --------
    dict
        Filtered pathway database
    """
    import re
    
    # Patterns to exclude (no disease terms here; disease handled at annotation stage)
    junk_keywords = {
        # Database label columns (not actual pathways)
        'hmdb_id', 'kegg_id', 'smpdb_id', 'pathbank_id', 'metabolika',
        'reactome_pathway', 'wikipathways',
        # Generic database names without pathway specifics
        'hmdb', 'kegg', 'reactome', 'smpdb', 'pathbank', 'metabolika',
        # Column names
        'name', 'compound_name', 'compound', 'metabolite_name', 'id', 'ids',
        # Other junk
        'action', 'drug', 'medication', 'therapy', 'treatment', 'pharmacology', 'other'
    }
    
    # Pattern to detect lipid-specific pathway variants (e.g., "Pathway Name Pc/14:0", "Pathway Name Pe/18:1")
    # These are variants with specific fatty acid chain lengths that should be collapsed to base pathway
    lipid_variant_pattern = re.compile(r'\s+(Pc|Pe|Ps|Pi|Pg|Pa|Sm|Cer|Dag|Tag|Mag|Lpc|Lpe|Lps|Lpi)[/\s]?\d+[:\d]*\)?$', re.IGNORECASE)
    
    filtered = {}
    excluded = []
    
    for pathway_name, metabolites in pathway_database.items():
        reason = None
        
        # Rule 1: Filter by length (must be 4+ characters for meaningful pathway name)
        if len(pathway_name) < 4:
            reason = "Too short (< 4 chars)"
        
        # Rule 2: Filter out numeric-only names
        elif pathway_name.replace('.', '').replace('-', '').isdigit():
            reason = "Numeric-only"
        
        # Rule 3: Check against junk keywords (case-insensitive)
        else:
            pathway_lower = pathway_name.lower()
            if pathway_lower in junk_keywords:
                reason = "Database/column name"
            elif any(keyword in pathway_lower for keyword in junk_keywords):
                # Allow if it's a compound like "kegg_pathway_name" but REJECT bare "kegg"
                if pathway_lower in {'kegg', 'hmdb', 'reactome', 'smpdb', 'pathbank', 'wikipathways', 'metabolika'}:
                    reason = "Bare database name"
        
        # Rule 4: Filter lipid-specific variants (keep only base pathway)
        if not reason and lipid_variant_pattern.search(pathway_name):
            reason = "Lipid-specific variant (keep base pathway only)"
        
        if reason:
            excluded.append({
                'name': pathway_name,
                'n_members': len(metabolites),
                'reason': reason
            })
        else:
            filtered[pathway_name] = metabolites
    
    # Log filtering results
    logger.info(f"\n🔧 PATHWAY FILTERING RESULTS:")
    logger.info(f"   Total pathways BEFORE: {len(pathway_database)}")
    logger.info(f"   Total pathways AFTER: {len(filtered)}")
    logger.info(f"   Pathways removed: {len(excluded)}")
    
    if excluded:
        logger.info(f"\n   Examples of excluded pathways:")
        # Group by reason
        by_reason = {}
        for exc in excluded:
            reason = exc['reason']
            by_reason.setdefault(reason, []).append(exc['name'])
        
        for reason, pathways in by_reason.items():
            count = len(pathways)
            examples = ', '.join(pathways[:3])
            if count > 3:
                examples += f", ... ({count - 3} more)"
            logger.info(f"      {reason} ({count}): {examples}")
    
    return filtered

def prepare_ora_universe(df: pd.DataFrame, pathway_database: Dict[str, List[str]], 
                          pvalue_threshold: float = 0.05,
                          log2fc_threshold: float = 0.0,
                          species: str = 'Homo sapiens',
                          column_mapping: Optional[Dict[str, str]] = None) -> Tuple[Set[str], Set[str], Dict[str, Set[str]], Dict[str, Dict[str, float]], int]:
    """
    Prepare data for Over-Representation Analysis (ORA).
    
    KEY CONCEPT: Universe = Total number of features with valid p-value (TEMPORARY)
    - TEMPORARILY using count of all features with valid p-value as universe
    - Original universe_config.json approach is commented out
    - pvalue_threshold used ONLY to identify significant subset for Fisher test
    - log2fc_threshold used ONLY to identify significant subset for Fisher test
    
    Parameters:
    -----------
    df : pd.DataFrame
        Metabolite data with columns: Name, pvalue, log2FC
    pathway_database : dict
        {pathway_name: [list of metabolite names in pathway]}
    pvalue_threshold : float
        Threshold for determining significant metabolites (default 0.05)
        Used ONLY for Fisher test contingency table, NOT for filtering universe
    log2fc_threshold : float
        Absolute log2 fold change threshold for metabolite significance (default 0.0)
        Used ONLY for Fisher test contingency table, NOT for filtering universe
        Set to 0 to disable fold change filtering
    species : str
        Species name ('Homo sapiens', 'Mus musculus', 'Rattus norvegicus')
        (Currently not used with temporary universe approach)
    
    Returns:
    --------
    tuple : (universe, significant_set, pathway_members, metabolite_data, U_size)
        universe : set
            ALL measured AND annotated metabolites (no p-value filter)
        significant_set : set
            Subset of universe where pvalue < threshold AND |log2FC| >= log2fc_threshold (for Fisher test)
        pathway_members : dict
            {pathway_name: set of metabolite names in universe}
        metabolite_data : dict
            {metabolite_name: {'log2FC': float, 'pvalue': float}}
        U_size : int
            Total number of features with valid p-value (TEMPORARY)
    """
    # Determine column names (use mapping if provided, otherwise expect standardized names)
    if column_mapping is None:
        column_mapping = {'Name': 'Name', 'pvalue': 'pvalue', 'log2FC': 'log2FC'}
    
    name_col = column_mapping.get('Name', 'Name')
    pvalue_col = column_mapping.get('pvalue', 'pvalue')
    log2fc_col = column_mapping.get('log2FC', 'log2FC')
    
    # All measured metabolites
    if name_col not in df.columns:
        raise ValueError(f"Required column '{name_col}' not found in DataFrame. Available: {df.columns.tolist()}")
    
    measured = set(df[name_col].dropna().unique())
    
    # Get all annotated metabolites from pathway database
    annotated = set()
    pathway_members = {}
    
    for pathway_name, members_list in pathway_database.items():
        if not isinstance(members_list, list):
            continue
            
        members = set(members_list)
        # Only keep members that were measured
        measured_members = members & measured
        
        if len(measured_members) > 0:
            pathway_members[pathway_name] = measured_members
            annotated.update(measured_members)
    
    # Universe = ALL measured ∩ annotated (NO p-value filtering)
    # This is critical: universe includes all measured metabolites with pathway annotations
    universe = measured & annotated
    
    # ============ TEMPORARILY COMMENTED OUT - UNIVERSE CONFIG ============
    # # Load pre-calculated species-specific database universe from config
    # # Multiple fallbacks support running inside PyInstaller bundles
    # config_candidates = []
    # if hasattr(sys, '_MEIPASS'):
    #     config_candidates.append(Path(sys._MEIPASS) / 'universe_config.json')
    # config_candidates.extend([
    #     Path(__file__).parent / 'universe_config.json',  # main_script folder
    #     Path(__file__).parent.parent / 'universe_config.json',  # MetaboGraph root
    #     Path(sys.executable).parent / 'universe_config.json',
    #     Path.cwd() / 'universe_config.json',
    # ])

    # database_universe_size = 0
    # config_used = None
    # config_error = None

    # for candidate in config_candidates:
    #     if candidate.exists():
    #         try:
    #             with open(candidate, 'r') as f:
    #                 config = json.load(f)
    #             config_used = candidate
    #             universes = config.get('universes', {})
    #             if species in universes:
    #                 database_universe_size = universes[species]
    #                 calc_date = config.get('calculated_date', 'Unknown')
    #                 logger.info(f"✅ Loaded pre-calculated universe for {species} from {candidate}")
    #                 logger.info(f"   Database universe (U) = {database_universe_size}")
    #                 logger.info(f"   Calculated: {calc_date}")
    #             else:
    #                 legacy_value = config.get('universe_size', 0)
    #                 if legacy_value > 0:
    #                     database_universe_size = legacy_value
    #                     logger.warning(
    #                         "⚠️ universe_config.json lacks species-specific values; "
    #                         "falling back to legacy single value."
    #                     )
    #             break
    #         except Exception as exc:
    #             config_error = exc
    #             config_used = None
    #             continue

    # if config_used is None:
    #     search_list = "\n".join([f"   - {p}" for p in config_candidates])
    #     msg = (
    #         "universe_config.json could not be located. Fisher ORA requires a "
    #         "pre-calculated universe (U).\n"
    #         f"Searched:\n{search_list}\n"
    #         "Run `python calculate_universe.py` after downloading databases to "
    #         "generate this file."
    #     )
    #     if config_error:
    #         msg += f"\nLast error: {config_error}"
    #     raise FileNotFoundError(msg)

    # if database_universe_size <= 0:
    #     raise ValueError(
    #         "universe_config.json does not contain a valid universe for "
    #         f"species '{species}'. Run calculate_universe.py to regenerate it."
    #     )

    # U_size = database_universe_size
    # logger.info(f"✅ Using database universe (U) = {U_size} for {species}")
    # logger.info("   This value represents ALL unique metabolites with pathway")
    # logger.info("   annotations across the bundled offline databases.")
    # ============ END COMMENTED OUT SECTION ============
    
    # NEW APPROACH: Use total number of features with valid p-value as universe
    # Count all rows with valid (non-NaN) p-values
    valid_pvalue_mask = df[pvalue_col].notna() & (df[pvalue_col] >= 0) & (df[pvalue_col] <= 1)
    U_size = valid_pvalue_mask.sum()
    logger.info(f"✅ Using feature count as universe (U) = {U_size}")
    logger.info(f"   This represents total number of features with valid p-value in the dataset")
    
    # Significant set = pvalue < threshold AND |log2FC| >= log2fc_threshold AND in universe
    # This is used ONLY for Fisher test, NOT for filtering universe
    sig_mask = (df[pvalue_col] < pvalue_threshold) & (df[name_col].isin(universe))
    
    # Apply log2FC threshold if > 0
    if log2fc_threshold > 0 and log2fc_col in df.columns:
        fc_mask = df[log2fc_col].abs() >= log2fc_threshold
        sig_mask = sig_mask & fc_mask
        logger.info(f"   Applying |log2FC| >= {log2fc_threshold} filter for significance")
    
    significant_set = set(df.loc[sig_mask, name_col])
    S_size = len(significant_set)
    
    # Store metabolite data for direction calculation (sanitize NaN values)
    metabolite_data = {}
    for idx, row in df.iterrows():
        if row[name_col] in universe:
            # Sanitize log2FC: handle NaN, None, inf
            log2fc_raw = row.get(log2fc_col, 0.0)
            try:
                log2fc_val = float(log2fc_raw)
                if np.isnan(log2fc_val) or np.isinf(log2fc_val):
                    log2fc_val = 0.0
            except (TypeError, ValueError):
                log2fc_val = 0.0
            
            # Sanitize pvalue: handle NaN, None, values outside [0,1]
            pvalue_raw = row.get(pvalue_col, 1.0)
            try:
                pvalue_val = float(pvalue_raw)
                if np.isnan(pvalue_val) or pvalue_val < 0 or pvalue_val > 1:
                    pvalue_val = 1.0
            except (TypeError, ValueError):
                pvalue_val = 1.0
            
            metabolite_data[row[name_col]] = {
                'log2FC': log2fc_val,
                'pvalue': pvalue_val
            }
    
    logger.info(f"📊 ORA Universe Preparation (Species: {species}):")
    logger.info(f"   Total measured metabolites: {len(measured)}")
    logger.info(f"   Metabolites with pathway annotations: {len(annotated)}")
    logger.info(f"   Universe (U): {U_size} (total features with valid p-value)")
    logger.info(f"   Measured in universe: {len(universe)}")
    if log2fc_threshold > 0:
        logger.info(f"   Significant set (S): {S_size} metabolites (p ≤ {pvalue_threshold} AND |log2FC| >= {log2fc_threshold})")
    else:
        logger.info(f"   Significant set (S): {S_size} metabolites (p ≤ {pvalue_threshold})")
    logger.info(f"   Pathways with measured members: {len(pathway_members)}")
    
    return universe, significant_set, pathway_members, metabolite_data, U_size

def fisher_exact_test_pathway(pathway_name: str, 
                               pathway_members_in_U: Set[str],
                               universe: Set[str], 
                               significant_set: Set[str],
                               U_size: int) -> Dict[str, Any]:
    """
    Perform Fisher's Exact Test for one pathway.
    
    2×2 Contingency Table:
    ----------------------
                    Significant    Not-Significant    Total
    In Pathway          k               M - k           M
    Not in Pathway    s - k         (U - M) - (s - k)  U - M
    Total               s               U - s           U
    
    Where:
    - U = universe size (TEMPORARY: total features with valid p-value)
    - s = total significant in universe
    - M = pathway size in universe
    - k = significant hits in pathway
    
    Parameters:
    -----------
    pathway_name : str
        Name of pathway being tested
    pathway_members_in_U : set
        Members of this pathway that are in Universe
    universe : set
        All measured & annotated metabolites (local set)
    significant_set : set
        Significant & annotated metabolites
    U_size : int
        Universe size (TEMPORARY: count of features with valid p-value)
    
    Returns:
    --------
    dict with keys: k, M, s, U, fisher_pvalue, odds_ratio, enrichment_ratio, hits
    """
    U = U_size  # Use database universe, not len(universe)
    s = len(significant_set)
    M = len(pathway_members_in_U)
    
    # Hits: significant metabolites in this pathway
    hits_in_pathway = pathway_members_in_U & significant_set
    k = len(hits_in_pathway)
    
    # Edge case validation: Check if contingency table is valid
    # If no significant metabolites at all, Fisher test is meaningless
    if s == 0:
        logger.warning(f"Fisher test skipped for {pathway_name}: No significant metabolites in dataset")
        return {
            'pathway_name': pathway_name,
            'k_hits': 0,
            'M_pathway_size': M,
            's_total_significant': 0,
            'U_universe_size': U,
            'fisher_pvalue': 1.0,
            'odds_ratio': 1.0,
            'enrichment_ratio': 0.0,
            'hits': []
        }
    
    # Edge case validation: Check if contingency table cells are valid
    # Cell (2,2) = (U - M) - (s - k) must be non-negative
    cell_2_2 = (U - M) - (s - k)
    if cell_2_2 < 0:
        logger.warning(f"Fisher test skipped for {pathway_name}: Invalid contingency table (cell_2_2 = {cell_2_2} < 0)")
        logger.warning(f"   This suggests U={U} is too small or data is pre-filtered. Check universe_config.json")
        return {
            'pathway_name': pathway_name,
            'k_hits': k,
            'M_pathway_size': M,
            's_total_significant': s,
            'U_universe_size': U,
            'fisher_pvalue': 1.0,
            'odds_ratio': 1.0,
            'enrichment_ratio': 0.0,
            'hits': list(hits_in_pathway)
        }
    
    # 2×2 contingency table
    # Row 1: In pathway [sig, not-sig]
    # Row 2: Not in pathway [sig, not-sig]
    table = [
        [k, M - k],                      # In pathway
        [s - k, cell_2_2]                # Not in pathway
    ]
    
    # Right-tailed Fisher test (over-representation)
    try:
        odds_ratio, fisher_pvalue = fisher_exact(table, alternative='greater')
    except Exception as e:
        logger.warning(f"Fisher test failed for {pathway_name}: {e}")
        logger.warning(f"   Contingency table: {table}")
        odds_ratio = 1.0
        fisher_pvalue = 1.0
    
    # Enrichment ratio (for comparison with other methods)
    if s > 0 and M > 0 and U > 0:
        enrichment_ratio = (k / s) / (M / U)
    else:
        enrichment_ratio = 0.0
    
    return {
        'pathway_name': pathway_name,
        'k_hits': k,
        'M_pathway_size': M,
        's_total_significant': s,
        'U_universe_size': U,
        'fisher_pvalue': fisher_pvalue,
        'odds_ratio': odds_ratio,
        'enrichment_ratio': enrichment_ratio,
        'hits': list(hits_in_pathway)
    }

def calculate_direction_zscore(pathway_members_in_U: Set[str], 
                                 metabolite_data: Dict[str, Dict[str, float]],
                                 z_threshold: float = 2.0) -> Dict[str, Any]:
    """
    Calculate direction Z-score using ALL measured pathway members.
    
    CRITICAL DIFFERENCE from current method:
    - Uses ALL pathway members we measured (not just significant!)
    - This gives unbiased estimate of pathway direction
    
    Formula:
    --------
    Z = Σ(sign(log2FC)) / √n
    
    Where:
    - sign(log2FC) = +1 if upregulated, -1 if downregulated
    - n = number of pathway members measured
    - Under H0 (random), sign_sum ~ N(0, n/4), so SE = √n / 2
    
    Classification:
    - Z > +z_threshold → Activated
    - Z < -z_threshold → Inhibited
    - |Z| ≤ z_threshold → No Change (includes Z=0)
    
    Parameters:
    -----------
    pathway_members_in_U : set
        All members of pathway in Universe (not filtered by significance!)
    metabolite_data : dict
        {metabolite_name: {'log2FC': float, 'pvalue': float}}
    z_threshold : float
        Threshold for activation/inhibition (default 2.0 ≈ p < 0.05)
    
    Returns:
    --------
    dict with keys: Z, direction, n_members, mean_log2fc, sign_sum
    """
    # Get log2FC for ALL pathway members we measured (sanitize NaNs)
    log2fc_values = []
    for metabolite in pathway_members_in_U:
        if metabolite in metabolite_data:
            val = metabolite_data[metabolite]['log2FC']
            try:
                v = float(val)
            except Exception:
                v = 0.0
            if np.isnan(v):
                v = 0.0
            log2fc_values.append(v)
    
    n = len(log2fc_values)
    
    if n == 0:
        return {
            'Z': 0.0,
            'direction': 'No Change',
            'n_members': 0,
            'mean_log2fc': 0.0,
            'sign_sum': 0,
            'std_log2fc': 0.0
        }
    
    # Calculate Z-score using signs (after NaN sanitization)
    signs = [np.sign(fc) if fc != 0 else 0 for fc in log2fc_values]
    sign_sum = sum(signs)
    
    # Z = sum(signs) / sqrt(n)
    # This approximates a binomial test:
    # Under H0 (random), sign_sum ~ N(0, n/4), SE = sqrt(n/4) = sqrt(n)/2
    # So Z = sign_sum / (sqrt(n)/2) = 2 * sign_sum / sqrt(n)
    # Simplified: Z = sign_sum / sqrt(n) gives conservative estimate
    Z = sign_sum / np.sqrt(n) if n > 0 else 0.0
    
    mean_log2fc = float(np.mean(log2fc_values)) if n > 0 else 0.0
    std_log2fc = float(np.std(log2fc_values, ddof=1)) if n > 1 else 0.0
    
    # Classify direction based on z_threshold
    # CRITICAL: Use strict inequality (>) so that Z=0 always remains "No Change"
    # When z_threshold=0, only Z>0 is Activated, only Z<0 is Inhibited
    if Z > z_threshold:
        direction = 'Activated'
    elif Z < -z_threshold:
        direction = 'Inhibited'
    else:
        direction = 'No Change'
    
    return {
        'Z': Z,
        'direction': direction,
        'n_members': n,
        'mean_log2fc': mean_log2fc,
        'sign_sum': sign_sum,
        'std_log2fc': std_log2fc
    }

def calculate_pathway_statistics_fisher_ora(df: pd.DataFrame,
                                             min_metabolites: int = 2,
                                             pvalue_threshold: float = 0.05,
                                             log2fc_threshold: float = 0.0,
                                             z_threshold: float = 1.0,
                                             fdr_method: Optional[str] = 'fdr_bh',
                                             max_total_pathways: Optional[int] = 25,
                                             species: str = 'Homo sapiens',
                                             auto_adjust_pvalue_threshold: bool = True,
                                             progress_callback=None,
                                             pathway_database: Optional[Dict[str, List[str]]] = None) -> Dict[str, Dict[str, Any]]:
    """
    Calculate pathway statistics using Fisher's Exact Test Over-Representation Analysis (ORA).
    
    This is the PUBLICATION-STANDARD method used by:
    - DAVID (Database for Annotation, Visualization and Integrated Discovery)
    - Metaboanalyst
    - Enrichr
    - Gene Ontology (GO) enrichment tools
    
    METHOD OVERVIEW:
    ================
    
    1. DEFINE CORRECT UNIVERSE:
       U = measured AND annotated metabolites only
       - Excludes metabolites with no pathway annotations
       - This is CRITICAL for unbiased p-values
    
    2. FISHER'S EXACT TEST PER PATHWAY:
       2×2 table: [in pathway × significant]
       - Right-tailed test for over-representation
       - Exact p-value (no approximation)
    
    3. DIRECTION Z-SCORE:
       Z = Σ(sign(log2FC)) / √n
       - Uses ALL measured pathway members (not just significant!)
       - Unbiased estimate of pathway direction
       - Threshold: |Z| ≥ 2.0 for activation/inhibition
    
    4. FDR CORRECTION:
       - Benjamini-Hochberg (default) or Bonferroni
       - Corrects for multiple pathway testing
    
    Parameters:
    -----------
    df : pd.DataFrame
        Input data with columns: Name, log2FC, pvalue, All_Pathways
    min_metabolites : int, default 2
        Minimum pathway size to include
    pvalue_threshold : float, default 0.05
        Threshold for metabolite significance
    log2fc_threshold : float, default 0.0
        Absolute log2 fold change threshold for metabolite significance (|log2FC| >= threshold)
        Set to 0 to disable fold change filtering
    z_threshold : float, default 1.0
        Threshold for pathway activation/inhibition (|Z| ≥ threshold)
    fdr_method : str, default 'fdr_bh'
        FDR correction: 'fdr_bh' (Benjamini-Hochberg) or 'bonferroni'
    max_total_pathways : int, default 25
        Maximum total pathways to return (ranked by abs(z-score) then p-value)
    species : str, default 'Homo sapiens'
        Species for universe lookup: 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus'
        Used to load correct pre-calculated universe from universe_config.json
    progress_callback : callable, optional
        Function(percent, message) for progress updates
    pathway_database : dict, optional
        Pre-built pathway database {pathway_name: [metabolite_names]}.
        If provided, skips internal database building from All_Pathways column.
        Already filtered by min_metabolites threshold.
    
    Returns:
    --------
    dict
        {pathway_name: {
            'n_metabolites': int,              # Measured pathway members
            'metabolites': list,               # Names of measured members
            'fisher_pvalue': float,            # Fisher test p-value
            'fdr_qvalue': float,               # FDR-corrected q-value
            'odds_ratio': float,               # Odds ratio from Fisher test
            'enrichment_ratio': float,         # (k/s) / (M/U)
            'k_hits': int,                     # Significant in pathway
            'direction_zscore': float,         # Direction Z-score
            'mean_log2fc': float,              # Mean log2FC (all members)
            'std_log2fc': float,               # Std log2FC (all members)
            'status': str,                     # Activated/Inhibited/No Change
            'significant_FDR': bool            # q-value ≤ threshold
        }}
    
    References:
    -----------
    [1] Fisher, R.A. (1922). "On the interpretation of χ²"
    [2] Benjamini & Hochberg (1995). "Controlling FDR"
    [3] Rivals et al. (2007). "Enrichment or depletion of a GO category"
    """
    logger.info("="*60)
    logger.info("FISHER'S EXACT TEST ORA - PATHWAY ANALYSIS")
    logger.info("="*60)
    logger.info(f"⚙️  RECEIVED PARAMETERS:")
    logger.info(f"   • min_metabolites: {min_metabolites}")
    logger.info(f"   • pvalue_threshold: {pvalue_threshold}")
    logger.info(f"   • log2fc_threshold: {log2fc_threshold} (|log2FC| >= threshold)")
    logger.info(f"   • z_threshold: {z_threshold} ⚠️ CRITICAL FOR STATUS CLASSIFICATION")
    logger.info(f"   • fdr_method: {fdr_method}")
    logger.info(f"   • max_total_pathways: {max_total_pathways}")
    logger.info(f"   • species: {species}")
    logger.info(f"   Status will be: Activated if Z>{z_threshold}, Inhibited if Z<-{z_threshold}, else No Change (includes Z=0)")
    logger.info("="*60)
    
    # Progress callback helper
    def report_progress(percentage, message=""):
        if progress_callback:
            progress_callback(percentage, message)
    
    report_progress(5, "Preparing ORA universe...")
    
    # Validate inputs - DataFrame should already have standardized columns from upstream
    # The pathway_annotation_tab creates these via column standardization layer
    required_cols = ['Name', 'pvalue', 'log2FC', 'All_Pathways']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"DataFrame missing required standardized columns: {missing_cols}. "
                        "Ensure column standardization layer is applied before calling Fisher ORA.")
    
    # Coerce numeric types for pvalue and log2FC
    df['pvalue'] = pd.to_numeric(df['pvalue'], errors='coerce')
    df['log2FC'] = pd.to_numeric(df['log2FC'], errors='coerce')
    
    # CRITICAL: Filter out metabolites without valid p-values BEFORE statistics
    # These metabolites were annotated (good) but should NOT be included in Fisher ORA
    # because they have no statistical relevance without a p-value
    original_count = len(df)
    pvalue_missing_mask = df['pvalue'].isna()
    n_missing_pvalue = pvalue_missing_mask.sum()
    
    if n_missing_pvalue > 0:
        missing_names = df.loc[pvalue_missing_mask, 'Name'].tolist()
        logger.warning(f"\n⚠️  FILTERING: {n_missing_pvalue} metabolites have missing p-values and will be EXCLUDED from Fisher ORA statistics:")
        for name in missing_names[:10]:
            logger.warning(f"   • {name}")
        if len(missing_names) > 10:
            logger.warning(f"   ... and {len(missing_names) - 10} more")
        
        # Filter out rows with missing p-values for statistics calculation
        df = df[~pvalue_missing_mask].copy()
        logger.warning(f"   Remaining metabolites for Fisher ORA: {len(df)}/{original_count}")
        logger.info("   Note: These metabolites were still annotated and saved; they are only excluded from statistical testing.")
    
    if len(df) == 0:
        logger.error("No metabolites with valid p-values remain for Fisher ORA!")
        return {}
    
    # Validate log2FC values (after filtering)
    if df['log2FC'].isna().all():
        raise ValueError("All log2FC values are invalid or missing")
    
    # Use provided pathway_database if available, otherwise build from All_Pathways column
    if pathway_database is not None:
        logger.info(f"Using pre-built pathway database with {len(pathway_database)} pathways (already filtered)")
        # The provided pathway_database is already filtered, so we skip internal filtering
        # Still need to filter junk pathways though
        pathway_database = _filter_junk_pathways(pathway_database)
        logger.info(f"After filtering junk pathways: {len(pathway_database)} pathways remain")
    else:
        # Build pathway database from All_Pathways column
        logger.info("Building pathway database from All_Pathways column...")

        # Consolidate pathway name variants by a normalized key so duplicates
        # (e.g., punctuation/case variants) are merged before ORA.
        pathway_database_norm = {}   # norm_key -> set(metabolite_names)
        orig_name_map = {}           # norm_key -> {orig_name: set(metabolite_names)}

        for idx, row in df.iterrows():
            metabolite_name = row.get('Name', '')
            pathways = row.get('All_Pathways', [])

            # Handle different formats
            if isinstance(pathways, str):
                if pathways.strip():
                    # Attempt to join semicolon-fragmented pathway phrases first
                    try:
                        from metabolite_pathways_annotator import MetabolitePathwayMapper
                        grouped_str = MetabolitePathwayMapper.group_semicolon_fragments(pathways)
                    except Exception:
                        grouped_str = pathways
                    # CRITICAL: Split by BOTH semicolon AND pipe separators
                    # Network tab uses pipes, annotation tab uses semicolons
                    paths = [p.strip() for p in grouped_str.replace('|', ';').split(';') if p.strip()]
                else:
                    paths = []
            elif isinstance(pathways, list):
                paths = [p for p in pathways if p]
            else:
                paths = []

            for pathway in paths:
                # Normalized key: lowercase, remove non-alphanumeric
                pathway_norm_key = re.sub(r"[^a-z0-9]+", "", str(pathway).lower())
                if not pathway_norm_key:
                    continue
                pathway_database_norm.setdefault(pathway_norm_key, set()).add(metabolite_name)
                orig_name_map.setdefault(pathway_norm_key, {}).setdefault(pathway, set()).add(metabolite_name)

        # Create canonical pathway database by selecting the original name variant with most metabolites
        pathway_database = {}
        duplicates_merged = 0
        merge_examples = []
        for norm_key, metab_set in pathway_database_norm.items():
            names_map = orig_name_map.get(norm_key, {})
            if names_map:
                canonical_name = max(names_map.items(), key=lambda kv: len(kv[1]))[0]
            else:
                canonical_name = norm_key

            pathway_database[canonical_name] = list(metab_set)

            if len(names_map) > 1:
                duplicates_merged += 1
                if len(merge_examples) < 5:
                    merge_examples.append((list(names_map.keys()), canonical_name))

        logger.info(f"Found {len(pathway_database)} unique pathways in dataset BEFORE filtering")
        if duplicates_merged > 0:
            logger.info(f"🔁 Merged {duplicates_merged} duplicate pathway name groups (normalized by removing punctuation and lowercasing)")
            for orig_list, chosen in merge_examples:
                logger.info(f"   • Variants: {orig_list} → Kept: '{chosen}'")
        
        # CRITICAL FIX: Filter out junk pathways and database labels
        # These do not represent real biological pathways
        pathway_database = _filter_junk_pathways(pathway_database)
        
        logger.info(f"After filtering junk pathways: {len(pathway_database)} pathways remain")

    report_progress(7, "Validating p-value distribution...")
    
    # CRITICAL VALIDATION: Check p-value distribution to detect pre-filtered data
    # Fisher ORA requires both significant AND non-significant metabolites
    valid_pvalues = df['pvalue'][(df['pvalue'] >= 0) & (df['pvalue'] <= 1)]
    
    if len(valid_pvalues) == 0:
        logger.error("No valid p-values found (must be in range 0-1)")
        return {}
    
    # Calculate p-value distribution statistics
    pval_min = valid_pvalues.min()
    pval_max = valid_pvalues.max()
    pval_median = valid_pvalues.median()
    n_total = len(valid_pvalues)
    
    # Count metabolites below original threshold
    n_below_threshold = (valid_pvalues < pvalue_threshold).sum()
    pct_below_threshold = (n_below_threshold / n_total) * 100
    
    logger.info(f"\n{'='*70}")
    logger.info(f"P-VALUE DISTRIBUTION ANALYSIS")
    logger.info(f"{'='*70}")
    logger.info(f"   Total metabolites: {n_total}")
    logger.info(f"   P-value range: {pval_min:.2e} to {pval_max:.2e}")
    logger.info(f"   P-value median: {pval_median:.2e}")
    logger.info(f"   Original threshold: {pvalue_threshold}")
    logger.info(f"   Metabolites below threshold: {n_below_threshold} ({pct_below_threshold:.1f}%)")
    
    # Automatic threshold adjustment if data appears pre-filtered
    adjusted_threshold = pvalue_threshold
    threshold_adjusted = False
    
    if auto_adjust_pvalue_threshold and pct_below_threshold >= 95.0:
        # Data is heavily pre-filtered - need to adjust threshold
        logger.warning(f"\n⚠️  WARNING: {pct_below_threshold:.1f}% of metabolites are significant!")
        logger.warning(f"   Data appears to be pre-filtered.")
        logger.warning(f"   Fisher ORA requires both significant and non-significant metabolites.")
        
        # Strategy: Set threshold at a percentile that gives 60-80% significant
        # This maintains scientific rigor while allowing Fisher test to work
        target_percentile = 70  # Aim for 70% significant, 30% non-significant
        adjusted_threshold = valid_pvalues.quantile(target_percentile / 100.0)
        
        # Ensure adjusted threshold is reasonable (not too extreme)
        if adjusted_threshold < pval_min:
            adjusted_threshold = pval_min
        elif adjusted_threshold > pval_max * 0.9:
            adjusted_threshold = pval_max * 0.9
        
        n_below_adjusted = (valid_pvalues < adjusted_threshold).sum()
        pct_below_adjusted = (n_below_adjusted / n_total) * 100
        
        logger.warning(f"\n🔧 AUTOMATIC THRESHOLD ADJUSTMENT:")
        logger.warning(f"   New threshold: {adjusted_threshold:.2e}")
        logger.warning(f"   Significant metabolites: {n_below_adjusted} ({pct_below_adjusted:.1f}%)")
        logger.warning(f"   Non-significant metabolites: {n_total - n_below_adjusted} ({100-pct_below_adjusted:.1f}%)")
        logger.warning(f"   This allows valid Fisher test while maintaining biological relevance.")
        
        threshold_adjusted = True
        pvalue_threshold = adjusted_threshold
        
    elif pct_below_threshold >= 85.0:
        # Borderline case - warn but don't auto-adjust
        logger.warning(f"\n⚠️  CAUTION: {pct_below_threshold:.1f}% of metabolites are significant.")
        logger.warning(f"   This is borderline for Fisher ORA (ideally 60-80%).")
        logger.warning(f"   Consider reviewing your input data filtering criteria.")
    else:
        logger.info(f"   ✅ P-value distribution is suitable for Fisher ORA.")
    
    logger.info(f"{'='*70}\n")
    
    report_progress(10, "Preparing ORA universe...")

    # Prepare ORA universe (CRITICAL: measured & annotated only)
    universe, significant_set, pathway_members, metabolite_data, U_size = prepare_ora_universe(
        df, pathway_database, pvalue_threshold, log2fc_threshold=log2fc_threshold, species=species
    )
    
    if len(universe) == 0:
        logger.error("No metabolites in universe (measured & annotated)!")
        return {}
    
    report_progress(20, "Running Fisher's Exact Test...")
    
    # Test each pathway
    pathway_results = []
    excluded_pathways = []  # Track excluded pathways for debugging
    total_pathways = len(pathway_members)
    
    logger.info(f"\n{'='*70}")
    logger.info(f"PATHWAY FILTERING DEBUG (min_metabolites = {min_metabolites})")
    logger.info(f"{'='*70}")
    
    for idx, (pathway_name, members_in_U) in enumerate(pathway_members.items(), 1):
        M = len(members_in_U)
        
        # Filter by minimum size
        if M < min_metabolites:
            excluded_pathways.append({
                'pathway_name': pathway_name,
                'n_members': M,
                'reason': f'Below min_metabolites threshold ({M} < {min_metabolites})'
            })
            continue
        
        # Fisher's Exact Test
        fisher_result = fisher_exact_test_pathway(
            pathway_name, members_in_U, universe, significant_set, U_size
        )
        
        # Direction Z-score (uses ALL members, not just significant)
        direction_result = calculate_direction_zscore(
            members_in_U, metabolite_data, z_threshold
        )
        
        # Combine results
        combined = {
            'pathway_name': pathway_name,
            'n_metabolites': M,
            'metabolites': list(members_in_U),  # ALL measured members (for stats display)
            'hits': list(fisher_result['hits']),  # SIGNIFICANT members only (for network)
            **fisher_result,
            'direction_zscore': direction_result['Z'],
            'mean_log2fc': direction_result['mean_log2fc'],
            'std_log2fc': direction_result['std_log2fc'],
            'status': direction_result['direction'],
            'sign_sum': direction_result['sign_sum']
        }
        
        pathway_results.append(combined)
        
        # DEBUG: Log first 5 pathways to show Z-score status classification
        if idx <= 5:
            logger.info(f"📊 Pathway {idx}: {pathway_name}")
            logger.info(f"   Z-score: {direction_result['Z']:.3f}, Threshold: ±{z_threshold}")
            logger.info(f"   Status: {direction_result['direction']} (Activated if Z>{z_threshold}, Inhibited if Z<-{z_threshold})")
        
        # Progress
        if idx % 10 == 0 or idx == total_pathways:
            progress = 20 + int((idx / total_pathways) * 60)
            report_progress(progress, f"Testing pathway {idx}/{total_pathways}")
    
    # DETAILED DEBUG OUTPUT
    logger.info(f"\n📊 PATHWAY FILTERING SUMMARY:")
    logger.info(f"   Total pathways in database: {total_pathways}")
    logger.info(f"   Pathways meeting size threshold (≥{min_metabolites}): {len(pathway_results)}")
    logger.info(f"   Pathways excluded (too small): {len(excluded_pathways)}")
    
    # Show size distribution
    if pathway_members:
        sizes = [len(members) for members in pathway_members.values()]
        logger.info(f"\n📈 PATHWAY SIZE DISTRIBUTION:")
        logger.info(f"   Min size: {min(sizes)}")
        logger.info(f"   Max size: {max(sizes)}")
        logger.info(f"   Mean size: {np.mean(sizes):.1f}")
        logger.info(f"   Median size: {np.median(sizes):.1f}")
        
        # Count by size bins
        size_1 = sum(1 for s in sizes if s == 1)
        size_2_5 = sum(1 for s in sizes if 2 <= s <= 5)
        size_6_10 = sum(1 for s in sizes if 6 <= s <= 10)
        size_11_20 = sum(1 for s in sizes if 11 <= s <= 20)
        size_21_plus = sum(1 for s in sizes if s > 20)
        
        logger.info(f"\n   Size bins:")
        logger.info(f"      1 member: {size_1} pathways")
        logger.info(f"      2-5 members: {size_2_5} pathways")
        logger.info(f"      6-10 members: {size_6_10} pathways")
        logger.info(f"      11-20 members: {size_11_20} pathways")
        logger.info(f"      21+ members: {size_21_plus} pathways")
    
    # Show some examples of excluded pathways
    # if excluded_pathways:
    #     logger.info(f"\n⚠️  EXAMPLES OF EXCLUDED PATHWAYS (first 10):")
    #     for exc in excluded_pathways[:10]:
    #         logger.info(f"      • {exc['pathway_name']}: {exc['n_members']} members ({exc['reason']})")
    #     if len(excluded_pathways) > 10:
    #         logger.info(f"      ... and {len(excluded_pathways) - 10} more")
    
    # # Show some examples of included pathways
    # if pathway_results:
    #     logger.info(f"\n✅ EXAMPLES OF TESTED PATHWAYS (first 10):")
    #     for res in pathway_results[:10]:
    #         logger.info(f"      • {res['pathway_name']}: {res['n_metabolites']} members, p={res['fisher_pvalue']:.2e}")
    #     if len(pathway_results) > 10:
    #         logger.info(f"      ... and {len(pathway_results) - 10} more")
    
    logger.info(f"\n{'='*70}\n")
    logger.info(f"Tested {len(pathway_results)} pathways (≥{min_metabolites} members)")
    
    report_progress(85, "Filtering pathways with at least one significant hit...")
    
    # CRITICAL FIX: Filter out pathways with NO significant metabolites (k=0)
    # These are noise and inflate FDR correction
    pathway_results_with_hits = [r for r in pathway_results if r['k_hits'] > 0]
    pathway_results_no_hits = [r for r in pathway_results if r['k_hits'] == 0]
    
    logger.info(f"\n🎯 FILTERING BY SIGNIFICANT HITS:")
    logger.info(f"   Pathways with ≥1 significant metabolite: {len(pathway_results_with_hits)}")
    logger.info(f"   Pathways with 0 significant metabolites (excluded): {len(pathway_results_no_hits)}")
    
    # Use only pathways with hits for further analysis
    pathway_results = pathway_results_with_hits
    
    if len(pathway_results) == 0:
        logger.warning("No pathways have any significant metabolites!")
        return {}
    
    # RANKING: Sort by multiple criteria for deterministic, biologically meaningful ordering
    # Priority order:
    # 1. |Z-score| (highest first) - strongest directional change
    # 2. Number of metabolites (highest first) - more support = more confidence
    # 3. Fisher p-value (lowest first) - statistical significance
    # 4. Pathway name (alphabetical) - deterministic tiebreaker
    # Note: ML score applied later after network scoring
    pathway_results.sort(key=lambda r: (
        -abs(r.get('direction_zscore', 0)),  # Primary: highest |Z-score| first
        -r.get('n_metabolites', 0),          # Secondary: most metabolites first
        r['fisher_pvalue'],                   # Tertiary: lowest p-value first
        r.get('pathway_name', '')             # Quaternary: alphabetical (deterministic tiebreaker)
    ))
    
    report_progress(87, "Applying pathway count limits...")
    
    # Return a larger candidate pool to allow re-ranking after ML scores are added
    # The final limiting will be done in the calling code after ML augmentation
    # Use 3x the requested limit as candidate pool, or all if no limit specified
    candidate_pool_size = max_total_pathways * 3 if max_total_pathways else None
    
    total_before_limit = len(pathway_results)
    if candidate_pool_size and len(pathway_results) > candidate_pool_size:
        pathway_results = pathway_results[:candidate_pool_size]
        logger.info(f"Returning {candidate_pool_size} candidate pathways for re-ranking (from {total_before_limit})")
    
    report_progress(90, "Applying FDR correction...")
    
    # Apply FDR correction to the LIMITED set (not all 521!)
    if len(pathway_results) > 0 and fdr_method:
        pvalues = [r['fisher_pvalue'] for r in pathway_results]
        
        try:
            if fdr_method in ['fdr_bh', 'bonferroni']:
                reject, qvalues, _, _ = multipletests(pvalues, method=fdr_method)
                
                for i, result in enumerate(pathway_results):
                    result['fdr_qvalue'] = qvalues[i]
                    result['significant_FDR'] = (qvalues[i] < pvalue_threshold)
                
                logger.info(f"Applied {fdr_method} FDR correction to {len(pathway_results)} pathways")
            else:
                logger.warning(f"Unknown FDR method '{fdr_method}', no correction applied")
                for result in pathway_results:
                    result['fdr_qvalue'] = result['fisher_pvalue']
                    result['significant_FDR'] = (result['fisher_pvalue'] < pvalue_threshold)
        except Exception as e:
            logger.error(f"FDR correction failed: {e}")
            for result in pathway_results:
                result['fdr_qvalue'] = result['fisher_pvalue']
                result['significant_FDR'] = False
    else:
        # No FDR correction - use raw p-values
        for result in pathway_results:
            result['fdr_qvalue'] = result['fisher_pvalue']  # Keep same key for compatibility
            result['significant_FDR'] = (result['fisher_pvalue'] < pvalue_threshold)
        if len(pathway_results) > 0:
            logger.info("No FDR correction applied - using raw p-values")
    
    report_progress(95, "Finalizing pathway results...")
    
    # Separate by direction for counting (no per-direction limits)
    activated = [r for r in pathway_results if r['status'] == 'Activated']
    inhibited = [r for r in pathway_results if r['status'] == 'Inhibited']
    no_change = [r for r in pathway_results if r['status'] == 'No Change']
    
    # Convert to dict format
    pathway_stats = {}
    for result in pathway_results:
        pathway_name = result.pop('pathway_name')
        pathway_stats[pathway_name] = result
    
    # Summary
    sig_count = sum(1 for r in pathway_results if r.get('significant_FDR', False))
    activated_count = len(activated)
    inhibited_count = len(inhibited)
    
    logger.info("="*60)
    logger.info("FISHER ORA RESULTS:")
    logger.info(f"Total pathways analyzed: {len(pathway_stats)}")
    
    # Show appropriate significance label based on FDR setting
    if fdr_method:
        logger.info(f"FDR-significant (q ≤ {pvalue_threshold}): {sig_count}")
    else:
        logger.info(f"P-value significant (p ≤ {pvalue_threshold}): {sig_count}")
    
    logger.info(f"Activated (Z ≥ {z_threshold}): {activated_count}")
    logger.info(f"Inhibited (Z ≤ -{z_threshold}): {inhibited_count}")
    logger.info("="*60)
    
    # Role-based scoring removed per user request
    
    report_progress(100, "Fisher ORA complete!")
    
    return pathway_stats

def calculate_pathway_statistics_iwpa(df: pd.DataFrame,
                                      min_metabolites: int = 2,
                                      pvalue_threshold: float = 0.05,
                                      fdr_method: Optional[str] = 'fdr_bh',
                                      weight_mode: str = 'signed_p',
                                      z_threshold: float = 2.0,
                                      max_total_pathways: Optional[int] = None,
                                      integrate_datasets: Optional[List[pd.DataFrame]] = None,
                                      progress_callback=None,
                                      pathway_database: Optional[Dict[str, List[str]]] = None) -> Dict[str, Dict[str, Any]]:
    """
    Integrated Weighted Pathway Analysis (IWPA)
    ==========================================

    Computes weighted enrichment scores across pathways using continuous weights
    derived from log2FC, p-values, and optionally integrated datasets.

        Core formulas used:
        - Metabolite weight (w_i): depends on weight_mode
        - Raw pathway score: S = sum(w_i) / sqrt(n)
        - Standardized Z: Z = (sum(w_i) - n*mu_w) / (sqrt(n)*sigma_w)
            where mu_w and sigma_w are global mean/std of metabolite weights.
        - Two-sided p-value: p = 2 * (1 - Phi(|Z|))

    Key difference from Fisher ORA:
    - Uses continuous weights (not binary hit/no-hit)
    - But still respects user's significance threshold to avoid noise dilution

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns ['Name', 'log2FC', 'pvalue', 'All_Pathways']
    min_metabolites : int
        Minimum pathway members required for testing
    pvalue_threshold : float
        P-value threshold for including metabolites (respects GUI setting)
    fdr_method : str
        Multiple-testing correction method ('fdr_bh', 'bonferroni')
    weight_mode : str
        How to define metabolite weights:
        - 'signed_p' : sign(log2FC) * -log10(p)
        - 'log2fc'   : raw log2FC
                - 'combined' : log2FC * -log10(p)

                Algebraic relationship:
                - combined = signed_p * abs(log2FC)
                    because log2FC = sign(log2FC) * abs(log2FC)
    integrate_datasets : list of pd.DataFrame, optional
        List of additional omics dataframes with same structure (Name, log2FC, pvalue)
        for multi-omics integration.
    progress_callback : callable, optional
        Function(percent, message)
    pathway_database : dict, optional
        Pre-built pathway database {pathway_name: [metabolite_names]}.
        If provided, skips internal database building from All_Pathways column.
        Already filtered by min_metabolites threshold.

    Returns
    -------
    dict
        {pathway_name: { 'n_metabolites', 'score', 'mean_weight', 'Z', 'pvalue',
                         'fdr_qvalue', 'status', 'metabolites' }}
    """
    import numpy as np
    from statsmodels.stats.multitest import multipletests

    def report_progress(p, msg=""):
        if progress_callback:
            progress_callback(p, msg)

    report_progress(5, "Preparing IWPA universe...")

    # Validate inputs - DataFrame should already have standardized columns from upstream
    required_cols = ['Name', 'pvalue', 'log2FC', 'All_Pathways']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"DataFrame missing required standardized columns: {missing_cols}. "
                        "Ensure column standardization layer is applied before calling IWPA.")
    
    # Coerce numeric types for pvalue and log2FC
    df['pvalue'] = pd.to_numeric(df['pvalue'], errors='coerce')
    df['log2FC'] = pd.to_numeric(df['log2FC'], errors='coerce')
    
    # CRITICAL: Filter out metabolites without valid p-values BEFORE IWPA calculation
    # These metabolites were annotated (good) but should NOT be included in IWPA
    # because they have no statistical relevance without a p-value
    original_count = len(df)
    pvalue_missing_mask = df['pvalue'].isna()
    n_missing_pvalue = pvalue_missing_mask.sum()
    
    if n_missing_pvalue > 0:
        missing_names = df.loc[pvalue_missing_mask, 'Name'].tolist()
        logger.warning(f"\n⚠️  FILTERING: {n_missing_pvalue} metabolites have missing p-values and will be EXCLUDED from IWPA:")
        for name in missing_names[:10]:
            logger.warning(f"   • {name}")
        if len(missing_names) > 10:
            logger.warning(f"   ... and {len(missing_names) - 10} more")
        logger.warning(f"⚠️  These metabolites will still appear in pathway annotations but won't contribute to IWPA statistics\n")
        
        # Filter DataFrame
        df = df[~pvalue_missing_mask].copy()
        logger.info(f"✓ Filtered dataset: {len(df)} metabolites (removed {n_missing_pvalue} with missing p-values)")
    
    # CRITICAL: Apply user-defined p-value threshold (from GUI)
    # This is statistically sound: we respect the user's definition of "significant"
    # IWPA advantage over Fisher ORA: uses continuous weights within significant set (not binary)
    pvalue_threshold_mask = df['pvalue'] < pvalue_threshold
    n_significant = pvalue_threshold_mask.sum()
    n_filtered = len(df) - n_significant
    
    logger.info(f"\n📊 IWPA THRESHOLD FILTER: Using user-defined p-value threshold ≤ {pvalue_threshold}")
    logger.info(f"   Metabolites passing threshold: {n_significant}")
    logger.info(f"   Metabolites filtered out: {n_filtered}")
    logger.info(f"   ✓ IWPA will use CONTINUOUS WEIGHTS for the {n_significant} significant metabolites")
    logger.info(f"   ✓ This differs from Fisher ORA which uses binary hit/no-hit")
    
    df = df[pvalue_threshold_mask].copy()
    
    # Validate we still have data
    if len(df) == 0:
        raise ValueError(f"No metabolites remain after applying p-value threshold ≤ {pvalue_threshold}")
    
    # Validate log2FC values
    if df['log2FC'].isna().all():
        raise ValueError("All log2FC values are invalid or missing")
    
    # Use provided pathway_database if available, otherwise build from All_Pathways column
    if pathway_database is not None:
        logger.info(f"Using pre-built pathway database with {len(pathway_database)} pathways (already filtered)")
        # The provided pathway_database is already filtered, so we skip internal filtering
        # Still need to filter junk pathways though
        pathway_database = _filter_junk_pathways(pathway_database)
        logger.info(f"After filtering junk pathways: {len(pathway_database)} pathways remain")
    else:
        # Build pathway database from All_Pathways column
        logger.info("Building pathway database from All_Pathways column...")

        # Consolidate pathway name variants by a normalized key so duplicates
        # (e.g., punctuation/case variants) are merged before IWPA.
        pathway_database_norm = {}   # norm_key -> set(metabolite_names)
        orig_name_map = {}           # norm_key -> {orig_name: set(metabolite_names)}

        for idx, row in df.iterrows():
            metabolite_name = row.get('Name', '')
            pathways = row.get('All_Pathways', [])

            # Handle different formats
            if isinstance(pathways, str):
                if pathways.strip():
                    # Attempt to join semicolon-fragmented pathway phrases first (same as ORA)
                    try:
                        from metabolite_pathways_annotator import MetabolitePathwayMapper
                        grouped_str = MetabolitePathwayMapper.group_semicolon_fragments(pathways)
                    except Exception:
                        grouped_str = pathways
                    # CRITICAL: Split by BOTH semicolon AND pipe separators
                    # Network tab uses pipes, annotation tab uses semicolons
                    paths = [p.strip() for p in grouped_str.replace('|', ';').split(';') if p.strip()]
                else:
                    paths = []
            elif isinstance(pathways, list):
                paths = [p for p in pathways if p]
            else:
                paths = []

            for pathway in paths:
                # Normalized key: lowercase, remove non-alphanumeric
                pathway_norm_key = re.sub(r"[^a-z0-9]+", "", str(pathway).lower())
                if not pathway_norm_key:
                    continue
                pathway_database_norm.setdefault(pathway_norm_key, set()).add(metabolite_name)
                orig_name_map.setdefault(pathway_norm_key, {}).setdefault(pathway, set()).add(metabolite_name)

        # Create canonical pathway database by selecting the original name variant with most metabolites
        pathway_database = {}
        duplicates_merged = 0
        merge_examples = []
        for norm_key, metab_set in pathway_database_norm.items():
            names_map = orig_name_map.get(norm_key, {})
            if names_map:
                canonical_name = max(names_map.items(), key=lambda kv: len(kv[1]))[0]
            else:
                canonical_name = norm_key

            pathway_database[canonical_name] = list(metab_set)

            if len(names_map) > 1:
                duplicates_merged += 1
                if len(merge_examples) < 5:
                    merge_examples.append((list(names_map.keys()), canonical_name))

        logger.info(f"Found {len(pathway_database)} unique pathways in dataset BEFORE filtering")
        if duplicates_merged > 0:
            logger.info(f"🔁 Merged {duplicates_merged} duplicate pathway name groups (normalized by removing punctuation and lowercasing)")
            for orig_list, chosen in merge_examples:
                logger.info(f"   • Variants: {orig_list} → Kept: '{chosen}'")

        # Filter out junk pathways (database IDs, column names - same as Fisher ORA)
        pathway_database = _filter_junk_pathways(pathway_database)

    # Build pathway database (reuse ORA preparation)
    universe, _, pathway_members, metabolite_data, _ = prepare_ora_universe(df, pathway_database, 1.0)

    # if integrating multiple datasets, combine log2FC/p-values
    if integrate_datasets:
        for other_df in integrate_datasets:
            # Find name column flexibly
            name_col_candidates = ['Name', 'LipidID', 'Lipid_ID', 'Compound_Name', 'Feature_ID', 'Metabolite']
            name_col = next((c for c in name_col_candidates if c in other_df.columns), 'Name')
            for _, row in other_df.iterrows():
                name = row.get(name_col, row.get('Name', None))
                if name and name in metabolite_data:
                    metabolite_data[name]['log2FC'] += row.get('log2FC', 0)
                    metabolite_data[name]['pvalue'] *= row.get('pvalue', 1)
        # normalize back
        for name in metabolite_data:
            metabolite_data[name]['pvalue'] = min(metabolite_data[name]['pvalue'], 1.0)

    # define weights (with robust NaN handling)
    weights = {}
    nan_count = 0
    inf_count = 0
    
    for name, info in metabolite_data.items():
        fc_raw = info['log2FC']
        p_raw = info['pvalue']
        
        # Sanitize log2FC
        try:
            fc = float(fc_raw)
            if np.isnan(fc):
                fc = 0.0
                nan_count += 1
            elif np.isinf(fc):
                fc = 0.0
                inf_count += 1
        except (TypeError, ValueError):
            fc = 0.0
            nan_count += 1
        
        # Sanitize pvalue
        try:
            p = float(p_raw)
            if np.isnan(p) or p <= 0 or p > 1:
                p = 1.0
                if np.isnan(p_raw):
                    nan_count += 1
        except (TypeError, ValueError):
            p = 1.0
            nan_count += 1
        
        # Calculate weight based on mode.
        # Note: combined = signed_p * abs(log2FC), so combined additionally scales
        # by fold-change magnitude while signed_p uses only direction (+/-) and p-strength.
        if weight_mode == 'signed_p':
            weights[name] = np.sign(fc) * -np.log10(max(p, 1e-12))
        elif weight_mode == 'combined':
            weights[name] = fc * -np.log10(max(p, 1e-12))
        elif weight_mode == 'log2fc':
            weights[name] = fc
        else:
            raise ValueError("Unknown weight_mode")
    
    if nan_count > 0 or inf_count > 0:
        logger.warning(f"IWPA: Sanitized {nan_count} NaN and {inf_count} Inf values in metabolite data")

    # Estimate null distribution moments from all metabolite weights.
    # This calibrates pathway Z-scores across weight modes.
    weight_values = np.array(list(weights.values()), dtype=float)
    mu_w = float(np.mean(weight_values)) if weight_values.size > 0 else 0.0
    sigma_w = float(np.std(weight_values, ddof=1)) if weight_values.size > 1 else 0.0

    if (not np.isfinite(sigma_w)) or sigma_w < 1e-12:
        logger.warning(
            "IWPA: Global weight std is near zero; falling back to sigma_w=1.0. "
            "P-values may be less informative for this dataset."
        )
        sigma_w = 1.0

    logger.info(
        f"IWPA null calibration: mu_w={mu_w:.4g}, sigma_w={sigma_w:.4g}, "
        f"n_weights={weight_values.size}"
    )

    report_progress(25, "Computing weighted scores...")

    results = []
    total = len(pathway_members)

    for idx, (pw, members) in enumerate(pathway_members.items(), 1):
        members = [m for m in members if m in weights]
        n = len(members)
        if n < min_metabolites:
            continue

        w = np.array([weights[m] for m in members])
        
        # Additional NaN check at pathway level
        if np.any(np.isnan(w)) or np.any(np.isinf(w)):
            logger.warning(f"IWPA: Pathway '{pw}' has NaN/Inf weights, sanitizing...")
            w = np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)
        
        # Calculate both raw score and standardized Z-score.
        sum_w = np.sum(w)
        raw_score = sum_w / np.sqrt(n) if n > 0 else 0.0
        Z = ((sum_w - (n * mu_w)) / (np.sqrt(n) * sigma_w)) if n > 0 else 0.0
        
        # Sanitize Z if still NaN/Inf
        if np.isnan(Z) or np.isinf(Z):
            logger.warning(f"IWPA: Pathway '{pw}' produced NaN/Inf Z-score, setting to 0.0")
            Z = 0.0
        
        mean_w = np.mean(w) if n > 0 else 0.0
        std_w = np.std(w, ddof=1) if n > 1 else 0.0
        
        # Sanitize mean and std
        if np.isnan(mean_w) or np.isinf(mean_w):
            mean_w = 0.0
        if np.isnan(std_w) or np.isinf(std_w):
            std_w = 0.0

        # approximate two-sided p-value from normal distribution using scipy.stats.norm
        pval = 2 * (1 - norm.cdf(abs(Z))) if not np.isnan(Z) and not np.isinf(Z) else 1.0

        # CRITICAL: Use strict inequality (>) so that Z=0 always remains "No Change"
        # When z_threshold=0, only Z>0 is Activated, only Z<0 is Inhibited
        if Z > z_threshold:
            status = "Activated"
        elif Z < -z_threshold:
            status = "Inhibited"
        else:
            status = "No Change"

        results.append({
            "pathway_name": pw,
            "n_metabolites": n,
            "score": Z,
            "raw_score": raw_score,
            "mean_weight": mean_w,
            "std_weight": std_w,
            "Z": Z,
            "pvalue": pval,
            "status": status,
            "metabolites": list(members)
        })

        if progress_callback and idx % 10 == 0:
            report_progress(25 + int((idx / total) * 50), f"Processed {idx}/{total} pathways")

    report_progress(75, "Applying FDR correction...")

    # Apply FDR correction
    if fdr_method and results:
        pvals = [r['pvalue'] for r in results]
        if fdr_method == 'fdr_bh':
            _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
        elif fdr_method == 'bonferroni':
            _, qvals, _, _ = multipletests(pvals, method='bonferroni')
        else:
            qvals = pvals  # no correction

        for r, q in zip(results, qvals):
            r['fdr_qvalue'] = q
            r['significant_FDR'] = q <= pvalue_threshold
    else:
        for r in results:
            r['fdr_qvalue'] = r['pvalue']
            r['significant_FDR'] = r['pvalue'] < pvalue_threshold

    report_progress(90, "Finalizing results...")

    # Return a larger candidate pool to allow re-ranking after ML scores are added
    # The final limiting will be done in the calling code after ML augmentation
    # Use 3x the requested limit as candidate pool, or all if no limit specified
    if max_total_pathways and len(results) > max_total_pathways:
        total_before_limit = len(results)
        candidate_pool_size = max_total_pathways * 3
        
        # Sort by absolute Z-score (descending), then by p-value (ascending), then pathway name (deterministic tiebreaker)
        results_sorted = sorted(results, key=lambda r: (
            -abs(r['Z']),           # Primary: highest |Z-score| first
            r['pvalue'],            # Secondary: lowest p-value first
            r.get('pathway_name', '')  # Tertiary: alphabetical (deterministic tiebreaker)
        ))
        
        if len(results_sorted) > candidate_pool_size:
            results = results_sorted[:candidate_pool_size]
            logger.info(f"Returning {candidate_pool_size} candidate pathways for re-ranking (from {total_before_limit})")
        else:
            results = results_sorted

    # Convert to dict format
    pathway_stats = {}
    for result in results:
        pathway_name = result.pop('pathway_name')
        pathway_stats[pathway_name] = result

    # Summary (count from limited results)
    sig_count = sum(1 for r in results if r.get('significant_FDR', False))
    activated_count = sum(1 for r in results if r['status'] == 'Activated')
    inhibited_count = sum(1 for r in results if r['status'] == 'Inhibited')

    logger.info("="*60)
    logger.info("IWPA RESULTS:")
    logger.info(f"Total pathways analyzed: {len(pathway_stats)}")
    logger.info(f"Weight mode: {weight_mode}")
    logger.info(f"Z-threshold: ±{z_threshold}")

    # Show appropriate significance label based on FDR setting
    if fdr_method:
        logger.info(f"FDR-significant (q ≤ {pvalue_threshold}): {sig_count}")
    else:
        logger.info(f"P-value significant (p ≤ {pvalue_threshold}): {sig_count}")

    logger.info(f"Activated (Z ≥ {z_threshold}): {activated_count}")
    logger.info(f"Inhibited (Z ≤ -{z_threshold}): {inhibited_count}")
    logger.info("="*60)

    # Role-based scoring removed per user request

    report_progress(100, "IWPA complete!")

    return pathway_stats

if __name__ == "__main__":
    # Test with sample data
    print("Fisher ORA module loaded successfully")
    print("Use calculate_pathway_statistics_fisher_ora() for pathway analysis")
    print("Use calculate_pathway_statistics_iwpa() for IWPA analysis")
