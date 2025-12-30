#!/usr/bin/env python3
"""
Calculate Species-Specific Universe Sizes (U) for Fisher ORA from Database Files

This script calculates the total number of UNIQUE metabolites with pathway annotations
for each species (Homo sapiens, Mus musculus, Rattus norvegicus) across all offline databases.

Universe Definition:
- U = Unique metabolites across all databases with at least 1 pathway annotation
- Deduplication uses multiple identifiers (HMDB_ID, KEGG_ID, PubChem_CID, ChEBI_ID, Name, InChIKey)
- Only metabolites with pathway annotations are counted
- Calculated separately for each species

Author: Metabolite Analysis Tool
Date: November 15, 2025
"""

import pandas as pd
import os
import json
import logging
from pathlib import Path
from typing import Set, Dict, List, Tuple
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def normalize_id(value):
    """Normalize an ID value for comparison."""
    if pd.isna(value) or not value:
        return None
    value_str = str(value).strip().upper()
    if value_str in ['NAN', 'NONE', '', 'NULL']:
        return None
    return value_str


def normalize_name(name: str) -> str:
    """Normalize metabolite name for comparison."""
    if not isinstance(name, str) or not name or name.lower() == 'nan':
        return ''
    return name.strip().lower().replace('-', '').replace(' ', '')


def extract_metabolites_hmdb(df: pd.DataFrame, species: str = 'Homo sapiens') -> Dict[str, Set[str]]:
    """Extract metabolites from HMDB database (has Pathways column)."""
    metabolites = {}
    
    for idx, row in df.iterrows():
        # Check if has pathway
        pathways = row.get('Pathways')
        if pd.isna(pathways) or not pathways or str(pathways).strip() == '':
            continue
        
        # Collect identifiers
        identifiers = set()
        
        hmdb_id = normalize_id(row.get('HMDB'))
        if hmdb_id:
            identifiers.add(f"HMDB:{hmdb_id}")
        
        name = row.get('Name')
        if name and pd.notna(name):
            identifiers.add(f"NAME:{normalize_name(name)}")
        
        inchikey = normalize_id(row.get('InChIKey'))
        if inchikey:
            identifiers.add(f"INCHIKEY:{inchikey}")
        
        cas = normalize_id(row.get('CAS'))
        if cas:
            identifiers.add(f"CAS:{cas}")
        
        if identifiers:
            met_key = list(identifiers)[0]  # Primary key
            if met_key not in metabolites:
                metabolites[met_key] = identifiers
            else:
                metabolites[met_key].update(identifiers)
    
    logger.info(f"  HMDB: {len(metabolites)} metabolites with pathways")
    return metabolites


def extract_metabolites_smp(df: pd.DataFrame, species: str = 'Homo sapiens') -> Dict[str, Set[str]]:
    """Extract metabolites from merged_SMP_metabolites (has SMP_Pathways column)."""
    metabolites = {}
    
    for idx, row in df.iterrows():
        # Check if has pathway
        pathways = row.get('SMP_Pathways')
        if pd.isna(pathways) or not pathways or str(pathways).strip() == '':
            continue
        
        # Collect identifiers
        identifiers = set()
        
        hmdb_id = normalize_id(row.get('HMDB'))
        if hmdb_id:
            identifiers.add(f"HMDB:{hmdb_id}")
        
        kegg_id = normalize_id(row.get('KEGG'))
        if kegg_id:
            identifiers.add(f"KEGG:{kegg_id}")
        
        name = row.get('Name')
        if name and pd.notna(name):
            identifiers.add(f"NAME:{normalize_name(name)}")
        
        inchikey = normalize_id(row.get('InChIKey'))
        if inchikey:
            identifiers.add(f"INCHIKEY:{inchikey}")
        
        cas = normalize_id(row.get('CAS'))
        if cas:
            identifiers.add(f"CAS:{cas}")
        
        if identifiers:
            met_key = list(identifiers)[0]
            if met_key not in metabolites:
                metabolites[met_key] = identifiers
            else:
                metabolites[met_key].update(identifiers)
    
    logger.info(f"  SMP: {len(metabolites)} metabolites with pathways")
    return metabolites


def extract_metabolites_pathbank(df: pd.DataFrame, species: str) -> Dict[str, Set[str]]:
    """Extract metabolites from PathBank (each row = metabolite-pathway pair, has Species column)."""
    # Filter by species first
    df_species = df[df['Species'] == species].copy()
    
    metabolites = {}
    
    for idx, row in df_species.iterrows():
        # Every row has a pathway (that's the structure)
        identifiers = set()
        
        hmdb_id = normalize_id(row.get('HMDB ID'))
        if hmdb_id:
            identifiers.add(f"HMDB:{hmdb_id}")
        
        kegg_id = normalize_id(row.get('KEGG ID'))
        if kegg_id:
            identifiers.add(f"KEGG:{kegg_id}")
        
        chebi_id = normalize_id(row.get('ChEBI ID'))
        if chebi_id:
            identifiers.add(f"CHEBI:{chebi_id}")
        
        name = row.get('Metabolite Name')
        if name and pd.notna(name):
            identifiers.add(f"NAME:{normalize_name(name)}")
        
        inchikey = normalize_id(row.get('InChI Key'))
        if inchikey:
            identifiers.add(f"INCHIKEY:{inchikey}")
        
        cas = normalize_id(row.get('CAS'))
        if cas:
            identifiers.add(f"CAS:{cas}")
        
        if identifiers:
            met_key = list(identifiers)[0]
            if met_key not in metabolites:
                metabolites[met_key] = identifiers
            else:
                metabolites[met_key].update(identifiers)
    
    logger.info(f"  PathBank ({species}): {len(metabolites)} metabolites with pathways")
    return metabolites


def extract_metabolites_wikipathways(df: pd.DataFrame, species: str) -> Dict[str, Set[str]]:
    """Extract metabolites from WikiPathways (has Pathways column and Organism column)."""
    # WikiPathways files are already species-specific, but check Organism column
    df_species = df[df['Organism'] == species].copy() if 'Organism' in df.columns else df.copy()
    
    metabolites = {}
    
    for idx, row in df_species.iterrows():
        # Check if has pathway
        pathways = row.get('Pathways')
        if pd.isna(pathways) or not pathways or str(pathways).strip() == '':
            continue
        
        # Collect identifiers
        identifiers = set()
        
        hmdb_id = normalize_id(row.get('HMDB_ID'))
        if hmdb_id:
            identifiers.add(f"HMDB:{hmdb_id}")
        
        kegg_id = normalize_id(row.get('KEGG'))
        if kegg_id:
            identifiers.add(f"KEGG:{kegg_id}")
        
        chebi_id = normalize_id(row.get('ChEBI'))
        if chebi_id:
            identifiers.add(f"CHEBI:{chebi_id}")
        
        pubchem_id = normalize_id(row.get('PubChem'))
        if pubchem_id:
            identifiers.add(f"PUBCHEM:{pubchem_id}")
        
        name = row.get('name')
        if name and pd.notna(name):
            identifiers.add(f"NAME:{normalize_name(name)}")
        
        inchikey = normalize_id(row.get('InChIKey'))
        if inchikey:
            identifiers.add(f"INCHIKEY:{inchikey}")
        
        cas = normalize_id(row.get('CAS'))
        if cas:
            identifiers.add(f"CAS:{cas}")
        
        if identifiers:
            met_key = list(identifiers)[0]
            if met_key not in metabolites:
                metabolites[met_key] = identifiers
            else:
                metabolites[met_key].update(identifiers)
    
    logger.info(f"  WikiPathways ({species}): {len(metabolites)} metabolites with pathways")
    return metabolites


def deduplicate_metabolites(metabolite_info: Dict[str, Set[str]]) -> Set[str]:
    """
    Deduplicate metabolites by exact match only (UNION approach).
    
    Keep ALL metabolites from all databases.
    Only remove metabolites that are EXACTLY the same across databases.
    
    Matching criteria: metabolite must have EXACT same identifiers (not just one match)
    This gives us a TRUE union: unique_count >= max(database_sizes)
    """
    logger.info("\n" + "="*70)
    logger.info("DEDUPLICATION PHASE - UNION Approach")
    logger.info("="*70)
    logger.info(f"Starting with {len(metabolite_info):,} metabolite entries...")
    logger.info("Strategy: Keep ALL metabolites, remove only exact duplicates")
    
    # Convert sets to frozensets for hashability (to use as dict keys)
    # This allows us to group metabolites with identical identifier sets
    identifier_sets_map = {}  # frozenset of identifiers -> list of metabolite keys
    
    logger.info("\nGrouping metabolites by identical identifier sets...")
    processed = 0
    
    for met_key, identifiers in metabolite_info.items():
        processed += 1
        if processed % 10000 == 0:
            logger.info(f"  Progress: {processed:,}/{len(metabolite_info):,}")
        
        # Create frozen set for comparison
        id_set = frozenset(identifiers)
        
        if id_set not in identifier_sets_map:
            identifier_sets_map[id_set] = []
        identifier_sets_map[id_set].append(met_key)
    
    logger.info(f"  Found {len(identifier_sets_map):,} unique identifier combinations")
    
    # For each unique identifier set, keep only ONE representative
    logger.info("\nSelecting representatives from duplicate groups...")
    unique_metabolites = set()
    duplicates_removed = 0
    
    for id_set, met_keys in identifier_sets_map.items():
        # Keep first metabolite, rest are duplicates
        unique_metabolites.add(met_keys[0])
        
        if len(met_keys) > 1:
            duplicates_removed += len(met_keys) - 1
    
    logger.info(f"  Completed selection")
    
    logger.info("\n" + "="*70)
    logger.info("DEDUPLICATION RESULTS")
    logger.info("="*70)
    logger.info(f"  Before: {len(metabolite_info):,} metabolites")
    logger.info(f"  After: {len(unique_metabolites):,} unique metabolites")
    logger.info(f"  Removed: {duplicates_removed:,} exact duplicates ({100*duplicates_removed/len(metabolite_info):.1f}%)")
    logger.info(f"  Expected: Unique >= max(database_sizes)")
    
    return unique_metabolites


def calculate_species_universe(database_folder: str, species: str) -> int:
    """
    Calculate universe size for a specific species.
    
    Parameters:
    -----------
    database_folder : str
        Path to folder containing database files
    species : str
        Species name ('Homo sapiens', 'Mus musculus', or 'Rattus norvegicus')
    
    Returns:
    --------
    int
        Universe size for this species
    """
    logger.info(f"\n{'='*70}")
    logger.info(f"CALCULATING UNIVERSE FOR: {species}")
    logger.info(f"{'='*70}")
    
    db_path = Path(database_folder)
    all_metabolites = {}
    raw_totals = {}  # Track raw counts before deduplication
    
    # 1. HMDB (species-agnostic - add to ALL species)
    hmdb_file = db_path / "hmdb_database.feather"
    if hmdb_file.exists():
        logger.info("Processing HMDB database...")
        df = pd.read_feather(hmdb_file)
        mets = extract_metabolites_hmdb(df, species)
        raw_totals['HMDB'] = len(mets)
        all_metabolites.update(mets)
    
    # 2. SMP (species-agnostic - add to ALL species)
    smp_file = db_path / "merged_SMP_metabolites.feather"
    if smp_file.exists():
        logger.info("Processing SMP database...")
        df = pd.read_feather(smp_file)
        mets = extract_metabolites_smp(df, species)
        raw_totals['SMP'] = len(mets)
        all_metabolites.update(mets)
    
    # 3. PathBank (has species column)
    pathbank_file = db_path / "pathbank_selected.feather"
    if pathbank_file.exists():
        logger.info("Processing PathBank database...")
        df = pd.read_feather(pathbank_file)
        mets = extract_metabolites_pathbank(df, species)
        raw_totals['PathBank'] = len(mets)
        all_metabolites.update(mets)
    
    # 4. WikiPathways (species-specific files)
    species_file_map = {
        'Homo sapiens': 'wikipathways_homo_sapiens.feather',
        'Mus musculus': 'wikipathways_mus_musculus.feather',
        'Rattus norvegicus': 'wikipathways_rattus_norvegicus.feather'
    }
    
    wiki_file = db_path / species_file_map.get(species, '')
    if wiki_file.exists():
        logger.info(f"Processing WikiPathways database...")
        df = pd.read_feather(wiki_file)
        mets = extract_metabolites_wikipathways(df, species)
        raw_totals['WikiPathways'] = len(mets)
        all_metabolites.update(mets)
    
    # Calculate raw sum and overlap
    raw_sum = sum(raw_totals.values())
    unique_before_dedup = len(all_metabolites)
    duplicates_before_dedup = raw_sum - unique_before_dedup
    
    logger.info(f"\nDatabase Summary (BEFORE deduplication):")
    for db_name, count in raw_totals.items():
        logger.info(f"  {db_name:15s}: {count:6,} metabolites")
    logger.info(f"\nRaw sum (all databases):     {raw_sum:,} metabolites")
    logger.info(f"After removing duplicates:   {unique_before_dedup:,} metabolites")
    logger.info(f"Duplicates (same name/ID):   {duplicates_before_dedup:,} ({100*duplicates_before_dedup/raw_sum:.1f}%)")
    
    # Deduplicate
    unique_metabolites = deduplicate_metabolites(all_metabolites)
    universe_size = len(unique_metabolites)
    
    logger.info(f"Universe size (U) for {species}: {universe_size}")
    
    return universe_size


def calculate_all_universes(database_folder: str = None) -> Dict[str, int]:
    """
    Calculate universe sizes for all three species.
    
    Parameters:
    -----------
    database_folder : str, optional
        Path to folder containing database files.
        If None, looks for 'Databases' folder in parent directory (MetaboGraph root)
    
    Returns:
    --------
    dict
        {species_name: universe_size}
    """
    logger.info("="*70)
    logger.info("METABOLITE UNIVERSE CALCULATOR - MULTI-SPECIES")
    logger.info("="*70)
    
    # Default to parent directory's Databases folder
    if database_folder is None:
        # Get parent folder (MetaboGraph root) from this script location
        script_dir = Path(__file__).parent
        db_path = script_dir.parent / "Databases"
    else:
        db_path = Path(database_folder)
    
    logger.info(f"Looking for database files in: {db_path}")
    
    if not db_path.exists():
        logger.error(f"Database folder not found: {db_path}")
        return {}
    
    species_list = ['Homo sapiens', 'Mus musculus', 'Rattus norvegicus']
    universes = {}
    
    for species in species_list:
        universe_size = calculate_species_universe(db_path, species)
        universes[species] = universe_size
    
    logger.info("\n" + "="*70)
    logger.info("FINAL RESULTS - ALL SPECIES")
    logger.info("="*70)
    for species, size in universes.items():
        logger.info(f"  {species:20s}: U = {size:6d} metabolites")
    
    return universes


def save_universe_config(universes: Dict[str, int], config_file: str = "universe_config.json"):
    """
    Save universe sizes to configuration file.
    
    Parameters:
    -----------
    universes : dict
        {species_name: universe_size}
    config_file : str
        Path to configuration file
    """
    script_dir = Path(__file__).parent
    config_path = script_dir / config_file
    
    # Also save to parent directory (MetaboGraph root) for easier access
    root_config_path = script_dir.parent / config_file
    
    config = {
        'universes': universes,
        'calculated_date': pd.Timestamp.now().isoformat(),
        'description': 'Species-specific universe sizes (U) for Fisher ORA',
        'note': 'Automatically recalculated when databases are updated'
    }
    
    # Save to main_script folder
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)
    
    # Save to root folder
    with open(root_config_path, 'w') as f:
        json.dump(config, f, indent=2)
    
    logger.info(f"\n{'='*70}")
    logger.info("CONFIGURATION SAVED")
    logger.info(f"{'='*70}")
    logger.info(f"File: {config_path}")
    logger.info(f"File: {root_config_path}")
    for species, size in universes.items():
        logger.info(f"  {species}: U = {size}")


if __name__ == "__main__":
    print("\n" + "="*70)
    print("METABOLITE UNIVERSE CALCULATOR - MULTI-SPECIES")
    print("="*70)
    print("\nCalculating universe (U) for all three species:")
    print("  - Homo sapiens")
    print("  - Mus musculus")
    print("  - Rattus norvegicus\n")
    
    # Calculate universes for all species
    universes = calculate_all_universes()
    
    if universes:
        # Save to config file
        save_universe_config(universes)
        
        print("\n" + "="*70)
        print("SUCCESS")
        print("="*70)
        print("Universe sizes calculated and saved to universe_config.json")
        print("\nThese values will be used for Fisher ORA analysis.")
        print("Rerun this script whenever you update your databases.")
        sys.exit(0)  # Success exit code
    else:
        print("\n" + "="*70)
        print("ERROR")
        print("="*70)
        print("Failed to calculate universe sizes.")
        print("Please check that database files exist in the 'Databases' folder.")
        sys.exit(1)  # Error exit code
        sys.exit(1)  # Error exit code
