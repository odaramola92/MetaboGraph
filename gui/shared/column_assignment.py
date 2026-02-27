"""
Unified Column Assignment Dialog for all tabs.

This module provides a centralized column verification and assignment system
that can be used across all tabs (Data Cleaning, Statistics, Visualization,
Pathway Annotation, ID Annotation, etc.) with tab-specific requirements.

Features:
- Two-column display (Column | Assignment) with easy user correction
- Dropdown assignments with all available options
- Tab-specific column requirements
- Auto-detection matching the actual algorithms used in:
  * Data Cleaner (feature columns, numeric samples, classifications)
  * Statistics Analysis (feature_columns_canonical, sample detection)
  * Visualization & Pathway Analysis (statistical columns detection)
  * ID Annotation (Feature ID detection, multi-sheet support)
- Auto-calculation of missing columns (log2FC, -log10_pvalue, etc.)
- Persistent column assignments during session
- Multi-sheet support for ID Annotation tab
"""

import tkinter as tk
from tkinter import ttk, messagebox
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Any, Set
import logging
import re
import threading

from gui.shared.utils import is_statistics_metadata_col

logger = logging.getLogger(__name__)

# ===== SHARED DETECTION CONSTANTS =====
# These come from the actual analysis modules to ensure consistency

COMMON_ID_COLUMNS = {
    'Name', 'HMDB_ID', 'KEGG_ID', 'PubChem_CID', 'ChEBI_ID', 'LipidMaps_ID',
    'InChIKey', 'InChI', 'SMILES', 'CAS', 'Formula', 'Endogenous_Source'
}

FEATURE_COLUMNS_CANONICAL = [
    'Name', 'Name_Key', 'Formula', 'Molecular_Formula', 'Molecular Formula', 'MW', 'Molecular_Weight', 'ppm',
    'Reference Ion', 'Reference_Ion', 'MS2', 'm/z', 'RT', 'RT [min]', 'Area (Max.)', 'Polarity', 'MS2_Purity', 'MS2 Purity [%]',
    'LipidMaps_ID', 'PubChem_CID', 'KEGG_ID', 'HMDB_ID', 'ChEBI_ID', 'CAS', 'SMILES', 'InChI', 'InChIKey', 'IUPAC_Name',
    'Super_Class', 'Class', 'Class_name', 'Sub_Class', 'Endogenous_Source', 'Metabolika Pathways', 'BioCyc Pathways'
]

# ===== COLUMN TYPE DEFINITIONS =====

COLUMN_TYPES = {
    'Feature ID': {
        'description': 'Main identifier (Name, Metabolite, Compound ID)',
        'required_by': ['data_cleaning', 'statistics', 'visualization', 'pathway', 'id_annotation'],
        'examples': ['Name', 'Metabolite', 'LipidID', 'Compound'],
    },
    'HMDB ID': {
        'description': 'HMDB database identifier',
        'required_by': [],
        'examples': ['HMDB_ID', 'HMDB', 'hmdb_id'],
        'category': 'Identifier',
    },
    'KEGG ID': {
        'description': 'KEGG database identifier',
        'required_by': [],
        'examples': ['KEGG_ID', 'KEGG', 'kegg_id'],
        'category': 'Identifier',
    },
    'ChEBI ID': {
        'description': 'ChEBI (Chemical Entities of Biological Interest) identifier for Reactome',
        'required_by': [],
        'examples': ['ChEBI_ID', 'ChEBI', 'chebi_id', 'CHEBI'],
        'category': 'Identifier',
    },
    'PubChem CID': {
        'description': 'PubChem compound identifier',
        'required_by': [],
        'examples': ['PubChem_CID', 'PubChem CID', 'pubchem_cid', 'pubchem'],
        'category': 'Identifier',
    },
    'LipidMaps ID': {
        'description': 'LipidMaps structural database identifier',
        'required_by': [],
        'examples': ['LipidMaps_ID', 'LipidMaps ID', 'lipidmaps_id', 'LMID'],
        'category': 'Identifier',
    },
    'InChIKey': {
        'description': 'InChIKey chemical structure identifier',
        'required_by': [],
        'examples': ['InChIKey', 'InChI_Key', 'inchikey'],
        'category': 'Identifier',
    },
    'P-Value': {
        'description': 'Statistical p-value (raw or adjusted)',
        'required_by': ['statistics', 'visualization', 'pathway'],
        'examples': ['pvalue', 'p_value', 'padj', 'p_adj'],
        'category': 'Statistics',
    },
    'Fold Change': {
        'description': 'Fold change (will be converted to log2 if needed)',
        'required_by': [],
        'examples': ['FC', 'fold_change', 'foldchange'],
        'category': 'Statistics',
    },
    'Log2 Fold Change': {
        'description': 'Log2-transformed fold change (preferred over FC)',
        'required_by': ['statistics', 'visualization', 'pathway'],
        'examples': ['log2FC', 'log2_fc', 'logFC', 'LogFC'],
        'category': 'Statistics',
    },
    '-Log10 P-Value': {
        'description': 'Negative log10 p-value (for visualization)',
        'required_by': [],
        'examples': ['neglog10', 'neg_log10', '-log10(pvalue)'],
        'category': 'Statistics',
    },
    'Class': {
        'description': 'Metabolite/Lipid classification',
        'required_by': [],
        'examples': ['Class', 'Metabolite_Class', 'Lipid_Class'],
        'category': 'Classification',
    },
    'Class_name': {
        'description': 'Lipid class name (for class-based network compression)',
        'required_by': [],
        'examples': ['Class_name', 'ClassName', 'Lipid_Class_Name'],
        'category': 'Classification',
    },
    'Super Class': {
        'description': 'Metabolite super class',
        'required_by': [],
        'examples': ['Super_Class', 'SuperClass'],
        'category': 'Classification',
    },
    'Sample Column': {
        'description': 'Numeric sample/intensity data column',
        'required_by': ['statistics', 'visualization'],
        'examples': ['Sample1', 'Intensity', 'Abundance'],
        'category': 'Sample Data',
    },
    'Sample Column': {
        'description': 'Mark columns as sample data (numeric intensity/abundance values). Use this to identify which columns contain your biological sample measurements.',
        'required_by': [],
        'examples': ['Any numeric column with sample names'],
        'category': 'Sample Data',
    },
    'Feature Column': {
        'description': 'Mark columns as feature metadata (Name, Formula, Class, RT, m/z, etc.). Use this to identify which columns contain metabolite/lipid annotation information.',
        'required_by': [],
        'examples': ['Name', 'Formula', 'Class', 'RT', 'm/z', 'MW'],
        'category': 'Feature Data',
    },
    'Metabolika Pathways': {
        'description': 'Metabolika pathway annotations',
        'required_by': [],
        'examples': ['Metabolika Pathways', 'Metabolika'],
        'category': 'Pathway',
    },
    'BioCyc Pathways': {
        'description': 'BioCyc pathway annotations',
        'required_by': [],
        'examples': ['BioCyc Pathways', 'BioCyc'],
        'category': 'Pathway',
    },
    'Group Area': {
        'description': 'Group Area sample columns (metabolite data)',
        'required_by': [],
        'examples': ['Group Area: LF_pos', 'Group Area: Sample1'],
        'category': 'Sample Data',
    },
    'Area': {
        'description': 'Area sample columns (lipid data)',
        'required_by': [],
        'examples': ['Area[DU143_3D_neg_1]', 'Area(Sample1)'],
        'category': 'Sample Data',
    },
    'Grade': {
        'description': 'Grade quality score columns (lipid data)',
        'required_by': [],
        'examples': ['Grade[DU143_3D_neg_1]', 'Grade[s1-1]'],
        'category': 'Sample Data',
    },
    'Rej': {
        'description': 'Rejection flag column (lipid data - True/False)',
        'required_by': [],
        'examples': ['Rej', 'Reject', 'Rejected'],
        'category': 'Lipid Data',
    },
    'LipidID': {
        'description': 'Lipid identifier (required for lipid data)',
        'required_by': ['lipid_cleaning'],
        'examples': ['LipidID', 'Lipid_ID', 'LipidName'],
        'category': 'Lipid Data',
    },
    'AdductIon': {
        'description': 'Adduct ion type for lipid filtering',
        'required_by': [],
        'examples': ['AdductIon', 'Adduct', 'Ion'],
        'category': 'Lipid Data',
    },
    'ObsMz': {
        'description': 'Calculated m/z values (lipid data - CalcMz)',
        'required_by': [],
        'examples': ['CalcMz', 'calcmz'],
        'category': 'Lipid Data',
    },
    'ObsRt': {
        'description': 'Base retention time (lipid data - BaseRt)',
        'required_by': [],
        'examples': ['BaseRt', 'basert'],
        'category': 'Lipid Data',
    },
    'Delta(PPM)': {
        'description': 'PPM difference values (lipid data)',
        'required_by': [],
        'examples': ['Delta(PPM)[Sample1]', 'PPM[Sample1]'],
        'category': 'Lipid Data',
    },
    'Class_name': {
        'description': 'Lipid class name (recommended for lipid ID annotation)',
        'required_by': [],
        'examples': ['Class_name', 'ClassName', 'Class Name', 'Lipid Class Name'],
        'category': 'Lipid Data',
    },
    'QC Sample': {
        'description': 'Quality Control (QC) sample columns - used for LOESS drift correction and monitoring',
        'required_by': [],
        'examples': ['QC', 'QC_1', 'QC_pool', 'Pool', 'Pooled_QC', 'Quality_Control'],
        'category': 'Sample Data',
    },
    'Internal Standard': {
        'description': 'Internal Standard (IS) feature - reference compound added at fixed concentration for IS normalization',
        'required_by': [],
        'examples': ['Internal_Standard', 'IS', 'ISTD', 'Standard', 'Ref_Standard'],
        'category': 'Feature Data',
    },
    # Data Cleaning specific columns
    'Annot. DeltaMass [ppm]': {
        'description': 'Annotation delta mass in PPM (data cleaning)',
        'required_by': [],
        'examples': ['Annot. DeltaMass [ppm]', 'Delta Mass', 'DeltaMass'],
        'category': 'Feature Data',
    },
    'Calc. MW': {
        'description': 'Calculated molecular weight (data cleaning)',
        'required_by': [],
        'examples': ['Calc. MW', 'Calculated MW', 'MW', 'Molecular Weight'],
        'category': 'Feature Data',
    },
    'm/z': {
        'description': 'Mass to charge ratio',
        'required_by': [],
        'examples': ['m/z', 'mz', 'Mass/Charge', 'Precursor m/z'],
        'category': 'Feature Data',
    },
    'MS2 Purity [%]': {
        'description': 'MS2 purity percentage (data cleaning)',
        'required_by': [],
        'examples': ['MS2 Purity [%]', 'MS2Purity', 'Purity'],
        'category': 'Feature Data',
    },
}

# ===== TAB-SPECIFIC REQUIREMENTS =====

TAB_REQUIREMENTS = {
    'data_cleaning': {
        'required': ['Feature ID'],  # Feature ID = Name column (critical)
        'essential': ['Reference Ion', 'Area (Max.)', 'RT [min]', 'Polarity', 'MS2'],  # Essential columns users should confirm
        'optional': ['Formula', 'Annot. DeltaMass [ppm]', 'Calc. MW', 'm/z', 'Reference Ion', 'RT [min]', 
                     'Area (Max.)', 'Metabolika Pathways', 'BioCyc Pathways', 'Polarity', 'MS2', 'MS2 Purity [%]', 
                     'Class', 'Super Class', 'Group Area'],
        'title': 'Data Cleaning Column Assignment',
        'description': 'Verify columns for data cleaning. Confirm essential columns before proceeding.',
        'detect_samples': False,
        'multi_sheet': False,
        'detect_group_area': True,  # Enable Group Area detection
        'column_mapping': {
            'Feature ID': 'Name',  # Feature ID is the metabolite Name
        }
    },
    'statistics': {
        'required': ['Feature ID', 'P-Value', 'Log2 Fold Change'],
        'optional': ['Fold Change', 'Class', 'Super Class', 'HMDB ID', 'KEGG ID'],
        'title': 'Statistics Column Assignment',
        'description': 'Verify columns for statistical analysis. Sample columns will be auto-detected from numeric data.',
        'detect_samples': True,
        'multi_sheet': False,
        'can_calculate': ['Log2 Fold Change'],
    },
    'statistics_metabolite': {
        'required': ['Feature ID'],  # Only Feature ID is required - no P-Value or Log2FC
        'optional': ['Sample Column', 'Feature Column'],  # Only sample/feature identifiers
        'title': 'Metabolite Statistics - Column Verification',
        'description': 'Verify Feature ID and sample columns for metabolite statistical analysis. Feature columns = metadata/annotation columns (Name, Formula, Class, etc.). Sample columns = numeric intensity columns for biological samples.',
        'detect_samples': True,  # Auto-detect numeric sample columns
        'multi_sheet': False,
        'can_calculate': [],  # No calculations - raw data only
        'column_mapping': {
            'Feature ID': 'Name',  # Feature ID is the metabolite Name
        }
    },
    'statistics_lipid': {
        'required': ['LipidID'],  # Only LipidID required
        'optional': ['Class', 'Class_name', 'Sample Column', 'Feature Column'],  # Class and Class_name are optional but important for classification
        'title': 'Lipid Statistics - Column Verification',
        'description': 'Verify LipidID and sample columns for lipid statistical analysis. Feature columns = lipid metadata (LipidID, Class, Class_name, CalcMz, etc.). Sample columns = numeric intensity columns for biological samples.',
        'detect_samples': True,  # Auto-detect numeric sample columns
        'multi_sheet': False,
        'can_calculate': [],  # No calculations - raw data only
        'column_mapping': {
            'LipidID': 'LipidID',  # Ensure LipidID is preserved
        }
    },
    'visualization': {
        'required': ['Feature ID', 'P-Value', 'Log2 Fold Change'],
        'optional': ['Log2 Fold Change', 'Fold Change', '-Log10 P-Value', 'Class', 'Super Class', 'HMDB ID', 'KEGG ID'],
        'title': 'Visualization Column Assignment',
        'description': 'Verify columns for data visualization. Sample columns will be auto-detected.',
        'detect_samples': True,
        'multi_sheet': False,
        'can_calculate': ['Log2 Fold Change', '-Log10 P-Value'],
    },
    'pathway': {
        'required': ['Feature ID', 'P-Value', 'Log2 Fold Change'],
        'optional': ['Fold Change', 'Class', 'Class_name', 'Super Class'],
        'hardcoded': ['HMDB_Pathways', 'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways', 'Metabolika_Pathways', 'Reactome_Pathways', 'KEGG_Pathways', 'All_Pathways_Display'],
        'title': 'Network Analysis Column Verification',
        'description': 'Verify comparison columns for network analysis. Pathway columns are auto-assigned.',
        'detect_samples': False,
        'multi_sheet': False,
        'can_calculate': ['Log2 Fold Change'],
    },
    'pathway_annotation': {
        'required': ['Feature ID'],
        'optional': ['Formula', 'HMDB ID', 'KEGG ID', 'PubChem CID', 'ChEBI ID', 
                    'CAS', 'SMILES', 'InChI', 'InChIKey', 'LipidMaps ID', 
                    'Class', 'Super Class', 'Sub Class', 'Metabolika Pathways'],
        'title': 'Pathway Annotation Column Assignment',
        'description': 'Verify columns for pathway annotation. Required: Feature ID (metabolite name). Optional: ID columns (HMDB ID, KEGG ID, PubChem CID, ChEBI ID, CAS, SMILES, InChI, InChIKey) and Metabolika Pathways (if from data cleaner).',
        'detect_samples': False,
        'multi_sheet': True,
        'can_calculate': [],
    },
    'lipid_pathway_annotation': {
        'required': ['Feature ID'],
        'recommended': ['Class_name'],
        'optional': ['Class', 'LipidMaps ID', 'HMDB ID', 'KEGG ID', 'PubChem CID', 'ChEBI ID',
                    'CAS', 'SMILES', 'InChI', 'InChIKey', 'Formula'],
        'title': 'Lipid Pathway Annotation Column Assignment',
        'description': 'Verify columns for lipid pathway annotation. Required: Feature ID (LipidID). Recommended: Class_name (for KEGG pathway lookup). Lipid mode performs ID annotation first, then proceeds to lipid pathway annotation.',
        'detect_samples': False,
        'multi_sheet': True,
        'can_calculate': [],
    },
    'id_annotation': {
        'required': ['Feature ID'],
        'optional': ['Formula'],
        'title': 'ID Annotation Column Assignment',
        'description': 'Verify columns for Custom ID Search. Required: Feature ID (metabolite name). Optional: Formula for additional validation.',
        'detect_samples': False,
        'multi_sheet': True,  # Enable multi-sheet support
        'column_mapping': {
            'Feature ID': 'Name',  # Feature ID is the metabolite Name for ID annotation
        }
    },
    'lipid_cleaning': {
        'required': ['LipidID'],  # Use LipidID directly, no abstraction
        'optional': ['Class', 'AdductIon', 'Area', 'Grade', 'Rej', 'ObsMz', 'ObsRt'],
        'title': 'Lipid Data Cleaning Column Assignment',
        'description': 'Verify columns for lipid data cleaning. Required: LipidID. Optional: Area[...] columns (sample data), Grade[...] columns (quality scores), Rej (rejection flag), ObsMz, ObsRt, Class, AdductIon.',
        'detect_samples': False,
        'multi_sheet': False,
        'detect_lipid_area': True,  # Enable Area[ ] column detection for lipids
    },
    'lipid_class': {
        'required': ['Class'],  # Class is the identifier for class sheets
        'optional': ['Class_name', 'Sample Column', 'Feature Column'],  # Same structure as statistics_lipid
        'title': 'Lipid Class Sheet - Column Verification',
        'description': 'Verify Class and sample columns for lipid class statistical analysis. Feature columns = class metadata (Class, Class_name, etc.). Sample columns = numeric intensity columns aggregated by class.',
        'detect_samples': True,  # Auto-detect numeric sample columns
        'multi_sheet': False,
        'can_calculate': [],  # No calculations - raw data only
        'detect_lipid_area': True,  # Enable Area[ ] column detection for class sheets
        'column_mapping': {
            'Class': 'Class',  # Ensure Class is preserved (equivalent to LipidID for class sheets)
        }
    },
    'lipid_id_annotation': {
        'required': ['LipidID'],
        'optional': ['Class', 'Class_name'],
        'title': 'Lipid ID Annotation Column Assignment',
        'description': 'Verify columns for Custom Lipid ID Search. Required: LipidID. Optional: Class. Recommended: Class_name.',
        'detect_samples': False,
        'multi_sheet': True,
    },
    'lipid_custom': {
        'required': ['Feature ID', 'P-Value', 'Log2 Fold Change'],
        'optional': ['HMDB ID', 'KEGG ID', 'InChIKey', 'ChEBI ID', 'Fold Change', 'Class', 'Super Class', 'LipidID', 'Class_name'],
        'title': 'Lipid Pathway Analysis Column Assignment',
        'description': 'Verify columns for lipid pathway annotation. Required: Feature ID (LipidID or Class), P-Value (p-value or p_adj), Log2 Fold Change. Sample columns will be auto-detected.',
        'detect_samples': True,
        'multi_sheet': False,
        'can_calculate': ['Log2 Fold Change'],
    },
}


# ===== COLUMN DETECTION & CALCULATION =====

class ColumnDetector:
    """Auto-detect and classify columns matching actual analysis module logic."""
    
    @staticmethod
    def normalize(text: str) -> str:
        """Normalize text for pattern matching.
        
        Removes underscores, dashes, spaces, dots, brackets, parentheses, and slashes
        to allow flexible matching of column names like 'Annot. DeltaMass [ppm]' or 'm/z'.
        """
        if not isinstance(text, str):
            text = str(text)
        # Remove common separators and special characters
        text = text.lower()
        for char in ['_', '-', ' ', '.', '[', ']', '(', ')', '/', '%']:
            text = text.replace(char, '')
        return text
    
    @staticmethod
    def is_numeric_convertible(series: pd.Series, min_ratio: float = 0.95) -> bool:
        """Check if a column is mostly numeric (matching statistics analysis logic)."""
        if pd.api.types.is_numeric_dtype(series):
            return True
        
        non_null = series.dropna()
        if len(non_null) == 0:
            return False
        
        convertible = 0
        for v in non_null.head(200):
            try:
                float(str(v).strip())
                convertible += 1
            except (ValueError, TypeError):
                pass
        
        ratio = convertible / max(1, len(non_null.head(200)))
        return ratio >= min_ratio
    
    @staticmethod
    def detect_column(col_name: str, col_type: str) -> bool:
        """Check if a column matches a given type using regex patterns.
        
        NOTE: All patterns are applied to NORMALIZED column names where underscores, 
        dashes, and spaces have been removed. E.g., 'adj_p' → 'adjp', 'log2_fc' → 'log2fc'
        """
        patterns = {
            # Feature ID: Prioritize exact 'name' match, exclude columns ending with 'id' 
            # (like metaboliteid, compoundid) to avoid false positives
            'Feature ID': [r'^name$', r'lipidid', r'^metabolite$', r'^compound$', r'^feature$', r'metabolitename', r'compoundname', r'featurename'],
            'HMDB ID': [r'hmdbiid', r'hmdbid', r'hmdb'],  # hmdb_id → hmdbiid, hmdb_id → hmdbid
            'KEGG ID': [r'keggid', r'kegg'],  # kegg_id → keggid
            'ChEBI ID': [r'chebiid', r'chebi'],  # chebi_id → chebiid
            'PubChem CID': [r'pubchemcid', r'pubchem'],  # pubchem_cid → pubchemcid
            'CAS': [r'^cas$', r'casnumber', r'casrn'],  # cas, cas_number → casnumber
            'SMILES': [r'^smiles$', r'smilesstring'],  # smiles
            'InChI': [r'^inchi$'],  # inchi (exact match to avoid matching inchikey)
            'InChIKey': [r'inchikey'],  # inchi_key → inchikey
            'LipidMaps ID': [r'lipidmapsid', r'lipidmaps', r'lmid'],  # lipidmaps_id → lipidmapsid
            'Sub Class': [r'subclass'],  # sub_class → subclass
            'P-Value': [r'adjp', r'padj', r'fdr', r'qvalue', r'pvalue', r'^p$'],  # adj_p → adjp, p_adj → padj
            'Fold Change': [r'fc$', r'^fc$'],  # Ends with fc (e.g., pdvscontrolfc) or exactly fc
            'Log2 Fold Change': [r'log2fc', r'logfc', r'logfoldch'],  # log2_fc → log2fc, log_fold_ch → logfoldch
            '-Log10 P-Value': [r'neglog', r'negadj', r'log10p', r'log10adj'],  # neg_log → neglog, neg_adj_p → negadjp
            'Class': [r'^class$', r'lipidclass$', r'metaboliteclass$'],  # Exactly 'class', 'lipidclass', or 'metaboliteclass'
            'Class_name': [r'^classname$', r'^class_name$', r'lipidclassname'],  # class_name, classname, lipid_class_name
            'Super Class': [r'superclass'],  # super_class → superclass
            'Formula': [r'^formula$', r'molecularformula', r'molformula'],  # molecular_formula → molecularformula
            'Metabolika Pathways': [r'^metabolikapathway', r'metabolikapathway$'],  # Must match exactly or contain word
            'BioCyc Pathways': [r'^biocycpathway', r'biocycpathway$'],  # Must match exactly or contain word
            'Group Area': [r'grouparea', r'group_area'],  # group_area → grouparea
            'Reference Ion': [r'referenceion', r'refion', r'reference'],  # Reference Ion, Ref Ion
            'Area (Max.)': [r'areamax', r'area\(max', r'maxarea'],  # Area (Max.), AreaMax, Area(Max)
            'RT [min]': [r'rt\[min', r'rtmin', r'^rt$'],  # RT [min], RT[min], RTmin, RT
            'Polarity': [r'polarity', r'^pol$'],  # Polarity, Pol
            'MS2': [r'^ms2$', r'msms'],  # MS2, MSMS (exact match for ms2 to avoid matching other columns)
            'MS2 Purity [%]': [r'ms2purity', r'purity'],  # MS2 Purity [%], MS2Purity, Purity
            'LipidID': [r'lipidid', r'^lipid$'],  # LipidID or Lipid
            'Area': [r'^area\['],  # Lipid area columns: Must start exactly with Area[ (no prefix like OrgArea or MeanArea)
            'Grade': [r'^grade\['],  # Lipid grade columns: Quality scores for filtering
            'Rej': [r'^rej$', r'reject', r'rejected'],  # Rejection flag column
            'AdductIon': [r'adduct', r'^ion$'],  # AdductIon, Adduct, or Ion
            'ObsMz': [r'^calcmz', r'calcmz'],  # CalcMz column
            'ObsRt': [r'^basert', r'basert'],  # BaseRt column
            # Data cleaning specific columns
            'Annot. DeltaMass [ppm]': [r'annotdeltamass', r'deltamass.*ppm', r'deltappm', r'ppm'],  # Annot. DeltaMass [ppm]
            'Calc. MW': [r'calcmw', r'calculatedmw', r'molecularweight', r'^mw$'],  # Calc. MW, Calculated MW, MW
            'm/z': [r'^mz$', r'm/z', r'masstocharge', r'precursormz'],  # m/z, mz, mass to charge
            #'Delta(PPM)': [r'ppm', r'delta.*ppm', r'deltappm'],  # PPM, Delta(PPM), Delta_PPM
        }
        
        if col_type not in patterns:
            return False
        
        normalized = ColumnDetector.normalize(col_name)
        for pattern in patterns[col_type]:
            if re.search(pattern, normalized):
                return True
        return False

    @staticmethod
    def select_best_pvalue_column(candidates: List[str]) -> Optional[str]:
        """Choose the most informative p-value column from candidates using priority rules."""
        if not candidates:
            return None
        
        # CRITICAL: Filter out -log10 transformed p-value columns FIRST
        # These are NOT actual p-values and will break Fisher ORA analysis
        valid_candidates = []
        for col in candidates:
            normalized = re.sub(r'[^a-z0-9]', '', str(col).lower())
            # Exclude any column with log transformation indicators
            if any(keyword in normalized for keyword in ['neglog10', 'log10', 'log10p', 'neglog']):
                continue
            valid_candidates.append(col)
        
        if not valid_candidates:
            return None

        def _priority(col_name: str) -> tuple:
            normalized = re.sub(r'[^a-z0-9]', '', str(col_name).lower())
            if any(keyword in normalized for keyword in ['adjp', 'padj', 'fdr', 'qvalue', 'adjustedp']):
                return (0, len(normalized))  # Highest priority: adjusted p-values
            if re.search(r'pvalue$', normalized) or normalized == 'p':
                return (1, len(normalized))  # Raw p-values next
            if 'fpvalue' in normalized or normalized.startswith('fstat'):
                return (2, len(normalized))  # F-statistic p-values last
            return (3, len(normalized))  # Unknown patterns fall back here

        return min(valid_candidates, key=_priority)
    
    
    @staticmethod
    def detect_sample_columns(df: pd.DataFrame, min_numeric_ratio: float = 0.95, exclude_cols: Optional[List[str]] = None) -> List[str]:
        """
        Detect sample/intensity columns using statistics module logic.
        
        Returns numeric columns that are NOT in the known feature set or metadata.
        Excludes group labels, metadata columns, pairwise stats columns, and any verified statistics columns.
        
        Args:
            df: DataFrame to analyze
            min_numeric_ratio: Threshold for numeric conversion
            exclude_cols: Additional columns to exclude (e.g., verified statistics columns like 'abd', 'cdf', 'gft')
        """
        sample_cols: List[str] = []
        
        # Columns to exclude
        exclude_set = set(exclude_cols) if exclude_cols else set()
        
        # Patterns for columns to exclude (metadata, group labels, pairwise stats)
        exclude_patterns = [
            r'^(PD|Control|TBI|Disease|Treatment|Group|Sample_Group)$',  # Group/condition labels
            r'_vs_',  # Pairwise comparison columns
            r'adj_p|p_value|pvalue|log2fc|logfc|fc\d+$',  # Stats columns
            r'neg_log10|neglog10|-log10',  # Transformed stats
        ]
        
        for col in df.columns:
            if is_statistics_metadata_col(col):
                continue
            col_lower = col.lower()
            
            # Skip if in explicit exclude list
            if col in exclude_set:
                continue
            
            # Skip known feature/metadata columns
            if any(col.lower() == c.lower() for c in COMMON_ID_COLUMNS) or any(col.lower() == c.lower() for c in FEATURE_COLUMNS_CANONICAL):
                continue
            
            # Skip if matches exclude patterns
            if any(re.search(pattern, col_lower) for pattern in exclude_patterns):
                continue
            
            # Check if column is numeric or numeric-convertible
            if ColumnDetector.is_numeric_convertible(df[col], min_numeric_ratio):
                sample_cols.append(col)
        
        return sample_cols
    
    @staticmethod
    def detect_feature_and_sample_columns(df: pd.DataFrame, min_numeric_ratio: float = 0.95) -> Tuple[List[str], List[str]]:
        """
        Detect feature and sample columns from dataframe.
        
        Returns:
            Tuple of (feature_columns_present, sample_columns)
            
        Matches the logic from metabolite_statistics_analysis.detect_feature_and_sample_columns
        """
        feature_present = [c for c in FEATURE_COLUMNS_CANONICAL if c in df.columns]
        sample_cols = ColumnDetector.detect_sample_columns(df, min_numeric_ratio)
        
        return feature_present, sample_cols
    
    @staticmethod
    def auto_detect_all(df: pd.DataFrame) -> Dict[str, Optional[str]]:
        """Auto-detect all known column types in dataframe, including numeric features like FC1, FC2, etc."""
        detected = {}
        numeric_features = {}  # Collect all FC* columns
        qc_columns = []  # Collect QC sample columns
        is_feature = None  # Internal Standard feature
        
        for col in df.columns:
            col_lower = col.lower()
            
            # Check if it's a numeric feature column (FC1, FC2, etc. - used in volcano plots)
            if re.match(r'^fc\d+$', col_lower) and ColumnDetector.is_numeric_convertible(df[col]):
                numeric_features[col] = col
                continue
            
            # Check for QC sample columns (contains 'qc', 'quality', 'pool' in name)
            qc_keywords = ['qc', 'quality', 'pool', 'blank', 'pooled']
            if any(keyword in col_lower for keyword in qc_keywords):
                if ColumnDetector.is_numeric_convertible(df[col]):
                    qc_columns.append(col)
                    continue
            
            # Standard column type detection
            for col_type in COLUMN_TYPES.keys():
                if ColumnDetector.detect_column(col, col_type):
                    if col_type not in detected:
                        detected[col_type] = col
                    break
        
        # Check for Internal Standard feature (row-based, not column-based)
        # Look for features with 'internal', 'IS', 'standard' in the name column
        if 'Feature ID' in detected:
            feature_col = detected['Feature ID']
            if feature_col in df.columns:
                for idx, val in enumerate(df[feature_col]):
                    val_str = str(val).strip().lower()
                    # Check for various IS naming patterns
                    is_patterns = [
                        r'\bis\b',  # standalone 'IS'
                        'internal',  # contains 'internal'
                        'standard',  # contains 'standard'
                        'istd',  # ISTD abbreviation
                    ]
                    if any(re.search(pattern, val_str) for pattern in is_patterns):
                        is_feature = str(df[feature_col].iloc[idx])
                        break
        
        # If we found numeric feature columns, store them for UI display
        if numeric_features:
            detected['_numeric_features'] = numeric_features  # Special key for UI to display FC1, FC2, etc.
        
        # Store detected QC columns
        if qc_columns:
            detected['_qc_columns'] = qc_columns  # List of QC sample columns
        
        # Store detected Internal Standard feature
        if is_feature:
            detected['Internal Standard'] = is_feature
        
        return detected
    
    @staticmethod
    def detect_group_area_columns_by_pattern(df: pd.DataFrame, pattern: str) -> List[str]:
        """Detect all Group Area columns matching given text pattern(s).
        
        Args:
            df: DataFrame to search
            pattern: Text pattern(s) to match (case-insensitive, supports comma-separated patterns)
            
        Returns:
            List of column names matching any of the patterns
        """
        if not pattern or not pattern.strip():
            return []
        
        # EXCLUDED columns - these are NEVER sample/Group Area columns
        # These are metadata columns that should not be assigned as Group Area
        EXCLUDED_COLUMN_PREFIXES = [
            'calc', 'formula', 'name', 'annotation', 'annot', 'm/z', 'mz', 'mass',
            'rt', 'retention', 'reference', 'polarity', 'ms2', 'ms1', 'adduct',
            'pathway', 'metabolika', 'kegg', 'hmdb', 'pubchem', 'inchi', 'smiles',
            'feature', 'compound', 'id', 'class', 'category', 'type', 'ion',
            'delta', 'ppm', 'score', 'match', 'library', 'fragmentation', 'area (max',
            'area(max', 'area max', 'checked', 'fill', 'molecular', 'neutral',
            'charge', 'isotope', 'peak', 'intensity', 'height'
        ]
        
        EXCLUDED_EXACT_COLUMNS = [
            'name', 'formula', 'm/z', 'rt', 'rt [min]', 'rt[min]', 'polarity', 'ms2',
            'reference ion', 'area (max.)', 'area(max)', 'area max', 'checked',
            'fill %', 'molecular weight', 'mw', 'calc. mw', 'calc mw', 'calculated mw'
        ]
        
        # Support comma-separated patterns
        patterns = [p.strip().lower() for p in pattern.split(',') if p.strip()]
        matching_cols = []
        
        for col in df.columns:
            col_lower = str(col).lower().strip()
            
            # Skip if column is in excluded exact matches
            if col_lower in EXCLUDED_EXACT_COLUMNS:
                continue
            
            # Skip if column starts with any excluded prefix
            is_excluded = False
            for excl in EXCLUDED_COLUMN_PREFIXES:
                if col_lower.startswith(excl):
                    is_excluded = True
                    break
            
            if is_excluded:
                continue
            
            # Check if column matches any of the user's patterns
            for p in patterns:
                # For very short patterns (1-2 chars), require exact word boundary match
                # to avoid matching "C" to "Calc. MW"
                if len(p) <= 2:
                    # Check for pattern at start followed by non-letter (word boundary)
                    # e.g., "C_sample" or "C1" or "C " matches, but "Calc" does not
                    import re
                    # Pattern should be at start and followed by non-letter or end of string
                    if re.match(rf'^{re.escape(p)}(?:[^a-zA-Z]|$)', col_lower):
                        # Verify it's numeric (sample data)
                        if ColumnDetector.is_numeric_convertible(df[col]):
                            matching_cols.append(col)
                            break
                else:
                    # For longer patterns, use startswith
                    if col_lower.startswith(p):
                        # Verify it's numeric (sample data)
                        if ColumnDetector.is_numeric_convertible(df[col]):
                            matching_cols.append(col)
                            break  # Only add once even if multiple patterns match
        
        return matching_cols
    
    @staticmethod
    def count_group_area_columns(df: pd.DataFrame, pattern: str = None) -> int:
        """Count the number of Group Area columns.
        
        Args:
            df: DataFrame to analyze
            pattern: Optional text pattern to match
            
        Returns:
            Number of matching Group Area columns
        """
        if pattern:
            return len(ColumnDetector.detect_group_area_columns_by_pattern(df, pattern))
        
        # Default: look for 'Group Area:' or 'group area'
        count = 0
        for col in df.columns:
            if 'group area' in col.lower() and ColumnDetector.is_numeric_convertible(df[col]):
                count += 1
        
        return count


class ColumnCalculator:
    """Calculate missing columns from existing data."""
    
    @staticmethod
    def calculate_log2fc(df: pd.DataFrame, fc_col: str) -> pd.Series:
        """Calculate log2 fold change from fold change."""
        fc_values = pd.to_numeric(df[fc_col], errors='coerce')
        result = np.log2(fc_values)
        # Replace inf values with NaN
        result = pd.Series(result, index=df.index)
        result = result.replace([np.inf, -np.inf], np.nan)
        return result
    
    @staticmethod
    def calculate_neglog10_pvalue(df: pd.DataFrame, pvalue_col: str) -> pd.Series:
        """Calculate -log10 p-value from p-value."""
        pvalue = pd.to_numeric(df[pvalue_col], errors='coerce')
        # Handle p-values of 0
        pvalue = pvalue.replace(0, 1e-300)
        result = -np.log10(pvalue)
        return pd.Series(result, index=df.index)
    
    @staticmethod
    def calculate_missing_columns(
        df: pd.DataFrame,
        assignments: Dict[str, Optional[str]],
        can_calculate: Optional[List[str]] = None
    ) -> Tuple[pd.DataFrame, Dict[str, Optional[str]]]:
        """
        Calculate missing columns if possible.
        
        Returns:
            Tuple of (modified_df, new_assignments)
        """
        if can_calculate is None:
            can_calculate = []
        
        new_assignments = assignments.copy()
        modified_df = df.copy()
        
        # Calculate log2FC from FC if needed
        if 'Log2 Fold Change' in can_calculate and not assignments.get('Log2 Fold Change'):
            fc_col = assignments.get('Fold Change')
            if fc_col and fc_col in df.columns:
                try:
                    new_col_name = 'log2FC_calculated'
                    modified_df[new_col_name] = ColumnCalculator.calculate_log2fc(modified_df, fc_col)
                    new_assignments['Log2 Fold Change'] = new_col_name
                    logger.info(f"Calculated log2FC from {fc_col}")
                except Exception as e:
                    logger.warning(f"Failed to calculate log2FC: {e}")
        
        # Calculate -log10 p-value if needed
        if '-Log10 P-Value' in can_calculate and not assignments.get('-Log10 P-Value'):
            pval_col = assignments.get('P-Value')
            if pval_col and pval_col in df.columns:
                try:
                    new_col_name = 'neglog10_pvalue_calculated'
                    modified_df[new_col_name] = ColumnCalculator.calculate_neglog10_pvalue(modified_df, pval_col)
                    new_assignments['-Log10 P-Value'] = new_col_name
                    logger.info(f"Calculated -log10 p-value from {pval_col}")
                except Exception as e:
                    logger.warning(f"Failed to calculate -log10 p-value: {e}")
        
        return modified_df, new_assignments


# ===== MULTI-COMPARISON COLUMN MATCHING =====

class ComparisonColumnMatcher:
    """
    Identify which columns are required/optional for each pairwise comparison.
    
    When users have multiple comparisons (e.g., Control vs PD, Control vs TBI, PD vs TBI),
    they may have different columns for each comparison. This matcher helps identify:
    - Which columns contain FC/pvalue/log2fc for each comparison
    - Which comparisons have complete data vs missing data
    - Per-comparison column requirements
    """
    
    @staticmethod
    def find_comparison_columns(
        df: pd.DataFrame,
        groups: Optional[List[str]] = None
    ) -> Dict[str, Dict[str, Optional[str]]]:
        """
        Find all pairwise comparison columns and their types.
        
        Returns:
            Dict mapping comparison name (e.g., 'Control_vs_PD') to column mappings:
            {
                'Control_vs_PD': {
                    'log2fc': 'Control_vs_PD_log2FC',
                    'pvalue': 'Control_vs_PD_pvalue',
                    '-log10p': 'Control_vs_PD_neg_log10_adj_p',
                    ...
                },
                'Control_vs_TBI': {...},
                ...
            }
        """
        comparisons: Dict[str, Dict[str, Optional[str]]] = {}
        
        # Pattern to match comparison columns: Group1_vs_Group2_*
        comparison_pattern = r'([A-Za-z0-9_]+)_vs_([A-Za-z0-9_]+)_(.+)$'
        
        for col in df.columns:
            match = re.match(comparison_pattern, col)
            if not match:
                continue
            
            group1, group2, stat_type = match.groups()
            comp_name = f"{group1}_vs_{group2}"
            
            if comp_name not in comparisons:
                comparisons[comp_name] = {}
            
            # Classify stat_type
            stat_normalized = ColumnDetector.normalize(stat_type)
            
            if re.search(r'log2fc|logfc', stat_normalized):
                comparisons[comp_name]['log2fc'] = col
            elif re.search(r'^fc$|fold', stat_normalized):
                comparisons[comp_name]['fc'] = col
            elif re.search(r'adjp|padj|fdr', stat_normalized):
                comparisons[comp_name]['pvalue_adj'] = col
            elif re.search(r'^p$|pvalue', stat_normalized):
                comparisons[comp_name]['pvalue'] = col
            elif re.search(r'neglog|negadj|log10p', stat_normalized):
                comparisons[comp_name]['-log10p'] = col
        
        return comparisons
    
    @staticmethod
    def validate_comparison_completeness(
        comparison_columns: Dict[str, Dict[str, Optional[str]]]
    ) -> Dict[str, Dict[str, bool]]:
        """
        Check if each comparison has required columns.
        
        Returns:
            Dict mapping comparison name to completeness dict:
            {
                'Control_vs_PD': {
                    'has_fc': True/False,          # has FC or log2FC
                    'has_pvalue': True/False,      # has raw or adjusted p-value
                    'complete': True/False,        # has both FC and p-value
                    'priority': 'high'/'medium'/'low'
                },
                ...
            }
        """
        completeness: Dict[str, Dict[str, bool]] = {}
        
        for comp_name, cols in comparison_columns.items():
            has_any_fc = 'log2fc' in cols or 'fc' in cols
            has_any_pvalue = 'pvalue_adj' in cols or 'pvalue' in cols or '-log10p' in cols
            
            completeness[comp_name] = {
                'has_fc': has_any_fc,
                'has_pvalue': has_any_pvalue,
                'complete': has_any_fc and has_any_pvalue,
                'priority': 'high' if (has_any_fc and has_any_pvalue) else 'medium' if (has_any_fc or has_any_pvalue) else 'low',
            }
        
        return completeness
    
    @staticmethod
    def get_required_columns_per_comparison(
        comparison_columns: Dict[str, Dict[str, Optional[str]]]
    ) -> Dict[str, List[str]]:
        """
        Get list of required/missing columns for each comparison.
        
        Returns:
            Dict mapping comparison to list of missing items:
            {
                'Control_vs_PD': [],           # Complete
                'Control_vs_TBI': ['log2FC'],  # Missing log2FC
                'PD_vs_TBI': ['p-value'],      # Missing p-value
            }
        """
        missing: Dict[str, List[str]] = {}
        
        for comp_name, cols in comparison_columns.items():
            missing_list: List[str] = []
            
            if 'log2fc' not in cols and 'fc' not in cols:
                missing_list.append('Fold Change (FC or log2FC)')
            
            if 'pvalue_adj' not in cols and 'pvalue' not in cols and '-log10p' not in cols:
                missing_list.append('P-Value')
            
            missing[comp_name] = missing_list
        
        return missing


# ===== MAIN ASSIGNMENT DIALOG =====

class ColumnAssignmentDialog(tk.Toplevel):
    """
    Unified column assignment dialog for all tabs with multi-sheet and sample detection support.
    
    Features:
    - Shows all columns and allows user to assign column types via dropdown
    - Multi-sheet support for ID Annotation (lets users select which sheet to use)
    - Auto-detection of sample columns for Statistics/Visualization/Pathway
    - Auto-calculation of missing columns (log2FC, -log10_pvalue)
    - Pre-fills with detected assignments
    """
    
    def __init__(
        self,
        parent: tk.Tk,
        df: pd.DataFrame,
        tab_type: str = 'general',
        auto_calculate: bool = True,
        existing_assignments: Optional[Dict[str, Optional[str]]] = None,
        excel_file_path: Optional[str] = None,  # For multi-sheet support
        group_area_pattern: Optional[str] = None,  # Custom pattern for Group Area detection
        detected_sample_cols: Optional[List[str]] = None,  # Pre-detected sample columns
        allow_skip: bool = False,  # When True, show a Skip button and allow caller to skip this sheet
    ):
        """
        Initialize column assignment dialog.
        
        Args:
            parent: Parent window
            df: DataFrame to assign columns for
            tab_type: Type of tab ('data_cleaning', 'statistics', 'visualization', 'pathway', 'id_annotation')
            auto_calculate: Whether to calculate missing columns
            existing_assignments: Previously confirmed assignments to pre-fill
            excel_file_path: Path to Excel file (for multi-sheet support in ID annotation)
            group_area_pattern: Optional custom pattern for Group Area column detection
            detected_sample_cols: Optional pre-detected sample columns to use instead of auto-detecting
        """
        super().__init__(parent)
        
        self.df = df
        self.tab_type = tab_type
        self.auto_calculate = auto_calculate
        self.result = None
        self.modified_df = None
        self.excel_file_path = excel_file_path
        self.group_area_pattern = group_area_pattern
        self.detected_sample_cols_input = detected_sample_cols  # Store pre-detected columns
        self.allow_skip = allow_skip
        self.existing_assignments = existing_assignments  # Store for pre-filling dialog
        self.group_area_auto_pattern = ""  # Will be set during auto-detection
        self.available_sheets = []
        self.selected_sheet = None
        self.detection_complete = False
        
        # Get requirements for this tab
        self.tab_config = TAB_REQUIREMENTS.get(tab_type, TAB_REQUIREMENTS['data_cleaning'])
        
        # Initialize these as empty - will be filled during detection
        self.detected = {}
        self.sample_cols = []
        
        # Column type variables (will store selected column for each type)
        self.assignments: Dict[str, tk.StringVar] = {}
        
        # Setup UI immediately (with loading indicator)
        self.title(self.tab_config['title'])
        self.geometry('900x750')  # Wider and taller for better view
        self.transient(parent)
        self.grab_set()
        
        self._setup_ui_with_loading()
        
        # Start detection in background after window is displayed
        self.after(100, self._perform_detection)
    
    def _perform_detection(self):
        """Perform column detection in a background thread after UI is shown."""
        # Start detection in a thread to avoid blocking UI
        detection_thread = threading.Thread(target=self._detection_worker, daemon=True)
        detection_thread.start()
    
    def _detection_worker(self):
        """Worker thread for column detection."""
        try:
            # Auto-detect columns using the real logic
            self.detected = ColumnDetector.auto_detect_all(self.df)
            
            # Auto-detect Group Area columns - those already containing "Group Area:" prefix
            if self.tab_config.get('detect_group_area', False):
                # Find columns that already have "Group Area:" prefix
                auto_group_area_cols = [col for col in self.df.columns 
                                       if isinstance(col, str) and 'Group Area:' in col]
                
                if auto_group_area_cols:
                    # Store them for auto-assignment
                    self.auto_assign_group_area_cols = auto_group_area_cols
                    self.group_area_auto_pattern = 'Group Area:'
                    self.sample_cols = auto_group_area_cols
                    logger.info(f"Auto-detected {len(auto_group_area_cols)} columns with 'Group Area:' prefix")
            
            # Auto-detect Area columns for lipids - those starting EXACTLY with "Area[" (no prefix)
            if self.tab_config.get('detect_lipid_area', False):
                # Find columns that start EXACTLY with "Area[" pattern (excludes OrgArea[, MeanArea[, etc.)
                auto_lipid_area_cols = [col for col in self.df.columns 
                                       if isinstance(col, str) and col.startswith('Area[') and col.endswith(']')]
                
                # Also detect Grade columns for lipid quality filtering
                auto_lipid_grade_cols = [col for col in self.df.columns 
                                        if isinstance(col, str) and col.startswith('Grade[') and col.endswith(']')]
                
                if auto_lipid_area_cols:
                    # Store them for auto-assignment
                    self.auto_assign_lipid_area_cols = auto_lipid_area_cols
                    self.auto_assign_lipid_grade_cols = auto_lipid_grade_cols
                    self.lipid_area_auto_pattern = 'Area['
                    self.sample_cols = auto_lipid_area_cols
                    logger.info(f"Auto-detected {len(auto_lipid_area_cols)} columns starting with 'Area[' pattern")
                    logger.info(f"Auto-detected {len(auto_lipid_grade_cols)} columns starting with 'Grade[' pattern")
            
            # Detect sample columns if needed
            if self.tab_config.get('detect_samples', False):
                # Use pre-detected sample columns if provided, otherwise auto-detect
                if self.detected_sample_cols_input is not None:
                    self.sample_cols = self.detected_sample_cols_input
                    logger.info(f"✅ Using pre-detected {len(self.sample_cols)} sample columns: {self.sample_cols[:5]}...")
                    print(f"DEBUG: Pre-detected sample columns count: {len(self.sample_cols)}")
                    print(f"DEBUG: First 5 sample columns: {self.sample_cols[:5]}")
                else:
                    _, self.sample_cols = ColumnDetector.detect_feature_and_sample_columns(self.df)
                    logger.info(f"⚠️ Auto-detected {len(self.sample_cols)} sample columns (no pre-detected input)")
                    print(f"DEBUG: Auto-detected sample columns count: {len(self.sample_cols)}")
                    print(f"DEBUG: First 5 sample columns: {self.sample_cols[:5]}")
            
            # Load available sheets if multi-sheet support enabled
            if self.tab_config.get('multi_sheet', False) and self.excel_file_path:
                self._load_available_sheets()
            
            # Mark detection complete and update UI in main thread
            self.detection_complete = True
            self.after(0, self._finish_loading)
            
        except Exception as e:
            logger.error(f"Error during column detection: {e}")
            self.detection_complete = True
            self.after(0, lambda: messagebox.showerror("Detection Error", f"Error detecting columns: {str(e)}"))
    
    def _finish_loading(self):
        """Finish loading process and populate columns (called in main thread)."""
        try:
            # Hide loading and populate columns
            self.loading_frame.pack_forget()
            self._populate_column_assignments()
        except Exception as e:
            logger.error(f"Error finishing loading: {e}")
            messagebox.showerror("Error", f"Error loading columns: {str(e)}")
    
    def _load_available_sheets(self):
        """Load sheet names from Excel file for ID annotation multi-sheet support."""
        try:
            if self.excel_file_path and self.excel_file_path.endswith(('.xls', '.xlsx')):
                xl_file = pd.ExcelFile(self.excel_file_path)
                self.available_sheets = [str(s) for s in xl_file.sheet_names]  # Ensure all are strings
                # Pre-select "Combined" sheet if available, otherwise first sheet
                if 'Combined' in self.available_sheets:
                    self.selected_sheet = 'Combined'
                else:
                    self.selected_sheet = self.available_sheets[0] if self.available_sheets else None
        except Exception as e:
            logger.warning(f"Could not load sheets from Excel file: {e}")
            self.available_sheets = []
    
    def _setup_ui_with_loading(self):
        """Build the dialog UI with a loading indicator."""
        # Main frame
        main_frame = ttk.Frame(self)
        main_frame.pack(fill='both', expand=True, padx=6, pady=4)
        
        # Compact title section
        title_frame = ttk.Frame(main_frame)
        title_frame.pack(fill='x', pady=(0, 3))
        
        title_label = ttk.Label(
            title_frame,
            text=self.tab_config['title'],
            font=('Arial', 10, 'bold')
        )
        title_label.pack(anchor='w')
        
        desc_label = ttk.Label(
            title_frame,
            text=self.tab_config['description'],
            font=('Arial', 8),
            foreground='gray'
        )
        desc_label.pack(anchor='w')
        
        # Required and Essential columns info (compact)
        info_frame = ttk.Frame(main_frame)
        info_frame.pack(fill='x', pady=(0, 3))
        
        required_text = ", ".join(self.tab_config['required'])
        ttk.Label(
            info_frame,
            text=f"✓ Required: {required_text}",
            font=('Arial', 8, 'bold'),
            foreground='#27ae60'
        ).pack(anchor='w')
        
        # Show recommended columns if defined
        if 'recommended' in self.tab_config:
            recommended_text = ", ".join(self.tab_config['recommended'])
            ttk.Label(
                info_frame,
                text=f"★ Recommended: {recommended_text}",
                font=('Arial', 8, 'bold'),
                foreground='#3498db'
            ).pack(anchor='w')
        
        # Show essential columns if defined
        if 'essential' in self.tab_config:
            essential_text = ", ".join(self.tab_config['essential'])
            ttk.Label(
                info_frame,
                text=f"⚠ Essential: {essential_text} (confirm assignments before proceeding)",
                font=('Arial', 8, 'bold'),
                foreground='#f39c12'
            ).pack(anchor='w')
        
        # Save reference to main_frame for later use
        self.main_frame = main_frame
        
        # Loading frame (will be replaced with content)
        self.loading_frame = ttk.Frame(main_frame)
        self.loading_frame.pack(fill='both', expand=True, pady=10)
        
        ttk.Label(
            self.loading_frame,
            text='🔄 Analyzing columns, please wait...',
            font=('Arial', 11, 'bold'),
            foreground='#3498db'
        ).pack(pady=20)
        
        ttk.Label(
            self.loading_frame,
            text='Detecting column types (Name, Formula, Pathways, etc.)',
            font=('Arial', 9),
            foreground='#7f8c8d'
        ).pack(pady=5)
        
        ttk.Label(
            self.loading_frame,
            text='This may take a few seconds for large files...',
            font=('Arial', 9, 'italic'),
            foreground='#95a5a6'
        ).pack(pady=5)
        
        # Create progress bar
        progress = ttk.Progressbar(self.loading_frame, mode='indeterminate', length=300)
        progress.pack(pady=10)
        progress.start()
        
        self.progress_bar = progress
    
    def _populate_column_assignments(self):
        """Populate the column assignment UI after detection is complete."""
        try:
            # Clear loading frame
            if hasattr(self, 'loading_frame'):
                for widget in self.loading_frame.winfo_children():
                    widget.destroy()
                self.loading_frame.pack_forget()
            
            # Setup the actual UI content
            self._setup_ui()
            
            # Update window to display new content
            self.update_idletasks()
            
        except Exception as e:
            logger.error(f"Error populating columns: {e}")
            messagebox.showerror("Error", f"Error populating columns: {str(e)}")
    
    def _setup_ui(self):
        """Build the dialog UI."""
        # Reuse the existing main frame created during loading so the
        # top info area remains fixed and only the lower section expands
        main_frame = getattr(self, 'main_frame', None)
        if main_frame is None or not isinstance(main_frame, ttk.Frame):
            # Fallback: create if missing (shouldn't normally happen)
            main_frame = ttk.Frame(self)
            main_frame.pack(fill='both', expand=True, padx=6, pady=4)
        
        # Multi-sheet selector for ID annotation
        if self.tab_config.get('multi_sheet', False) and self.available_sheets:
            sheet_frame = ttk.LabelFrame(main_frame, text='Select Sheet for ID Search', padding=5)
            sheet_frame.pack(fill='x', pady=(0, 10))
            
            sheet_value = str(self.selected_sheet) if self.selected_sheet else ''
            self.sheet_var = tk.StringVar(value=sheet_value)
            sheet_combo = ttk.Combobox(
                sheet_frame,
                textvariable=self.sheet_var,
                values=self.available_sheets,
                state='readonly',
                width=50
            )
            sheet_combo.pack(side='left', padx=5)
            
            ttk.Label(
                sheet_frame,
                text='(Typically "Combined" contains all metabolites for ID search)',
                font=('Arial', 8),
                foreground='gray'
            ).pack(side='left', padx=5)
        
        # Sample columns info (if detected) - more compact
        if self.sample_cols:
            sample_frame = ttk.LabelFrame(main_frame, text='Detected Sample Columns', padding=2)
            sample_frame.pack(fill='x', pady=(0, 3))
            
            sample_text = ', '.join(self.sample_cols[:5])
            if len(self.sample_cols) > 5:
                sample_text += f', ... +{len(self.sample_cols) - 5} more'
            
            ttk.Label(
                sample_frame,
                text=f"{len(self.sample_cols)} columns: {sample_text}",
                font=('Arial', 8),
                foreground='#3498db'
            ).pack(anchor='w', padx=3, pady=2)
        
        # Feature columns info (for statistics tabs) - show how many non-sample columns detected
        if self.tab_type in ['statistics_metabolite', 'statistics_lipid']:
            feature_cols = [col for col in self.df.columns if col not in self.sample_cols]
            if feature_cols:
                feature_frame = ttk.LabelFrame(main_frame, text='Detected Feature/Metadata Columns', padding=2)
                feature_frame.pack(fill='x', pady=(0, 3))
                
                feature_text = ', '.join(feature_cols[:5])
                if len(feature_cols) > 5:
                    feature_text += f', ... +{len(feature_cols) - 5} more'
                
                ttk.Label(
                    feature_frame,
                    text=f"{len(feature_cols)} columns: {feature_text}",
                    font=('Arial', 8),
                    foreground='#27ae60'
                ).pack(anchor='w', padx=3, pady=2)
        
        # Group Area pattern input (for data cleaning tab)
        if self.tab_config.get('detect_group_area', False):
            pattern_frame = ttk.LabelFrame(main_frame, text='🔍 Group Area Pattern (optional)', padding=5)
            pattern_frame.pack(fill='x', pady=(0, 10))
            
            ttk.Label(
                pattern_frame,
                text='Enter pattern(s) to search for Group Area columns (e.g., "Group Area", "LF_pos, Peak", etc.). Use commas to separate multiple patterns.',
                font=('Arial', 8),
                foreground='#7f8c8d'
            ).pack(anchor='w', padx=5, pady=(0, 3))
            
            # Pre-fill with auto-detected pattern if available
            auto_pattern = getattr(self, 'group_area_auto_pattern', 'Group Area')
            self.group_area_pattern_var = tk.StringVar(value=auto_pattern)
            pattern_entry = ttk.Entry(
                pattern_frame,
                textvariable=self.group_area_pattern_var,
                width=50
            )
            pattern_entry.pack(fill='x', padx=5, pady=(0, 3))
            
            def apply_pattern(*args):
                """Apply pattern and show matching Group Area columns"""
                pattern = self.group_area_pattern_var.get().strip()
                if pattern:
                    detected = ColumnDetector.detect_group_area_columns_by_pattern(self.df, pattern)
                    if detected:
                        status_text = f"Found {len(detected)} Group Area column(s): {', '.join(detected[:3])}"
                        if len(detected) > 3:
                            status_text += f" ... and {len(detected) - 3} more"
                        pattern_status.config(text=status_text, foreground='#27ae60')
                        assign_button.config(state='normal')
                    else:
                        pattern_status.config(text="No columns matched this pattern", foreground='#e74c3c')
                        assign_button.config(state='disabled')
                else:
                    pattern_status.config(text="", foreground='#7f8c8d')
                    assign_button.config(state='disabled')
            
            self.group_area_pattern_var.trace('w', apply_pattern)
            
            pattern_status = ttk.Label(
                pattern_frame,
                text="",
                font=('Arial', 8),
                foreground='#7f8c8d'
            )
            pattern_status.pack(anchor='w', padx=5)
            
            # Button frame for Assign Columns button
            button_frame = ttk.Frame(pattern_frame)
            button_frame.pack(fill='x', padx=5, pady=(5, 0))
            
            def assign_all_group_areas():
                """Assign and rename all detected Group Area columns with 'Group Area:' prefix"""
                pattern = self.group_area_pattern_var.get().strip()
                if not pattern:
                    messagebox.showwarning("No Pattern", "Please enter a pattern first")
                    return
                
                detected = ColumnDetector.detect_group_area_columns_by_pattern(self.df, pattern)
                if not detected:
                    messagebox.showwarning("No Matches", f"No columns found matching pattern '{pattern}'")
                    return
                
                # Rename all matching columns with 'Group Area:' prefix
                renamed_count = 0
                new_columns = {}
                renamed_mapping = {}  # Track old -> new column names
                
                for col in self.df.columns:
                    if col in detected:
                        # Add "Group Area:" prefix if not already present
                        if not col.startswith('Group Area:'):
                            new_col_name = f"Group Area: {col}"
                            new_columns[col] = new_col_name
                            renamed_mapping[col] = new_col_name
                            renamed_count += 1
                        else:
                            new_columns[col] = col
                            renamed_mapping[col] = col
                    else:
                        new_columns[col] = col
                
                if renamed_count > 0 or len(detected) > 0:
                    # Rename columns in dataframe
                    if renamed_count > 0:
                        self.df = self.df.rename(columns=new_columns)
                        logger.info(f"Assigned {renamed_count} Group Area columns with 'Group Area:' prefix")
                    
                    # Update assignments dictionary: create new entries for renamed columns
                    new_assignments = {}
                    for old_col, var in list(self.assignments.items()):
                        new_col = renamed_mapping.get(old_col, old_col)
                        if new_col != old_col:
                            # Column was renamed, create new assignment for new name
                            new_assignments[new_col] = tk.StringVar(value='Group Area')
                        else:
                            # Column not renamed, keep as is
                            new_assignments[old_col] = var
                    
                    self.assignments = new_assignments
                    
                    # Rebuild the UI with new column names and assignments
                    self._rebuild_column_rows()
                    
                    # Update area count validation
                    self._update_area_count()
                    
                    messagebox.showinfo(
                        "✅ Done",
                        f"Assigned {len(detected)} column(s) as Group Area"
                    )
                else:
                    messagebox.showinfo(
                        "ℹ️ Already Assigned",
                        "All Group Area columns already assigned"
                    )
            
            def reset_group_areas():
                """Reset Group Area assignments"""
                # Reset all Group Area dropdowns to Ignore
                reset_count = 0
                for excel_col, var in self.assignments.items():
                    if var.get() == 'Group Area':
                        var.set('Ignore')
                        reset_count += 1
                
                # Update area count validation
                self._update_area_count()
                
                messagebox.showinfo(
                    "✅ Reset",
                    f"Reset {reset_count} Group Area assignment(s)"
                )
            
            assign_button = ttk.Button(
                button_frame,
                text="✨ Assign Columns",
                command=assign_all_group_areas,
                state='disabled'
            )
            assign_button.pack(side='left', padx=5, pady=0)
            
            reset_button = ttk.Button(
                button_frame,
                text="🔄 Reset",
                command=reset_group_areas
            )
            reset_button.pack(side='left', padx=0)
            
            # Trigger initial pattern check if auto-pattern is set
            if auto_pattern:
                self.after(100, apply_pattern)
        
        # Lipid Area pattern input (for lipid data cleaning)
        if self.tab_config.get('detect_lipid_area', False):
            lipid_pattern_frame = ttk.LabelFrame(main_frame, text='🔍 Lipid Area Pattern (optional)', padding=5)
            lipid_pattern_frame.pack(fill='x', pady=(0, 10))
            
            ttk.Label(
                lipid_pattern_frame,
                text='Enter pattern(s) to search for Area columns (e.g., "Area", "Area_, Height", etc.). Use commas to separate multiple patterns.',
                font=('Arial', 8),
                foreground='#7f8c8d'
            ).pack(anchor='w', padx=5, pady=(0, 3))
            
            # Pre-fill with auto-detected pattern if available
            auto_lipid_pattern = getattr(self, 'lipid_area_auto_pattern', 'Area[')
            self.lipid_area_pattern_var = tk.StringVar(value=auto_lipid_pattern)
            lipid_pattern_entry = ttk.Entry(
                lipid_pattern_frame,
                textvariable=self.lipid_area_pattern_var,
                width=50
            )
            lipid_pattern_entry.pack(fill='x', padx=5, pady=(0, 3))
            
            def apply_lipid_pattern(*args):
                """Apply pattern and show matching Area columns for lipids"""
                pattern = self.lipid_area_pattern_var.get().strip()
                if pattern:
                    # Support comma-separated patterns
                    patterns = [p.strip() for p in pattern.split(',') if p.strip()]
                    detected_lipid = []
                    for col in self.df.columns:
                        if isinstance(col, str):
                            # Check if column starts with any pattern (case-insensitive)
                            for p in patterns:
                                if col.lower().startswith(p.lower()):
                                    detected_lipid.append(col)
                                    break
                    
                    # Also count Grade columns
                    detected_grade = [col for col in self.df.columns 
                                     if isinstance(col, str) and col.startswith('Grade[') and col.endswith(']')]
                    
                    if detected_lipid:
                        status_text = f"Found {len(detected_lipid)} Area column(s): {', '.join(detected_lipid[:3])}"
                        if len(detected_lipid) > 3:
                            status_text += f" ... and {len(detected_lipid) - 3} more"
                        # Add Grade column info
                        if detected_grade:
                            status_text += f"\nFound {len(detected_grade)} Grade column(s) for quality filtering"
                        lipid_pattern_status.config(text=status_text, foreground='#27ae60')
                        lipid_assign_button.config(state='normal')
                    else:
                        lipid_pattern_status.config(text="No columns matched this pattern", foreground='#e74c3c')
                        lipid_assign_button.config(state='disabled')
                else:
                    lipid_pattern_status.config(text="", foreground='#7f8c8d')
                    lipid_assign_button.config(state='disabled')
            
            self.lipid_area_pattern_var.trace('w', apply_lipid_pattern)
            
            lipid_pattern_status = ttk.Label(
                lipid_pattern_frame,
                text="",
                font=('Arial', 8),
                foreground='#7f8c8d'
            )
            lipid_pattern_status.pack(anchor='w', padx=5)
            
            # Button frame for Assign Columns button
            lipid_button_frame = ttk.Frame(lipid_pattern_frame)
            lipid_button_frame.pack(fill='x', padx=5, pady=(5, 0))
            
            def assign_all_lipid_areas():
                """Assign all detected Area columns for lipids"""
                pattern = self.lipid_area_pattern_var.get().strip()
                if not pattern:
                    messagebox.showwarning("No Pattern", "Please enter a pattern first")
                    return
                
                # Support comma-separated patterns
                patterns = [p.strip() for p in pattern.split(',') if p.strip()]
                detected_lipid = []
                for col in self.df.columns:
                    if isinstance(col, str):
                        # Check if column starts with any pattern (case-insensitive)
                        for p in patterns:
                            if col.lower().startswith(p.lower()):
                                detected_lipid.append(col)
                                break
                
                if not detected_lipid:
                    messagebox.showwarning("No Matches", f"No columns found matching pattern(s) '{pattern}'")
                    return
                
                # Update assignments dictionary: mark all as Area
                for col in detected_lipid:
                    if col in self.assignments:
                        self.assignments[col].set('Area')
                
                # Rebuild the UI to reflect the changes
                self._rebuild_column_rows()
                
                # Update area count validation
                self._update_area_count()
                
                messagebox.showinfo(
                    "✅ Done",
                    f"Assigned {len(detected_lipid)} column(s) as Area"
                )
            
            def reset_lipid_areas():
                """Reset Area assignments for lipids"""
                # Reset all Area dropdowns to Ignore
                reset_count = 0
                for excel_col, var in self.assignments.items():
                    if var.get() == 'Area':
                        var.set('Ignore')
                        reset_count += 1
                
                # Update area count validation
                self._update_area_count()
                
                messagebox.showinfo(
                    "✅ Reset",
                    f"Reset {reset_count} Area assignment(s)"
                )
            
            lipid_assign_button = ttk.Button(
                lipid_button_frame,
                text="✨ Assign Columns",
                command=assign_all_lipid_areas,
                state='disabled'
            )
            lipid_assign_button.pack(side='left', padx=5, pady=0)
            
            lipid_reset_button = ttk.Button(
                lipid_button_frame,
                text="🔄 Reset",
                command=reset_lipid_areas
            )
            lipid_reset_button.pack(side='left', padx=0)
            
            # Trigger initial pattern check if auto-pattern is set
            if auto_lipid_pattern:
                self.after(100, apply_lipid_pattern)
        
        # === GRADE PATTERN ASSIGNMENT (for lipid cleaning) ===
        if self.tab_type == 'lipid_cleaning':
            grade_pattern_frame = ttk.LabelFrame(
                main_frame,
                text='📊 Grade Pattern Assignment (optional)',
                padding=5
            )
            grade_pattern_frame.pack(fill='x', pady=(0, 5))
            
            ttk.Label(
                grade_pattern_frame,
                text='Enter pattern to search for Grade columns (e.g., "Grade["):',
                font=('Arial', 8)
            ).pack(anchor='w', padx=5, pady=(0, 3))
            
            # Pre-fill with auto-detected pattern if available
            auto_grade_pattern = 'Grade[' if getattr(self, 'auto_assign_lipid_grade_cols', []) else ''
            self.lipid_grade_pattern_var = tk.StringVar(value=auto_grade_pattern)
            grade_pattern_entry = ttk.Entry(
                grade_pattern_frame,
                textvariable=self.lipid_grade_pattern_var,
                width=50
            )
            grade_pattern_entry.pack(fill='x', padx=5, pady=(0, 3))
            
            def apply_grade_pattern(*args):
                """Apply pattern and show matching Grade columns for lipids"""
                pattern = self.lipid_grade_pattern_var.get().strip()
                if pattern:
                    # Support comma-separated patterns
                    patterns = [p.strip() for p in pattern.split(',') if p.strip()]
                    detected_grade = []
                    for col in self.df.columns:
                        if isinstance(col, str):
                            # Check if column starts with any pattern (case-insensitive)
                            for p in patterns:
                                if col.lower().startswith(p.lower()):
                                    detected_grade.append(col)
                                    break
                    
                    if detected_grade:
                        status_text = f"Found {len(detected_grade)} Grade column(s): {', '.join(detected_grade[:3])}"
                        if len(detected_grade) > 3:
                            status_text += f" ... and {len(detected_grade) - 3} more"
                        grade_pattern_status.config(text=status_text, foreground='#27ae60')
                        grade_assign_button.config(state='normal')
                    else:
                        grade_pattern_status.config(text="No columns matched this pattern", foreground='#e74c3c')
                        grade_assign_button.config(state='disabled')
                else:
                    grade_pattern_status.config(text="", foreground='#7f8c8d')
                    grade_assign_button.config(state='disabled')
            
            self.lipid_grade_pattern_var.trace('w', apply_grade_pattern)
            
            grade_pattern_status = ttk.Label(
                grade_pattern_frame,
                text="",
                font=('Arial', 8),
                foreground='#7f8c8d'
            )
            grade_pattern_status.pack(anchor='w', padx=5)
            
            # Button frame for Assign Columns button
            grade_button_frame = ttk.Frame(grade_pattern_frame)
            grade_button_frame.pack(fill='x', padx=5, pady=(5, 0))
            
            def assign_all_lipid_grades():
                """Assign all detected Grade columns for lipids"""
                pattern = self.lipid_grade_pattern_var.get().strip()
                if not pattern:
                    messagebox.showwarning("No Pattern", "Please enter a pattern first")
                    return
                
                # Support comma-separated patterns
                patterns = [p.strip() for p in pattern.split(',') if p.strip()]
                detected_grade = []
                for col in self.df.columns:
                    if isinstance(col, str):
                        # Check if column starts with any pattern (case-insensitive)
                        for p in patterns:
                            if col.lower().startswith(p.lower()):
                                detected_grade.append(col)
                                break
                
                if not detected_grade:
                    messagebox.showwarning("No Matches", f"No columns found matching pattern(s) '{pattern}'")
                    return
                
                # Update assignments dictionary: mark all as Grade
                for col in detected_grade:
                    if col in self.assignments:
                        self.assignments[col].set('Grade')
                
                # Rebuild the UI to reflect the changes
                self._rebuild_column_rows()
                
                messagebox.showinfo(
                    "✅ Done",
                    f"Assigned {len(detected_grade)} column(s) as Grade"
                )
            
            def reset_lipid_grades():
                """Reset Grade assignments for lipids"""
                # Reset all Grade dropdowns to Ignore
                reset_count = 0
                for excel_col, var in self.assignments.items():
                    if var.get() == 'Grade':
                        var.set('Ignore')
                        reset_count += 1
                
                messagebox.showinfo(
                    "✅ Reset",
                    f"Reset {reset_count} Grade assignment(s)"
                )
            
            grade_assign_button = ttk.Button(
                grade_button_frame,
                text="✨ Assign Columns",
                command=assign_all_lipid_grades,
                state='disabled'
            )
            grade_assign_button.pack(side='left', padx=5, pady=0)
            
            grade_reset_button = ttk.Button(
                grade_button_frame,
                text="🔄 Reset",
                command=reset_lipid_grades
            )
            grade_reset_button.pack(side='left', padx=0)
            
            # Trigger initial pattern check if auto-pattern is set
            if auto_grade_pattern:
                self.after(100, apply_grade_pattern)
        
        # ========== SAMPLE/AREA COLUMN COUNT VALIDATION SECTION ==========
        # For data_cleaning: validate Group Area columns
        # For lipid_cleaning: validate Area columns
        if self.tab_config.get('detect_group_area', False) or self.tab_config.get('detect_lipid_area', False):
            validation_frame = ttk.LabelFrame(main_frame, text='🔒 Sample/Area Column Validation (REQUIRED)', padding=5)
            validation_frame.pack(fill='x', pady=(0, 10))
            
            # Determine column type name based on tab type
            if self.tab_config.get('detect_group_area', False):
                col_type_name = "Group Area"
                self._area_col_type = 'Group Area'
            else:
                col_type_name = "Area"
                self._area_col_type = 'Area'
            
            # Calculate initial count of assigned columns from current assignments
            # This ensures the count reflects auto-detected columns immediately
            initial_count = 0
            for col, var in self.assignments.items():
                if var.get() == col_type_name:
                    initial_count += 1
            # Fallback to sample_cols if no assignments yet
            if initial_count == 0 and self.sample_cols:
                initial_count = len(self.sample_cols)
            self._detected_area_count = initial_count
            
            # Validation row
            validation_row = ttk.Frame(validation_frame)
            validation_row.pack(fill='x', padx=3, pady=(2, 3))
            
            ttk.Label(
                validation_row,
                text=f'Expected {col_type_name} columns:',
                font=('Arial', 8, 'bold')
            ).pack(side='left', padx=(0, 5))
            
            # Entry for expected count
            self._expected_area_count_var = tk.StringVar(value='')
            expected_entry = ttk.Entry(
                validation_row,
                textvariable=self._expected_area_count_var,
                width=8,
                font=('Arial', 10)
            )
            expected_entry.pack(side='left', padx=(0, 10))
            
            # Dynamic detected count display (updates in real-time as user changes assignments)
            self._area_count_display = ttk.Label(
                validation_row,
                text=f'Assigned: {self._detected_area_count}',
                font=('Arial', 9, 'bold'),
                foreground='#3498db'
            )
            self._area_count_display.pack(side='left')
            
            # Validation status label - starts as PENDING (user must enter count)
            self._area_validation_status = ttk.Label(
                validation_frame,
                text='⏳ Enter expected count to validate',
                font=('Arial', 8, 'bold'),
                foreground='#e67e22'
            )
            self._area_validation_status.pack(anchor='w', padx=3, pady=(0, 2))
            
            # Validation function
            def validate_area_count(*args):
                """Check if entered count matches current assigned count - REQUIRED for confirmation"""
                try:
                    entered = self._expected_area_count_var.get().strip()
                    if not entered:
                        # REQUIRED: User must enter expected count
                        self._area_validation_status.config(
                            text='⏳ Enter expected count to validate',
                            foreground='#e67e22'
                        )
                        self._area_count_valid = False  # REQUIRED - not valid if empty
                        return
                    
                    entered_count = int(entered)
                    current_assigned = self._detected_area_count
                    
                    if entered_count == current_assigned:
                        self._area_validation_status.config(
                            text=f'✅ Match! {entered_count} = {current_assigned}',
                            foreground='#27ae60'
                        )
                        self._area_count_valid = True
                    else:
                        self._area_validation_status.config(
                            text=f'❌ Mismatch! Expected {entered_count}, Assigned {current_assigned}',
                            foreground='#e74c3c'
                        )
                        self._area_count_valid = False
                except ValueError:
                    self._area_validation_status.config(
                        text='❌ Enter valid number',
                        foreground='#e74c3c'
                    )
                    self._area_count_valid = False
            
            # Store validation function for use by assignment change handler
            self._validate_area_count_func = validate_area_count
            
            # Initialize validation state - REQUIRED (not valid until user enters matching count)
            self._area_count_valid = False
            
            # Bind validation to entry changes
            self._expected_area_count_var.trace('w', validate_area_count)
        
        # Button frame at TOP RIGHT - compact
        btn_frame = tk.Frame(main_frame, bg='#f0f0f0')
        btn_frame.pack(side='top', fill='x', pady=(0, 5), anchor='e')
        
        # Status label (left side)
        self.status_label = ttk.Label(
            btn_frame,
            text='',
            foreground='green'
        )
        self.status_label.pack(side='left', padx=5)
        
        # Optional Skip button (for workflows where a sheet/polarity is optional)
        if self.allow_skip:
            skip_btn = tk.Button(
                btn_frame,
                text='⏭ Skip',
                command=self._on_skip,
                bg='#f1c40f',
                fg='black',
                font=('Arial', 11, 'bold'),
                padx=16,
                pady=8,
                relief='raised',
                cursor='hand2'
            )
            skip_btn.pack(side='right', padx=5)
        
        # Cancel button (red) - pack from right so it's rightmost
        cancel_btn = tk.Button(
            btn_frame,
            text='❌ Cancel',
            command=self._on_cancel,
            bg='#e74c3c',
            fg='white',
            font=('Arial', 11, 'bold'),
            padx=20,
            pady=8,
            relief='raised',
            cursor='hand2'
        )
        cancel_btn.pack(side='right', padx=5)
        
        # Confirm button (prominent green) - pack from right
        confirm_btn = tk.Button(
            btn_frame,
            text='✅ Confirm',
            command=self._on_confirm,
            bg='#27ae60',
            fg='white',
            font=('Arial', 11, 'bold'),
            padx=20,
            pady=8,
            relief='raised',
            cursor='hand2'
        )
        confirm_btn.pack(side='right', padx=5)
        
        # Column assignments frame with scrollbar (pack AFTER buttons) - maximize space
        scroll_frame = ttk.LabelFrame(main_frame, text='Column Assignment', padding=2)
        scroll_frame.pack(fill='both', expand=True, pady=(0, 2))

        # Create scrollbar
        scrollbar = ttk.Scrollbar(scroll_frame)
        scrollbar.grid(row=0, column=1, sticky='ns')

        # Create canvas for scrolling
        # Match canvas background to dialog background (avoid ttk cget background errors)
        try:
            parent_bg = self.cget('background')
        except Exception:
            parent_bg = '#f0f0f0'
        canvas = tk.Canvas(
            scroll_frame,
            yscrollcommand=scrollbar.set,
            highlightthickness=0,
            bg=parent_bg,
            relief='flat',
            borderwidth=0,
        )
        scrollbar.config(command=canvas.yview)
        canvas.grid(row=0, column=0, sticky='nsew')
        scroll_frame.grid_rowconfigure(0, weight=1)
        scroll_frame.grid_columnconfigure(0, weight=1)

        # Create scrollable frame
        scrollable_frame = ttk.Frame(canvas)
        canvas.create_window((0, 0), window=scrollable_frame, anchor='nw')
        
        # Store references for rebuilding later
        self.scrollable_frame = scrollable_frame
        self.scroll_canvas = canvas

        # Header
        header_frame = ttk.Frame(scrollable_frame)
        header_frame.pack(fill='x', pady=(0, 4))
        self.header_frame = header_frame
        ttk.Label(header_frame, text='Excel Column', font=('Arial', 9, 'bold'), width=30).pack(side='left', padx=4)
        ttk.Label(header_frame, text='Assign To ↓', font=('Arial', 9, 'bold'), width=24).pack(side='left', padx=4)
        
        # Populate rows with interactive dropdowns
        self._populate_interactive_rows(scrollable_frame)
        
        # Update scroll region
        def on_frame_configure(event=None):
            canvas.configure(scrollregion=canvas.bbox('all'))
        scrollable_frame.bind('<Configure>', on_frame_configure)
    
    def _populate_interactive_rows(self, parent_frame):
        """Populate interactive dropdown rows for each Excel column."""
        # Build list of available column types for this tab (required + recommended + optional)
        available_types = ['Ignore'] + self.tab_config['required'] + self.tab_config.get('recommended', []) + self.tab_config.get('optional', [])
        hardcoded_types = self.tab_config.get('hardcoded', [])
        
        # Get columns that should be auto-assigned as Group Area (those with "Group Area:" prefix)
        auto_group_area_cols = getattr(self, 'auto_assign_group_area_cols', [])
        
        # Get columns that should be auto-assigned as Area (lipid Area[ ] columns)
        auto_lipid_area_cols = getattr(self, 'auto_assign_lipid_area_cols', [])
        
        # Get columns that should be auto-assigned as Grade (lipid Grade[ ] columns)
        auto_lipid_grade_cols = getattr(self, 'auto_assign_lipid_grade_cols', [])
        
        # Auto-detect hardcoded pathway columns (for network analysis)
        hardcoded_cols = {}
        if hardcoded_types:
            for col in self.df.columns:
                col_str = str(col)
                for hardcoded_type in hardcoded_types:
                    # Match exact or with underscore/space variations
                    if col_str == hardcoded_type or col_str.replace(' ', '_') == hardcoded_type.replace(' ', '_'):
                        hardcoded_cols[col] = hardcoded_type
                        break
        
        # Check if we already have assignments (rebuilding) - preserve them
        preserve_existing = hasattr(self, 'assignments') and self.assignments
        
        # PRIORITY: Check for exact 'Name' column first for Feature ID assignment
        # This ensures 'Name' is prioritized over 'metabolite_id', 'compound_id', etc.
        exact_name_col = None
        for excel_col in self.df.columns:
            if str(excel_col).strip().lower() == 'name':
                exact_name_col = excel_col
                break
        
        # First pass: auto-detect all columns
        all_detections = {}
        for excel_col in self.df.columns:
            # If column already has "Group Area:" prefix, mark it as Group Area
            if excel_col in auto_group_area_cols:
                all_detections[excel_col] = 'Group Area'
            # If column matches lipid Area pattern, mark it as Area
            elif excel_col in auto_lipid_area_cols:
                all_detections[excel_col] = 'Area'
            # If column matches lipid Grade pattern, mark it as Grade
            elif excel_col in auto_lipid_grade_cols:
                all_detections[excel_col] = 'Grade'
            # For statistics tabs, auto-detect sample and feature columns
            elif self.tab_type in ['statistics_metabolite', 'statistics_lipid']:
                # Check if this column is in our detected sample columns
                if excel_col in self.sample_cols:
                    all_detections[excel_col] = 'Sample Column'
                # Otherwise, try to detect specific feature type (LipidID, Class, Class_name, etc.)
                else:
                    # Use the specific type detection to identify known metadata columns
                    detected_type = self._auto_detect_column_type(excel_col)
                    # If a specific type was detected (not 'Ignore'), use it
                    if detected_type != 'Ignore':
                        all_detections[excel_col] = detected_type
                    else:
                        # Otherwise, it's a generic feature/metadata column
                        all_detections[excel_col] = 'Feature Column'
            else:
                all_detections[excel_col] = self._auto_detect_column_type(excel_col)
        
        # Identify preferred P-Value column before resolving duplicates
        pvalue_candidates = [col for col, detected in all_detections.items() if detected == 'P-Value']
        preferred_pvalue_col = ColumnDetector.select_best_pvalue_column(pvalue_candidates)

        # Second pass: remove duplicates, keep only highest priority
        assigned_types = {}
        final_assignments = {}
        
        # PRIORITY 1: Use existing_assignments from previous session if available
        if self.existing_assignments:
            # Convert existing assignments (col_type: excel_col) to (excel_col: col_type)
            for col_type, col_name in self.existing_assignments.items():
                if col_name and col_name in self.df.columns:
                    if isinstance(col_name, list):
                        # Multiple columns assigned to same type
                        for c in col_name:
                            if c in self.df.columns:
                                final_assignments[c] = col_type
                    else:
                        final_assignments[col_name] = col_type
                        assigned_types[col_type] = col_name
        
        # If we found an exact 'Name' column, pre-assign it as Feature ID (unless already assigned)
        if exact_name_col and 'Feature ID' in (self.tab_config['required'] + self.tab_config.get('optional', [])):
            if 'Feature ID' not in assigned_types:
                assigned_types['Feature ID'] = exact_name_col
        
        for excel_col in self.df.columns:
            # Skip if already assigned from existing_assignments
            if excel_col in final_assignments:
                continue
                
            # If rebuilding and this column already has an assignment, preserve it
            if preserve_existing and excel_col in self.assignments:
                final_assignments[excel_col] = self.assignments[excel_col].get()
                continue
            
            # If this is the exact 'Name' column, assign it as Feature ID
            if excel_col == exact_name_col and 'Feature ID' not in final_assignments.values():
                final_assignments[excel_col] = 'Feature ID'
                continue
                
            detected = all_detections[excel_col]
            
            # Group Area, Area, Grade, Sample Column, and Feature Column can have multiple columns
            if detected in ['Group Area', 'Area', 'Grade', 'Sample Column', 'Feature Column']:
                final_assignments[excel_col] = detected
            elif detected == 'P-Value':
                # Only assign the preferred P-Value column (adjusted > raw > F-statistic)
                if preferred_pvalue_col is None and 'P-Value' not in assigned_types:
                    final_assignments[excel_col] = 'P-Value'
                    assigned_types['P-Value'] = excel_col
                elif excel_col == preferred_pvalue_col and 'P-Value' not in assigned_types:
                    final_assignments[excel_col] = 'P-Value'
                    assigned_types['P-Value'] = excel_col
                else:
                    final_assignments[excel_col] = 'Ignore'
            # If this type hasn't been assigned yet, use this column
            elif detected != 'Ignore' and detected not in assigned_types:
                final_assignments[excel_col] = detected
                assigned_types[detected] = excel_col
            else:
                # Already assigned or Ignore
                final_assignments[excel_col] = 'Ignore'
        
        # Create rows with interactive dropdowns
        for excel_col in self.df.columns:
            # Skip hardcoded columns - don't show them in the dialog
            if excel_col in hardcoded_cols:
                # Store assignment but don't create UI row
                var = tk.StringVar(value=hardcoded_cols[excel_col])
                self.assignments[excel_col] = var
                continue
            
            # Create a StringVar for this column (use preserved value if rebuilding)
            var = tk.StringVar(value=final_assignments.get(excel_col, 'Ignore'))
            self.assignments[excel_col] = var
            
            # Create row frame
            row_frame = ttk.Frame(parent_frame)
            row_frame.pack(fill='x', padx=10, pady=5)
            
            # Column name label
            ttk.Label(row_frame, text=excel_col, width=40, anchor='w').pack(side='left', padx=(5, 10))
            
            # Dropdown combobox - Interactive combobox per row
            combo = ttk.Combobox(
                row_frame,
                textvariable=var,
                values=available_types,
                state='readonly',
                width=28
            )
            combo.pack(side='left', padx=5)
            
            # Bind combobox change to update area count validation
            def on_combo_change(event=None, col=excel_col, v=var):
                """Update area count when assignment changes"""
                self._update_area_count()
            
            combo.bind('<<ComboboxSelected>>', on_combo_change)
            
            # Store reference for later
            if not hasattr(self, 'combo_widgets'):
                self.combo_widgets = {}
            self.combo_widgets[excel_col] = combo
    
    def _update_area_count(self):
        """Update the detected area/sample column count based on current assignments."""
        if not hasattr(self, '_detected_area_count'):
            return
        
        # Determine which column type to count
        col_type = getattr(self, '_area_col_type', 'Area')
        
        # Count columns assigned to this type
        count = 0
        for excel_col, var in self.assignments.items():
            if var.get() == col_type:
                count += 1
        
        # Update the count
        self._detected_area_count = count
        
        # Update display label if it exists
        if hasattr(self, '_area_count_display'):
            self._area_count_display.config(text=f'Assigned: {count}')
        
        # Re-validate if validation function exists
        if hasattr(self, '_validate_area_count_func'):
            self._validate_area_count_func()
    
    def _rebuild_column_rows(self):
        """Rebuild the column assignment rows after columns have been renamed."""
        # Find the scrollable frame containing the rows
        if not hasattr(self, 'scrollable_frame'):
            return
        
        # Clear existing rows (keep header)
        for widget in list(self.scrollable_frame.winfo_children()):
            if widget != getattr(self, 'header_frame', None):
                widget.destroy()
        
        # Repopulate with new assignments
        self._populate_interactive_rows(self.scrollable_frame)
        
        # Update scroll region
        self.scrollable_frame.update_idletasks()
    
    def _auto_detect_column_type(self, excel_col: str) -> str:
        """Auto-detect the best column type for an Excel column using proper detection patterns."""
        # Get available column types for this tab (required + recommended + optional)
        available_types = self.tab_config['required'] + self.tab_config.get('recommended', []) + self.tab_config.get('optional', [])
        
        # For statistics tabs, detect required and optional feature columns
        if self.tab_type in ['statistics_metabolite', 'statistics_lipid']:
            # Check for Feature ID (metabolite Name)
            if 'Feature ID' in available_types and ColumnDetector.detect_column(excel_col, 'Feature ID'):
                return 'Feature ID'
            # Check for LipidID (highest priority for lipid mode)
            if 'LipidID' in available_types and ColumnDetector.detect_column(excel_col, 'LipidID'):
                return 'LipidID'
            # Check for Class (important for lipid classification)
            if 'Class' in available_types and ColumnDetector.detect_column(excel_col, 'Class'):
                return 'Class'
            # Check for Class_name (important for lipid classification and pathway lookup)
            if 'Class_name' in available_types and ColumnDetector.detect_column(excel_col, 'Class_name'):
                return 'Class_name'
            # Check for other optional columns like Formula, Pathways, etc.
            for col_type in available_types:
                if col_type not in ['Sample Column', 'Feature Column']:  # Skip sample/feature identifiers
                    if ColumnDetector.detect_column(excel_col, col_type):
                        return col_type
            # Don't auto-detect Sample/Feature identifiers here - already done above
            return 'Ignore'
        
        # Check each available type using the proper regex detection
        for col_type in available_types:
            if ColumnDetector.detect_column(excel_col, col_type):
                return col_type
        
        return 'Ignore'
    
    def _on_confirm(self):
        """Confirm assignments and close."""
        # VALIDATION: Check area/sample column count validation (for data cleaning tabs)
        if hasattr(self, '_area_count_valid') and hasattr(self, '_detected_area_count'):
            # Get the entered value for error message
            entered = getattr(self, '_expected_area_count_var', tk.StringVar()).get().strip()
            col_type_name = getattr(self, '_area_col_type', 'Area')
            
            # MANDATORY: User MUST enter expected count
            if not entered:
                messagebox.showerror(
                    f"❌ {col_type_name} Column Count Required",
                    f"You must enter the expected number of {col_type_name} columns.\n\n"
                    f"Currently assigned: {self._detected_area_count} column(s)\n\n"
                    f"Please enter the expected count in the validation field\n"
                    f"to confirm your column assignments are correct."
                )
                return  # Don't close - user must fix
            
            # MANDATORY: Count must match
            if not self._area_count_valid:
                try:
                    entered_count = int(entered)
                    messagebox.showerror(
                        f"❌ {col_type_name} Column Count Mismatch",
                        f"Expected {col_type_name} columns ({entered_count}) does not match assigned ({self._detected_area_count}).\n\n"
                        f"This usually means:\n"
                        f"• Wrong columns are assigned as '{col_type_name}'\n"
                        f"• Some columns need to be manually changed\n"
                        f"• The pattern might need adjustment\n\n"
                        f"Please adjust column assignments in the list below,\n"
                        f"or use the pattern field to re-detect."
                    )
                except ValueError:
                    messagebox.showerror(
                        "❌ Invalid Input",
                        f"Please enter a valid number for the expected column count."
                    )
                return  # Don't close - user must fix
        
        # VALIDATION: MANDATORY - Group Area columns must be assigned for data cleaning
        if self.tab_config.get('detect_group_area', False):
            # Count how many columns are ACTUALLY assigned as Group Area (by dropdown value, not column name)
            # This respects user's manual changes to assignments
            group_area_count = sum(1 for var in self.assignments.values() if var.get() == 'Group Area')
            
            if group_area_count == 0:
                messagebox.showerror(
                    "❌ No Group Area Columns Assigned",
                    f"You must assign Group Area columns before proceeding.\n\n"
                    f"Data cleaning requires Group Area columns to identify sample data.\n"
                    f"Without them, the cleaning process will fail.\n\n"
                    f"💡 HOW TO FIX:\n"
                    f"1. Enter a pattern in the 'Group Area Pattern' field\n"
                    f"   (e.g., 'Group Area', 'LF_pos', 'Peak')\n"
                    f"2. Click the '✨ Assign Columns' button\n"
                    f"3. Verify the count matches your expected number\n"
                    f"4. Then click Confirm again"
                )
                return  # Don't close - user must fix
        
        # VALIDATION: MANDATORY - Area columns must be assigned for lipid cleaning
        if self.tab_config.get('detect_lipid_area', False):
            # Count how many columns are assigned as Area
            # Count how many columns are ACTUALLY assigned as Area (by dropdown value, not column name)
            area_count = sum(1 for var in self.assignments.values() if var.get() == 'Area')
            
            if area_count == 0:
                messagebox.showerror(
                    "❌ No Area Columns Assigned",
                    f"You must assign Area columns before proceeding.\n\n"
                    f"Lipid cleaning requires Area columns to identify sample data.\n"
                    f"Without them, the cleaning process will fail.\n\n"
                    f"💡 HOW TO FIX:\n"
                    f"1. Enter a pattern in the 'Lipid Area Pattern' field\n"
                    f"   (e.g., 'Area[')\n"
                    f"2. Click the '✨ Assign Columns' button\n"
                    f"3. Verify the count matches your expected number\n"
                    f"4. Then click Confirm again"
                )
                return  # Don't close - user must fix
        
        # If lipid area/grade patterns are provided and user didn't click Assign, apply them once here
        if self.tab_config.get('detect_lipid_area', False):
            # Area pattern
            try:
                if hasattr(self, 'lipid_area_pattern_var'):
                    patterns = [p.strip() for p in self.lipid_area_pattern_var.get().split(',') if p.strip()]
                    if patterns:
                        for col in self.df.columns:
                            if isinstance(col, str) and self.assignments.get(col, tk.StringVar(value='Ignore')).get() == 'Ignore':
                                for p in patterns:
                                    if col.lower().startswith(p.lower()):
                                        self.assignments[col].set('Area')
                                        break
                # Grade pattern
                if hasattr(self, 'lipid_grade_pattern_var'):
                    patterns = [p.strip() for p in self.lipid_grade_pattern_var.get().split(',') if p.strip()]
                    if patterns:
                        for col in self.df.columns:
                            if isinstance(col, str) and self.assignments.get(col, tk.StringVar(value='Ignore')).get() == 'Ignore':
                                for p in patterns:
                                    if col.lower().startswith(p.lower()):
                                        self.assignments[col].set('Grade')
                                        break
            except Exception:
                pass

        # Build assignments dict from tree
        assignments = {}
        group_area_assigned_cols = []  # Track which columns are assigned as Group Area
        
        for excel_col, var in self.assignments.items():
            selected = var.get()
            col_type = selected if selected != 'Ignore' else None
            assignments[excel_col] = col_type
            
            # Collect Group Area assignments
            if col_type == 'Group Area':
                group_area_assigned_cols.append(excel_col)
        
        # Apply trailing text removal if requested
        if hasattr(self, 'remove_trailing_text_var') and self.remove_trailing_text_var.get():
            # Common trailing patterns to remove
            trailing_patterns = [
                '_pos', '_positive', '_Pos', '_Positive', '_POS', '_POSITIVE',
                '_neg', '_negative', '_Neg', '_Negative', '_NEG', '_NEGATIVE',
                '.pos', '.positive', '.Pos', '.Positive', '.POS', '.POSITIVE',
                '.neg', '.negative', '.Neg', '.Negative', '.NEG', '.NEGATIVE',
                ' pos', ' positive', ' Pos', ' Positive', ' POS', ' POSITIVE',
                ' neg', ' negative', ' Neg', ' Negative', ' NEG', ' NEGATIVE'
            ]
            
            rename_dict = {}
            for col in self.df.columns:
                if isinstance(col, str):
                    new_col = col
                    for pattern in trailing_patterns:
                        if col.endswith(pattern):
                            new_col = col[:-len(pattern)]
                            break
                    if new_col != col:
                        rename_dict[col] = new_col
            
            if rename_dict:
                self.df = self.df.rename(columns=rename_dict)
                logger.info(f"Removed trailing text from {len(rename_dict)} columns: {list(rename_dict.keys())[:5]}")
                
                # Update assignments dict with new column names
                new_assignments = {}
                for old_col, col_type in assignments.items():
                    new_col = rename_dict.get(old_col, old_col)
                    new_assignments[new_col] = col_type
                assignments = new_assignments
                
                # Update sample columns list
                if self.sample_cols:
                    self.sample_cols = [rename_dict.get(col, col) for col in self.sample_cols]
        
        # IMPORTANT: Ensure all Group Area assigned columns have "Group Area:" prefix
        # If a column is assigned as Group Area but doesn't have the prefix, add it
        new_columns = {}
        for col in self.df.columns:
            if col in group_area_assigned_cols and not col.startswith('Group Area:'):
                new_col_name = f"Group Area: {col}"
                new_columns[col] = new_col_name
            else:
                new_columns[col] = col
        
        # Apply any needed renames
        if new_columns and any(old != new for old, new in new_columns.items()):
            self.df = self.df.rename(columns=new_columns)
            logger.info(f"Added 'Group Area:' prefix to {len([x for x in new_columns.values() if x.startswith('Group Area:')])} assigned columns")
            
            # Update assignments dictionary with new column names
            new_assignments = {}
            for old_col, col_type in assignments.items():
                new_col = new_columns.get(old_col, old_col)
                new_assignments[new_col] = col_type
            assignments = new_assignments
        
        # Reverse to: {column_type: excel_column} for consistency
        # NOTE: Group Area is handled via column renaming with "Group Area:" prefix,
        # so we don't include it in the assignments dict
        # For lipid Area/Grade/ObsMz/ObsRt/Delta(PPM), collect all matching columns
        # For statistics Sample/Feature Columns, also collect all matching columns
        reversed_assignments = {}
        pattern_columns = {
            'Area': [], 
            'Grade': [], 
            'ObsMz': [], 
            'ObsRt': [], 
            'Delta(PPM)': [],
            'Sample Column': [],
            'Feature Column': []
        }  # Columns with multiple matches
        
        for excel_col, col_type in assignments.items():
            if col_type and col_type != 'Group Area':  # Skip Group Area - it's handled by column renaming
                if col_type in pattern_columns:
                    # For pattern-based columns, store as list
                    pattern_columns[col_type].append(excel_col)
                else:
                    # For single-value columns, store directly
                    reversed_assignments[col_type] = excel_col
        
        # Now add pattern columns to reversed_assignments
        # If only one match, store as string; if multiple, store as list (for backward compatibility)
        for col_type, cols in pattern_columns.items():
            if cols:
                # Store all matches - even single items as list for consistency with renaming logic
                # Only add if there are actual columns (don't add None for empty lists)
                if len(cols) > 1:
                    reversed_assignments[col_type] = cols
                elif len(cols) == 1:
                    reversed_assignments[col_type] = cols[0]
                # If cols is empty, don't add to reversed_assignments at all

        # Safety filter: Remove any None or empty values that might have slipped through
        reversed_assignments = {k: v for k, v in reversed_assignments.items() if v is not None and v != ''}
        
        # Validate required columns
        missing = []
        
        for req_col in self.tab_config['required']:
            # Special case: For statistics_lipid, if LipidID is required but Class is assigned instead, accept it
            # This handles lipid class sheets where Class is the identifier instead of LipidID
            if req_col == 'LipidID' and req_col not in reversed_assignments:
                # Check if Class is assigned (lipid class sheet scenario)
                if 'Class' in reversed_assignments and reversed_assignments['Class']:
                    # Class is present, use it instead of LipidID
                    continue
            
            if req_col not in reversed_assignments:
                missing.append(req_col)
        
        if missing:
            messagebox.showerror(
                'Missing Required Columns',
                f"Please assign all required columns:\n{', '.join(missing)}"
            )
            return
        
        # Calculate missing columns if enabled
        modified_df = self.df
        if self.auto_calculate:
            can_calculate = self.tab_config.get('can_calculate', [])
            modified_df, reversed_assignments = ColumnCalculator.calculate_missing_columns(
                self.df,
                reversed_assignments,
                can_calculate
            )
        
        # Apply column_mapping from tab config (renames columns to expected names)
        # e.g., for id_annotation: 'Feature ID' -> 'Name'
        column_mapping_config = self.tab_config.get('column_mapping', {})
        if column_mapping_config:
            #print(f"📝 column_mapping_config: {column_mapping_config}")
            #print(f"📝 reversed_assignments before renaming: {reversed_assignments}")
            
            # Build rename dict: {user_col_name: target_name}
            rename_dict = {}
            for col_type, target_name in column_mapping_config.items():
                if col_type in reversed_assignments:
                    user_col = reversed_assignments[col_type]
                    if isinstance(user_col, str):
                        rename_dict[user_col] = target_name
                        print(f"   → Will rename '{user_col}' to '{target_name}'")
                else:
                    print(f"   ⚠️ '{col_type}' not in reversed_assignments (optional column?)")
            
            if rename_dict:
                modified_df = modified_df.rename(columns=rename_dict)
                #print(f"📝 Applied column mapping: {rename_dict}")
                #print(f"📝 Final dataframe columns: {list(modified_df.columns)}")

                # Update reversed_assignments to reflect the new names
                for col_type, target_name in column_mapping_config.items():
                    if col_type in reversed_assignments:
                        reversed_assignments[col_type] = target_name
        
        # Re-detect sample columns EXCLUDING verified statistics columns
        # This ensures that columns verified as log2FC, pvalue, etc. are not included in samples
        exclude_stats_cols = list(reversed_assignments.values()) if reversed_assignments else []
        exclude_stats_cols = [col for col in exclude_stats_cols if isinstance(col, str)]
        
        # For statistics tabs, extract sample columns from the assignments instead of re-detecting
        if self.tab_type in ['statistics_metabolite', 'statistics_lipid']:
            # Extract all columns that were assigned as 'Sample Column' from user selections
            user_selected_sample_cols = []
            for excel_col, var in self.assignments.items():
                if var.get() == 'Sample Column':
                    user_selected_sample_cols.append(excel_col)
            
            # Use user-selected sample columns from the dialog (respects user overrides)
            updated_sample_cols = user_selected_sample_cols if user_selected_sample_cols else self.sample_cols
            logger.info(f"Statistics tab: User selected {len(user_selected_sample_cols)} sample columns, returning {len(updated_sample_cols)} total")
        elif self.tab_config.get('detect_samples', False):
            updated_sample_cols = ColumnDetector.detect_sample_columns(modified_df, exclude_cols=exclude_stats_cols)
        elif self.tab_config.get('detect_group_area', False):
            # For metabolite cleaning with Group Area detection - use ACTUAL assignments, not column names
            # This respects user's manual changes (e.g., unassigning a wrongly detected column)
            updated_sample_cols = [col for col, var in self.assignments.items() if var.get() == 'Group Area']
        elif self.tab_config.get('detect_lipid_area', False):
            # For lipid cleaning: prefer user-assigned Area/Grade columns, fall back to pattern detection
            assigned_area_cols = [col for col, var in self.assignments.items() if var.get() == 'Area']
            assigned_grade_cols = [col for col, var in self.assignments.items() if var.get() == 'Grade']

            area_cols = [col for col in modified_df.columns
                        if isinstance(col, str) and col.startswith('Area[') and col.endswith(']')]
            grade_cols = [col for col in modified_df.columns
                         if isinstance(col, str) and col.startswith('Grade[') and col.endswith(']')]

            if assigned_area_cols or assigned_grade_cols:
                updated_sample_cols = assigned_area_cols + assigned_grade_cols
            else:
                updated_sample_cols = area_cols + grade_cols
        else:
            updated_sample_cols = self.sample_cols
        
        self.result = {
            'assignments': reversed_assignments,
            'dataframe': modified_df,
            'tab_type': self.tab_type,
            'sample_cols': updated_sample_cols,  # Use updated sample columns excluding verified stats
            'selected_sheet': getattr(self, 'sheet_var', None) and self.sheet_var.get(),  # Include selected sheet for ID annotation
            'validation_passed': getattr(self, '_area_count_valid', True),  # Validation status for caller
            'expected_count': getattr(self, '_expected_area_count_var', tk.StringVar()).get().strip() if hasattr(self, '_expected_area_count_var') else None,
            'assigned_count': getattr(self, '_detected_area_count', len(updated_sample_cols)),  # How many were assigned
        }
        self.modified_df = modified_df
        
        self.destroy()
       
    def _on_cancel(self):
        """Cancel and close."""
        self.result = None
        self.destroy()

    def _on_skip(self):
        """Skip this sheet/polarity but keep overall workflow running.

        The caller can detect a skip via result.get('skipped', False).
        """
        self.result = {
            'skipped': True,
            'tab_type': self.tab_type,
        }
        self.destroy()
    
    def wait_for_result(self) -> Optional[Dict[str, Any]]:
        """Wait for dialog to close and return result."""
        self.wait_window()
        return self.result


# ===== PUBLIC FUNCTIONS =====

def show_column_assignment_dialog(
    parent: tk.Tk,
    df: pd.DataFrame,
    tab_type: str = 'general',
    auto_calculate: bool = True,
    existing_assignments: Optional[Dict[str, Optional[str]]] = None,
    excel_file_path: Optional[str] = None,
    dialog_title: Optional[str] = None,
    group_area_pattern: Optional[str] = None,
    detected_sample_cols: Optional[List[str]] = None,
    allow_skip: bool = False,
) -> Optional[Dict[str, Any]]:
    """
    Show the unified column assignment dialog.
    
    Args:
        parent: Parent window
        df: DataFrame to assign columns for
        tab_type: Type of tab ('data_cleaning', 'statistics', 'visualization', 'pathway', 'id_annotation')
        auto_calculate: Whether to calculate missing columns
        existing_assignments: Previously confirmed assignments
        excel_file_path: Path to Excel file (for multi-sheet support)
        dialog_title: Optional custom dialog title (ignored - uses tab_type)
        group_area_pattern: Optional custom pattern for Group Area column detection
        detected_sample_cols: Optional pre-detected sample columns to use instead of auto-detecting
    
    Returns:
        Dict with 'assignments', 'dataframe', 'sample_cols', and 'selected_sheet' keys, or None if cancelled
    """
    # dialog_title parameter is accepted but not used (title is determined by tab_type)
    dialog = ColumnAssignmentDialog(
        parent,
        df,
        tab_type,
        auto_calculate,
        existing_assignments,
        excel_file_path,
        group_area_pattern,
        detected_sample_cols,
        allow_skip=allow_skip,
    )
    return dialog.wait_for_result()


__all__ = [
    'ColumnAssignmentDialog',
    'ColumnDetector',
    'ColumnCalculator',
    'ComparisonColumnMatcher',
    'show_column_assignment_dialog',
    'COLUMN_TYPES',
    'TAB_REQUIREMENTS',
    'COMMON_ID_COLUMNS',
    'FEATURE_COLUMNS_CANONICAL',
]

