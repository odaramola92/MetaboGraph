"""
Database Builder
Processes raw database files into optimized feather format for metabolite annotation.

This script converts:
1. hmdb_metabolites.xml → hmdb_database.feather + All_metabolites_synonyms_hmdb.feather
2. structures.sdf (LipidMaps) → lipidmap.feather
3. pathbank_all_metabolites.csv → pathbank_selected.feather (Human, Rat, Mouse only)
4. smpdb_metabolites.csv/*.csv → merged_SMP_metabolites.feather
5. WikiPathways GPML folders → wikipathways_homo_sapiens.feather (+ rat, mouse)

Usage:
    python database_builder.py --hmdb hmdb_metabolites.xml --lipidmaps structures.sdf
    python database_builder.py --pathbank pathbank_all_metabolites.csv
    python database_builder.py --smpdb smpdb_metabolites.csv
    python database_builder.py --wikipathways
    
Or run interactively:
    python database_builder.py
"""

import pandas as pd
import xml.etree.ElementTree as ET
import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Any
import logging
from datetime import datetime
import glob

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class HMDBDatabaseBuilder:
    """Build HMDB feather files from XML."""
    
    def __init__(self, xml_path: str):
        self.xml_path = xml_path
        self.namespace = {'hmdb': 'http://www.hmdb.ca'}
        
    def parse_hmdb_xml(self) -> pd.DataFrame:
        """
        Parse HMDB XML file and extract metabolite data.
        
        Returns:
            DataFrame with metabolite information
        """
        logger.info(f"Parsing HMDB XML: {self.xml_path}")
        
        if not os.path.exists(self.xml_path):
            raise FileNotFoundError(f"HMDB XML file not found: {self.xml_path}")
        
        metabolites = []
        
        try:
            # Parse XML iteratively to handle large files
            context = ET.iterparse(self.xml_path, events=('start', 'end'))
            context = iter(context)
            event, root = next(context)
            
            count = 0
            for event, elem in context:
                if event == 'end' and elem.tag == '{http://www.hmdb.ca}metabolite':
                    metabolite_data = self._extract_metabolite_data(elem)
                    if metabolite_data:
                        metabolites.append(metabolite_data)
                        count += 1
                        
                        if count % 1000 == 0:
                            logger.info(f"Processed {count} metabolites...")
                    
                    # Clear element to save memory
                    elem.clear()
                    root.clear()
            
            logger.info(f"✅ Successfully parsed {count} metabolites from HMDB XML")
            
        except Exception as e:
            logger.error(f"Error parsing HMDB XML: {e}")
            raise
        
        df = pd.DataFrame(metabolites)
        
        # Add computed columns
        if 'Name' in df.columns:
            df['Name_lc'] = df['Name'].str.lower()
        if 'Formula' in df.columns:
            df['Formula_norm'] = df['Formula'].str.replace(r'\s+', '', regex=True)
        
        return df
    
    def _extract_metabolite_data(self, metabolite_elem) -> Dict[str, Any]:
        """Extract data from a single metabolite XML element."""
        data = {}
        
        # Helper to safely get text
        def get_text(path: str, default='') -> str:
            elem = metabolite_elem.find(path, self.namespace)
            return elem.text.strip() if elem is not None and elem.text else default
        
        # Basic information
        data['HMDB'] = get_text('hmdb:accession')
        data['Name'] = get_text('hmdb:name')
        data['Formula'] = get_text('hmdb:chemical_formula')
        
        # Molecular properties
        try:
            data['Molecular_Weight'] = float(get_text('hmdb:average_molecular_weight') or 0)
        except:
            data['Molecular_Weight'] = None
        
        try:
            data['Monoisotopic_Mass'] = float(get_text('hmdb:monisotopic_molecular_weight') or 0)
        except:
            data['Monoisotopic_Mass'] = None
        
        # Taxonomy
        taxonomy = metabolite_elem.find('hmdb:taxonomy', self.namespace)
        if taxonomy is not None:
            data['Kingdom'] = get_text('hmdb:kingdom', default='')
            data['Direct_Parent'] = taxonomy.find('hmdb:direct_parent', self.namespace)
            data['Direct_Parent'] = data['Direct_Parent'].text if data['Direct_Parent'] is not None else ''
            data['Super_Class'] = taxonomy.find('hmdb:super_class', self.namespace)
            data['Super_Class'] = data['Super_Class'].text if data['Super_Class'] is not None else ''
            data['Class'] = taxonomy.find('hmdb:class', self.namespace)
            data['Class'] = data['Class'].text if data['Class'] is not None else ''
            data['Sub_Class'] = taxonomy.find('hmdb:sub_class', self.namespace)
            data['Sub_Class'] = data['Sub_Class'].text if data['Sub_Class'] is not None else ''
            data['Molecular_Framework'] = taxonomy.find('hmdb:molecular_framework', self.namespace)
            data['Molecular_Framework'] = data['Molecular_Framework'].text if data['Molecular_Framework'] is not None else ''
        
        # Synonyms (collect all)
        synonyms_elem = metabolite_elem.find('hmdb:synonyms', self.namespace)
        synonyms = []
        if synonyms_elem is not None:
            for syn in synonyms_elem.findall('hmdb:synonym', self.namespace):
                if syn.text:
                    synonyms.append(syn.text.strip())
        data['Synonyms'] = '|'.join(synonyms) if synonyms else ''
        
        # Get biological_properties (contains biospecimen, tissue, cellular, and pathways)
        bio_props = metabolite_elem.find('hmdb:biological_properties', self.namespace)
        
        # Biospecimen locations (inside biological_properties)
        if bio_props is not None:
            biospecimen = bio_props.find('hmdb:biospecimen_locations', self.namespace)
            if biospecimen is not None:
                locations = [loc.text.strip() for loc in biospecimen.findall('hmdb:biospecimen', self.namespace) if loc.text]
                data['In_Saliva'] = 'Yes' if 'Saliva' in locations else 'No'
                data['In_Blood'] = 'Yes' if 'Blood' in locations else 'No'
                data['In_Urine'] = 'Yes' if 'Urine' in locations else 'No'
                data['In_CSF'] = 'Yes' if 'Cerebrospinal Fluid (CSF)' in locations else 'No'
                data['In_Breast_Milk'] = 'Yes' if 'Breast Milk' in locations else 'No'
                data['In_Feces'] = 'Yes' if 'Feces' in locations else 'No'
                data['In_Sweat'] = 'Yes' if 'Sweat' in locations else 'No'
        
        # Tissue locations (inside biological_properties)
        if bio_props is not None:
            tissue_elem = bio_props.find('hmdb:tissue_locations', self.namespace)
            if tissue_elem is not None:
                tissues = [t.text.strip() for t in tissue_elem.findall('hmdb:tissue', self.namespace) if t.text]
                data['Tissues'] = '|'.join(tissues) if tissues else ''
        
        # Cellular locations (inside biological_properties)
        if bio_props is not None:
            cellular_elem = bio_props.find('hmdb:cellular_locations', self.namespace)
            if cellular_elem is not None:
                cellular = [c.text.strip() for c in cellular_elem.findall('hmdb:cellular', self.namespace) if c.text]
                data['Cellular_Locations'] = '|'.join(cellular) if cellular else ''
        
        # Pathways (inside biological_properties)
        pathways = []
        if bio_props is not None:
            pathways_elem = bio_props.find('hmdb:pathways', self.namespace)
            if pathways_elem is not None:
                for pw in pathways_elem.findall('hmdb:pathway', self.namespace):
                    pw_name = pw.find('hmdb:name', self.namespace)
                    if pw_name is not None and pw_name.text:
                        pathways.append(pw_name.text.strip())
        data['Pathways'] = '|'.join(pathways) if pathways else ''
        
        # Diseases
        diseases_elem = metabolite_elem.find('hmdb:diseases', self.namespace)
        if diseases_elem is not None:
            diseases = []
            for dis in diseases_elem.findall('hmdb:disease', self.namespace):
                dis_name = dis.find('hmdb:name', self.namespace)
                if dis_name is not None and dis_name.text:
                    diseases.append(dis_name.text.strip())
            data['Associated_Diseases'] = '|'.join(diseases) if diseases else ''
        
        # Protein associations
        protein_elem = metabolite_elem.find('hmdb:protein_associations', self.namespace)
        if protein_elem is not None:
            enzymes, enzyme_accessions, enzyme_genes = [], [], []
            transporters, transporter_genes, transporter_ids = [], [], []
            
            for protein in protein_elem.findall('hmdb:protein', self.namespace):
                p_type = protein.find('hmdb:protein_type', self.namespace)
                p_name = protein.find('hmdb:name', self.namespace)
                p_gene = protein.find('hmdb:gene_name', self.namespace)
                p_uniprot = protein.find('hmdb:uniprot_id', self.namespace)
                
                if p_type is not None and p_type.text:
                    if 'enzyme' in p_type.text.lower():
                        if p_name is not None and p_name.text:
                            enzymes.append(p_name.text.strip())
                        if p_uniprot is not None and p_uniprot.text:
                            enzyme_accessions.append(p_uniprot.text.strip())
                        if p_gene is not None and p_gene.text:
                            enzyme_genes.append(p_gene.text.strip())
                    elif 'transporter' in p_type.text.lower():
                        if p_name is not None and p_name.text:
                            transporters.append(p_name.text.strip())
                        if p_gene is not None and p_gene.text:
                            transporter_genes.append(p_gene.text.strip())
                        if p_uniprot is not None and p_uniprot.text:
                            transporter_ids.append(p_uniprot.text.strip())
            
            data['Enzymes'] = '|'.join(enzymes) if enzymes else ''
            data['Enzymes_Accession'] = '|'.join(enzyme_accessions) if enzyme_accessions else ''
            data['Enzyme_Gene_name'] = '|'.join(enzyme_genes) if enzyme_genes else ''
            data['Transporters'] = '|'.join(transporters) if transporters else ''
            data['Transporter_Gene_name'] = '|'.join(transporter_genes) if transporter_genes else ''
            data['Transporter_Uniprot_ID'] = '|'.join(transporter_ids) if transporter_ids else ''
        
        # Chemical identifiers
        data['IUPAC_Name'] = get_text('hmdb:iupac_name')
        data['CAS'] = get_text('hmdb:cas_registry_number')
        data['SMILES'] = get_text('hmdb:smiles')
        data['InChI'] = get_text('hmdb:inchi')
        data['InChIKey'] = get_text('hmdb:inchikey')
        
        # State
        data['State'] = get_text('hmdb:state')
        
        # Description
        data['Description'] = get_text('hmdb:description')
        
        # Endogenous/Exogenous status (from ontology)
        is_endogenous = False
        is_exogenous = False
        ontology = metabolite_elem.find('hmdb:ontology', self.namespace)
        if ontology is not None:
            # Search for Endogenous/Exogenous in ontology > Disposition > Source descendants
            for root_elem in ontology.findall('hmdb:root', self.namespace):
                term = root_elem.find('hmdb:term', self.namespace)
                if term is not None and term.text == 'Disposition':
                    descendants = root_elem.find('hmdb:descendants', self.namespace)
                    if descendants is not None:
                        for desc in descendants.findall('hmdb:descendant', self.namespace):
                            desc_term = desc.find('hmdb:term', self.namespace)
                            if desc_term is not None and desc_term.text == 'Source':
                                # Check its descendants for Endogenous/Exogenous
                                desc_descendants = desc.find('hmdb:descendants', self.namespace)
                                if desc_descendants is not None:
                                    for sub_desc in desc_descendants.findall('hmdb:descendant', self.namespace):
                                        sub_term = sub_desc.find('hmdb:term', self.namespace)
                                        if sub_term is not None:
                                            if sub_term.text == 'Endogenous':
                                                is_endogenous = True
                                            elif sub_term.text == 'Exogenous':
                                                is_exogenous = True
        
        data['Is_Endogenous'] = 'Yes' if is_endogenous else 'No'
        data['Is_Exogenous'] = 'Yes' if is_exogenous else 'No'
        
        return data
    
    def build_synonyms_table(self, hmdb_df: pd.DataFrame) -> pd.DataFrame:
        """
        Build expanded synonyms table from HMDB data.
        
        Args:
            hmdb_df: Main HMDB DataFrame
            
        Returns:
            DataFrame with HMDB ID and individual synonyms
        """
        logger.info("Building synonyms table...")
        
        synonym_rows = []
        
        for _, row in hmdb_df.iterrows():
            hmdb_id = row.get('HMDB', '')
            synonyms_str = row.get('Synonyms', '')
            
            if not hmdb_id or not synonyms_str:
                continue
            
            # Split synonyms and create individual rows
            synonyms = synonyms_str.split('|')
            for syn in synonyms:
                syn = syn.strip()
                if syn:
                    synonym_rows.append({
                        'HMDB': hmdb_id,
                        'Synonyms': syn
                    })
        
        synonyms_df = pd.DataFrame(synonym_rows)
        logger.info(f"✅ Created synonyms table with {len(synonyms_df)} entries")
        
        return synonyms_df
    
    def build(self, output_dir: str = '.') -> tuple:
        """
        Build both HMDB feather files.
        
        Args:
            output_dir: Directory to save feather files
            
        Returns:
            Tuple of (main_db_path, synonyms_path)
        """
        # Parse XML
        hmdb_df = self.parse_hmdb_xml()
        
        # Build synonyms table
        synonyms_df = self.build_synonyms_table(hmdb_df)
        
        # Save to feather
        main_path = os.path.join(output_dir, 'hmdb_database.feather')
        synonyms_path = os.path.join(output_dir, 'All_metabolites_synonyms_hmdb.feather')
        
        logger.info(f"Saving main database to: {main_path}")
        hmdb_df.to_feather(main_path)
        
        logger.info(f"Saving synonyms table to: {synonyms_path}")
        synonyms_df.to_feather(synonyms_path)
        
        logger.info(f"✅ HMDB databases created successfully!")
        logger.info(f"   • Main DB: {len(hmdb_df)} metabolites, {main_path}")
        logger.info(f"   • Synonyms: {len(synonyms_df)} entries, {synonyms_path}")
        
        return main_path, synonyms_path

class LipidMapsDatabaseBuilder:
    """Build LipidMaps feather file from SDF."""
    
    def __init__(self, sdf_path: str):
        self.sdf_path = sdf_path
    
    def parse_sdf(self) -> pd.DataFrame:
        """
        Parse LipidMaps SDF file.
        
        Returns:
            DataFrame with lipid information
        """
        logger.info(f"Parsing LipidMaps SDF: {self.sdf_path}")
        
        if not os.path.exists(self.sdf_path):
            raise FileNotFoundError(f"SDF file not found: {self.sdf_path}")
        
        lipids = []
        
        try:
            with open(self.sdf_path, 'r', encoding='utf-8', errors='ignore') as f:
                current_lipid = {}
                molblock_lines = []
                in_properties = False
                count = 0
                
                for line in f:
                    line = line.rstrip('\n')
                    
                    # Start of properties section
                    if line.startswith('> <'):
                        in_properties = True
                        # Extract property name
                        prop_name = line.strip('> <>').strip()
                        # Read next line for value
                        continue
                    
                    # End of molecule
                    if line.strip() == '$$$$':
                        if current_lipid:
                            lipids.append(current_lipid)
                            count += 1
                            
                            if count % 1000 == 0:
                                logger.info(f"Processed {count} lipids...")
                        
                        current_lipid = {}
                        molblock_lines = []
                        in_properties = False
                        continue
                    
                    # Property value
                    if in_properties:
                        # The line after property name is the value
                        if line and not line.startswith('> <'):
                            # Find the last property name added
                            if '> <' in line:
                                continue
                            # This is a value line
                            # Need to track which property we're reading
                            pass
                
                # Alternative: use RDKit if available
                try:
                    from rdkit import Chem
                    from rdkit.Chem import Descriptors, AllChem
                    
                    logger.info("Using RDKit for SDF parsing (more reliable)")
                    return self._parse_sdf_with_rdkit()
                    
                except ImportError:
                    logger.warning("RDKit not available, using basic parsing")
                    return self._parse_sdf_basic()
        
        except Exception as e:
            logger.error(f"Error parsing SDF: {e}")
            raise
    
    def _parse_sdf_basic(self) -> pd.DataFrame:
        """Basic SDF parser without RDKit."""
        logger.info("Using basic SDF parser (properties may be limited)")
        
        lipids = []
        
        with open(self.sdf_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            molecules = content.split('$$$$\n')
            
            for i, mol_text in enumerate(molecules):
                if not mol_text.strip():
                    continue
                
                lipid_data = {}
                lines = mol_text.split('\n')
                
                # Try to extract LM_ID from first line
                if lines:
                    lipid_data['MoleculeName'] = lines[0].strip()
                
                # Parse properties
                in_prop = False
                current_prop = None
                
                for line in lines:
                    if line.startswith('> <'):
                        current_prop = line.strip('> <>').strip()
                        in_prop = True
                    elif in_prop and line.strip() and not line.startswith('> <'):
                        lipid_data[current_prop] = line.strip()
                        in_prop = False
                
                if lipid_data:
                    lipids.append(lipid_data)
                
                if (i + 1) % 1000 == 0:
                    logger.info(f"Processed {i + 1} molecules...")
        
        df = pd.DataFrame(lipids)
        logger.info(f"✅ Parsed {len(df)} lipids from SDF")
        
        return df
    
    def _parse_sdf_with_rdkit(self) -> pd.DataFrame:
        """Parse SDF using RDKit for better reliability."""
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        
        lipids = []
        supplier = Chem.SDMolSupplier(self.sdf_path)
        
        count = 0
        for mol in supplier:
            if mol is None:
                continue
            
            lipid_data = {}
            
            # Get all properties
            for prop_name in mol.GetPropNames():
                lipid_data[prop_name] = mol.GetProp(prop_name)
            
            # Add calculated properties if not present
            if 'EXACT_MASS' not in lipid_data:
                lipid_data['EXACT_MASS'] = Descriptors.ExactMolWt(mol)
            
            if 'FORMULA' not in lipid_data:
                lipid_data['FORMULA'] = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
            
            # Generate SMILES if not present
            if 'SMILES' not in lipid_data:
                lipid_data['SMILES'] = Chem.MolToSmiles(mol)
            
            if 'canonical_smiles' not in lipid_data:
                lipid_data['canonical_smiles'] = Chem.MolToSmiles(mol, canonical=True)
            
            lipids.append(lipid_data)
            count += 1
            
            if count % 1000 == 0:
                logger.info(f"Processed {count} lipids...")
        
        df = pd.DataFrame(lipids)
        logger.info(f"✅ Parsed {len(df)} lipids using RDKit")
        
        return df
    
    def standardize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Standardize column names to match expected format."""
        
        # Column name mappings
        column_map = {
            'LM_ID': 'LM_ID',
            'LMID': 'LM_ID',
            'LIPIDMAPS_ID': 'LM_ID',
            'COMMON_NAME': 'NAME',
            'SYNONYMS': 'NAME',
            'PUBCHEM_CID': 'PUBCHEM_CID',
            'CHEBI_ID': 'CHEBI_ID',
        }
        
        # Rename columns
        for old_name, new_name in column_map.items():
            if old_name in df.columns and new_name not in df.columns:
                df = df.rename(columns={old_name: new_name})
        
        # Ensure required columns exist
        required_cols = ['LM_ID', 'NAME', 'CATEGORY', 'MAIN_CLASS', 'SUB_CLASS', 
                        'EXACT_MASS', 'FORMULA', 'INCHI_KEY', 'SMILES']
        
        for col in required_cols:
            if col not in df.columns:
                df[col] = ''
        
        return df
    
    def build(self, output_dir: str = '.') -> str:
        """
        Build LipidMaps feather file.
        
        Args:
            output_dir: Directory to save feather file
            
        Returns:
            Path to created feather file
        """
        # Parse SDF
        lipid_df = self.parse_sdf()
        
        # Standardize columns
        lipid_df = self.standardize_columns(lipid_df)
        
        # Save to feather
        output_path = os.path.join(output_dir, 'lipidmap.feather')
        
        logger.info(f"Saving LipidMaps database to: {output_path}")
        lipid_df.to_feather(output_path)
        
        logger.info(f"✅ LipidMaps database created successfully!")
        logger.info(f"   • {len(lipid_df)} lipids, {output_path}")
        
        return output_path

class PathBankDatabaseBuilder:
    """Build PathBank feather file from CSV, filtered for Human, Rat, and Mouse."""
    
    ALLOWED_SPECIES = ['Homo sapiens', 'Mus musculus', 'Rattus norvegicus']
    
    def __init__(self, csv_path: str):
        self.csv_path = csv_path
    
    def parse_and_filter_csv(self) -> pd.DataFrame:
        """
        Parse PathBank CSV and filter for relevant species.
        
        Returns:
            DataFrame with pathway-metabolite associations for Human, Rat, Mouse
        """
        logger.info(f"Parsing PathBank CSV: {self.csv_path}")
        
        if not os.path.exists(self.csv_path):
            raise FileNotFoundError(f"PathBank CSV file not found: {self.csv_path}")
        
        try:
            # Read CSV in chunks for efficiency
            logger.info("Reading PathBank data...")
            df = pd.read_csv(self.csv_path, dtype=str, low_memory=False)
            
            total_rows = len(df)
            logger.info(f"Total rows in file: {total_rows:,}")
            
            # Show species distribution before filtering
            if 'Species' in df.columns:
                logger.info("Species distribution (before filtering):")
                species_counts = df['Species'].value_counts()
                for species, count in species_counts.head(10).items():
                    logger.info(f"  • {species}: {count:,}")
            
            # Filter for allowed species
            logger.info(f"\nFiltering for: {', '.join(self.ALLOWED_SPECIES)}")
            df_filtered = df[df['Species'].isin(self.ALLOWED_SPECIES)].copy()
            
            filtered_rows = len(df_filtered)
            logger.info(f"Rows after filtering: {filtered_rows:,} ({filtered_rows/total_rows*100:.1f}%)")
            
            # Show species distribution after filtering
            logger.info("\nSpecies distribution (after filtering):")
            species_counts_filtered = df_filtered['Species'].value_counts()
            for species, count in species_counts_filtered.items():
                logger.info(f"  • {species}: {count:,}")
            
            # Ensure all required columns are present
            required_columns = [
                'PathBank ID', 'Pathway Name', 'Pathway Subject', 'Species',
                'Metabolite ID', 'Metabolite Name', 'HMDB ID', 'KEGG ID',
                'ChEBI ID', 'DrugBank ID', 'CAS', 'Formula', 'IUPAC',
                'SMILES', 'InChI', 'InChI Key'
            ]
            
            missing_cols = [col for col in required_columns if col not in df_filtered.columns]
            if missing_cols:
                logger.warning(f"Missing columns: {missing_cols}")
            
            logger.info(f"\n✅ Successfully parsed and filtered PathBank data")
            
            return df_filtered
            
        except Exception as e:
            logger.error(f"Error parsing PathBank CSV: {e}")
            raise
    
    def build(self, output_dir: str = '.') -> str:
        """
        Build PathBank feather file.
        
        Args:
            output_dir: Directory to save feather file
            
        Returns:
            Path to created feather file
        """
        # Parse and filter CSV
        pathbank_df = self.parse_and_filter_csv()
        
        # Save to feather
        output_path = os.path.join(output_dir, 'pathbank_selected.feather')
        
        logger.info(f"\nSaving PathBank database to: {output_path}")
        pathbank_df.to_feather(output_path)
        
        logger.info(f"✅ PathBank database created successfully!")
        logger.info(f"   • {len(pathbank_df):,} pathway-metabolite associations")
        logger.info(f"   • Species: {', '.join(self.ALLOWED_SPECIES)}")
        logger.info(f"   • Output: {output_path}")
        
        return output_path

class WikiPathwaysDatabaseBuilder:
    """Build WikiPathways organism-specific feather files from GPML archives."""
    
    # Organism configurations
    ORGANISMS = {
        "Homo sapiens": {
            "folder_pattern": "wikipathways-*-gpml-Homo_sapiens",
            "output_file": "wikipathways_homo_sapiens.feather"
        },
        "Rattus norvegicus": {
            "folder_pattern": "wikipathways-*-gpml-Rattus_norvegicus",
            "output_file": "wikipathways_rattus_norvegicus.feather"
        },
        "Mus musculus": {
            "folder_pattern": "wikipathways-*-gpml-Mus_musculus",
            "output_file": "wikipathways_mus_musculus.feather"
        }
    }
    
    def __init__(self, base_folder: str = "."):
        self.base_folder = base_folder
        self.namespace = {"gpml": "http://pathvisio.org/GPML/2013a"}
    
    def find_gpml_folder(self, pattern: str) -> str:
        """Find GPML folder matching the pattern."""
        matches = glob.glob(os.path.join(self.base_folder, pattern))
        if matches:
            return matches[0]
        return None
    
    def parse_gpml_files(self, gpml_folder: str, organism_name: str) -> pd.DataFrame:
        """
        Parse all GPML files in a folder and extract metabolite information.
        
        Args:
            gpml_folder: Path to folder containing GPML files
            organism_name: Name of the organism (e.g., "Homo sapiens")
        
        Returns:
            DataFrame with metabolite data and pathway annotations
        """
        if not os.path.exists(gpml_folder):
            logger.error(f"Folder not found: {gpml_folder}")
            return None
        
        logger.info(f"Processing {organism_name} from: {gpml_folder}")
        records = []
        
        # Get all GPML files
        gpml_files = [f for f in os.listdir(gpml_folder) if f.endswith(".gpml")]
        logger.info(f"Found {len(gpml_files)} GPML files")
        
        for file in gpml_files:
            file_path = os.path.join(gpml_folder, file)
            try:
                tree = ET.parse(file_path)
                root = tree.getroot()

                pathway_name = root.attrib.get("Name", "")
                pathway_id = root.attrib.get("Identifier", "")
                version = root.attrib.get("Version", "")
                organism = root.attrib.get("Organism", organism_name)
                data_source = root.attrib.get("Data-Source", "WikiPathways")

                # Extract DataNode elements (focus on metabolites)
                for node in root.findall(".//gpml:DataNode", self.namespace):
                    node_type = node.attrib.get("Type", "")
                    
                    if node_type.lower() == "metabolite":
                        xref = node.find("gpml:Xref", self.namespace)
                        xref_db = xref.attrib.get("Database", "") if xref is not None else ""
                        xref_id = xref.attrib.get("ID", "") if xref is not None else ""

                        records.append({
                            "TextLabel": node.attrib.get("TextLabel", ""),
                            "Type": node_type,
                            "Xref_Database": xref_db,
                            "Xref_ID": xref_id,
                            "Pathway_Name": pathway_name,
                            "Pathway_ID": pathway_id,
                            "Version": version,
                            "Organism": organism,
                            "Data_Source": data_source,
                            "Source_File": file
                        })
            except Exception as e:
                logger.warning(f"Error parsing {file}: {e}")
        
        if not records:
            logger.warning(f"No metabolite DataNode elements found")
            return None
        
        logger.info(f"Extracted {len(records):,} metabolite records from {len(gpml_files)} pathways")
        
        # Create DataFrame
        df = pd.DataFrame(records)
        
        # Group by metabolite name + Xref to collect all pathways
        logger.info("Grouping by metabolite and database identifiers...")
        grouped = (
            df.groupby(["TextLabel", "Xref_Database", "Xref_ID", "Organism"], dropna=False)
              .agg({
                  "Pathway_Name": lambda x: ", ".join(sorted(set(filter(None, x)))),
                  "Pathway_ID": lambda x: ", ".join(sorted(set(filter(None, x)))),
                  "Version": "first",
                  "Data_Source": "first",
                  "Type": "first"
              })
              .reset_index()
        )
        
        # Map each Xref_Database into its own column
        logger.info("Organizing database identifiers...")
        id_cols = ["TextLabel", "Xref_Database", "Xref_ID", "Pathway_Name", "Pathway_ID", "Organism"]
        
        final_grouped = (
            grouped[id_cols]
            .groupby(["Organism", "TextLabel"], dropna=False)
            .apply(lambda g: {
                db: ", ".join(sorted(set(g.loc[g["Xref_Database"] == db, "Xref_ID"])))
                for db in g["Xref_Database"].unique() if db
            } | {
                "Pathways": ", ".join(sorted(set(g["Pathway_Name"]))),
                "Pathway_IDs": ", ".join(sorted(set(g["Pathway_ID"])))
            })
            .reset_index()
        )
        
        # Expand dicts into columns
        expanded = final_grouped["TextLabel"].to_frame()
        expanded.columns = ["Metabolite"]
        expanded["Organism"] = final_grouped["Organism"]
        
        meta_df = pd.json_normalize(final_grouped[0])
        final_df = pd.concat([expanded, meta_df], axis=1)
        
        # Rename columns to match expected output
        rename_map = {
            "HMDB": "HMDB_ID",
            "PubChem": "PubChem",
            "KEGG": "KEGG",
            "InChIKey": "InChIKey",
            "ChEBI": "ChEBI",
            "Entrez Gene": "Entrez",
            "Wikidata": "Wikidata",
            "ChemSpider": "ChemSpider",
            "CAS": "CAS"
        }
        final_df = final_df.rename(columns=rename_map)
        
        # Ensure all required columns exist
        required_cols = ["HMDB_ID", "PubChem", "KEGG", "InChIKey", "ChEBI", "Entrez", 
                        "Wikidata", "ChemSpider", "CAS"]
        for col in required_cols:
            if col not in final_df.columns:
                final_df[col] = ""
        
        # Add 'name' column (same as Metabolite for compatibility)
        final_df["name"] = final_df["Metabolite"]
        
        # Reorder columns
        column_order = [
            "Organism", "name", "Metabolite", "HMDB_ID", "PubChem", "KEGG", 
            "InChIKey", "ChEBI", "Entrez", "Pathways", "Pathway_IDs", 
            "Wikidata", "ChemSpider", "CAS"
        ]
        
        final_columns = [col for col in column_order if col in final_df.columns]
        final_df = final_df[final_columns]
        
        logger.info(f"Created dataframe with {len(final_df):,} unique metabolites")
        
        return final_df
    
    def build(self, output_dir: str = ".") -> List[str]:
        """
        Build WikiPathways databases for all available organisms.
        
        Args:
            output_dir: Directory to save output files
        
        Returns:
            List of successfully created output files
        """
        logger.info("Building WikiPathways databases...")
        
        # Check which organism folders exist
        found_folders = []
        missing_folders = []
        
        for organism, config in self.ORGANISMS.items():
            folder = self.find_gpml_folder(config["folder_pattern"])
            if folder:
                found_folders.append((organism, folder, config["output_file"]))
            else:
                missing_folders.append((organism, config["folder_pattern"]))
        
        if missing_folders:
            logger.warning("Missing WikiPathways data for:")
            for organism, pattern in missing_folders:
                logger.warning(f"  - {organism} (looking for: {pattern})")
            logger.info("\n📥 Download from: https://zenodo.org/communities/wikipathways/records")
            logger.info("   Extract GPML archives to the working directory\n")
        
        if not found_folders:
            logger.error("No WikiPathways GPML folders found!")
            return []
        
        logger.info(f"Found {len(found_folders)} organism folder(s)")
        
        created_files = []
        
        for organism, folder, output_file in found_folders:
            try:
                df = self.parse_gpml_files(folder, organism)
                
                if df is not None and len(df) > 0:
                    output_path = os.path.join(output_dir, output_file)
                    df.to_feather(output_path)
                    logger.info(f"✅ Saved {organism} database: {output_path}")
                    logger.info(f"   {len(df):,} metabolites, {df.columns.tolist()}")
                    created_files.append(output_path)
            except Exception as e:
                logger.error(f"Error processing {organism}: {e}")
        
        return created_files


class SMPDBDatabaseBuilder:
    """Build SMPDB feather file by merging individual pathway CSV files."""
    
    def __init__(self, csv_folder: str):
        self.csv_folder = csv_folder
    
    def merge_csv_files(self) -> pd.DataFrame:
        """
        Merge all SMPDB pathway CSV files into a single DataFrame.
        
        Returns:
            DataFrame with all metabolites and their associated pathways
        """
        logger.info(f"Merging SMPDB CSV files from: {self.csv_folder}")
        
        if not os.path.exists(self.csv_folder):
            raise FileNotFoundError(f"SMPDB folder not found: {self.csv_folder}")
        
        # Get all CSV files
        csv_pattern = os.path.join(self.csv_folder, "SMP*.csv")
        csv_files = glob.glob(csv_pattern)
        
        if not csv_files:
            raise ValueError(f"No SMPDB CSV files found in: {self.csv_folder}")
        
        logger.info(f"Found {len(csv_files):,} SMPDB pathway files")
        
        # Dictionary to store metabolites with their pathways
        metabolite_pathways = {}
        
        try:
            for idx, csv_file in enumerate(csv_files, 1):
                try:
                    # Read pathway CSV
                    df = pd.read_csv(csv_file, dtype=str)
                    
                    # Get pathway ID from filename (e.g., SMP0000001)
                    pathway_id = os.path.basename(csv_file).replace('_metabolites.csv', '')
                    
                    # Get pathway name if available
                    pathway_name = df['Pathway Name'].iloc[0] if 'Pathway Name' in df.columns and len(df) > 0 else pathway_id
                    
                    # Process each metabolite in this pathway
                    for _, row in df.iterrows():
                        # Use metabolite name as key
                        met_name = row.get('Metabolite Name', '').strip()
                        if not met_name:
                            continue
                        
                        # Initialize metabolite entry if not exists
                        if met_name not in metabolite_pathways:
                            metabolite_pathways[met_name] = {
                                'Name': met_name,
                                'SMP_Pathways': [],
                                'HMDB': row.get('HMDB ID', ''),
                                'KEGG': row.get('KEGG ID', ''),
                                'CAS': row.get('CAS', ''),
                                'Formula': row.get('Formula', ''),
                                'IUPAC': row.get('IUPAC', ''),
                                'SMILES': row.get('SMILES', ''),
                                'InChI': row.get('InChI', ''),
                                'InChIKey': row.get('InChI Key', '')
                            }
                        
                        # Add pathway to this metabolite
                        metabolite_pathways[met_name]['SMP_Pathways'].append(f"{pathway_id}:{pathway_name}")
                        
                        # Update IDs if they were empty before
                        for id_field in ['HMDB', 'KEGG', 'CAS', 'Formula', 'IUPAC', 'SMILES', 'InChI', 'InChIKey']:
                            csv_col_map = {
                                'HMDB': 'HMDB ID',
                                'KEGG': 'KEGG ID',
                                'InChIKey': 'InChI Key'
                            }
                            csv_col = csv_col_map.get(id_field, id_field)
                            
                            if not metabolite_pathways[met_name][id_field] and csv_col in row:
                                value = row.get(csv_col, '').strip()
                                if value and value.lower() not in ['nan', 'n/a', '']:
                                    metabolite_pathways[met_name][id_field] = value
                    
                    if (idx % 1000 == 0) or (idx == len(csv_files)):
                        logger.info(f"Processed {idx:,}/{len(csv_files):,} files ({idx/len(csv_files)*100:.1f}%) - {len(metabolite_pathways):,} unique metabolites")
                
                except Exception as e:
                    logger.warning(f"Error processing {csv_file}: {e}")
                    continue
            
            # Convert to DataFrame
            logger.info("\nConverting to DataFrame...")
            records = []
            for met_name, data in metabolite_pathways.items():
                # Join pathways with pipe separator
                data['SMP_Pathways'] = '|'.join(data['SMP_Pathways'])
                records.append(data)
            
            df = pd.DataFrame(records)
            
            # Sort by metabolite name
            df = df.sort_values('Name').reset_index(drop=True)
            
            logger.info(f"✅ Successfully merged {len(csv_files):,} SMPDB files")
            logger.info(f"   • Total unique metabolites: {len(df):,}")
            logger.info(f"   • Metabolites with HMDB IDs: {df['HMDB'].notna().sum():,}")
            logger.info(f"   • Metabolites with KEGG IDs: {df['KEGG'].notna().sum():,}")
            
            return df
            
        except Exception as e:
            logger.error(f"Error merging SMPDB files: {e}")
            raise
    
    def build(self, output_dir: str = '.') -> str:
        """
        Build SMPDB merged feather file.
        
        Args:
            output_dir: Directory to save feather file
            
        Returns:
            Path to created feather file
        """
        # Merge all CSV files
        smpdb_df = self.merge_csv_files()
        
        # Save to feather
        output_path = os.path.join(output_dir, 'merged_SMP_metabolites.feather')
        
        logger.info(f"\nSaving SMPDB database to: {output_path}")
        smpdb_df.to_feather(output_path)
        
        logger.info(f"✅ SMPDB database created successfully!")
        logger.info(f"   • {len(smpdb_df):,} unique metabolites")
        logger.info(f"   • Output: {output_path}")
        
        return output_path

def interactive_build():
    """Interactive mode for building databases."""
    print("\n" + "="*80)
    print(" METABOLITE DATABASE BUILDER".center(80))
    print("="*80 + "\n")
    
    print("This tool will help you convert raw database files to feather format.\n")
    
    # HMDB
    print("─" * 80)
    print("📊 HMDB Database")
    print("─" * 80)
    hmdb_choice = input("\nDo you want to build HMDB databases? (y/n): ").strip().lower()
    
    if hmdb_choice == 'y':
        hmdb_xml = input("Enter path to HMDB XML file (e.g., hmdb_metabolites.xml): ").strip()
        
        if os.path.exists(hmdb_xml):
            try:
                builder = HMDBDatabaseBuilder(hmdb_xml)
                builder.build()
                print("\n✅ HMDB databases created successfully!\n")
            except Exception as e:
                print(f"\n❌ Error building HMDB databases: {e}\n")
        else:
            print(f"\n❌ File not found: {hmdb_xml}\n")
    
    # LipidMaps
    print("─" * 80)
    print("🧬 LipidMaps Database")
    print("─" * 80)
    lipid_choice = input("\nDo you want to build LipidMaps database? (y/n): ").strip().lower()
    
    if lipid_choice == 'y':
        lipid_sdf = input("Enter path to LipidMaps SDF file (e.g., structures.sdf): ").strip()
        
        if os.path.exists(lipid_sdf):
            try:
                builder = LipidMapsDatabaseBuilder(lipid_sdf)
                builder.build()
                print("\n✅ LipidMaps database created successfully!\n")
            except Exception as e:
                print(f"\n❌ Error building LipidMaps database: {e}\n")
        else:
            print(f"\n❌ File not found: {lipid_sdf}\n")
    
    # PathBank
    print("─" * 80)
    print("🧬 PathBank Database")
    print("─" * 80)
    print("Note: This will filter for Human, Rat, and Mouse pathways only.")
    pathbank_choice = input("\nDo you want to build PathBank database? (y/n): ").strip().lower()
    
    if pathbank_choice == 'y':
        pathbank_csv = input("Enter path to PathBank CSV file (e.g., pathbank_all_metabolites.csv): ").strip()
        
        if os.path.exists(pathbank_csv):
            try:
                builder = PathBankDatabaseBuilder(pathbank_csv)
                builder.build()
                print("\n✅ PathBank database created successfully!\n")
            except Exception as e:
                print(f"\n❌ Error building PathBank database: {e}\n")
        else:
            print(f"\n❌ File not found: {pathbank_csv}\n")
    
    # SMPDB
    print("─" * 80)
    print("🔬 SMPDB Database")
    print("─" * 80)
    print("Note: This will merge all individual pathway CSV files.")
    smpdb_choice = input("\nDo you want to build SMPDB database? (y/n): ").strip().lower()
    
    if smpdb_choice == 'y':
        smpdb_folder = input("Enter path to SMPDB CSV folder (e.g., smpdb_metabolites.csv): ").strip()
        
        if os.path.exists(smpdb_folder) and os.path.isdir(smpdb_folder):
            try:
                builder = SMPDBDatabaseBuilder(smpdb_folder)
                builder.build()
                print("\n✅ SMPDB database created successfully!\n")
            except Exception as e:
                print(f"\n❌ Error building SMPDB database: {e}\n")
        else:
            print(f"\n❌ Folder not found or not a directory: {smpdb_folder}\n")
    
    # WikiPathways
    print("─" * 80)
    print("🌐 WikiPathways Database (Organism-Specific)")
    print("─" * 80)
    print("Note: This will create separate databases for Human, Rat, and Mouse.")
    print("Download GPML archives from: https://zenodo.org/communities/wikipathways/records")
    wikipathways_choice = input("\nDo you want to build WikiPathways databases? (y/n): ").strip().lower()
    
    if wikipathways_choice == 'y':
        try:
            builder = WikiPathwaysDatabaseBuilder()
            created_files = builder.build()
            if created_files:
                print(f"\n✅ WikiPathways databases created: {len(created_files)} organism(s)\n")
                for file in created_files:
                    print(f"   • {file}")
            else:
                print("\n⚠️ No WikiPathways data found. Please download GPML archives.\n")
        except Exception as e:
            print(f"\n❌ Error building WikiPathways databases: {e}\n")
    
    print("="*80)
    print("\n✨ Database building complete! You can now use these files for annotation.\n")

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Build metabolite database feather files from raw XML/SDF/CSV sources',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Build all databases
  python database_builder.py --hmdb hmdb_metabolites.xml --lipidmaps structures.sdf \\
                             --pathbank pathbank_all_metabolites.csv \\
                             --smpdb smpdb_metabolites.csv \\
                             --wikipathways
  
  # Build only PathBank, SMPDB, and WikiPathways
  python database_builder.py --pathbank pathbank_all_metabolites.csv \\
                             --smpdb smpdb_metabolites.csv \\
                             --wikipathways
  
  # Build only WikiPathways (auto-detects GPML folders)
  python database_builder.py --wikipathways
  
  # Interactive mode
  python database_builder.py --interactive
        """
    )
    parser.add_argument('--hmdb', type=str, help='Path to HMDB XML file')
    parser.add_argument('--lipidmaps', type=str, help='Path to LipidMaps SDF file')
    parser.add_argument('--pathbank', type=str, help='Path to PathBank CSV file (will filter for Human/Rat/Mouse)')
    parser.add_argument('--smpdb', type=str, help='Path to SMPDB CSV folder containing individual pathway files')
    parser.add_argument('--wikipathways', action='store_true', help='Build WikiPathways databases (auto-detects GPML folders)')
    parser.add_argument('--output', type=str, default='.', help='Output directory for feather files')
    parser.add_argument('--interactive', action='store_true', help='Run in interactive mode')
    
    args = parser.parse_args()
    
    # Interactive mode if no arguments provided
    if args.interactive or (not args.hmdb and not args.lipidmaps and not args.pathbank and not args.smpdb and not args.wikipathways):
        interactive_build()
        return
    
    # Track successes and failures
    built_databases = []
    failed_databases = []
    
    # Build HMDB if provided
    if args.hmdb:
        try:
            logger.info("\n" + "="*80)
            logger.info(" Building HMDB Database")
            logger.info("="*80)
            builder = HMDBDatabaseBuilder(args.hmdb)
            builder.build(args.output)
            built_databases.append("HMDB")
        except Exception as e:
            logger.error(f"Failed to build HMDB databases: {e}")
            failed_databases.append(("HMDB", str(e)))
    
    # Build LipidMaps if provided
    if args.lipidmaps:
        try:
            logger.info("\n" + "="*80)
            logger.info(" Building LipidMaps Database")
            logger.info("="*80)
            builder = LipidMapsDatabaseBuilder(args.lipidmaps)
            builder.build(args.output)
            built_databases.append("LipidMaps")
        except Exception as e:
            logger.error(f"Failed to build LipidMaps database: {e}")
            failed_databases.append(("LipidMaps", str(e)))
    
    # Build PathBank if provided
    if args.pathbank:
        try:
            logger.info("\n" + "="*80)
            logger.info(" Building PathBank Database (Human/Rat/Mouse only)")
            logger.info("="*80)
            builder = PathBankDatabaseBuilder(args.pathbank)
            builder.build(args.output)
            built_databases.append("PathBank")
        except Exception as e:
            logger.error(f"Failed to build PathBank database: {e}")
            failed_databases.append(("PathBank", str(e)))
    
    # Build SMPDB if provided
    if args.smpdb:
        try:
            logger.info("\n" + "="*80)
            logger.info(" Building SMPDB Database")
            logger.info("="*80)
            builder = SMPDBDatabaseBuilder(args.smpdb)
            builder.build(args.output)
            built_databases.append("SMPDB")
        except Exception as e:
            logger.error(f"Failed to build SMPDB database: {e}")
            failed_databases.append(("SMPDB", str(e)))
    
    # Build WikiPathways if requested
    if args.wikipathways:
        try:
            logger.info("\n" + "="*80)
            logger.info(" Building WikiPathways Databases (Organism-Specific)")
            logger.info("="*80)
            builder = WikiPathwaysDatabaseBuilder(args.output)
            created_files = builder.build(args.output)
            if created_files:
                built_databases.append(f"WikiPathways ({len(created_files)} organisms)")
            else:
                logger.warning("No WikiPathways GPML folders found")
        except Exception as e:
            logger.error(f"Failed to build WikiPathways databases: {e}")
            failed_databases.append(("WikiPathways", str(e)))
    
    # Print summary
    print("\n" + "="*80)
    print(" BUILD SUMMARY".center(80))
    print("="*80)
    
    if built_databases:
        print(f"\n✅ Successfully built {len(built_databases)} database(s):")
        for db_name in built_databases:
            print(f"   • {db_name}")
    
    if failed_databases:
        print(f"\n❌ Failed to build {len(failed_databases)} database(s):")
        for db_name, error in failed_databases:
            print(f"   • {db_name}: {error}")
        sys.exit(1)
    
    if not built_databases and not failed_databases:
        print("\nNo databases were specified. Use --help for usage information.")
        sys.exit(1)
    
    print("\n" + "="*80 + "\n")
    logger.info("✅ All requested databases built successfully!")

if __name__ == "__main__":
    main()
