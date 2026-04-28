#!/usr/bin/env python3
"""
Data Cleaner for Metabolite Analysis - Clean Column Format
Extracts and cleans metabolite data with Metabolika Pathways and BioCyc Pathways columns
"""

import pandas as pd
import numpy as np
import os
import sys
import logging
import re
from collections import Counter
from pathlib import Path
from typing import Optional, Dict, Any, Callable, List


class CleanColumnDataCleaner:
    """Data cleaner for the clean column format with pathways"""

    def __init__(self, progress_callback: Optional[Callable] = None):
        """
        Initialize the data cleaner

        Args:
            progress_callback: Optional callback function to report progress messages
        """
        self.progress_callback = progress_callback or self._default_progress_callback
        
        # Set up logging
        logging.basicConfig(level=logging.INFO)
        
        # Combined reference ion list for user selection (all 10 ions)
        # This is the master list shown in the GUI with all ions checked by default
        self.REFERENCE_IONS_ALL = [
            '[M+H]+1',      # Positive
            '[M-H]-1',      # Negative
            '[M+2H]+2',     # Positive
            '[M-2H]-2',     # Negative
            '[M+H-H2O]+1',  # Positive
            '[M-H-H2O]-1',  # Negative
            '[M+H+MeOH]+1', # Positive
            '[M+FA+H]+1',   # Positive
            '[M-H-MeOH]-1', # Negative
            '[M+FA-H]-1'    # Negative
        ]
        
        # Ion order definitions (explicit charges) – user-specified
        # Positive mode (default priority order)
        self.positive_ion_order = [
            '[M+H]+1',
            '[M+2H]+2',
            '[M+H-H2O]+1',
            '[M+H+MeOH]+1',
            '[M+FA+H]+1'
        ]
        # Negative mode (default priority order)
        self.negative_ion_order = [
            '[M-H]-1',
            '[M-2H]-2',
            '[M-H-H2O]-1',
            '[M-H-MeOH]-1',
            '[M+FA-H]-1'
        ]
        
        # Lipid adduct order definitions
        self.lipid_positive_adducts = [
            'M+H',
            'M+2H', 
            'M+NH4',
            'M+H-H2O',
            'M+Na',
            'M+K',
            'M+Li'
        ]
        self.lipid_negative_adducts = [
            'M-H',
            'M-2H',
            'M-CH3',
            'M+HCOO',
            'M+CH3COO'
        ]

    def get_selected_ion_order(self, selected_ions: List[str], polarity: str = 'positive') -> List[str]:
        """
        Build ion order list based on user selection while preserving priority order
        
        Args:
            selected_ions: List of ions selected by user
            polarity: 'positive' or 'negative'
        
        Returns:
            List of ions in priority order (filtered to selected ions only)
        """
        # Reference order for each polarity
        pos_order = ['[M+H]+1', '[M+2H]+2', '[M+H-H2O]+1', '[M+H+MeOH]+1', '[M+FA+H]+1']
        neg_order = ['[M-H]-1', '[M-2H]-2', '[M-H-H2O]-1', '[M-H-MeOH]-1', '[M+FA-H]-1']
        
        reference = pos_order if polarity == 'positive' else neg_order
        
        # Filter to only selected ions, maintaining order
        return [ion for ion in reference if ion in selected_ions]
    
    def clean_column_names(self, df: pd.DataFrame, patterns_to_remove: List[str]) -> tuple:
        """
        Clean numeric column names by removing specified patterns
        
        Args:
            df: DataFrame with columns to clean
            patterns_to_remove: List of text patterns to remove
        
        Returns:
            Tuple of (cleaned DataFrame, dict mapping old->new names)
        """
        rename_map = {}
        
        for col in df.columns:
            clean_name = col
            for pattern in patterns_to_remove:
                if pattern:  # Skip empty patterns
                    clean_name = re.sub(re.escape(pattern), '', clean_name, flags=re.IGNORECASE)
            
            # Clean up extra spaces and special characters
            clean_name = ' '.join(clean_name.split())  # Remove multiple spaces
            clean_name = clean_name.strip()
            
            if clean_name != col and clean_name:  # Only rename if name changed and not empty
                rename_map[col] = clean_name
        
        # Apply renaming
        if rename_map:
            df_clean = df.rename(columns=rename_map)
            return df_clean, rename_map
        
        return df, rename_map

    def _normalize_kegg_id(self, value: str) -> str:
        """Normalize KEGG ID to allowed formats C##### or D##### (uppercase).
        Removes prefixes like 'cpd:', 'cpd+', 'dr:', 'dr+'. Returns '' if not matching.
        Also handles semicolon-separated lists like 'C00022; C00023'
        """
        try:
            if value is None:
                return ''
            s = str(value).strip()
            if not s:
                return ''
            
            # Handle semicolon-separated lists
            if ';' in s:
                parts = [part.strip() for part in s.split(';')]
                normalized_parts = []
                for part in parts:
                    normalized = self._normalize_single_kegg_id(part)
                    if normalized:
                        normalized_parts.append(normalized)
                return '; '.join(normalized_parts) if normalized_parts else ''
            
            return self._normalize_single_kegg_id(s)
        except Exception:
            return ''

    def _normalize_single_kegg_id(self, s: str) -> str:
        """Normalize a single KEGG ID"""
        s_lower = s.lower()
        for prefix in ('cpd:', 'cpd+', 'dr:', 'dr+'):
            if s_lower.startswith(prefix):
                s = s[len(prefix):]
                break
        s = s.strip().upper()
        import re
        return s if re.fullmatch(r'[CD]\d{5}', s) else s  # Return original if not matching

    def _normalize_ion(self, ion: str) -> str:
        """Normalize a Reference Ion string to explicit charge-count form.

        Rules:
        - Trim spaces
        - Ensure trailing charge includes +N/-N; if only + or -, assume 1
        - Leave bracketed content as-is
        Examples:
            '[M+H]+'   -> '[M+H]+1'
            '[M-H]-'   -> '[M-H]-1'
            '[M+2H]+2' -> '[M+2H]+2' (unchanged)
            'M+H'     -> '[M+H]+1'
        """
        try:
            if ion is None:
                return ''
            s = str(ion).strip()
            if not s:
                return s
            
            # Handle bracketed ions
            if s.startswith('[') and ']' in s:
                bracket_end = s.find(']') + 1
                base_ion = s[:bracket_end]
                charge_part = s[bracket_end:]
                
                # If no charge part, add +1
                if not charge_part:
                    return base_ion + '+1'
                # If ends with + or - followed by optional digits
                if charge_part.endswith('+'):
                    return base_ion + charge_part + '1'
                if charge_part.endswith('-'):
                    return base_ion + charge_part + '1'
                # Already has charge number
                return s
            else:
                # Handle non-bracketed ions like 'M+H'
                if '+' in s:
                    return f'[{s}]+1'
                elif '-' in s:
                    return f'[{s}]-1'
                # If it has a charge number already
                return s
        except Exception:
            return str(ion)

    def _default_progress_callback(self, message: str):
        """Default progress callback that prints to console"""
        try:
            print(message.strip())
        except UnicodeEncodeError:
            # Fall back to ASCII-safe output when the console encoding
            # cannot represent some characters (e.g., emoji on Windows cp1252)
            safe = message.strip().encode('ascii', errors='replace').decode('ascii')
            print(safe)

    def _log(self, message: str):
        """Log a message using the progress callback"""
        self.progress_callback(message)
    
    def _update_progress(self, current: int, total: int, message: str):
        """Update progress with current step and message"""
        progress_msg = f"[{current}/{total}] {message}"
        self._log(progress_msg)
        logging.info(progress_msg)

    def clean_data(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Main entry point for data cleaning - compatible with GUI interface

        Args:
            config: Configuration dictionary containing:
                - cleaning_mode: "combined" or "separate"
                - files: Dictionary with file paths
                - output: Dictionary with output settings
                - use_negative_mode: Boolean for negative mode processing
                - stage1_options: Optional dictionary of Stage 1 cleanup toggles
                - verified_dataframes: Optional dict with pre-loaded DataFrames {key: dataframe}
                - verify_columns_used: Boolean indicating if verify columns was used

        Returns:
            Dictionary with results and output file paths
        """
        try:
            self._log("🚀 Starting data cleaning process...\n")

            results = {
                'success': False,
                'files_created': [],
                'pos_result': None,
                'neg_result': None,
                'combined_result': None,
                'metabolite_ids_df': None,           # Combined metabolite IDs
                'pos_metabolite_ids_df': None,       # Positive metabolite IDs
                'neg_metabolite_ids_df': None,       # Negative metabolite IDs
                'pos_enhanced_df': None,             # Enhanced positive DataFrame
                'neg_enhanced_df': None,             # Enhanced negative DataFrame
                'combined_enhanced_df': None,        # Enhanced combined DataFrame
                'message': '',
                'error': None
            }

            stage1_options = config.get('stage1_options')
            column_assignments = config.get('column_assignments', {})  # Get verified column assignments
            pos_assignments = config.get('pos_column_assignments', column_assignments)  # Specific for positive
            neg_assignments = config.get('neg_column_assignments', column_assignments)  # Specific for negative
            verified_dataframes = config.get('verified_dataframes', {})  # Pre-loaded DataFrames
            verify_columns_used = config.get('verify_columns_used', False)

            # Handle different cleaning modes
            if config['cleaning_mode'] == "combined":
                # For combined mode, process the single file
                input_file = config['files']['combined_file']
                output_filename = f"{config['output']['combined_filename']}.xlsx"
                output_file = os.path.join(config['output']['directory'], output_filename)

                # Clean the data - pass pre-loaded dataframe if available from verify columns
                combined_df = verified_dataframes.get('combined') if verify_columns_used else None
                clean_results = self._clean_single_file(
                    input_file, 
                    output_file, 
                    stage1_options=stage1_options, 
                    column_assignments=column_assignments,
                    pre_loaded_df=combined_df  # Pass verified dataframe if available
                )

                if clean_results['success']:
                    results['combined_result'] = clean_results['data']
                    # Store metabolite IDs if available
                    if clean_results.get('metabolite_ids_df') is not None:
                        results['metabolite_ids_df'] = clean_results['metabolite_ids_df']
                    results['files_created'].append(output_file)
                    results['success'] = True
                    results['message'] = "Combined file processed successfully"

            else:
                # For separate mode, process positive file with pos-specific assignments
                pos_input = config['files']['pos_raw_file']
                # Stage 1: do not write pos file; keep in-memory only
                pos_df = verified_dataframes.get('positive') if verify_columns_used else None  # Get verified positive dataframe
                pos_results = self._clean_single_file(
                    pos_input,
                    None,
                    mzcloud_file=config['files'].get('pos_metabolika_file'),
                    stage1_options=stage1_options,
                    column_assignments=pos_assignments,  # Use positive-specific assignments
                    pre_loaded_df=pos_df  # Pass verified dataframe if available
                )

                if pos_results['success']:
                    results['pos_result'] = pos_results['data']
                    # Store positive metabolite IDs if available
                    if pos_results.get('metabolite_ids_df') is not None:
                        results['pos_metabolite_ids_df'] = pos_results['metabolite_ids_df']

                # Process negative file if requested with neg-specific assignments
                if config.get('use_negative_mode', False) and 'neg_raw_file' in config['files']:
                    neg_input = config['files']['neg_raw_file']
                    # Stage 1: do not write neg file; keep in-memory only
                    neg_df = verified_dataframes.get('negative') if verify_columns_used else None  # Get verified negative dataframe
                    neg_results = self._clean_single_file(
                        neg_input,
                        None,
                        mzcloud_file=config['files'].get('neg_metabolika_file'),
                        stage1_options=stage1_options,
                        column_assignments=neg_assignments,  # Use negative-specific assignments
                        pre_loaded_df=neg_df  # Pass verified dataframe if available
                    )

                    if neg_results['success']:
                        results['neg_result'] = neg_results['data']
                        # Store negative metabolite IDs if available
                        if neg_results.get('metabolite_ids_df') is not None:
                            results['neg_metabolite_ids_df'] = neg_results['metabolite_ids_df']

                # Create combined result (simplified) when processing separate mode
                if results['pos_result'] is not None or results['neg_result'] is not None:
                    # Build pos/neg subsets (Name + Formula) when available
                    subsets = []
                    if results['pos_result'] is not None:
                        pos_subset = results['pos_result'][['Name', 'Formula']].copy() if 'Formula' in results['pos_result'].columns else results['pos_result'][['Name']].copy()
                        subsets.append(pos_subset)
                    if results['neg_result'] is not None:
                        neg_subset = results['neg_result'][['Name', 'Formula']].copy() if 'Formula' in results['neg_result'].columns else results['neg_result'][['Name']].copy()
                        subsets.append(neg_subset)

                    # Concatenate available subsets and deduplicate on Name
                    if subsets:
                        combined_subset = pd.concat(subsets, ignore_index=True)
                        combined_subset = combined_subset.drop_duplicates(subset=['Name'], keep='first').reset_index(drop=True)
                    else:
                        combined_subset = pd.DataFrame(columns=['Name', 'Formula'])

                    # Expose the simplified two-column Combined subset in results for downstream consumers
                    results['combined_subset'] = combined_subset

                    # Also create a full combined_data (full rows) for backward compatibility
                    full_parts = [df for df in (results.get('pos_result'), results.get('neg_result')) if df is not None]
                    if full_parts:
                        combined_data = pd.concat(full_parts, ignore_index=True)
                        # Deduplicate combined data by Name to avoid redundant ID annotation searches
                        original_count = len(combined_data)
                        if 'Name' in combined_data.columns:
                            combined_data = combined_data.drop_duplicates(subset=['Name'], keep='first').reset_index(drop=True)
                            dedup_count = original_count - len(combined_data)
                            if dedup_count > 0:
                                self._log(f"   🔍 Deduplicated combined data: removed {dedup_count} duplicate metabolite(s)\n")
                                self._log(f"   Unique metabolites for ID annotation: {len(combined_data)}\n")
                        results['combined_result'] = combined_data
                    else:
                        results['combined_result'] = combined_subset

                    # Save combined file with the simplified Combined sheet as FIRST sheet
                    combined_filename = f"{config['output']['combined_filename']}.xlsx"
                    combined_output = os.path.join(config['output']['directory'], combined_filename)

                    with pd.ExcelWriter(combined_output, engine='openpyxl') as writer:
                        # Write simplified Combined sheet first
                        combined_subset.to_excel(writer, sheet_name='Combined', index=False)

                        # Then write full Positive/Negative sheets if present
                        if results['pos_result'] is not None:
                            results['pos_result'].to_excel(writer, sheet_name='Positive', index=False)
                        if results['neg_result'] is not None:
                            results['neg_result'].to_excel(writer, sheet_name='Negative', index=False)

                        # Optionally include combined Metabolites_IDS aggregated from provided mzcloud files
                        # If mzcloud files were provided in config, try to include an aggregated Metabolites_IDS sheet
                        mzcloud_parts = []
                        if config['files'].get('pos_metabolika_file'):
                            try:
                                mpos = pd.read_excel(config['files']['pos_metabolika_file'])
                                mzcloud_parts.append(mpos)
                            except Exception:
                                pass
                        if config['files'].get('neg_metabolika_file'):
                            try:
                                mneg = pd.read_excel(config['files']['neg_metabolika_file'])
                                mzcloud_parts.append(mneg)
                            except Exception:
                                pass

                        if mzcloud_parts:
                            try:
                                mz_all = pd.concat(mzcloud_parts, ignore_index=True)
                                mz_clean = self._clean_metabolika_ids_df(mz_all)
                                mz_clean.to_excel(writer, sheet_name='Metabolites_IDS', index=False)
                                # Store combined metabolite IDs in memory
                                results['metabolite_ids_df'] = mz_clean
                            except Exception:
                                pass

                    results['files_created'].append(combined_output)

                if results['pos_result'] is not None or results['neg_result'] is not None:
                    results['success'] = True
                    results['message'] = "Processing completed successfully"

            if not results['success']:
                results['message'] = "No data was processed successfully"
                self._log("❌ No data was processed successfully\n")

            return results

        except Exception as e:
            self._log(f"❌ Error: {str(e)}\n")
            return {
                'success': False,
                'files_created': [],
                'message': f"Data cleaning failed: {str(e)}",
                'error': str(e)
            }

    def _clean_single_file(self, input_file: str, output_file: Optional[str], mzcloud_file: Optional[str] = None, stage1_options: Optional[Dict[str, Any]] = None, column_assignments: Optional[Dict[str, str]] = None, pre_loaded_df: Optional[pd.DataFrame] = None) -> Dict[str, Any]:
        """Clean a single file and return results

        Args:
            input_file: path to raw input excel
            output_file: path to save cleaned excel
            mzcloud_file: optional path to mzCloud/metabolika excel containing ID columns to clean and write to Metabolites_IDS
            column_assignments: optional dictionary mapping expected column names to actual columns in the file
            pre_loaded_df: optional pre-loaded DataFrame (from verify columns step) to use instead of reading from file
        """
        try:
            # Ensure column_assignments is a dict
            if column_assignments is not None and not isinstance(column_assignments, dict):
                self._log(f"⚠️ Warning: column_assignments is {type(column_assignments)}, converting to dict\n")
                column_assignments = {}
            
            # Use pre-loaded dataframe if available, otherwise read from file
            if pre_loaded_df is not None and isinstance(pre_loaded_df, pd.DataFrame):
                self._log(f"📂 Using pre-loaded verified dataframe\n")
                df = pre_loaded_df.copy()  # Use copy to avoid modifying original
            else:
                # Validate input file
                if not input_file or input_file is None:
                    raise ValueError("input_file cannot be None or empty")
                
                # Read the input file
                self._log(f"📂 Loading input file: {input_file}\n")
                df = pd.read_excel(input_file, sheet_name=0)

            self._log(f"📊 Input data: {len(df)} rows, {len(df.columns)} columns\n")
            
            # Apply column assignments (rename columns according to verified mapping)
            # Only rename if the actual column name differs from what backend expects
            # For data cleaning, backend expects: Name, Formula, Metabolika Pathways, BioCyc Pathways
            if column_assignments:
                self._log(f"📋 Processing column assignments...\n")
                # column_assignments format: {'expected_role': 'actual_column_in_file'}
                # Example: {'Feature ID': 'Name', 'Formula': 'Formula', 'Metabolika Pathways': '# Metabolika Pathways'}
                # Map expected roles to backend column names
                role_to_backend = {
                    'Feature ID': 'Name',
                    'Metabolika Pathways': 'Metabolika Pathways',
                    'BioCyc Pathways': 'BioCyc Pathways',
                    'Formula': 'Formula',
                    'Class': 'Class',
                    'Super Class': 'Super_Class'
                }
                
                rename_map = {}
                for role, actual_col in column_assignments.items():
                    if actual_col is not None and actual_col:  # Ensure actual_col is not empty
                        backend_col = role_to_backend.get(role, role)
                        # Only rename if actual column is different from what backend expects
                        if actual_col != backend_col:
                            # Verify the column exists in the dataframe before adding to rename map
                            if actual_col in df.columns:
                                rename_map[actual_col] = backend_col
                            else:
                                self._log(f"   ⚠️ Warning: Column '{actual_col}' not found in file\n")
                
                if rename_map:
                    self._log(f"   Renaming columns: {rename_map}\n")
                    df = df.rename(columns=rename_map)

            # Clean the data
            df_clean = self._clean_dataframe(df, options=stage1_options)

            # Drop rows where all numeric values are zero before any downstream use
            group_area_missing = False
            if hasattr(self, '_remove_zero_metabolites'):
                # For metabolites with Group Area columns
                df_clean_result = self._remove_zero_metabolites(df_clean)
                # Check if it was a tuple (df, flag) or just df
                if isinstance(df_clean_result, tuple):
                    df_clean, group_area_missing = df_clean_result
                else:
                    df_clean = df_clean_result
                    # Check if Group Area columns exist
                    group_area_missing = not any('Group Area:' in str(col) for col in df_clean.columns)
            elif hasattr(self, '_drop_all_zero_rows'):
                df_clean = self._drop_all_zero_rows(df_clean)

            # Prepare optional Metabolites_IDS dataframe if mzcloud/metabolika file provided
            metabolite_ids_df = None
            if mzcloud_file and mzcloud_file.strip() and mzcloud_file != 'None':
                try:
                    self._log(f"  📊 Loading mzCloud/metabolika IDs file: {mzcloud_file}\n")
                    if os.path.exists(mzcloud_file):
                        df_ids = pd.read_excel(mzcloud_file)
                        metabolite_ids_df = self._clean_metabolika_ids_df(df_ids)
                        self._log(f"  📊 Metabolites_IDS entries: {len(metabolite_ids_df)}\n")
                    else:
                        self._log(f"  ⚠️ mzCloud/metabolika file not found: {mzcloud_file}\n")
                except Exception as e:
                    self._log(f"  ⚠️ Error reading mzCloud/metabolika file: {e}\n")
                    metabolite_ids_df = None

            # Save the cleaned data (this will create Cleaned_Data, Pathway_Summary and optional Metabolites_IDS)
            if output_file:
                self._save_cleaned_data(df_clean, output_file, metabolite_ids_df=metabolite_ids_df)

            return {
                'success': True,
                'data': df_clean,
                'metabolite_ids_df': metabolite_ids_df,  # Store metabolite IDs DataFrame in memory
                'input_rows': len(df),
                'output_rows': len(df_clean),
                'group_area_missing': group_area_missing  # Flag if Group Area columns were not found
            }

        except Exception as e:
            self._log(f"❌ Error processing file {input_file}: {str(e)}\n")
            return {
                'success': False,
                'data': None,
                'error': str(e)
            }

    def _clean_dataframe(self, df: pd.DataFrame, options: Optional[Dict[str, Any]] = None) -> pd.DataFrame:
        """Clean the DataFrame by removing unneeded columns and filtering rows"""
        try:
            self._log("🧹 Cleaning DataFrame...\n")

            # Make a copy to avoid modifying original
            df_clean = df.copy()
            options = options or {}
            
            # Debug: Check for duplicate columns at start
            duplicate_cols = df_clean.columns[df_clean.columns.duplicated()].tolist()
            if duplicate_cols:
                self._log(f"   ⚠️ DEBUG: Found duplicate column names at start: {duplicate_cols}\n")
                self._log(f"   📊 DEBUG: All columns: {list(df_clean.columns)}\n")

            # Remove completely empty columns
            df_clean = df_clean.dropna(axis=1, how="all")
            self._log(f"   Removed empty columns: {len(df.columns) - len(df_clean.columns)}\n")

            # Define columns to remove (unneeded metadata)
            columns_to_remove = [
                'Checked', 'Tags', 'Comments'
            ]

            # Remove unneeded columns if they exist
            existing_remove_cols = [col for col in columns_to_remove if col in df_clean.columns]
            if existing_remove_cols:
                df_clean = df_clean.drop(columns=existing_remove_cols)
                self._log(f"   Removed unneeded columns: {existing_remove_cols}\n")

            # Filter out invalid rows
            original_rows = len(df_clean)

            # Remove rows where Name is NaN or empty
            if 'Name' in df_clean.columns:
                try:
                    df_clean = df_clean.dropna(subset=['Name'])
                    # Strip whitespace from Name column and remove empty strings
                    name_col = df_clean['Name'].astype(str)
                    if isinstance(name_col, pd.Series):
                        df_clean = df_clean[name_col.str.strip() != '']
                    self._log(f"   Removed rows with empty Name: {original_rows - len(df_clean)}\n")
                except Exception as e:
                    self._log(f"   ⚠️ Error filtering Name column: {str(e)}\n")
            else:
                self._log(f"   ⚠️ Name column not found in DataFrame. Available columns: {list(df_clean.columns)[:5]}...\n")

            # Remove rows that appear to be headers or pathway names
            if 'Name' in df_clean.columns:
                header_values = ['Pathway Name', 'Structure', 'Checked', 'Tags', 'Comments']
                rows_before_header_filter = len(df_clean)
                mask = ~df_clean['Name'].isin(header_values)
                df_clean = df_clean[mask].reset_index(drop=True)
                self._log(f"   Removed header-like rows: {rows_before_header_filter - len(df_clean)}\n")

            # Ensure essential columns are present
            essential_columns = [
                'Name', 'Formula', 'Metabolika Pathways', 'BioCyc Pathways'
            ]

            missing_essential = [col for col in essential_columns if col not in df_clean.columns]
            if missing_essential:
                self._log(f"⚠️ Warning: Missing essential columns: {missing_essential}\n")

            # Debug: Check for duplicate columns before pathway processing
            duplicate_cols = df_clean.columns[df_clean.columns.duplicated()].tolist()
            if duplicate_cols:
                self._log(f"   ⚠️ DEBUG: Found duplicate columns before pathway processing: {duplicate_cols}\n")
                self._log(f"   📊 DEBUG: Column names and counts:\n")
                for col in set(duplicate_cols):
                    count = (df_clean.columns == col).sum()
                    self._log(f"      - '{col}': {count} occurrences\n")
                # Remove all but first occurrence of duplicates
                self._log(f"   🔧 DEBUG: Removing duplicate columns...\n")
                df_clean = df_clean.loc[:, ~df_clean.columns.duplicated(keep='first')]
                self._log(f"   ✅ DEBUG: Duplicate columns removed. Remaining columns: {len(df_clean.columns)}\n")

            # Clean up pathway columns - replace NaN with empty strings
            pathway_columns = ['Metabolika Pathways', 'BioCyc Pathways']
            for col in pathway_columns:
                if col in df_clean.columns:
                    try:
                        col_data = df_clean[col]
                        # Ensure it's a Series before calling string methods
                        if isinstance(col_data, pd.Series):
                            # Safely convert to string, handling NaN and other types
                            # Use .loc to avoid reindexing issues with non-unique indices
                            df_clean.loc[:, col] = col_data.fillna('').astype(str).str.strip()
                        else:
                            # If not a Series, convert to Series first with integer index
                            temp_series = pd.Series(col_data, index=range(len(df_clean))).fillna('').astype(str).str.strip()
                            df_clean.loc[:, col] = temp_series.values
                    except Exception as e:
                        self._log(f"   ⚠️ Error processing '{col}': {str(e)}\n")
                        # As fallback, try to convert column to string without strip
                        try:
                            df_clean.loc[:, col] = df_clean[col].astype(str).fillna('')
                        except Exception as e2:
                            self._log(f"   ⚠️ Could not convert '{col}' to string: {str(e2)}\n")

            # Remove metabolites with all-zero Group Area values
            result = self._remove_zero_metabolites(df_clean)
            # Unpack tuple: (dataframe, group_area_missing_flag)
            if isinstance(result, tuple):
                df_clean, group_area_missing_internal = result
            else:
                df_clean = result
                group_area_missing_internal = False

            # Sort by Name for consistency
            # if 'Name' in df_clean.columns:
            #     df_clean = df_clean.sort_values('Name').reset_index(drop=True)

            self._log(f"   Final cleaned data: {len(df_clean)} rows, {len(df_clean.columns)} columns\n")

            return df_clean

        except Exception as e:
            self._log(f"❌ Error cleaning DataFrame: {str(e)}\n")
            raise

    def _remove_zero_metabolites(self, df: pd.DataFrame) -> tuple:
        """Remove metabolites that have zeros in ALL 'Group Area:' columns.
        
        This function identifies columns containing 'Group Area:' in their name
        and removes rows where ALL such columns contain zero, NaN, or empty values.
        
        Args:
            df: DataFrame to filter
            
        Returns:
            Tuple of (Filtered DataFrame with zero-only metabolites removed, group_area_missing_flag)
        """
        try:
            self._log("🔍 Removing metabolites with all-zero Group Area values...\n")
            
            # Identify columns containing 'Group Area:' in their name
            group_area_cols = [col for col in df.columns if isinstance(col, str) and 'Group Area:' in col]
            
            if not group_area_cols:
                self._log("\n   ⚠️⚠️⚠️ CRITICAL ERROR ⚠️⚠️⚠️\n")
                self._log("   No 'Group Area:' columns found!\n")
                self._log("   This indicates the 'Group Area' column was not properly assigned.\n")
                self._log("   \n")
                self._log("   ACTION REQUIRED:\n")
                self._log("   Please use the 'Verify Columns' button to explicitly assign\n")
                self._log("   the Group Area column before running data cleaning again.\n")
                self._log("   \n")
                # Return original df and True flag indicating missing Group Area
                return df, True
            
            self._log(f"   Found {len(group_area_cols)} 'Group Area:' columns\n")
            
            # Store original row count
            original_rows = len(df)
            
            # Create a copy to avoid modifying original
            df_filtered = df.copy()
            
            # Reset index to ensure unique index values (prevents reindexing errors)
            df_filtered = df_filtered.reset_index(drop=True)
            
            # Convert Group Area columns to numeric, coercing errors to NaN
            for col in group_area_cols:
                df_filtered[col] = pd.to_numeric(df_filtered[col], errors='coerce')
            
            # Create mask: keep rows where at least ONE Group Area column has a non-zero value
            # A row should be kept if any Group Area column is not (zero or NaN)
            mask = pd.Series([False] * len(df_filtered), index=df_filtered.index)
            
            for col in group_area_cols:
                # Mark rows that have non-zero, non-NaN values in this column
                mask |= (df_filtered[col].notna()) & (df_filtered[col] != 0)
            
            # Apply the mask to keep only rows with at least one non-zero value
            df_filtered = df_filtered[mask].reset_index(drop=True)
            
            removed_count = original_rows - len(df_filtered)
            
            if removed_count > 0:
                self._log(f"   ✅ Removed {removed_count} metabolite(s) with all-zero Group Area values\n")
                self._log(f"   Retained {len(df_filtered)} metabolites with valid data\n")
            else:
                self._log(f"   ✅ No metabolites with all-zero Group Area values found\n")
            
            # Return tuple: (dataframe, group_area_missing_flag)
            # group_area_missing is False because we found Group Area columns
            return df_filtered, False
            
        except Exception as e:
            self._log(f"❌ Error removing zero metabolites: {str(e)}\n")
            # Return original dataframe with True flag if error occurs
            return df, True

    def _save_cleaned_data(self, df: pd.DataFrame, output_file: str, metabolite_ids_df: Optional[pd.DataFrame] = None):
        """Save the cleaned data to Excel file"""
        try:
            # Validate output_file parameter
            if not output_file or output_file is None:
                raise ValueError("output_file cannot be None or empty")
            
            self._log(f"💾 Saving cleaned data to: {output_file}\n")

            # Create output directory if it doesn't exist
            output_dir = os.path.dirname(output_file)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)

            # Save to Excel
            # Before saving: normalize Reference Ion and coerce Area to numeric
            # Reset index FIRST to ensure unique indices for all operations
            df = df.reset_index(drop=True)
            
            # Remove any duplicate columns that may still exist
            duplicate_cols = df.columns[df.columns.duplicated()].tolist()
            if duplicate_cols:
                self._log(f"   ⚠️ DEBUG: Found duplicate columns in save method: {duplicate_cols}\n")
                df = df.loc[:, ~df.columns.duplicated(keep='first')]
                self._log(f"   ✅ DEBUG: Removed duplicates. Columns now: {len(df.columns)}\n")
            
            try:
                if 'Reference Ion' in df.columns:
                    df.loc[:, 'Reference Ion'] = df['Reference Ion'].astype(str).str.strip().apply(self._normalize_ion)
            except Exception as e:
                self._log(f"   ⚠️ Warning: Could not normalize Reference Ion: {e}\n")
                pass
            try:
                if 'Area (Max.)' in df.columns:
                    df.loc[:, 'Area (Max.)'] = pd.to_numeric(df['Area (Max.)'], errors='coerce')
            except Exception as e:
                self._log(f"   ⚠️ Warning: Could not coerce Area to numeric: {e}\n")
                pass
            
            # Clean up pathway columns one more time before saving
            for col in ['Metabolika Pathways', 'BioCyc Pathways']:
                if col in df.columns:
                    try:
                        df.loc[:, col] = df[col].astype(str).fillna('').str.strip()
                    except Exception:
                        pass

            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                df.to_excel(writer, sheet_name='Cleaned_Data', index=False)

                # Create a summary sheet with pathway information
                if 'Metabolika Pathways' in df.columns or 'BioCyc Pathways' in df.columns:
                    pathway_summary = self._create_pathway_summary(df)
                    if not pathway_summary.empty:
                        pathway_summary.to_excel(writer, sheet_name='Pathway_Summary', index=False)

                # If caller provides a metabolite IDs dataframe, write it to Metabolites_IDS sheet
                # (columns already cleaned and deduplicated by _clean_metabolika_ids_df)
                if metabolite_ids_df is not None and not metabolite_ids_df.empty:
                    metabolite_ids_df.to_excel(writer, sheet_name='Metabolites_IDS', index=False)

            self._log(f"   ✅ Saved successfully\n")

            self._log(f"   ✅ Saved successfully\n")

        except Exception as e:
            self._log(f"❌ Error saving data: {str(e)}\n")
            raise

    def _create_pathway_summary(self, df: pd.DataFrame) -> pd.DataFrame:
        """Create a summary of pathway information"""
        try:
            summary_data = []
            
            # Ensure we have a clean copy with proper index
            df = df.reset_index(drop=True)

            # Count metabolites with pathways
            if 'Metabolika Pathways' in df.columns:
                try:
                    pathways_col = df['Metabolika Pathways']
                    # Handle different data types safely
                    if isinstance(pathways_col, pd.Series):
                        metabolika_count = pathways_col.astype(str).str.strip().ne('').sum()
                    else:
                        # If not a Series, try to convert first
                        metabolika_count = pd.Series(pathways_col).astype(str).str.strip().ne('').sum()
                except Exception as e:
                    self._log(f"   ⚠️ Error counting Metabolika Pathways: {e}\n")
                    metabolika_count = 0
                
                summary_data.append({
                    'Pathway_Type': 'Metabolika Pathways',
                    'Metabolites_With_Pathways': metabolika_count,
                    'Total_Metabolites': len(df),
                    'Percentage': f"{(metabolika_count / len(df) * 100):.1f}%"
                })

            if 'BioCyc Pathways' in df.columns:
                try:
                    pathways_col = df['BioCyc Pathways']
                    # Handle different data types safely
                    if isinstance(pathways_col, pd.Series):
                        biocyc_count = pathways_col.astype(str).str.strip().ne('').sum()
                    else:
                        # If not a Series, try to convert first
                        biocyc_count = pd.Series(pathways_col).astype(str).str.strip().ne('').sum()
                except Exception as e:
                    self._log(f"   ⚠️ Error counting BioCyc Pathways: {e}\n")
                    biocyc_count = 0
                
                summary_data.append({
                    'Pathway_Type': 'BioCyc Pathways',
                    'Metabolites_With_Pathways': biocyc_count,
                    'Total_Metabolites': len(df),
                    'Percentage': f"{(biocyc_count / len(df) * 100):.1f}%"
                })

            return pd.DataFrame(summary_data)

        except Exception as e:
            self._log(f"⚠️ Error creating pathway summary: {str(e)}\n")
            return pd.DataFrame()

    def _clean_metabolika_ids_df(self, df_ids: pd.DataFrame) -> pd.DataFrame:
        """Clean metabolika/mzCloud ID dataframe for Metabolites_IDS sheet

        Steps:
        - Remove 'Tags' and 'Checked' columns if present
        - Normalize 'Name' into 'Name_key' (lowercase stripped)
        - Sort by 'HMDB ID' if present
        - Remove duplicates based on 'Name_key' keeping first
        """
        try:
            df = df_ids.copy()
            # Ensure index is reset and unique to prevent reindexing errors
            df = df.reset_index(drop=True)

            # Remove empty columns
            df = df.dropna(axis=1, how='all')

            # Remove tags and checked
            for col in ['Tags', 'Checked']:
                if col in df.columns:
                    df = df.drop(columns=[col])

            # Ensure Name column exists
            if 'Name' not in df.columns:
                return pd.DataFrame()

            # Create Name_key
            df['Name_key'] = df['Name'].astype(str).str.strip().str.lower()

            # Sort by HMDB ID if present
            if 'HMDB ID' in df.columns:
                try:
                    df = df.sort_values('HMDB ID').reset_index(drop=True)
                except Exception:
                    df = df.reset_index(drop=True)

            # Drop duplicates by Name_key - ensure index is reset after this operation
            df = df.drop_duplicates(subset=['Name_key'], keep='first')
            df = df.reset_index(drop=True)

            # Normalize column names to a standard set for merging later
            rename_map = {
                'HMDB ID': 'HMDB_ID',
                'PubChem CID': 'PubChem_CID',
                'KEGG ID': 'KEGG_ID',
                'CAS Number': 'CAS',
                'IUPAC Name': 'IUPAC_Name',
            }
            # Apply renames for any present columns
            existing_renames = {k: v for k, v in rename_map.items() if k in df.columns}
            if existing_renames:
                df = df.rename(columns=existing_renames)

            # Clean KEGG IDs
            if 'KEGG_ID' in df.columns:
                df['KEGG_ID'] = df['KEGG_ID'].apply(self._normalize_kegg_id)

            # Build final ordered subset. Preserve 'Compound Class' as-is if present
            desired_order = ['Name', 'Formula', 'IUPAC_Name', 'KEGG_ID', 'CAS', 'SMILES', 'InChI', 'InChIKey', 'HMDB_ID', 'PubChem_CID', 'Compound Class', 'Name_key']
            final_cols = [c for c in desired_order if c in df.columns]
            # Ensure Name and Name_key exist
            if 'Name' in df.columns and 'Name_key' in df.columns:
                df = df[final_cols] if final_cols else df[['Name', 'Name_key']]
            else:
                df = df

            return df

        except Exception as e:
            self._log(f"⚠️ Error cleaning metabolika ids df: {e}\n")
            return pd.DataFrame()

    def load_and_clean_data(self, 
                           file_path: str, 
                           group_patterns: List[str],
                           ppm_threshold: float = 10.0,
                           rt_threshold: float = 4.0,
                           ion_mode: str = "positive",
                           ms2_filter_mode: str = 'none',
                           selected_ions: Optional[List[str]] = None,
                           column_clean_patterns: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Load and clean metabolite data from Excel file
        
        Args:
            file_path: Path to Excel file
            ppm_threshold: Maximum allowed ppm difference
            rt_threshold: RT deviation threshold for deduplication
            ion_mode: "positive" or "negative"
            ms2_filter_mode: 'none' | 'preferred_only' | 'preferred_or_nonpreferred'
            selected_ions: Optional list of user-selected reference ions to keep
            column_clean_patterns: Optional list of text patterns to remove from column names
            
        Returns:
            Cleaned DataFrame
        """
        logger = logging.getLogger(__name__)
        logger.info(f"Loading data from {file_path}")
        self._update_progress(1, 10, f"Loading data from {os.path.basename(file_path)}")
        
        # Load data
        df = pd.read_excel(file_path)
        
        # Initial column selection
        base_columns = ["Name", "Formula", "Calc. MW", "Annot. DeltaMass [ppm]",
                       "m/z", "RT [min]", "Area (Max.)", "MS2", "Reference Ion"]
        
        # Add group columns
        available_columns = list(df.columns)
        selected_columns = base_columns.copy()
        
        for pattern in group_patterns:
            matching_cols = [col for col in available_columns if pattern in col]
            selected_columns.extend(matching_cols)
        
        # Filter columns that exist
        existing_columns = [col for col in selected_columns if col in available_columns]
        df = df[existing_columns]
        
        self._update_progress(2, 10, "Performing initial filtering")
        
        # Initial filtering
        df = df[
            (df['Name'].notna()) & 
            (df['Name'] != "") &
            (df['Formula'].notna()) & 
            (df['Formula'] != "") &
            (~df['Name'].str.contains("Similar", na=False))
        ]
        
        # PPM filtering
        df['Annot. DeltaMass [ppm]'] = df['Annot. DeltaMass [ppm]'].abs()
        df = df[df['Annot. DeltaMass [ppm]'] < ppm_threshold]

        # Optional MS2 filtering based on mode
        if 'MS2' in df.columns and ms2_filter_mode != 'none':
            before = len(df)
            df['MS2'] = df['MS2'].astype(str)
            if ms2_filter_mode == 'preferred_only':
                df = df[df['MS2'].str.contains('DDA for preferred ion', case=False, na=False)]
            elif ms2_filter_mode == 'preferred_or_nonpreferred':
                df = df[df['MS2'].str.contains('DDA for preferred ion', case=False, na=False) |
                        df['MS2'].str.contains('DDA for non-preferred ion', case=False, na=False)]
            after = len(df)
            self._update_progress(3, 10, f"MS2 filter ({ms2_filter_mode}): kept {after}/{before} rows")
        
        # Create Name_Key for deduplication
        df['Name_Key'] = df['Name'].str.lower().str.strip()
        
        self._update_progress(3, 10, "Filtering reference ions")
        
        # Ion filtering - use selected ions if provided, otherwise use default order
        if selected_ions:
            # Use user-selected ions, but maintain priority order
            ion_order = self.get_selected_ion_order(selected_ions, ion_mode)
        else:
            # Use default order
            ion_order = self.positive_ion_order if ion_mode == "positive" else self.negative_ion_order
        
        df['Reference Ion'] = df['Reference Ion'].astype(str).str.strip().apply(self._normalize_ion)
        df = df[df['Reference Ion'].isin(ion_order)]
        
        self._update_progress(4, 10, "Cleaning column names")
        
        # Clean column names
        df.columns = df.columns.str.replace(r'RT \[min\]', 'RT', regex=True)
        df.columns = df.columns.str.replace(r'Annot. DeltaMass \[ppm\]', 'ppm', regex=True)
        df.columns = df.columns.str.replace(r'Calc. MW', 'MW', regex=True)
        
        # Apply column name cleaning if patterns provided
        if column_clean_patterns:
            df, rename_map = self.clean_column_names(df, column_clean_patterns)
            if rename_map:
                self._update_progress(4, 10, f"Cleaned {len(rename_map)} column names")
        
        # Clean column names
        df.columns = df.columns.str.replace(r'RT \[min\]', 'RT', regex=True)
        df.columns = df.columns.str.replace(r'Annot. DeltaMass \[ppm\]', 'ppm', regex=True)
        df.columns = df.columns.str.replace(r'Calc. MW', 'MW', regex=True)
        
        # Ensure numeric area for sorting
        df['Area (Max.)'] = pd.to_numeric(df['Area (Max.)'], errors='coerce')

        self._update_progress(5, 10, "Ranking reference ions")

        # Ion ranking and two-step stable ordering to implement the user's
        # requested behavior: first sort by Area (desc) using a stable sort
        # (preserves relative order), then sort by ion priority using a
        # stable sort so that ion-order becomes the primary grouping while
        # Area ordering is preserved inside each ion group.
        df['Area (Max.)'] = pd.to_numeric(df['Area (Max.)'], errors='coerce')
        try:
            # Step 1: Area descending (stable)
            df = df.sort_values('Area (Max.)', ascending=False, kind='mergesort')
            # Step 2: ion priority (stable) — this makes ion order primary
            df['ion_rank'] = pd.Categorical(df['Reference Ion'], categories=ion_order, ordered=True)
            df = df.sort_values('ion_rank', kind='mergesort')
            df.drop('ion_rank', axis=1, inplace=True)
        except Exception:
            # Fallback to previous combined sort if mergesort fails
            df['ion_rank'] = pd.Categorical(df['Reference Ion'], categories=ion_order, ordered=True)
            df = df.sort_values(['Area (Max.)', 'ion_rank'], ascending=[False, True])
            df.drop('ion_rank', axis=1, inplace=True)
        # After sorting by Area and ion priority, select a representative row per
        # Name_Key (the first row after the sort), then perform the normal
        # group-and-summarize to restore summed sample columns. Finally, merge
        # the representative Reference Ion/Area back into the aggregated result
        # so the metadata reflects the chosen representative.
        self._update_progress(6, 10, "Selecting representative rows and restoring sums")

        # Ensure Name_Key exists
        if 'Name_Key' not in df.columns and 'Name' in df.columns:
            df['Name_Key'] = df['Name'].astype(str).str.lower().str.strip()

        rep_df = None
        if 'Name_Key' in df.columns:
            try:
                # df is already sorted by Area desc then ion priority, so head(1) gives the representative
                rep_df = df.groupby('Name_Key', sort=False).head(1)[['Name_Key', 'Reference Ion', 'Area (Max.)']].copy()
                # Ensure unique Name_Key values with reset index and explicit deduplication
                rep_df = rep_df.drop_duplicates(subset=['Name_Key'], keep='first').reset_index(drop=True)
                # Set Name_Key as index temporarily to ensure unique values for merge
                rep_df = rep_df.set_index('Name_Key', drop=False).reset_index(drop=True)
            except Exception as e:
                self._log(f"⚠️ Warning: Could not extract representative rows: {e}\n")
                rep_df = None

        # Now aggregate/sum numeric sample columns using the existing helper
        agg_df = self._group_and_summarize(df)

        # Merge representative metadata back in (prefer rep_df values where present)
        if rep_df is not None and 'Name_Key' in agg_df.columns and not rep_df.empty:
            try:
                # Only merge the representative Reference Ion so that the
                # aggregated Area (which is now a SUM) is preserved.
                # Ensure rep_df has unique Name_Key values before merging
                rep_merge = rep_df[['Name_Key', 'Reference Ion']].drop_duplicates(subset=['Name_Key'], keep='first').reset_index(drop=True)
                agg_df = agg_df.merge(rep_merge, on='Name_Key', how='left', suffixes=('', '_rep'))
                if 'Reference Ion_rep' in agg_df.columns:
                    agg_df['Reference Ion'] = agg_df['Reference Ion_rep'].fillna(agg_df['Reference Ion'])
                    agg_df = agg_df.drop(columns=['Reference Ion_rep'])
            except Exception as e:
                self._log(f"⚠️ Warning: Could not merge representative ions: {e}\n")

        # Do not re-sort the aggregated DataFrame here; preserve ordering
        # determined earlier (representative selection + aggregation).
        # Normalize Reference Ion and coerce Area to numeric for consistency.
        if 'Reference Ion' in agg_df.columns:
            try:
                agg_df['Reference Ion'] = agg_df['Reference Ion'].astype(str).str.strip().apply(self._normalize_ion)
            except Exception:
                pass
        if 'Area (Max.)' in agg_df.columns:
            try:
                agg_df['Area (Max.)'] = pd.to_numeric(agg_df['Area (Max.)'], errors='coerce')
            except Exception:
                pass

        df = agg_df
        self._update_progress(8, 10, "Final cleanup")
        
        logger.info(f"Cleaned data: {len(df)} metabolites")
        self._update_progress(10, 10, f"Cleaning complete: {len(df)} metabolites")
        
        return df

    def _filter_polarity_columns(self, df: pd.DataFrame, keep_polarity: str) -> pd.DataFrame:
        """
        Step 7: Filter out cross-polarity sample columns.
        Remove negative columns from positive data and vice versa.
        This function detects sample columns using 'Area[' or 'Area' prefix and
        removes columns whose names include polarity suffixes or tokens that
        belong to the opposite polarity as requested in keep_polarity.
        """
        self._log(f"🧹 Step 7: Filtering polarity-specific columns ({keep_polarity} mode)...\n")

        if df.empty:
            self._log(f"   No {keep_polarity} lipids to filter\n")
            return df

        # Find sample columns (Area[...] or Area_ style)
        sample_cols = [c for c in df.columns if isinstance(c, str) and ('Area[' in c or c.startswith('Area'))]

        if keep_polarity == 'positive':
            cols_to_drop = [col for col in sample_cols if ('_neg' in col.lower() or 'neg' in col.lower() or 'negative' in col.lower() or '-neg' in col.lower())]
            opposite = 'negative'
        else:
            cols_to_drop = [col for col in sample_cols if ('_pos' in col.lower() or '_positive' in col.lower() or 'pos' in col.lower() or 'positive' in col.lower())]
            opposite = 'positive'

        if cols_to_drop:
            df = df.drop(columns=cols_to_drop)
            self._log(f"   ✅ Removed {len(cols_to_drop)} {opposite} mode columns\n")
            retained_count = len([c for c in df.columns if isinstance(c, str) and ('Area[' in c or c.startswith('Area'))])
            self._log(f"   ✅ Retained {retained_count} {keep_polarity} mode columns\n")
        else:
            self._log(f"   ✅ No {opposite} mode columns found to remove\n")

        return df
    
    def _deduplicate_by_rt(self, df: pd.DataFrame, rt_threshold: float) -> pd.DataFrame:
        """
        Deduplicate metabolites based on RT deviation
        """
        # Calculate RT statistics by Name_Key
        rt_stats = df.groupby('Name_Key')['RT'].agg(['mean', 'count']).reset_index()
        rt_stats.columns = ['Name_Key', 'RT_mean', 'RT_count']
        
        # Merge back
        df = df.merge(rt_stats, on='Name_Key', how='left')
        
        # Calculate RT deviation
        df['RT_deviation'] = df['RT'] - df['RT_mean']
        df['RT_abs_dev'] = df['RT_deviation'].abs()
        
        # Filter by RT threshold
        df_filtered = df[df['RT_deviation'].abs() < rt_threshold].copy()
        
        # Drop temporary columns for final result
        df_filtered.drop(['RT_deviation', 'RT_abs_dev'], axis=1, inplace=True)
        
        return df_filtered

    def _group_and_summarize(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Group by Name_Key and summarize data
        """
        logger = logging.getLogger(__name__)
        
        # Metadata (do not sum)
        base_meta = ["Name_Key", "Name", "Formula", "MW", "ppm", "Reference Ion",
                     "MS2", "m/z", "RT", "RT_mean", "Area (Max.)"]
        # Optional metadata if present
        optional_meta = [c for c in ["Metabolika Pathways", "BioCyc Pathways", "Polarity", "MS2_Purity"]
                         if c in df.columns]
        # Keep only metadata columns that actually exist in df
        metadata_cols = [c for c in base_meta if c in df.columns] + optional_meta
        # Deduplicate while preserving order
        seen = set()
        metadata_cols = [c for c in metadata_cols if not (c in seen or seen.add(c))]

        # Candidate sample columns = everything else (except Name_Key)
        candidate_cols = [c for c in df.columns if c not in metadata_cols and c != "Name_Key"]

        # Keep only truly numeric sample columns (coerce and keep if any numeric exists)
        sample_cols = []
        for c in candidate_cols:
            coerced = pd.to_numeric(df[c], errors="coerce")
            if coerced.notna().any():
                df[c] = coerced.fillna(0)
                sample_cols.append(c)

        if not sample_cols:
            logger.info("No numeric sample columns found to sum.")
        else:
            logger.info(f"Summing {len(sample_cols)} sample columns.")

        # Note: representative Reference Ion selection is handled by callers
        # (load_and_clean_data or enhanced_stage2_cleaning) which may choose the
        # representative row based on the desired two-step sorting behavior.
        # Here we only perform numeric aggregation and leave any representative
        # metadata merge to the caller. This ensures Area aggregation can be
        # computed as the SUM of group members while callers may inject the
        # preferred Reference Ion separately.
        rep_ion_df = None

        # Build aggregation: metadata -> first, samples -> sum
        agg_meta = {c: "first" for c in metadata_cols if c != "Name_Key" and c in df.columns}
        # If desired, you can change these:
        if "Area (Max.)" in df.columns:
            # Aggregate the Area across grouped rows by summing so the final
            # aggregated Area reflects the total signal for that metabolite.
            agg_meta["Area (Max.)"] = "sum"
        if "RT" in df.columns:
            agg_meta["RT"] = "mean"

        agg_samples = {c: "sum" for c in sample_cols}

        # Single-pass groupby/agg
        result_df = (
            df.groupby("Name_Key", as_index=False)
              .agg({**agg_meta, **agg_samples})
        )

        # Ensure Name_Key is unique after aggregation (safety check)
        if 'Name_Key' in result_df.columns:
            result_df = result_df.drop_duplicates(subset=['Name_Key'], keep='first').reset_index(drop=True)

        # Note: do not overwrite Area or other aggregated numeric metadata here;
        # callers may merge representative Reference Ion after aggregation if
        # they want the Reference Ion to reflect a particular chosen row.

        # Final column order
        final_columns = ["Name", "Name_Key", "Formula", "MW", "ppm", "Reference Ion",
                         "MS2", "m/z", "RT", "RT_mean", "Area (Max.)"] + optional_meta + sample_cols
        existing_final_cols = [c for c in final_columns if c in result_df.columns]
        
        # Do not perform an automatic final sort here; callers control ordering.
        if 'Area (Max.)' in result_df.columns:
            try:
                result_df['Area (Max.)'] = pd.to_numeric(result_df['Area (Max.)'], errors='coerce')
            except Exception:
                pass
        
        return result_df[existing_final_cols]

    def enhanced_stage2_cleaning(self, 
                                df: pd.DataFrame,
                                ppm_threshold: float = 10.0,
                                rt_threshold: float = 4.0,
                                ion_mode: str = "positive",
                                ms2_filter_mode: str = 'none',
                                progress_callback=None,
                                file_path: Optional[str] = None,
                                selected_ions: Optional[List[str]] = None,
                                column_clean_patterns: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Enhanced Stage 2 cleaning for already-cleaned data from Stage 1
        Focuses on column filtering, deduplication, and sorting without group patterns
        
        IMPORTANT: Stage 2 should ONLY work with in-memory DataFrames from Stage 1.
        Excel reading should only happen in Stage 1 (clean_data method).
        
        Args:
            df: Already-cleaned DataFrame from Stage 1 (REQUIRED)
            ppm_threshold: Maximum allowed ppm difference
            rt_threshold: RT deviation threshold for deduplication
            ion_mode: "positive" or "negative"
            ms2_filter_mode: 'none' | 'preferred_only' | 'preferred_or_nonpreferred'
            progress_callback: Optional progress callback function
            file_path: [DEPRECATED] Path to file - only for backward compatibility
            selected_ions: Optional list of user-selected reference ions to keep
            column_clean_patterns: Optional list of text patterns to remove from column names
            
        Returns:
            Enhanced cleaned DataFrame
        """
        logger = logging.getLogger(__name__)
        
        # Validate the provided DataFrame (required for Stage 2)
        if df is None:
            if file_path:
                logger.warning(f"Stage 2 should use DataFrames, not files. Falling back to file: {file_path}")
                if not os.path.exists(file_path):
                    raise ValueError(f"File path does not exist: {file_path}")
                df = pd.read_excel(file_path)
                logger.info(f"Loaded {len(df)} rows, {len(df.columns)} columns from file")
            else:
                raise ValueError("Stage 2 requires a DataFrame from Stage 1. Neither df nor file_path provided.")
        else:
            # Validate the provided DataFrame
            if not isinstance(df, pd.DataFrame):
                raise ValueError(f"Expected DataFrame, got {type(df)}")
            if df.empty:
                logger.warning("Received empty DataFrame for enhanced stage 2 cleaning")
                return pd.DataFrame()  # Return empty DataFrame for empty input
            df = df.copy()
            logger.info(f"Stage 2: Processing in-memory DataFrame: {len(df)} rows, {len(df.columns)} columns")
        
        if progress_callback:
            progress_callback("[1/8] Processing cleaned data from Stage 1")
        original_row_count = len(df)
        last_nonempty_df = df.copy()
        
        if progress_callback:
            progress_callback(f"[2/8] Filtering essential feature columns")
        
        # Step 1: Keep only essential feature columns
        essential_columns = [
            "Name", "Formula", "Annot. DeltaMass [ppm]", "Calc. MW", "m/z", 
            "Reference Ion", "RT [min]", "Area (Max.)", "Metabolika Pathways", 
            "BioCyc Pathways", "Polarity", "MS2", "MS2 Purity [%]"
        ]
        
        # Find existing essential columns (handle missing columns gracefully)
        existing_essential = [col for col in essential_columns if col in df.columns]
        logger.info(f"Found {len(existing_essential)} essential columns out of {len(essential_columns)}")
        
        # Get all column names for numeric column detection
        all_columns = list(df.columns)
        
        if progress_callback:
            progress_callback(f"[3/8] Removing annotation and match columns")
        
        # Step 2: Remove unwanted columns containing specific strings
        unwanted_patterns = ["Peak Rating", "Annot. Source:", "#", "Match"]
        columns_to_remove = []
        
        for col in all_columns:
            for pattern in unwanted_patterns:
                if pattern in str(col):
                    columns_to_remove.append(col)
                    break
        
        # Remove unwanted columns
        remaining_columns = [col for col in all_columns if col not in columns_to_remove]
        logger.info(f"Removed {len(columns_to_remove)} annotation/match columns")
        
        if progress_callback:
            progress_callback(f"[4/8] Identifying numeric sample columns")
        
        # Step 3: Identify numeric columns (sample data columns)
        numeric_columns = []
        for col in remaining_columns:
            if col not in existing_essential:
                # Check if column contains numeric data
                try:
                    pd.to_numeric(df[col], errors='coerce')
                    numeric_columns.append(col)
                except:
                    pass
        
        # Final column selection: essential + numeric
        final_columns = existing_essential + numeric_columns
        df = df[final_columns]
        
        logger.info(f"Final dataset: {len(existing_essential)} feature columns + {len(numeric_columns)} sample columns")
        
        if progress_callback:
            progress_callback(f"[5/8] Applying filtering and cleaning")
        
        # Initialize warnings list for skipped steps
        skipped_steps = []
        
        # Step 4: Basic filtering
        # Remove rows with missing essential data
        if 'Name' in df.columns:
            before = len(df)
            df = df[df['Name'].notna() & (df['Name'].astype(str).str.strip() != "")]
            after = len(df)
            logger.info(f"Name filter: kept {after}/{before} rows ({after/max(before,1):.1%})")
            if after == 0:
                logger.warning("All rows removed by Name filter – reverting and skipping this filter.")
                df = last_nonempty_df.copy()
            else:
                last_nonempty_df = df.copy()
        if 'Formula' in df.columns:
            before = len(df)
            df = df[df['Formula'].notna() & (df['Formula'].astype(str).str.strip() != "")]
            after = len(df)
            logger.info(f"Formula filter: kept {after}/{before} rows ({after/max(before,1):.1%})")
            if after == 0:
                logger.warning("All rows removed by Formula filter – reverting and skipping this filter.")
                df = last_nonempty_df.copy()
            else:
                last_nonempty_df = df.copy()
        else:
            skipped_steps.append(("Formula check", "Formula column not found - skipped data validation"))
        
        # PPM filtering if column exists
        if 'Annot. DeltaMass [ppm]' in df.columns:
            before = len(df)
            df['Annot. DeltaMass [ppm]'] = pd.to_numeric(df['Annot. DeltaMass [ppm]'], errors='coerce')
            valid_mask = df['Annot. DeltaMass [ppm]'].notna()
            numeric_count = valid_mask.sum()
            if numeric_count == 0:
                logger.warning("PPM filter skipped: no numeric values present in 'Annot. DeltaMass [ppm]'.")
            else:
                ppm_filtered = df[valid_mask & (df['Annot. DeltaMass [ppm]'].abs() < ppm_threshold)]
                after = len(ppm_filtered)
                if after == 0:
                    # Determine distribution to help user
                    try:
                        max_ppm = df['Annot. DeltaMass [ppm]'].abs().max()
                        min_ppm = df['Annot. DeltaMass [ppm]'].abs().min()
                        logger.warning(
                            f"PPM filter would remove all {before} rows (threshold={ppm_threshold}). "
                            f"Observed abs(ppm) range: {min_ppm} – {max_ppm}. Skipping filter.")
                    except Exception:
                        logger.warning("PPM filter would remove all rows; skipping, could not compute range.")
                    # Skip - keep previous df
                else:
                    logger.info(f"PPM filter: kept {after}/{before} rows ({after/max(before,1):.1%}) with abs(ppm) < {ppm_threshold}")
                    df = ppm_filtered
                    last_nonempty_df = df.copy()
        else:
            skipped_steps.append(("PPM filter", "'Annot. DeltaMass [ppm]' column not found - skipped PPM threshold filtering"))

        # Optional MS2 filtering based on mode
        if 'MS2' in df.columns and ms2_filter_mode != 'none':
            before = len(df)
            df['MS2'] = df['MS2'].astype(str)
            if ms2_filter_mode == 'preferred_only':
                ms2_filtered = df[df['MS2'].str.contains('DDA for preferred ion', case=False, na=False)]
            elif ms2_filter_mode == 'preferred_or_nonpreferred':
                ms2_filtered = df[df['MS2'].str.contains('DDA for preferred ion', case=False, na=False) |
                                   df['MS2'].str.contains('DDA for non-preferred ion', case=False, na=False)]
            else:
                ms2_filtered = df
            after = len(ms2_filtered)
            logger.info(f"MS2 filter ({ms2_filter_mode}): kept {after}/{before} rows")
            df = ms2_filtered
            last_nonempty_df = df.copy()
        elif ms2_filter_mode != 'none':
            skipped_steps.append(("MS2 filter", "'MS2' column not found - skipped MS2 data filtering"))
        
        # Ion filtering (compulsory) if Reference Ion column exists
        if 'Reference Ion' in df.columns:
            # Use selected ions if provided, otherwise use default order
            if selected_ions:
                ion_order = self.get_selected_ion_order(selected_ions, ion_mode)
            else:
                ion_order = self.positive_ion_order if ion_mode == "positive" else self.negative_ion_order
            
            df['Reference Ion'] = df['Reference Ion'].astype(str).str.strip().apply(self._normalize_ion)
            before = len(df)
            df = df[df['Reference Ion'].isin(ion_order)]
            after = len(df)
            logger.info(f"Ion filter (compulsory): kept {after}/{before} rows; allowed ions: {ion_order}")
            last_nonempty_df = df.copy()
        else:
            skipped_steps.append(("Reference Ion filter", "'Reference Ion' column not found - skipped ion filtering"))
        
        if progress_callback:
            progress_callback(f"[6/8] Cleaning column names and sorting")
        
        # Step 5: Clean column names for simplicity
        column_renames = {
            'RT [min]': 'RT',
            'Annot. DeltaMass [ppm]': 'ppm',
            'Calc. MW': 'MW',
            'MS2 Purity [%]': 'MS2_Purity'
        }
        df = df.rename(columns=column_renames)
        
        # Apply column name cleaning if patterns provided
        if column_clean_patterns:
            df, rename_map = self.clean_column_names(df, column_clean_patterns)
            if rename_map and progress_callback:
                progress_callback(f"[6/8] Cleaned {len(rename_map)} column names")
        
        # Two-step stable sorting: first Area desc, then ion priority (stable)
        if 'Reference Ion' in df.columns and 'Area (Max.)' in df.columns:
            # Use the same ion order determined above for sorting
            if selected_ions:
                ion_order = self.get_selected_ion_order(selected_ions, ion_mode)
            else:
                ion_order = self.positive_ion_order if ion_mode == "positive" else self.negative_ion_order
            
            df['Area (Max.)'] = pd.to_numeric(df['Area (Max.)'], errors='coerce')
            try:
                df = df.sort_values('Area (Max.)', ascending=False, kind='mergesort')
                df['ion_rank'] = pd.Categorical(df['Reference Ion'], categories=ion_order, ordered=True)
                df = df.sort_values('ion_rank', kind='mergesort')
                df.drop('ion_rank', axis=1, inplace=True)
            except Exception:
                df['ion_rank'] = pd.Categorical(df['Reference Ion'], categories=ion_order, ordered=True)
                df = df.sort_values(['Area (Max.)', 'ion_rank'], ascending=[False, True])
                df.drop('ion_rank', axis=1, inplace=True)
        elif 'Area (Max.)' in df.columns:
            df['Area (Max.)'] = pd.to_numeric(df['Area (Max.)'], errors='coerce')
            df = df.sort_values('Area (Max.)', ascending=False)
        else:
            skipped_steps.append(("Area-based sorting", "'Area (Max.)' column not found - deduplication without intensity-based prioritization"))
        
        if progress_callback:
            progress_callback(f"[7/8] Performing RT-based deduplication")
        
        # Step 6: RT-based deduplication
        if 'Name' in df.columns and 'RT' in df.columns:
            before = len(df)
            dedup_df = self._enhanced_deduplicate_by_rt(df, rt_threshold)
            after = len(dedup_df)
            if after == 0:
                logger.warning(f"RT-based deduplication removed all rows (threshold={rt_threshold}); reverting.")
            else:
                logger.info(f"RT-based deduplication: kept {after}/{before} rows ({after/max(before,1):.1%})")
                df = dedup_df
                last_nonempty_df = df.copy()
        elif 'Name' in df.columns:
            skipped_steps.append(("RT-based deduplication", "'RT [min]' column not found - skipped retention time-based duplicate merging"))
        
        if progress_callback:
            progress_callback(f"[8/8] Final cleanup and saving")
        
        # Step 7: Final grouping & summarization (sum numeric sample columns)
        if 'Name' in df.columns:
            before = len(df)
            df['Name_Key'] = df['Name'].astype(str).str.lower().str.strip()
            df = self._group_and_summarize(df)
            after = len(df)
            logger.info(f"Final grouping: {after} unique names (from {before} rows)")

        # Reset index
        df = df.reset_index(drop=True)

        # No fallback to non-empty state here to respect compulsory filters
        logger.info(f"Enhanced cleaning complete: {len(df)} metabolites")

        # Return dict with data and warnings about skipped steps
        return {
            'data': df,
            'skipped_steps': skipped_steps
        }
    
    def _enhanced_deduplicate_by_rt(self, df: pd.DataFrame, rt_threshold: float) -> pd.DataFrame:
        """
        Enhanced RT-based deduplication for stage 2 cleaning
        """
        if 'Name' not in df.columns or 'RT' not in df.columns:
            return df
        
        # Create Name_Key for grouping
        df = df.copy()
        df['Name_Key'] = df['Name'].astype(str).str.lower().str.strip()

        # Calculate RT statistics by Name_Key
        rt_stats = df.groupby('Name_Key', as_index=False)['RT'].agg(RT_mean='mean', RT_count='count')

        # Remove any pre-existing RT_mean/RT_count to avoid suffixing on merge
        for tmp_col in ['RT_mean', 'RT_count']:
            if tmp_col in df.columns:
                df = df.drop(columns=[tmp_col])

        # Merge back
        df = df.merge(rt_stats, on='Name_Key', how='left')

        # Calculate RT deviation (guard against missing RT_mean)
        df['RT_deviation'] = df['RT'] - df['RT_mean']
        
        # Filter by RT threshold
        df_filtered = df[df['RT_deviation'].abs() < rt_threshold].copy()
        
        # Drop temporary columns
        columns_to_drop = ['Name_Key', 'RT_deviation', 'RT_mean', 'RT_count']
        existing_temp_cols = [col for col in columns_to_drop if col in df_filtered.columns]
        df_filtered = df_filtered.drop(columns=existing_temp_cols)
        
        return df_filtered

    def store_enhanced_dataframes(self, results: Dict[str, Any], 
                                 pos_enhanced_df: Optional[pd.DataFrame] = None,
                                 neg_enhanced_df: Optional[pd.DataFrame] = None,
                                 combined_enhanced_df: Optional[pd.DataFrame] = None) -> Dict[str, Any]:
        """
        Store enhanced DataFrames in the results dictionary for memory access
        
        Args:
            results: The results dictionary to update
            pos_enhanced_df: Enhanced positive DataFrame
            neg_enhanced_df: Enhanced negative DataFrame  
            combined_enhanced_df: Enhanced combined DataFrame
            
        Returns:
            Updated results dictionary with enhanced DataFrames stored
        """
        if pos_enhanced_df is not None:
            results['pos_enhanced_df'] = pos_enhanced_df
            self._log(f"📊 Stored positive enhanced DataFrame: {len(pos_enhanced_df)} rows\n")
        if neg_enhanced_df is not None:
            results['neg_enhanced_df'] = neg_enhanced_df
            self._log(f"📊 Stored negative enhanced DataFrame: {len(neg_enhanced_df)} rows\n")
        if combined_enhanced_df is not None:
            results['combined_enhanced_df'] = combined_enhanced_df
            self._log(f"📊 Stored combined enhanced DataFrame: {len(combined_enhanced_df)} rows\n")
        return results

    # ==================== LIPID CLEANING FUNCTIONS ====================
    
    def clean_lipid_data(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Main entry point for lipid data cleaning
        
        Args:
            config: Configuration dictionary containing:
                - input_file: Path to LipidSearch export file
                - output_file: Path for output Excel file
                - adduct_filter: List of adducts to keep (optional)
                - use_default_adducts: Boolean to use default priority adducts
                - column_assignments: Optional dictionary mapping expected column names to actual columns
                
        Returns:
            Dictionary with results and output file paths
        """
        try:
            self._log("🚀 Starting lipid data cleaning process...\n")
            
            results = {
                'success': False,
                'files_created': [],
                'pos_lipids': None,
                'neg_lipids': None,
                'combined_lipids': None,
                'pos_class_df': None,
                'neg_class_df': None,
                'message': '',
                'error': None,
                'stats': {}
            }
            
            mode = config.get('mode', 'combined')
            output_file = config.get('output_file')
            
            if not output_file:
                raise ValueError("Output file path is required")
            
            # Check if verified dataframes were provided (to avoid reloading from Excel)
            verified_dataframes = config.get('verified_dataframes', {})
            
            # Handle separate mode - process positive and negative files INDEPENDENTLY
            if mode == 'separate':
                self._log(f"\n🔄 Processing in SEPARATE mode - positive and negative handled independently\n")
                
                # Initialize as None - will only process if file is provided
                pos_cleaned = None
                neg_cleaned = None
                
                # ===== PROCESS POSITIVE FILE (if provided) =====
                pos_file = config.get('positive_file')
                if pos_file:
                    self._log(f"\n📂 ===== PROCESSING POSITIVE FILE =====\n")
                    if 'positive' in verified_dataframes and verified_dataframes['positive'] is not None:
                        self._log(f"📂 Using verified POSITIVE lipid dataframe\n")
                        pos_df = verified_dataframes['positive'].copy()
                        self._log(f"📊 Loaded {len(pos_df)} rows, {len(pos_df.columns)} columns\n")
                    else:
                        self._log(f"📂 Loading POSITIVE lipid data from: {pos_file}\n")
                        pos_df = pd.read_excel(pos_file)
                        self._log(f"📊 Loaded {len(pos_df)} rows, {len(pos_df.columns)} columns\n")
                    
                    # Get positive assignments
                    positive_assignments = config.get('positive_assignments', {})
                    self._log(f"   DEBUG: Positive assignments: {positive_assignments}\n")
                    
                    # Create config with positive-specific sample mapping
                    pos_config = config.copy()
                    pos_config['sample_mapping'] = config.get('positive_sample_mapping', None)
                    
                    # Clean positive file with its assignments
                    pos_cleaned = self._clean_single_lipid_file(
                        pos_df, 
                        positive_assignments, 
                        pos_config,
                        polarity='positive'
                    )
                else:
                    self._log(f"\n⚠️ No positive file provided - skipping positive processing\n")
                
                # ===== PROCESS NEGATIVE FILE (if provided) =====
                neg_file = config.get('negative_file')
                if neg_file:
                    self._log(f"\n📂 ===== PROCESSING NEGATIVE FILE =====\n")
                    if 'negative' in verified_dataframes and verified_dataframes['negative'] is not None:
                        self._log(f"📂 Using verified NEGATIVE lipid dataframe\n")
                        neg_df = verified_dataframes['negative'].copy()
                        self._log(f"📊 Loaded {len(neg_df)} rows, {len(neg_df.columns)} columns\n")
                    else:
                        self._log(f"📂 Loading NEGATIVE lipid data from: {neg_file}\n")
                        neg_df = pd.read_excel(neg_file)
                        self._log(f"📊 Loaded {len(neg_df)} rows, {len(neg_df.columns)} columns\n")
                    
                    # Get negative assignments
                    negative_assignments = config.get('negative_assignments', {})
                    self._log(f"   DEBUG: Negative assignments: {negative_assignments}\n")
                    
                    # Create config with negative-specific sample mapping
                    neg_config = config.copy()
                    neg_config['sample_mapping'] = config.get('negative_sample_mapping', None)
                    
                    # Clean negative file with its assignments
                    neg_cleaned = self._clean_single_lipid_file(
                        neg_df, 
                        negative_assignments, 
                        neg_config,
                        polarity='negative'
                    )
                else:
                    self._log(f"\n⚠️ No negative file provided - skipping negative processing\n")
                
                # Separate mode is complete - skip to results generation
                self._log(f"\n✅ Separate mode processing complete\n")
                if pos_cleaned is not None:
                    self._log(f"   Positive: {len(pos_cleaned)} lipids\n")
                else:
                    self._log(f"   Positive: Not processed (no file provided)\n")
                if neg_cleaned is not None:
                    self._log(f"   Negative: {len(neg_cleaned)} lipids\n")
                else:
                    self._log(f"   Negative: Not processed (no file provided)\n")
                
                # Set statistics for separate mode
                results['stats'] = {
                    'positive_after_cleanup': len(pos_cleaned) if pos_cleaned is not None else 0,
                    'negative_after_cleanup': len(neg_cleaned) if neg_cleaned is not None else 0
                }
                
                # Set column_assignments for downstream checks (separate mode already processed)
                column_assignments = {}  # Empty for separate mode since already handled
                
            else:
                # ===== COMBINED MODE =====
                self._log(f"\n🔄 Processing in COMBINED mode\n")
                
                # Check if we have a pre-loaded/verified dataframe
                if 'combined' in verified_dataframes and verified_dataframes['combined'] is not None:
                    self._log(f"📂 Using verified lipid dataframe\n")
                    df = verified_dataframes['combined'].copy()
                    self._log(f"📊 Loaded {len(df)} rows, {len(df.columns)} columns\n")
                else:
                    input_file = config.get('input_file')
                    if not input_file:
                        raise ValueError("Input file path is required")
                    self._log(f"📂 Loading lipid data from: {input_file}\n")
                    df = pd.read_excel(input_file)
                    self._log(f"📊 Loaded {len(df)} rows, {len(df.columns)} columns\n")
                
                column_assignments = config.get('column_assignments', {})
                
                # Process combined mode using verified assignments
                if column_assignments:
                    # Use the new clean_single_lipid_file method for combined mode too
                    df_processed = self._clean_single_lipid_file(
                        df,
                        column_assignments,
                        config,
                        polarity='combined'  # Will be separated later
                    )
                    
                    # In COMBINED mode, save the full processed dataframe as combined_lipids
                    # This will be saved as the primary sheet
                    combined_lipids = df_processed.copy()
                    
                    # Also separate into positive and negative based on Polarity column for optional separate sheets
                    if 'Polarity' in df_processed.columns:
                        pos_cleaned = df_processed[df_processed['Polarity'] == '+'].copy()
                        neg_cleaned = df_processed[df_processed['Polarity'] == '-'].copy()
                    else:
                        # If no polarity column, treat all as positive
                        self._log("   ⚠️ No Polarity column found, treating all as positive\n")
                        pos_cleaned = df_processed.copy()
                        neg_cleaned = pd.DataFrame()  # Empty negative
                else:
                    # No assignments - try old hardcoded pattern method (fallback)
                    self._log("   ⚠️ No verified assignments found, using fallback pattern matching\n")
            
            # OLD METHOD CODE - only runs if no assignments provided in combined mode
            # Skip this entirely for separate mode (column_assignments will be empty dict)
            if mode == 'combined' and not column_assignments:
                self._log(f"📋 Processing column assignments...\n")
                self._log(f"   DEBUG: Received assignments: {column_assignments}\n")
                # For lipids: Add prefix format to user-assigned columns when not auto-detected
                # e.g., Groupz[PC3_3D_pos_3] -> Area[PC3_3D_pos_3]
                # This allows downstream code to find columns by pattern (e.g., startswith('Area['))
                rename_map = {}
                
                for expected_col, actual_col_or_list in column_assignments.items():
                    # Handle both single values and lists of values
                    actual_cols = actual_col_or_list if isinstance(actual_col_or_list, list) else [actual_col_or_list] if actual_col_or_list else []
                    
                    for actual_col in actual_cols:
                        if actual_col is not None and actual_col != expected_col:
                            # Check if this is a pattern-based column that needs prefix
                            # Pattern columns are: Area, CalcMz, BaseRt, Delta(PPM)
                            if expected_col in ['Area', 'CalcMz', 'BaseRt', 'Delta(PPM)']:
                                # Extract the bracket portion and reconstruct with expected prefix
                                # e.g., Groupz[PC3_3D_pos_3] -> extract [PC3_3D_pos_3] -> Area[PC3_3D_pos_3]
                                if '[' in actual_col and ']' in actual_col:
                                    bracket_portion = actual_col[actual_col.index('['):]
                                    new_col_name = expected_col + bracket_portion
                                    rename_map[actual_col] = new_col_name
                            else:
                                # For non-pattern columns (LipidID, Class, etc.), rename directly
                                rename_map[actual_col] = expected_col
                
                if rename_map:
                    self._log(f"   Renaming columns: {rename_map}\n")
                    df = df.rename(columns=rename_map)
            
                original_count = len(df)
                
                # Step 1: Filter and keep only required columns
                df_clean = self._filter_lipid_columns(df)
            
                # Check if Area columns were found (critical check for lipid data)
                area_cols = [col for col in df_clean.columns if isinstance(col, str) and col.startswith('Area[') and col.endswith(']')]
                area_columns_missing = False
                if not area_cols:
                    area_columns_missing = True
                    self._log("\n   ⚠️⚠️⚠️ CRITICAL ERROR ⚠️⚠️⚠️\n")
                    self._log("   No 'Area[' columns found!\n")
                    self._log("   This indicates the Area/sample columns were not properly detected.\n")
                    self._log("   \\n")
                    self._log("   ACTION REQUIRED:\n")
                    self._log("   Please use the 'Verify Columns' button to explicitly assign\n")
                    self._log("   the Area columns before running lipid cleaning again.\n")
                    self._log("   \\n")
                
                # Step 2: Merge CalcMz, BaseRt, and Delta(PPM) columns
                df_clean = self._merge_lipid_observation_columns(df_clean)
                
                # Step 3: Filter adducts
                adduct_filter = config.get('adduct_filter', None)
                df_clean = self._filter_lipid_adducts(df_clean, adduct_filter)
                
                # Step 4: Clean LipidID and create Polarity column
                df_clean = self._clean_lipid_id_and_polarity(df_clean)
                
                # Step 5: Separate positive and negative
                pos_df, neg_df = self._separate_lipid_polarity(df_clean)
                
                # Step 6: Clean both dataframes (sum duplicates, keep first occurrence)
                pos_cleaned = self._clean_duplicate_lipids(pos_df, mode='positive')
                neg_cleaned = self._clean_duplicate_lipids(neg_df, mode='negative')
                
                # Step 7: Filter polarity-specific columns (remove cross-polarity columns)
                pos_cleaned = self._filter_polarity_columns(pos_cleaned, 'positive')
                neg_cleaned = self._filter_polarity_columns(neg_cleaned, 'negative')

                # Step 8: Keep only raw Area[...] sample columns and rename (Area[Sample] -> Sample)
                # Also apply column name cleaning if patterns provided
                column_clean_patterns = config.get('column_clean_patterns', None)
                sample_mapping = config.get('sample_mapping', None)
                pos_cleaned = self._reduce_to_area_sample_columns(pos_cleaned, mode='positive', column_clean_patterns=column_clean_patterns, sample_mapping=sample_mapping)
                neg_cleaned = self._reduce_to_area_sample_columns(neg_cleaned, mode='negative', column_clean_patterns=column_clean_patterns, sample_mapping=sample_mapping)
                
                # Collect statistics for combined mode fallback
                results['stats'] = {
                        'original_count': original_count,
                        'after_filtering': len(df_clean),
                        'positive_before_cleanup': len(pos_df),
                        'positive_after_cleanup': len(pos_cleaned),
                        'negative_before_cleanup': len(neg_df),
                        'negative_after_cleanup': len(neg_cleaned)
                    }
            
            # ===== FINAL RESULTS GENERATION (both separate and combined modes) =====
            self._log(f"\n📊 Generating final results...\n")
            
            # Log statistics after split
            self._log(f"\n📊 Splitting into positive and negative...\n")
            self._log(f"   Positive lipids: {len(pos_cleaned) if pos_cleaned is not None else 0}\n")
            self._log(f"   Negative lipids: {len(neg_cleaned) if neg_cleaned is not None else 0}\n")
            
            # Note: Grade filtering was already applied earlier (before sample mapping)
            # No need to apply it again here
            
            # Final cleanup: drop rows where ALL Area columns are zero
            self._log(f"\n🧹 Removing lipid rows with all-zero Area values...\n")
            if pos_cleaned is not None:
                before = len(pos_cleaned)
                pos_cleaned = self._drop_all_zero_rows(pos_cleaned)
                after = len(pos_cleaned)
                if before != after:
                    self._log(f"   Positive: {before} → {after} lipids (removed {before - after} all-zero rows)\n")
            if neg_cleaned is not None:
                before = len(neg_cleaned)
                neg_cleaned = self._drop_all_zero_rows(neg_cleaned)
                after = len(neg_cleaned)
                if before != after:
                    self._log(f"   Negative: {before} → {after} lipids (removed {before - after} all-zero rows)\n")
            
            # Defer class abundance table generation until after column pruning
            pos_class_df = None
            neg_class_df = None

            # Build combined (two-column) sheet for ID annotation (LipidID + Class)
            combined_lipids = None
            try:
                frames_for_combined = []
                for _df in [pos_cleaned, neg_cleaned]:
                    if _df is not None and isinstance(_df, pd.DataFrame) and not _df.empty and 'LipidID' in _df.columns:
                        subset_cols = [c for c in ['LipidID', 'Class'] if c in _df.columns]
                        frames_for_combined.append(_df[subset_cols].copy())
                if frames_for_combined:
                    combined_lipids = pd.concat(frames_for_combined, ignore_index=True)
                    # Drop duplicate LipidID keeping first
                    combined_lipids = combined_lipids.drop_duplicates(subset=['LipidID'])
                    self._log(f"   🧬 Combined lipid sheet rows: {len(combined_lipids)}\n")
            except Exception as e:
                self._log(f"   ⚠️ Could not build combined lipid sheet: {e}\n")
                combined_lipids = None

            # Step 9: Add Class_name annotation from lipid_class_annotation.txt
            self._log(f"\n🏷️  Adding lipid class name annotations...\n")
            annotation_df = self._load_lipid_class_annotation()
            
            # Apply annotation to all dataframes that have a Class column
            if pos_cleaned is not None:
                pos_cleaned = self._annotate_lipid_class_names(pos_cleaned, annotation_df)
            if neg_cleaned is not None:
                neg_cleaned = self._annotate_lipid_class_names(neg_cleaned, annotation_df)
            if combined_lipids is not None:
                combined_lipids = self._annotate_lipid_class_names(combined_lipids, annotation_df)
            if pos_class_df is not None:
                pos_class_df = self._annotate_lipid_class_names(pos_class_df, annotation_df)
            if neg_class_df is not None:
                neg_class_df = self._annotate_lipid_class_names(neg_class_df, annotation_df)
            
            # Step 10: Reorder columns to place Class_name after Class
            self._log(f"\n📋 Reordering columns...\n")
            if pos_cleaned is not None:
                pos_cleaned = self._reorder_lipid_columns(pos_cleaned)
            if neg_cleaned is not None:
                neg_cleaned = self._reorder_lipid_columns(neg_cleaned)
            if combined_lipids is not None:
                combined_lipids = self._reorder_lipid_columns(combined_lipids)
            # Class tables are already ordered correctly by _generate_lipid_class_tables

            # Step 11: Remove polarity-specific columns (0 columns) from pos and neg
            # Keep only positive samples in pos sheet and negative samples in neg sheet
            self._log(f"\n🔄 Filtering polarity-specific columns before export...\n")
            if pos_cleaned is not None:
                pos_cleaned = self._filter_polarity_columns(pos_cleaned, 'positive')
            if neg_cleaned is not None:
                neg_cleaned = self._filter_polarity_columns(neg_cleaned, 'negative')

            # Step: Remove all-zero Area columns (columns where all rows are 0)
            self._log(f"\n🧹 Removing all-zero Area columns...\n")
            if pos_cleaned is not None:
                pos_cleaned_before = len(pos_cleaned.columns)
                pos_cleaned = self._drop_all_zero_columns(pos_cleaned)
                pos_cleaned_after = len(pos_cleaned.columns)
                if pos_cleaned_before != pos_cleaned_after:
                    self._log(f"   Positive: {pos_cleaned_before} → {pos_cleaned_after} cols (removed {pos_cleaned_before - pos_cleaned_after} all-zero Area columns)\n")
            if neg_cleaned is not None:
                neg_cleaned_before = len(neg_cleaned.columns)
                neg_cleaned = self._drop_all_zero_columns(neg_cleaned)
                neg_cleaned_after = len(neg_cleaned.columns)
                if neg_cleaned_before != neg_cleaned_after:
                    self._log(f"   Negative: {neg_cleaned_before} → {neg_cleaned_after} cols (removed {neg_cleaned_before - neg_cleaned_after} all-zero Area columns)\n")

            # Now generate class abundance tables using fully cleaned pos/neg data
            if pos_cleaned is not None or neg_cleaned is not None:
                pos_class_df, neg_class_df = self._generate_lipid_class_tables(pos_cleaned, neg_cleaned)
                if pos_class_df is not None:
                    self._log(f"   🧪 Generated positive class table internally: {len(pos_class_df)} rows\n")
                else:
                    self._log("   ⚠️ No positive class table generated (missing 'Class' or sample columns)\n")
                if neg_class_df is not None:
                    self._log(f"   🧪 Generated negative class table internally: {len(neg_class_df)} rows\n")
                else:
                    self._log("   ⚠️ No negative class table generated (missing 'Class' or sample columns)\n")

            # Save results to Excel with combined sheet first (and optional class sheets)
            # Get sample_mapping based on mode:
            # - Separate mode: Files were processed independently with polarity-specific mappings already applied during cleaning
            #   However, we need to pass the mappings here for the save step since we moved mapping to save time
            # - Combined mode: Use the single sample_mapping from config
            if mode == 'separate':
                # Keep polarity mappings separate to avoid key collisions (s1-1 etc.)
                # that overwrite positive names with negative names.
                pos_mapping = config.get('positive_sample_mapping', None)
                neg_mapping = config.get('negative_sample_mapping', None)
                sample_mapping = None
                self._log(f"   🔍 DEBUG - Separate mode sample mapping:\n")
                self._log(f"      Positive mapping: {pos_mapping}\n")
                self._log(f"      Negative mapping: {neg_mapping}\n")
            else:
                # Combined mode
                sample_mapping = config.get('sample_mapping', None)
                pos_mapping = None
                neg_mapping = None
                self._log(f"   🔍 DEBUG - Combined mode sample mapping: {sample_mapping}\n")
            
            column_clean_patterns = config.get('column_clean_patterns', None)
            self._save_lipid_data(
                pos_cleaned,
                neg_cleaned,
                output_file,
                combined_df=combined_lipids,
                pos_class_df=pos_class_df,
                neg_class_df=neg_class_df,
                sample_mapping=sample_mapping,
                pos_sample_mapping=pos_mapping,
                neg_sample_mapping=neg_mapping,
                column_clean_patterns=column_clean_patterns
            )
            
            results['success'] = True
            results['pos_lipids'] = pos_cleaned
            results['neg_lipids'] = neg_cleaned
            results['combined_lipids'] = combined_lipids
            results['pos_class_df'] = pos_class_df
            results['neg_class_df'] = neg_class_df
            results['files_created'].append(output_file)
            results['message'] = "Lipid data cleaning completed successfully"
            
            self._log(f"\n✅ Lipid cleaning completed successfully!\n")
            self._log(f"💾 Output saved to: {output_file}\n")
            
            return results
            
        except Exception as e:
            self._log(f"❌ Error during lipid cleaning: {str(e)}\n")
            return {
                'success': False,
                'files_created': [],
                'message': f"Lipid cleaning failed: {str(e)}",
                'error': str(e)
            }
    
    def _load_lipid_class_annotation(self) -> Optional[pd.DataFrame]:
        """
        Load lipid class annotation database from lipid_class_annotation.txt
        
        Returns:
            DataFrame with Short_Annotation and Class columns, or None if file not found
        """
        try:
            # Try to find the annotation file - check multiple locations for both development and PyInstaller
            base_dir = os.path.dirname(os.path.abspath(__file__))
            
            possible_paths = []
            
            # Check if running as PyInstaller executable
            if getattr(sys, 'frozen', False):
                # Running as PyInstaller executable - check bundled resources first
                if hasattr(sys, '_MEIPASS'):
                    # PyInstaller bundled resources location
                    possible_paths.append(os.path.join(sys._MEIPASS, 'Databases', 'lipid_class_annotation.txt'))
                
                # Also check next to executable
                exe_dir = os.path.dirname(sys.executable)
                possible_paths.extend([
                    os.path.join(exe_dir, 'Databases', 'lipid_class_annotation.txt'),  # Databases next to exe
                    os.path.join(exe_dir, 'lipid_class_annotation.txt'),  # Same dir as exe
                ])
            
            # Standard development/runtime paths (works for both Python script and PyInstaller fallback)
            possible_paths.extend([
                os.path.join(base_dir, '..', 'Databases', 'lipid_class_annotation.txt'),  # ../Databases/ folder (MAIN LOCATION)
                os.path.join(base_dir, 'lipid_class_annotation.txt'),  # main_script/ folder (same as this file)
                os.path.join(base_dir, '..', 'lipid_class_annotation.txt'),  # Project root
                os.path.join(base_dir, 'Databases', 'lipid_class_annotation.txt'),  # main_script/Databases subfolder
            ])
            
            annotation_file = None
            for path in possible_paths:
                if os.path.exists(path):
                    annotation_file = path
                    break
            
            if not annotation_file:
                self._log(f"   ℹ️  Lipid class annotation file not found in: {possible_paths}\n")
                self._log(f"   ℹ️  Skipping class name annotation (file missing)\n")
                return None
            
            # Read the tab-delimited file
            annotation_df = pd.read_csv(annotation_file, sep='\t')
            
            # Check if required columns exist
            required_cols = ['Short_Annotation', 'Class']
            if not all(col in annotation_df.columns for col in required_cols):
                self._log(f"   ⚠️  Warning: Annotation file missing required columns: {required_cols}\n")
                return None
            
            # Keep only the columns we need
            annotation_df = annotation_df[required_cols].copy()
            
            # Strip whitespace from Short_Annotation column to prevent matching issues
            annotation_df['Short_Annotation'] = annotation_df['Short_Annotation'].astype(str).str.strip()
            annotation_df['Class'] = annotation_df['Class'].astype(str).str.strip()
            
            # Remove duplicates, keeping first occurrence
            annotation_df = annotation_df.drop_duplicates(subset=['Short_Annotation'], keep='first')
            
            self._log(f"   ✅ Loaded {len(annotation_df)} lipid class annotations\n")
            
            return annotation_df
            
        except Exception as e:
            self._log(f"   ⚠️  Error loading lipid class annotation: {e}\n")
            return None
    
    def _annotate_lipid_class_names(self, df: pd.DataFrame, annotation_df: Optional[pd.DataFrame]) -> pd.DataFrame:
        """
        Add Class_name column to lipid dataframe by mapping Class to annotation database
        
        Args:
            df: Lipid DataFrame with 'Class' column
            annotation_df: Annotation DataFrame with 'Short_Annotation' and 'Class' columns
            
        Returns:
            DataFrame with added 'Class_name' column
        """
        if df is None or df.empty:
            return df
        
        if 'Class' not in df.columns:
            self._log(f"   ⚠️  Warning: 'Class' column not found, cannot add Class_name\n")
            return df
        
        # If no annotation database, add Class_name as NA
        if annotation_df is None or annotation_df.empty:
            df['Class_name'] = 'NA'
            return df
        
        # Strip whitespace from Class column in the lipid data to match annotation database
        df_class_stripped = df['Class'].astype(str).str.strip()
        
        # Create a mapping dictionary: Short_Annotation -> Class
        class_name_map = dict(zip(annotation_df['Short_Annotation'], annotation_df['Class']))
        
        # Map the stripped Class column to Class_name, filling unmatched with 'NA'
        df['Class_name'] = df_class_stripped.map(class_name_map).fillna('NA')
        
        # Count matches
        matched = (df['Class_name'] != 'NA').sum()
        total = len(df)
        
        self._log(f"   📋 Mapped {matched}/{total} lipids to class names ({matched/total*100:.1f}%)\n")
        
        # Log any unmatched classes for debugging
        unmatched_classes = df_class_stripped[df['Class_name'] == 'NA'].unique()
        if len(unmatched_classes) > 0:
            self._log(f"   ⚠️  Unmatched classes: {', '.join(unmatched_classes[:10])}")
            if len(unmatched_classes) > 10:
                self._log(f" ... and {len(unmatched_classes) - 10} more\n")
            else:
                self._log(f"\n")
        
        return df
    
    def _reorder_lipid_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Reorder columns in lipid dataframe to place Class_name immediately after Class
        
        Args:
            df: Lipid DataFrame
            
        Returns:
            DataFrame with reordered columns
        """
        if df is None or df.empty:
            return df
        
        # If Class_name doesn't exist, return as-is
        if 'Class_name' not in df.columns:
            return df
        
        # Define desired column order for feature columns
        feature_order = [
            'LipidID', 'Class', 'Class_name', 'LipidGroup', 'Charge', 'CalcMz', 'BaseRt', 
            'SubClass', 'AdductIon', 'IonFormula', 'MolStructure', 
            'PPM_diff', 'Polarity'
        ]
        
        # Get existing feature columns in the desired order
        ordered_feature_cols = [col for col in feature_order if col in df.columns]
        
        # Get all other columns (sample columns)
        other_cols = [col for col in df.columns if col not in ordered_feature_cols]
        
        # Combine: feature columns first (in order), then sample columns
        final_column_order = ordered_feature_cols + other_cols
        
        return df[final_column_order]
    
    def _clean_single_lipid_file(self, df: pd.DataFrame, assignments: dict, config: dict, polarity: str) -> pd.DataFrame:
        """
        Clean a single lipid file using verified column assignments.
        NO hardcoded patterns - uses exactly what user verified.
        
        Args:
            df: Input dataframe
            assignments: Column assignments from verification {expected_name: actual_column(s)}
            config: Full config dict for adduct_filter, column_clean_patterns, etc.
            polarity: 'positive' or 'negative'
            
        Returns:
            Cleaned dataframe with standardized columns
        """
        self._log(f"\n🧽 Processing {polarity} file with verified assignments...\n")
        
        # Extract verified columns from assignments
        # assignments structure: {'LipidID': 'LipidID', 'Area': ['Area_sample1', 'Area_sample2'], 'ObsMz': 'CalcMz', ...}
        
        # Required feature columns that must exist or be mapped
        required_feature_cols = ['LipidID', 'Class']
        
        # Build list of columns to keep
        cols_to_keep = []
        area_columns = []
        mz_column = None
        rt_column = None
        lipid_id_column = None
        class_column = None
        adduct_ion_column = None
        other_feature_columns = {}
        
        for expected_col, actual_col_or_list in assignments.items():
            # Skip if assignment is None or empty
            if actual_col_or_list is None:
                continue
                
            actual_cols = actual_col_or_list if isinstance(actual_col_or_list, list) else [actual_col_or_list] if actual_col_or_list else []
            
            for actual_col in actual_cols:
                if actual_col and actual_col in df.columns:
                    cols_to_keep.append(actual_col)
                    
                    # Track special columns
                    if expected_col == 'Area':
                        area_columns.append(actual_col)
                    elif expected_col in ['ObsMz', 'CalcMz']:
                        mz_column = actual_col
                    elif expected_col in ['ObsRt', 'BaseRt']:
                        rt_column = actual_col
                    elif expected_col == 'LipidID':
                        lipid_id_column = actual_col
                    elif expected_col == 'Class':
                        class_column = actual_col
                    elif expected_col == 'AdductIon':
                        adduct_ion_column = actual_col
                    else:
                        # Track other feature columns (LipidGroup, SubClass, etc.)
                        other_feature_columns[expected_col] = actual_col
        
        # Enhanced debug output showing ALL essential columns
        self._log(f"\n   📋 VERIFIED COLUMN ASSIGNMENTS:\n")
        self._log(f"   {'✅' if lipid_id_column else '⚠️ '} LipidID: {lipid_id_column if lipid_id_column else '(not assigned)'}\n")
        self._log(f"   {'✅' if class_column else '⚠️ '} Class: {class_column if class_column else '(not assigned)'}\n")
        self._log(f"   {'✅' if adduct_ion_column else '⚠️ '} AdductIon: {adduct_ion_column if adduct_ion_column else '(not assigned - adduct filtering will be skipped)'}\n")
        self._log(f"   {'✅' if mz_column else '⚠️ '} Mz column: {mz_column if mz_column else '(not assigned)'}\n")
        self._log(f"   {'✅' if rt_column else '⚠️ '} Rt column: {rt_column if rt_column else '(not assigned)'}\n")
        self._log(f"   ✅ Area columns: {len(area_columns)}\n")
        if other_feature_columns:
            self._log(f"   ✅ Other features: {', '.join(other_feature_columns.keys())}\n")
        self._log(f"\n")

        if not area_columns:
            raise ValueError("No Area columns found in verified assignments")
        
        # Keep only verified columns
        df_clean = df[cols_to_keep].copy()
        
        # CRITICAL: Force Area columns to numeric and clean problematic values
        self._log(f"   🔢 Converting Area columns to numeric and cleaning values...\n")
        problem_strings = ['NA', 'NaN', 'NP', 'ND', 'P', 'na', 'nan', 'np', 'nd', 'p']
        total_replacements = 0
        
        for col in area_columns:
            if col in df_clean.columns:
                # Replace problematic strings with NaN first
                if df_clean[col].dtype == 'object' or df_clean[col].dtype == 'string':
                    mask = df_clean[col].astype(str).str.strip().isin(problem_strings)
                    replacements = mask.sum()
                    if replacements > 0:
                        total_replacements += replacements
                        df_clean.loc[mask, col] = np.nan
                
                # Force to numeric
                df_clean[col] = pd.to_numeric(df_clean[col], errors='coerce')
        
        if total_replacements > 0:
            self._log(f"   ✅ Cleaned {total_replacements} problematic values (NA, ND, NP) from Area columns\n")
        self._log(f"   ✅ Forced {len(area_columns)} Area columns to numeric\n")

        # Standardize column names for downstream processing
        rename_map = {}
        for expected_col, actual_col_or_list in assignments.items():
            # Skip if assignment is None or empty
            if actual_col_or_list is None:
                continue
                
            actual_cols = actual_col_or_list if isinstance(actual_col_or_list, list) else [actual_col_or_list] if actual_col_or_list else []

            # For Grade: normalize to Grade[<token>] so grade filtering picks them up even if not already in bracket form
            if expected_col == 'Grade' and actual_cols:
                for actual_col in actual_cols:
                    if actual_col and actual_col in df_clean.columns:
                        token = re.sub(r'(?i)^grade', '', str(actual_col)).strip()
                        token = token.strip(' _[]()') or str(actual_col).strip()
                        token = token.replace(' ', '')
                        new_name = f"Grade[{token}]"
                        rename_map[actual_col] = new_name
                continue

            # For single-value columns (not Area/Grade), rename to standard name
            if expected_col not in ['Area', 'Delta(PPM)', 'Grade'] and len(actual_cols) == 1:
                actual_col = actual_cols[0]
                if actual_col and actual_col in df_clean.columns and actual_col != expected_col:
                    # Map old names to new standard names
                    if expected_col in ['ObsMz', 'CalcMz']:
                        rename_map[actual_col] = 'CalcMz'
                    elif expected_col in ['ObsRt', 'BaseRt']:
                              rename_map[actual_col] = 'BaseRt'                
                    else:
                        rename_map[actual_col] = expected_col
        if rename_map:
            self._log(f"\ud83d\udcdd Standardizing column names: {rename_map}\\n")
        df_clean = df_clean.rename(columns=rename_map)

        # Filter adducts if specified
        adduct_filter = config.get('adduct_filter')
        if adduct_filter and 'AdductIon' in df_clean.columns:
            df_clean = self._filter_lipid_adducts(df_clean, adduct_filter)
        elif adduct_filter:
            self._log(f"   ⚠️  Adduct filtering requested but AdductIon column not available - skipping adduct filtering\n")

        # Clean LipidID and add Polarity column (skip if polarity is already known from separate mode)
        if polarity in ['positive', 'negative']:
            # Polarity is already known - clean LipidID by removing ions and add explicit Polarity column
            self._log(f"🧹 Step 4: Cleaning LipidID and removing ions (polarity already known: {polarity})...\n")
            if 'LipidID' in df_clean.columns:
                # Use the unified ion removal function (same as combined mode)
                df_clean['LipidID'] = df_clean['LipidID'].apply(self._remove_lipid_ions)
                self._log(f"   📝 Sample cleaned IDs: {', '.join(df_clean['LipidID'].head().tolist())}\n")
            # Add explicit polarity marker
            df_clean['Polarity'] = 'pos' if polarity == 'positive' else 'neg'
        else:
            # Combined mode - need to extract polarity from LipidID
            df_clean = self._clean_lipid_id_and_polarity(df_clean)

        # Clean duplicate lipids (sum duplicates) - use verified Area columns, no pattern matching needed
        df_clean = self._clean_duplicate_lipids_verified(df_clean, mode=polarity)

        # EARLY: Apply grade-based filtering BEFORE sample name mapping to avoid duplicate collisions
        grade_threshold = config.get('grade_threshold', '')
        filter_rej = config.get('filter_rej', False)
        grade_filter_remove_rows = config.get('grade_filter_remove_rows', False)
        grade_filter_row_keep_pct = config.get('grade_filter_row_keep_pct', 80)
        try:
            df_clean = self._apply_grade_filtering(
                df_clean,
                grade_threshold,
                polarity,
                filter_rej,
                remove_rows_by_grade=grade_filter_remove_rows,
                min_keep_percentage=grade_filter_row_keep_pct,
            )
        except Exception as _:
            # Continue even if grade filtering fails; better to proceed than abort
            pass

        # Drop Grade columns immediately after filtering; they are no longer needed downstream
        grade_cols_early = [c for c in df_clean.columns if isinstance(c, str) and c.startswith('Grade[') and c.endswith(']')]
        if grade_cols_early:
            df_clean = df_clean.drop(columns=grade_cols_early)
            self._log(f"   🗑️ Dropped {len(grade_cols_early)} Grade columns after early filtering\n")

        # Keep Area columns with original names (no File ID mapping in combined mode)
        # This avoids duplicate column name issues after splitting into pos/neg
        # File ID mapping can be applied later if needed, or skipped entirely
        column_clean_patterns = config.get('column_clean_patterns')
        # Get sample_mapping from config (will be None if not provided)
        sample_mapping = config.get('sample_mapping', None)
        df_clean = self._reduce_to_area_sample_columns_verified(
            df_clean,
            area_columns=area_columns,
            mode=polarity,
            column_clean_patterns=column_clean_patterns,
            sample_mapping=sample_mapping  # Use sample_mapping from config (separate mode uses polarity-specific mappings)
        )

        # Drop rows where ALL numeric sample columns are zero
        df_clean = self._drop_all_zero_rows(df_clean)

        self._log(f"\u2705 {polarity.capitalize()} file cleaned: {len(df_clean)} lipids\n")
        return df_clean

    def _filter_lipid_columns(self, df: pd.DataFrame) -> pd.DataFrame:

        """
        Step 1: Filter and keep only required columns from LipidSearch export
        """
        self._log("🔍 Step 1: Filtering required columns...\n")
        
        # Define required feature columns
        feature_cols = [
            'LipidID', 'LipidGroup', 'Charge', 'CalcMz', 'BaseRt', 
            'Class', 'SubClass', 'AdductIon', 'IonFormula', 'MolStructure'
        ]
        
        # Find existing feature columns
        existing_feature_cols = [col for col in feature_cols if col in df.columns]
        missing_feature_cols = [col for col in feature_cols if col not in df.columns]
        
        if missing_feature_cols:
            self._log(f"   ⚠️  Warning: Missing feature columns: {missing_feature_cols}\n")
        
        # Find columns containing pattern indicators
        calc_mz_cols = [col for col in df.columns if 'CalcMz[' in col]
        base_rt_cols = [col for col in df.columns if 'BaseRt[' in col]
        delta_ppm_cols = [col for col in df.columns if 'Delta(PPM)[' in col]
        area_cols = [col for col in df.columns if 'Area[' in col and 'Area[s1]' not in col]
        grade_cols = [col for col in df.columns if 'Grade[' in col]
        #height_cols = [col for col in df.columns if 'Height[' in col]
        
        self._log(f"   Found {len(calc_mz_cols)} CalcMz columns\n")
        self._log(f"   Found {len(base_rt_cols)} BaseRt columns\n")
        self._log(f"   Found {len(delta_ppm_cols)} Delta(PPM) columns\n")
        self._log(f"   Found {len(area_cols)} Area columns\n")
        self._log(f"   Found {len(grade_cols)} Grade columns\n")
        #self._log(f"   Found {len(height_cols)} Height columns\n")
        
        # Check if at least one numeric column type exists
        if not area_cols:
            raise ValueError("Error: Area columns are missing")
            
        # Combine all columns to keep
        cols_to_keep = (existing_feature_cols + calc_mz_cols + base_rt_cols + 
                       delta_ppm_cols + area_cols + grade_cols)
        
        # Remove duplicates while preserving order
        cols_to_keep = list(dict.fromkeys(cols_to_keep))
        
        # Filter dataframe
        df_filtered = df[cols_to_keep].copy()
        
        self._log(f"   ✅ Kept {len(cols_to_keep)} columns out of {len(df.columns)}\n")
        
        # Replace problematic string values ONLY in Area[...] columns
        # Important: Do NOT touch Grade[...] columns here, to preserve grade values like 'A','B','C','P'
        self._log(f"   🧹 Cleaning problematic Area values (NA, NaN, NP, ND, P) -> NaN...\n")
        problem_strings_area = ['NA', 'NaN', 'NP', 'ND', 'P', 'na', 'nan', 'np', 'nd', 'p']
        
        # Identify Area columns in the filtered dataframe
        area_cols_for_clean = [col for col in df_filtered.columns if isinstance(col, str) and col.startswith('Area[') and col.endswith(']')]
        total_replacements = 0
        for col in area_cols_for_clean:
            if df_filtered[col].dtype == 'object' or df_filtered[col].dtype == 'string':
                mask = df_filtered[col].astype(str).str.strip().isin(problem_strings_area)
                replacements = mask.sum()
                if replacements > 0:
                    total_replacements += replacements
                    df_filtered.loc[mask, col] = np.nan
        
        if total_replacements > 0:
            self._log(f"   ✅ Replaced {total_replacements} problematic Area values with NaN\n")
        else:
            self._log(f"   ✅ No problematic Area string values found\n")
        
        # CRITICAL: Force all Area columns to numeric (convert scientific notation, text numbers, etc.)
        self._log(f"   🔢 Converting all Area columns to numeric...\n")
        area_cols = [col for col in df_filtered.columns if isinstance(col, str) and col.startswith('Area[') and col.endswith(']')]
        if area_cols:
            for col in area_cols:
                # Convert to numeric, coercing all non-numeric to NaN
                df_filtered[col] = pd.to_numeric(df_filtered[col], errors='coerce')
            self._log(f"   ✅ Forced {len(area_cols)} Area columns to numeric (all non-numeric values → NaN)\n")
        
        return df_filtered
    
    def _merge_lipid_observation_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Step 2: Merge all CalcMz, BaseRt, and Delta(PPM) columns into single columns
        Keep only first numeric value and drop the sample-specific columns
        """
        self._log("🔄 Step 2: Merging observation columns...\n")
        
        # Find observation columns
        calc_mz_cols = [col for col in df.columns if 'CalcMz[' in col]
        base_rt_cols = [col for col in df.columns if 'BaseRt[' in col]
        delta_ppm_cols = [col for col in df.columns if 'Delta(PPM)[' in col]
        
        # For each row, take the first non-null value across all CalcMz columns
        if calc_mz_cols:
            df['CalcMz'] = df[calc_mz_cols].apply(
                lambda row: next((val for val in row if pd.notna(val)), np.nan), axis=1
            )
            df = df.drop(columns=calc_mz_cols)
            self._log(f"   ✅ Merged {len(calc_mz_cols)} CalcMz columns into 'CalcMz'\n")
        
        # For BaseRt columns
        if base_rt_cols:
            df['BaseRt'] = df[base_rt_cols].apply(
                lambda row: next((val for val in row if pd.notna(val)), np.nan), axis=1
            )
            df = df.drop(columns=base_rt_cols)
            self._log(f"   ✅ Merged {len(base_rt_cols)} BaseRt columns into 'BaseRt'\n")
        
        # For Delta(PPM) columns
        if delta_ppm_cols:
            df['PPM_diff'] = df[delta_ppm_cols].apply(
                lambda row: next((val for val in row if pd.notna(val)), np.nan), axis=1
            )
            df = df.drop(columns=delta_ppm_cols)
            self._log(f"   ✅ Merged {len(delta_ppm_cols)} Delta(PPM) columns into 'PPM_diff'\n")
        
        return df
    
    def _filter_lipid_adducts(self, df: pd.DataFrame, custom_adducts: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Step 3: Filter adducts based on custom selection or default priority
        If custom_adducts is provided and not empty, use those.
        Otherwise, use default priority adducts.
        """
        self._log("🎯 Step 3: Filtering adducts...\n")
        
        if 'AdductIon' not in df.columns:
            self._log("   ⚠️  Warning: AdductIon column not found, skipping adduct filtering\n")
            return df
        
        original_count = len(df)
        
        if custom_adducts and len(custom_adducts) > 0:
            # Use custom adduct list
            df_filtered = df[df['AdductIon'].isin(custom_adducts)].copy()
            self._log(f"   Using custom adducts: {custom_adducts}\n")
        else:
            # Use default priority adducts
            #default_adducts = ['M+H', 'M-H', 'M+2H', 'M-2H', 'M+NH4', 'M+H-H2O']
            df_filtered = df.copy()
            self._log(f"   No adducts applied, keeping all rows\n")
        
        filtered_count = len(df_filtered)
        self._log(f"   ✅ Kept {filtered_count} out of {original_count} rows ({filtered_count/original_count*100:.1f}%)\n")
        
        return df_filtered
    
    def _remove_lipid_ions(self, lipid_id: str) -> str:
        """
        Remove all adduct ions from a LipidID string.
        Used for both separated and combined modes.
        
        Examples:
            BiotinylPE 27:1-H        -> BiotinylPE 27:1
            BiotinylPE 29:1+H        -> BiotinylPE 29:1
            BiotinylPE 42:11+2H      -> BiotinylPE 42:11
            PA(38:10)+H              -> PA(38:10)
            Ch-D7+NH4                -> Ch-D7
            Cer(t35:0)+H-H2O         -> Cer(t35:0)
            CL(O-26:3)-2H            -> CL(O-26:3)
        """
        if pd.isna(lipid_id):
            return lipid_id
        
        lipid_str = str(lipid_id).strip()
        
        # Check if there's a closing bracket ')'
        has_bracket = ')' in lipid_str
        
        if has_bracket:
            # Find position of the last closing bracket
            last_bracket_pos = lipid_str.rfind(')')
            # Return everything up to and including the last closing bracket
            return lipid_str[:last_bracket_pos + 1].strip()
        else:
            # No bracket - lipid like "Ch-D7+H", "BiotinylPE 27:1-H", "Cer(t35:0)+H-H2O", "PC(34:1)+HCOO"
            # Robust removal: strip any trailing chain of known adduct tokens
            adduct_tokens = [
                'H-H2O', '2H', 'NH4', 'HCOO', 'CH3COO', 'Na', 'K', 'Li', 'H', 'CH3'
            ]
            # Build regex matching one or more trailing adduct tokens each prefixed by + or -
            # Example matches: +H, -H, +H-H2O, +H+Na, -2H, +HCOO
            token_pattern = r'(?:[+-](?:' + '|'.join(adduct_tokens) + '))+'
            regex = re.compile(token_pattern + r'$')
            base = regex.sub('', lipid_str).strip()
            return base
    
    def _clean_lipid_id_and_polarity(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Step 4: Clean LipidID column and create Polarity column
        Remove adduct ions from LipidID and extract polarity sign
        
        Handles both bracketed formats: PA(38:10)+H  -> PA(38:10), polarity='+'
        And non-bracketed formats:    Ch-D7+H       -> Ch-D7, polarity='+'
        
                Polarity determination precedence:
                    1. Use AdductIon column if present (leading 'M+' = positive, 'M-' = negative)
                    2. Otherwise parse adduct substring after base ID and look for known tokens.
        
                Known positive adduct tokens (any present => positive):
                    +H, +2H, +NH4, +H-H2O, +Na, +K, +Li
                Known negative adduct tokens (only used if no positive tokens detected):
                    -H, -2H, -H-H2O, -CH3, +HCOO, +CH3COO
        
                Notes:
                    - Complex multi-part adducts like '+H-H2O' must remain classified as positive.
                    - We avoid using the last sign; instead we tokenize the whole adduct string.
                    - Removal of adduct from LipidID strips everything from the first '+' or '-' that begins an adduct token.
        """
        self._log("🧹 Step 4: Cleaning LipidID and extracting polarity...\n")
        
        if 'LipidID' not in df.columns:
            self._log("   ⚠️  Warning: LipidID column not found\n")
            return df
        
        # Define known tokens for tokenized adduct parsing (fallback when AdductIon absent)
        positive_tokens = {"+H", "+2H", "+NH4", "+H-H2O", "+Na", "+K", "+Li"}
        negative_tokens = {"-H", "-2H", "-H-H2O", "-CH3", "+HCOO", "+CH3COO"}

        def row_clean(row):
            lipid_id = row.get('LipidID')
            if pd.isna(lipid_id):
                return pd.Series([lipid_id, None])
            lipid_str = str(lipid_id).strip()
            polarity = None
            # 1. Use AdductIon column when present
            if 'AdductIon' in row and pd.notna(row.get('AdductIon')):
                adduct_val = str(row.get('AdductIon')).strip()
                # Normalize adduct string for matching (remove brackets and trailing charge numbers)
                adduct_core = adduct_val.replace('[', '').replace(']', '')
                # Remove trailing charge info like '+1' or '-1'
                adduct_core = re.sub(r"[+-]\d+$", "", adduct_core)
                # Stronger detection: check against explicit adduct lists
                if any(ad in adduct_core for ad in self.lipid_positive_adducts):
                    polarity = '+'
                elif any(ad in adduct_core for ad in self.lipid_negative_adducts):
                    polarity = '-'
                else:
                    # Fallback: use first sign after 'M' in string
                    if 'M+' in adduct_core:
                        polarity = '+'
                    elif 'M-' in adduct_core:
                        polarity = '-'

            # 2. Determine base LipidID and adduct remainder
            has_bracket = ')' in lipid_str
            if has_bracket:
                last_bracket_pos = lipid_str.rfind(')')
                clean_id = lipid_str[:last_bracket_pos + 1]
                remainder = lipid_str[last_bracket_pos + 1:]
            else:
                # Use same robust stripping as _remove_lipid_ions to find base ID
                adduct_tokens = [
                    'H-H2O', '2H', 'NH4', 'HCOO', 'CH3COO', 'Na', 'K', 'Li', 'H', 'CH3'
                ]
                token_pattern = r'(?:[+-](?:' + '|'.join(adduct_tokens) + '))+'
                regex = re.compile(token_pattern + r'$')
                # Identify remainder before stripping for polarity token search
                m = regex.search(lipid_str)
                if m:
                    clean_id = lipid_str[:m.start()]
                    remainder = lipid_str[m.start():]
                else:
                    clean_id = lipid_str
                    remainder = ''

            # 3. Token classification when polarity not determined yet
            if polarity is None and remainder:
                tokens = remainder.replace('+', ' +').replace('-', ' -').split()
                if any(t in positive_tokens for t in tokens):
                    polarity = '+'
                elif any(t in negative_tokens for t in tokens):
                    polarity = '-'

            return pd.Series([clean_id.strip(), polarity])

        # Apply cleaning function to each row to obtain clean LipidID and Polarity
        df[['LipidID_clean', 'Polarity']] = df.apply(row_clean, axis=1)
        # Replace original LipidID with cleaned version and drop temporary column
        df['LipidID'] = df['LipidID_clean']
        df = df.drop(columns=['LipidID_clean'])

        polarity_counts = df['Polarity'].value_counts(dropna=True)
        polarity_dict = {str(k): int(v) for k, v in polarity_counts.items()}
        self._log(f"   ✅ Extracted polarity: {polarity_dict}\n")
        sample_ids = df['LipidID'].head(10).tolist()
        self._log(f"   📝 Sample cleaned IDs: {', '.join(map(str, sample_ids[:5]))}\n")
        return df

    def _reduce_to_area_sample_columns_verified(self, df: pd.DataFrame, area_columns: list, mode: str, column_clean_patterns: Optional[List[str]] = None, sample_mapping: Optional[Dict[str, str]] = None) -> pd.DataFrame:
        """
        Verified version: Uses exact area column names from user verification.
        Also handles Grade columns if present.
        NO pattern searching - uses what was verified.
        
        Args:
            df: DataFrame to process
            area_columns: Exact list of verified Area column names
            mode: 'positive' or 'negative'
            column_clean_patterns: Patterns to remove from column names
            sample_mapping: Dict mapping DataID to File Name
        """
        if df.empty:
            return df
        
        self._log(f"🧪 Step: Reducing to verified Area sample columns ({mode})...\n")
        
        # Feature columns to keep
        feature_candidates = [
            'LipidID','Class','Class_name','LipidGroup','Charge','CalcMz','BaseRt','SubClass',
            'AdductIon','IonFormula','MolStructure','Polarity'
        ]
        id_cols = [c for c in feature_candidates if c in df.columns]
        
        # Use VERIFIED area columns (they might be Area_sample or Area[sample] or anything else)
        verified_area_cols = [c for c in area_columns if c in df.columns]
        
        # Find corresponding Grade columns (NEW)
        # Note: In the current pipeline, Grade columns may already have been dropped after early filtering.
        grade_cols = [c for c in df.columns if isinstance(c, str) and c.startswith('Grade[') and c.endswith(']')]
        
        if not verified_area_cols:
            self._log("   ⚠️ No verified Area columns found in dataframe\n")
            return df
        
        self._log(f"   ✅ Keeping {len(verified_area_cols)} verified Area columns and {len(grade_cols)} Grade columns\n")
        
        # Build reduced dataframe including Grade columns
        reduced = df[id_cols + verified_area_cols + grade_cols].copy()
        
        # NOTE: Sample mapping is NOT applied here - it's applied as the final step before saving
        # This allows all other processing (deduplication, grade filtering, etc.) to work correctly
        # with the original column names from the user's file
        
        # Apply column name cleaning if patterns provided (but exclude Grade and Area columns)
        if column_clean_patterns:
            # Column name cleaning should NOT apply to Area[] and Grade[] columns
            # as they need to keep their bracketed format for downstream processing
            area_grade_col_set = set(grade_cols) | set(verified_area_cols)
            
            # Only clean non-Area, non-Grade columns (like feature column names)
            cols_to_clean = [c for c in reduced.columns if c not in area_grade_col_set]
            
            if cols_to_clean:
                # Create a subset dataframe without Area/Grade columns for cleaning
                df_without_area_grade = reduced[cols_to_clean].copy()
                df_without_area_grade, clean_map = self.clean_column_names(df_without_area_grade, column_clean_patterns)
                
                if clean_map:
                    # Apply the cleaned names back to the main dataframe
                    reduced = reduced.rename(columns=clean_map)
                    self._log(f"   🧹 Cleaned {len(clean_map)} feature column names (preserved {len(verified_area_cols)} Area + {len(grade_cols)} Grade columns)\n")
        
        self._log(f"   ✅ Kept {len(id_cols)} feature cols, {len(verified_area_cols)} sample cols, {len(grade_cols)} grade cols\n")
        return reduced


    def _reduce_to_area_sample_columns(self, df: pd.DataFrame, mode: str, column_clean_patterns: Optional[List[str]] = None, sample_mapping: Optional[Dict[str, str]] = None) -> pd.DataFrame:
       
        """Step 8: Retain only Area[...] sample intensity columns, Grade[...] columns, and essential ID columns.
        - Drop Height[...] and any summary/stat columns (OrgMeanArea, MeanArea, etc.)
        - Apply sample DataID -> File Name mapping if provided (to BOTH Area and Grade columns)
        - Rename Area[Sample] -> Sample
        - Keep Grade[Sample] columns unchanged (for grade-based filtering)
        - Keep LipidID (and Class if present)
        - Optionally clean column names based on provided patterns
        
        Args:
            df: DataFrame to process
            mode: 'positive' or 'negative'
            column_clean_patterns: Patterns to remove from column names
            sample_mapping: Dict mapping DataID to File Name (e.g., {'s1-1': '231_1'})
        """
        try:
            if df.empty:
                return df
            self._log(f"🧪 Step 8: Reducing to Area sample columns ({mode})...\n")
            import re
            # Identify feature / metadata columns to retain (broader than previous id-only approach)
            feature_candidates = [
                'LipidID','Class','Class_name','LipidGroup','Charge','CalcMz','BaseRt','SubClass',
                'AdductIon','IonFormula','MolStructure','Polarity'
            ]
            id_cols = [c for c in feature_candidates if c in df.columns]
            
            # Area sample columns only
            area_cols = [c for c in df.columns if isinstance(c, str) and c.startswith('Area[') and c.endswith(']')]
            
            # Grade columns (NEW - keep these for grade-based filtering)
            grade_cols = [c for c in df.columns if isinstance(c, str) and c.startswith('Grade[') and c.endswith(']')]
            
            reduced = df[id_cols + area_cols + grade_cols].copy()
            
            # Rename Area columns (Area[sample] -> sample)
            pattern_area = re.compile(r'^Area\[(.+)\]$')
            rename_map = {}
            for col in area_cols:
                m = pattern_area.match(col)
                if m:
                    new_name = m.group(1)
                    # Avoid collision with id columns
                    if new_name in reduced.columns:
                        new_name = f"sample_{new_name}"
                    rename_map[col] = new_name
            if rename_map:
                reduced = reduced.rename(columns=rename_map)
            
            # Note: Grade columns stay as Grade[sample] for correlation in _apply_grade_filtering
            
            # Apply column name cleaning if patterns provided (but exclude Grade columns)
            if column_clean_patterns:
                # Only clean non-Grade columns to preserve Grade[sample] format
                grade_col_set = set(grade_cols)
                
                # Temporarily separate Grade columns
                cols_to_clean = [c for c in reduced.columns if c not in grade_col_set]
                
                # Create a subset dataframe without Grade columns for cleaning
                df_without_grades = reduced[cols_to_clean].copy()
                df_without_grades, clean_map = self.clean_column_names(df_without_grades, column_clean_patterns)
                
                if clean_map:
                    # Apply the cleaned names back to the main dataframe
                    reduced = reduced.rename(columns=clean_map)
                    self._log(f"   🧹 Cleaned {len(clean_map)} {mode} column names (preserved {len(grade_cols)} Grade columns)\n")
            
            self._log(f"   ✅ Kept {len(id_cols)} feature cols, {len(area_cols)} sample cols, {len(grade_cols)} grade cols\n")
            return reduced
        except Exception as e:
            self._log(f"   ❌ Error in Step 8 reduction: {e}\n")
            return df
    
    def _separate_lipid_polarity(self, df: pd.DataFrame) -> tuple:
        """
        Step 5: Separate dataframe into positive and negative based on Polarity
        """
        self._log("📊 Step 5: Separating by polarity...\n")
        
        if 'Polarity' not in df.columns:
            self._log("   ⚠️  Warning: Polarity column not found, cannot separate\n")
            return df.copy(), pd.DataFrame()
        
        pos_df = df[df['Polarity'] == '+'].copy()
        neg_df = df[df['Polarity'] == '-'].copy()
        
        self._log(f"   ✅ Positive lipids: {len(pos_df)}\n")
        self._log(f"   ✅ Negative lipids: {len(neg_df)}\n")
        
        return pos_df, neg_df
    
    def _clean_duplicate_lipids_verified(self, df: pd.DataFrame, mode: str = 'positive') -> pd.DataFrame:
        """
        Clean duplicates with VERIFIED Area columns (no pattern searching needed)
        Sort by adduct priority before grouping
        Used when Area columns are already verified and known
        """
        self._log(f"🧹 Step 6: Cleaning duplicate lipids ({mode} mode)...\n")
        
        if df.empty:
            self._log(f"   No {mode} lipids to clean\n")
            return df
        
        original_count = len(df)
        self._log(f"   Starting with {original_count} rows\n")
        
        # Get adduct priority order based on mode
        self._log(f"   Getting adduct priority order for {mode} mode\n")
        if mode == 'positive':
            adduct_order = self.lipid_positive_adducts
        else:
            adduct_order = self.lipid_negative_adducts
        
        # Create adduct rank for sorting (if AdductIon column exists)
        self._log(f"   Checking for AdductIon column\n")
        if 'AdductIon' in df.columns:
            self._log(f"   Creating adduct ranking\n")
            def get_adduct_rank(x):
                if pd.isna(x) or not isinstance(x, str):
                    return len(adduct_order)
                return adduct_order.index(x) if x in adduct_order else len(adduct_order)
            
            df['adduct_rank'] = df['AdductIon'].apply(get_adduct_rank)
            self._log(f"   ✅ Adduct ranking created successfully\n")
            
            # Sort by LipidID and adduct rank
            self._log(f"   Sorting by LipidID and adduct rank\n")
            df = df.sort_values(['LipidID', 'adduct_rank'])
            self._log(f"   ✅ Sorting completed successfully\n")
        
        # Identify numeric columns - already verified as numeric Area columns
        self._log(f"   Identifying verified numeric sample columns\n")
        feature_cols = {'LipidID','Class','Class_name','LipidGroup','Charge','CalcMz','BaseRt','SubClass',
                       'AdductIon','IonFormula','MolStructure','ObsMz','ObsRt','Polarity','adduct_rank'}
        numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c]) and c not in feature_cols]
        self._log(f"   ✅ Found {len(numeric_cols)} verified numeric columns\n")
        
        # Metadata columns (take first occurrence)
        metadata_cols = [col for col in df.columns if col not in numeric_cols and col != 'LipidID' and col != 'adduct_rank']
        self._log(f"   ✅ Found {len(metadata_cols)} metadata columns\n")
        
        # Group by LipidID
        self._log(f"   Starting groupby operations\n")
        if numeric_cols:
            agg_dict = {col: 'first' for col in metadata_cols}
            agg_dict.update({col: 'sum' for col in numeric_cols})
            
            self._log(f"   Summing {len(numeric_cols)} numeric columns, keeping first for {len(metadata_cols)} metadata columns\n")
            df_grouped = df.groupby('LipidID', as_index=False).agg(agg_dict)
        else:
            self._log(f"   No numeric columns found, taking first occurrence only\n")
            df_grouped = df.groupby('LipidID', as_index=False).first()
        
        self._log(f"   ✅ Groupby completed successfully\n")
        
        # Remove adduct_rank if it exists
        if 'adduct_rank' in df_grouped.columns:
            self._log(f"   Cleaning up adduct_rank column\n")
            df_grouped = df_grouped.drop(columns=['adduct_rank'])
            self._log(f"   ✅ Removed adduct_rank column\n")
        
        # Check for remaining duplicates
        duplicates_after = df_grouped['LipidID'].duplicated().sum()
        if duplicates_after > 0:
            self._log(f"   ⚠️  Warning: {duplicates_after} duplicate LipidIDs remain after grouping\n")
        
        final_count = len(df_grouped)
        self._log(f"   ✅ Cleaned: {original_count} → {final_count} lipids\n")
        
        return df_grouped
    
    def _clean_duplicate_lipids(self, df: pd.DataFrame, mode: str = 'positive') -> pd.DataFrame:
        """
        Step 6: Clean duplicates by summing numeric columns and keeping first occurrence
        Sort by adduct priority before grouping
        """
        self._log(f"🧹 Step 6: Cleaning duplicate lipids ({mode} mode)...\n")
        
        try:
            if df.empty:
                self._log(f"   No {mode} lipids to clean\n")
                return df
            
            original_count = len(df)
            self._log(f"   Starting with {original_count} rows\n")
        
            # Get adduct priority order based on mode
            self._log(f"   Getting adduct priority order for {mode} mode\n")
            if mode == 'positive':
                adduct_order = self.lipid_positive_adducts
            else:
                adduct_order = self.lipid_negative_adducts
            
            # Create adduct rank for sorting (if AdductIon column exists)
            self._log(f"   Checking for AdductIon column\n")
            if 'AdductIon' in df.columns:
                self._log(f"   Creating adduct ranking\n")
                try:
                    def get_adduct_rank(x):
                        # Handle NaN and non-string values
                        if pd.isna(x) or not isinstance(x, str):
                            return len(adduct_order)
                        return adduct_order.index(x) if x in adduct_order else len(adduct_order)
                    
                    df['adduct_rank'] = df['AdductIon'].apply(get_adduct_rank)
                    self._log(f"   ✅ Adduct ranking created successfully\n")
                except Exception as e:
                    self._log(f"   ❌ Error creating adduct ranking: {str(e)}\n")
                    import traceback
                    self._log(f"   Traceback: {traceback.format_exc()}\n")
                    raise
                
                try:
                    # Sort by LipidID and adduct rank (so preferred adducts come first)
                    self._log(f"   Sorting by LipidID and adduct rank\n")
                    df = df.sort_values(['LipidID', 'adduct_rank'])
                    self._log(f"   ✅ Sorting completed successfully\n")
                except Exception as e:
                    self._log(f"   ❌ Error sorting dataframe: {str(e)}\n")
                    import traceback
                    self._log(f"   Traceback: {traceback.format_exc()}\n")
                    raise
        
            # Identify numeric columns (Area columns)
            self._log(f"   Identifying numeric columns\n")
            try:
                # First identify columns by name pattern
                potential_numeric_cols = [col for col in df.columns if isinstance(col, str) and 'Area[' in col and col.endswith(']')]
                self._log(f"   Found {len(potential_numeric_cols)} potential numeric columns by name pattern\n")
                
                # Now filter to only those that actually contain numeric data
                numeric_cols = []
                for col in potential_numeric_cols:
                    try:
                        # Check if column can be converted to numeric (allowing for NaN)
                        if pd.api.types.is_numeric_dtype(df[col]):
                            numeric_cols.append(col)
                        else:
                            # Try converting to numeric to see if it's convertible
                            test_conversion = pd.to_numeric(df[col], errors='coerce')
                            # If most values are numeric (not NaN after conversion), include it
                            if test_conversion.notna().sum() > 0:
                                self._log(f"   ⚠️  Column '{col}' contains mixed data, converting to numeric\n")
                                df[col] = test_conversion
                                numeric_cols.append(col)
                            else:
                                self._log(f"   ⚠️  Column '{col}' appears to be non-numeric, treating as metadata\n")
                    except Exception as col_error:
                        self._log(f"   ⚠️  Error checking column '{col}': {str(col_error)}, treating as metadata\n")
                
                self._log(f"   ✅ Verified {len(numeric_cols)} truly numeric columns\n")
            except Exception as e:
                self._log(f"   ❌ Error identifying numeric columns: {str(e)}\n")
                import traceback
                self._log(f"   Traceback: {traceback.format_exc()}\n")
                raise
            
            # Identify metadata columns to keep (take first occurrence)
            self._log(f"   Identifying metadata columns\n")
            try:
                metadata_cols = [col for col in df.columns if col not in numeric_cols and col != 'LipidID']
                self._log(f"   ✅ Found {len(metadata_cols)} metadata columns\n")
            except Exception as e:
                self._log(f"   ❌ Error identifying metadata columns: {str(e)}\n")
                import traceback
                self._log(f"   Traceback: {traceback.format_exc()}\n")
                raise
        
            # Group by LipidID
            self._log(f"   Starting groupby operations\n")
            if numeric_cols:
                self._log(f"   Found {len(numeric_cols)} numeric columns to sum\n")
                try:
                    # Create aggregation dictionary
                    self._log(f"   Creating aggregation dictionary\n")
                    agg_dict = {col: 'first' for col in metadata_cols if col in df.columns and isinstance(col, str)}
                    agg_dict.update({col: 'sum' for col in numeric_cols if col in df.columns and isinstance(col, str)})
                    
                    # Remove adduct_rank from aggregation if it exists
                    if 'adduct_rank' in agg_dict:
                        del agg_dict['adduct_rank']
                    
                    self._log(f"   ✅ Aggregation dictionary created with {len(agg_dict)} columns\n")
                except Exception as e:
                    self._log(f"   ❌ Error creating aggregation dictionary: {str(e)}\n")
                    import traceback
                    self._log(f"   Traceback: {traceback.format_exc()}\n")
                    raise
                
                try:
                    # Check for any NaN values in LipidID before grouping
                    self._log(f"   Checking for NaN LipidID values\n")
                    nan_lipids = df['LipidID'].isna().sum()
                    if nan_lipids > 0:
                        self._log(f"   ⚠️  Warning: Found {nan_lipids} rows with NaN LipidID, filling with 'Unknown'\n")
                        df['LipidID'] = df['LipidID'].fillna('Unknown')
                    
                    self._log(f"   Performing groupby and aggregation\n")
                    df_clean = df.groupby('LipidID', as_index=False).agg(agg_dict)
                    self._log(f"   ✅ Groupby completed successfully\n")
                except Exception as e:
                    self._log(f"   ❌ Error during groupby/aggregation: {str(e)}\n")
                    import traceback
                    self._log(f"   Traceback: {traceback.format_exc()}\n")
                    raise
            else:
                self._log(f"   No numeric columns found, taking first occurrence only\n")
                try:
                    # No numeric columns to sum, just take first occurrence
                    df_clean = df.groupby('LipidID', as_index=False).first()
                    self._log(f"   ✅ First occurrence groupby completed\n")
                except Exception as e:
                    self._log(f"   ❌ Error during first occurrence groupby: {str(e)}\n")
                    import traceback
                    self._log(f"   Traceback: {traceback.format_exc()}\n")
                    raise
        
            # Remove adduct_rank column if it exists
            try:
                self._log(f"   Cleaning up adduct_rank column\n")
                if 'adduct_rank' in df_clean.columns:
                    df_clean = df_clean.drop(columns=['adduct_rank'])
                    self._log(f"   ✅ Removed adduct_rank column\n")
            except Exception as e:
                self._log(f"   ❌ Error removing adduct_rank column: {str(e)}\n")
                import traceback
                self._log(f"   Traceback: {traceback.format_exc()}\n")
                raise
            
            # Check for remaining duplicates
            try:
                self._log(f"   Checking for remaining duplicates\n")
                duplicates_remain = df_clean['LipidID'].duplicated().sum()
                
                if duplicates_remain > 0:
                    self._log(f"   ⚠️  Warning: {duplicates_remain} duplicate LipidIDs remain after cleanup\n")
                
                cleaned_count = len(df_clean)
                self._log(f"   ✅ Cleaned: {original_count} → {cleaned_count} lipids\n")
            except Exception as e:
                self._log(f"   ❌ Error checking duplicates: {str(e)}\n")
                import traceback
                self._log(f"   Traceback: {traceback.format_exc()}\n")
                raise
            
            return df_clean
            
        except Exception as e:
            self._log(f"❌ Error in _clean_duplicate_lipids: {str(e)}\n")
            import traceback
            self._log(f"Full traceback: {traceback.format_exc()}\n")
            raise
    
    def _generate_lipid_class_tables(self, pos_df: Optional[pd.DataFrame], neg_df: Optional[pd.DataFrame]) -> tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
        """Aggregate lipid class tables internally from positive/negative lipid DataFrames.

        For each polarity DataFrame:
          - Requires a 'Class' column
          - Identifies sample intensity columns (numeric columns not in feature set)
          - Groups by Class and sums sample columns
          - Preserves Class_name mapping from the source dataframe
          - Drops rows where all summed sample values are zero
          - Returns (pos_class_df, neg_class_df) or None where not applicable
        """
        feature_cols = {
            'LipidID','Class','Class_name','LipidGroup','Charge','CalcMz','BaseRt','SubClass',
            'AdductIon','IonFormula','MolStructure','ObsMz','ObsRt','Polarity'
        }

        def build(df: Optional[pd.DataFrame]) -> Optional[pd.DataFrame]:
            if df is None or df.empty or 'Class' not in df.columns:
                return None
            
            # Sample columns: look for Area[...] columns first
            sample_cols = [c for c in df.columns if isinstance(c, str) and c.startswith('Area[') and c.endswith(']')]
            
            if not sample_cols:
                # Fallback: numeric columns not in features and not Grade
                numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
                sample_cols = [c for c in numeric_cols if c not in feature_cols and not c.startswith('Grade[')]
            
            if not sample_cols:
                return None
            
            # If Class_name exists, preserve it by taking first occurrence per Class
            groupby_cols = ['Class']
            agg_dict = {col: 'sum' for col in sample_cols}
            
            if 'Class_name' in df.columns:
                # Get the first Class_name for each Class group
                agg_dict['Class_name'] = 'first'
            
            grouped = df.groupby('Class', as_index=False).agg(agg_dict)
            
            # Drop rows with all-zero sample values
            mask_all_zero = (grouped[sample_cols].fillna(0) == 0).all(axis=1)
            grouped = grouped.loc[~mask_all_zero]
            
            # Reorder columns: Class, Class_name (if exists), then sample columns
            if 'Class_name' in grouped.columns:
                column_order = ['Class', 'Class_name'] + sample_cols
            else:
                column_order = ['Class'] + sample_cols
            
            grouped = grouped[column_order]
            
            self._log(f"   ✅ Grouped {len(df)} → {len(grouped)} lipids by class\n")
            return grouped if not grouped.empty else None

        return build(pos_df), build(neg_df)

    def _apply_grade_filtering(
        self,
        df: pd.DataFrame,
        grade_threshold: str,
        polarity: str,
        filter_rej: bool = False,
        remove_rows_by_grade: bool = False,
        min_keep_percentage: float = 80.0,
    ) -> pd.DataFrame:
        """
        Apply grade-based filtering using one of two modes:
        1) Default mode: set Area values to 0 for poor quality grades.
        2) Alternative mode: remove full rows when acceptable-grade coverage is below threshold.
        Optionally filter by Rej column to remove rejected rows.
        
        Args:
            df: DataFrame with Area and Grade columns
            grade_threshold: Threshold string (e.g., 'A,B,C'). Empty string skips grade filtering.
            polarity: 'positive' or 'negative' for logging
            filter_rej: If True, keep only rows where Rej column is False
            remove_rows_by_grade: If True, remove rows below acceptable-grade coverage threshold
            min_keep_percentage: Minimum percentage of mapped samples per row with acceptable grades
            
        Returns:
            DataFrame with grade filtering applied and/or Rej filtering
        """
        try:
            # Check if we should skip both filters
            skip_grade = not grade_threshold or grade_threshold.strip() == ''
            
            if skip_grade and not filter_rej:
                self._log(f"   No grade threshold specified and Rej filtering disabled - skipping quality filtering\n")
                return df
            
            # Apply Rej column filtering first (removes entire rows)
            if filter_rej:
                if 'Rej' in df.columns:
                    initial_rows = len(df)
                    # Keep only rows where Rej is False (remove True)
                    # Handle various representations: 'False', False, 'false', 0, etc.
                    df = df[df['Rej'].astype(str).str.upper() != 'TRUE']
                    removed_rows = initial_rows - len(df)
                    self._log(f"\n🔍 Rej Column Filtering ({polarity})...\n")
                    self._log(f"   Initial rows: {initial_rows:,}\n")
                    self._log(f"   Removed rows (Rej=True): {removed_rows:,}\n")
                    self._log(f"   Remaining rows: {len(df):,}\n\n")
                else:
                    self._log(f"   ⚠️ Rej column filtering enabled but 'Rej' column not found - skipping Rej filtering\n")
            
            # Skip grade filtering if no threshold specified
            if skip_grade:
                self._log(f"   No grade threshold specified - skipping grade filtering\n")
                return df
            
            # Find all Grade columns (pattern: Grade[sample_name])
            grade_cols = [col for col in df.columns if isinstance(col, str) and col.startswith('Grade[') and col.endswith(']')]
            
            if not grade_cols:
                self._log(f"   No Grade columns found for {polarity} data - skipping grade filtering\n")
                return df
            
            # Find all sample columns
            # After sample mapping, Area columns are named: Area[sample_name]
            # So sample columns are either:
            # 1. Columns starting with 'Area[' and ending with ']'
            # 2. Numeric columns that are not features and not Grade
            feature_cols = {
                'LipidID','Class','Class_name','LipidGroup','Charge','CalcMz','BaseRt','SubClass',
                'AdductIon','IonFormula','MolStructure','ObsMz','ObsRt','Polarity', 'TotalScore', 
                'FAKey', 'TotalGrade'
            }
            
            # Find Area columns (these are the sample intensity columns after mapping)
            area_cols = [c for c in df.columns if isinstance(c, str) and c.startswith('Area[') and c.endswith(']')]
            
            if not area_cols:
                # Fallback: use numeric columns not in features or Grade
                numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
                sample_cols = [c for c in numeric_cols if c not in feature_cols and not c.startswith('Grade[')]
            else:
                sample_cols = area_cols
            
            self._log(f"\n📊 Step: Grade-based Quality Filtering ({polarity})...\n")
            self._log(f"   Starting with {len(df)} lipid rows\n")
            self._log(f"   Found {len(grade_cols)} Grade columns and {len(sample_cols)} sample columns\n")
            
            # Build mapping from Grade columns to sample columns
            # Extract sample names from Grade[sample] and match to Area[sample]
            grade_to_sample_map = {}
            
            # If using Area columns, extract sample names from both Grade and Area
            if area_cols:
                # Create mapping from Grade[sample] -> Area[sample]
                for grade_col in grade_cols:
                    # Extract sample name: Grade[DU145 3D 1] -> DU145 3D 1
                    if '[' in grade_col and ']' in grade_col:
                        sample_name = grade_col[grade_col.index('[')+1:grade_col.rindex(']')]
                        # Find matching Area column
                        matching_area = f"Area[{sample_name}]"
                        if matching_area in sample_cols:
                            grade_to_sample_map[grade_col] = matching_area
                        else:
                            self._log(f"   ⚠️ No matching sample column found for {grade_col}\n")
            else:
                # Fallback: match by sample name extraction
                for grade_col in grade_cols:
                    if '[' in grade_col and ']' in grade_col:
                        sample_name = grade_col[grade_col.index('[')+1:grade_col.rindex(']')]
                        matching_samples = [col for col in sample_cols if sample_name == col or sample_name in col or col in sample_name]
                        if matching_samples:
                            grade_to_sample_map[grade_col] = matching_samples[0]
                        else:
                            self._log(f"   ⚠️ No matching sample column found for {grade_col}\n")
            
            # Report mapping coverage; continue with mapped subset even if extras are present
            if not grade_to_sample_map:
                self._log(f"\n⚠️ No Grade→Sample mappings could be resolved; skipping grade filtering\n\n")
                return df
            unmatched_grades = [g for g in grade_cols if g not in grade_to_sample_map]
            unmatched_samples = [s for s in sample_cols if s not in grade_to_sample_map.values()]
            self._log(f"   ✓ Mapped {len(grade_to_sample_map)}/{len(grade_cols)} Grade→Sample pairs\n")
            if unmatched_grades:
                self._log(f"   ⚠️ Unmatched Grade columns (ignored): {', '.join(unmatched_grades[:5])}{'...' if len(unmatched_grades) > 5 else ''}\n")
            if unmatched_samples:
                self._log(f"   ⚠️ Sample columns without Grade (left unfiltered): {', '.join(unmatched_samples[:5])}{'...' if len(unmatched_samples) > 5 else ''}\n")
            
            # Parse threshold - comma-separated grades (e.g., "A,B,C" or "A,B,C,D" or "A,B,C,P")
            grade_threshold_clean = grade_threshold.strip().upper()
            if ',' in grade_threshold_clean:
                # Comma-separated format: "A,B,C"
                acceptable_grades = [g.strip() for g in grade_threshold_clean.split(',') if g.strip()]
            else:
                # Fallback for legacy format: if user enters just grades without comma
                acceptable_grades = [c for c in grade_threshold_clean if c.isalpha()]
            
            # Remove duplicates and sort
            acceptable_grades = sorted(list(set(acceptable_grades)))
            
            if not acceptable_grades:
                # Default if parsing failed
                acceptable_grades = ['A', 'B', 'C']
                self._log(f"   Warning: Could not parse grade threshold '{grade_threshold}', using default: A, B, C\n")
            
            # Parse and validate row-coverage threshold
            try:
                min_keep_percentage = float(str(min_keep_percentage).replace('%', '').strip())
            except Exception:
                min_keep_percentage = 80.0
            if min_keep_percentage < 0:
                min_keep_percentage = 0.0
            if min_keep_percentage > 100:
                min_keep_percentage = 100.0

            self._log(f"   Grade threshold input: {grade_threshold}\n")
            self._log(f"   Acceptable grades: {', '.join(acceptable_grades)}\n")
            if remove_rows_by_grade:
                self._log(f"   Mode: REMOVE ROWS below {min_keep_percentage:.1f}% acceptable-grade coverage\n\n")
            else:
                self._log(f"   Mode: SET AREA TO 0 for poor grades\n")
                self._log(f"   Poor grades (will be set to 0): All other grades\n\n")
            
            # DEBUG: Show sample grade values
            self._log(f"   📋 DEBUG: Inspecting Grade column values...\n")
            for grade_col in list(grade_to_sample_map.keys())[:3]:  # Show first 3
                sample_values = df[grade_col].unique()[:5]
                self._log(f"      {grade_col}: {sample_values}\n")
            
            if remove_rows_by_grade:
                # Row-level mode: keep rows meeting the acceptable-grade percentage threshold.
                total_pairs = len(grade_to_sample_map)
                if total_pairs == 0:
                    self._log("   ⚠️ No Grade→Sample pairs available for row-level filtering; skipping\n")
                    return df

                acceptable_counts = pd.Series(0, index=df.index, dtype='int64')
                for grade_col in grade_to_sample_map.keys():
                    acceptable_mask = df[grade_col].astype(str).str.upper().isin(acceptable_grades)
                    acceptable_counts = acceptable_counts + acceptable_mask.astype('int64')

                acceptable_pct = (acceptable_counts / float(total_pairs)) * 100.0
                keep_mask = acceptable_pct >= min_keep_percentage

                rows_before = len(df)
                rows_kept = int(keep_mask.sum())
                rows_removed = rows_before - rows_kept

                self._log("   📉 Row-Level Grade Coverage Results:\n")
                self._log(f"      Total rows before: {rows_before:,}\n")
                self._log(f"      Required acceptable-grade coverage: {min_keep_percentage:.1f}%\n")
                self._log(f"      Grade/sample pairs considered per row: {total_pairs}\n")
                self._log(f"      Rows kept: {rows_kept:,}\n")
                self._log(f"      Rows removed: {rows_removed:,}\n")

                if rows_before > 0:
                    self._log(f"      Retention: {(rows_kept / rows_before) * 100:.1f}%\n")

                df = df.loc[keep_mask].copy()
            else:
                # Value-level mode (default): set Area values to zero for poor grades.
                rows_with_nonzero_before = (df[sample_cols] != 0).any(axis=1).sum()
                total_nonzero_values_before = (df[sample_cols] != 0).sum().sum()

                total_changed = 0
                sample_stats = []

                for grade_col, sample_col in grade_to_sample_map.items():
                    poor_grade_mask = ~df[grade_col].astype(str).str.upper().isin(acceptable_grades)
                    values_to_change = ((df[sample_col] != 0) & poor_grade_mask).sum()

                    if values_to_change > 0:
                        df.loc[poor_grade_mask, sample_col] = 0
                        total_changed += values_to_change
                        sample_stats.append(f"      • {sample_col}: {values_to_change} values → 0")

                total_nonzero_values_after = (df[sample_cols] != 0).sum().sum()
                rows_with_all_zeros_after = (df[sample_cols] == 0).all(axis=1).sum()
                rows_still_nonzero = (df[sample_cols] != 0).any(axis=1).sum()

                self._log(f"   📉 Filtering Results:\n")
                self._log(f"      Total Area values set to 0 (due to poor grades): {total_changed:,}\n")
                self._log(f"      Non-zero values before: {total_nonzero_values_before:,}\n")
                self._log(f"      Non-zero values after:  {total_nonzero_values_after:,}\n")
                if total_nonzero_values_before > 0:
                    reduction_pct = ((total_nonzero_values_before - total_nonzero_values_after) / total_nonzero_values_before * 100)
                else:
                    reduction_pct = 0.0
                self._log(f"      Reduction: {total_nonzero_values_before - total_nonzero_values_after:,} values ({reduction_pct:.1f}%)\n\n")

                self._log(f"   📊 Row Impact:\n")
                self._log(f"      Rows with non-zero values before: {rows_with_nonzero_before:,}\n")
                self._log(f"      Rows with non-zero values after:  {rows_still_nonzero:,}\n")
                self._log(f"      Rows now all-zero (will be removed later): {rows_with_all_zeros_after:,}\n\n")

                if sample_stats:
                    self._log(f"   📝 Per-Sample Breakdown:\n")
                    for stat in sample_stats[:10]:
                        self._log(f"{stat}\n")
                    if len(sample_stats) > 10:
                        self._log(f"      ... and {len(sample_stats) - 10} more samples\n")
            
            self._log(f"\n✅ Grade filtering complete!\n\n")
            
        except Exception as e:
            self._log(f"\n⚠️ Error during grade filtering: {e}\n")
            import traceback
            self._log(f"{traceback.format_exc()}\n")
        
        return df
    
    def _drop_all_zero_rows(self, df: Optional[pd.DataFrame]) -> Optional[pd.DataFrame]:
        """Remove rows where ALL numeric columns are zero (ignoring non-numeric columns).

        This is applied before returning cleaned DataFrames to avoid retaining
        empty signal rows.
        """
        if df is None or df.empty:
            return df
        try:
            # Identify feature (non-intensity) columns to exclude from zero check
            feature_cols = {
                'LipidID','Class','Class_name','LipidGroup','Charge','CalcMz','BaseRt','SubClass',
                'AdductIon','IonFormula','MolStructure','ObsMz','ObsRt','Polarity'
            }
            # Numeric columns
            numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
            # Sample intensity columns = numeric columns not in feature_cols
            sample_cols = [c for c in numeric_cols if c not in feature_cols]
            if not sample_cols:
                return df  # nothing to filter on
            
            rows_before = len(df)
            mask_all_zero = (df[sample_cols].fillna(0) == 0).all(axis=1)
            removed = int(mask_all_zero.sum())
            
            # Only log if we actually removed rows (silent otherwise)
            # Note: This removes ONLY rows where ALL Area columns are zero
            
            return df.loc[~mask_all_zero].copy()
        except Exception as e:
            self._log(f"   ⚠️ Error in _drop_all_zero_rows: {e}\n")
            import traceback
            self._log(f"{traceback.format_exc()}\n")
            return df

    def _drop_all_zero_columns(self, df: Optional[pd.DataFrame]) -> Optional[pd.DataFrame]:
        """Remove columns where ALL rows are zero (for Area/sample columns only).
        
        This removes sample columns that have no signal in any lipid,
        keeping only columns with at least some non-zero values.
        """
        if df is None or df.empty:
            self._log(f"   ⚠️ DataFrame is None or empty\n")
            return df
        try:
            # Identify feature (non-intensity) columns to preserve
            feature_cols = {
                'LipidID','Class','Class_name','LipidGroup','Charge','CalcMz','BaseRt','SubClass',
                'AdductIon','IonFormula','MolStructure','ObsMz','ObsRt','Polarity'
            }
            
            # Sample columns are Area[...] columns
            area_cols = [c for c in df.columns if isinstance(c, str) and c.startswith('Area[') and c.endswith(']')]
            
            self._log(f"      Found {len(area_cols)} Area columns\n")
            
            if not area_cols:
                # Fallback: use numeric columns not in feature_cols
                numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
                area_cols = [c for c in numeric_cols if c not in feature_cols and not c.startswith('Grade[')]
                self._log(f"      Fallback: Found {len(area_cols)} numeric columns\n")
            
            if not area_cols:
                self._log(f"      No Area columns found - skipping analysis\n")
                return df  # nothing to filter on
            
            # Calculate sum for each Area column and identify all-zero columns
            self._log(f"      Analyzing {len(area_cols)} Area columns...\n")
            
            cols_to_drop = []
            cols_low_signal = []
            
            for col in area_cols:
                col_sum = df[col].fillna(0).sum()
                
                # Check if ALL values are zero
                if col_sum == 0:
                    cols_to_drop.append(col)
                # Flag columns with very low signal (< 1000)
                elif col_sum < 1000:
                    cols_low_signal.append((col, col_sum))
            
            # Report low signal columns
            if cols_low_signal:
                self._log(f"      ⚠️ {len(cols_low_signal)} Area columns with sum < 1000 (low signal):\n")
                for col, col_sum in sorted(cols_low_signal, key=lambda x: x[1]):
                    self._log(f"         {col}: {col_sum:.0f}\n")
            
            # Report all-zero columns
            self._log(f"      📊 Area columns before: {len(area_cols)}\n")
            self._log(f"      📊 All-zero Area columns found: {len(cols_to_drop)}\n")
            
            if cols_to_drop:
                cols_after = len(area_cols) - len(cols_to_drop)
                self._log(f"      ✅ Removing {len(cols_to_drop)} all-zero Area columns\n")
                self._log(f"      📊 Area columns after: {cols_after}\n")
                
                if len(cols_to_drop) <= 5:
                    for col in cols_to_drop:
                        self._log(f"         - {col}\n")
                else:
                    for col in cols_to_drop[:5]:
                        self._log(f"         - {col}\n")
                    self._log(f"         ... and {len(cols_to_drop) - 5} more\n")
                
                return df.drop(columns=cols_to_drop)
            else:
                self._log(f"      ✅ No all-zero Area columns found - keeping all {len(area_cols)} columns\n")
                return df
        except Exception as e:
            self._log(f"      ⚠️ Error in all-zero column analysis: {e}\n")
            import traceback
            self._log(f"      {traceback.format_exc()}\n")
            return df
            self._log(f"   ⚠️ Error removing all-zero columns: {e}\n")
            return df

    def _save_lipid_data(self, pos_df: pd.DataFrame, neg_df: pd.DataFrame, output_file: str,
                          combined_df: Optional[pd.DataFrame] = None,
                          pos_class_df: Optional[pd.DataFrame] = None,
                          neg_class_df: Optional[pd.DataFrame] = None,
                          sample_mapping: Optional[Dict[str, str]] = None,
                          pos_sample_mapping: Optional[Dict[str, str]] = None,
                          neg_sample_mapping: Optional[Dict[str, str]] = None,
                          column_clean_patterns: Optional[List[str]] = None):
        """Save cleaned lipid (and class) data to an Excel file.
        
        Applies File ID mapping and column cleaning patterns at save time if provided.

        Sheets written (in order):
          1. Combined_Lipids (if provided and non-empty; columns: LipidID + Class)
          2. Positive_Lipids
          3. Negative_Lipids
          4. Positive_Lipid_Class (optional)
          5. Negative_Lipid_Class (optional)
        """
        self._log(f"💾 Saving lipid data to: {output_file}\n")
        self._log(f"   🔍 DEBUG - SAVE FUNCTION PARAMETERS:\n")
        self._log(f"      sample_mapping type: {type(sample_mapping)}\n")
        self._log(f"      sample_mapping value: {sample_mapping}\n")
        self._log(f"      pos_sample_mapping type: {type(pos_sample_mapping)}\n")
        self._log(f"      pos_sample_mapping value: {pos_sample_mapping}\n")
        self._log(f"      neg_sample_mapping type: {type(neg_sample_mapping)}\n")
        self._log(f"      neg_sample_mapping value: {neg_sample_mapping}\n")
        self._log(f"      column_clean_patterns: {column_clean_patterns}\n")

        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Check if at least one sheet will be written
        has_valid_sheet = (
            (isinstance(combined_df, pd.DataFrame) and not combined_df.empty) or
            (isinstance(pos_df, pd.DataFrame) and not pos_df.empty) or
            (isinstance(neg_df, pd.DataFrame) and not neg_df.empty) or
            (isinstance(pos_class_df, pd.DataFrame) and not pos_class_df.empty) or
            (isinstance(neg_class_df, pd.DataFrame) and not neg_class_df.empty)
        )
        
        if not has_valid_sheet:
            self._log("   ❌ ERROR: No valid data to save (all dataframes are empty)\n")
            return
        
        # Apply File ID mapping and column cleaning patterns if provided
        # This is the FINAL step before saving - simple text replacement in column names
        self._log(f"\n📝 ========== FINAL COLUMN NAME TRANSFORMATIONS ==========\n")
        self._log(f"   DEBUG: sample_mapping type: {type(sample_mapping)}\n")
        self._log(f"   DEBUG: sample_mapping value: {sample_mapping}\n")
        self._log(f"   DEBUG: column_clean_patterns: {column_clean_patterns}\n")
        
        if sample_mapping or pos_sample_mapping or neg_sample_mapping or column_clean_patterns:
            self._log(f"\n   ✅ At least one transformation is enabled\n")
            if sample_mapping:
                self._log(f"   🔄 FILE ID MAPPING ENABLED: {len(sample_mapping)} replacements\n")
                self._log(f"   🔄 Full mapping dictionary:\n")
                for data_id, file_name in sample_mapping.items():
                    self._log(f"      {data_id} : {file_name}\n")

                # Detect one-to-many collisions in mapping values (different DataIDs -> same sample name)
                value_to_keys = {}
                for data_id, mapped_name in sample_mapping.items():
                    value_to_keys.setdefault(str(mapped_name), []).append(str(data_id))
                mapping_collisions = {
                    mapped_name: data_ids
                    for mapped_name, data_ids in value_to_keys.items()
                    if len(data_ids) > 1
                }
                if mapping_collisions:
                    self._log("   ❌ FILE ID discrepancy detected: duplicate mapped sample names in File ID mapping\n")
                    for mapped_name, data_ids in list(mapping_collisions.items())[:10]:
                        self._log(f"      • {mapped_name} <= {', '.join(data_ids)}\n")
                    if len(mapping_collisions) > 10:
                        self._log(f"      ... and {len(mapping_collisions) - 10} more collisions\n")
                    raise ValueError(
                        "File ID discrepancy detected: multiple DataIDs map to the same sample name. "
                        "Please ensure each DataID maps to a unique sample name in your File ID mapping file."
                    )
            elif pos_sample_mapping or neg_sample_mapping:
                self._log("   🔄 FILE ID MAPPING ENABLED (separate mode):\n")
                self._log(f"      Positive mapping entries: {len(pos_sample_mapping) if pos_sample_mapping else 0}\n")
                self._log(f"      Negative mapping entries: {len(neg_sample_mapping) if neg_sample_mapping else 0}\n")
            else:
                self._log(f"   ⚠️ NO FILE ID MAPPING (sample_mapping is None or empty)\n")
            
            if column_clean_patterns:
                self._log(f"   🧹 COLUMN CLEANING ENABLED: {column_clean_patterns}\n")
            else:
                self._log(f"   ⚠️ NO COLUMN CLEANING (patterns is None or empty)\n")
            
            def apply_mapping_and_cleaning_to_df(df, current_mapping=None):
                """
                Apply column name transformations with EXACT word matching:
                1. FIRST: Apply column cleaning patterns to extract pure DataID
                2. THEN: Check if cleaned name EXACTLY matches a DataID in mapping
                3. If exact match: Replace with File Name
                4. If no match: Keep cleaned name and warn user
                
                Example: Area[S1-1] → (clean) → S1-1 → (exact match) → LC_HCtl_1
                """
                if df is None or df.empty:
                    self._log(f"   ⏭️ Skipping empty/None dataframe\n")
                    return df

                active_mapping = current_mapping if current_mapping is not None else sample_mapping
                
                self._log(f"\n   🔍 Processing dataframe with {len(df.columns)} columns:\n")
                self._log(f"      Sample columns: {list(df.columns)[:5]}...\n")
                
                rename_map = {}
                unmatched_columns = []

                def clean_name(raw_name):
                    """Apply user-provided text cleanup patterns to a column/sample name."""
                    name = str(raw_name)
                    if column_clean_patterns:
                        for pattern in column_clean_patterns:
                            if pattern:
                                name = name.replace(pattern, '')
                    name = ' '.join(name.split()).strip()
                    name = name.strip('_').strip()
                    return name
                
                for col in df.columns:
                    col_str = str(col)
                    new_name = col_str
                    
                    # STEP 1: Apply column cleaning patterns FIRST to extract pure DataID
                    if column_clean_patterns:
                        new_name = clean_name(new_name)
                    
                    # STEP 2: Now check if cleaned name EXACTLY matches a DataID in active mapping
                    if active_mapping and new_name in active_mapping:
                        file_name = active_mapping[new_name]
                        final_name = clean_name(file_name) if column_clean_patterns else str(file_name)
                        self._log(
                            f"      ✅ EXACT MATCH: {col_str} → (clean) → {new_name} "
                            f"→ (map) → {file_name} → (final) → {final_name}\n"
                        )
                        rename_map[col] = final_name
                    elif column_clean_patterns and new_name != col_str:
                        # Column was cleaned but no mapping match found
                        self._log(f"      🧹 CLEANED ONLY: {col_str} → {new_name} (no mapping match)\n")
                        if active_mapping:
                            # Warn that this cleaned column didn't match any mapping
                            unmatched_columns.append((col_str, new_name))
                        rename_map[col] = new_name
                    # else: no cleaning patterns or no change, keep original
                
                # Report unmatched columns that user may need to fix
                if unmatched_columns:
                    self._log(f"\n      ⚠️  WARNING: {len(unmatched_columns)} cleaned column(s) did NOT match any File ID mapping:\n")
                    for orig, cleaned in unmatched_columns[:10]:  # Show first 10
                        self._log(f"         • {orig} → {cleaned} (not found in mapping)\n")
                    if len(unmatched_columns) > 10:
                        self._log(f"         ... and {len(unmatched_columns) - 10} more\n")
                    self._log(f"      💡 TIP: Update 'Column Name Cleaning' patterns to properly extract DataIDs\n")
                    self._log(f"      💡 TIP: Or add these DataIDs to your File ID mapping sheet\n")
                
                if rename_map:
                    self._log(f"   ✅ RENAMED {len(rename_map)} COLUMNS in this sheet\n")
                    df = df.rename(columns=rename_map)
                else:
                    self._log(f"   ℹ️ NO COLUMN RENAMES (no patterns applied or matches found)\n")

                # Validate transformed names to avoid ambiguous duplicate sample columns.
                # Duplicates are usually caused by mismatches between file IDs and sample names.
                str_cols = [str(c) for c in df.columns]
                name_counts = Counter(str_cols)
                duplicate_names = [name for name, count in name_counts.items() if count > 1]
                if duplicate_names:
                    self._log("\n      ❌ File ID / sample name discrepancy detected after column transformation:\n")
                    for dup_name in duplicate_names[:10]:
                        self._log(f"         • Duplicate sample name: {dup_name} (appears {name_counts[dup_name]} times)\n")
                    if len(duplicate_names) > 10:
                        self._log(f"         ... and {len(duplicate_names) - 10} more duplicate names\n")
                    self._log("      💡 Check File ID mapping and column cleaning patterns so each sample name is unique.\n")
                    raise ValueError(
                        "File ID and sample name discrepancy detected: duplicate sample names were produced "
                        "after applying mapping/cleaning. Please verify your File ID mapping sheet and "
                        "column-cleaning patterns."
                    )
                return df
            
            pos_df = apply_mapping_and_cleaning_to_df(pos_df, current_mapping=pos_sample_mapping)
            neg_df = apply_mapping_and_cleaning_to_df(neg_df, current_mapping=neg_sample_mapping)
            pos_class_df = apply_mapping_and_cleaning_to_df(pos_class_df, current_mapping=pos_sample_mapping)
            neg_class_df = apply_mapping_and_cleaning_to_df(neg_class_df, current_mapping=neg_sample_mapping)
            combined_df = apply_mapping_and_cleaning_to_df(combined_df, current_mapping=sample_mapping)
            
            self._log(f"   ✅ All column transformations complete\n")

        # Post-mapping cleanup for lipid exports:
        # - If both positive and negative sample columns are present, keep only
        #   the correct polarity for each sheet based on '_pos'/'_neg' suffix.
        # - Strip '_pos'/'_neg' suffixes for cleaner final column names.
        # - If duplicate sample column names remain (e.g., because of mapping
        #   collisions), collapse numeric duplicates by summing.
        def _drop_columns_with_suffix(df: Optional[pd.DataFrame], suffix_to_drop: str, label: str) -> Optional[pd.DataFrame]:
            if df is None or df.empty:
                return df
            suffix = suffix_to_drop.lower()
            cols_to_drop = [
                c for c in df.columns
                if isinstance(c, str) and c.lower().endswith(suffix)
            ]
            if cols_to_drop:
                self._log(f"   🧽 Dropping {len(cols_to_drop)} columns ending with '{suffix_to_drop}' ({label})\n")
                return df.drop(columns=cols_to_drop)
            return df

        def _strip_suffix(df: Optional[pd.DataFrame], suffix_to_strip: str) -> Optional[pd.DataFrame]:
            if df is None or df.empty:
                return df
            suffix = suffix_to_strip
            new_cols = []
            for c in df.columns:
                if isinstance(c, str) and c.endswith(suffix):
                    new_cols.append(c[: -len(suffix)])
                else:
                    new_cols.append(c)
            df = df.copy()
            df.columns = new_cols
            return df

        def _collapse_duplicate_columns(df: Optional[pd.DataFrame], label: str) -> Optional[pd.DataFrame]:
            if df is None or df.empty:
                return df
            if not getattr(df.columns, 'duplicated', None):
                return df

            dup_mask = df.columns.duplicated(keep=False)
            if not dup_mask.any():
                return df

            dup_names = [c for c in df.columns[dup_mask] if isinstance(c, str)]
            # Preserve order of first occurrence
            seen = set()
            dup_names = [c for c in dup_names if not (c in seen or seen.add(c))]

            self._log(f"   ⚠️ Duplicate column names detected ({label}): {len(dup_names)}\n")

            out = pd.DataFrame(index=df.index)
            for col_name in df.columns:
                if col_name in out.columns:
                    continue
                if col_name in dup_names:
                    # Selecting by duplicate label may return a DataFrame (not Series).
                    selected = df.loc[:, col_name]
                    if isinstance(selected, pd.Series):
                        series_list = [selected]
                    else:
                        series_list = [selected.iloc[:, i] for i in range(selected.shape[1])]

                    # If any of the duplicates are numeric, sum numeric values across them.
                    numeric_summed = None
                    for series_raw in series_list:
                        series = pd.to_numeric(series_raw, errors='coerce')
                        if numeric_summed is None:
                            numeric_summed = series.fillna(0)
                        else:
                            numeric_summed = numeric_summed.add(series.fillna(0), fill_value=0)
                    if numeric_summed is not None:
                        out[col_name] = numeric_summed
                    else:
                        out[col_name] = series_list[0]
                else:
                    out[col_name] = df[col_name]
            return out

        # Keep polarity-specific mapped columns if both exist
        pos_df = _drop_columns_with_suffix(pos_df, '_neg', label='POS')
        neg_df = _drop_columns_with_suffix(neg_df, '_pos', label='NEG')
        pos_class_df = _drop_columns_with_suffix(pos_class_df, '_neg', label='POS_CLASS')
        neg_class_df = _drop_columns_with_suffix(neg_class_df, '_pos', label='NEG_CLASS')

        # Strip polarity suffixes for final sheet readability
        pos_df = _strip_suffix(pos_df, '_pos')
        neg_df = _strip_suffix(neg_df, '_neg')
        pos_class_df = _strip_suffix(pos_class_df, '_pos')
        neg_class_df = _strip_suffix(neg_class_df, '_neg')

        # Collapse any duplicate sample columns caused by mapping collisions
        pos_df = _collapse_duplicate_columns(pos_df, label='POS')
        neg_df = _collapse_duplicate_columns(neg_df, label='NEG')
        pos_class_df = _collapse_duplicate_columns(pos_class_df, label='POS_CLASS')
        neg_class_df = _collapse_duplicate_columns(neg_class_df, label='NEG_CLASS')

        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Combined sheet first for downstream ID annotation
            if isinstance(combined_df, pd.DataFrame) and not combined_df.empty:
                # Remove Grade columns before saving
                grade_cols = [c for c in combined_df.columns if isinstance(c, str) and 'Grade[' in c]
                if grade_cols:
                    combined_df = combined_df.drop(columns=grade_cols)
                    self._log(f"   🗑️ Removed {len(grade_cols)} Grade columns from combined data\n")
                combined_df.to_excel(writer, sheet_name='Combined_Lipids', index=False)
                self._log(f"   ✅ Saved combined lipid sheet ({len(combined_df)} rows)\n")
            elif combined_df is not None:
                self._log("   ⚠️ Combined lipid dataframe provided but empty – sheet omitted\n")

            if isinstance(pos_df, pd.DataFrame) and not pos_df.empty:
                # Final pass: remove all-zero Area columns before saving
                before_cols = len(pos_df.columns)
                pos_df = self._drop_all_zero_columns(pos_df)
                after_cols = len(pos_df.columns)
                if after_cols != before_cols:
                    self._log(f"   🧹 Final column cleanup (POS): {before_cols} → {after_cols} (removed all-zero Area columns)\n")

                # Remove Grade columns before saving
                grade_cols = [c for c in pos_df.columns if isinstance(c, str) and 'Grade[' in c]
                if grade_cols:
                    pos_df = pos_df.drop(columns=grade_cols)
                    self._log(f"   🗑️ Removed {len(grade_cols)} Grade columns from positive data\n")
                pos_df.to_excel(writer, sheet_name='Positive_Lipids', index=False)
                self._log(f"   ✅ Saved {len(pos_df)} positive lipids ({len(pos_df.columns)} cols)\n")
            else:
                self._log("   ⚠️ No positive lipids to save\n")

            if isinstance(neg_df, pd.DataFrame) and not neg_df.empty:
                # Final pass: remove all-zero Area columns before saving
                before_cols = len(neg_df.columns)
                neg_df = self._drop_all_zero_columns(neg_df)
                after_cols = len(neg_df.columns)
                if after_cols != before_cols:
                    self._log(f"   🧹 Final column cleanup (NEG): {before_cols} → {after_cols} (removed all-zero Area columns)\n")

                # Remove Grade columns before saving
                grade_cols = [c for c in neg_df.columns if isinstance(c, str) and 'Grade[' in c]
                if grade_cols:
                    neg_df = neg_df.drop(columns=grade_cols)
                    self._log(f"   🗑️ Removed {len(grade_cols)} Grade columns from negative data\n")
                neg_df.to_excel(writer, sheet_name='Negative_Lipids', index=False)
                self._log(f"   ✅ Saved {len(neg_df)} negative lipids ({len(neg_df.columns)} cols)\n")
            else:
                self._log("   ⚠️ No negative lipids to save\n")

            if isinstance(pos_class_df, pd.DataFrame) and not pos_class_df.empty:
                pos_class_df.to_excel(writer, sheet_name='Positive_Lipid_Class', index=False)
                self._log(f"   ✅ Saved positive class sheet ({len(pos_class_df)} rows)\n")
            elif pos_class_df is not None:
                self._log("   ⚠️ Positive class provided but empty – sheet omitted\n")

            if isinstance(neg_class_df, pd.DataFrame) and not neg_class_df.empty:
                neg_class_df.to_excel(writer, sheet_name='Negative_Lipid_Class', index=False)
                self._log(f"   ✅ Saved negative class sheet ({len(neg_class_df)} rows)\n")
            elif neg_class_df is not None:
                self._log("   ⚠️ Negative class provided but empty – sheet omitted\n")

        self._log("   ✅ File saved successfully!\n")


def main():
    """Main function for command line usage"""
    import argparse

    parser = argparse.ArgumentParser(description='Clean metabolite data with pathways')
    parser.add_argument('input_file', help='Input Excel file path')
    parser.add_argument('output_file', help='Output Excel file path')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')

    args = parser.parse_args()

    # Create cleaner instance
    cleaner = CleanColumnDataCleaner()

    # For backward compatibility, create a simple config
    config = {
        'cleaning_mode': 'combined',
        'files': {'combined_file': args.input_file},
        'output': {
            'directory': os.path.dirname(args.output_file) or os.getcwd(),
            'combined_filename': os.path.splitext(os.path.basename(args.output_file))[0]
        }
    }

    # Clean the data
    results = cleaner.clean_data(config)

    if results['success']:
        print(f"\n✅ Success! Cleaned data saved to: {args.output_file}")
        if 'combined_result' in results and results['combined_result'] is not None:
            print(f"   Rows: {len(results['combined_result'])}")
        print(f"   Files created: {', '.join(results['files_created'])}")
    else:
        print(f"\n❌ Failed: {results['message']}")
        exit(1)


if __name__ == "__main__":
    main()