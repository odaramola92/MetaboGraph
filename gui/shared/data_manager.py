"""Centralized data manager for sharing data between tabs"""
from typing import Optional, Dict, Any, List
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class DataManager:
    """
    Centralized data storage for all tabs.
    
    This class manages all shared data that needs to persist between tabs:
    - DataFrames (raw data, cleaned data, annotated data, pathway data)
    - File paths (for auto-loading between tabs)
    - Generic memory store (for any additional data)
    
    Usage:
        data_manager = DataManager()
        data_manager.set_cleaned_data(df)
        
        # Later, in another tab:
        df = data_manager.get_cleaned_data()
    """
    
    def __init__(self):
        """Initialize data storage"""
        logger.info("Initializing DataManager")
        
        # ========== DataFrames ==========
        self.current_data: Optional[pd.DataFrame] = None
        self.cleaned_data: Optional[pd.DataFrame] = None
        self.annotated_data: Optional[pd.DataFrame] = None
        self.pathway_data: Optional[pd.DataFrame] = None
        
        # ========== File Paths (for auto-loading between tabs) ==========
        self.cleaned_excel_path: Optional[str] = None
        self.id_annotated_excel_path: Optional[str] = None
        self.pathway_annotated_excel_path: Optional[str] = None
        self.cleaned_lipid_excel_path: Optional[str] = None
        
        # ========== Lipid-specific data ==========
        self.lipid_input_file: Optional[str] = None
        self.lipid_output_file: Optional[str] = None
        self.lipid_id_annotation_file: Optional[str] = None
        
        # ========== Generic memory store ==========
        # For any additional data not covered by specific attributes
        self.memory_store: Dict[str, Any] = {}
        
        # ========== Processing state ==========
        self.is_processing = False
        self.current_operation: Optional[str] = None
    
    # ========== Current Data ==========
    
    def set_current_data(self, df: pd.DataFrame):
        """Store current/raw data"""
        self.current_data = df
        logger.info(f"Stored current data: {df.shape[0]} rows × {df.shape[1]} cols")
    
    def get_current_data(self) -> Optional[pd.DataFrame]:
        """Retrieve current/raw data"""
        return self.current_data
    
    def has_current_data(self) -> bool:
        """Check if current data is available"""
        return self.current_data is not None
    
    # ========== Cleaned Data ==========
    
    def set_cleaned_data(self, df: pd.DataFrame):
        """Store cleaned data"""
        self.cleaned_data = df
        logger.info(f"Stored cleaned data: {df.shape[0]} rows × {df.shape[1]} cols")
    
    def get_cleaned_data(self) -> Optional[pd.DataFrame]:
        """Retrieve cleaned data"""
        return self.cleaned_data
    
    def has_cleaned_data(self) -> bool:
        """Check if cleaned data is available"""
        return self.cleaned_data is not None
    
    # ========== Annotated Data ==========
    
    def set_annotated_data(self, df: pd.DataFrame):
        """Store ID-annotated data"""
        self.annotated_data = df
        logger.info(f"Stored annotated data: {df.shape[0]} rows × {df.shape[1]} cols")
    
    def get_annotated_data(self) -> Optional[pd.DataFrame]:
        """Retrieve ID-annotated data"""
        return self.annotated_data
    
    def has_annotated_data(self) -> bool:
        """Check if annotated data is available"""
        return self.annotated_data is not None
    
    # ========== Pathway Data ==========
    
    def set_pathway_data(self, df: pd.DataFrame):
        """Store pathway-annotated data"""
        self.pathway_data = df
        logger.info(f"Stored pathway data: {df.shape[0]} rows × {df.shape[1]} cols")
    
    def get_pathway_data(self) -> Optional[pd.DataFrame]:
        """Retrieve pathway-annotated data"""
        return self.pathway_data
    
    def has_pathway_data(self) -> bool:
        """Check if pathway data is available"""
        return self.pathway_data is not None
    
    # ========== File Paths ==========
    
    def set_cleaned_excel_path(self, path: str):
        """Store path to cleaned data Excel file"""
        self.cleaned_excel_path = path
        logger.info(f"Stored cleaned Excel path: {path}")
    
    def get_cleaned_excel_path(self) -> Optional[str]:
        """Get path to cleaned data Excel file"""
        return self.cleaned_excel_path
    
    def set_id_annotated_excel_path(self, path: str):
        """Store path to ID-annotated Excel file"""
        self.id_annotated_excel_path = path
        logger.info(f"Stored ID annotated Excel path: {path}")
    
    def get_id_annotated_excel_path(self) -> Optional[str]:
        """Get path to ID-annotated Excel file"""
        return self.id_annotated_excel_path
    
    def set_pathway_annotated_excel_path(self, path: str):
        """Store path to pathway-annotated Excel file"""
        self.pathway_annotated_excel_path = path
        logger.info(f"Stored pathway annotated Excel path: {path}")
    
    def get_pathway_annotated_excel_path(self) -> Optional[str]:
        """Get path to pathway-annotated Excel file"""
        return self.pathway_annotated_excel_path
    
    def set_cleaned_lipid_excel_path(self, path: str):
        """Store path to cleaned lipid data Excel file"""
        self.cleaned_lipid_excel_path = path
        logger.info(f"Stored cleaned lipid Excel path: {path}")
    
    def get_cleaned_lipid_excel_path(self) -> Optional[str]:
        """Get path to cleaned lipid data Excel file"""
        return self.cleaned_lipid_excel_path
    
    # ========== Lipid-specific Data ==========
    
    def set_lipid_input_file(self, path: str):
        """Store lipid input file path"""
        self.lipid_input_file = path
    
    def get_lipid_input_file(self) -> Optional[str]:
        """Get lipid input file path"""
        return self.lipid_input_file
    
    def set_lipid_output_file(self, path: str):
        """Store lipid output file path"""
        self.lipid_output_file = path
    
    def get_lipid_output_file(self) -> Optional[str]:
        """Get lipid output file path"""
        return self.lipid_output_file
    
    def set_lipid_id_annotation_file(self, path: str):
        """Store lipid ID annotation file path"""
        self.lipid_id_annotation_file = path
    
    def get_lipid_id_annotation_file(self) -> Optional[str]:
        """Get lipid ID annotation file path"""
        return self.lipid_id_annotation_file
    
    # ========== Generic Memory Store ==========
    
    def set_value(self, key: str, value: Any):
        """Store arbitrary value in memory store"""
        self.memory_store[key] = value
        logger.debug(f"Stored memory value: {key}")
    
    def get_value(self, key: str, default=None) -> Any:
        """Retrieve arbitrary value from memory store"""
        return self.memory_store.get(key, default)
    
    def has_value(self, key: str) -> bool:
        """Check if value exists in memory store"""
        return key in self.memory_store
    
    def clear_value(self, key: str):
        """Remove value from memory store"""
        if key in self.memory_store:
            del self.memory_store[key]
            logger.debug(f"Cleared memory value: {key}")
    
    def get_all_values(self) -> Dict[str, Any]:
        """Get all values from memory store"""
        return self.memory_store.copy()
    
    # ========== Processing State ==========
    
    def set_processing(self, is_processing: bool, operation: Optional[str] = None):
        """Set processing state"""
        self.is_processing = is_processing
        self.current_operation = operation
        if is_processing:
            logger.info(f"Processing started: {operation}")
        else:
            logger.info(f"Processing completed: {operation}")
    
    def is_currently_processing(self) -> bool:
        """Check if currently processing"""
        return self.is_processing
    
    def get_current_operation(self) -> Optional[str]:
        """Get current operation name"""
        return self.current_operation
    
    # ========== Utilities ==========
    
    def clear_all(self):
        """Clear all stored data (use with caution!)"""
        logger.warning("Clearing all data in DataManager")
        self.current_data = None
        self.cleaned_data = None
        self.annotated_data = None
        self.pathway_data = None
        self.cleaned_excel_path = None
        self.id_annotated_excel_path = None
        self.pathway_annotated_excel_path = None
        self.cleaned_lipid_excel_path = None
        self.lipid_input_file = None
        self.lipid_output_file = None
        self.lipid_id_annotation_file = None
        self.memory_store.clear()
    
    def get_status(self) -> Dict[str, Any]:
        """Get status of all stored data"""
        return {
            "has_current_data": self.has_current_data(),
            "has_cleaned_data": self.has_cleaned_data(),
            "has_annotated_data": self.has_annotated_data(),
            "has_pathway_data": self.has_pathway_data(),
            "cleaned_excel_path": self.cleaned_excel_path,
            "id_annotated_excel_path": self.id_annotated_excel_path,
            "pathway_annotated_excel_path": self.pathway_annotated_excel_path,
            "cleaned_lipid_excel_path": self.cleaned_lipid_excel_path,
            "memory_store_keys": list(self.memory_store.keys()),
            "is_processing": self.is_processing,
            "current_operation": self.current_operation,
        }
    
    def print_status(self):
        """Print status of all stored data (useful for debugging)"""
        status = self.get_status()
        print("\n========== DataManager Status ==========")
        for key, value in status.items():
            print(f"{key}: {value}")
        print("=========================================\n")
