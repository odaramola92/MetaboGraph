"""
Factor Mapping Configuration Manager

Handles saving, loading, and managing Two-Way ANOVA factor configurations.
Allows users to define and reuse factor mappings across different datasets.
"""

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class FactorMappingConfig:
    """Manages user-defined factor mappings for Two-Way ANOVA."""
    
    def __init__(self, config_path: str = "factor_mapping_config.json"):
        """Initialize factor mapping configuration manager.
        
        Args:
            config_path: Path to the configuration JSON file
        """
        self.config_path = config_path
        self.config = self._load_config()
    
    def _load_config(self) -> dict:
        """Load configuration from JSON file."""
        if os.path.exists(self.config_path):
            try:
                with open(self.config_path, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError):
                return self._create_default_config()
        return self._create_default_config()
    
    def _create_default_config(self) -> dict:
        """Create default configuration structure."""
        return {
            "version": "1.0",
            "description": "Two-Way ANOVA Factor Configuration Storage",
            "saved_mappings": {}
        }
    
    def _save_config(self):
        """Save configuration to JSON file."""
        try:
            os.makedirs(os.path.dirname(self.config_path) or '.', exist_ok=True)
            with open(self.config_path, 'w') as f:
                json.dump(self.config, f, indent=2)
        except IOError as e:
            print(f"Error saving config: {e}")
    
    def save_mapping(self, name: str, factor_a_name: str, factor_b_name: str,
                    factor_a_groups: Dict[str, str], factor_b_groups: Dict[str, str],
                    description: str = "", notes: str = "") -> bool:
        """Save a factor mapping configuration.
        
        Args:
            name: Identifier for this mapping
            factor_a_name: Display name for Factor A (e.g., "Injury Status")
            factor_b_name: Display name for Factor B (e.g., "Heat Condition")
            factor_a_groups: Dict mapping group names to Factor A levels
            factor_b_groups: Dict mapping group names to Factor B levels
            description: Brief description of the mapping
            notes: Additional notes about the design
            
        Returns:
            True if saved successfully, False otherwise
        """
        try:
            # Validate that all groups are consistent
            groups_a = set(factor_a_groups.keys())
            groups_b = set(factor_b_groups.keys())
            if groups_a != groups_b:
                print("Error: Factor A and B must have same groups")
                return False
            
            # Check if mapping is balanced
            is_balanced = self._check_balanced_design(factor_a_groups, factor_b_groups)
            
            mapping = {
                "description": description,
                "created": datetime.now().isoformat(),
                "last_modified": datetime.now().isoformat(),
                "factorA": {
                    "name": factor_a_name,
                    "groups": factor_a_groups
                },
                "factorB": {
                    "name": factor_b_name,
                    "groups": factor_b_groups
                },
                "design_type": "factorial",
                "balanced": is_balanced,
                "notes": notes
            }
            
            self.config["saved_mappings"][name] = mapping
            self._save_config()
            return True
        except Exception as e:
            print(f"Error saving mapping: {e}")
            return False
    
    def load_mapping(self, name: str) -> Optional[Dict]:
        """Load a saved factor mapping configuration.
        
        Args:
            name: Identifier of the mapping to load
            
        Returns:
            Mapping dict or None if not found
        """
        return self.config["saved_mappings"].get(name)
    
    def list_mappings(self) -> List[str]:
        """Get list of all saved mapping names."""
        return list(self.config["saved_mappings"].keys())
    
    def delete_mapping(self, name: str) -> bool:
        """Delete a saved mapping.
        
        Args:
            name: Identifier of mapping to delete
            
        Returns:
            True if deleted successfully
        """
        if name in self.config["saved_mappings"]:
            del self.config["saved_mappings"][name]
            self._save_config()
            return True
        return False
    
    def update_mapping(self, name: str, factor_a_name: str = None, 
                      factor_b_name: str = None, factor_a_groups: Dict = None,
                      factor_b_groups: Dict = None, notes: str = None) -> bool:
        """Update an existing mapping.
        
        Args:
            name: Identifier of mapping to update
            factor_a_name: New name for Factor A (optional)
            factor_b_name: New name for Factor B (optional)
            factor_a_groups: New Factor A groups mapping (optional)
            factor_b_groups: New Factor B groups mapping (optional)
            notes: New notes (optional)
            
        Returns:
            True if updated successfully
        """
        if name not in self.config["saved_mappings"]:
            return False
        
        mapping = self.config["saved_mappings"][name]
        
        if factor_a_name:
            mapping["factorA"]["name"] = factor_a_name
        if factor_b_name:
            mapping["factorB"]["name"] = factor_b_name
        if factor_a_groups:
            mapping["factorA"]["groups"] = factor_a_groups
        if factor_b_groups:
            mapping["factorB"]["groups"] = factor_b_groups
        if notes:
            mapping["notes"] = notes
        
        mapping["last_modified"] = datetime.now().isoformat()
        mapping["balanced"] = self._check_balanced_design(
            factor_a_groups or mapping["factorA"]["groups"],
            factor_b_groups or mapping["factorB"]["groups"]
        )
        
        self._save_config()
        return True
    
    @staticmethod
    def _check_balanced_design(factor_a_groups: Dict[str, str], 
                              factor_b_groups: Dict[str, str]) -> bool:
        """Check if the design is balanced (equal replicates per cell).
        
        Args:
            factor_a_groups: Factor A group assignments
            factor_b_groups: Factor B group assignments
            
        Returns:
            True if balanced design, False otherwise
        """
        from collections import Counter
        
        # Count occurrences of each (factorA, factorB) combination
        cell_counts = Counter(
            (factor_a_groups[g], factor_b_groups[g]) 
            for g in factor_a_groups.keys()
        )
        
        # Balanced if all cells have same count
        counts = list(cell_counts.values())
        return len(set(counts)) == 1 if counts else False
    
    def get_design_summary(self, factor_a_groups: Dict[str, str],
                          factor_b_groups: Dict[str, str]) -> Dict:
        """Generate design summary including cell counts and balance info.
        
        Args:
            factor_a_groups: Factor A group assignments
            factor_b_groups: Factor B group assignments
            
        Returns:
            Dict with design summary
        """
        from collections import Counter
        
        cell_counts = Counter(
            (factor_a_groups[g], factor_b_groups[g])
            for g in factor_a_groups.keys()
        )
        
        factor_a_levels = set(factor_a_groups.values())
        factor_b_levels = set(factor_b_groups.values())
        
        return {
            "factorA_levels": sorted(factor_a_levels),
            "factorB_levels": sorted(factor_b_levels),
            "cell_counts": dict(cell_counts),
            "total_samples": len(factor_a_groups),
            "balanced": len(set(cell_counts.values())) == 1,
            "min_reps": min(cell_counts.values()) if cell_counts else 0,
            "max_reps": max(cell_counts.values()) if cell_counts else 0,
        }
