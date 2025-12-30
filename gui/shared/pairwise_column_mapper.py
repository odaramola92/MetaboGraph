"""
Pairwise Column Mapper for Visualization Tab

Handles auto-detection and user verification of pairwise stats columns
(log2FC, p-value, fold-change) for selected group pairs.
"""

import tkinter as tk
from tkinter import ttk, messagebox
import logging
import pandas as pd
from typing import Optional, Dict, List, Tuple

logger = logging.getLogger(__name__)


class PairwiseColumnMapper:
    """Maps group pairs to their corresponding stats columns with user verification."""
    
    def __init__(self, df: pd.DataFrame, groups: List[str]):
        """
        Initialize mapper with data and groups.
        
        Args:
            df: Complete metabolite dataframe
            groups: List of group names (e.g., ['Control', 'TBI'])
        """
        self.df = df
        self.groups = groups
        self.mappings: Dict[Tuple[str, str], Dict[str, Optional[str]]] = {}
        # Key format: (g1, g2) -> {'log2FC': 'Control_vs_TBI_log2FC', 'pvalue': 'Control_vs_TBI_adj_p', ...}
    
    def auto_detect_columns(self) -> Dict[Tuple[str, str], Dict[str, Optional[str]]]:
        """
        Auto-detect pairwise stats columns for all group pairs.
        
        Returns dict mapping (g1, g2) -> {'log2FC': col, 'pvalue': col}
        """
        from itertools import combinations
        
        results = {}
        
        for g1, g2 in combinations(self.groups, 2):
            detected = self._detect_for_pair(g1, g2)
            results[(g1, g2)] = detected
        
        return results
    
    def _detect_for_pair(self, g1: str, g2: str) -> Dict[str, Optional[str]]:
        """
        Detect columns for a specific group pair.
        
        Returns: {'log2FC': col or None, 'pvalue': col or None}
        """
        detection = {'log2FC': None, 'pvalue': None}
        
        # Try both orderings
        prefixes_to_try = [f"{g1}_vs_{g2}", f"{g2}_vs_{g1}"]
        
        for prefix in prefixes_to_try:
            # Look for log2FC
            if detection['log2FC'] is None:
                fc_col = f"{prefix}_log2FC"
                if fc_col in self.df.columns:
                    detection['log2FC'] = fc_col
            
            # Look for p-value (prefer adj_p, skip neg_log10)
            if detection['pvalue'] is None:
                for p_col in [f"{prefix}_adj_p", f"{prefix}_p_value", f"{prefix}_pvalue", f"{prefix}_p"]:
                    if p_col in self.df.columns and 'neg_log10' not in p_col.lower():
                        detection['pvalue'] = p_col
                        break
            
            # If both found, stop searching
            if detection['log2FC'] and detection['pvalue']:
                break
        
        return detection
    
    def show_verification_dialog(self, parent_window) -> Optional[Dict]:
        """
        Show GUI dialog for user to verify/adjust auto-detected columns.
        
        Returns: Dict mapping (g1, g2) -> {'log2FC': col, 'pvalue': col}, or None if cancelled
        """
        dialog = PairwiseColumnDialog(parent_window, self.df, self.groups, self.auto_detect_columns())
        result = dialog.wait_for_result()
        if result:
            self.mappings = result
        return result


class PairwiseColumnDialog(tk.Toplevel):
    """Dialog for verifying pairwise column mappings."""
    
    def __init__(self, parent, df: pd.DataFrame, groups: List[str], auto_detected: Dict):
        """
        Initialize dialog.
        
        Args:
            parent: Parent window
            df: Dataframe
            groups: List of groups
            auto_detected: Auto-detected columns from mapper
        """
        super().__init__(parent)
        self.title("Pairwise Column Verification")
        self.geometry("900x500")
        self.resizable(True, True)
        
        self.df = df
        self.groups = groups
        self.auto_detected = auto_detected
        self.result = None
        
        # Get all available *_log2FC and *_adj_p columns for dropdowns
        self.available_fc_cols = sorted([c for c in df.columns if '_log2FC' in c or '_FC' in c])
        self.available_p_cols = sorted([c for c in df.columns if ('_adj_p' in c or '_p_value' in c or '_pvalue' in c) and 'neg_log10' not in c.lower()])
        
        # Build UI
        self._build_ui()
        
        # Make modal
        self.transient(parent)
        self.grab_set()
    
    def _build_ui(self):
        """Build the dialog UI."""
        # Title frame
        title_frame = ttk.Frame(self)
        title_frame.pack(fill='x', padx=10, pady=10)
        
        ttk.Label(title_frame, text="Pairwise Stats Column Mapping", font=('Arial', 12, 'bold')).pack()
        ttk.Label(title_frame, text="Verify or adjust the auto-detected columns for each group comparison:").pack()
        
        # Main content frame with scrollbar
        content_frame = ttk.Frame(self)
        content_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        canvas = tk.Canvas(content_frame)
        scrollbar = ttk.Scrollbar(content_frame, orient='vertical', command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Create rows for each pair
        self.column_vars = {}  # Store OptionMenu vars for retrieval
        
        from itertools import combinations
        pairs = list(combinations(self.groups, 2))
        
        for pair_idx, (g1, g2) in enumerate(pairs):
            pair_key = (g1, g2)
            detected = self.auto_detected.get(pair_key, {})
            
            # Pair label
            pair_label = ttk.Label(scrollable_frame, text=f"{g1} vs {g2}", font=('Arial', 10, 'bold'))
            pair_label.grid(row=pair_idx, column=0, sticky='w', padx=5, pady=10)
            
            # log2FC selector
            fc_label = ttk.Label(scrollable_frame, text="log2FC:")
            fc_label.grid(row=pair_idx, column=1, sticky='e', padx=5)
            
            fc_var = tk.StringVar(value=detected.get('log2FC') or '')
            fc_menu = ttk.Combobox(scrollable_frame, textvariable=fc_var, values=self.available_fc_cols, width=40, state='readonly')
            fc_menu.grid(row=pair_idx, column=2, sticky='ew', padx=5)
            self.column_vars[(pair_key, 'log2FC')] = fc_var
            
            # p-value selector
            p_label = ttk.Label(scrollable_frame, text="p-value:")
            p_label.grid(row=pair_idx, column=3, sticky='e', padx=5)
            
            p_var = tk.StringVar(value=detected.get('pvalue') or '')
            p_menu = ttk.Combobox(scrollable_frame, textvariable=p_var, values=self.available_p_cols, width=40, state='readonly')
            p_menu.grid(row=pair_idx, column=4, sticky='ew', padx=5)
            self.column_vars[(pair_key, 'pvalue')] = p_var
            
            # Status indicator
            if detected.get('log2FC') and detected.get('pvalue'):
                status_label = ttk.Label(scrollable_frame, text="✓ Auto-detected", foreground='green')
            else:
                status_label = ttk.Label(scrollable_frame, text="⚠ Manual selection needed", foreground='orange')
            status_label.grid(row=pair_idx, column=5, padx=5)
        
        # Configure grid weights for responsiveness
        scrollable_frame.columnconfigure(2, weight=1)
        scrollable_frame.columnconfigure(4, weight=1)
        
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')
        
        # Button frame
        button_frame = ttk.Frame(self)
        button_frame.pack(fill='x', padx=10, pady=10)
        
        ttk.Button(button_frame, text="✓ Confirm", command=self._on_confirm).pack(side='left', padx=5)
        ttk.Button(button_frame, text="✗ Cancel", command=self._on_cancel).pack(side='left', padx=5)
    
    def _on_confirm(self):
        """Validate and confirm selections."""
        from itertools import combinations
        
        result = {}
        pairs = list(combinations(self.groups, 2))
        
        for g1, g2 in pairs:
            pair_key = (g1, g2)
            fc_col = self.column_vars.get((pair_key, 'log2FC'), tk.StringVar()).get()
            p_col = self.column_vars.get((pair_key, 'pvalue'), tk.StringVar()).get()
            
            result[pair_key] = {'log2FC': fc_col or None, 'pvalue': p_col or None}
        
        # Validate that at least one pair has both columns
        has_valid_pair = any(v.get('log2FC') and v.get('pvalue') for v in result.values())
        
        if not has_valid_pair:
            messagebox.showwarning("Validation", "At least one comparison must have both log2FC and p-value columns selected.")
            return
        
        self.result = result
        self.destroy()
    
    def _on_cancel(self):
        """Cancel and close."""
        self.result = None
        self.destroy()
    
    def wait_for_result(self):
        """Wait for dialog to close and return result."""
        self.wm_deiconify()
        self.update_idletasks()
        self.wait_window()
        return self.result


__all__ = ['PairwiseColumnMapper', 'PairwiseColumnDialog']
