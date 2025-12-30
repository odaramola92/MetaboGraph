"""Shared utility functions for all tabs"""
from __future__ import annotations
import time
import logging
import tkinter as tk
import os
import json
import re
import sys
import shutil
from pathlib import Path
from typing import Dict, Any, Optional, Set


logger = logging.getLogger(__name__)


def find_database_paths():
    """
    Get prioritized list of paths to search for database files.
    
    Returns:
        list: List of directory paths to search, in priority order
    """
    if getattr(sys, 'frozen', False):
        # Running as PyInstaller executable
        exe_dir = os.path.dirname(sys.executable)
        return [
            os.path.join(exe_dir, "Databases"),  # Databases subfolder next to exe (PREFERRED)
            exe_dir,  # Same directory as executable
            os.path.join(os.getcwd(), "Databases"),  # Working directory Databases subfolder
            os.getcwd(),  # Working directory
        ]
    else:
        # Running as Python script
        return [
            os.path.join(os.getcwd(), "Databases"),  # Databases subfolder
            os.getcwd(),  # Working directory
            "Databases",  # Relative path
            ".",  # Current directory
        ]


def find_database(filename):
    """
    Find a database file in standard locations.
    
    Args:
        filename: Name of the database file (e.g., 'hmdb_database.feather')
    
    Returns:
        str: Full path to the database file, or None if not found
    """
    for search_path in find_database_paths():
        full_path = os.path.join(search_path, filename)
        if os.path.exists(full_path):
            return full_path
    return None


def get_user_config_dir(app_name: str = "MetaboliteAnnotationTool") -> Path:
    """Return a writable directory for storing user config/state files."""
    if sys.platform.startswith("win"):
        base_dir = Path(os.getenv("APPDATA", Path.home()))
    elif sys.platform == "darwin":
        base_dir = Path.home() / "Library" / "Application Support"
    else:
        base_dir = Path(os.getenv("XDG_CONFIG_HOME", Path.home() / ".config"))
    target = base_dir / app_name
    target.mkdir(parents=True, exist_ok=True)
    return target


def resolve_runtime_config_path(
    filename: str,
    template_dirs: Optional[Set[Path | str]] = None,
    app_name: str = "MetaboliteAnnotationTool",
) -> str:
    """Return path to writable config, cloning template if needed."""
    destination = get_user_config_dir(app_name) / filename
    if destination.exists():
        return str(destination)

    search_dirs = []
    if template_dirs:
        for directory in template_dirs:
            search_dirs.append(Path(directory))
    search_dirs.append(Path(__file__).resolve().parent)

    for directory in search_dirs:
        candidate = directory / filename
        if candidate.exists():
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(candidate, destination)
            return str(destination)

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.touch(exist_ok=True)
    return str(destination)


class TimerManager:
    """Manages timer, progress, and ETA calculations for long-running operations"""
    
    def __init__(self):
        """Initialize timer manager"""
        self.start_time = None
        self.current_step = 0
        self.total_steps = 100
        self.paused_time = 0
        self.pause_start = None
    
    def start(self, total_steps=100):
        """Start the timer"""
        self.start_time = time.time()
        self.current_step = 0
        self.total_steps = total_steps
        self.paused_time = 0
        self.pause_start = None
    
    def update_progress(self, current_step):
        """Update current step in progress"""
        self.current_step = current_step
    
    def pause(self):
        """Pause the timer"""
        if self.pause_start is None:
            self.pause_start = time.time()
    
    def resume(self):
        """Resume the timer"""
        if self.pause_start is not None:
            self.paused_time += time.time() - self.pause_start
            self.pause_start = None
    
    def get_elapsed_time(self):
        """Get elapsed time in seconds (excluding paused time)"""
        if self.start_time is None:
            return 0
        
        current_pause_offset = 0
        if self.pause_start is not None:
            current_pause_offset = time.time() - self.pause_start
        
        return time.time() - self.start_time - self.paused_time - current_pause_offset
    
    def get_elapsed_formatted(self):
        """Get elapsed time as MM:SS string"""
        elapsed = self.get_elapsed_time()
        minutes = int(elapsed // 60)
        seconds = int(elapsed % 60)
        return f"{minutes:02d}:{seconds:02d}"
    
    def get_eta_formatted(self):
        """Get ETA as MM:SS string"""
        if self.current_step <= 0 or self.total_steps <= 0:
            return "--:--"
        
        elapsed = self.get_elapsed_time()
        rate = elapsed / self.current_step
        remaining_steps = self.total_steps - self.current_step
        eta_seconds = remaining_steps * rate
        
        if eta_seconds < 0:
            eta_seconds = 0
        
        minutes = int(eta_seconds // 60)
        seconds = int(eta_seconds % 60)
        return f"{minutes:02d}:{seconds:02d}"
    
    def get_progress_percentage(self):
        """Get progress as percentage (0-100)"""
        if self.total_steps <= 0:
            return 0
        return int((self.current_step / self.total_steps) * 100)


class SessionManager:
    """Manages saving and loading GUI session state"""
    
    def __init__(self, app_name="MetaboliteAnnotation"):
        """Initialize session manager"""
        self.app_name = app_name
        self.session_dir = Path.home() / f".{app_name}_sessions"
        self.session_dir.mkdir(exist_ok=True)
        self.session_file = self.session_dir / "last_session.json"
    
    def save_session(self, session_data: Dict[str, Any]):
        """Save session data to file"""
        try:
            with open(self.session_file, 'w') as f:
                json.dump(session_data, f, indent=2, default=str)
            logger.info(f"Session saved to {self.session_file}")
            return True
        except Exception as e:
            logger.error(f"Failed to save session: {e}")
            return False
    
    def load_session(self) -> Optional[Dict[str, Any]]:
        """Load session data from file"""
        try:
            if self.session_file.exists():
                with open(self.session_file, 'r') as f:
                    data = json.load(f)
                logger.info(f"Session loaded from {self.session_file}")
                return data
        except Exception as e:
            logger.error(f"Failed to load session: {e}")
        return None
    
    def clear_session(self):
        """Clear saved session"""
        try:
            if self.session_file.exists():
                self.session_file.unlink()
            logger.info("Session cleared")
            return True
        except Exception as e:
            logger.error(f"Failed to clear session: {e}")
            return False


class ProgressUpdater:
    """Helper class to update progress bars and labels in a thread-safe way"""
    
    @staticmethod
    def update_progress(root_window, progress_var, percentage, message_label=None, message_text=""):
        """Update progress bar and message from background thread"""
        try:
            # Schedule GUI update in main thread
            root_window.after(0, lambda: ProgressUpdater._do_update(
                progress_var, percentage, message_label, message_text
            ))
        except Exception as e:
            logger.error(f"Failed to update progress: {e}")
    
    @staticmethod
    def _do_update(progress_var, percentage, message_label, message_text):
        """Actually perform the GUI update (must be called from main thread)"""
        try:
            progress_var.set(percentage)
            if message_label is not None and message_text:
                message_label.config(text=message_text)
        except tk.TclError:
            # GUI closed, ignore
            pass
        except Exception as e:
            logger.error(f"Error in progress update: {e}")


def format_file_size(size_bytes):
    """Format file size in human-readable format"""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} TB"


def validate_file_path(file_path, must_exist=True):
    """Validate file path"""
    try:
        path = Path(file_path)
        
        if must_exist and not path.exists():
            return False, f"File does not exist: {file_path}"
        
        if path.exists() and not path.is_file():
            return False, f"Path is not a file: {file_path}"
        
        return True, "Valid"
    except Exception as e:
        return False, str(e)


def validate_directory_path(dir_path):
    """Validate directory path"""
    try:
        path = Path(dir_path)
        
        if not path.exists():
            return False, f"Directory does not exist: {dir_path}"
        
        if not path.is_dir():
            return False, f"Path is not a directory: {dir_path}"
        
        # Check write permissions
        test_file = path / ".write_test"
        test_file.touch()
        test_file.unlink()
        
        return True, "Valid"
    except PermissionError:
        return False, f"No write permission: {dir_path}"
    except Exception as e:
        return False, str(e)


"""Shared normalization and canonical feature utilities for both statistics and visualization.

Centralizes logic to avoid drift between modules.
"""

# Canonical normalized lipid feature names (lowercase, non-alnum removed)
LIPID_FEATURE_CANONICAL: Set[str] = {
    'lipidid', 'class', 'lipidgroup', 'charge', 'calcmz', 'basert', 'subclass',
    'adduction', 'ionformula', 'molstructure', 'obsmz', 'obsrt', 'ppmdiff', 'polarity'
}

STAT_META_TOKENS = (
    'stat_', '_stat', 'statistic',
    'adj_p', 'adjp', 'padj',
    '_pvalue', '_p_value', '_p',
    'neg_log10'
)

def is_statistics_metadata_col(col_name: str) -> bool:
    """Return True if column name represents a statistical summary (not a raw sample)."""
    try:
        lowered = str(col_name).lower().strip()
    except Exception:
        lowered = str(col_name)
    if not lowered:
        return False
    if lowered.startswith('stat'):
        return True
    if lowered.startswith('adj') and ('p' in lowered or lowered.startswith('adjp')):
        return True
    if lowered.startswith('neg_log10'):
        return True
    
    # Exclude _pos and _neg suffixes before checking for _p token
    # These are polarity indicators, not p-value columns
    if lowered.endswith('_pos') or lowered.endswith('_neg'):
        return False
    if lowered.endswith('.pos') or lowered.endswith('.neg'):
        return False
    
    for token in STAT_META_TOKENS:
        if token in lowered:
            return True
    return False

def normalize_col(col_name: str) -> str:
    """Normalize a column name for robust matching.

    Lowercase and strip all non-alphanumeric characters so that variants like
    'CalcMz', 'Calc_MZ', 'calc m/z' all map to 'calcmz'.
    """
    try:
        s = str(col_name).lower()
        return re.sub(r'[^a-z0-9]', '', s)
    except Exception:
        return str(col_name)

def is_lipid_feature(col_name: str) -> bool:
    """Return True if the given column name corresponds to a canonical lipid feature field."""
    return normalize_col(col_name) in LIPID_FEATURE_CANONICAL


def identify_sample_columns(df) -> list[str]:
    """
    Identify sample columns in a DataFrame.
    
    Samples are typically numeric columns that are NOT lipid feature fields.
    Returns a list of column names that appear to be sample data.
    """
    sample_cols = []
    for col in df.columns:
        # Skip lipid feature or statistics metadata columns
        if is_lipid_feature(col) or is_statistics_metadata_col(col):
            continue
        
        # Check if column contains numeric data
        try:
            # Try to convert to numeric
            import pandas as pd
            numeric_col = pd.to_numeric(df[col], errors='coerce')
            # If most values are numeric, it's likely a sample column
            if numeric_col.notna().sum() / len(df) > 0.8:
                sample_cols.append(col)
        except Exception:
            pass
    
    return sample_cols


def identify_numeric_columns(df) -> list[str]:
    """
    Identify numeric columns in a DataFrame.
    
    Returns a list of column names that contain numeric data
    (excluding those that are clearly sample identifiers).
    """
    numeric_cols = []
    for col in df.columns:
        if is_statistics_metadata_col(col):
            continue
        try:
            import pandas as pd
            # Try to convert to numeric
            numeric_col = pd.to_numeric(df[col], errors='coerce')
            # If most values are numeric (>80%), it's numeric
            if numeric_col.notna().sum() / len(df) > 0.8:
                numeric_cols.append(col)
        except Exception:
            pass
    
    return numeric_cols


# ===== Group Management Utilities =====
# Shared utilities for sample grouping across Statistics and Visualization tabs

def build_sample_to_group_mapping(sample_group_vars: dict, group_definitions: dict) -> dict[str, str]:
    """
    Build a mapping from sample column names to group labels.
    
    Args:
        sample_group_vars: Dict of {sample_col_name: tk.StringVar with group assignment}
        group_definitions: Dict of {group_id: group_label}
    
    Returns:
        Dict of {sample_col_name: group_label}
    """
    mapping = {}
    if not sample_group_vars:
        return mapping
    
    # Build reverse lookup: group_id -> group_label
    current_labels = {gid: group_definitions.get(gid, gid) for gid in group_definitions}
    label_set = set(current_labels.values())
    
    for col_name, group_var in sample_group_vars.items():
        try:
            sel = group_var.get() if hasattr(group_var, 'get') else str(group_var)
            
            # If selection is an internal ID, convert to its label
            if sel in current_labels:
                mapping[col_name] = current_labels[sel]
            elif sel in label_set:
                mapping[col_name] = sel
            else:
                # Fallback to first defined label if something unexpected
                mapping[col_name] = next(iter(label_set)) if label_set else sel
        except Exception:
            pass
    
    return mapping


def get_group_sample_counts(sample_to_group: dict[str, str]) -> dict[str, int]:
    """
    Count how many samples are assigned to each group.
    
    Args:
        sample_to_group: Dict of {sample_col_name: group_label}
    
    Returns:
        Dict of {group_label: count}
    """
    from collections import Counter
    return dict(Counter(sample_to_group.values()))


def filter_groups_by_min_samples(
    sample_to_group: dict[str, str], 
    min_required: int
) -> tuple[dict[str, str], dict[str, int], dict[str, int]]:
    """
    Filter out groups that don't meet the minimum sample requirement.
    
    Args:
        sample_to_group: Dict of {sample_col_name: group_label}
        min_required: Minimum number of samples required per group
    
    Returns:
        Tuple of (filtered_mapping, kept_counts, removed_counts)
    """
    from collections import Counter
    
    # Count samples per group
    group_counts = Counter(sample_to_group.values())
    
    # Identify groups that meet the minimum
    valid_groups = {group for group, count in group_counts.items() if count >= min_required}
    
    # Filter mapping
    filtered_mapping = {
        sample: group 
        for sample, group in sample_to_group.items() 
        if group in valid_groups
    }
    
    # Track kept and removed counts
    kept_counts = {g: c for g, c in group_counts.items() if g in valid_groups}
    removed_counts = {g: c for g, c in group_counts.items() if g not in valid_groups}
    
    return filtered_mapping, kept_counts, removed_counts


def split_pathways(pathway_str: str) -> list:
    """
    Split pathway string by both ';' and '|' separators.
    
    Handles cases like:
    - "Path1; Path2; Path3"
    - "Path1|Path2|Path3"
    - "Path1; Path2|Path3"
    
    Args:
        pathway_str: String containing pathways separated by ; or |
        
    Returns:
        list: Individual pathway names (stripped and non-empty)
    """
    if not pathway_str or not isinstance(pathway_str, str):
        return []
    # Split by both ; and | characters
    pathways = re.split(r'[;|]', pathway_str)
    # Strip whitespace and filter empty strings
    return [p.strip() for p in pathways if p.strip()]


__all__ = [
    'TimerManager',
    'SessionManager',
    'ProgressUpdater',
    'validate_directory_path',
    'LIPID_FEATURE_CANONICAL',
    'normalize_col',
    'is_lipid_feature',
    'identify_sample_columns',
    'identify_numeric_columns',
    'build_sample_to_group_mapping',
    'get_group_sample_counts',
    'filter_groups_by_min_samples',
    'split_pathways'
]
