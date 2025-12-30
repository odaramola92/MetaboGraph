"""Shared utilities and base classes for all tabs"""

from .utils import (
    TimerManager,
    SessionManager,
    ProgressUpdater,
    validate_file_path,
    validate_directory_path,
    LIPID_FEATURE_CANONICAL,
    normalize_col,
    is_lipid_feature,
    identify_sample_columns,
    identify_numeric_columns,
    build_sample_to_group_mapping,
    get_group_sample_counts,
    filter_groups_by_min_samples,
    split_pathways,
)

from .column_assignment import (
    show_column_assignment_dialog,
    ColumnDetector,
    ColumnCalculator,
    ColumnAssignmentDialog,
    COLUMN_TYPES,
    TAB_REQUIREMENTS,
)

__all__ = [
    # Utils
    'TimerManager',
    'SessionManager',
    'ProgressUpdater',
    'validate_file_path',
    'validate_directory_path',
    'LIPID_FEATURE_CANONICAL',
    'normalize_col',
    'is_lipid_feature',
    'identify_sample_columns',
    'identify_numeric_columns',
    'build_sample_to_group_mapping',
    'get_group_sample_counts',
    'filter_groups_by_min_samples',
    'split_pathways',
    # Column Assignment (New Unified System)
    'show_column_assignment_dialog',
    'ColumnDetector',
    'ColumnCalculator',
    'ColumnAssignmentDialog',
    'COLUMN_TYPES',
    'TAB_REQUIREMENTS',
]
