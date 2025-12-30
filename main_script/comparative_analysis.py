"""Comparative analysis utilities for MetaboGraph.

This module operates on multiple *annotated* datasets (same omics type across
studies or conditions) and produces summary tables that can be visualized in
Tkinter (heatmaps, barplots, etc.).

Each input dataset is expected to be a pandas DataFrame with at least:
- a feature identifier column (standardized to 'Name' or 'Feature_ID')
- a p-value column (standardized to 'pvalue')
- a log2 fold-change column (standardized to 'log2FC')
- one or more pathway columns OR a combined 'All_Pathways' column

The Multi-Omics/Comparative GUI tab is responsible for column-detection and
standardization before calling into this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class ComparativeDataset:
    """Lightweight wrapper for a single comparative dataset.

    Attributes
    ----------
    name:
        Human-readable label for the dataset (e.g. "Study1_Metabolomics").
    df:
        Annotated metabolite DataFrame.
    feature_col:
        Column name for feature identifier (typically 'Name').
    pvalue_col:
        Column name for p-values (typically 'pvalue').
    log2fc_col:
        Column name for log2 fold-change (typically 'log2FC').
    pathway_cols:
        List of pathway-related columns for this dataset. If empty, we fall
        back to any global 'All_Pathways' column if present.
    """

    name: str
    df: pd.DataFrame
    feature_col: str = "Name"
    pvalue_col: str = "pvalue"
    log2fc_col: str = "log2FC"
    pathway_cols: Optional[List[str]] = None

    def effective_pathway_columns(self) -> List[str]:
        """Return the list of pathway columns actually present in df.

        If explicit pathway_cols were provided, intersect them with df.columns.
        Otherwise, try to infer from df (including 'All_Pathways').
        """
        cols: List[str] = []
        if self.pathway_cols:
            cols = [c for c in self.pathway_cols if c in self.df.columns]
        if not cols:
            # Prefer explicit All_Pathways if present
            if "All_Pathways" in self.df.columns:
                cols = ["All_Pathways"]
            else:
                cols = [c for c in self.df.columns if "pathway" in c.lower()]
        return cols


def _standardize_pvalues(df: pd.DataFrame, p_col: str) -> pd.Series:
    """Return a numeric p-value Series with NaNs where conversion fails."""
    s = pd.to_numeric(df.get(p_col, np.nan), errors="coerce")
    return s


def build_feature_overlap_matrix(
    datasets: Sequence[ComparativeDataset],
    pvalue_threshold: float = 0.05,
) -> pd.DataFrame:
    """Build a binary feature-overlap matrix across datasets.

    Parameters
    ----------
    datasets:
        Iterable of ComparativeDataset objects.
    pvalue_threshold:
        P-value cutoff defining "significant" features in each dataset.

    Returns
    -------
    pd.DataFrame
        Index = feature identifiers (union across datasets).
        Columns = dataset names.
        Values = 1 if feature is significant (p < threshold) in that dataset,
        else 0. Features that never pass the threshold are dropped.
    """
    if not datasets:
        return pd.DataFrame()

    matrices = []
    for ds in datasets:
        df = ds.df
        if df is None or df.empty:
            continue
        if ds.feature_col not in df.columns:
            logger.warning("Dataset %s missing feature column %s", ds.name, ds.feature_col)
            continue
        pvals = _standardize_pvalues(df, ds.pvalue_col)
        signif_mask = pvals < float(pvalue_threshold)
        features = df.loc[signif_mask, ds.feature_col].astype(str)
        if features.empty:
            continue
        mat = pd.DataFrame(
            {ds.name: 1}, index=features
        )
        matrices.append(mat)

    if not matrices:
        return pd.DataFrame()

    overlap = pd.concat(matrices, axis=1).fillna(0).astype(int)
    # Drop features that are never significant in any dataset
    overlap = overlap.loc[overlap.sum(axis=1) > 0]
    # Ensure sorted index for stable display
    overlap = overlap.sort_index()
    logger.info("Comparative feature overlap matrix built with %d features and %d datasets", overlap.shape[0], overlap.shape[1])
    return overlap


def _extract_pathways_from_row(row: pd.Series, pathway_cols: Sequence[str]) -> List[str]:
    """Extract pathway names from a row given candidate columns.

    Splits on common separators (';', '|') and strips whitespace.
    """
    pathways: List[str] = []
    for col in pathway_cols:
        if col not in row.index:
            continue
        val = row.get(col)
        if val is None or (isinstance(val, float) and pd.isna(val)):
            continue
        s = str(val).strip()
        if not s or s.lower() == "nan":
            continue
        for part in str(s).replace("|", ";").split(";"):
            p = part.strip()
            if p:
                pathways.append(p)
    return pathways


def build_pathway_comparison_table(
    datasets: Sequence[ComparativeDataset],
    pvalue_threshold: float = 0.05,
) -> pd.DataFrame:
    """Summarize pathway involvement across datasets.

    For each dataset, we:
    - filter to significant features (p < threshold),
    - collect all pathways those features belong to,
    - count #significant features per pathway
      and compute the mean log2FC of those features.

    The returned table is suitable for a heatmap / dot plot view.

    Returns
    -------
    pd.DataFrame
        Multi-index columns:
            (dataset, 'count') and (dataset, 'mean_log2fc') for each dataset.
        Index:
            pathway names (union across datasets).
    """
    if not datasets:
        return pd.DataFrame()

    per_dataset_frames: List[pd.DataFrame] = []

    for ds in datasets:
        df = ds.df
        if df is None or df.empty:
            continue
        if ds.feature_col not in df.columns:
            logger.warning("Dataset %s missing feature column %s", ds.name, ds.feature_col)
            continue

        pvals = _standardize_pvalues(df, ds.pvalue_col)
        signif_mask = pvals < float(pvalue_threshold)
        if not signif_mask.any():
            logger.info("Dataset %s has no significant features at p < %g", ds.name, pvalue_threshold)
            continue

        work = df.loc[signif_mask].copy()
        # Ensure log2FC is numeric
        work[ds.log2fc_col] = pd.to_numeric(work.get(ds.log2fc_col, np.nan), errors="coerce")

        pathway_cols = ds.effective_pathway_columns()
        if not pathway_cols:
            logger.info("Dataset %s has no pathway columns; skipping pathway summary", ds.name)
            continue

        records: Dict[str, List[float]] = {}

        for _, row in work.iterrows():
            log2fc_val = row.get(ds.log2fc_col)
            try:
                log2fc_float = float(log2fc_val) if not pd.isna(log2fc_val) else np.nan
            except Exception:
                log2fc_float = np.nan

            pathways = _extract_pathways_from_row(row, pathway_cols)
            for p_name in pathways:
                key = str(p_name).strip()
                if not key:
                    continue
                if key not in records:
                    records[key] = []
                records[key].append(log2fc_float)

        if not records:
            logger.info("Dataset %s produced no pathway records after filtering", ds.name)
            continue

        rows = []
        for p_name, fc_values in records.items():
            arr = np.array([v for v in fc_values if not np.isnan(v)], dtype=float)
            count = int(len(fc_values))
            mean_fc = float(arr.mean()) if arr.size > 0 else np.nan
            rows.append({
                "Pathway": p_name,
                f"{ds.name}__count": count,
                f"{ds.name}__mean_log2fc": mean_fc,
            })

        if not rows:
            continue

        ds_table = pd.DataFrame(rows).set_index("Pathway")
        per_dataset_frames.append(ds_table)

    if not per_dataset_frames:
        return pd.DataFrame()

    # Outer-join all dataset-specific tables on pathway name
    combined = per_dataset_frames[0]
    for tbl in per_dataset_frames[1:]:
        combined = combined.join(tbl, how="outer")

    # Sort index for stable display
    combined = combined.sort_index()

    # Optionally, build a column MultiIndex: (dataset, metric)
    tuples: List[Tuple[str, str]] = []
    for col in combined.columns:
        if "__" in col:
            ds_name, metric = col.split("__", 1)
        else:
            ds_name, metric = "", col
        tuples.append((ds_name, metric))
    multi_cols = pd.MultiIndex.from_tuples(tuples, names=["dataset", "metric"])
    combined.columns = multi_cols

    logger.info(
        "Comparative pathway table built with %d pathways and %d datasets",
        combined.shape[0], len({d for d, _ in multi_cols}),
    )
    return combined


def summarize_comparative_results(
    datasets: Sequence[ComparativeDataset],
    pvalue_threshold: float = 0.05,
) -> Dict[str, pd.DataFrame]:
    """Compute all core comparative tables in one call.

    Returns a dict with keys:
    - 'features': binary feature overlap matrix
    - 'pathways': pathway comparison table
    Additional tables can be added later (upstream, diseases, etc.).
    """
    try:
        feature_overlap = build_feature_overlap_matrix(datasets, pvalue_threshold=pvalue_threshold)
        pathway_table = build_pathway_comparison_table(datasets, pvalue_threshold=pvalue_threshold)
        return {
            "features": feature_overlap,
            "pathways": pathway_table,
        }
    except Exception as e:
        logger.error("Error computing comparative results: %s", e, exc_info=True)
        return {}
