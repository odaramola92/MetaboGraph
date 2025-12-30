"""Machine learning enhancements for pathway-level analysis and networks.

This module adds publication-ready ML capabilities that complement the existing
ORA (Fisher) and IWPA scoring pipelines without changing their behavior.

Features
--------
- Graph-based diffusion (Random Walk with Restart) to propagate metabolite
  signals onto pathways for robust pathway activity scores.
- Network centrality metrics (PageRank, betweenness) for pathway prioritization.
- Compact API to augment existing pathway_stats with ML-derived scores.

Design
------
- Uses only core dependencies: numpy, pandas, networkx.
- Builds a bipartite graph between metabolites and pathways from your filtered
  `pathway_stats` and your original `df`.
- Seeds the graph with metabolite weights derived from log2FC/p-values and
  computes steady-state scores via RWR; returns per-pathway diffusion scores.

Usage
-----
>>> from metabolite_pathway_network import calculate_pathway_statistics
>>> from ml_pathway_models import augment_with_ml_scores
>>> pathway_stats = calculate_pathway_statistics(df)
>>> augmented = augment_with_ml_scores(df, pathway_stats)

"""
from __future__ import annotations

from typing import Dict, Iterable, Tuple
import numpy as np
import pandas as pd
import networkx as nx


def _build_bipartite_from_stats(
    df: pd.DataFrame,
    pathway_stats: Dict[str, dict],
) -> Tuple[nx.Graph, list[str], list[str]]:
    """Construct a metabolite–pathway bipartite graph from filtered stats.

    Nodes
    - metabolites: label `node_type='metabolite'`
    - pathways:    label `node_type='pathway'`

    Edges
    - metabolite -- pathway if metabolite listed under pathway in `pathway_stats`

    Returns (G, metabolite_nodes, pathway_nodes)
    """
    G = nx.Graph()

    # Collect metabolites per pathway from stats (supports alternate key names)
    def _extract_mets(stats_obj: dict) -> list[str]:
        if not isinstance(stats_obj, dict):
            return []
        for key in (
            'metabolites', 'Metabolites', 'metabolites_list', 'Metabolites_List',
            'Metabolites_in_Pathway', 'Metabolite_Names'
        ):
            if key in stats_obj and stats_obj[key] is not None:
                v = stats_obj[key]
                if isinstance(v, (list, tuple, set)):
                    return [str(x).strip() for x in v if str(x).strip()]
                if isinstance(v, str):
                    parts = [p.strip() for p in re.split(r"[|,;]", v) if p and p.strip()]
                    return parts
        return []

    import re as _re
    re = _re  # local alias for split

    # Add pathway nodes and cache metabolites
    pathway_nodes: list[str] = []
    pathway_to_mets: Dict[str, list[str]] = {}
    for p_name, stats in pathway_stats.items():
        if p_name.startswith('_'):
            # Skip metadata keys like _validation_report
            continue
        mets = _extract_mets(stats)
        if not mets:
            continue
        pathway_nodes.append(p_name)
        pathway_to_mets[p_name] = mets
        G.add_node(p_name, node_type='pathway')

    # Add metabolite nodes and edges
    metabolite_nodes: set[str] = set()
    for p_name, mets in pathway_to_mets.items():
        for m in mets:
            if not G.has_node(m):
                G.add_node(m, node_type='metabolite')
            metabolite_nodes.add(m)
            G.add_edge(m, p_name)

    return G, sorted(list(metabolite_nodes)), sorted(list(pathway_nodes))


def _seed_weights_from_df(
    df: pd.DataFrame,
    metabolites: Iterable[str],
    *,
    mode: str = 'signed_p',
    p_floor: float = 1e-16,
) -> np.ndarray:
    """Create a seed weight vector over metabolite nodes.

    Modes
    - 'signed_p': sign(log2FC) * -log10(p)
    - 'combined': log2FC * -log10(p)
    - 'log2fc':   raw log2FC
    """
    idx = {m: i for i, m in enumerate(metabolites)}
    s = np.zeros(len(idx), dtype=float)

    # Build lookup for speed
    df_idx = df.set_index('Name', drop=False)
    for m, i in idx.items():
        if m not in df_idx.index:
            s[i] = 0.0
            continue
        row = df_idx.loc[m]
        
        # Sanitize log2FC (handle NaN, None, inf)
        fc_raw = row.get('log2FC', 0.0)
        try:
            fc = float(fc_raw) if fc_raw is not None else 0.0
            if np.isnan(fc) or np.isinf(fc):
                fc = 0.0
        except (TypeError, ValueError):
            fc = 0.0
        
        # Sanitize pvalue (handle NaN, None, invalid range)
        pv_raw = row.get('pvalue', 1.0)
        try:
            pv = float(pv_raw) if pv_raw is not None else 1.0
            if np.isnan(pv) or pv <= 0 or pv > 1:
                pv = 1.0
        except (TypeError, ValueError):
            pv = 1.0
        pv = max(min(pv, 1.0), p_floor)

        if mode == 'signed_p':
            weight = np.sign(fc) * -np.log10(pv)
        elif mode == 'combined':
            weight = fc * -np.log10(pv)
        elif mode == 'log2fc':
            weight = fc
        else:
            raise ValueError(f"Unknown mode '{mode}'")
        
        # Final NaN/Inf check on computed weight
        if np.isnan(weight) or np.isinf(weight):
            s[i] = 0.0
        else:
            s[i] = weight

    # Normalize to unit L1 to prevent scale issues
    norm = np.sum(np.abs(s))
    if norm > 0:
        s = s / norm
    return s


def random_walk_with_restart(
    G: nx.Graph,
    metabolite_nodes: list[str],
    pathway_nodes: list[str],
    seed: np.ndarray,
    *,
    restart_prob: float = 0.3,
    max_iter: int = 200,
    tol: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray]:
    """Run Random Walk with Restart on the bipartite graph.

    Returns (p_metabolites, p_pathways) steady-state probabilities.

    Implementation details
    - Constructs adjacency over all nodes in order [metabolites..., pathways...]
    - Column-normalizes to a stochastic matrix W
    - Iterates: p_{t+1} = (1-r) * W @ p_t + r * p0, with p0 seeded on metabolite nodes
    """
    all_nodes = metabolite_nodes + pathway_nodes
    idx = {n: i for i, n in enumerate(all_nodes)}

    n = len(all_nodes)
    A = np.zeros((n, n), dtype=float)
    for u, v in G.edges():
        iu, iv = idx[u], idx[v]
        A[iu, iv] = 1.0
        A[iv, iu] = 1.0

    # Column-normalize (handle isolates defensively)
    col_sums = A.sum(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        W = np.divide(A, col_sums, where=col_sums != 0)

    # Initial distribution p0: mass only on metabolite nodes
    # Use absolute values to ensure non-negative probabilities in RWR
    p0 = np.zeros(n, dtype=float)
    p0[: len(metabolite_nodes)] = np.abs(seed)

    p = p0.copy()
    for _ in range(max_iter):
        p_next = (1.0 - restart_prob) * (W @ p) + restart_prob * p0
        if np.linalg.norm(p_next - p, 1) < tol:
            p = p_next
            break
        p = p_next

    return p[: len(metabolite_nodes)], p[len(metabolite_nodes) :]


def compute_diffusion_pathway_scores(
    df: pd.DataFrame,
    pathway_stats: Dict[str, dict],
    *,
    weight_mode: str = 'signed_p',
    restart_prob: float = 0.3,
    max_iter: int = 200,
    tol: float = 1e-8,
) -> pd.Series:
    """Compute RWR diffusion scores for pathway nodes.

    Returns a pd.Series indexed by pathway name with higher absolute values
    indicating stronger support from neighboring metabolite signals.
    """
    G, met_nodes, path_nodes = _build_bipartite_from_stats(df, pathway_stats)
    if not met_nodes or not path_nodes:
        return pd.Series(dtype=float)

    seed = _seed_weights_from_df(df, met_nodes, mode=weight_mode)
    _, p_path = random_walk_with_restart(
        G, met_nodes, path_nodes, seed,
        restart_prob=restart_prob, max_iter=max_iter, tol=tol
    )
    return pd.Series(p_path, index=path_nodes, name='ml_rwr_score').sort_values(ascending=False)


def compute_pathway_centrality(
    df: pd.DataFrame,
    pathway_stats: Dict[str, dict],
) -> pd.DataFrame:
    """Compute basic centrality metrics for pathway nodes.

    Returns DataFrame with columns: 'pagerank', 'betweenness'.
    """
    G, met_nodes, path_nodes = _build_bipartite_from_stats(df, pathway_stats)
    if not path_nodes:
        return pd.DataFrame(columns=['pagerank', 'betweenness'])

    pr = nx.pagerank(G, alpha=0.85)
    btw = nx.betweenness_centrality(G, normalized=True)
    data = {
        'pagerank': [pr.get(p, 0.0) for p in path_nodes],
        'betweenness': [btw.get(p, 0.0) for p in path_nodes],
    }
    return pd.DataFrame(data, index=path_nodes)


def augment_with_ml_scores(
    df: pd.DataFrame,
    pathway_stats: Dict[str, dict],
    *,
    weight_mode: str = 'signed_p',
    restart_prob: float = 0.3,
    max_iter: int = 200,
    tol: float = 1e-8,
) -> Dict[str, dict]:
    """Augment a pathway_stats dict with ML-derived scores (in-place copy).

    Adds the following keys per pathway:
    - 'ml_rwr_score': diffusion score from RWR (unitless, >= 0)
    - 'ml_pagerank':   PageRank centrality in the bipartite graph
    - 'ml_betweenness': betweenness centrality in the bipartite graph

    Returns a new dict (does not mutate input).
    """
    rwr = compute_diffusion_pathway_scores(
        df, pathway_stats, weight_mode=weight_mode, restart_prob=restart_prob, max_iter=max_iter, tol=tol
    )
    cent = compute_pathway_centrality(df, pathway_stats)

    enriched: Dict[str, dict] = {}
    for p, stats in pathway_stats.items():
        if p.startswith('_'):
            # pass through metadata entries unchanged
            enriched[p] = stats
            continue
        new_stats = dict(stats)
        new_stats['ml_rwr_score'] = float(rwr.get(p, 0.0))
        if p in cent.index:
            new_stats['ml_pagerank'] = float(cent.loc[p, 'pagerank'])
            new_stats['ml_betweenness'] = float(cent.loc[p, 'betweenness'])
        else:
            new_stats['ml_pagerank'] = 0.0
            new_stats['ml_betweenness'] = 0.0
        enriched[p] = new_stats

    return enriched


__all__ = [
    'compute_diffusion_pathway_scores',
    'compute_pathway_centrality',
    'augment_with_ml_scores',
    'compute_disease_diffusion_scores',
    'compute_disease_centrality',
]
 
# ===== Diseases ML utilities =====
from typing import Dict as _Dict, List as _List

def _build_metabolite_disease_graph(
    disease_to_metabolites: _Dict[str, _List[str]]
) -> Tuple[nx.Graph, list[str], list[str]]:
    """Construct bipartite graph between metabolites and diseases.

    Returns (G, metabolite_nodes, disease_nodes).
    """
    G = nx.Graph()
    disease_nodes: list[str] = []
    metabolite_nodes: set[str] = set()
    for disease, mets in disease_to_metabolites.items():
        if not disease or not mets:
            continue
        disease_nodes.append(disease)
        G.add_node(disease, node_type='disease')
        for m in mets:
            if not m:
                continue
            if not G.has_node(m):
                G.add_node(m, node_type='metabolite')
            metabolite_nodes.add(m)
            G.add_edge(m, disease)
    return G, sorted(list(metabolite_nodes)), sorted(disease_nodes)


def compute_disease_diffusion_scores(
    df: pd.DataFrame,
    disease_to_metabolites: _Dict[str, _List[str]],
    *,
    weight_mode: str = 'signed_p',
    restart_prob: float = 0.3,
    max_iter: int = 200,
    tol: float = 1e-8,
) -> pd.Series:
    """Compute RWR diffusion scores for diseases using metabolite evidence.

    Returns a pd.Series indexed by disease name.
    """
    G, met_nodes, dis_nodes = _build_metabolite_disease_graph(disease_to_metabolites)
    if not met_nodes or not dis_nodes:
        return pd.Series(dtype=float)

    seed = _seed_weights_from_df(df, met_nodes, mode=weight_mode)
    _, p_dis = random_walk_with_restart(
        G, met_nodes, dis_nodes, seed,
        restart_prob=restart_prob, max_iter=max_iter, tol=tol
    )
    return pd.Series(p_dis, index=dis_nodes, name='ml_disease_rwr_score').sort_values(ascending=False)


def compute_disease_centrality(
    disease_to_metabolites: _Dict[str, _List[str]],
) -> pd.DataFrame:
    """Compute PageRank and betweenness centrality for disease nodes."""
    G, met_nodes, dis_nodes = _build_metabolite_disease_graph(disease_to_metabolites)
    if not dis_nodes:
        return pd.DataFrame(columns=['pagerank', 'betweenness'])

    pr = nx.pagerank(G, alpha=0.85)
    btw = nx.betweenness_centrality(G, normalized=True)
    data = {
        'pagerank': [pr.get(d, 0.0) for d in dis_nodes],
        'betweenness': [btw.get(d, 0.0) for d in dis_nodes],
    }
    return pd.DataFrame(data, index=dis_nodes)

# ===== Upstream (Enzymes/Transporters) ML utilities =====
def _build_metabolite_upstream_graph(
    upstream_to_metabolites: _Dict[str, _List[str]]
) -> Tuple[nx.Graph, list[str], list[str]]:
    """Construct bipartite graph between metabolites and upstream regulators."""
    G = nx.Graph()
    upstream_nodes: list[str] = []
    metabolite_nodes: set[str] = set()
    for upstream, mets in upstream_to_metabolites.items():
        if not upstream or not mets:
            continue
        upstream_nodes.append(upstream)
        G.add_node(upstream, node_type='upstream')
        for m in mets:
            if not m:
                continue
            if not G.has_node(m):
                G.add_node(m, node_type='metabolite')
            metabolite_nodes.add(m)
            G.add_edge(m, upstream)
    return G, sorted(list(metabolite_nodes)), sorted(upstream_nodes)

def compute_upstream_diffusion_scores(
    df: pd.DataFrame,
    upstream_to_metabolites: _Dict[str, _List[str]],
    *,
    weight_mode: str = 'signed_p',
    restart_prob: float = 0.3,
    max_iter: int = 200,
    tol: float = 1e-8,
) -> pd.Series:
    """Compute RWR diffusion scores for upstream regulators."""
    G, met_nodes, up_nodes = _build_metabolite_upstream_graph(upstream_to_metabolites)
    if not met_nodes or not up_nodes:
        return pd.Series(dtype=float)
    seed = _seed_weights_from_df(df, met_nodes, mode=weight_mode)
    _, p_up = random_walk_with_restart(
        G, met_nodes, up_nodes, seed,
        restart_prob=restart_prob, max_iter=max_iter, tol=tol
    )
    return pd.Series(p_up, index=up_nodes, name='ml_upstream_rwr_score').sort_values(ascending=False)

def compute_upstream_centrality(
    upstream_to_metabolites: _Dict[str, _List[str]],
) -> pd.DataFrame:
    """Compute PageRank and betweenness centrality for upstream nodes."""
    G, met_nodes, up_nodes = _build_metabolite_upstream_graph(upstream_to_metabolites)
    if not up_nodes:
        return pd.DataFrame(columns=['pagerank', 'betweenness'])
    pr = nx.pagerank(G, alpha=0.85)
    btw = nx.betweenness_centrality(G, normalized=True)
    data = {
        'pagerank': [pr.get(u, 0.0) for u in up_nodes],
        'betweenness': [btw.get(u, 0.0) for u in up_nodes],
    }
    return pd.DataFrame(data, index=up_nodes)

