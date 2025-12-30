"""
Interactive Network Viewer for Compound Pathway Analysis
============================================================
Provides an embedded interactive network visualization similar to MetaboAnalyst's EnrichNet.

Uses Dash + Cytoscape.js for web-based interactivity that can be:
1. Opened in a browser window (best experience)
2. Embedded in tkinter via webview (requires pywebview)

Features:
- Interactive zoom/pan
- Node selection with details panel
- Multiple layout algorithms (force-directed, hierarchical, circular)
- Node size based on significance/hits
- Color coding by p-value/FDR
- Edge weights showing pathway connections
- Export to PNG/SVG

Requirements:
    pip install dash dash-cytoscape pandas networkx

Usage:
    from interactive_network_viewer import PathwayNetworkViewer
    
    viewer = PathwayNetworkViewer(pathway_data_df)
    viewer.run()  # Opens in browser
    # OR
    viewer.run_embedded(parent_window)  # Embeds in tkinter
"""

import pandas as pd
import networkx as nx
import threading
import webbrowser
import socket
import time
from typing import Optional, Dict, List, Any, Callable

# Dash and Cytoscape imports
try:
    import dash
    from dash import html, dcc, Input, Output, State, callback_context
    import dash_cytoscape as cyto
    DASH_AVAILABLE = True
except ImportError:
    DASH_AVAILABLE = False
    print("⚠️ dash and dash-cytoscape not installed. Install with: pip install dash dash-cytoscape")

# Optional: pywebview for tkinter embedding
try:
    import webview
    WEBVIEW_AVAILABLE = True
except ImportError:
    WEBVIEW_AVAILABLE = False


def find_free_port(start=8050, max_attempts=100):
    """Find an available port starting from 'start'."""
    for port in range(start, start + max_attempts):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('127.0.0.1', port))
                return port
        except OSError:
            continue
    raise RuntimeError(f"Could not find free port in range {start}-{start+max_attempts}")


class PathwayNetworkViewer:
    """
    Interactive pathway network viewer using Dash + Cytoscape.js
    
    Provides two visualization modes:
    1. Simple Mode: Pathway-only network colored by p-value
    2. Cytoscape Mode: Full network with compounds, activation/inhibition coloring
    
    Features:
    - Force-directed and hierarchical layouts
    - Node sizing by hits/significance
    - Color coding by p-value (simple) or activation/inhibition (cytoscape)
    - Click to select and view details
    - Zoom, pan, and export
    """
    
    # Color scales for p-value visualization (Simple Mode)
    PVALUE_COLORS = {
        'very_sig': '#FF4500',      # Red-orange (p < 0.001)
        'highly_sig': '#FF8C00',    # Dark orange (p < 0.01)
        'significant': '#FFD700',   # Gold (p < 0.05)
        'marginal': '#FFFFE0',      # Light yellow (p < 0.1)
        'not_sig': '#D3D3D3'        # Light gray (p >= 0.1)
    }
    
    # Cytoscape-style colors for activation/inhibition
    CYTOSCAPE_COLORS = {
        'activated': '#FF0000',      # Red for activated pathways
        'inhibited': '#3498DB',      # Blue for inhibited pathways  
        'neutral': '#BDC3C7',        # Gray for neutral
        'upregulated': '#FF8C00',    # Orange for upregulated compounds
        'downregulated': '#228B22',  # Green for downregulated compounds
        'no_change': '#CCCCCC',      # Gray for no change
        'upstream': '#8e44ad',       # Purple for enzymes/transporters
        'disease': '#8e44ad',        # Purple for diseases
        'association': '#888888',    # Gray for disease/upstream associations
    }
    
    # Available layout algorithms
    LAYOUTS = {
        'cose': 'Force-Directed (CoSE)',
        'cola': 'Force-Directed (Cola)',
        'dagre': 'Hierarchical (Dagre)',
        'circle': 'Circular',
        'concentric': 'Concentric',
        'grid': 'Grid',
        'breadthfirst': 'Breadth-First',
        'random': 'Random'
    }
    
    def __init__(
        self,
        pathway_data: pd.DataFrame = None,
        metabolite_data: pd.DataFrame = None,
        upstream_data: pd.DataFrame = None,
        disease_data: pd.DataFrame = None,
        pathway_stats: Dict = None,
        title: str = "Pathway Enrichment Network",
        mode: str = 'simple',        # 'simple' or 'cytoscape'
        show_edges: bool = True,
        min_edge_weight: float = 0.1,
        node_size_by: str = 'hits',  # 'hits', 'pvalue', or 'fdr'
        color_by: str = 'pvalue',    # 'pvalue' or 'fdr' (simple mode only)
        focused_mode: bool = True   # True = only pathway-connected compounds; False = all compounds
    ):
        """
        Initialize the network viewer.
        
        Args:
            pathway_data: DataFrame with columns like ['Pathway', 'Hits', 'P-Value', 'FDR', 'Compounds']
            metabolite_data: DataFrame with compound info ['Name', 'log2FC', 'pvalue'] (for cytoscape mode)
            upstream_data: DataFrame with upstream regulator info (for cytoscape mode)
            disease_data: DataFrame with disease info (for cytoscape mode)
            pathway_stats: Optional dict with pathway statistics including z_score, status, etc.
            title: Title for the network visualization
            mode: 'simple' (pathway-only, p-value colors) or 'cytoscape' (full network, activation colors)
            show_edges: Whether to show edges between pathways
            min_edge_weight: Minimum edge weight to display (0-1)
            node_size_by: What to base node size on
            color_by: What to base node color on (simple mode)
            focused_mode: If True, only show compounds connected to pathways; if False, show all compounds
        """
        self.pathway_data = pathway_data
        self.metabolite_data = metabolite_data
        self.upstream_data = upstream_data
        self.disease_data = disease_data
        self.pathway_stats = pathway_stats or {}
        self.title = title
        self.mode = mode
        self.show_edges = show_edges
        self.min_edge_weight = min_edge_weight
        self.node_size_by = node_size_by
        self.color_by = color_by
        self.focused_mode = focused_mode  # Focused mode from parameter
        
        self.app = None
        self.port = None
        self.server_thread = None
        self._shutdown_event = threading.Event()
        
        # Build network graph
        self.G = nx.Graph()
        self.elements = []
        
        # Debug logging
        print(f"\n🔍 PathwayNetworkViewer initialization:")
        print(f"   Mode: {self.mode}")
        print(f"   Pathway data: {len(pathway_data) if pathway_data is not None else 0} rows")
        print(f"   Compound data: {len(metabolite_data) if metabolite_data is not None else 0} rows")
        if metabolite_data is not None and not metabolite_data.empty:
            print(f"   Compound columns: {metabolite_data.columns.tolist()}")
        
        if pathway_data is not None:
            self._build_network()
    
    def _get_node_color(self, pvalue: float) -> str:
        """Get node color based on p-value."""
        if pvalue < 0.001:
            return self.PVALUE_COLORS['very_sig']
        elif pvalue < 0.01:
            return self.PVALUE_COLORS['highly_sig']
        elif pvalue < 0.05:
            return self.PVALUE_COLORS['significant']
        elif pvalue < 0.1:
            return self.PVALUE_COLORS['marginal']
        else:
            return self.PVALUE_COLORS['not_sig']
    
    def _get_node_size(self, hits: int, pvalue: float) -> int:
        """Calculate node size based on hits or significance."""
        if self.node_size_by == 'hits':
            # Scale by number of hits (compounds in pathway)
            return max(30, min(100, 20 + hits * 6))
        elif self.node_size_by == 'pvalue':
            # Scale by -log10(pvalue)
            import math
            log_p = -math.log10(max(pvalue, 1e-50))
            return max(30, min(100, 15 + log_p * 6))
        else:
            return 50  # Default size
    
    def _build_network(self):
        """Build network elements from pathway data."""
        if self.pathway_data is None or self.pathway_data.empty:
            return
        
        if self.mode == 'full_network':
            self._build_full_network()
        else:
            self._build_simple_network()
    
    def _build_simple_network(self):
        """Build simple pathway-only network with p-value coloring."""
        df = self.pathway_data.copy()
        
        # Detect column names (flexible naming)
        name_col = self._find_column(df, ['Pathway', 'Name', 'pathway', 'name', 'Pathway_Name'])
        hits_col = self._find_column(df, ['Hits', 'hits', '# Metabolites', '# Compounds', 'Count', 'count'])
        pval_col = self._find_column(df, ['P-Value', 'pvalue', 'p_value', 'P_Value', 'pval', 'Pvalue'])
        fdr_col = self._find_column(df, ['FDR', 'fdr', 'q_value', 'adj_pvalue', 'padj', 'BH_FDR'])
        zscore_col = self._find_column(df, ['Z-Score', 'z_score', 'zscore', 'Z_Score'])
        status_col = self._find_column(df, ['Status', 'status'])
        met_col = self._find_column(df, ['Metabolites', 'metabolites', 'Compounds', 'Members', 'Genes'])
        
        if name_col is None:
            print("⚠️ Could not find pathway name column")
            return
        
        # Build nodes
        nodes = []
        for idx, row in df.iterrows():
            pathway_name = str(row[name_col])
            hits = int(row[hits_col]) if hits_col and pd.notna(row.get(hits_col)) else 1
            pvalue = float(row[pval_col]) if pval_col and pd.notna(row.get(pval_col)) else 1.0
            fdr = float(row[fdr_col]) if fdr_col and pd.notna(row.get(fdr_col)) else 1.0
            z_score = float(row[zscore_col]) if zscore_col and pd.notna(row.get(zscore_col)) else 0.0
            status = str(row[status_col]) if status_col and pd.notna(row.get(status_col)) else 'N/A'
            metabolites = str(row[met_col]) if met_col and pd.notna(row.get(met_col)) else ""
            
            # Determine color
            color_val = pvalue if self.color_by == 'pvalue' else fdr
            node_color = self._get_node_color(color_val)
            
            # Determine size
            node_size = self._get_node_size(hits, pvalue)
            
            # Create node
            node = {
                'data': {
                    'id': pathway_name,
                    'label': pathway_name,
                    'full_name': pathway_name,
                    'hits': hits,
                    'pvalue': f"{pvalue:.2e}",
                    'fdr': f"{fdr:.2e}",
                    'z_score': f"{z_score:.2f}",
                    'status': status,
                    'metabolites': metabolites,
                    'size': node_size,
                    'color': node_color,
                    'node_type': 'pathway',
                    'shape': 'ellipse',  # Circles for simple mode
                }
            }
            nodes.append(node)
            
            # Add to NetworkX graph for edge calculation
            self.G.add_node(pathway_name, hits=hits, pvalue=pvalue, fdr=fdr, z_score=z_score, status=status, metabolites=metabolites, node_type='pathway')
        
        # Build edges based on shared compounds
        edges = []
        if self.show_edges and met_col:
            pathway_metabolites = {}
            for idx, row in df.iterrows():
                pathway_name = str(row[name_col])
                mets = str(row.get(met_col, ""))
                if mets and mets != 'nan':
                    # Parse metabolite list (comma, semicolon, or pipe separated)
                    met_list = set()
                    for sep in [';', ',', '|']:
                        if sep in mets:
                            met_list = set(m.strip() for m in mets.split(sep) if m.strip())
                            break
                    if not met_list:
                        met_list = {mets.strip()}
                    pathway_metabolites[pathway_name] = met_list
            
            # Calculate Jaccard similarity for edges
            pathways = list(pathway_metabolites.keys())
            for i, p1 in enumerate(pathways):
                for j, p2 in enumerate(pathways[i+1:], i+1):
                    mets1 = pathway_metabolites[p1]
                    mets2 = pathway_metabolites[p2]
                    
                    intersection = len(mets1 & mets2)
                    union = len(mets1 | mets2)
                    
                    if union > 0 and intersection > 0:
                        jaccard = intersection / union
                        if jaccard >= self.min_edge_weight:
                            edge = {
                                'data': {
                                    'source': p1,
                                    'target': p2,
                                    'weight': jaccard,
                                    'shared': intersection,
                                }
                            }
                            edges.append(edge)
                            self.G.add_edge(p1, p2, weight=jaccard, shared=intersection)
        
        self.elements = nodes + edges
        print(f"✅ Built simple network with {len(nodes)} pathways and {len(edges)} connections")
    
    def _build_full_network(self):
        """Build full pathway-metabolite network with activation/inhibition coloring.
        
        Supports two modes:
        - Focused Mode (focused_mode=True): Only show metabolites connected to pathways.
          Upstream/diseases only connect to pathway-linked metabolites.
        - General Mode (focused_mode=False): Show all metabolites including secondary
          connections from upstream/disease that may not have pathway links.
        """
        print(f"\n🔶 Building FULL NETWORK (pathways + metabolites + upstream + diseases)...")
        print(f"   Focused Mode: {self.focused_mode}")
        nodes = []
        edges = []
        
        df_pathways = self.pathway_data.copy()
        df_metabolites = self.metabolite_data.copy() if self.metabolite_data is not None else None
        
        print(f"   Pathway DataFrame: {len(df_pathways)} rows")
        print(f"   Metabolite DataFrame: {len(df_metabolites) if df_metabolites is not None else 0} rows")
        print(f"   df_metabolites is None: {df_metabolites is None}")
        print(f"   df_metabolites is empty: {df_metabolites.empty if df_metabolites is not None else 'N/A'}")
        
        # Find pathway columns
        pathway_name_col = self._find_column(df_pathways, ['Pathway', 'Name', 'pathway'])
        if pathway_name_col is None:
            print("⚠️ Could not find pathway name column")
            return
        
        # First pass: collect all metabolites linked to pathways
        pathway_metabolites = set()  # Metabolites that have pathway connections
        met_col = self._find_column(df_pathways, ['Metabolites', 'metabolites', 'Compounds'])
        if met_col:
            for idx, row in df_pathways.iterrows():
                mets_str = str(row.get(met_col, ""))
                if mets_str and mets_str != 'nan':
                    for sep in [';', ',', '|']:
                        if sep in mets_str:
                            pathway_metabolites.update([m.strip() for m in mets_str.split(sep) if m.strip()])
                            break
                    else:
                        pathway_metabolites.add(mets_str.strip())
        
        print(f"   📊 Found {len(pathway_metabolites)} metabolites linked to pathways")
        
        # Build metabolite nodes if available
        metabolite_map = {}
        if df_metabolites is not None and not df_metabolites.empty:
            met_name_col = self._find_column(df_metabolites, ['Name', 'name', 'Metabolite', 'metabolite'])
            log2fc_col = self._find_column(df_metabolites, ['log2FC', 'log2_FC', 'Log2FC', 'log2FoldChange'])
            pval_col = self._find_column(df_metabolites, ['pvalue', 'p_value', 'P-Value'])
            
            print(f"   📊 Metabolite data columns: {df_metabolites.columns.tolist()}")
            print(f"   📊 Found columns - Name: {met_name_col}, log2FC: {log2fc_col}, pvalue: {pval_col}")
            
            if met_name_col and log2fc_col:
                for idx, row in df_metabolites.iterrows():
                    met_name = str(row[met_name_col])
                    log2fc = float(row[log2fc_col]) if pd.notna(row.get(log2fc_col)) else 0.0
                    pval = float(row[pval_col]) if pval_col and pd.notna(row.get(pval_col)) else 1.0
                    
                    # Determine regulation
                    if log2fc > 0:
                        regulation = 'upregulated'
                        color = self.CYTOSCAPE_COLORS['upregulated']
                    elif log2fc < 0:
                        regulation = 'downregulated'
                        color = self.CYTOSCAPE_COLORS['downregulated']
                    else:
                        regulation = 'no_change'
                        color = self.CYTOSCAPE_COLORS['no_change']
                    
                    # Store in map for reference
                    metabolite_map[met_name] = {
                        'log2fc': log2fc,
                        'pval': pval,
                        'regulation': regulation,
                        'color': color,
                        'has_pathway': met_name in pathway_metabolites,
                    }
                    
                    # In focused mode, only add metabolite nodes if they're connected to a pathway
                    # OR if they are significant (not gray/no_change)
                    # In general mode, add all metabolites
                    should_add_node = True
                    if self.focused_mode:
                        # Only add if connected to pathway
                        should_add_node = met_name in pathway_metabolites
                    
                    if should_add_node:
                        node = {
                            'data': {
                                'id': f'met_{met_name}',
                                'label': met_name,
                                'full_name': met_name,
                                'log2fc': f"{log2fc:.2f}",
                                'pvalue': f"{pval:.2e}",
                                'regulation': regulation,
                                'size': 30,
                                'color': color,
                                'node_type': 'metabolite',
                                'shape': 'ellipse',
                                'has_pathway': met_name in pathway_metabolites,
                            }
                        }
                        nodes.append(node)
                        self.G.add_node(f'met_{met_name}', node_type='metabolite', log2fc=log2fc, regulation=regulation)
                
                print(f"   ✓ Created {len([n for n in nodes if n['data']['node_type']=='metabolite'])} metabolite nodes")
            else:
                print(f"   ⚠️ Missing required columns for metabolites: name={met_name_col}, log2fc={log2fc_col}")
        
        # Build set of metabolite node IDs we've created
        created_met_nodes = set(n['data']['id'].replace('met_', '') for n in nodes if n['data']['node_type'] == 'metabolite')
        
        # Build pathway nodes with activation/inhibition coloring
        zscore_col = self._find_column(df_pathways, ['Z-Score', 'z_score', 'zscore', 'Z_Score'])
        status_col = self._find_column(df_pathways, ['Status', 'status'])
        hits_col = self._find_column(df_pathways, ['Hits', 'hits', '# Metabolites', '# Compounds', 'Count', 'count'])
        pval_col = self._find_column(df_pathways, ['P-Value', 'pvalue', 'p_value', 'P_Value', 'pval', 'Pvalue'])
        
        for idx, row in df_pathways.iterrows():
            pathway_name = str(row[pathway_name_col])
            
            # Get z-score and status from DataFrame columns
            z_score = float(row[zscore_col]) if zscore_col and pd.notna(row.get(zscore_col)) else 0.0
            status = str(row[status_col]) if status_col and pd.notna(row.get(status_col)) else 'No Change'
            
            # Determine color based on status (not z-score alone)
            # Map status text to colors - status takes priority over z-score
            status_lower = str(status).lower()
            if 'activat' in status_lower:
                color = self.CYTOSCAPE_COLORS['activated']
            elif 'inhibit' in status_lower:
                color = self.CYTOSCAPE_COLORS['inhibited']
            else:
                # "No Change" or any other status gets neutral gray
                color = self.CYTOSCAPE_COLORS['neutral']
            
            hits = int(row[hits_col]) if hits_col and pd.notna(row.get(hits_col)) else 1
            pvalue = float(row[pval_col]) if pval_col and pd.notna(row.get(pval_col)) else 1.0
            
            node = {
                'data': {
                    'id': f'path_{pathway_name}',
                    'label': pathway_name,
                    'full_name': pathway_name,
                    'hits': hits,
                    'status': status,
                    'z_score': f"{z_score:.2f}" if z_score is not None else 'N/A',
                    'pvalue': f"{pvalue:.2e}",
                    'size': self._get_node_size(hits, pvalue),
                    'color': color,
                    'node_type': 'pathway',
                    'shape': 'diamond',
                }
            }
            nodes.append(node)
            self.G.add_node(f'path_{pathway_name}', node_type='pathway', z_score=z_score, status=status)
            
            # Connect pathway to metabolites
            if met_col:
                mets_str = str(row.get(met_col, ""))
                if mets_str and mets_str != 'nan':
                    # Parse metabolite list
                    met_list = []
                    for sep in [';', ',', '|']:
                        if sep in mets_str:
                            met_list = [m.strip() for m in mets_str.split(sep) if m.strip()]
                            break
                    if not met_list:
                        met_list = [mets_str.strip()]
                    
                    # Create edges from pathway to metabolites
                    for met in met_list:
                        # Only connect if metabolite node exists and is SIGNIFICANT (not gray)
                        if met in created_met_nodes and met in metabolite_map:
                            met_info = metabolite_map[met]
                            # Skip gray (no_change) metabolites - don't connect them
                            if met_info['regulation'] == 'no_change':
                                continue
                                
                            # Determine edge color based on regulation - use RED for upregulated
                            if met_info['regulation'] == 'upregulated':
                                edge_color = self.CYTOSCAPE_COLORS['activated']  # Red
                            elif met_info['regulation'] == 'downregulated':
                                edge_color = self.CYTOSCAPE_COLORS['inhibited']  # Blue
                            else:
                                edge_color = '#808080'  # Gray
                            
                            edge = {
                                'data': {
                                    'source': f'path_{pathway_name}',
                                    'target': f'met_{met}',
                                    'weight': 1.0,
                                    'color': edge_color,
                                }
                            }
                            edges.append(edge)
                            self.G.add_edge(f'path_{pathway_name}', f'met_{met}')
        
        # Build upstream regulator nodes (enzymes/transporters) if provided
        print(f"\n   🔍 Checking upstream_data:")
        print(f"      upstream_data is None: {self.upstream_data is None}")
        if self.upstream_data is not None:
            print(f"      upstream_data type: {type(self.upstream_data)}")
            print(f"      upstream_data empty: {self.upstream_data.empty if hasattr(self.upstream_data, 'empty') else 'N/A'}")
            if hasattr(self.upstream_data, 'columns'):
                print(f"      upstream_data columns: {self.upstream_data.columns.tolist()}")
                print(f"      upstream_data shape: {self.upstream_data.shape}")
        
        if self.upstream_data is not None and not self.upstream_data.empty:
            name_col = self._find_column(self.upstream_data, ['name', 'Name', 'gene', 'Gene'])
            type_col = self._find_column(self.upstream_data, ['type', 'Type'])
            met_col_up = self._find_column(self.upstream_data, ['metabolites', 'Metabolites'])
            
            if name_col:
                for idx, row in self.upstream_data.iterrows():
                    upstream_name = str(row[name_col])
                    upstream_type = str(row.get(type_col, 'enzyme')) if type_col else 'enzyme'
                    mets = row.get(met_col_up, []) if met_col_up else []
                    
                    if isinstance(mets, str):
                        mets = [m.strip() for m in mets.split(';') if m.strip()]
                    
                    # Filter metabolites based on mode
                    # In focused mode: only connect to metabolites that have pathway connections AND are significant
                    # In general mode: connect to all metabolites that exist as nodes AND are significant
                    valid_mets = []
                    for met in mets:
                        if met not in metabolite_map:
                            continue
                        met_info = metabolite_map[met]
                        # Skip gray (no_change) metabolites
                        if met_info['regulation'] == 'no_change':
                            continue
                        if self.focused_mode:
                            # Only include if metabolite has pathway connection
                            if met_info.get('has_pathway', False):
                                valid_mets.append(met)
                        else:
                            # General mode: include all significant metabolites
                            valid_mets.append(met)
                    
                    # Only add upstream node if it has valid metabolite connections
                    if not valid_mets:
                        continue
                    
                    node = {
                        'data': {
                            'id': f'upstream_{upstream_name}',
                            'label': upstream_name,
                            'full_name': upstream_name,
                            'type': upstream_type,
                            'n_metabolites': len(valid_mets),
                            'size': 35,
                            'color': self.CYTOSCAPE_COLORS['upstream'],
                            'node_type': 'upstream',
                            'shape': 'triangle',
                        }
                    }
                    nodes.append(node)
                    self.G.add_node(f'upstream_{upstream_name}', node_type='upstream')
                    
                    # Connect upstream to metabolites with gray (association) edges
                    for met in valid_mets:
                        if met in created_met_nodes:
                            edge = {
                                'data': {
                                    'source': f'upstream_{upstream_name}',
                                    'target': f'met_{met}',
                                    'weight': 1.0,
                                    'color': self.CYTOSCAPE_COLORS['association'],  # Gray for associations
                                }
                            }
                            edges.append(edge)
                            self.G.add_edge(f'upstream_{upstream_name}', f'met_{met}')
        
        # Build disease nodes if provided
        print(f"\n   🔍 Checking disease_data:")
        print(f"      disease_data is None: {self.disease_data is None}")
        if self.disease_data is not None:
            print(f"      disease_data type: {type(self.disease_data)}")
            print(f"      disease_data empty: {self.disease_data.empty if hasattr(self.disease_data, 'empty') else 'N/A'}")
            if hasattr(self.disease_data, 'columns'):
                print(f"      disease_data columns: {self.disease_data.columns.tolist()}")
                print(f"      disease_data shape: {self.disease_data.shape}")
        
        if self.disease_data is not None and not self.disease_data.empty:
            name_col = self._find_column(self.disease_data, ['name', 'Name', 'disease', 'Disease'])
            met_col_dis = self._find_column(self.disease_data, ['metabolites', 'Metabolites'])
            
            if name_col:
                for idx, row in self.disease_data.iterrows():
                    disease_name = str(row[name_col])
                    mets = row.get(met_col_dis, []) if met_col_dis else []
                    
                    if isinstance(mets, str):
                        mets = [m.strip() for m in mets.split(';') if m.strip()]
                    
                    # Filter metabolites based on mode (same logic as upstream)
                    valid_mets = []
                    for met in mets:
                        if met not in metabolite_map:
                            continue
                        met_info = metabolite_map[met]
                        # Skip gray (no_change) metabolites
                        if met_info['regulation'] == 'no_change':
                            continue
                        if self.focused_mode:
                            # Only include if metabolite has pathway connection
                            if met_info.get('has_pathway', False):
                                valid_mets.append(met)
                        else:
                            # General mode: include all significant metabolites
                            valid_mets.append(met)
                    
                    # Only add disease node if it has valid metabolite connections
                    if not valid_mets:
                        continue
                    
                    node = {
                        'data': {
                            'id': f'disease_{disease_name}',
                            'label': disease_name,
                            'full_name': disease_name,
                            'n_metabolites': len(valid_mets),
                            'size': 38,
                            'color': self.CYTOSCAPE_COLORS['disease'],
                            'node_type': 'disease',
                            'shape': 'hexagon',
                        }
                    }
                    nodes.append(node)
                    self.G.add_node(f'disease_{disease_name}', node_type='disease')
                    
                    # Connect disease to metabolites with gray (association) edges
                    for met in valid_mets:
                        if met in created_met_nodes:
                            edge = {
                                'data': {
                                    'source': f'disease_{disease_name}',
                                    'target': f'met_{met}',
                                    'weight': 1.0,
                                    'color': self.CYTOSCAPE_COLORS['association'],  # Gray for associations
                                }
                            }
                            edges.append(edge)
                            self.G.add_edge(f'disease_{disease_name}', f'met_{met}')
        
        # Remove orphan gray metabolite nodes (nodes with no connections)
        connected_mets = set()
        for edge in edges:
            src = edge['data']['source']
            tgt = edge['data']['target']
            if src.startswith('met_'):
                connected_mets.add(src)
            if tgt.startswith('met_'):
                connected_mets.add(tgt)
        
        # Filter out unconnected metabolite nodes
        nodes = [n for n in nodes if n['data']['node_type'] != 'metabolite' or n['data']['id'] in connected_mets]
        
        self.elements = nodes + edges
        node_counts = {
            'pathway': len([n for n in nodes if n['data']['node_type']=='pathway']),
            'metabolite': len([n for n in nodes if n['data']['node_type']=='metabolite']),
            'upstream': len([n for n in nodes if n['data']['node_type']=='upstream']),
            'disease': len([n for n in nodes if n['data']['node_type']=='disease']),
        }
        print(f"✅ Built full network with {node_counts['pathway']} pathways, "
              f"{node_counts['metabolite']} metabolites, {node_counts['upstream']} upstream, "
              f"{node_counts['disease']} diseases, and {len(edges)} connections")
        print(f"   Mode: {'Focused' if self.focused_mode else 'General'}")
    
    def _find_column(self, df, candidates):
        """Find first matching column from list of candidates."""
        for col in candidates:
            if col in df.columns:
                return col
        return None
    
    def _get_legend_html(self, mode=None):
        """Get legend HTML based on current mode."""
        mode = mode if mode is not None else self.mode
        if mode == 'full_network':
            # Cytoscape mode legend - activation/inhibition
            legend_items = [
                html.P("Pathways:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                html.Div([
                    html.Span("◆", style={'color': self.CYTOSCAPE_COLORS['activated'], 'fontSize': '20px'}),
                    html.Span(" Activated", style={'marginLeft': '5px'})
                ]),
                html.Div([
                    html.Span("◆", style={'color': self.CYTOSCAPE_COLORS['inhibited'], 'fontSize': '20px'}),
                    html.Span(" Inhibited", style={'marginLeft': '5px'})
                ]),
                html.Div([
                    html.Span("◆", style={'color': self.CYTOSCAPE_COLORS['neutral'], 'fontSize': '20px'}),
                    html.Span(" Neutral", style={'marginLeft': '5px'})
                ]),
                html.Br(),
                html.P("Metabolites:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                html.Div([
                    html.Span("●", style={'color': self.CYTOSCAPE_COLORS['upregulated'], 'fontSize': '20px'}),
                    html.Span(" Upregulated", style={'marginLeft': '5px'})
                ]),
                html.Div([
                    html.Span("●", style={'color': self.CYTOSCAPE_COLORS['downregulated'], 'fontSize': '20px'}),
                    html.Span(" Downregulated", style={'marginLeft': '5px'})
                ]),
            ]
            
            # Add upstream/disease legend if data is present
            if self.upstream_data is not None and not self.upstream_data.empty:
                legend_items.extend([
                    html.Br(),
                    html.P("Upstream:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                    html.Div([
                        html.Span("▲", style={'color': self.CYTOSCAPE_COLORS['upstream'], 'fontSize': '20px'}),
                        html.Span(" Regulators", style={'marginLeft': '5px'})
                    ]),
                ])
            
            if self.disease_data is not None and not self.disease_data.empty:
                legend_items.extend([
                    html.Br(),
                    html.P("Diseases:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                    html.Div([
                        html.Span("⬢", style={'color': self.CYTOSCAPE_COLORS['disease'], 'fontSize': '20px', 
                                              'textShadow': f"0 0 0 {self.CYTOSCAPE_COLORS['disease']}"}),
                        html.Span(" Disease Nodes", style={'marginLeft': '5px'})
                    ]),
                ])
            
            # Add edge legend
            legend_items.extend([
                html.Br(),
                html.P("Connections:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                html.Div([
                    html.Span("⋯⋯", style={'color': self.CYTOSCAPE_COLORS['activated'], 'fontSize': '16px', 'fontWeight': 'bold'}),
                    html.Span(" Activating", style={'marginLeft': '5px', 'fontSize': '11px'})
                ]),
                html.Div([
                    html.Span("⋯⋯", style={'color': self.CYTOSCAPE_COLORS['inhibited'], 'fontSize': '16px', 'fontWeight': 'bold'}),
                    html.Span(" Inhibiting", style={'marginLeft': '5px', 'fontSize': '11px'})
                ]),
                html.Div([
                    html.Span("⋯⋯", style={'color': self.CYTOSCAPE_COLORS['association'], 'fontSize': '16px', 'fontWeight': 'bold'}),
                    html.Span(" Association", style={'marginLeft': '5px', 'fontSize': '11px'})
                ]),
            ])
            
            return html.Div(legend_items)
        else:
            # Simple mode legend - p-value with circles
            return html.Div([
                html.P("P-value:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                html.Div([
                    html.Span("●", style={'color': self.PVALUE_COLORS['very_sig'], 'fontSize': '20px'}),
                    html.Span(" p < 0.001", style={'marginLeft': '5px'})
                ]),
                html.Div([
                    html.Span("●", style={'color': self.PVALUE_COLORS['highly_sig'], 'fontSize': '20px'}),
                    html.Span(" p < 0.01", style={'marginLeft': '5px'})
                ]),
                html.Div([
                    html.Span("●", style={'color': self.PVALUE_COLORS['significant'], 'fontSize': '20px'}),
                    html.Span(" p < 0.05", style={'marginLeft': '5px'})
                ]),
                html.Div([
                    html.Span("●", style={'color': self.PVALUE_COLORS['marginal'], 'fontSize': '20px'}),
                    html.Span(" p < 0.1", style={'marginLeft': '5px'})
                ]),
                html.Div([
                    html.Span("●", style={'color': self.PVALUE_COLORS['not_sig'], 'fontSize': '20px'}),
                    html.Span(" p ≥ 0.1", style={'marginLeft': '5px'})
                ]),
            ])
    
    def set_data(self, pathway_data: pd.DataFrame):
        """Update the pathway data and rebuild the network."""
        self.pathway_data = pathway_data
        self.G = nx.Graph()
        self.elements = []
        self._build_network()
    
    def _create_app(self):
        """Create the Dash application."""
        if not DASH_AVAILABLE:
            raise ImportError("Dash is not installed. Install with: pip install dash dash-cytoscape")
        
        # Load extra layouts
        cyto.load_extra_layouts()
        
        # Create Dash app
        self.app = dash.Dash(__name__, title=self.title)
        
        # Stylesheet for Cytoscape - supports both simple and cytoscape modes
        default_stylesheet = [
            # Pathway nodes (diamond in cytoscape, circle in simple)
            {
                'selector': 'node[node_type="pathway"]',
                'style': {
                    'label': 'data(label)',
                    'width': 'data(size)',
                    'height': 'data(size)',
                    'background-color': 'data(color)',
                    'shape': 'data(shape)',  # Shape from node data
                    'border-width': 0,
                    'border-color': 'transparent',
                    'font-size': '16px',
                    'font-weight': 'bold',
                    'text-valign': 'center',
                    'text-halign': 'center',
                    'text-wrap': 'wrap',
                    'text-max-width': '200px',
                    'color': '#fff',
                    'text-outline-color': '#000',
                    'text-outline-width': 1,
                }
            },
            # Metabolite nodes (circle shape) - white text with black outline for visibility on any background
            {
                'selector': 'node[node_type="metabolite"]',
                'style': {
                    'label': 'data(label)',
                    'width': 'data(size)',
                    'height': 'data(size)',
                    'background-color': 'data(color)',
                    'shape': 'ellipse',
                    'border-width': 0,
                    'border-color': 'transparent',
                    'font-size': '13px',
                    'font-weight': 'bold',
                    'text-valign': 'center',
                    'text-halign': 'center',
                    'text-wrap': 'wrap',
                    'text-max-width': '150px',
                    'color': '#fff',
                    'text-outline-color': '#000',
                    'text-outline-width': 1,
                }
            },
            # Upstream nodes (triangle shape)
            {
                'selector': 'node[node_type=\"upstream\"]',
                'style': {
                    'label': 'data(label)',
                    'width': 'data(size)',
                    'height': 'data(size)',
                    'background-color': 'data(color)',
                    'shape': 'triangle',
                    'border-width': 2,
                    'border-color': '#333',
                    'font-size': '14px',
                    'font-weight': 'bold',
                    'text-valign': 'center',
                    'text-halign': 'center',
                    'text-wrap': 'wrap',
                    'text-max-width': '150px',
                    'color': '#fff',
                }
            },
            # Disease nodes (hexagon shape)
            {
                'selector': 'node[node_type=\"disease\"]',
                'style': {
                    'label': 'data(label)',
                    'width': 'data(size)',
                    'height': 'data(size)',
                    'background-color': 'data(color)',
                    'shape': 'hexagon',
                    'border-width': 2,
                    'border-color': '#333',
                    'font-size': '14px',
                    'font-weight': 'bold',
                    'text-valign': 'center',
                    'text-halign': 'center',
                    'text-wrap': 'wrap',
                    'text-max-width': '150px',
                    'color': '#fff',
                }
            },
            # Selected node
            {
                'selector': 'node:selected',
                'style': {
                    'border-width': 4,
                    'border-color': '#FFD700',
                    'overlay-color': '#FFD700',
                    'overlay-opacity': 0.3,
                }
            },
            # Edge styles - default gray dashed
            {
                'selector': 'edge',
                'style': {
                    'width': 4,
                    'line-color': '#666',
                    'opacity': 0.8,
                    'curve-style': 'bezier',
                    'line-style': 'dashed',
                    'line-dash-pattern': [8, 8],
                }
            },
            # Edges with color (full network mode - activation/inhibition/association)
            {
                'selector': 'edge[color]',
                'style': {
                    'width': 4,
                    'line-color': 'data(color)',
                    'opacity': 0.85,
                    'curve-style': 'bezier',
                    'line-style': 'dashed',
                    'line-dash-pattern': [8, 8],
                }
            },
            # Selected edge
            {
                'selector': 'edge:selected',
                'style': {
                    'line-color': '#FFD700',
                    'opacity': 1,
                    'width': 5,
                }
            },
        ]
        
        # Layout
        self.app.layout = html.Div([
            # Header
            html.Div([
                html.H2(self.title, style={'margin': '10px', 'color': '#333'}),
            ], style={'backgroundColor': '#f8f9fa', 'borderBottom': '1px solid #ddd'}),
            
            # Main content
            html.Div([
                # Left panel - Controls
                html.Div([
                    html.H4("Controls", style={'marginTop': '10px'}),
                    
                    # Layout selector
                    html.Label("Layout Algorithm:"),
                    dcc.Dropdown(
                        id='layout-dropdown',
                        options=[{'label': v, 'value': k} for k, v in self.LAYOUTS.items()],
                        value='cose',
                        clearable=False,
                        style={'marginBottom': '15px'}
                    ),
                    
                    # Network Mode Switcher (NEW)
                    html.Label("Network Mode:", style={'fontWeight': 'bold'}),
                    dcc.RadioItems(
                        id='network-mode',
                        options=[
                            {'label': '🔵 Pathway Similarity', 'value': 'simple'},
                            {'label': '🔶 Pathway-Metabolite Network', 'value': 'full_network'},
                        ],
                        value=self.mode,  # Start with current mode
                        style={'marginBottom': '10px'}
                    ),
                    
                    # Focused Mode checkbox (only visible in full_network mode)
                    html.Div([
                        dcc.Checklist(
                            id='focused-mode-checkbox',
                            options=[{'label': ' 🎯 Focused Mode', 'value': 'focused'}],
                            value=['focused'] if self.focused_mode else [],  # Default checked
                            style={'marginBottom': '5px', 'fontWeight': 'bold'}
                        ),
                        html.Div([
                            html.Span("ℹ️ ", style={'fontSize': '11px'}),
                            html.Span("Focused: Show only pathway-linked metabolites", 
                                     style={'fontSize': '10px', 'color': '#666'})
                        ], style={'marginLeft': '22px', 'marginBottom': '10px'}),
                    ], id='focused-mode-container', style={'display': 'block' if self.mode == 'full_network' else 'none'}),
                    
                    html.Hr(),
                    
                    # Background color
                    html.Label("Background:"),
                    dcc.Dropdown(
                        id='background-dropdown',
                        options=[
                            {'label': 'Dark', 'value': '#1a1a2e'},
                            {'label': 'Light', 'value': '#ffffff'},
                            {'label': 'Gray', 'value': '#2d2d2d'},
                        ],
                        value='#1a1a2e',
                        clearable=False,
                        style={'marginBottom': '15px'}
                    ),
                    
                    # Node size slider
                    html.Label("Node Size:"),
                    dcc.Slider(
                        id='node-size-slider',
                        min=0.2,
                        max=2.0,
                        step=0.1,
                        value=1.0,
                        marks={0.2: '0.2x', 0.5: '0.5x', 1.0: '1x', 1.5: '1.5x', 2.0: '2x'},
                        tooltip={'placement': 'bottom', 'always_visible': False},
                    ),
                    html.Br(),
                    
                    # Font size slider
                    html.Label("Font Size:"),
                    dcc.Slider(
                        id='font-size-slider',
                        min=8,
                        max=40,
                        step=1,
                        value=14,
                        marks={8: '8px', 16: '16px', 24: '24px', 32: '32px', 40: '40px'},
                        tooltip={'placement': 'bottom', 'always_visible': False},
                    ),
                    
                    html.Hr(),
                    
                    # Always show with metabolites - removed view mode selector
                    
                    # Filter by significance
                    html.Label("P-value Filter:"),
                    dcc.Slider(
                        id='pvalue-slider',
                        min=-5,
                        max=0,
                        step=0.5,
                        value=0,
                        marks={-5: '1e-5', -4: '1e-4', -3: '0.001', -2: '0.01', -1: '0.1', 0: '1.0'},
                        tooltip={'placement': 'bottom', 'always_visible': False},
                    ),
                    
                    html.Hr(),
                    
                    # Export section with actual download button
                    html.Label("Export Network:", style={'fontWeight': 'bold'}),
                    html.Button(
                        '📷 Download PNG',
                        id='btn-get-png',
                        n_clicks=0,
                        style={
                            'width': '100%',
                            'backgroundColor': '#27ae60',
                            'color': 'white',
                            'border': 'none',
                            'padding': '8px',
                            'borderRadius': '4px',
                            'cursor': 'pointer',
                            'marginBottom': '5px',
                            'fontWeight': 'bold',
                        }
                    ),
                    html.Button(
                        '📄 Download SVG',
                        id='btn-get-svg',
                        n_clicks=0,
                        style={
                            'width': '100%',
                            'backgroundColor': '#3498db',
                            'color': 'white',
                            'border': 'none',
                            'padding': '8px',
                            'borderRadius': '4px',
                            'cursor': 'pointer',
                            'marginBottom': '10px',
                            'fontWeight': 'bold',
                        }
                    ),
                    dcc.Download(id='download-image'),
                    html.Div(id='export-status', style={'fontSize': '10px', 'color': '#666', 'marginBottom': '10px'}),
                    
                    html.Hr(),
                    
                    # Selected node info (legend removed - now only in floating panel)
                    html.H4("Selected Node"),
                    html.Div(id='node-info', children="Click a node to see details"),
                    
                ], style={
                    'width': '250px',
                    'padding': '10px',
                    'backgroundColor': '#f8f9fa',
                    'borderRight': '1px solid #ddd',
                    'overflowY': 'auto',
                    'height': 'calc(100vh - 60px)',
                }),
                
                # Right panel - Network visualization with floating legend
                html.Div([
                    cyto.Cytoscape(
                        id='cytoscape-network',
                        elements=self.elements,
                        layout={'name': 'cose', 'animate': True, 'animationDuration': 500},
                        style={'width': '100%', 'height': 'calc(100vh - 60px)'},
                        stylesheet=default_stylesheet,
                        boxSelectionEnabled=True,
                        autoungrabify=False,
                        userZoomingEnabled=True,
                        userPanningEnabled=True,
                        minZoom=0.1,
                        maxZoom=3,
                    ),
                    # Floating legend overlay (top-right corner)
                    html.Div([
                        html.Div(id='floating-legend-content', children=self._get_legend_html(),
                                style={'backgroundColor': 'rgba(248, 249, 250, 0.95)',
                                       'padding': '10px',
                                       'borderRadius': '8px',
                                       'border': '1px solid #ddd',
                                       'fontSize': '11px',
                                       'maxWidth': '200px',
                                       'boxShadow': '0 2px 8px rgba(0,0,0,0.15)'})
                    ], style={'position': 'absolute',
                             'top': '15px',
                             'right': '15px',
                             'zIndex': '1000',
                             'pointerEvents': 'none'}),  # Allow clicks to pass through to network
                ], style={'flex': '1', 'position': 'relative'}),
                
            ], style={'display': 'flex', 'height': 'calc(100vh - 60px)'}),
            
        ], style={'fontFamily': 'Arial, sans-serif'})
        
        # Callbacks
        @self.app.callback(
            Output('cytoscape-network', 'layout'),
            Input('layout-dropdown', 'value')
        )
        def update_layout(layout_name):
            return {'name': layout_name, 'animate': True, 'animationDuration': 500}
        
        @self.app.callback(
            Output('cytoscape-network', 'style'),
            Input('background-dropdown', 'value')
        )
        def update_background(bg_color):
            return {'width': '100%', 'height': 'calc(100vh - 60px)', 'backgroundColor': bg_color}
        
        @self.app.callback(
            Output('node-info', 'children'),
            Input('cytoscape-network', 'tapNodeData')
        )
        def display_node_info(data):
            if data is None:
                return "Click a node to see details"
            
            node_type = data.get('node_type', 'pathway')
            
            if node_type == 'pathway':
                # Display pathway info
                info = [
                    html.P([html.B("Type: "), "Pathway"]),
                    html.P([html.B("Name: "), data.get('full_name', 'N/A')]),
                    html.P([html.B("Hits: "), str(data.get('hits', 'N/A'))]),
                ]
                
                if self.mode == 'full_network':
                    info.extend([
                        html.P([html.B("Status: "), data.get('status', 'N/A')]),
                        html.P([html.B("Z-score: "), data.get('z_score', 'N/A')]),
                    ])
                
                info.extend([
                    html.P([html.B("P-value: "), data.get('pvalue', 'N/A')]),
                ])
                
                if 'fdr' in data:
                    info.append(html.P([html.B("FDR: "), data.get('fdr', 'N/A')]))
                
                # Only show metabolites list in simple mode (in full network mode they're visible as nodes)
                if self.mode == 'simple':
                    info.extend([
                        html.Hr(),
                        html.P([html.B("Metabolites:")]),
                        html.P(data.get('metabolites', 'N/A'), style={'fontSize': '11px', 'wordWrap': 'break-word'}),
                    ])
                
                return html.Div(info)
                
            elif node_type == 'metabolite':
                # Display metabolite info
                return html.Div([
                    html.P([html.B("Type: "), "Compound"]),
                    html.P([html.B("Name: "), data.get('full_name', 'N/A')]),
                    html.P([html.B("Log2FC: "), data.get('log2fc', 'N/A')]),
                    html.P([html.B("Regulation: "), data.get('regulation', 'N/A')]),
                    html.P([html.B("P-value: "), data.get('pvalue', 'N/A')]),
                ])
            elif node_type == 'upstream':
                # Display upstream regulator info
                return html.Div([
                    html.P([html.B("Type: "), "Upstream Regulator"]),
                    html.P([html.B("Name: "), data.get('full_name', 'N/A')]),
                    html.P([html.B("Category: "), data.get('type', 'N/A')]),
                    html.P([html.B("Metabolites: "), str(data.get('n_metabolites', 0))]),
                ])
            elif node_type == 'disease':
                # Display disease info
                return html.Div([
                    html.P([html.B("Type: "), "Disease"]),
                    html.P([html.B("Name: "), data.get('full_name', 'N/A')]),
                    html.P([html.B("Metabolites: "), str(data.get('n_metabolites', 0))]),
                ])
            else:
                return html.Div([html.P("Node details not available")])
        
        @self.app.callback(
            Output('export-status', 'children'),
            Input('btn-get-png', 'n_clicks'),
            Input('btn-get-svg', 'n_clicks'),
            prevent_initial_call=True
        )
        def export_network(png_clicks, svg_clicks):
            """Trigger image generation."""
            ctx = callback_context
            if not ctx.triggered:
                return ""
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]
            if button_id == 'btn-get-png':
                return "📸 Downloading PNG image..."
            elif button_id == 'btn-get-svg':
                return "📄 Downloading SVG image..."
            return ""
        
        # Combined callback for PNG/SVG download using generateImage
        @self.app.callback(
            Output('cytoscape-network', 'generateImage'),
            [Input('btn-get-png', 'n_clicks'),
             Input('btn-get-svg', 'n_clicks')],
            State('background-dropdown', 'value'),
            prevent_initial_call=True
        )
        def download_image(png_clicks, svg_clicks, bg_color):
            """Generate PNG or SVG image for download based on which button was clicked."""
            ctx = callback_context
            if not ctx.triggered:
                return {}
            
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]
            
            if button_id == 'btn-get-png':
                return {
                    'type': 'png',
                    'action': 'download',
                    'filename': 'pathway_network',
                    'options': {
                        'bg': bg_color or '#1a1a2e',
                        'full': True,
                        'scale': 2
                    }
                }
            elif button_id == 'btn-get-svg':
                return {
                    'type': 'svg',
                    'action': 'download',
                    'filename': 'pathway_network',
                    'options': {
                        'full': True,
                        'scale': 2
                    }
                }
            return {}
        
        # Callback for updating stylesheet based on node size and font size
        @self.app.callback(
            Output('cytoscape-network', 'stylesheet'),
            [Input('node-size-slider', 'value'),
             Input('font-size-slider', 'value'),
             Input('background-dropdown', 'value')]
        )
        def update_stylesheet(node_scale, font_size, bg_color):
            """Update node sizes and font colors based on background."""
            # Determine font color based on background
            # Light backgrounds (#ffffff) get dark text, dark backgrounds get light text
            is_light_bg = bg_color in ('#ffffff', '#FFFFFF', 'white')
            
            # Font colors - all change based on background
            pathway_font_color = '#000000' if is_light_bg else '#ffffff'
            metabolite_font_color = '#000000' if is_light_bg else '#ffffff'
            upstream_disease_font_color = '#000000' if is_light_bg else '#ffffff'
            
            # Text outline for readability - contrast with font color
            # Light bg: black text with no/minimal outline
            # Dark bg: white text with black outline
            text_outline_color = '#ffffff' if is_light_bg else '#000000'
            text_outline_width = 0 if is_light_bg else 1
            metabolite_outline_color = '#ffffff' if is_light_bg else '#000000'
            metabolite_outline_width = 0 if is_light_bg else 1
            # For upstream/disease nodes
            upstream_disease_outline_color = '#ffffff' if is_light_bg else '#000000'
            upstream_disease_outline_width = 0 if is_light_bg else 1
            
            return [
                # Pathway nodes (diamond in cytoscape, circle in simple) - no border
                {
                    'selector': 'node[node_type="pathway"]',
                    'style': {
                        'label': 'data(label)',
                        'width': f'mapData(size, 30, 100, {30 * node_scale}, {100 * node_scale})',
                        'height': f'mapData(size, 30, 100, {30 * node_scale}, {100 * node_scale})',
                        'background-color': 'data(color)',
                        'shape': 'data(shape)',
                        'border-width': 0,
                        'border-color': 'transparent',
                        'font-size': f'{font_size}px',
                        'font-weight': 'bold',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'text-wrap': 'wrap',
                        'text-max-width': '200px',
                        'color': pathway_font_color,
                        'text-outline-color': text_outline_color,
                        'text-outline-width': text_outline_width,
                    }
                },
                # Metabolite nodes (circle shape) - dynamic text color based on background
                {
                    'selector': 'node[node_type="metabolite"]',
                    'style': {
                        'label': 'data(label)',
                        'width': f'{30 * node_scale}',
                        'height': f'{30 * node_scale}',
                        'background-color': 'data(color)',
                        'shape': 'ellipse',
                        'border-width': 0,
                        'border-color': 'transparent',
                        'font-size': f'{max(font_size - 3, 8)}px',
                        'font-weight': 'bold',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'text-wrap': 'wrap',
                        'text-max-width': '150px',
                        'color': metabolite_font_color,
                        'text-outline-color': metabolite_outline_color,
                        'text-outline-width': metabolite_outline_width,
                    }
                },
                # Upstream nodes (triangle shape) - no border
                {
                    'selector': 'node[node_type="upstream"]',
                    'style': {
                        'label': 'data(label)',
                        'width': f'{35 * node_scale}',
                        'height': f'{35 * node_scale}',
                        'background-color': 'data(color)',
                        'shape': 'triangle',
                        'border-width': 0,
                        'border-color': 'transparent',
                        'font-size': f'{font_size}px',
                        'font-weight': 'bold',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'text-wrap': 'wrap',
                        'text-max-width': '150px',
                        'color': upstream_disease_font_color,
                        'text-outline-color': upstream_disease_outline_color,
                        'text-outline-width': upstream_disease_outline_width,
                    }
                },
                # Disease nodes (hexagon shape) - no border
                {
                    'selector': 'node[node_type="disease"]',
                    'style': {
                        'label': 'data(label)',
                        'width': f'{38 * node_scale}',
                        'height': f'{38 * node_scale}',
                        'background-color': 'data(color)',
                        'shape': 'hexagon',
                        'border-width': 0,
                        'border-color': 'transparent',
                        'font-size': f'{font_size}px',
                        'font-weight': 'bold',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'text-wrap': 'wrap',
                        'text-max-width': '150px',
                        'color': upstream_disease_font_color,
                        'text-outline-color': upstream_disease_outline_color,
                        'text-outline-width': upstream_disease_outline_width,
                    }
                },
                # Selected node
                {
                    'selector': 'node:selected',
                    'style': {
                        'border-width': 4,
                        'border-color': '#FFD700',
                        'overlay-color': '#FFD700',
                        'overlay-opacity': 0.3,
                    }
                },
                # Edge styles - default gray dashed
                {
                    'selector': 'edge',
                    'style': {
                        'width': 4,
                        'line-color': '#666',
                        'opacity': 0.8,
                        'curve-style': 'bezier',
                        'line-style': 'dashed',
                        'line-dash-pattern': [8, 8],
                    }
                },
                # Edges with color (full network mode - activation/inhibition/association)
                {
                    'selector': 'edge[color]',
                    'style': {
                        'width': 4,
                        'line-color': 'data(color)',
                        'opacity': 0.85,
                        'curve-style': 'bezier',
                        'line-style': 'dashed',
                        'line-dash-pattern': [8, 8],
                    }
                },
                # Selected edge
                {
                    'selector': 'edge:selected',
                    'style': {
                        'line-color': '#FFD700',
                        'opacity': 1,
                        'width': 5,
                    }
                },
            ]
        
        @self.app.callback(
            Output('floating-legend-content', 'children'),
            Input('network-mode', 'value')
        )
        def update_floating_legend(network_mode):
            """Update floating legend when network mode changes."""
            # Don't update self.mode here - let filter_elements handle it
            # Otherwise the mode will already be changed when filter_elements checks it
            return self._get_legend_html(mode=network_mode)
        
        # View mode selector removed - always using 'with_mets' mode
        
        @self.app.callback(
            Output('focused-mode-container', 'style'),
            Input('network-mode', 'value')
        )
        def toggle_focused_mode_visibility(network_mode):
            """Show/hide focused mode checkbox based on network mode."""
            if network_mode == 'full_network':
                return {'display': 'block'}
            else:
                return {'display': 'none'}
        
        @self.app.callback(
            Output('cytoscape-network', 'elements'),
            [Input('pvalue-slider', 'value'),
             Input('network-mode', 'value'),
             Input('focused-mode-checkbox', 'value')]  # Add focused mode
        )
        def filter_elements(log_pvalue_threshold, network_mode, focused_mode_value):
            """Filter nodes by p-value threshold, network mode, and focused mode."""
            print(f"\n📊 filter_elements callback triggered:")
            print(f"   Current mode: {self.mode}")
            print(f"   Requested network_mode: {network_mode}")
            print(f"   focused_mode_value: {focused_mode_value}")
            print(f"   Elements before filtering: {len(self.elements)}")
            
            # Determine if focused mode is enabled
            new_focused_mode = 'focused' in (focused_mode_value or [])
            
            # Rebuild network if mode changed or focused mode changed
            need_rebuild = False
            if network_mode != self.mode:
                print(f"   🔄 Mode changed from '{self.mode}' to '{network_mode}'")
                self.mode = network_mode
                need_rebuild = True
            if new_focused_mode != self.focused_mode:
                print(f"   🔄 Focused mode changed from {self.focused_mode} to {new_focused_mode}")
                self.focused_mode = new_focused_mode
                need_rebuild = True
            
            if need_rebuild:
                print(f"   🔄 Rebuilding network...")
                self._build_network()  # Rebuild with new mode
                print(f"   Elements after rebuild: {len(self.elements)}")
            
            # Now filter by p-value and view mode
            pvalue_threshold = 10 ** log_pvalue_threshold
            print(f"   P-value threshold: {pvalue_threshold}")
            
            filtered = []
            visible_nodes = set()
            
            # Filter nodes
            for elem in self.elements:
                if 'source' not in elem['data']:  # It's a node
                    node_type = elem['data'].get('node_type', 'pathway')
                    
                    # Always show all node types (metabolites, upstream, diseases) with pathways
                    # Filter by p-value (only for pathways)
                    if node_type == 'pathway':
                        try:
                            pval = float(elem['data'].get('pvalue', '1').replace('e', 'E').replace('E', 'e'))
                        except:
                            pval = 1.0
                        
                        if pval <= pvalue_threshold:
                            filtered.append(elem)
                            visible_nodes.add(elem['data']['id'])
                    else:
                        # Metabolites, upstream, and disease nodes always visible if view mode allows
                        filtered.append(elem)
                        visible_nodes.add(elem['data']['id'])
            
            # Add edges only between visible nodes
            for elem in self.elements:
                if 'source' in elem['data']:  # It's an edge
                    if elem['data']['source'] in visible_nodes and elem['data']['target'] in visible_nodes:
                        filtered.append(elem)
            
            node_types_in_filtered = {}
            for elem in filtered:
                if 'source' not in elem['data']:
                    node_type = elem['data'].get('node_type', 'unknown')
                    node_types_in_filtered[node_type] = node_types_in_filtered.get(node_type, 0) + 1
            
            print(f"   ✅ Returning {len(filtered)} elements to display:")
            print(f"      Node types: {node_types_in_filtered}")
            print(f"      Edges: {len([e for e in filtered if 'source' in e['data']])}")
            
            return filtered
    
    def run(self, port: int = None, open_browser: bool = True, debug: bool = False):
        """
        Run the network viewer in a web browser.
        
        Args:
            port: Port to run on (auto-finds free port if None)
            open_browser: Whether to automatically open browser
            debug: Enable Dash debug mode
        """
        if not DASH_AVAILABLE:
            print("❌ Dash is not installed. Install with: pip install dash dash-cytoscape")
            return
        
        self._create_app()
        self.port = port or find_free_port()
        
        print(f"\n🌐 Starting Pathway Network Viewer...")
        print(f"   Open in browser: http://127.0.0.1:{self.port}")
        print(f"   Press Ctrl+C to stop\n")
        
        if open_browser:
            threading.Timer(1.5, lambda: webbrowser.open(f'http://127.0.0.1:{self.port}')).start()
        
        self.app.run(debug=debug, port=self.port, host='127.0.0.1')
    
    def run_threaded(self, port: int = None, open_browser: bool = True) -> int:
        """
        Run the network viewer in a background thread.
        
        Returns:
            The port number the server is running on
        """
        if not DASH_AVAILABLE:
            print("❌ Dash is not installed. Install with: pip install dash dash-cytoscape")
            return None
        
        self._create_app()
        self.port = port or find_free_port()
        
        def run_server():
            import logging
            log = logging.getLogger('werkzeug')
            log.setLevel(logging.ERROR)
            self.app.run(debug=False, port=self.port, host='127.0.0.1', use_reloader=False)
        
        self.server_thread = threading.Thread(target=run_server, daemon=True)
        self.server_thread.start()
        
        time.sleep(1)  # Give server time to start
        
        print(f"🌐 Network viewer running at: http://127.0.0.1:{self.port}")
        
        if open_browser:
            webbrowser.open(f'http://127.0.0.1:{self.port}')
        
        return self.port
    
    def run_embedded(self, parent_window=None, width: int = 1200, height: int = 800):
        """
        Run the network viewer embedded in a window (requires pywebview).
        
        Args:
            parent_window: Parent tkinter window (optional)
            width: Window width
            height: Window height
        """
        if not WEBVIEW_AVAILABLE:
            print("⚠️ pywebview not installed. Falling back to browser mode.")
            print("   Install with: pip install pywebview")
            self.run(open_browser=True)
            return
        
        port = self.run_threaded(open_browser=False)
        
        if port:
            webview.create_window(
                self.title,
                f'http://127.0.0.1:{port}',
                width=width,
                height=height,
            )
            webview.start()


def create_demo_data() -> pd.DataFrame:
    """Create demo pathway data for testing."""
    import random
    
    pathways = [
        ("Phenylalanine, tyrosine and tryptophan biosynthesis", 2),
        ("Fatty acid biosynthesis", 3),
        ("Lysine degradation", 3),
        ("Purine metabolism", 12),
        ("Arginine biosynthesis", 4),
        ("Arginine and proline metabolism", 7),
        ("Amino sugar and nucleotide sugar metabolism", 3),
        ("Alanine, aspartate and glutamate metabolism", 4),
        ("Tyrosine metabolism", 2),
        ("Vitamin B6 metabolism", 1),
        ("Glycerolipid metabolism", 1),
        ("Fatty acid elongation", 1),
        ("Ubiquinone biosynthesis", 2),
        ("Glyoxylate metabolism", 2),
        ("Starch and sucrose metabolism", 2),
    ]
    
    data = []
    for name, hits in pathways:
        pval = random.uniform(1e-13, 1e-5) if hits > 3 else random.uniform(1e-9, 1e-6)
        fdr = pval * random.uniform(1, 10)
        mets = "; ".join([f"Met{i}" for i in range(hits)])
        data.append({
            'Pathway': name,
            'Hits': hits,
            'P-Value': pval,
            'FDR': fdr,
            'Metabolites': mets,
        })
    
    return pd.DataFrame(data)


# Simple tkinter integration example
def show_network_in_tkinter(pathway_data: pd.DataFrame, parent=None):
    """
    Show network in a new window from tkinter.
    Opens in system browser with dark theme like MetaboAnalyst.
    
    Args:
        pathway_data: DataFrame with pathway enrichment results
        parent: Optional tkinter parent window
    """
    viewer = PathwayNetworkViewer(
        pathway_data=pathway_data,
        title="Pathway Enrichment Network",
        show_edges=True,
        min_edge_weight=0.05,
    )
    
    # Run in background thread and open browser
    viewer.run_threaded(open_browser=True)
    
    return viewer


if __name__ == "__main__":
    # Demo mode
    print("=" * 60)
    print("Interactive Pathway Network Viewer - Demo")
    print("=" * 60)
    
    if not DASH_AVAILABLE:
        print("\n❌ Required packages not installed!")
        print("   Run: pip install dash dash-cytoscape")
        exit(1)
    
    # Create demo data
    demo_df = create_demo_data()
    print(f"\nDemo data created with {len(demo_df)} pathways")
    print(demo_df[['Pathway', 'Hits', 'P-Value']].head(10))
    
    # Launch viewer
    viewer = PathwayNetworkViewer(
        pathway_data=demo_df,
        title="Demo: Pathway Enrichment Network",
    )
    viewer.run(open_browser=True)
