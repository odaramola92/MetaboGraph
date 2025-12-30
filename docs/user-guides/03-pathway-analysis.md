# Pathway Analysis - User Guide

## Overview
MetaboGraph provides **two tabs** for pathway analysis:
1. **Pathway Annotation Tab**: Maps metabolites/lipids to biological pathways
2. **Pathway Network Tab**: Enrichment analysis, network visualization, and Cytoscape integration

This guide covers both tabs as they work together in the analysis workflow.

---

# Part 1: Pathway Annotation Tab

## Purpose
The Pathway Annotation tab maps ID-annotated metabolites to biological pathways using multiple pathway databases (SMPDB, WikiPathways, PathBank).

## Key Features
- **Data Mode Selection**: Metabolite or Lipid
- **Multi-Database Mapping**: SMPDB, WikiPathways, PathBank
- **Species-Specific**: Human, Mouse, Rat
- **Project-Based Output**: Auto-saves with project name and timestamp

## Required Inputs
- **ID-Annotated Excel File**: Output from ID Annotation tab
- **Working Folder**: Where output files are saved
- **Project Name**: Used for output filename

## Setup Before Running

### 1. Project Setup (Required)
1. **Select Working Folder**: Click "📁 Browse" to choose output directory
2. **Enter Project Name**: Used in output filename
3. ✅ Green checkmark appears when configured

### 2. Data Mode Selection
- **Metabolite Mode**: Requires Feature ID column; uses HMDB, KEGG, PubChem, ChEBI, CAS, SMILES, InChI, InChIKey for mapping
- **Lipid Mode**: Uses LipidMaps ID and lipid class for mapping

### 3. Input File Selection
- Auto-loaded from ID Annotation tab OR
- Browse to select custom ID-annotated file
- File must have ID columns (HMDB_ID, KEGG_ID, etc.)

### 4. P-Value & Log2FC Assignment
For filtering and network analysis:
1. Click **"🔍 Verify P-Value"** - Select your significance column
2. Click **"🔍 Verify Log2FC"** - Select your fold change column

## Workflow

1. **Configure Project** (folder + name)
2. **Select Data Mode** (Metabolite or Lipid)
3. **Load Input File** (ID-annotated data)
4. **Verify Columns** (P-Value, Log2FC)
5. **Click "🚀 Start Pathway Annotation"**
6. **Monitor Progress** in log panel
7. **Auto-exports** to Network Tab when complete

## Output File
**Filename**: `[ProjectName]_pathway_annotated_[timestamp].xlsx`

**Contains**:
- All original columns
- HMDB_Pathways, KEGG_Pathways, SMPDB_Pathways, WikiPathways
- Disease associations
- Upstream regulators (enzymes)
- Merged All_Pathways column

---

# Part 2: Pathway Network Tab

## Purpose
The Pathway Network tab performs **Fisher's Over-Representation Analysis (ORA)**, visualizes pathway enrichment, and integrates with Cytoscape for network analysis.

## Key Features
- **Fisher's Exact Test**: Over-representation analysis with FDR correction
- **Interactive Network Viewer**: Dash-based visualization (MetaboAnalyst-like)
- **Multiple Pathway Selection**: Pathways, Upstream Regulators, Diseases
- **Cytoscape Integration**: Export networks for advanced analysis

## Data Loading

### Option 1: Auto-load from Pathway Annotation
After running pathway annotation, data is automatically pushed to Network tab.

### Option 2: Manual Import
1. Click **"📂 Browse & Import"**
2. Select pathway-annotated Excel file
3. System auto-detects columns

### Column Verification
1. Click **"🔍 Verify Columns"**
2. Confirm/assign:
   - **Feature ID**: Metabolite name column
   - **P-Value**: Significance column (prefers adjusted p-values)
   - **Log2 Fold Change**: Direction/magnitude of change
   - **Class** / **Class_name**: For lipid data

## Enrichment Analysis (Fisher's ORA)

### Configure Parameters
| Parameter | Description | Default |
|-----------|-------------|---------|
| P-Value Threshold | Significance cutoff | 0.05 |
| Log2FC Threshold | Fold change cutoff | 1.0 |
| Direction | Up/Down/Both | Both |
| Universe | Reference metabolite set | Auto-detect |

### Run Enrichment
1. Set filtering thresholds
2. Click **"📊 Run Fisher's Enrichment"**
3. Results populate in tables:
   - **Pathways Tab**: Enriched biological pathways
   - **Upstream Tab**: Enzyme regulators
   - **Diseases Tab**: Disease associations

### Enrichment Results Columns
| Column | Description |
|--------|-------------|
| Pathway | Pathway name |
| Hits | Metabolites in your data that map to pathway |
| Total | Total metabolites in pathway (database) |
| P-Value | Fisher's exact test p-value |
| FDR | Benjamini-Hochberg corrected p-value |
| Enrichment Ratio | Observed/Expected |

## Interactive Selection

### Selection Trees
Three tabs with hierarchical selection:
1. **Pathways**: Biological pathways
2. **Upstream**: Enzyme regulators
3. **Diseases**: Disease associations

### Selection Controls
- ☑️ **Check** pathways to include in network
- **Auto-Connect**: When enabled, checking a pathway highlights connected metabolites
- **Statistics**: Shows selected/total counts

### Quick Selection Actions
- **Select All**: Check all pathways
- **Deselect All**: Uncheck all
- **Select by Threshold**: Auto-select based on criteria

## Network Generation

### Build Network
1. Select pathways/upstream/diseases using checkboxes
2. Click **"🌐 Generate Network"**
3. Network generates with selected elements

### Network Export Options
| Format | Description | Use Case |
|--------|-------------|----------|
| JSON | Node/edge data | Web viewers |
| GraphML | XML format | Cytoscape |
| TSV | Tab-separated | Spreadsheets |
| PNG | Image | Presentations |

## Interactive Viewer (Dash)

### Launch Viewer
1. Select pathways
2. Click **"🔬 Launch Interactive Viewer"**
3. Browser opens with network visualization

### Viewer Features
- **Zoom/Pan**: Navigate large networks
- **Node Selection**: Click nodes for details
- **Layout Options**: Force-directed, hierarchical
- **Color Coding**: By fold change, significance, or category
- **Export**: Download as image

## Cytoscape Integration

### Export for Cytoscape
1. Generate network
2. Click **"🔄 Export to Cytoscape"**
3. Files saved: `.graphml` or `.json`
4. Open in Cytoscape for advanced analysis

### Cytoscape Format Options
- **GraphML**: Full network with attributes
- **SIF**: Simple interaction format
- **JSON**: Cytoscape.js compatible

## Output Files

### Network Export
**Folder**: Output folder specified in settings

**Files Generated**:
- `[ProjectName]_network.graphml` - Network file
- `[ProjectName]_nodes.tsv` - Node attributes
- `[ProjectName]_edges.tsv` - Edge list
- `[ProjectName]_enrichment_results.xlsx` - Enrichment tables

### Selected Pathways Export
- `Selected_Pathways` sheet added to output Excel
- Contains: Pathway name, P-value, FDR, Hits, metabolite list

---

## Workflow Summary

```
ID-Annotated Data
       ↓
┌──────────────────────┐
│ Pathway Annotation   │
│ Tab                  │
│ • Map to pathways    │
│ • Multi-database     │
│ • Add pathway cols   │
└──────────────────────┘
       ↓
Pathway-Annotated Data
       ↓
┌──────────────────────┐
│ Pathway Network Tab  │
│ • Fisher's ORA       │
│ • Enrichment results │
│ • Interactive viewer │
│ • Cytoscape export   │
└──────────────────────┘
       ↓
Network Visualization + Enrichment Results
```

---

## Tips & Best Practices

### Before Pathway Analysis
- ✅ Complete ID Annotation first (good coverage improves results)
- ✅ Verify P-Value and Log2FC columns are correct
- ✅ Set working folder for organized output
- ✅ Use meaningful project names

### Enrichment Analysis
- ✅ Use FDR-corrected p-values for multiple testing
- ✅ Start with standard thresholds (p<0.05, |FC|>1.0)
- ✅ Consider biological relevance alongside statistics
- ✅ Pathways with more hits are generally more reliable

### Network Visualization
- ✅ Select relevant pathways (not all significant ones)
- ✅ Use auto-connect to see metabolite relationships
- ✅ Export to Cytoscape for publication-quality figures
- ✅ Color by fold change direction for biological interpretation

### Interpretation
- ✅ Focus on pathways with multiple metabolite hits
- ✅ Consider pathway coverage (hits/total)
- ✅ Look for consistent patterns (related pathways enriched together)
- ✅ Validate key findings with literature

---

## Troubleshooting

### No Pathways Found
- **Cause**: Poor ID annotation or wrong data mode
- **Solution**: Check ID Annotation coverage, verify column assignments

### Low Enrichment Hits
- **Cause**: Stringent thresholds or small dataset
- **Solution**: Relax P-value/FC thresholds, check universe size

### Network Too Dense
- **Cause**: Too many pathways selected
- **Solution**: Select fewer pathways, use top enriched only

### Interactive Viewer Won't Launch
- **Cause**: Dash not installed or port conflict
- **Solution**: Check Dash installation, try different port

### Cytoscape Export Fails
- **Cause**: Invalid characters in pathway names
- **Solution**: System auto-sanitizes; check output folder permissions

---

## Related Tabs
- **ID Annotation Tab**: Provides input data with database IDs
- **Comparative Analysis Tab**: Compare pathway results across datasets
- **Multi-Omics Tab**: Merge datasets before network analysis

---

## Support
For questions or issues, refer to the Help tab or submit an issue on GitHub.
- Filter by fold-change threshold
- Show only top pathways
- Adjust edge confidence threshold
- Use hierarchical layout

## Advanced Features

### Custom Universe
```
1. Prepare universe file (list of all measurable metabolites)
2. Include HMDB IDs or KEGG IDs
3. Load via "Custom Universe" option
4. System uses this for background
```

### Pathway Hierarchy Navigation
- View super-pathways (broad categories)
- Drill down to sub-pathways (specific processes)
- Collapse/expand pathway trees
- Filter by hierarchy level

### Metabolite Set Analysis
- Import custom metabolite sets
- Run enrichment on multiple sets simultaneously
- Compare set overlaps
- Venn diagrams for set intersections

### Export for External Tools
**Compatible Formats**:
- **Metaboanalyst**: Upload enrichment results
- **Cytoscape**: Import network (GraphML, SIF)
- **Ingenuity Pathway Analysis (IPA)**: ID mapping
- **KEGG Mapper**: Pathway visualization
- **Reactome**: Pathway analysis

## Next Steps

After Pathway Analysis:
1. **Interpret Biology**: Connect enriched pathways to phenotype
2. **Literature Review**: Validate findings with publications
3. **Hypothesis Generation**: Identify targets for follow-up
4. **Comparative Analysis**: Compare across experimental conditions
5. **Integration**: Combine with transcriptomics/proteomics data (Multiomics tab)

## Additional Resources
- [Pathway Analysis Best Practices](./pathway-analysis-best-practices.md)
- [Statistical Methods Explained](./statistical-methods.md)
- [Database Coverage and Updates](./database-info.md)
- [Publication Guidelines](./publication-guidelines.md)

## Support
For pathway analysis questions:
- [FAQ - Pathway Analysis](../faq.md#pathway-analysis)
- [Interpretation Guide](./interpretation-guide.md)
- GitHub Issues: [github-link]
- Email: support@metabograph.org
