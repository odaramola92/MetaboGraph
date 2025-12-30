# Comparative Analysis Tab - User Guide

## Overview
The Comparative Analysis tab allows you to compare **pathway-level outputs** across multiple datasets. It takes the pathway enrichment results exported from the Network Analysis tab and visualizes how pathways compare across different experiments, conditions, or time points.

## Key Purpose
This tab is **NOT** for comparing raw metabolite data or running statistical tests like t-tests. Instead, it compares **pathway enrichment results** from multiple studies/conditions to identify:
- Pathways consistently enriched across datasets
- Differences in pathway activation between conditions
- Overlap of significant pathways

## Features

### 1. Upload Pathway Exports
- Upload pathway export Excel files from the Network Analysis tab
- Each file should contain the "Selected_Pathways" sheet with enrichment results
- Required columns: Pathway name, P-value (and optionally Z-score, Mean_Log2FC)

### 2. Visualization Types

#### Z-Score Heatmap (Default)
- **Color Scale**: Blue (inhibited) → Gray (neutral) → Red (activated)
- **Rows**: Pathways sorted by max absolute Z-score
- **Columns**: Each uploaded dataset
- **White cells**: Pathway absent in that dataset
- Shows pathway activation/inhibition direction

#### P-Value Heatmap
- **Color Scale**: Yellow-Orange-Red (-log10 p-value)
- **Higher intensity**: More significant p-values
- Shows statistical significance across datasets

#### Pathway Overlap Bar Chart
- Shows count of significant pathways per dataset
- Compares pathway coverage across conditions

#### Venn Diagram
- Shows overlap of significant pathways between 2-3 datasets
- Identifies shared vs. unique pathways

### 3. Settings

| Setting | Description | Default |
|---------|-------------|---------|
| P-value threshold | Significance cutoff | 0.05 |
| Top N pathways | Maximum pathways to display | 50 |
| Pathway filter | All or Common pathways only | All |
| Figure size | Width x Height in inches | 10 x 8 |

## Workflow

### Step 1: Export from Network Tab
First, for each condition/experiment you want to compare:
1. Complete pathway analysis in the Pathway Network tab
2. Export the pathway results (Excel workbook with "Selected_Pathways" sheet)

### Step 2: Upload Pathway Exports
1. Click **"Add Workbook(s)"**
2. Select one or more pathway export files
3. Give each dataset a descriptive name (e.g., "Control", "Treatment", "24h")
4. Files appear in the list showing rows and columns

### Step 3: Configure Analysis
1. Set **P-value threshold** (default 0.05)
2. Choose **Visualization type** (Z-Score Heatmap recommended)
3. Set **Top N pathways** to limit display
4. Optionally filter to **Common Pathways Only**

### Step 4: Run Analysis
1. Click **"Run Comparative Analysis"**
2. System aggregates pathways across all datasets
3. Calculates overlap statistics
4. Generates visualizations

### Step 5: Interpret Results
- **Z-Score Heatmap**: 
  - Red = pathway activated (positive Z-score)
  - Blue = pathway inhibited (negative Z-score)
  - Gray = neutral (Z-score ≈ 0)
  - White = pathway not present in dataset
- **Common pathways** appear at top, sorted by significance
- **Log panel** shows statistics: significant pathways per dataset, overlap counts

### Step 6: Export
- **Export Excel**: Summary table with all pathway statistics
- **Save PNG**: Raster image of current visualization
- **Save SVG**: Vector image for publications

## Expected Input Format

The uploaded pathway exports should have columns like:
```
| Pathway | Combined_Pvalue | Z_Score | Mean_Log2FC | N_Metabolites |
|---------|-----------------|---------|-------------|---------------|
| Glycolysis | 0.001 | 2.45 | 0.89 | 12 |
| TCA Cycle | 0.023 | 1.87 | 0.45 | 8 |
```

**Column Detection**:
- Pathway name: "pathway", "name" (auto-detected)
- P-value: "combined_pvalue", "pvalue", "p_value"
- Z-score: "z_score", "zscore" (optional)
- Fold change: "mean_log2fc", "log2fc" (optional)

## Example Use Cases

### Example 1: Compare Disease vs. Control
```
Datasets:
  - Control_pathways.xlsx (from healthy samples)
  - Disease_pathways.xlsx (from patient samples)

Analysis:
  - Z-Score heatmap shows which pathways are activated (red) or 
    inhibited (blue) in disease vs. control
  - Identify disease-specific pathway signatures
```

### Example 2: Time-Course Study
```
Datasets:
  - 0h_pathways.xlsx
  - 6h_pathways.xlsx  
  - 24h_pathways.xlsx

Analysis:
  - Track pathway activation over time
  - Identify early vs. late response pathways
```

### Example 3: Multi-Tissue Comparison
```
Datasets:
  - Liver_pathways.xlsx
  - Muscle_pathways.xlsx
  - Adipose_pathways.xlsx

Analysis:
  - Identify tissue-specific pathway enrichment
  - Find systemic pathways (significant in all tissues)
```

## Log Output

The analysis log shows:
```
📊 Pathway Statistics:
   • Dataset1: 45/120 significant (p ≤ 0.05)
   • Dataset2: 38/105 significant (p ≤ 0.05)

🔗 Pathways present in ALL datasets: 85
🔗 Pathways significant in ALL datasets: 12
```

## Tips & Best Practices

### Data Preparation
- ✅ Use consistent pathway databases across all datasets
- ✅ Use same organism for all comparisons
- ✅ Apply same enrichment method (ORA/GSEA)
- ✅ Export complete results (not just significant)

### Visualization
- ✅ Use Z-Score heatmap for direction of change
- ✅ Use P-Value heatmap for significance only
- ✅ Limit to top 30-50 pathways for readability
- ✅ Use "Common Pathways Only" to focus on shared biology

### Interpretation
- ✅ Look for consistent patterns across datasets
- ✅ Pay attention to pathway direction (activated vs. inhibited)
- ✅ White cells indicate missing data, not zero
- ✅ Consider biological context when interpreting

## Troubleshooting

### "No pathway-like column found"
- Ensure your export has a "Pathway" or similar column
- Check column names for typos

### "No Z-Score data available"
- Falls back to P-value heatmap
- Z-scores require enrichment analysis that generates Z-scores

### No significant pathways in overlap
- Check p-value threshold (try increasing to 0.1)
- Verify datasets are comparable
- Consider using "All Pathways" filter instead of "Common"

## Next Steps

After comparative analysis:
1. **Identify key pathways** consistently enriched across conditions
2. **Follow up** with targeted metabolite analysis
3. **Generate publication figures** using SVG export
4. **Document findings** in the exported Excel summary
