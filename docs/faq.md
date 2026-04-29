# Frequently Asked Questions (FAQ)

## Table of Contents
- [General Questions](#general-questions)
- [Installation & Setup](#installation--setup)
- [Data Requirements](#data-requirements)
- [Analysis & Features](#analysis--features)
- [Databases](#databases)
- [Results & Interpretation](#results--interpretation)
- [Troubleshooting](#troubleshooting)
- [Citing & Contributing](#citing--contributing)

---

## General Questions

### What is MetaboGraph?

MetaboGraph is a comprehensive desktop application for metabolomics and lipidomics data analysis. It provides tools for:
- Data cleaning and preprocessing
- Metabolite/lipid identification and annotation
- Pathway enrichment analysis
- Network visualization
- Comparative analysis across datasets
- Multi-omics integration

### Who is MetaboGraph for?

- **Metabolomics researchers** analyzing LC-MS data
- **Lipidomics researchers** working with lipid profiles
- **Systems biologists** integrating omics datasets
- **Students** learning metabolomics data analysis
- **Core facilities** providing analysis services

No programming experience required!

### Is MetaboGraph free?

Yes! MetaboGraph is open-source software licensed under the MIT License. You can use it freely for academic and commercial purposes.

Note: Some databases used by MetaboGraph have their own licenses. See [Database Licenses](user-guides/04-database-setup.md) for details.

### What platforms does MetaboGraph support?

Currently:
- ✅ **Windows 10/11** (fully tested, standalone executable available)
- ⚠️ **macOS/Linux** (source installation only, limited testing)

Future versions may include native macOS/Linux executables.

### How does MetaboGraph compare to other tools?

| Feature | MetaboGraph | MetaboAnalyst | XCMS Online | Compound Discoverer |
|---------|-------------|---------------|-------------|---------------------|
| **Offline** | ✅ | ❌ | ❌ | ✅ |
| **Free** | ✅ | ✅ | ✅ | ❌ ($$$) |
| **Lipid Annotation** | ✅ | ⚠️ | ❌ | ✅ |
| **Network Analysis** | ✅ | ⚠️ | ❌ | ⚠️ |
| **Multi-Omics** | ✅ | ⚠️ | ❌ | ⚠️ |
| **No File Upload** | ✅ | ❌ | ❌ | ✅ |
| **API Annotation** | ✅ | ❌ | ❌ | ✅ |

---

## Installation & Setup

### Do I need Python installed?

**For standalone executable:** No! The `.exe` version includes everything needed.

**For source installation:** Yes, Python 3.8 or higher required.

### How much disk space do I need?

- **Application:** ~250 MB
- **Databases:** ~3 GB
- **Working space:** 2-5 GB (for temporary files during analysis)
- **Total recommended:** 10+ GB free space

### How long does database setup take?

First-time setup: **15-30 minutes** depending on internet speed.

- Download: 5-10 minutes (fast connection)
- Processing: 10-20 minutes (one-time only)

Subsequent launches use cached databases instantly.

### Can I skip database download?

You can skip databases you won't use:
- Skip **PathBank/SMPDB** if you only want ID annotation (no pathway analysis)
- Skip **LipidMaps** if you're only analyzing metabolites (not lipids)

At minimum, download **HMDB** for basic metabolite annotation.

### Can I use my own databases?

Not directly in current version. Future versions may support custom database import.

Workaround: Manually edit `.feather` files in `Databases/` folder (advanced users only).

---

## Data Requirements

### What file formats are supported?

**Input:**
- Excel (`.xlsx`) - **primary format**
- CSV (`.csv`) - limited support
- Compound Discoverer exports (`.xlsx`)

**Output:**
- Excel (`.xlsx`) with multiple sheets
- CSV (`.csv`) for raw data
- PNG/SVG for plots

### What columns are required in my data?

**For Data Cleaning:**
- Compound names/IDs
- m/z values
- Retention time
- Optional: Intensities, ion types

**For ID Annotation:**
- Compound names (for searching) OR
- m/z values + retention time (for matching)

**For Pathway Analysis:**
- At least one ID column: HMDB ID, KEGG ID, or PubChem CID

See [User Guides](user-guides/) for detailed column requirements.

### Can I analyze Compound Discoverer output directly?

Yes! The Data Cleaning tab is specifically designed for Compound Discoverer exports:
1. Load `.xlsx` export directly
2. Select positive/negative mode
3. Process and export
4. Use cleaned data for annotation

### What organisms are supported?

- **Metabolite databases:** All organisms (HMDB, KEGG are multi-species)
- **Pathway databases:**
  - Human (*Homo sapiens*) - full support
  - Mouse (*Mus musculus*) - full support
  - Rat (*Rattus norvegicus*) - full support
  - Other species - limited pathway data

### How many metabolites can I analyze?

Tested with up to **10,000 features** per dataset. Performance depends on your system:
- **8 GB RAM:** Up to 2,000-3,000 features comfortably
- **16 GB RAM:** Up to 5,000-8,000 features
- **32+ GB RAM:** 10,000+ features

For larger datasets, consider filtering or batch processing.

---

## Analysis & Features

### What metabolite identification methods are used?

MetaboGraph uses a **multi-tier approach**:

1. **Primary:** API queries (PubChem, KEGG) for latest data
2. **Fallback:** Local HMDB database for offline/fast lookup
3. **Lipid-specific:** LipidMaps classification system

Advantages:
- Always up-to-date (via APIs)
- Works offline (via local databases)
- Handles lipids specially (class/subclass annotation)

### How accurate is the annotation?

Accuracy depends on data quality and search parameters:
- **Exact mass match** (<5 ppm): ~80-90% correct for common metabolites
- **Mass + name match:** >95% correct
- **Manual verification recommended** for critical identifications

False positives are possible - always validate important results manually.

### What pathway analysis methods are available?

1. **Over-Representation Analysis (ORA):**
   - Fisher's exact test
   - Identifies enriched pathways
   - Fast, straightforward interpretation

2. **Network Analysis:**
   - Pathway-pathway interactions
   - Visualize relationships
   - Community detection

3. **Comparative Analysis:**
   - Compare multiple datasets
   - Z-score heatmaps
   - Venn diagrams

Future versions may include GSEA (Gene Set Enrichment Analysis).

### Can I do statistical analysis?

MetaboGraph focuses on **annotation and pathway analysis**, not statistical testing.

For statistics (t-tests, ANOVA, PCA), use:
- MetaboAnalyst (online)
- R/Bioconductor
- Python (scipy, statsmodels)
- Excel/GraphPad Prism

Import your statistically analyzed data into MetaboGraph for pathway annotation.

### What visualizations are available?

- **Heatmaps:** Z-score, P-value, fold-change
- **Network graphs:** Pathway networks with interactive layout
- **Venn diagrams:** Dataset overlap (2-3 sets)
- **Bar charts:** Enrichment scores, pathway counts
- **Tables:** Searchable, sortable results

All plots exportable as high-resolution PNG or SVG.

---

## Databases

### Why are databases not included?

**File size:** Databases total ~3 GB, making the download too large.

**Licensing:** Some databases have redistribution restrictions.

**Updates:** Separate download allows easy updates without reinstalling app.

### How often should I update databases?

Recommended: **Every 6-12 months**

Databases are updated regularly by providers:
- **HMDB:** Annual major updates
- **KEGG:** Monthly updates
- **WikiPathways:** Weekly updates

To update: Re-run "Download All Databases" in Database Setup tab.

### Can I use older database versions?

Yes, but not recommended. Older versions may:
- Miss newly discovered metabolites
- Have outdated pathway information
- Contain errors that have been fixed

Keep databases updated for best results.

### What if a metabolite isn't in the database?

Options:
1. **Try different databases** - some metabolites only in KEGG, not HMDB
2. **Update databases** - newly discovered metabolites added regularly
3. **Check spelling** - typos prevent matching
4. **Try synonyms** - metabolites have multiple names
5. **Use m/z matching** - more permissive than name search

Some metabolites genuinely unknown - this is expected in untargeted metabolomics.

### Are database downloads secure?

Yes:
- Downloaded from official sources (HMDB, KEGG, etc.)
- HTTPS encrypted connections
- File integrity checks after download
- No executable code in database files

---

## Results & Interpretation

### What is a good P-value threshold?

**Standard:** P < 0.05 (with FDR correction)

**Exploratory:** P < 0.10 may be acceptable

**Conservative:** P < 0.01 for high confidence

Always apply multiple testing correction (Benjamini-Hochberg FDR) when testing many pathways.

### How do I interpret enrichment scores?

**Fold enrichment** = (Observed / Expected)
- **> 1:** More metabolites than expected (enriched)
- **< 1:** Fewer metabolites than expected (depleted)
- **~1:** No enrichment

Higher scores = stronger enrichment, but consider P-value too!

**Example:** Fold enrichment = 3.5, P = 0.001
- Pathway has 3.5× more significant metabolites than expected by chance
- Very high confidence (P < 0.001)

### What does "Universe" mean in pathway analysis?

**Universe** = all metabolites tested in your experiment.

Used for background correction in enrichment analysis. MetaboGraph automatically:
1. Detects all metabolites in your dataset
2. Maps them to pathways
3. Uses as background for statistical testing

Proper universe prevents false positives.

### Why do I get different results than MetaboAnalyst?

Possible reasons:
1. **Different databases** - versions/sources differ
2. **Different algorithms** - ORA implementation details vary
3. **Different ID mapping** - how names/IDs are matched
4. **Different universe** - background set differs

Both results can be valid - use biological knowledge to interpret.

### How many pathways should I report?

**Typical:** Top 10-20 significant pathways tell coherent story

**Maximum:** Top 50 for comprehensive analysis

**Avoid:** Reporting 100+ pathways - likely includes noise

Focus on pathways that:
- Have strong enrichment (high fold-change)
- Are statistically significant (P < 0.05)
- Make biological sense for your experiment

---

## Troubleshooting

### The application won't start. What do I do?

See [Troubleshooting Guide](troubleshooting.md) for detailed solutions.

Quick checks:
1. Run as Administrator
2. Check antivirus isn't blocking
3. Verify all files extracted
4. Check `logs/` folder for errors

### Analysis is taking forever. Is it frozen?

**Check if active:**
- Task Manager shows CPU usage
- Progress bar is updating
- Status messages changing

**Typical times:**
- ID annotation (1000 metabolites): 5-15 minutes
- Pathway annotation: 2-5 minutes
- Network generation: 1-3 minutes

If truly frozen (>30 min, no progress), restart application.

### I got an error message. What does it mean?

Check [Troubleshooting Guide](troubleshooting.md) → Common Error Messages section.

Most errors indicate:
- File format issues → Re-save as `.xlsx`
- Missing columns → Check required columns
- Database issues → Verify databases loaded
- Memory issues → Close other apps, use smaller dataset

### My results look wrong. How do I validate?

**Validation steps:**
1. **Spot check:** Manually verify a few metabolite annotations on HMDB.ca
2. **Known metabolites:** Include standards/controls, verify they annotate correctly
3. **Biological sense:** Do enriched pathways match your experimental design?
4. **Cross-reference:** Compare with literature, other tools
5. **Repeat analysis:** Use slightly different parameters, check consistency

**Red flags:**
- 0 pathways found (likely database/ID issue)
- All pathways significant (too lenient threshold)
- Completely unrelated pathways (poor data quality)

---

## Citing & Contributing

### How do I cite MetaboGraph?

**Software:**
```bibtex
@software{metabograph2024,
  author = {Your Name},
  title = {MetaboGraph: A Comprehensive Tool for Metabolomics Data Analysis},
  year = {2024},
  url = {https://github.com/yourusername/MetaboGraph},
  doi = {10.5281/zenodo.XXXXXXX}
}
```

**Manuscript:**
> Your Name et al. (2024). MetaboGraph: A Comprehensive Tool for Metabolomics Data Analysis. *Journal Name*, Volume(Issue), pages. DOI: 10.XXXX/XXXXXX

**Also cite databases used:**
- HMDB, KEGG, LipidMaps, PathBank, SMPDB, WikiPathways
- See [Database Setup Guide](user-guides/04-database-setup.md) for citations

### Can I contribute to MetaboGraph?

Yes! Contributions welcome:
- **Bug reports:** [GitHub Issues](https://github.com/yourusername/MetaboGraph/issues)
- **Feature requests:** Submit enhancement proposals
- **Code contributions:** Fork, improve, submit pull request
- **Documentation:** Help improve guides and examples
- **Testing:** Test on different systems, report results

See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines.

### Can I use MetaboGraph for commercial purposes?

Yes! MIT License allows commercial use. Requirements:
- Include original copyright notice
- Include MIT License text

No obligation to share modifications (but encouraged!).

### Can I modify the code?

Yes! You can:
- Modify for your needs
- Distribute modifications
- Use in proprietary software

Encouraged to contribute improvements back to community (optional).

---

## Additional Resources

### Where can I learn more?

- **User Guides:** [docs/user-guides/](user-guides/) - Detailed tutorials for each feature
- **Installation:** [docs/installation.md](installation.md) - Setup instructions
- **Troubleshooting:** [docs/troubleshooting.md](troubleshooting.md) - Fix common issues
- **GitHub:** [Issues, discussions, updates](https://github.com/yourusername/MetaboGraph)

### Example workflows?

See `examples/` folder for:
- Sample data files
- Step-by-step workflows
- Expected outputs

### Is there a video tutorial?

Coming soon! Check GitHub repository for updates.

### Can I get one-on-one help?

For research collaborations or consulting:
- Email: support@metabograph.org
- Mention your institution and project

For general support, use GitHub Issues (helps entire community).

---

## Still have questions?

**Ask the community:**
- GitHub Discussions: [github.com/yourusername/MetaboGraph/discussions](https://github.com/yourusername/MetaboGraph/discussions)
- GitHub Issues: [github.com/yourusername/MetaboGraph/issues](https://github.com/yourusername/MetaboGraph/issues)

**Contact directly:**
- Email: support@metabograph.org
- Include: MetaboGraph version, OS, detailed question

**Check documentation:**
- User guides, troubleshooting, and installation guides cover most topics

---

**Last Updated:** December 2024  
**Version:** 1.0.1
