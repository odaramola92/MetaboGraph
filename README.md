# MetaboGraph

**MetaboGraph** is a Framework for Metabolomics and Lipidomics Annotation and Pathway Network Analysis that provides end-to-end workflow for metabolite identification, pathway enrichment, net analysis, and multi-omics integration. Designed for both beginner and advanced users, MetaboGraph combines powerful computational methods with an intuitive graphical interface.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)]()

## 🌟 Features

### 📊 Data Cleaning & Preprocessing
- Import data from Excel formats (Compound Discoverer, LipidSearch exports)
- Separate processing for metabolomics and lipidomics data
- Handle positive and negative ionization modes (combined or separate files)
- Reference ion filtering with customizable ion selection
- PPM and RT threshold-based duplicate detection
- Column name cleaning with customizable patterns
- MZCloud file integration for enriched annotations
- MS2 filtering options (DDA preferred ion)

### 🔬 ID Annotation
- **Multi-database annotation**: HMDB, PubChem, KEGG, LipidMaps, ChEBI, CAS
- **API-first approach**: PubChem/KEGG APIs with offline HMDB fallback
- **Cross-reference retrieval**: Automatic ID mapping across databases
- **Standard and Custom modes**: Flexible file structure handling
- **Lipid-specific annotation** with LipidMaps integration
- Re-filter existing results without re-running annotation
- Parallel processing with configurable workers

### 🧬 Pathway Analysis
- Map metabolites to biological pathways
- **Databases**: SMPDB, WikiPathways, PathBank, HMDB, KEGG, REACTOME
- **Fisher's Exact Test**: Over-representation analysis with FDR correction
- Interactive pathway network viewer (Dash-based, MetaboAnalyst-like)
- Cytoscape export for publication-quality networks
- Project-based output with auto-naming
- Upstream regulators and disease associations
- Organism-specific (Human, Mouse, Rat)

### 🔀 Comparative Analysis
- Compare pathway enrichment results across multiple datasets
- Z-score heatmaps with clustering (blue-gray-red scale)
- P-value heatmaps
- Pathway overlap bar charts
- Venn diagrams for pathway comparison (2-5 datasets)
- Filter by Z-score and P-value thresholds
- Export visualizations as PNG/PDF

### 🧩 Multi-Omics Integration
- Merge multiple pathway-annotated datasets (metabolite + lipid)
- Column validation (Feature ID, P-value, Log2FC)
- Smart merging with consolidation rules
- Automatic push to Network Tab for unified analysis
- Combine positive/negative mode results
- Source tracking for data origin

### 🗄️ Database Management
- Download and process metabolite databases
- HMDB, LipidMaps, PathBank, SMPDB, WikiPathways
- Automatic format conversion to efficient `.feather` format
- Database validation and update checking
- Local storage for offline analysis

## 🚀 Quick Start

### Installation

#### Option 1: Standalone Executable (Windows)
1. Download the latest release from [Releases](https://github.com/odaramola92/MetaboGraph/releases)
2. Extract the ZIP file
3. Run the executable (named `MetaboGraph_v<version>.exe` or similar – the version matches the release tag)
4. Set up databases (see [Database Setup Guide](docs/user-guides/04-database-setup.md))

**See [INSTALLATION.md](INSTALLATION.md) for detailed build and deployment instructions.**

#### Option 2: Python Installation
```bash
# Clone the repository
git clone https://github.com/odaramola92/MetaboGraph.git
cd MetaboGraph

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run the application
python metabograph.py
```

For comprehensive build and deployment instructions, see [INSTALLATION.md](INSTALLATION.md).

> **Note for non‑Python users:** you do **not** need to install Python to use MetaboGraph. Simply download the
> latest release (see the “Standalone Executable” option above), extract the ZIP and run the included executable (it will be named `MetaboGraph_v<version>.exe`). If you build the executable yourself with PyInstaller, you
> can ZIP the `MetaboGraph` folder and host it (e.g. on Massive or any file server) so that
> other users can download and unzip without any Python installation. This is the recommended
> distribution method for colleagues who are unfamiliar with Python or who cannot install it.

### System Requirements
- **OS**: Windows 10/11, macOS 10.14+, or Linux
- **Python**: 3.8 or higher (tested on 3.8-3.13)
- **RAM**: 8 GB minimum, 16 GB recommended
- **Disk Space**: 5 GB (includes databases)
- **Internet**: Required for database downloads and online API queries

## 🟢 Application Status

**Operational:** All 5 main tabs load successfully without errors. See [STARTUP_STATUS.md](STARTUP_STATUS.md) for verification details.

### First-Time Setup

1. **Download Databases**
   - Go to Database Setup tab
   - Download HMDB (required)
   - Download other databases as needed
   - Process each database (converts to efficient format)

2. **Verify Installation**
   - Check Database Status shows all databases found
   - Load sample data in Data Cleaning tab
   - Test annotation with a few metabolites

3. **Start Your Analysis**
   - Follow the [User Guides](docs/user-guides/) for each tab

## 📖 Documentation

### User Guides
- [Data Cleaning Tab](docs/user-guides/01-data-cleaning.md)
- [ID Annotation Tab](docs/user-guides/02-id-annotation.md)
- [Pathway Analysis Tab](docs/user-guides/03-pathway-analysis.md)
- [Database Setup Tab](docs/user-guides/04-database-setup.md)
- [Comparative Analysis Tab](docs/user-guides/05-comparative-analysis.md)
- [Multi-Omics Integration Tab](docs/user-guides/06-multiomics-integration.md)

### Additional Documentation
- [Installation & Build Guide](INSTALLATION.md)
- [Diagnostics Guide](DIAGNOSTICS.md)
- [Startup Status Report](STARTUP_STATUS.md)
- [Pre-Publication Checklist](PRE_PUBLICATION_CHECKLIST.md)
- [Troubleshooting](docs/troubleshooting.md)
- [FAQ](docs/faq.md)
- [API Documentation](docs/api-documentation.md)
- [Contributing Guidelines](CONTRIBUTING.md)

## 💡 Example Workflow

```
1. Data Cleaning
   └── Load raw LC-MS data
   └── Remove duplicates (m/z ±0.005, RT ±0.1 min)
   └── Handle missing values
   └── Apply log2 transformation
   └── Export cleaned data

2. ID Annotation
   └── Load cleaned data
   └── Annotate using HMDB + PubChem
   └── Review confidence scores
   └── Export annotated dataset

3. Pathway Analysis
   └── Load annotated metabolites
   └── Map to pathways (PathBank + WikiPathways)
   └── Run Fisher's exact test
   └── Visualize enriched pathways
   └── Generate pathway network

4. Comparative Analysis
   └── Define Control vs. Treatment groups
   └── Run t-test with FDR correction
   └── Create volcano plot
   └── Export differential metabolites

5. Multi-Omics Integration (Optional)
   └── Load metabolomics + transcriptomics
   └── Calculate correlations
   └── Build integrated network
   └── Joint pathway analysis
```

## 🔧 Development

### Project Structure
```
MetaboGraph/
├── 📄 Configuration & Setup
│   ├── metabograph.py              # Main entry point
│   ├── metabograph.spec            # PyInstaller spec file
│   ├── setup.py                    # Package setup
│   ├── requirements.txt            # Python dependencies
│   ├── run_metabograph.bat         # Windows batch launcher
│   └── run_metabograph.sh          # Unix shell launcher
├── 📋 Documentation & Metadata
│   ├── README.md                   # This file
│   ├── LICENSE                     # MIT License
│   ├── CITATION.cff                # Citation metadata
│   ├── CONTRIBUTING.md             # Contribution guidelines
│   ├── INSTALLATION.md             # Installation & build guide
│   ├── DIAGNOSTICS.md              # Diagnostics & troubleshooting
│   ├── STARTUP_STATUS.md           # Startup verification report
│   ├── UNICODE_FIX_COMPLETE.md     # Unicode fix documentation
│   ├── PRE_PUBLICATION_CHECKLIST.md # Pre-publication tasks
│   └── id_annotation_prefs.json    # ID annotation preferences
├── gui/                            # GUI modules
│   ├── __init__.py
│   ├── main.py                     # Main application window
│   ├── shared/                     # Shared utilities
│   │   ├── __init__.py
│   │   ├── data_manager.py         # Data management layer
│   │   ├── base_tab.py             # Base class for all tabs
│   │   ├── column_assignment.py    # Column mapping interface
│   │   ├── pairwise_column_mapper.py
│   │   └── utils.py                # Common utilities
│   └── tabs/                       # Tab implementations
│       ├── __init__.py
│       ├── data_cleaning_tab.py    # Data preprocessing tab
│       ├── id_annotation_tab.py    # Metabolite ID annotation
│       ├── pathway_analysis_parent_tab.py  # Pathway analysis parent
│       ├── pathway_annotation_tab.py       # Pathway annotation
│       ├── pathway_network_tab.py          # Network visualization
│       ├── comparative_analysis_tab.py     # Multi-dataset comparison
│       ├── multiomics_analysis_tab.py      # Multi-omics integration
│       ├── database_setup_tab.py           # Database management
│       ├── help_tab.py                     # Help & documentation
│       └── [obsolete files] pathway_annotation_tab_old.py
├── main_script/                    # Core analysis modules
│   ├── metabolite_data_cleaner.py  # Data cleaning engine
│   ├── metabolite_ID_annotator.py  # ID annotation engine
│   ├── metabolite_pathways_annotator.py     # Pathway mapping
│   ├── metabolite_pathway_network.py        # Network construction
│   ├── fisher_ora_pathway_analysis.py       # Pathway enrichment
│   ├── database_builder.py         # Database processing
│   ├── comparative_analysis.py     # Comparative analysis
│   ├── interactive_network_viewer.py        # Network visualization
│   ├── lipid_search.py             # Lipid-specific utilities
│   ├── ml_pathway_models.py        # Machine learning models
│   ├── calculate_universe.py       # Background universe calculation
│   ├── factor_mapping_manager.py   # Factor management
│   ├── generate_license.py         # License generation
│   ├── runtime_license_hook.py     # License hooks
│   ├── help.py                     # Help content
│   ├── run_gui.py                  # GUI launcher
│   ├── lipid_class_annotation.txt  # Lipid classification reference
│   └── universe_config.json        # Universe configuration
├── Databases/                      # Database storage
│   ├── hmdb_database.feather       # HMDB metabolite database
│   ├── lipidmap.feather            # LipidMaps database
│   ├── pathbank_selected.feather   # PathBank pathways
│   ├── All_metabolites_synonyms_hmdb.feather
│   ├── merged_SMP_metabolites.feather      # Merged SMP data
│   ├── wikipathways_homo_sapiens.feather   # WikiPathways - Human
│   ├── wikipathways_mus_musculus.feather   # WikiPathways - Mouse
│   ├── wikipathways_rattus_norvegicus.feather  # WikiPathways - Rat
│   └── lipid_class_annotation.txt  # Lipid classifications
├── docs/                           # Documentation
│   ├── faq.md                      # Frequently asked questions
│   ├── troubleshooting.md          # Troubleshooting guide
│   ├── installation.md             # Installation instructions
│   ├── images/                     # Screenshots and diagrams
│   └── user-guides/                # Tab-specific user guides
│       ├── 01-data-cleaning.md
│       ├── 02-id-annotation.md
│       ├── 03-pathway-analysis.md
│       ├── 04-database-setup.md
│       ├── 05-comparative-analysis.md
│       └── 06-multiomics-integration.md
├── build/                          # PyInstaller build artifacts
├── dist/                           # Distribution executables
├── logs/                           # Runtime logs
└── [Cache files] *.pkl, *.cache    # Database/API caches
```

### Building from Source

#### Create Standalone Executable
```bash
# Install PyInstaller
pip install pyinstaller

# Build executable
pyinstaller metabograph.spec

# Output: dist/MetaboGraph/<app_name>.exe (e.g. MetaboGraph_v1.0.0.exe)
```

### Running Tests
```bash
# Install test dependencies
pip install pytest pytest-cov

# Run tests
pytest tests/

# With coverage
pytest --cov=main_script --cov=gui tests/
```

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Ways to Contribute
- 🐛 Report bugs
- 💡 Suggest new features
- 📝 Improve documentation
- 🔧 Submit pull requests
- 🌍 Add support for new databases
- 🧪 Add test cases

### Development Setup
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests if applicable
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## 📊 Supported Databases

| Database | Purpose | Version | Size | Required |
|----------|---------|---------|------|----------|
| [HMDB](https://hmdb.ca) | Metabolite annotation | 5.0 | ~150 MB | ✅ Yes |
| [LipidMaps](https://lipidmaps.org) | Lipid annotation | 2024 | ~50 MB | For lipids |
| [PathBank](https://pathbank.org) | Pathway mapping | 2024 | ~30 MB | For pathways |
| [SMPDB](https://smpdb.ca) | Pathway enrichment | 2024 | ~20 MB | For pathways |
| [WikiPathways](https://wikipathways.org) | Curated pathways | 2024 | ~30 MB | For pathways |
| PubChem | Online queries | API | - | Optional |
| KEGG | Online queries | API | - | Optional |

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citation

If you use MetaboGraph in your research, please cite the associated publication and this:

```bibtex
@software{daramola_metabograph_2026,
  author  = {Daramola, Oluwatosin and Mechref, Yehia},
  title   = {MetaboGraph: A Framework for Metabolomics and Lipidomics Annotation and Pathway Network Analysis},
  year    = {2026},
  version = {1.0},
  url     = {https://github.com/odaramola92/MetaboGraph}
}
```

Also see [CITATION.cff](CITATION.cff) for structured citation information.

## 🙏 Acknowledgments

- **Database Providers**: HMDB, LipidMaps, PathBank, SMPDB, WikiPathways
- **Libraries**: Pandas, NumPy, SciPy, Scikit-learn, NetworkX, Plotly, Dash
- **Contributors**: See [CONTRIBUTORS.md](CONTRIBUTORS.md)

## 📞 Support & Troubleshooting

- **Documentation**: [docs/](docs/)
- **Diagnostics**: Check [DIAGNOSTICS.md](DIAGNOSTICS.md) for startup issues and log file locations
- **Issues**: [GitHub Issues](https://github.com/odaramola92/MetaboGraph/issues)
- **Discussions**: [GitHub Discussions](https://github.com/odaramola92/MetaboGraph/discussions)
- **Log Files**: Located at `%TEMP%\metabograph_*.log` for debugging
- **Email**: odaramol@ttu.edu

## ⚠️ Disclaimer

MetaboGraph is provided "as is" without warranty of any kind. The database content and annotations are derived from third-party sources and may contain errors. Always validate critical findings experimentally.

## 📈 Statistics

![GitHub stars](https://img.shields.io/github/stars/odaramola92/MetaboGraph?style=social)
![GitHub forks](https://img.shields.io/github/forks/odaramola92/MetaboGraph?style=social)
![GitHub issues](https://img.shields.io/github/issues/odaramola92/MetaboGraph)
![GitHub pull requests](https://img.shields.io/github/issues-pr/odaramola92/MetaboGraph)

---

**Made with ❤️ for the metabolomics community**
