# MetaboGraph

**MetaboGraph** is a comprehensive metabolomics and lipidomics data analysis platform that provides end-to-end workflow for metabolite identification, pathway enrichment, and multi-omics integration. Designed for both beginner and advanced users, MetaboGraph combines powerful computational methods with an intuitive graphical interface.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.1+](https://img.shields.io/badge/python-3.1+-blue.svg)](https://www.python.org/downloads/)
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

## 📸 Screenshots

see 'docs/images'

## 🚀 Quick Start

### Installation

#### Option 1: Standalone Executable (Windows)
1. Download the latest release from [Releases](link-to-releases)
2. Extract the ZIP file
3. Run `MetaboliteAnnotationTool.exe`
4. Set up databases (see [Database Setup Guide](docs/user-guides/04-database-setup.md))

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

### System Requirements
- **OS**: Windows 10/11, macOS 10.14+, or Linux
- **Python**: 3.1 or higher (developed on 3.13)
- **RAM**: 8 GB minimum, 16 GB recommended
- **Disk Space**: 5 GB (includes databases)
- **Internet**: Required for database downloads and online API queries

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
- [Installation Guide](docs/installation.md)
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
├── metabograph.py           # Main entry point
├── metabograph.spec         # PyInstaller spec file
├── requirements.txt         # Python dependencies
├── setup.py                 # Package setup
├── gui/                     # GUI modules
│   ├── main.py              # Main application window
│   ├── shared/              # Shared utilities
│   │   ├── data_manager.py
│   │   ├── utils.py
│   │   └── column_assignment.py
│   └── tabs/                # Tab implementations
│       ├── data_cleaning_tab.py
│       ├── id_annotation_tab.py
│       ├── pathway_analysis_parent_tab.py
│       ├── database_setup_tab.py
│       └── ...
├── main_script/             # Core analysis modules
│   ├── metabolite_data_cleaner.py
│   ├── metabolite_ID_annotator.py
│   ├── metabolite_pathways_annotator.py
│   ├── fisher_ora_pathway_analysis.py
│   ├── database_builder.py
│   └── ...
├── Databases/               # Database files (user-provided)
├── docs/                    # Documentation
│   ├── user-guides/
│   └── images/
└── build/                   # Build artifacts (generated)
```

### Building from Source

#### Create Standalone Executable
```bash
# Install PyInstaller
pip install pyinstaller

# Build executable
pyinstaller metabograph.spec

# Output: dist/MetaboGraph/MetaboliteAnnotationTool.exe
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

If you use MetaboGraph in your research, please cite:

```bibtex
@software{metabograph2024,
  author = {Your Name},
  title = {MetaboGraph: A Comprehensive Platform for Metabolomics Data Analysis},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/odaramola92/MetaboGraph}
}
```

Also see [CITATION.cff](CITATION.cff) for structured citation information.

## 🙏 Acknowledgments

- **Database Providers**: HMDB, LipidMaps, PathBank, SMPDB, WikiPathways
- **Libraries**: Pandas, NumPy, SciPy, Scikit-learn, NetworkX, Plotly, Dash
- **Contributors**: See [CONTRIBUTORS.md](CONTRIBUTORS.md)

## 📞 Support

- **Documentation**: [docs/](docs/)
- **Issues**: [GitHub Issues](https://github.com/odaramola92/MetaboGraph/issues)
- **Discussions**: [GitHub Discussions](https://github.com/odaramola92/MetaboGraph/discussions)
- **Email**: support@metabograph.org

## 🗺️ Roadmap

### Current Version (v1.0)
- ✅ Data cleaning and preprocessing
- ✅ Multi-database ID annotation
- ✅ Pathway enrichment analysis
- ✅ Comparative analysis
- ✅ Multi-omics integration
- ✅ Database management

## ⚠️ Disclaimer

MetaboGraph is provided "as is" without warranty of any kind. The database content and annotations are derived from third-party sources and may contain errors. Always validate critical findings experimentally.

## 📈 Statistics

![GitHub stars](https://img.shields.io/github/stars/odaramola92/MetaboGraph?style=social)
![GitHub forks](https://img.shields.io/github/forks/odaramola92/MetaboGraph?style=social)
![GitHub issues](https://img.shields.io/github/issues/odaramola92/MetaboGraph)
![GitHub pull requests](https://img.shields.io/github/issues-pr/odaramola92/MetaboGraph)

---

**Made with ❤️ for the metabolomics community**

[Documentation](docs/)
