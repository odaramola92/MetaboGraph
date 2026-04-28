# Contributing to MetaboGraph

First off, thank you for considering contributing to MetaboGraph! It's people like you that make MetaboGraph such a great tool for the metabolomics community.

## Code of Conduct

This project and everyone participating in it is governed by our [Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code.

## How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check the [existing issues](https://github.com/yourusername/MetaboGraph/issues) to avoid duplicates.

**When reporting a bug, please include:**
- **Clear descriptive title**
- **Detailed steps to reproduce** the bug
- **Expected behavior** vs **actual behavior**
- **Screenshots** if applicable
- **System information**:
  - OS (Windows 10, macOS 12, Ubuntu 22.04, etc.)
  - Python version (if using Python installation)
  - MetaboGraph version
- **Error messages** or log files
- **Sample data** (if possible and not confidential)

**Bug Report Template:**
```markdown
### Description
[Clear description of the bug]

### Steps to Reproduce
1. Go to '...'
2. Click on '...'
3. Enter '...'
4. See error

### Expected Behavior
[What you expected to happen]

### Actual Behavior
[What actually happened]

### Environment
- OS: Windows 11
- Python: 3.13 (or 3.8+)
- MetaboGraph: v1.0.0

### Additional Context
[Screenshots, error logs, etc.]
```

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues.

**When suggesting an enhancement, please include:**
- **Clear title** and detailed description
- **Use case**: Why is this enhancement useful?
- **Examples**: How would it work?
- **Alternatives**: Other ways to achieve the same goal?
- **Mockups**: If it's a UI enhancement

### Pull Requests

**Process:**
1. Fork the repo and create your branch from `main`
2. Make your changes
3. Add tests if applicable
4. Update documentation
5. Ensure tests pass
6. Submit pull request

**Pull Request Guidelines:**
- Follow the existing code style
- Write clear commit messages
- Reference related issues
- Update documentation
- Add tests for new features
- Keep changes focused (one feature per PR)

## Development Setup

### Prerequisites
- Python 3.8 or higher (developed on Python 3.13)
- Git
- Virtual environment tool (venv, conda)

### Setup Steps
```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/MetaboGraph.git
cd MetaboGraph

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
pip install -r requirements-dev.txt  # Development dependencies

# Install in editable mode
pip install -e .

# Run tests
pytest tests/
```

### Project Structure
```
MetaboGraph/
├── metabograph.py          # Main entry point
├── gui/                    # GUI modules
│   ├── main.py             # Main application window
│   ├── tabs/               # Tab implementations
│   │   ├── data_cleaning_tab.py
│   │   ├── id_annotation_tab.py
│   │   ├── pathway_annotation_tab.py
│   │   ├── pathway_network_tab.py
│   │   ├── comparative_analysis_tab.py
│   │   ├── multiomics_analysis_tab.py
│   │   ├── database_setup_tab.py
│   │   └── help_tab.py
│   └── shared/             # Shared utilities
│       ├── data_manager.py
│       ├── utils.py
│       ├── base_tab.py
│       └── column_assignment.py
├── main_script/            # Core analysis modules
│   ├── metabolite_data_cleaner.py
│   ├── metabolite_ID_annotator.py
│   ├── metabolite_pathways_annotator.py
│   ├── metabolite_pathway_network.py
│   ├── fisher_ora_pathway_analysis.py
│   ├── comparative_analysis.py
│   └── database_builder.py
├── docs/                   # Documentation
│   └── user-guides/        # Per-tab user guides
├── Databases/              # Database files (.feather)
├── requirements.txt        # Python dependencies
├── setup.py                # Package setup
└── metabograph.spec        # PyInstaller build spec
```

## Coding Standards

### Python Style Guide
- Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/)
- Use meaningful variable and function names
- Maximum line length: 120 characters
- Use type hints where appropriate
- Document all public functions and classes

### Documentation Standards
```python
def annotate_metabolites(data: pd.DataFrame, 
                         mz_tolerance: float = 0.005,
                         confidence_threshold: str = "medium") -> pd.DataFrame:
    """
    Annotate metabolites using multiple databases.
    
    Args:
        data: DataFrame with metabolite features
        mz_tolerance: m/z matching tolerance in Daltons
        confidence_threshold: Minimum confidence level ("high", "medium", "low")
    
    Returns:
        DataFrame with added annotation columns
    
    Raises:
        ValueError: If data is missing required columns
        FileNotFoundError: If databases are not found
    
    Example:
        >>> data = pd.read_excel("metabolites.xlsx")
        >>> annotated = annotate_metabolites(data, mz_tolerance=0.01)
    """
    # Implementation
    pass
```

### Testing Standards
- Write tests for new features
- Maintain or improve code coverage
- Use pytest for testing
- Follow Arrange-Act-Assert pattern

```python
def test_remove_duplicates():
    # Arrange
    data = create_test_data_with_duplicates()
    
    # Act
    result = remove_duplicates(data, mz_tolerance=0.005, rt_tolerance=0.1)
    
    # Assert
    assert len(result) < len(data)
    assert result['Name'].is_unique
```

### Commit Message Guidelines
Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <subject>

<body>

<footer>
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, no logic change)
- `refactor`: Code refactoring
- `test`: Adding or updating tests
- `chore`: Maintenance tasks

**Examples:**
```
feat(annotation): add ChEBI database integration

- Implement ChEBI API client
- Add ChEBI ID to annotation output
- Update tests for ChEBI support

Closes #123
```

```
fix(cleaning): handle negative RT values correctly

Previously, negative retention times caused crashes.
Now they are flagged as invalid and user is notified.

Fixes #456
```

## Adding New Features

### New Database Integration
1. **Create database module** in `main_script/`
2. **Implement functions**:
   - `download_database()`: Download raw data
   - `process_database()`: Convert to .feather format
   - `load_database()`: Read into memory
   - `search_database()`: Query metabolites
3. **Add to database_builder.py**
4. **Update GUI** in `database_setup_tab.py`
5. **Add tests** in `tests/test_databases.py`
6. **Update documentation** in `docs/user-guides/04-database-setup.md`

### New Analysis Method
1. **Create module** in `main_script/`
2. **Implement core logic** with proper error handling
3. **Add GUI interface** in appropriate tab
4. **Write tests** with sample data
5. **Document usage** in user guide
6. **Add example workflow** to README

### New Visualization
1. **Choose library** (matplotlib, plotly, seaborn)
2. **Implement plotting function** with customization options
3. **Add export options** (PNG, SVG, PDF)
4. **Make interactive** if using plotly
5. **Add to appropriate tab**
6. **Document parameters** and usage

## Testing

### Running Tests
```bash
# Run all tests
pytest tests/

# Run specific test file
pytest tests/test_annotation.py

# Run with coverage
pytest --cov=main_script --cov=gui tests/

# Run with verbose output
pytest -v tests/

# Run and stop at first failure
pytest -x tests/
```

### Writing Tests
```python
# tests/test_data_cleaner.py
import pytest
import pandas as pd
from main_script.metabolite_data_cleaner import DataCleaner

@pytest.fixture
def sample_data():
    """Fixture providing sample metabolite data"""
    return pd.DataFrame({
        'Name': ['Glucose', 'Lactate', 'Pyruvate'],
        'mz': [180.0634, 90.0320, 88.0164],
        'RT': [5.2, 3.1, 2.8],
        'Intensity': [10000, 5000, 3000]
    })

def test_remove_duplicates(sample_data):
    """Test duplicate removal functionality"""
    # Add a duplicate
    duplicate = sample_data.iloc[0:1].copy()
    data_with_dup = pd.concat([sample_data, duplicate], ignore_index=True)
    
    cleaner = DataCleaner()
    result = cleaner.remove_duplicates(data_with_dup, mz_tol=0.005, rt_tol=0.1)
    
    assert len(result) == len(sample_data)
    assert result['Name'].is_unique

def test_handle_missing_values(sample_data):
    """Test missing value imputation"""
    # Introduce missing value
    data = sample_data.copy()
    data.loc[0, 'Intensity'] = None
    
    cleaner = DataCleaner()
    result = cleaner.fill_missing(data, method='mean')
    
    assert result['Intensity'].isna().sum() == 0
```

## Documentation

### Updating Documentation
- **User Guides**: Add/update docs in `docs/user-guides/`
- **API Docs**: Update docstrings in code
- **README**: Update main README.md for new features
- **Changelog**: Add entry to CHANGELOG.md

### Documentation Style
- Use clear, concise language
- Include code examples
- Add screenshots for UI features
- Link to related documentation
- Use consistent formatting (Markdown)

## Release Process

### Version Numbering
We use [Semantic Versioning](https://semver.org/):
- **MAJOR**: Incompatible API changes
- **MINOR**: New features (backward compatible)
- **PATCH**: Bug fixes (backward compatible)

### Release Checklist
- [ ] Update version in `setup.py` and `metabograph.py`
- [ ] Update CHANGELOG.md
- [ ] Run all tests
- [ ] Update documentation
- [ ] Build executable
- [ ] Test executable on clean machine
- [ ] Create GitHub release
- [ ] Update website (if applicable)

## Questions?

- **General Questions**: [GitHub Discussions](https://github.com/yourusername/MetaboGraph/discussions)
- **Bug Reports**: [GitHub Issues](https://github.com/yourusername/MetaboGraph/issues)
- **Security Issues**: Email security@metabograph.org (do not create public issue)
- **Other**: Email support@metabograph.org

## Recognition

Contributors will be:
- Listed in [CONTRIBUTORS.md](CONTRIBUTORS.md)
- Acknowledged in release notes
- Mentioned in documentation (for major contributions)

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to MetaboGraph! 🎉
