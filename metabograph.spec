# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller Spec File for MetaboGraph
======================================
This spec file builds a standalone executable for MetaboGraph - Metabolite Annotation & Pathway Analysis Tool.

To build:
    pyinstaller metabograph.spec

Output:
    dist/MetaboGraph/MetaboGraph.exe
"""

block_cipher = None

import sys
import os
from PyInstaller.utils.hooks import collect_data_files, collect_submodules, collect_dynamic_libs

# Get the current directory
spec_root = os.path.abspath(SPECPATH)

# Collect all data files from packages that need them
datas = []
datas += collect_data_files('sklearn')
datas += collect_data_files('scipy')
datas += collect_data_files('matplotlib')
datas += collect_data_files('seaborn')
datas += collect_data_files('pandas')
datas += collect_data_files('plotly')
datas += collect_data_files('statsmodels')
datas += collect_data_files('networkx')
datas += collect_data_files('kaleido')
datas += collect_data_files('pyarrow')
# Dash and dash_cytoscape need their package files
datas += collect_data_files('dash')
datas += collect_data_files('dash_cytoscape')
datas += collect_data_files('dash_core_components')
datas += collect_data_files('dash_html_components')
datas += collect_data_files('dash_table')

# Bundle critical JSON/config assets only (if they exist)
# NOTE: Database .feather files are NOT bundled - users must provide them in the working directory
import os

# Add config files if they exist
config_files = [
    ('main_script/universe_config.json', '.'),
    ('main_script/factor_mapping_config.json', '.'),
    ('main_script/id_annotation_prefs.json', '.'),
]

for src, dest in config_files:
    full_path = os.path.join(spec_root, src)
    if os.path.exists(full_path):
        datas.append((src, dest))
    else:
        print(f"Note: {src} not found - will be generated at runtime if needed")

# Add lipid class annotation if it exists
lipid_annotation_path = os.path.join(spec_root, 'Databases', 'lipid_class_annotation.txt')
if os.path.exists(lipid_annotation_path):
    datas.append(('Databases/lipid_class_annotation.txt', 'Databases'))

# Collect hidden imports - all modules used by the application
hiddenimports = [
    # Core Python libraries
    'tkinter',
    'tkinter.ttk',
    'tkinter.filedialog',
    'tkinter.messagebox',
    'tkinter.scrolledtext',
    'tkinter.simpledialog',
    'abc',
    'argparse',
    'io',
    
    # Data processing
    'pandas',
    'pandas.io.excel._openpyxl',
    'pandas.io.excel._xlsxwriter',
    'pandas.io.formats.style',
    'pandas.io.json._json',
    'pandas.core.computation.expressions',
    'numpy',
    'numpy.core._methods',
    'numpy.lib.format',
    'numpy.fft',
    'numpy.random',
    'numpy.linalg',
    
    # Datetime
    'datetime',
    'time',
    
    # Excel support
    'openpyxl',
    'openpyxl.cell._writer',
    'openpyxl.styles.stylesheet',
    'xlsxwriter',
    'xlrd',
    'pyarrow',
    'pyarrow.feather',
    'pyarrow.parquet',
    
    # Scientific computing
    'scipy',
    'scipy._lib',
    'scipy._lib.messagestream',
    'scipy._cyutility',
    'scipy.stats',
    'scipy.stats._stats_py',
    'scipy.stats.distributions',
    'scipy.stats.contingency',
    'scipy.spatial',
    'scipy.spatial.distance',
    'scipy.special',
    'scipy.special.cython_special',
    'scipy.integrate',
    'scipy.sparse',
    'scipy.sparse._sparsetools',
    'scipy.sparse.csgraph',
    'scipy.sparse.csgraph._tools',
    
    # Machine learning
    'sklearn',
    'sklearn.base',
    'sklearn._cyutility',
    'sklearn.decomposition',
    'sklearn.decomposition._pca',
    'sklearn.preprocessing',
    'sklearn.preprocessing._data',
    'sklearn.metrics',
    'sklearn.metrics._classification',
    'sklearn.metrics._ranking',
    'sklearn.metrics.pairwise',
    'sklearn.cluster',
    'sklearn.cluster._agglomerative',
    'sklearn.utils',
    'sklearn.utils._typedefs',
    'sklearn.utils._heap',
    'sklearn.utils._sorting',
    'sklearn.utils._cython_blas',
    'sklearn.utils.validation',
    'sklearn.utils._isfinite',
    'sklearn.utils._chunking',
    'sklearn.utils._param_validation',
    'sklearn.neighbors',
    'sklearn.neighbors._partition_nodes',
    'sklearn.externals',
    'sklearn.externals.array_api_compat',
    'sklearn.externals.array_api_compat.numpy',
    'sklearn.externals.array_api_compat.numpy.fft',
    'sklearn.tree',
    'sklearn.ensemble',
    'sklearn.linear_model',
    
    # Statistics
    'statsmodels',
    'statsmodels.stats',
    'statsmodels.stats.multicomp',
    'statsmodels.stats.multitest',
    'statsmodels.stats.libqsturng',
    'statsmodels.stats.libqsturng.qsturng_',
    
    # Visualization
    'matplotlib',
    'matplotlib.pyplot',
    'matplotlib.backends',
    'matplotlib.backends.backend_tkagg',
    'matplotlib.backends.backend_agg',
    'matplotlib.figure',
    'matplotlib.colors',
    'matplotlib.patches',
    'mpl_toolkits',
    'mpl_toolkits.mplot3d',
    'mpl_toolkits.mplot3d.axes3d',
    'mpl_toolkits.mplot3d.art3d',
    'seaborn',
    'seaborn.matrix',
    'seaborn.categorical',
    'seaborn.distributions',
    'seaborn.regression',
    
    # Plotly (optional interactive plots)
    'plotly',
    'plotly.express',
    'plotly.graph_objects',
    'plotly.graph_objs',
    'plotly.subplots',
    'plotly.io',
    'plotly.offline',
    
    # Dash and Cytoscape for interactive networks
    'dash',
    'dash.dependencies',
    'dash.development',
    'dash_cytoscape',
    'dash_core_components',
    'dash_html_components',
    'dash_table',
    'flask',
    'flask_compress',
    'werkzeug',
    
    # Network analysis
    'networkx',
    'networkx.algorithms',
    'networkx.algorithms.community',
    'networkx.drawing',
    'networkx.drawing.layout',
    'networkx.drawing.spring_layout',
    'networkx.generators',
    'networkx.utils',
    'networkx.classes',
    
    # Image/Plot rendering for network export
    'kaleido',
    'kaleido.scopes',
    'kaleido.scopes.plotly',
    
    # HTTP requests
    'requests',
    'requests.packages',
    'requests.packages.urllib3',
    'urllib3',
    'urllib',
    'urllib.request',
    'urllib.error',
    'urllib.parse',
    'certifi',
    'charset_normalizer',
    
    # Web and networking
    'socket',
    'webbrowser',
    
    # XML parsing
    'xml',
    'xml.etree',
    'xml.etree.ElementTree',
    'lxml',
    'lxml.etree',
    
    # Crypto (if using license generation)
    'cryptography',
    'cryptography.fernet',
    
    # Utilities
    'pickle',
    'json',
    'logging',
    'logging.config',
    'logging.handlers',
    're',
    'threading',
    'concurrent.futures',
    'pathlib',
    'tempfile',
    'traceback',
    'dataclasses',
    'typing',
    'typing_extensions',
    'collections',
    'collections.abc',
    'collections.defaultdict',
    'itertools',
    'functools',
    'operator',
    'copy',
    'warnings',
    'weakref',
    'queue',
    'enum',
    'uuid',
    'glob',
    'shutil',
    'subprocess',
    'platform',
    
    # Encoding
    'encodings',
    'encodings.utf_8',
    'encodings.cp1252',
    'encodings.ascii',
    
    # Compression (needed by PyArrow)
    'zlib',
    'gzip',
    'bz2',
    'lzma',
    
    # Math
    'math',
    'cmath',
    'decimal',
    'fractions',
    
    # Platform specific
    'ctypes',
    'ctypes.wintypes',
    
    # Image processing
    'PIL',
    'PIL.Image',
    
    # User custom modules (ensure they're included)
    'gui',
    'gui.main',
    'gui.tabs',
    'gui.tabs.pathway_analysis_parent_tab',
    'gui.tabs.pathway_network_tab',
    'gui.tabs.pathway_annotation_tab',
    'gui.tabs.id_annotation_tab',
    'gui.tabs.data_cleaning_tab',
    'gui.tabs.database_setup_tab',
    'gui.tabs.help_tab',
    'gui.tabs.multiomics_analysis_tab',
    'gui.tabs.comparative_analysis_tab',
    'gui.shared',
    'gui.shared.base_tab',
    'gui.shared.data_manager',
    'gui.shared.utils',
    'gui.shared.column_assignment',
    'gui.shared.pairwise_column_mapper',
    'main_script',
    'main_script.metabolite_ID_annotator',
    'main_script.metabolite_pathways_annotator',
    'main_script.metabolite_pathway_network',
    'main_script.metabolite_data_cleaner',
    'main_script.fisher_ora_pathway_analysis',
    'main_script.comparative_analysis',
    'main_script.ml_pathway_models',
    'main_script.lipid_search',
    'main_script.database_builder',
    'main_script.calculate_universe',
    'main_script.factor_mapping_manager',
    'main_script.interactive_network_viewer',
    'main_script.help',
]

# Collect all submodules for key packages - this ensures all Cython extensions are included
hiddenimports += collect_submodules('sklearn')  # Critical: includes all sklearn Cython extensions
hiddenimports += collect_submodules('scipy.stats')
hiddenimports += collect_submodules('scipy.special')
hiddenimports += collect_submodules('scipy.sparse')
hiddenimports += collect_submodules('statsmodels.stats')
hiddenimports += collect_submodules('matplotlib.backends')
hiddenimports += collect_submodules('networkx')
hiddenimports += collect_submodules('plotly')
hiddenimports += collect_submodules('pandas.io')
hiddenimports += collect_submodules('numpy.core')
hiddenimports += collect_submodules('numpy.lib')

# Add user custom analysis modules
hiddenimports += [
    'main_script.metabolite_pathway_network',
    'main_script.metabolite_ID_annotator',
    'main_script.metabolite_data_cleaner',
    'main_script.fisher_ora_pathway_analysis',
    'main_script.factor_mapping_manager',
    'main_script.metabolite_pathways_annotator',
    'main_script.database_builder',
    # 'main_script.calculate_universe',  # Commented out - universe now calculated from user's metabolite list
    'main_script.lipid_search',
    'main_script.interactive_network_viewer',
    'main_script.help',
]

# Collect binaries from packages with compiled extensions
from PyInstaller.utils.hooks import collect_dynamic_libs
binaries = []
binaries += collect_dynamic_libs('sklearn')
binaries += collect_dynamic_libs('scipy')
binaries += collect_dynamic_libs('numpy')

# Analysis of main script and all dependencies
a = Analysis(
    ['metabograph.py'],
    pathex=[spec_root],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=['main_script/runtime_license_hook.py'],
    excludes=[
        'IPython',
        'jupyter',
        'notebook',
        'PyQt5',
        'PyQt6',
        'PySide2',
        'PySide6',
        'wx',
        'pytest',
        'sphinx',
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

# Collect PYZ archive
pyz = PYZ(
    a.pure,
    a.zipped_data,
    cipher=block_cipher
)

# Create executable
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='MetaboliteAnnotationTool',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,  # DISABLED: UPX can corrupt scipy/numpy DLLs causing crashes on other machines
    console=True,  # Set to False for GUI app (no console window)
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='logo.ico',
    version=None,  # Add version info file if desired
)

# Collect all files into distribution directory
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,  # DISABLED: UPX can corrupt scipy/numpy DLLs causing crashes on other machines
    upx_exclude=[],
    name='MetaboGraph'
)

# Create empty Databases folder in the output directory
import os
databases_dir = os.path.join(DISTPATH, 'MetaboGraph', 'Databases')
os.makedirs(databases_dir, exist_ok=True)
print(f"Created Databases directory: {databases_dir}")
print("Note: Users must copy database .feather files to this folder after installation.")
