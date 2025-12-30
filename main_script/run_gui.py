#!/usr/bin/env python3
"""
Entry point for the MetaboGraph GUI application.

This script:
1. Properly configures Python paths
2. Checks for required database files
3. Launches the MetaboGraph GUI
4. Handles cleanup on exit
"""
import sys
import os

# Determine the correct base directory for PyInstaller or normal execution
if getattr(sys, 'frozen', False):
    # Running as compiled PyInstaller executable
    # Use the directory where the exe is located
    application_path = os.path.dirname(sys.executable)
else:
    # Running as normal Python script
    application_path = os.path.dirname(os.path.abspath(__file__))

# Change working directory to the application path
# This ensures relative paths work correctly when double-clicking the exe
os.chdir(application_path)

# Add the script directory to Python path
sys.path.insert(0, application_path)

# Check for required database files in working directory or Databases subfolder
required_databases = [
    'hmdb_database.feather',
    'lipidmap.feather',
    'merged_SMP_metabolites.feather',
    'pathbank_selected.feather',
    'wikipathways_homo_sapiens.feather',
    'wikipathways_mus_musculus.feather',
    'wikipathways_rattus_norvegicus.feather',
]

# Debug: Show what directory we're checking
print(f"Application path: {application_path}")
print(f"Current working directory: {os.getcwd()}")
print(f"Checking for database files...")

# Check multiple possible locations
possible_paths = [
    os.getcwd(),  # Working directory
    os.path.join(os.getcwd(), "Databases"),  # Databases subfolder
]

missing_databases = []
found_databases = []
database_location = None

for db in required_databases:
    found = False
    for search_path in possible_paths:
        db_path = os.path.join(search_path, db)
        if os.path.exists(db_path):
            found_databases.append(db)
            database_location = search_path
            print(f"  ✓ Found: {db}")
            found = True
            break
    
    if not found:
        missing_databases.append(db)
        print(f"  ✗ Missing: {db}")

if missing_databases:
    warning_msg = f"""
╔════════════════════════════════════════════════════════════════╗
║         DATABASE FILES NOT FOUND - SETUP REQUIRED              ║
╚════════════════════════════════════════════════════════════════╝

Working Directory: {os.getcwd()}
Database Location: {database_location if database_location else 'Not found'}

Found {len(found_databases)} of {len(required_databases)} database files.

Missing files:
"""
    for db in missing_databases:
        warning_msg += f"  ❌ {db}\n"
    
    warning_msg += """
INSTRUCTIONS:
1. The application will start, but some features will not work until databases are set up
2. Go to the 'Database' tab in the application
3. Follow the instructions to download and generate the required database files
4. Once databases are ready, features will be available

The application is starting...
"""
    print(warning_msg)
else:
    print(f"\n✓ All {len(found_databases)} database files found!\n")

# Now we can import and run the GUI
if __name__ == "__main__":
    try:
        from gui.main import main
        main()
    except Exception as e:
        print(f"Error launching GUI: {e}")
        import traceback
        traceback.print_exc()
        # Keep console open so user can see the error
        input("\nPress Enter to exit...")
        sys.exit(1)
