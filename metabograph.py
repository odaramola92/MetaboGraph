#!/usr/bin/env python3
"""
MetaboGraph - Quick Launch Script

This is a convenience script to launch MetaboGraph from the root directory.
Simply run: python metabograph.py
"""

import sys
import os

# Add the directory containing this script to Python path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
os.chdir(script_dir)

# Import and run the main GUI
if __name__ == "__main__":
    try:
        from gui.main import main
        main()
    except Exception as e:
        print(f"Error launching MetaboGraph: {e}")
        import traceback
        traceback.print_exc()
        input("\nPress Enter to exit...")
        sys.exit(1)
