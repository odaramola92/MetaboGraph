#!/usr/bin/env python3
"""
MetaboGraph - Quick Launch Script

This is a convenience script to launch MetaboGraph from the root directory.
Simply run: python metabograph.py
"""

import sys
import os
import tempfile
from datetime import datetime

# Add the directory containing this script to Python path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
os.chdir(script_dir)

# Setup file-based logging for diagnostics (useful for frozen executables)
log_file = None
try:
    # Try to write log to temp directory with timestamp
    temp_dir = tempfile.gettempdir()
    log_file = os.path.join(temp_dir, f"metabograph_startup_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    log_handle = open(log_file, 'w')
    
    def log_diagnostic(message: str):
        """Write diagnostic message to both console and log file"""
        timestamp = datetime.now().strftime('%H:%M:%S.%f')[:-3]
        msg = f"[{timestamp}] {message}"
        try:
            print(msg, flush=True)
        except UnicodeEncodeError:
            # Fallback for console that can't encode UTF-8
            print(msg.encode('ascii', 'ignore').decode('ascii'), flush=True)
        try:
            log_handle.write(msg + '\n')
            log_handle.flush()
        except Exception:
            pass
    
except Exception as e:
    def log_diagnostic(message: str):
        """Fallback logging if file opening fails"""
        print(f"[LOG] {message}", flush=True)


# Import and run the main GUI
if __name__ == "__main__":
    try:
        log_diagnostic("=" * 80)
        log_diagnostic("MetaboGraph Startup Sequence")
        log_diagnostic("=" * 80)
        log_diagnostic(f"Python version: {sys.version}")
        log_diagnostic(f"Working directory: {os.getcwd()}")
        log_diagnostic(f"Script directory: {script_dir}")
        
        # Test basic imports with diagnostic logging
        log_diagnostic("Testing Python standard library imports...")
        import tkinter
        log_diagnostic("[OK] tkinter imported successfully")
        
        import pandas as pd
        log_diagnostic("[OK] pandas imported successfully")
        
        import numpy as np
        log_diagnostic("[OK] numpy imported successfully")
        
        log_diagnostic("Importing MetaboGraph GUI main module...")
        from gui.main import main
        log_diagnostic("[OK] gui.main imported successfully")
        
        log_diagnostic("Launching main application...")
        main()
        log_diagnostic("[OK] Application exited normally")
        
    except Exception as e:
        import traceback
        
        error_msg = f"Error launching MetaboGraph: {e}\n\n"
        print("\n" + "=" * 80)
        print(error_msg)
        print("=" * 80)
        
        # Log full traceback
        log_diagnostic("\n" + "=" * 80)
        log_diagnostic("[ERROR] EXCEPTION ENCOUNTERED:")
        log_diagnostic("=" * 80)
        traceback_text = traceback.format_exc()
        for line in traceback_text.split('\n'):
            if line:
                log_diagnostic(line)
        
        # Try to show error dialog for frozen executable
        try:
            import tkinter as tk
            from tkinter import messagebox
            
            root = tk.Tk()
            root.withdraw()
            
            dialog_msg = f"MetaboGraph Failed to Start\n\n{error_msg}"
            if log_file:
                dialog_msg += f"\nDiagnostic log: {log_file}"
            
            messagebox.showerror("MetaboGraph Error", dialog_msg)
            root.destroy()
        except Exception as dialog_error:
            log_diagnostic(f"Failed to show error dialog: {dialog_error}")
        
        # Wait for input if running from console
        try:
            input("\nPress Enter to exit...")
        except Exception:
            pass
        
        sys.exit(1)
    
    finally:
        # Close log file
        try:
            if log_file and log_handle:
                log_diagnostic("Closing diagnostic log")
                log_handle.close()
        except Exception:
            pass
