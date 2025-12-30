# runtime_license_hook.py
# PyInstaller runtime hook: it runs before your bundled app's main code.
# It must be present at build time and path resolved in .spec runtime_hooks.

# License check functionality - currently disabled as license_check module is not implemented
# To enable license checking in the future:
# 1. Create a license_check.py module with verify_license() function
# 2. Uncomment the code below

# try:
#     import license_check
#     license_check.verify_license()
# except Exception as e:
#     # If any unexpected error happens, attempt to show an error and exit.
#     try:
#         import tkinter as tk
#         from tkinter import messagebox
#         root = tk.Tk()
#         root.withdraw()
#         messagebox.showerror("Startup Error", f"License check failed: {e}")
#         root.destroy()
#     except Exception:
#         pass
#     import sys
#     sys.exit(1)

# For now, just pass through without license check
pass
