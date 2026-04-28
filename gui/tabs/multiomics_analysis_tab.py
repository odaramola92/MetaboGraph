"""
Multi-Omics Integration Tab

This tab focuses exclusively on merging multiple pathway-ready datasets (for example,
metabolite + lipid outputs) into a single network-ready file that can be pushed directly
to the Pathway Network tab.
"""
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, scrolledtext, simpledialog
import logging
import os
import pandas as pd
import threading
from datetime import datetime

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# Import shared components
from gui.shared.base_tab import BaseTab, _setup_global_styles
from gui.shared.column_assignment import ColumnDetector

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MultiOmicsAnalysisTab(BaseTab):
    """Multi-Omics/Comparative Analysis Tab - Integrate or compare multiple datasets"""
    
    def __init__(self, parent, data_manager):
        """Initialize Multi-Omics Analysis tab"""
        super().__init__(parent, data_manager)
        
        # Setup global styles
        _setup_global_styles()
        
        # Get root window for dialogs
        self.root = parent.winfo_toplevel()
        
        # Store reference to parent (for accessing network tab later)
        self.parent_notebook = parent
        
        # Data storage
        self.uploaded_files = []  # List of dicts: {'name': str, 'path': str, 'df': DataFrame}
        self.merged_data = None
        
        # Create UI
        self.setup_ui()
        logger.info("Multi-Omics Analysis Tab initialized")
    
    def setup_ui(self):
        """Build the simplified multi-omics integration UI."""
        main_frame = tk.Frame(self.frame, bg='#f0f0f0')
        main_frame.pack(fill='both', expand=True, padx=10, pady=10)

        title_label = tk.Label(
            main_frame,
            text="🧬 Multi-Omics Integration",
            font=('Arial', 14, 'bold'),
            bg='#f0f0f0',
            fg='#2c3e50'
        )
        title_label.pack(pady=(0, 15))

        columns_frame = tk.Frame(main_frame, bg='#f0f0f0')
        columns_frame.pack(fill='both', expand=True)
        columns_frame.grid_columnconfigure(0, weight=1)
        columns_frame.grid_columnconfigure(1, weight=1)
        columns_frame.grid_rowconfigure(0, weight=1)

        left_frame = tk.LabelFrame(
            columns_frame,
            text="⚙️ Upload & Configuration",
            font=('Arial', 11, 'bold'),
            bg='#f0f0f0',
            fg='#2c3e50'
        )
        left_frame.grid(row=0, column=0, sticky='nsew', padx=(0, 5))

        left_canvas = tk.Canvas(left_frame, bg='#f0f0f0')
        left_scrollbar = ttk.Scrollbar(left_frame, orient="vertical", command=left_canvas.yview)
        left_content = tk.Frame(left_canvas, bg='#f0f0f0')

        left_content.bind(
            "<Configure>",
            lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all"))
        )

        canvas_window = left_canvas.create_window((0, 0), window=left_content, anchor="nw")
        left_canvas.configure(yscrollcommand=left_scrollbar.set)

        left_canvas.pack(side="left", fill="both", expand=True, padx=10, pady=10)
        left_scrollbar.pack(side="right", fill="y")

        def configure_scroll_region(event):
            left_canvas.configure(scrollregion=left_canvas.bbox("all"))
            left_canvas.itemconfig(canvas_window, width=event.width)

        left_canvas.bind('<Configure>', configure_scroll_region)

        def on_mousewheel(event):
            left_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        left_canvas.bind("<MouseWheel>", on_mousewheel)
        left_content.bind("<MouseWheel>", on_mousewheel)

        # overview_frame = tk.LabelFrame(
        #     left_content,
        #     text="🛠️ Workflow Overview",
        #     font=('Arial', 10, 'bold'),
        #     bg='#f0f0f0',
        #     fg='#2c3e50'
        # )
        # overview_frame.pack(fill='x', padx=10, pady=(10, 0))

        # overview_text = (
        #     "1. Upload two or more pathway-annotated datasets (e.g., metabolite + lipid).\n"
        #     "2. Confirm required columns via the detector dialog.\n"
        #     "3. Merge the validated datasets into a single Excel file and push it to the Network Analysis tab."
        # )
        # tk.Label(
        #     overview_frame,
        #     text=overview_text,
        #     bg='#f0f0f0',
        #     justify='left',
        #     wraplength=420,
        #     font=('Arial', 9)
        # ).pack(anchor='w', padx=10, pady=10)

        # === Upload Section ===
        upload_frame = tk.LabelFrame(
            left_content,
            text="📁 Upload Datasets",
            font=('Arial', 10, 'bold'),
            bg='#f0f0f0',
            fg='#2c3e50'
        )
        upload_frame.pack(fill='x', padx=10, pady=10)
        
        # Notice about unlimited uploads
        notice_frame = tk.Frame(upload_frame, bg='#e3f2fd', relief='solid', borderwidth=1)
        notice_frame.pack(fill='x', padx=10, pady=(10, 5))
        
        tk.Label(
            notice_frame,
            text="ℹ️ NOTICE: You can upload unlimited files for analysis",
            bg='#e3f2fd',
            font=('Arial', 9, 'bold'),
            fg='#1976d2'
        ).pack(pady=5)
        
        tk.Label(
            notice_frame,
            text="Each file should be a pathway-annotated dataset with:\n"
                 "• Feature ID column\n"
                 "• P-value column\n"
                 "• Log2FC column\n"
                 "• Pathway columns (optional but recommended)",
            bg='#e3f2fd',
            font=('Arial', 8),
            fg='#424242',
            justify='left'
        ).pack(pady=(0, 5), padx=10)
        
        # Upload button
        upload_btn_frame = tk.Frame(upload_frame, bg='#f0f0f0')
        upload_btn_frame.pack(fill='x', padx=10, pady=5)
        
        self.upload_button = tk.Button(
            upload_btn_frame,
            text="📤 Add Files",
            command=self._upload_files,
            bg='#2196f3',
            fg='white',
            font=('Arial', 9, 'bold'),
            relief='raised',
            borderwidth=2,
            padx=20,
            pady=5
        )
        self.upload_button.pack(side='left', padx=5)
        
        self.clear_files_button = tk.Button(
            upload_btn_frame,
            text="🗑️ Clear All",
            command=self._clear_files,
            bg='#f44336',
            fg='white',
            font=('Arial', 9, 'bold'),
            relief='raised',
            borderwidth=2,
            padx=20,
            pady=5,
            state='disabled'
        )
        self.clear_files_button.pack(side='left', padx=5)
        
        # Files list
        files_list_frame = tk.Frame(upload_frame, bg='#f0f0f0')
        files_list_frame.pack(fill='both', expand=True, padx=10, pady=(5, 10))
        
        tk.Label(
            files_list_frame,
            text="Uploaded Files:",
            bg='#f0f0f0',
            font=('Arial', 9, 'bold')
        ).pack(anchor='w', pady=(0, 5))
        
        # Treeview for files
        files_tree_frame = tk.Frame(files_list_frame, bg='#f0f0f0')
        files_tree_frame.pack(fill='both', expand=True)
        
        files_scrollbar = ttk.Scrollbar(files_tree_frame, orient="vertical")
        
        self.files_tree = ttk.Treeview(
            files_tree_frame,
            columns=('Dataset', 'Rows', 'Columns'),
            show='headings',
            height=6,
            yscrollcommand=files_scrollbar.set
        )
        files_scrollbar.config(command=self.files_tree.yview)
        
        self.files_tree.heading('Dataset', text='Dataset')
        self.files_tree.heading('Rows', text='Rows')
        self.files_tree.heading('Columns', text='Columns')
        
        self.files_tree.column('Dataset', width=200)
        self.files_tree.column('Rows', width=80, anchor='center')
        self.files_tree.column('Columns', width=80, anchor='center')
        
        self.files_tree.pack(side='left', fill='both', expand=True)
        files_scrollbar.pack(side='right', fill='y')
        
        # Remove selected file button and Rename Data button
        manage_btn_frame = tk.Frame(upload_frame, bg='#f0f0f0')
        manage_btn_frame.pack(fill='x', padx=10, pady=(0, 10))
        
        self.remove_file_button = tk.Button(
            manage_btn_frame,
            text="➖ Remove Selected",
            command=self._remove_selected_file,
            bg='#ff9800',
            fg='white',
            font=('Arial', 8, 'bold'),
            state='disabled',
            padx=8,
            pady=2
        )
        self.remove_file_button.pack(side='left', padx=(0, 5))
        
        self.rename_file_button = tk.Button(
            manage_btn_frame,
            text="✏️ Rename Data",
            command=self._rename_selected_file,
            bg="#4caf50",
            fg='white',
            font=('Arial', 8, 'bold'),
            state='disabled',
            padx=8,
            pady=2
        )
        self.rename_file_button.pack(side='left')
        
        process_frame = tk.Frame(left_content, bg='#f0f0f0')
        process_frame.pack(fill='x', padx=10, pady=20)

        self.process_button = tk.Button(
            process_frame,
            text="⚡ Process & Validate Data",
            command=self._process_data,
            bg='#4caf50',
            fg='white',
            font=('Arial', 11, 'bold'),
            relief='raised',
            borderwidth=3,
            padx=30,
            pady=10,
            state='disabled'
        )
        self.process_button.pack()
        
        # === 4. Merge & Export Section (Multi-Omic only) ===
        self.merge_frame = tk.LabelFrame(
            left_content,
            text="🔗 Merge & Export",
            font=('Arial', 10, 'bold'),
            bg='#f0f0f0',
            fg='#2c3e50'
        )
        self.merge_frame.pack(fill='x', padx=10, pady=10)
        
        tk.Label(
            self.merge_frame,
            text="After validation, merge datasets into a single file for network analysis",
            bg='#f0f0f0',
            font=('Arial', 8, 'italic'),
            fg='#7f8c8d',
            wraplength=400,
            justify='left'
        ).pack(padx=10, pady=(10, 5))
        
        self.merge_button = tk.Button(
            self.merge_frame,
            text="🔗 Merge & Export to Network Tab",
            command=self._merge_and_export,
            bg='#9c27b0',
            fg='white',
            font=('Arial', 10, 'bold'),
            relief='raised',
            borderwidth=2,
            padx=20,
            pady=8,
            state='disabled'
        )
        self.merge_button.pack(pady=(5, 10))
        
        # ============ RIGHT COLUMN: Progress Log ============
        right_frame = tk.LabelFrame(
            columns_frame,
            text="📋 Analysis Log",
            font=('Arial', 11, 'bold'),
            bg='#f0f0f0',
            fg='#2c3e50'
        )
        right_frame.grid(row=0, column=1, sticky='nsew', padx=(5, 0))
        
        self.log_text = scrolledtext.ScrolledText(
            right_frame,
            wrap=tk.WORD,
            font=('Courier New', 9),
            bg='#ffffff',
            fg='#000000',
            state='disabled',
            height=22
        )
        self.log_text.pack(fill='both', expand=True, padx=10, pady=10)

        self._log("Welcome to the Multi-Omics Integration workspace! 🎉")
        self._log("Upload at least two pathway-annotated datasets to get started.\n")
    
    def _log(self, message):
        """Add message to log"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_text.config(state='normal')
        self.log_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.log_text.see(tk.END)
        self.log_text.config(state='disabled')
    
    def _upload_files(self):
        """Upload multiple files"""
        try:
            file_paths = filedialog.askopenfilenames(
                title="Select Pathway-Annotated Files",
                filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
            )
            
            if not file_paths:
                return
            
            self._log(f"\n📂 Loading {len(file_paths)} file(s)...")
            
            for file_path in file_paths:
                filename = os.path.basename(file_path)
                
                # Check if already uploaded
                if any(f['path'] == file_path for f in self.uploaded_files):
                    self._log(f"⚠️ {filename} already uploaded, skipping")
                    continue
                
                try:
                    # Load Excel file
                    df = pd.read_excel(file_path)

                    label = self._prompt_dataset_label(filename)
                    if label is None:
                        self._log(f"⚠️ Skipped {filename} (no dataset name provided)")
                        continue
                    
                    # Store file info
                    file_info = {
                        'name': label,
                        'original_name': filename,
                        'path': file_path,
                        'df': df,
                        'validated': False
                    }
                    self.uploaded_files.append(file_info)
                    
                    # Add to treeview
                    self.files_tree.insert(
                        '',
                        'end',
                        values=(label, len(df), len(df.columns))
                    )
                    
                    self._log(f"✅ Loaded {label} ({filename}): {len(df)} rows, {len(df.columns)} columns")
                    
                except Exception as e:
                    self._log(f"❌ Failed to load {filename}: {str(e)}")
                    logger.error(f"Error loading {filename}: {e}", exc_info=True)
            
            if len(self.uploaded_files) > 0:
                self.clear_files_button.config(state='normal')
                self.remove_file_button.config(state='normal')
                self.rename_file_button.config(state='normal')

                if len(self.uploaded_files) < 2:
                    self._log("\n⚠️ Multi-omics merging requires at least 2 datasets (e.g., metabolite + lipid)")
                    self.process_button.config(state='disabled')
                else:
                    self.process_button.config(state='normal')
            
        except Exception as e:
            error_msg = f"Error uploading files: {str(e)}"
            self._log(f"❌ {error_msg}")
            messagebox.showerror("Upload Error", error_msg)
            logger.error(error_msg, exc_info=True)
    
    def _clear_files(self):
        """Clear all uploaded files"""
        if messagebox.askyesno("Clear Files", "Remove all uploaded files?"):
            self.uploaded_files.clear()
            self.files_tree.delete(*self.files_tree.get_children())
            self.clear_files_button.config(state='disabled')
            self.remove_file_button.config(state='disabled')
            self.rename_file_button.config(state='disabled')
            self.process_button.config(state='disabled')
            self.merge_button.config(state='disabled')
            self._log("\n🗑️ All files cleared")
    
    def _remove_selected_file(self):
        """Remove selected file from list"""
        selected = self.files_tree.selection()
        if not selected:
            messagebox.showwarning("No Selection", "Please select a file to remove")
            return
        
        # Get filename from treeview
        item = selected[0]
        values = self.files_tree.item(item, 'values')
        filename = values[0]
        
        # Remove from list
        self.uploaded_files = [f for f in self.uploaded_files if f['name'] != filename]
        
        # Remove from treeview
        self.files_tree.delete(item)
        
        self._log(f"➖ Removed: {filename}")
        
        # Update button states
        if len(self.uploaded_files) == 0:
            self.clear_files_button.config(state='disabled')
            self.remove_file_button.config(state='disabled')
            self.rename_file_button.config(state='disabled')
            self.process_button.config(state='disabled')
            self.merge_button.config(state='disabled')

    def _rename_selected_file(self):
        """Rename selected file in the list"""
        selected = self.files_tree.selection()
        if not selected:
            messagebox.showwarning("No Selection", "Please select a file to rename")
            return
        
        # Get current name from treeview
        item = selected[0]
        values = self.files_tree.item(item, 'values')
        old_name = values[0]
        
        # Get existing names excluding the current one
        existing = {f['name'] for f in self.uploaded_files if f['name'] != old_name}
        
        while True:
            new_name = simpledialog.askstring(
                "Rename Dataset",
                f"Enter a new name for '{old_name}':",
                initialvalue=old_name,
                parent=self.frame
            )
            
            if new_name is None:
                return
            
            new_name = new_name.strip()
            if not new_name:
                messagebox.showwarning("Invalid Name", "Dataset name cannot be empty.")
                continue
            
            if new_name in existing:
                messagebox.showwarning(
                    "Duplicate Name",
                    "A dataset with this name already exists. Please choose another."
                )
                continue
            
            break
        
        # Update the record
        for f in self.uploaded_files:
            if f['name'] == old_name:
                f['name'] = new_name
                break
        
        # Update treeview
        self.files_tree.item(item, values=(new_name, values[1], values[2]))
        self._log(f"✏️ Renamed: {old_name} → {new_name}")

    def _prompt_dataset_label(self, default_name):
        """Prompt the user for a custom dataset label"""
        existing_labels = {f['name'] for f in self.uploaded_files}
        while True:
            label = simpledialog.askstring(
                "Dataset Name",
                "Enter a name for this dataset:",
                initialvalue=default_name,
                parent=self.frame
            )

            if label is None:
                return None

            label = label.strip()
            if not label:
                messagebox.showwarning("Invalid Name", "Dataset name cannot be empty.")
                continue

            if label in existing_labels:
                messagebox.showwarning(
                    "Duplicate Name",
                    "A dataset with this name already exists. Please choose another."
                )
                continue

            return label
    
    def _process_data(self):
        """Process and validate uploaded data"""
        if len(self.uploaded_files) < 2:
            messagebox.showerror(
                "Insufficient Data",
                "Multi-omics merging requires at least 2 datasets (e.g., metabolite + lipid)"
            )
            return
        
        self._log("\n🔄 Starting data validation...")
        
        # Run validation in background thread
        thread = threading.Thread(target=self._validate_data_thread, daemon=True)
        thread.start()
    
    def _validate_data_thread(self):
        """Background thread for data validation"""
        from gui.shared.column_assignment import show_column_assignment_dialog
        
        try:
            required_cols = ['Feature ID', 'P-Value', 'Log2 Fold Change']
            all_valid = True
            
            for idx, file_info in enumerate(self.uploaded_files):
                filename = file_info['name']
                df = file_info['df']
                
                self.root.after(0, lambda f=filename: self._log(f"\n📊 Validating {f}..."))
                
                # Auto-detect required columns
                detected = {}
                
                # Feature ID
                feature_col = None
                for col in df.columns:
                    if ColumnDetector.detect_column(col, 'Feature ID'):
                        feature_col = col
                        break
                
                if feature_col:
                    detected['Feature ID'] = feature_col
                    self.root.after(0, lambda c=feature_col: self._log(f"  ✓ Feature ID: {c}"))
                else:
                    self.root.after(0, lambda: self._log(f"  ❌ Feature ID column not found"))
                
                # P-Value
                pvalue_candidates = []
                for col in df.columns:
                    if ColumnDetector.detect_column(col, 'P-Value'):
                        pvalue_candidates.append(col)
                
                pvalue_col = ColumnDetector.select_best_pvalue_column(pvalue_candidates) if pvalue_candidates else None
                
                if pvalue_col:
                    detected['P-Value'] = pvalue_col
                    priority_note = "(adjusted)" if any(keyword in pvalue_col.lower().replace('_', '').replace('-', '') for keyword in ['adjp', 'padj', 'fdr', 'qvalue', 'adjustedp']) else ''
                    priority_note = f" {priority_note}" if priority_note else ''
                    self.root.after(0, lambda c=pvalue_col, n=priority_note: self._log(f"  ✓ P-Value: {c}{n}"))
                    
                    # Note if multiple p-value columns exist
                    if len(pvalue_candidates) > 1:
                        self.root.after(0, lambda count=len(pvalue_candidates): self._log(f"    ℹ️ Note: {count} p-value columns detected, selected best match"))
                else:
                    self.root.after(0, lambda: self._log(f"  ❌ P-Value column not found"))
                
                # Log2 Fold Change
                log2fc_col = None
                log2fc_candidates = []
                for col in df.columns:
                    if ColumnDetector.detect_column(col, 'Log2 Fold Change'):
                        log2fc_candidates.append(col)
                        if log2fc_col is None:
                            log2fc_col = col
                
                if log2fc_col:
                    detected['Log2 Fold Change'] = log2fc_col
                    self.root.after(0, lambda c=log2fc_col: self._log(f"  ✓ Log2FC: {c}"))
                    
                    # Note if multiple log2FC columns exist
                    if len(log2fc_candidates) > 1:
                        self.root.after(0, lambda count=len(log2fc_candidates): self._log(f"    ℹ️ Note: {count} log2FC columns detected, selected first match"))
                else:
                    self.root.after(0, lambda: self._log(f"  ❌ Log2FC column not found"))
                
                # Check for pathway columns
                pathway_cols = [col for col in df.columns if 'pathway' in col.lower()]
                if pathway_cols:
                    self.root.after(0, lambda p=pathway_cols: self._log(f"  ✓ Found {len(p)} pathway column(s)"))
                else:
                    self.root.after(0, lambda: self._log(f"  ⚠️ No pathway columns found (optional)"))
                
                # Store detected columns temporarily
                file_info['detected_columns'] = detected
                
                # Show column assignment dialog for user confirmation/editing
                self.root.after(0, lambda: self._log(f"\n⏳ Please confirm column assignments for {filename}..."))
                
                # Show dialog on main thread
                result = [None]  # Use list to store result from main thread
                dialog_complete = [False]
                
                def show_dialog():
                    try:
                        res = show_column_assignment_dialog(
                            parent=self.root,
                            df=df,
                            tab_type='pathway',  # Requires Feature ID, P-Value, Log2FC
                            auto_calculate=False,
                            existing_assignments=detected,  # Pre-populate with detected values
                            excel_file_path=file_info['path']
                        )
                        result[0] = res
                    except Exception as e:
                        logger.error(f"Dialog error: {e}", exc_info=True)
                        result[0] = None
                    finally:
                        dialog_complete[0] = True
                
                self.root.after(0, show_dialog)
                
                # Wait for dialog to complete
                import time
                while not dialog_complete[0]:
                    time.sleep(0.1)
                
                if result[0]:
                    # User confirmed/edited assignments
                    confirmed = result[0]['assignments']
                    file_info['detected_columns'] = confirmed
                    file_info['validated'] = True
                    
                    self.root.after(0, lambda: self._log(f"✅ Column assignments confirmed for {filename}"))
                    self.root.after(0, lambda c=confirmed: self._log(f"    Feature ID: {c.get('Feature ID', 'N/A')}"))
                    self.root.after(0, lambda c=confirmed: self._log(f"    P-Value: {c.get('P-Value', 'N/A')}"))
                    self.root.after(0, lambda c=confirmed: self._log(f"    Log2FC: {c.get('Log2 Fold Change', 'N/A')}"))
                else:
                    # User cancelled
                    self.root.after(0, lambda: self._log(f"❌ Column assignment cancelled for {filename}"))
                    all_valid = False
                    file_info['validated'] = False
            
            # Update UI based on validation results
            validated_count = sum(1 for f in self.uploaded_files if f.get('validated', False))
            
            if all_valid and validated_count == len(self.uploaded_files):
                self.root.after(0, lambda: self._log(f"\n✅ All {validated_count} file(s) validated successfully!"))
                self.root.after(0, lambda: self._log("📝 You can now proceed with merging or analysis"))
                self.root.after(0, lambda: self.merge_button.config(state='normal'))
                self.root.after(0, lambda: messagebox.showinfo(
                    "Validation Complete",
                    f"✅ All {validated_count} dataset(s) validated!\n\n"
                    "Column assignments have been confirmed.\n\n"
                    "Click 'Merge & Export to Network Tab' to create the unified dataset."
                ))
            elif validated_count > 0:
                self.root.after(0, lambda: self._log(f"\n⚠️ Only {validated_count} of {len(self.uploaded_files)} files validated"))
                self.root.after(0, lambda: messagebox.showwarning(
                    "Partial Validation",
                    f"Only {validated_count} of {len(self.uploaded_files)} files were validated.\n\n"
                    "Some column assignments were cancelled.\n\n"
                    "Remove invalid files or try validation again."
                ))
            else:
                self.root.after(0, lambda: self._log("\n❌ Validation cancelled or failed"))
                self.root.after(0, lambda: messagebox.showerror(
                    "Validation Failed",
                    "No files were successfully validated.\n\n"
                    "Column assignments were cancelled or required columns are missing."
                ))
        
        except Exception as e:
            error_msg = f"Validation error: {str(e)}"
            self.root.after(0, lambda: self._log(f"❌ {error_msg}"))
            logger.error(error_msg, exc_info=True)
    
    def _merge_and_export(self):
        """Merge datasets and export to network tab"""
        self._log("\n🔗 Starting merge process...")
        
        # Run merge in background thread
        thread = threading.Thread(target=self._merge_thread, daemon=True)
        thread.start()

    def _merge_thread(self):
        """Background thread for merging data"""
        try:
            if len(self.uploaded_files) < 2:
                self.root.after(0, lambda: messagebox.showerror(
                    "Insufficient Data",
                    "Need at least 2 datasets to merge"
                ))
                return
            
            self.root.after(0, lambda: self._log("📊 Preparing datasets for merge..."))
            
            # Prepare datasets with standardized column names
            standardized_dfs = []
            
            for file_info in self.uploaded_files:
                if not file_info.get('validated', False):
                    continue
                
                df = file_info['df'].copy()
                detected = file_info['detected_columns']
                filename = file_info['name']
                
                # Rename to standard names
                rename_map = {}
                if 'Feature ID' in detected:
                    rename_map[detected['Feature ID']] = 'Feature_ID'
                if 'P-Value' in detected:
                    rename_map[detected['P-Value']] = 'pvalue'
                if 'Log2 Fold Change' in detected:
                    rename_map[detected['Log2 Fold Change']] = 'log2FC'
                
                df.rename(columns=rename_map, inplace=True)
                
                # Add source column to track origin
                source_name = filename
                df['Source'] = source_name
                
                # Get pathway columns
                pathway_cols = [col for col in df.columns if 'pathway' in col.lower()]
                
                # Get upstream columns (including enzymes and transporters)
                upstream_cols = []
                for col in df.columns:
                    col_lower = col.lower()
                    if any(keyword in col_lower for keyword in ['upstream', 'enzyme', 'transporter', 'regulator']):
                        upstream_cols.append(col)
                
                # Get disease columns
                disease_cols = [col for col in df.columns if 'disease' in col.lower()]
                
                # Check for Class_name column (for lipid class compression in network analysis)
                class_name_col = None
                for col in df.columns:
                    if col.lower() == 'class_name':
                        class_name_col = col
                        break
                
                # Keep only necessary columns
                keep_cols = ['Feature_ID', 'pvalue', 'log2FC', 'Source']
                keep_cols.extend(pathway_cols)
                keep_cols.extend(upstream_cols)
                keep_cols.extend(disease_cols)
                
                # Add Class_name if available (needed for network class compression)
                if class_name_col:
                    keep_cols.append(class_name_col)
                
                # Filter to existing columns
                keep_cols = [col for col in keep_cols if col in df.columns]
                df_filtered = df[keep_cols]
                
                standardized_dfs.append(df_filtered)
                
                # Log what was prepared, including Class_name status
                class_status = " (Class_name ✓)" if class_name_col and class_name_col in df_filtered.columns else ""
                self.root.after(0, lambda s=source_name, r=len(df_filtered), p=len(pathway_cols), u=len(upstream_cols), d=len(disease_cols), cs=class_status: 
                               self._log(f"  ✓ Prepared {s}: {r} rows, {p} pathway cols, {u} upstream cols, {d} disease cols{cs}"))
            
            if len(standardized_dfs) < 2:
                self.root.after(0, lambda: messagebox.showerror(
                    "Error",
                    "Not enough validated datasets to merge"
                ))
                return
            
            # Merge datasets
            self.root.after(0, lambda: self._log("\n🔄 Merging datasets..."))
            
            # Use outer join to keep all features from both datasets
            merged = standardized_dfs[0]
            
            for i, df in enumerate(standardized_dfs[1:], 1):
                self.root.after(0, lambda idx=i: self._log(f"  Merging dataset {idx+1}..."))
                
                # Merge on Feature_ID
                merged = pd.merge(
                    merged,
                    df,
                    on='Feature_ID',
                    how='outer',
                    suffixes=('', f'_{i}')
                )
            
            # Consolidate duplicate columns
            self.root.after(0, lambda: self._log("\n🔧 Consolidating columns..."))
            
            # For pvalue, keep the most significant (minimum) value
            if 'pvalue' in merged.columns:
                pvalue_cols = [col for col in merged.columns if col.startswith('pvalue')]
                if len(pvalue_cols) > 1:
                    # Use min to get most significant p-value, ignoring NaN
                    merged['pvalue'] = merged[pvalue_cols].min(axis=1, skipna=True)
                    # Drop suffixed versions
                    merged.drop([col for col in pvalue_cols if col != 'pvalue'], axis=1, inplace=True)
            
            # For log2FC, keep the one with largest absolute value
            if 'log2FC' in merged.columns:
                log2fc_cols = [col for col in merged.columns if col.startswith('log2FC')]
                if len(log2fc_cols) > 1:
                    # For each row, find the value with largest absolute value across all log2FC columns
                    def get_max_abs_value(row):
                        values = row[log2fc_cols]
                        # Drop NaN values
                        valid_values = values.dropna()
                        if len(valid_values) == 0:
                            return pd.NA
                        # Return value with largest absolute magnitude
                        abs_values = valid_values.abs()
                        max_idx = abs_values.idxmax()
                        return valid_values[max_idx]
                    
                    merged['log2FC'] = merged.apply(get_max_abs_value, axis=1)
                    # Drop suffixed versions
                    merged.drop([col for col in log2fc_cols if col != 'log2FC'], axis=1, inplace=True)
            
            # Consolidate Source column
            source_cols = [col for col in merged.columns if col.startswith('Source')]
            if len(source_cols) > 1:
                # Combine sources, removing empty strings
                merged['Source'] = merged[source_cols].fillna('').apply(
                    lambda x: ' + '.join([str(v) for v in x if str(v).strip() and str(v) != 'nan']),
                    axis=1
                )
                merged.drop([col for col in source_cols if col != 'Source'], axis=1, inplace=True)
            
            # Helper to strip only numeric merge suffixes (e.g. "HMDB_Pathways_1" -> "HMDB_Pathways")
            def _get_base_name(col_name: str) -> str:
                if '_' in col_name:
                    left, right = col_name.rsplit('_', 1)
                    if right.isdigit():
                        return left
                return col_name

            # Consolidate pathway columns
            pathway_base_names = set()
            for col in merged.columns:
                if 'pathway' in col.lower():
                    base_name = _get_base_name(col)
                    pathway_base_names.add(base_name)
            
            for base_name in pathway_base_names:
                cols = [col for col in merged.columns if col.startswith(base_name)]
                if len(cols) > 1:
                    # Concatenate pathway values with semicolon, removing duplicates
                    def merge_pathways(row):
                        all_pathways = []
                        for col in cols:
                            val = str(row[col])
                            if val and val.strip() and val != 'nan':
                                # Split by semicolon in case there are multiple pathways
                                pathways = [p.strip() for p in val.split(';') if p.strip()]
                                all_pathways.extend(pathways)
                        # Remove duplicates while preserving order
                        unique_pathways = []
                        seen = set()
                        for p in all_pathways:
                            if p not in seen:
                                unique_pathways.append(p)
                                seen.add(p)
                        return ';'.join(unique_pathways) if unique_pathways else ''
                    
                    merged[base_name] = merged.apply(merge_pathways, axis=1)
                    # Drop suffixed versions
                    merged.drop([col for col in cols if col != base_name], axis=1, inplace=True)
            
            # Consolidate upstream-related columns (upstream, enzyme, transporter, regulator)
            upstream_base_names = set()
            for col in merged.columns:
                col_lower = col.lower()
                if any(keyword in col_lower for keyword in ['upstream', 'enzyme', 'transporter', 'regulator']):
                    base_name = _get_base_name(col)
                    upstream_base_names.add(base_name)

            for base_name in upstream_base_names:
                cols = [col for col in merged.columns if col.startswith(base_name)]
                if len(cols) > 1:
                    # Concatenate upstream annotations with semicolon, removing duplicates
                    def merge_upstream(row):
                        all_items = []
                        for col in cols:
                            val = str(row[col])
                            if val and val.strip() and val != 'nan':
                                parts = [p.strip() for p in val.split(';') if p.strip()]
                                all_items.extend(parts)
                        unique_items = []
                        seen = set()
                        for item in all_items:
                            if item not in seen:
                                unique_items.append(item)
                                seen.add(item)
                        return ';'.join(unique_items) if unique_items else ''

                    merged[base_name] = merged.apply(merge_upstream, axis=1)
                    merged.drop([col for col in cols if col != base_name], axis=1, inplace=True)

            # Consolidate disease-related columns
            disease_base_names = set()
            for col in merged.columns:
                if 'disease' in col.lower():
                    base_name = _get_base_name(col)
                    disease_base_names.add(base_name)

            for base_name in disease_base_names:
                cols = [col for col in merged.columns if col.startswith(base_name)]
                if len(cols) > 1:
                    # Concatenate disease annotations with semicolon, removing duplicates
                    def merge_disease(row):
                        all_items = []
                        for col in cols:
                            val = str(row[col])
                            if val and val.strip() and val != 'nan':
                                parts = [p.strip() for p in val.split(';') if p.strip()]
                                all_items.extend(parts)
                        unique_items = []
                        seen = set()
                        for item in all_items:
                            if item not in seen:
                                unique_items.append(item)
                                seen.add(item)
                        return ';'.join(unique_items) if unique_items else ''

                    merged[base_name] = merged.apply(merge_disease, axis=1)
                    merged.drop([col for col in cols if col != base_name], axis=1, inplace=True)
            
            # Consolidate Class_name column (for lipid class compression)
            class_name_cols = [col for col in merged.columns if col.lower() == 'class_name' or col.startswith('Class_name')]
            if len(class_name_cols) > 1:
                # For Class_name, take the first non-empty value (they should be consistent for same feature)
                def merge_class_name(row):
                    for col in class_name_cols:
                        val = str(row[col])
                        if val and val.strip() and val != 'nan' and val.lower() != 'none':
                            return val
                    return ''
                
                merged['Class_name'] = merged.apply(merge_class_name, axis=1)
                # Drop suffixed versions
                merged.drop([col for col in class_name_cols if col != 'Class_name'], axis=1, inplace=True)
                self.root.after(0, lambda: self._log("  ✓ Consolidated Class_name column"))
            elif len(class_name_cols) == 1 and class_name_cols[0] != 'Class_name':
                # Rename single Class_name column to standard name
                merged.rename(columns={class_name_cols[0]: 'Class_name'}, inplace=True)
                self.root.after(0, lambda: self._log("  ✓ Standardized Class_name column"))
            
            # Create All_Pathways column (required by Network Analysis tab)
            self.root.after(0, lambda: self._log("\n🔗 Creating combined pathway column..."))
            
            # Identify all pathway-related columns (excluding disease, upstream, and pathway metrics)
            all_pathway_cols = []
            for col in merged.columns:
                col_lower = col.lower()
                # Include pathway columns but exclude disease, upstream, and pathway summary metrics
                if ('pathway' in col_lower or col in ['HMDB', 'KEGG', 'SMPDB', 'Reactome', 'WikiPathways', 
                                                        'PathBank', 'Metabolika', 'wikipathways']) \
                   and 'disease' not in col_lower and 'upstream' not in col_lower \
                   and not col_lower.startswith('pathway_sources') and not col_lower.startswith('total_pathways'):
                    all_pathway_cols.append(col)
            
            if all_pathway_cols:
                self.root.after(0, lambda cols=all_pathway_cols: self._log(f"  📝 Combining pathways from: {', '.join(cols)}"))
                
                def combine_all_pathways(row):
                    """Combine all pathways from different database columns"""
                    all_pathways = []
                    
                    # Database labels to exclude (these are not pathway names)
                    database_labels = {
                        'hmdb', 'kegg', 'smpdb', 'reactome', 'wikipathways', 'pathbank', 
                        'metabolika', 'hmdb_id', 'kegg_id', 'pathway', 'pathways',
                        'hmdb, pathbank, smpdb', 'hmdb, kegg', 'kegg, smpdb'
                    }
                    
                    # Keywords to exclude from pathway names
                    excluded_keywords = [
                        'plant', 'plants', 'arabidopsis', 'rice', 'maize', 'wheat',
                        'indole-3-acetate', 'auxin', 'auxin conjugate', 'iaa biosynthesis',
                        'action', 'toxicity', 'drug', 'medication', 'therapy', 
                        'treatment', 'pharmacology', 'yeast', 'bacteria', 'bacterial'
                    ]
                    
                    # Disease/disorder keywords that can be overridden by pathway keywords
                    disease_keywords = ['disease', 'disorder', 'syndrome', 'deficiency']
                    
                    # Pathway keywords that indicate it's a valid pathway (override disease exclusion)
                    pathway_keywords = [
                        'metabolism', 'metabolic', 'biosynthesis', 'degradation',
                        'pathway', 'cycle', 'synthesis', 'oxidation', 'reduction',
                        'transport', 'signaling', 'activation'
                    ]
                    
                    for col in all_pathway_cols:
                        val = str(row[col])
                        if val and val.strip() and val != 'nan' and val.lower() != 'none':
                            # Split by semicolon in case there are multiple pathways
                            pathways = [p.strip() for p in val.split(';') if p.strip()]
                            
                            for p in pathways:
                                # Skip if it's a database label or too short
                                p_lower = p.lower().replace(' ', '').replace(',', '')
                                if p_lower not in database_labels and len(p) > 3:
                                    # Additional check: exclude if it's just a comma-separated list of DB names
                                    parts = [part.strip().lower() for part in p.split(',')]
                                    is_db_list = all(part in database_labels for part in parts)
                                    
                                    # Check if pathway contains excluded keywords
                                    p_check = p.lower()
                                    
                                    # Check for basic excluded keywords (always exclude)
                                    has_basic_excluded = any(keyword in p_check for keyword in excluded_keywords)
                                    
                                    # Check for disease keywords
                                    has_disease = any(keyword in p_check for keyword in disease_keywords)
                                    
                                    # Check for pathway keywords (override disease exclusion)
                                    has_pathway_keyword = any(keyword in p_check for keyword in pathway_keywords)
                                    
                                    # Exclude if:
                                    # 1. It's a database list, OR
                                    # 2. Has basic excluded keywords, OR
                                    # 3. Has disease keywords BUT NOT pathway keywords
                                    should_exclude = is_db_list or has_basic_excluded or (has_disease and not has_pathway_keyword)
                                    
                                    if not should_exclude:
                                        all_pathways.append(p)
                    
                    # Remove duplicates while preserving order
                    unique_pathways = []
                    seen = set()
                    for p in all_pathways:
                        if p not in seen:
                            unique_pathways.append(p)
                            seen.add(p)
                    
                    return ';'.join(unique_pathways) if unique_pathways else ''
                
                merged['All_Pathways'] = merged.apply(combine_all_pathways, axis=1)
                self.root.after(0, lambda: self._log(f"  ✅ All_Pathways column created"))
            else:
                self.root.after(0, lambda: self._log(f"  ⚠️ No pathway columns found to combine"))
            
            # After generic disease consolidation, ensure Associated_Diseases includes Reactome_Disease (no duplicates)
            if 'Associated_Diseases' in merged.columns or 'Reactome_Disease' in merged.columns:
                def _merge_disease_population(row):
                    vals = []
                    for col in ['Associated_Diseases', 'Reactome_Disease']:
                        if col in merged.columns:
                            v = row.get(col, '')
                            if v is None:
                                continue
                            s = str(v)
                            if not s or s.strip().lower() in ('nan', 'none'):
                                continue
                            # Split on common separators
                            parts = []
                            for sep in ['|', ';']:
                                if sep in s:
                                    for p in s.split(sep):
                                        p = p.strip()
                                        if p:
                                            parts.append(p)
                                    break
                            if not parts:
                                parts = [s.strip()]
                            vals.extend(parts)
                    # De-duplicate while preserving order
                    seen = set()
                    out = []
                    for d in vals:
                        if d not in seen:
                            seen.add(d)
                            out.append(d)
                    return ' | '.join(out) if out else ''

                merged['Associated_Diseases'] = merged.apply(_merge_disease_population, axis=1)
            
            # Remove rows without p-values (no statistics)
            if 'pvalue' in merged.columns:
                before_rows = len(merged)
                merged = merged[merged['pvalue'].notna()]
                removed_rows = before_rows - len(merged)
                if removed_rows > 0:
                    self.root.after(0, lambda r=removed_rows: self._log(f"  🧹 Removed {r} rows without p-value"))

            # Remove rows with missing Feature_ID
            merged.dropna(subset=['Feature_ID'], inplace=True)
            
            # Final cleanup of pathway summary metrics across datasets
            # Consolidate pathway_sources columns
            ps_cols = [c for c in merged.columns if c.startswith('pathway_sources')]
            if len(ps_cols) > 1:
                def _merge_pathway_sources(row):
                    tokens = []
                    for c in ps_cols:
                        v = row.get(c, '')
                        if v is None:
                            continue
                        s = str(v).strip()
                        if not s or s.lower() in ('nan', 'none'):
                            continue
                        # Split on common separators
                        parts = []
                        for sep in [',', ';', '|']:
                            if sep in s:
                                for p in s.split(sep):
                                    p = p.strip()
                                    if p:
                                        parts.append(p)
                                break
                        if not parts:
                            parts = [s]
                        tokens.extend(parts)
                    seen = set()
                    out = []
                    for t in tokens:
                        if t not in seen:
                            seen.add(t)
                            out.append(t)
                    return ', '.join(out) if out else ''

                merged['pathway_sources'] = merged.apply(_merge_pathway_sources, axis=1)
                merged.drop([c for c in ps_cols if c != 'pathway_sources'], axis=1, inplace=True)
            elif len(ps_cols) == 1 and ps_cols[0] != 'pathway_sources':
                merged.rename(columns={ps_cols[0]: 'pathway_sources'}, inplace=True)

            # Consolidate total_pathways columns (take max across datasets)
            tp_cols = [c for c in merged.columns if c.startswith('total_pathways')]
            if len(tp_cols) > 1:
                merged['total_pathways'] = merged[tp_cols].max(axis=1, skipna=True)
                merged.drop([c for c in tp_cols if c != 'total_pathways'], axis=1, inplace=True)
            elif len(tp_cols) == 1 and tp_cols[0] != 'total_pathways':
                merged.rename(columns={tp_cols[0]: 'total_pathways'}, inplace=True)

            # Rename back to expected format for network tab
            merged.rename(columns={'Feature_ID': 'Name'}, inplace=True)

            # Reorder and ensure expected columns exist (fill missing as empty/NaN)
            expected_cols = [
                'Name', 'Class_name', 'pvalue', 'log2FC',
                'HMDB_Pathways', 'Enzymes_Accession', 'Enzyme_Gene_name',
                'Associated_Diseases', 'Transporters', 'Transporter_Gene_name', 'Transporter_Uniprot_ID',
                'PathBank_Pathways', 'SMPDB_Pathways', 'WikiPathways', 'Metabolika_Pathways',
                'Reactome_Pathways', 'Reactome_Disease', 'KEGG_Pathways',
                'All_Pathways', 'pathway_sources', 'total_pathways'
            ]

            for col in expected_cols:
                if col not in merged.columns:
                    # Use empty string for text cols and NaN for numeric-ish
                    if col in ['pvalue', 'log2FC', 'total_pathways']:
                        merged[col] = pd.NA
                    else:
                        merged[col] = ''

            # Place expected columns first, then any extras
            ordered_cols = [c for c in expected_cols if c in merged.columns]
            remaining_cols = [c for c in merged.columns if c not in ordered_cols]
            merged = merged[ordered_cols + remaining_cols]

            self.merged_data = merged
            
            # Check if Class_name is in final dataset
            has_class_name = 'Class_name' in merged.columns
            
            self.root.after(0, lambda: self._log(f"\n✅ Merge complete!"))
            self.root.after(0, lambda: self._log(f"📊 Merged dataset: {len(merged)} rows, {len(merged.columns)} columns"))
            
            if has_class_name:
                non_empty_class = merged['Class_name'].notna().sum()
                self.root.after(0, lambda n=non_empty_class: self._log(f"🏷️  Class_name preserved: {n}/{len(merged)} rows have class annotations"))
            
            # Export to Excel (must be called on main thread)
            self.root.after(100, self._export_merged_data)
            
        except Exception as e:
            error_msg = f"Merge error: {str(e)}"
            self.root.after(0, lambda: self._log(f"❌ {error_msg}"))
            self.root.after(0, lambda: messagebox.showerror("Merge Error", error_msg))
            logger.error(error_msg, exc_info=True)

    def _export_merged_data(self):
        """Export merged data to Excel and load to network tab"""
        try:
            # Ask user where to save
            save_path = filedialog.asksaveasfilename(
                title="Save Merged Dataset",
                defaultextension=".xlsx",
                filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")],
                initialfile="MultiOmic_Merged.xlsx"
            )
            
            if not save_path:
                self._log("❌ Export cancelled")
                return
            
            # Save to Excel
            self._log(f"\n💾 Saving to: {os.path.basename(save_path)}")
            self.merged_data.to_excel(save_path, index=False)
            self._log("✅ File saved successfully!")
            
            # Try to load into network tab
            try:
                # Use stored reference to network tab (set by parent)
                if hasattr(self, 'network_tab_instance') and self.network_tab_instance:
                    network_tab = self.network_tab_instance
                    
                    self._log("\n🔄 Loading merged data to Network Analysis tab...")

                    # Load from the saved file path (not in-memory DataFrame) so auto-load
                    # and manual Browse import follow the exact same parsing path/results.
                    network_tab.load_annotated_data(save_path)
                    
                    self._log("✅ Data loaded to Network Analysis tab!")
                    
                    # Switch to network tab (selection by frame to avoid index assumptions)
                    try:
                        network_frame = self.network_tab_instance.get_frame()
                        self.parent_notebook.select(network_frame)
                    except Exception:
                        logger.debug("Could not auto-select the network tab frame", exc_info=True)
                    
                    messagebox.showinfo(
                        "Success",
                        f"✅ Merged dataset created and loaded!\n\n"
                        f"File: {os.path.basename(save_path)}\n"
                        f"Rows: {len(self.merged_data)}\n"
                        f"Columns: {len(self.merged_data.columns)}\n\n"
                        f"The data has been loaded to the Network Analysis tab.\n"
                        f"Please verify columns and proceed with pathway enrichment."
                    )
                else:
                    self._log("⚠️ Could not auto-load to Network tab")
                    messagebox.showinfo(
                        "Success",
                        f"✅ Merged dataset created!\n\n"
                        f"File: {os.path.basename(save_path)}\n"
                        f"Rows: {len(self.merged_data)}\n"
                        f"Columns: {len(self.merged_data.columns)}\n\n"
                        f"Please import this file in the Network Analysis tab."
                    )
            
            except Exception as e:
                logger.error(f"Error loading to network tab: {e}", exc_info=True)
                self._log(f"⚠️ Could not auto-load to Network tab: {str(e)}")
                messagebox.showinfo(
                    "Partial Success",
                    f"✅ Merged dataset saved to:\n{save_path}\n\n"
                    f"Please manually import this file in the Network Analysis tab."
                )
        
        except Exception as e:
            error_msg = f"Export error: {str(e)}"
            self._log(f"❌ {error_msg}")
            messagebox.showerror("Export Error", error_msg)
            logger.error(error_msg, exc_info=True)
