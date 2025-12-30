"""
Comparative Analysis Tab

Allows users to compare pathway-level outputs across multiple datasets.
Uses Z-Score based heatmap visualization with blue (negative) to red (positive)
color scale, and provides multiple visualization options.
"""
from __future__ import annotations

import logging
import os
from typing import List, Optional

import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext, simpledialog

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

from gui.shared.base_tab import BaseTab, _setup_global_styles

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ComparativeAnalysisTab(BaseTab):
    """Standalone GUI tab for pathway-level comparative analysis."""

    def __init__(self, parent, data_manager):
        super().__init__(parent, data_manager)
        _setup_global_styles()

        self.root = parent.winfo_toplevel()

        self.uploaded_files: List[dict] = []
        self.results_table: Optional[pd.DataFrame] = None
        self.heatmap_matrix: Optional[pd.DataFrame] = None
        self._current_figure: Optional[Figure] = None
        self._heatmap_canvas: Optional[FigureCanvasTkAgg] = None
        
        # Store all generated figures
        self._all_figures: dict = {}  # {viz_type: Figure}
        self._all_matrices: dict = {}  # {viz_type: DataFrame}

        # Variables for settings
        self.pvalue_threshold = tk.StringVar(value="0.05")
        self.viz_type = tk.StringVar(value="Z-Score Heatmap")
        self.fig_width = tk.StringVar(value="10")
        self.fig_height = tk.StringVar(value="8")
        self.top_n_pathways = tk.StringVar(value="50")
        self.pathway_filter = tk.StringVar(value="All Pathways")

        self.setup_ui()
        self._log("Welcome to Comparative Analysis! Upload pathway exports to begin.\n")

    # ------------------------------------------------------------------
    # UI assembly
    # ------------------------------------------------------------------
    def setup_ui(self):
        main_frame = tk.Frame(self.frame, bg="#f0f0f0")
        main_frame.pack(fill="both", expand=True, padx=10, pady=10)

        title = tk.Label(
            main_frame,
            text="📊 Comparative Pathway Analysis",
            font=("Arial", 14, "bold"),
            bg="#f0f0f0",
            fg="#2c3e50",
        )
        title.pack(pady=(0, 15))

        columns_frame = tk.Frame(main_frame, bg="#f0f0f0")
        columns_frame.pack(fill="both", expand=True)
        columns_frame.grid_columnconfigure(0, weight=1)
        columns_frame.grid_columnconfigure(1, weight=2)
        columns_frame.grid_rowconfigure(0, weight=1)

        # ---------------- Left column ----------------
        left_frame = tk.LabelFrame(
            columns_frame,
            text="⚙️ Inputs & Settings",
            font=("Arial", 11, "bold"),
            bg="#f0f0f0",
            fg="#2c3e50",
        )
        left_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 5))

        left_canvas = tk.Canvas(left_frame, bg="#f0f0f0")
        left_scroll = ttk.Scrollbar(left_frame, orient="vertical", command=left_canvas.yview)
        left_content = tk.Frame(left_canvas, bg="#f0f0f0")

        left_content.bind(
            "<Configure>",
            lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all")),
        )
        canvas_window = left_canvas.create_window((0, 0), window=left_content, anchor="nw")
        left_canvas.configure(yscrollcommand=left_scroll.set)

        left_canvas.pack(side="left", fill="both", expand=True, padx=10, pady=10)
        left_scroll.pack(side="right", fill="y")

        def _sync_width(event):
            left_canvas.itemconfig(canvas_window, width=event.width)

        left_canvas.bind("<Configure>", _sync_width)

        def _on_mousewheel(event):
            left_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        left_canvas.bind("<MouseWheel>", _on_mousewheel)
        left_content.bind("<MouseWheel>", _on_mousewheel)

        overview = tk.LabelFrame(
            left_content,
            text="🛠️ Workflow",
            font=("Arial", 10, "bold"),
            bg="#f0f0f0",
            fg="#2c3e50",
        )
        overview.pack(fill="x", padx=10, pady=(5, 0))
        tk.Label(
            overview,
            text=(
                "1. Export pathway results for each comparison using the Network Analysis tab.\n"
                "2. Upload each workbook here.\n"
                "3. Run analysis to compare pathways across datasets using Z-Score heatmap."
            ),
            font=("Arial", 9),
            justify="left",
            wraplength=380,
            bg="#f0f0f0",
        ).pack(anchor="w", padx=10, pady=10)

        # === Upload Section ===
        upload_frame = tk.LabelFrame(
            left_content,
            text="📁 Upload Pathway Exports",
            font=("Arial", 10, "bold"),
            bg="#f0f0f0",
            fg="#2c3e50",
        )
        upload_frame.pack(fill="x", padx=10, pady=10)

        btn_row = tk.Frame(upload_frame, bg="#f0f0f0")
        btn_row.pack(fill="x", padx=10, pady=(10, 5))

        tk.Button(
            btn_row,
            text="📤 Add Workbook(s)",
            command=self._upload_files,
            bg="#2196f3",
            fg="white",
            font=("Arial", 9, "bold"),
            relief="raised",
            padx=12,
            pady=4,
        ).pack(side="left", padx=(0, 5))

        self.clear_button = tk.Button(
            btn_row,
            text="🗑️ Clear All",
            command=self._clear_files,
            bg="#f44336",
            fg="white",
            font=("Arial", 9, "bold"),
            relief="raised",
            state="disabled",
            padx=12,
            pady=4,
        )
        self.clear_button.pack(side="left")

        files_frame = tk.Frame(upload_frame, bg="#f0f0f0")
        files_frame.pack(fill="both", expand=True, padx=10, pady=(0, 5))

        self.files_tree = ttk.Treeview(
            files_frame,
            columns=("Dataset", "Rows", "Columns"),
            show="headings",
            height=5,
        )
        self.files_tree.heading("Dataset", text="Dataset")
        self.files_tree.heading("Rows", text="Rows")
        self.files_tree.heading("Columns", text="Columns")
        self.files_tree.column("Dataset", width=180)
        self.files_tree.column("Rows", width=60, anchor="center")
        self.files_tree.column("Columns", width=60, anchor="center")

        tree_scroll = ttk.Scrollbar(files_frame, orient="vertical", command=self.files_tree.yview)
        self.files_tree.configure(yscrollcommand=tree_scroll.set)
        self.files_tree.pack(side="left", fill="both", expand=True)
        tree_scroll.pack(side="right", fill="y")

        # === Remove/Rename buttons ===
        manage_btn_row = tk.Frame(upload_frame, bg="#f0f0f0")
        manage_btn_row.pack(fill="x", padx=10, pady=(0, 10))

        self.remove_button = tk.Button(
            manage_btn_row,
            text="➖ Remove Selected",
            command=self._remove_selected,
            bg="#ff9800",
            fg="white",
            font=("Arial", 8, "bold"),
            relief="raised",
            state="disabled",
            padx=10,
            pady=3,
        )
        self.remove_button.pack(side="left", padx=(0, 5))

        self.rename_button = tk.Button(
            manage_btn_row,
            text="✏️ Rename Data",
            command=self._rename_selected,
            bg="#607d8b",
            fg="white",
            font=("Arial", 8, "bold"),
            relief="raised",
            state="disabled",
            padx=10,
            pady=3,
        )
        self.rename_button.pack(side="left")

        # === Settings Section ===
        settings_frame = tk.LabelFrame(
            left_content,
            text="🎯 Analysis Settings",
            font=("Arial", 10, "bold"),
            bg="#f0f0f0",
            fg="#2c3e50",
        )
        settings_frame.pack(fill="x", padx=10, pady=10)

        # P-value threshold
        pval_row = tk.Frame(settings_frame, bg="#f0f0f0")
        pval_row.pack(fill="x", padx=10, pady=(10, 5))
        tk.Label(pval_row, text="P-value threshold:", bg="#f0f0f0", font=("Arial", 9)).pack(side="left")
        tk.Entry(pval_row, textvariable=self.pvalue_threshold, width=8).pack(side="left", padx=(5, 0))

        # Visualization type dropdown
        viz_row = tk.Frame(settings_frame, bg="#f0f0f0")
        viz_row.pack(fill="x", padx=10, pady=5)
        tk.Label(viz_row, text="Visualization:", bg="#f0f0f0", font=("Arial", 9)).pack(side="left")
        self.viz_dropdown = ttk.Combobox(
            viz_row,
            textvariable=self.viz_type,
            values=["Z-Score Heatmap", "P-Value Heatmap", "Pathway Overlap Bar", "Venn Diagram"],
            state="readonly",
            width=18,
        )
        self.viz_dropdown.pack(side="left", padx=(5, 0))
        self.viz_dropdown.bind("<<ComboboxSelected>>", self._on_viz_changed)

        # Top N pathways control
        topn_row = tk.Frame(settings_frame, bg="#f0f0f0")
        topn_row.pack(fill="x", padx=10, pady=5)
        tk.Label(topn_row, text="Top N pathways:", bg="#f0f0f0", font=("Arial", 9)).pack(side="left")
        tk.Entry(topn_row, textvariable=self.top_n_pathways, width=6).pack(side="left", padx=(5, 0))
        tk.Label(topn_row, text="(for heatmaps)", bg="#f0f0f0", font=("Arial", 8), fg="#666").pack(side="left", padx=(5, 0))

        # Pathway filter dropdown
        filter_row = tk.Frame(settings_frame, bg="#f0f0f0")
        filter_row.pack(fill="x", padx=10, pady=5)
        tk.Label(filter_row, text="Pathway filter:", bg="#f0f0f0", font=("Arial", 9)).pack(side="left")
        self.filter_dropdown = ttk.Combobox(
            filter_row,
            textvariable=self.pathway_filter,
            values=["All Pathways", "Common Pathways Only"],
            state="readonly",
            width=18,
        )
        self.filter_dropdown.pack(side="left", padx=(5, 0))

        # Figure size controls
        size_row = tk.Frame(settings_frame, bg="#f0f0f0")
        size_row.pack(fill="x", padx=10, pady=5)
        tk.Label(size_row, text="Figure Size:", bg="#f0f0f0", font=("Arial", 9)).pack(side="left")
        tk.Label(size_row, text="W:", bg="#f0f0f0", font=("Arial", 8)).pack(side="left", padx=(10, 0))
        tk.Entry(size_row, textvariable=self.fig_width, width=4).pack(side="left", padx=(2, 5))
        tk.Label(size_row, text="H:", bg="#f0f0f0", font=("Arial", 8)).pack(side="left")
        tk.Entry(size_row, textvariable=self.fig_height, width=4).pack(side="left", padx=(2, 0))

        # === Action Buttons ===
        action_frame = tk.Frame(left_content, bg="#f0f0f0")
        action_frame.pack(fill="x", padx=10, pady=10)

        self.run_button = tk.Button(
            action_frame,
            text="⚡ Run Comparative Analysis",
            command=self._run_analysis,
            bg="#4caf50",
            fg="white",
            font=("Arial", 11, "bold"),
            relief="raised",
            borderwidth=3,
            state="disabled",
            padx=20,
            pady=8,
        )
        self.run_button.pack(fill="x")

        # === Export Section ===
        export_frame = tk.LabelFrame(
            left_content,
            text="💾 Export Options",
            font=("Arial", 10, "bold"),
            bg="#f0f0f0",
            fg="#2c3e50",
        )
        export_frame.pack(fill="x", padx=10, pady=10)

        export_btn_row = tk.Frame(export_frame, bg="#f0f0f0")
        export_btn_row.pack(fill="x", padx=10, pady=10)

        self.export_excel_btn = tk.Button(
            export_btn_row,
            text="📊 Export Excel",
            command=self._export_summary,
            bg="#9c27b0",
            fg="white",
            font=("Arial", 9, "bold"),
            relief="raised",
            state="disabled",
            padx=10,
            pady=4,
        )
        self.export_excel_btn.pack(side="left", padx=(0, 5))

        self.export_png_btn = tk.Button(
            export_btn_row,
            text="🖼️ Save PNG",
            command=lambda: self._export_figure("png"),
            bg="#00bcd4",
            fg="white",
            font=("Arial", 9, "bold"),
            relief="raised",
            state="disabled",
            padx=10,
            pady=4,
        )
        self.export_png_btn.pack(side="left", padx=(0, 5))

        self.export_svg_btn = tk.Button(
            export_btn_row,
            text="📐 Save SVG",
            command=lambda: self._export_figure("svg"),
            bg="#795548",
            fg="white",
            font=("Arial", 9, "bold"),
            relief="raised",
            state="disabled",
            padx=10,
            pady=4,
        )
        self.export_svg_btn.pack(side="left")

        # ---------------- Right column ----------------
        right_frame = tk.LabelFrame(
            columns_frame,
            text="📋 Analysis Log & Visualization",
            font=("Arial", 11, "bold"),
            bg="#f0f0f0",
            fg="#4caf50",
        )
        right_frame.grid(row=0, column=1, sticky="nsew", padx=(5, 0))
        right_frame.grid_rowconfigure(0, weight=0)
        right_frame.grid_rowconfigure(1, weight=1)
        right_frame.grid_columnconfigure(0, weight=1)

        self.log_text = scrolledtext.ScrolledText(
            right_frame,
            wrap=tk.WORD,
            font=("Courier New", 9),
            bg="#ffffff",
            fg="#000000",
            state="disabled",
            height=6,
        )
        self.log_text.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=10, pady=(10, 5))

        # Plot area with scrollbars
        self.plot_container = tk.Frame(right_frame, bg="#ffffff")
        self.plot_container.grid(row=1, column=0, sticky="nsew", padx=(10, 0), pady=(0, 10))

        self.plot_canvas = tk.Canvas(self.plot_container, bg="#ffffff", highlightthickness=1, highlightbackground="#cccccc")
        
        self.plot_scrollbar_y = ttk.Scrollbar(right_frame, orient="vertical", command=self.plot_canvas.yview)
        self.plot_scrollbar_x = ttk.Scrollbar(self.plot_container, orient="horizontal", command=self.plot_canvas.xview)
        
        self.plot_canvas.configure(yscrollcommand=self.plot_scrollbar_y.set, xscrollcommand=self.plot_scrollbar_x.set)
        
        self.plot_scrollbar_x.pack(side="bottom", fill="x")
        self.plot_canvas.pack(side="left", fill="both", expand=True)
        self.plot_scrollbar_y.grid(row=1, column=1, sticky="ns", pady=(0, 10))

        self.plot_inner = tk.Frame(self.plot_canvas, bg="#ffffff")
        self._plot_window = self.plot_canvas.create_window((0, 0), window=self.plot_inner, anchor="nw")

        def _update_plot_scroll(event):
            self.plot_canvas.configure(scrollregion=self.plot_canvas.bbox("all"))

        self.plot_inner.bind("<Configure>", _update_plot_scroll)

        def _on_plot_mousewheel(event):
            self.plot_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        self.plot_canvas.bind("<MouseWheel>", _on_plot_mousewheel)
        self.plot_inner.bind("<MouseWheel>", _on_plot_mousewheel)

    # ------------------------------------------------------------------
    # Logging helper
    # ------------------------------------------------------------------
    def _log(self, message: str):
        timestamp = pd.Timestamp.now().strftime("%H:%M:%S")
        self.log_text.configure(state="normal")
        self.log_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.log_text.see(tk.END)
        self.log_text.configure(state="disabled")

    # ------------------------------------------------------------------
    # File handling
    # ------------------------------------------------------------------
    def _upload_files(self):
        file_paths = filedialog.askopenfilenames(
            title="Select comparative pathway exports",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")],
        )
        if not file_paths:
            return

        new_files = 0
        for path in file_paths:
            name = os.path.basename(path)
            if any(item["path"] == path for item in self.uploaded_files):
                self._log(f"⚠️ {name} already loaded; skipping")
                continue

            try:
                df = self._load_pathway_sheet(path)
                if df is None or df.empty:
                    self._log(f"❌ {name}: could not find a pathway sheet with data")
                    continue

                label = self._prompt_dataset_label(name)
                if label is None:
                    self._log(f"⚠️ Skipped {name} (no dataset name provided)")
                    continue

                record = {"name": label, "original_name": name, "path": path, "df": df}
                self.uploaded_files.append(record)
                self.files_tree.insert("", "end", values=(label, len(df), len(df.columns)))
                new_files += 1
                self._log(f"✅ Loaded {label} ({name})")
            except Exception as exc:
                self._log(f"❌ Failed to read {name}: {exc}")
                logger.error("Comparative upload error", exc_info=True)

        self._update_button_states()

    def _prompt_dataset_label(self, default_name: str) -> Optional[str]:
        existing = {item['name'] for item in self.uploaded_files}
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

            if label in existing:
                messagebox.showwarning(
                    "Duplicate Name",
                    "A dataset with this name already exists. Please choose another."
                )
                continue

            return label

    def _remove_selected(self):
        """Remove selected dataset from the list"""
        selected = self.files_tree.selection()
        if not selected:
            messagebox.showwarning("No Selection", "Please select a dataset to remove.")
            return

        item = selected[0]
        values = self.files_tree.item(item, 'values')
        dataset_name = values[0]

        self.uploaded_files = [f for f in self.uploaded_files if f['name'] != dataset_name]
        self.files_tree.delete(item)
        self._log(f"➖ Removed: {dataset_name}")
        self._update_button_states()

    def _rename_selected(self):
        """Rename selected dataset"""
        selected = self.files_tree.selection()
        if not selected:
            messagebox.showwarning("No Selection", "Please select a dataset to rename.")
            return

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

    def _clear_files(self):
        if not self.uploaded_files:
            return
        if not messagebox.askyesno("Clear", "Remove all uploaded comparative datasets?"):
            return
        self.uploaded_files.clear()
        self.files_tree.delete(*self.files_tree.get_children())
        self.results_table = None
        self.heatmap_matrix = None
        self._clear_plot()
        self._update_button_states()
        self._log("🗑️ Cleared all datasets")

    def _update_button_states(self):
        """Update button states based on current data"""
        has_files = len(self.uploaded_files) > 0
        has_multiple = len(self.uploaded_files) >= 2
        has_results = self.results_table is not None

        self.clear_button.config(state="normal" if has_files else "disabled")
        self.remove_button.config(state="normal" if has_files else "disabled")
        self.rename_button.config(state="normal" if has_files else "disabled")
        self.run_button.config(state="normal" if has_multiple else "disabled")
        self.export_excel_btn.config(state="normal" if has_results else "disabled")
        self.export_png_btn.config(state="normal" if self._current_figure else "disabled")
        self.export_svg_btn.config(state="normal" if self._current_figure else "disabled")

    def _load_pathway_sheet(self, path: str) -> Optional[pd.DataFrame]:
        """Load the Selected_Pathways sheet (or the first sheet) from a workbook."""
        xl = pd.ExcelFile(path)
        sheet_name = "Selected_Pathways" if "Selected_Pathways" in xl.sheet_names else xl.sheet_names[0]
        df = xl.parse(sheet_name)
        df.columns = [str(col).strip() for col in df.columns]
        return df

    # ------------------------------------------------------------------
    # Analysis pipeline
    # ------------------------------------------------------------------
    def _run_analysis(self):
        if len(self.uploaded_files) < 2:
            messagebox.showerror("Need more data", "Upload at least two pathway exports to compare.")
            return

        try:
            p_thr = float(self.pvalue_threshold.get())
            if not 0 < p_thr < 1:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid threshold", "Enter a valid p-value threshold between 0 and 1.")
            return

        try:
            self._log("\n🔄 Running comparative analysis...")
            combined = None
            datasets = []

            for info in self.uploaded_files:
                df = info["df"].copy()
                path_col = self._detect_column(df, ["pathway", "name"], fallback=df.columns[0])
                if path_col is None:
                    self._log(f"⚠️ {info['name']} has no pathway-like column; skipping")
                    continue

                df[path_col] = df[path_col].astype(str).str.strip()
                df = df[df[path_col].astype(bool)]
                df.rename(columns={path_col: "Pathway"}, inplace=True)

                p_col = self._detect_column(df, ["combined_pvalue", "pvalue", "p_value", "p-value"])
                z_col = self._detect_column(df, ["z_score", "zscore", "z-score"], optional=True)
                fc_col = self._detect_column(df, ["mean_log2fc", "log2fc"], optional=True)
                count_col = self._detect_column(df, ["n_metabolites", "k_hits", "count", "compoundn"], optional=True)
                status_col = self._detect_column(df, ["status"], optional=True)

                subset = df[["Pathway"]].copy()
                if p_col:
                    subset[f"{info['name']}__pvalue"] = pd.to_numeric(df[p_col], errors="coerce")
                if z_col:
                    subset[f"{info['name']}__zscore"] = pd.to_numeric(df[z_col], errors="coerce")
                if fc_col:
                    subset[f"{info['name']}__mean_log2fc"] = pd.to_numeric(df[fc_col], errors="coerce")
                if count_col:
                    subset[f"{info['name']}__count"] = pd.to_numeric(df[count_col], errors="coerce")
                if status_col:
                    subset[f"{info['name']}__status"] = df[status_col].astype(str)

                subset.set_index("Pathway", inplace=True)
                datasets.append((info["name"], subset))

                if combined is None:
                    combined = subset
                else:
                    combined = combined.join(subset, how="outer")

            if not datasets or combined is None or combined.empty:
                self._log("❌ No usable pathway tables were found across the uploaded files.")
                return

            self.results_table = combined
            self._log(f"✅ Aggregated {combined.shape[0]} pathways across {len(datasets)} datasets")

            # Count pathways per dataset and overlap
            self._analyze_overlap(combined, datasets, p_thr)

            # Build and render visualization based on selected type
            viz_type = self.viz_type.get()
            self._render_visualization(datasets, combined, p_thr, viz_type)

            self._update_button_states()

        except Exception as exc:
            logger.error("Comparative analysis failed", exc_info=True)
            messagebox.showerror("Error", f"Comparative analysis failed: {exc}")
            self._log(f"❌ Error: {exc}")

    def _on_viz_changed(self, event=None):
        """Handle visualization dropdown changes"""
        if not self._all_figures:
            return  # No analysis run yet
        
        viz_type = self.viz_type.get()
        self._switch_to_visualization(viz_type)
    
    def _switch_to_visualization(self, viz_type):
        """Switch to display a specific visualization from stored figures"""
        if viz_type not in self._all_figures:
            self._log(f"⚠️ {viz_type} not available")
            return
        
        self._clear_plot()
        self._current_figure = self._all_figures[viz_type]
        if viz_type in self._all_matrices:
            self.heatmap_matrix = self._all_matrices[viz_type]
        
        self._heatmap_canvas = FigureCanvasTkAgg(self._current_figure, master=self.plot_inner)
        self._heatmap_canvas.draw()
        self._heatmap_canvas.get_tk_widget().pack(fill="both", expand=True)
        self._log(f"✅ Displaying {viz_type}")
        self._update_button_states()

    def _analyze_overlap(self, combined, datasets, p_thr):
        """Analyze and log pathway overlap statistics"""
        # Count significant pathways per dataset
        self._log("\n📊 Pathway Statistics:")
        for name, _ in datasets:
            p_col = f"{name}__pvalue"
            if p_col in combined.columns:
                sig_count = (combined[p_col] <= p_thr).sum()
                total = combined[p_col].notna().sum()
                self._log(f"   • {name}: {sig_count}/{total} significant (p ≤ {p_thr})")

        # Find common pathways (present in all datasets)
        common_pathways = set(combined.index)
        for name, _ in datasets:
            p_col = f"{name}__pvalue"
            if p_col in combined.columns:
                present = set(combined[combined[p_col].notna()].index)
                common_pathways &= present

        self._log(f"\n🔗 Pathways present in ALL datasets: {len(common_pathways)}")

        # Find commonly significant
        sig_in_all = None
        for name, _ in datasets:
            p_col = f"{name}__pvalue"
            if p_col in combined.columns:
                sig_mask = combined[p_col] <= p_thr
                sig_in_all = sig_mask if sig_in_all is None else (sig_in_all & sig_mask)

        if sig_in_all is not None:
            common_sig = sig_in_all.sum()
            self._log(f"🔗 Pathways significant in ALL datasets: {common_sig}")

    def _detect_column(self, df: pd.DataFrame, keywords, fallback: Optional[str] = None, optional: bool = False):
        for col in df.columns:
            normalized = col.lower().replace(" ", "").replace("-", "").replace("_", "")
            for kw in keywords:
                kw_normalized = kw.replace("_", "").replace("-", "")
                if kw_normalized in normalized:
                    return col
        if fallback and fallback in df.columns and not optional:
            return fallback
        return None

    def _render_visualization(self, datasets, combined, p_thr, viz_type):
        """Generate all visualizations and display the selected one"""
        self._clear_plot()
        self._all_figures.clear()
        self._all_matrices.clear()

        try:
            width = float(self.fig_width.get())
            height = float(self.fig_height.get())
        except ValueError:
            width, height = 10, 8

        try:
            top_n = int(self.top_n_pathways.get())
            if top_n < 1:
                top_n = 50
        except ValueError:
            top_n = 50

        common_only = self.pathway_filter.get() == "Common Pathways Only"
        
        # Generate ALL visualizations
        self._log("\n🎨 Generating all visualizations...")
        self._render_zscore_heatmap(datasets, combined, width, height, top_n, common_only, store=True)
        self._render_pvalue_heatmap(datasets, combined, width, height, store=True)
        self._render_overlap_bar(datasets, combined, p_thr, width, height, store=True)
        self._render_venn_diagram(datasets, combined, p_thr, width, height, store=True)
        
        # Display the selected one
        self._switch_to_visualization(viz_type)

    def _render_zscore_heatmap(self, datasets, combined, width, height, top_n=50, common_only=False, store=False):
        """Render Z-Score heatmap with blue-white-red color scale"""
        # Build Z-score matrix
        zscore_matrix = pd.DataFrame(index=combined.index)
        has_data = False

        for name, _ in datasets:
            z_col = f"{name}__zscore"
            if z_col in combined.columns:
                zscore_matrix[name] = combined[z_col]
                has_data = True
            else:
                zscore_matrix[name] = np.nan

        if not has_data:
            if not store:
                self._log("⚠️ No Z-Score data available. Falling back to p-value heatmap.")
                self._render_pvalue_heatmap(datasets, combined, width, height, store=False)
            return

        # Filter to rows with at least one non-NaN value
        zscore_matrix = zscore_matrix.dropna(how='all')
        if zscore_matrix.empty:
            self._log("⚠️ No valid Z-Score data to display.")
            return

        # Identify common pathways (present in ALL datasets)
        common_mask = zscore_matrix.notna().all(axis=1)
        common_pathways = zscore_matrix[common_mask].copy()
        unique_pathways = zscore_matrix[~common_mask].copy()

        # Apply common pathways filter if selected
        if common_only:
            if common_pathways.empty:
                self._log("⚠️ No pathways found in all datasets.")
                return
            zscore_matrix = common_pathways
            self._log(f"📊 Filtering to {len(zscore_matrix)} common pathways")
        else:
            # Sort common pathways by max absolute Z-score (from any dataset)
            if not common_pathways.empty:
                common_pathways['_sort'] = common_pathways.abs().max(axis=1)
                common_pathways = common_pathways.sort_values('_sort', ascending=False).drop('_sort', axis=1)
            
            # Sort unique pathways by max absolute Z-score
            if not unique_pathways.empty:
                unique_pathways['_sort'] = unique_pathways.abs().max(axis=1)
                unique_pathways = unique_pathways.sort_values('_sort', ascending=False).drop('_sort', axis=1)
            
            # Concatenate: common pathways first, then unique
            zscore_matrix = pd.concat([common_pathways, unique_pathways])
            self._log(f"📊 {len(common_pathways)} common + {len(unique_pathways)} unique pathways")

        # Limit to top pathways for readability
        display_n = min(top_n, len(zscore_matrix))
        display_matrix = zscore_matrix.head(display_n)

        self.heatmap_matrix = display_matrix

        # Create custom blue-gray-red colormap (0=gray, not white)
        colors = ['#2166ac', '#67a9cf', '#d1e5f0', '#d0d0d0', '#fddbc7', '#ef8a62', '#b2182b']
        cmap = LinearSegmentedColormap.from_list('zscore_cmap', colors, N=256)
        cmap.set_bad(color='white')  # NaN (absent pathways) = white

        # Determine symmetric color limits
        vmax = max(abs(display_matrix.min().min()), abs(display_matrix.max().max()))
        vmax = max(vmax, 0.1)  # Ensure some range

        fig = Figure(figsize=(width, height))
        ax = fig.add_subplot(111)

        # Convert to masked array to properly handle NaN
        masked_data = np.ma.masked_invalid(display_matrix.values)
        im = ax.imshow(masked_data, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)

        ax.set_xticks(range(len(display_matrix.columns)))
        ax.set_xticklabels(display_matrix.columns, rotation=45, ha="right", fontsize=10, fontweight='bold')
        ax.set_yticks(range(len(display_matrix.index)))
        ax.set_yticklabels(display_matrix.index, fontsize=9, fontweight='bold')
        ax.set_title("Z-Score Heatmap: Pathway Activation Comparison", fontsize=12, fontweight='bold')
        ax.tick_params(axis='both', which='major', labelsize=9)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Z-Score (Blue=Inhibited, Gray=0, Red=Activated)", fontsize=10, fontweight='bold')
        cbar.ax.tick_params(labelsize=9)

        fig.tight_layout()
        if store:
            self._all_figures["Z-Score Heatmap"] = fig
            self._all_matrices["Z-Score Heatmap"] = display_matrix
        else:
            self._display_figure(fig)
            self._log(f"✅ Rendered Z-Score heatmap ({display_n} pathways)")

    def _render_pvalue_heatmap(self, datasets, combined, width, height, store=False):
        """Render -log10(p) heatmap"""
        heatmap = pd.DataFrame(index=combined.index)
        for name, _ in datasets:
            p_col = f"{name}__pvalue"
            if p_col in combined.columns:
                series = combined[p_col].astype(float)
                # Keep NaN for absent pathways, compute -log10 for present ones
                heatmap[name] = series.apply(lambda x: -np.log10(x) if pd.notna(x) and x > 0 else np.nan)
            else:
                heatmap[name] = np.nan

        heatmap = heatmap.replace([np.inf, -np.inf], np.nan).dropna(how='all')
        
        top_n = min(50, len(heatmap))
        heatmap = heatmap.sort_values(by=list(heatmap.columns), ascending=False, na_position='last').head(top_n)
        self.heatmap_matrix = heatmap

        fig = Figure(figsize=(width, height))
        ax = fig.add_subplot(111)
        
        # Create colormap with white for NaN (absent pathways)
        from matplotlib import cm
        cmap = cm.get_cmap("YlOrRd").copy()
        cmap.set_bad(color='white')
        
        # Convert to masked array to handle NaN
        masked_data = np.ma.masked_invalid(heatmap.values)
        im = ax.imshow(masked_data, aspect="auto", cmap=cmap)

        ax.set_xticks(range(len(heatmap.columns)))
        ax.set_xticklabels(heatmap.columns, rotation=45, ha="right", fontsize=10, fontweight='bold')
        ax.set_yticks(range(len(heatmap.index)))
        ax.set_yticklabels(heatmap.index, fontsize=9, fontweight='bold')
        ax.set_title("-log10(p-value) Heatmap", fontsize=12, fontweight='bold')
        ax.tick_params(axis='both', which='major', labelsize=9)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("-log10(p)", fontsize=10, fontweight='bold')
        cbar.ax.tick_params(labelsize=9)
        fig.tight_layout()
        if store:
            self._all_figures["P-Value Heatmap"] = fig
            self._all_matrices["P-Value Heatmap"] = heatmap
        else:
            self._display_figure(fig)
            self._log("✅ Rendered P-value heatmap")

    def _render_overlap_bar(self, datasets, combined, p_thr, width, height, store=False):
        """Render bar chart showing pathway overlap"""
        dataset_names = [name for name, _ in datasets]
        
        # Calculate counts
        unique_counts = []
        shared_counts = []
        
        for name, _ in datasets:
            p_col = f"{name}__pvalue"
            if p_col in combined.columns:
                sig_in_this = combined[p_col] <= p_thr
                
                # Check if significant in other datasets
                sig_in_others = pd.Series(False, index=combined.index)
                for other_name, _ in datasets:
                    if other_name != name:
                        other_p_col = f"{other_name}__pvalue"
                        if other_p_col in combined.columns:
                            sig_in_others |= (combined[other_p_col] <= p_thr)
                
                unique = (sig_in_this & ~sig_in_others).sum()
                shared = (sig_in_this & sig_in_others).sum()
                unique_counts.append(unique)
                shared_counts.append(shared)
            else:
                unique_counts.append(0)
                shared_counts.append(0)

        fig = Figure(figsize=(width, height))
        ax = fig.add_subplot(111)

        x = np.arange(len(dataset_names))
        bar_width = 0.35

        bars1 = ax.bar(x - bar_width/2, unique_counts, bar_width, label='Unique', color='#2196f3')
        bars2 = ax.bar(x + bar_width/2, shared_counts, bar_width, label='Shared with others', color='#ff9800')

        ax.set_ylabel('Number of Significant Pathways', fontsize=10, fontweight='bold')
        ax.set_title(f'Pathway Overlap Analysis (p ≤ {p_thr})', fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(dataset_names, rotation=45, ha='right', fontsize=10, fontweight='bold')
        ax.tick_params(axis='both', which='major', labelsize=9)
        ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0, fontsize=9)

        # Add value labels
        for bar in bars1:
            height = bar.get_height()
            ax.annotate(f'{int(height)}', xy=(bar.get_x() + bar.get_width()/2, height),
                       xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9, fontweight='bold')
        for bar in bars2:
            height = bar.get_height()
            ax.annotate(f'{int(height)}', xy=(bar.get_x() + bar.get_width()/2, height),
                       xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9, fontweight='bold')

        fig.tight_layout(rect=[0, 0, 0.85, 1])
        if store:
            self._all_figures["Pathway Overlap Bar"] = fig
        else:
            self._display_figure(fig)
            self._log("✅ Rendered pathway overlap bar chart")

    def _render_status_comparison(self, datasets, combined, width, height, top_n=50, common_only=False):
        """Render comparison of pathway activation status"""
        # Build status matrix
        status_matrix = pd.DataFrame(index=combined.index)
        has_status = False

        for name, _ in datasets:
            status_col = f"{name}__status"
            if status_col in combined.columns:
                # Convert to numeric: Activated=1, Inhibited=-1, else=0
                status_series = combined[status_col].str.lower()
                numeric_status = status_series.map({
                    'activated': 1, 'active': 1, 'up': 1,
                    'inhibited': -1, 'inactive': -1, 'down': -1
                }).fillna(0)
                status_matrix[name] = numeric_status
                has_status = True
            else:
                # Try to infer from Z-score
                z_col = f"{name}__zscore"
                if z_col in combined.columns:
                    z_vals = combined[z_col]
                    status_matrix[name] = np.sign(z_vals).fillna(0)
                    has_status = True

        if not has_status:
            self._log("⚠️ No status data available.")
            return

        status_matrix = status_matrix.dropna(how='all')
        if status_matrix.empty:
            self._log("⚠️ No valid status data to display.")
            return

        # Apply common pathways filter if selected
        if common_only:
            status_matrix = status_matrix.dropna(how='any')
            if status_matrix.empty:
                self._log("⚠️ No pathways found in all datasets.")
                return
            self._log(f"📊 Filtering to {len(status_matrix)} common pathways")

        # Filter to pathways with any non-zero status
        status_matrix = status_matrix[(status_matrix != 0).any(axis=1)]
        display_n = min(top_n, len(status_matrix))
        display_matrix = status_matrix.head(display_n)

        # Custom colormap: blue (-1) -> white (0) -> red (1)
        colors = ['#2166ac', '#ffffff', '#b2182b']
        cmap = LinearSegmentedColormap.from_list('status_cmap', colors, N=3)

        fig = Figure(figsize=(width, height))
        ax = fig.add_subplot(111)

        im = ax.imshow(display_matrix.values, aspect="auto", cmap=cmap, vmin=-1, vmax=1)

        ax.set_xticks(range(len(display_matrix.columns)))
        ax.set_xticklabels(display_matrix.columns, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(len(display_matrix.index)))
        ax.set_yticklabels(display_matrix.index, fontsize=8)
        ax.set_title("Pathway Status Comparison", fontsize=11, fontweight='bold')

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, ticks=[-1, 0, 1])
        cbar.ax.set_yticklabels(['Inhibited', 'Neutral', 'Activated'])

        fig.tight_layout()
        self._display_figure(fig)
        self._log(f"✅ Rendered status comparison ({display_n} pathways)")

    def _render_venn_diagram(self, datasets, combined, p_thr, width, height, store=False):
        """Render Venn diagram showing pathway overlap between datasets"""
        if len(datasets) < 2 or len(datasets) > 3:
            self._log("⚠️ Venn diagram supports 2-3 datasets only.")
            if len(datasets) > 3:
                self._log("💡 Using first 3 datasets for Venn diagram.")
                datasets = datasets[:3]
            else:
                return

        # Get significant pathways for each dataset
        sig_sets = []
        set_labels = []
        
        for name, _ in datasets[:3]:
            p_col = f"{name}__pvalue"
            if p_col in combined.columns:
                sig_pathways = set(combined[combined[p_col] <= p_thr].index)
                sig_sets.append(sig_pathways)
                set_labels.append(name)
            else:
                sig_sets.append(set())
                set_labels.append(name)

        fig = Figure(figsize=(width, height))
        ax = fig.add_subplot(111)

        if len(sig_sets) == 2:
            # Calculate actual counts for equal-sized circles
            venn = venn2(sig_sets, set_labels=set_labels, ax=ax, set_colors=('#2196f3', '#ff9800'), alpha=0.6)
            # Color the overlap region
            if venn.get_patch_by_id('11'):
                venn.get_patch_by_id('11').set_color('#4caf50')
                venn.get_patch_by_id('11').set_alpha(0.6)
            # Make set labels bold
            # Make set labels bold
            for text in venn.set_labels:
                if text:
                    text.set_fontsize(10)
                    text.set_fontweight('bold')
            # Make subset labels bold
            for text in venn.subset_labels:
                if text:
                    text.set_fontsize(9)
                    text.set_fontweight('bold')
        else:  # 3 datasets
            venn = venn3(sig_sets, set_labels=set_labels, ax=ax, set_colors=('#2196f3', '#ff9800', '#9c27b0'), alpha=0.6)
            # Color the overlap regions
            colors = {'110': '#4caf50', '101': '#00bcd4', '011': '#ff5722', '111': '#607d8b'}
            for region_id, color in colors.items():
                patch = venn.get_patch_by_id(region_id)
                if patch:
                    patch.set_color(color)
                    patch.set_alpha(0.6)
            # Make set labels bold
            for text in venn.set_labels:
                if text:
                    text.set_fontsize(10)
                    text.set_fontweight('bold')
            # Make subset labels bold
            for text in venn.subset_labels:
                if text:
                    text.set_fontsize(9)
                    text.set_fontweight('bold')

        ax.set_title(f'Significant Pathway Overlap (p ≤ {p_thr})', fontsize=12, fontweight='bold')

        # Calculate and log overlap statistics
        if len(sig_sets) >= 2:
            overlap_2 = sig_sets[0] & sig_sets[1]
            self._log(f"🔗 Overlap ({set_labels[0].split(chr(10))[0]} ∩ {set_labels[1].split(chr(10))[0]}): {len(overlap_2)} pathways")
        if len(sig_sets) == 3:
            overlap_all = sig_sets[0] & sig_sets[1] & sig_sets[2]
            self._log(f"🔗 Overlap (all 3): {len(overlap_all)} pathways")

        fig.tight_layout()
        if store:
            self._all_figures["Venn Diagram"] = fig
        else:
            self._display_figure(fig)
            self._log("✅ Rendered Venn diagram")

    def _display_figure(self, fig: Figure):
        """Display a matplotlib figure in the plot area"""
        self._current_figure = fig
        self._heatmap_canvas = FigureCanvasTkAgg(fig, master=self.plot_inner)
        self._heatmap_canvas.draw()
        self._heatmap_canvas.get_tk_widget().pack(fill="both", expand=True)
        self._update_button_states()

    def _clear_plot(self):
        if self._heatmap_canvas is not None:
            try:
                self._heatmap_canvas.get_tk_widget().destroy()
            except Exception:
                pass
            self._heatmap_canvas = None
        self._current_figure = None
        if hasattr(self, "plot_canvas"):
            self.plot_canvas.yview_moveto(0)
            self.plot_canvas.xview_moveto(0)

    # ------------------------------------------------------------------
    # Export helpers
    # ------------------------------------------------------------------
    def _export_summary(self):
        if self.results_table is None or self.results_table.empty:
            messagebox.showinfo("No data", "Run the comparative analysis before exporting.")
            return

        path = filedialog.asksaveasfilename(
            title="Save comparative summary",
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")],
            initialfile="Comparative_Pathway_Summary.xlsx",
        )
        if not path:
            return

        try:
            with pd.ExcelWriter(path, engine="openpyxl") as writer:
                self.results_table.reset_index().to_excel(writer, sheet_name="Combined_Table", index=False)
                
                # Export all matrices from all visualizations
                for viz_name, matrix in self._all_matrices.items():
                    safe_sheet_name = viz_name.replace(" ", "_").replace("-", "")[:31]  # Excel sheet name limit
                    matrix.reset_index().to_excel(writer, sheet_name=safe_sheet_name, index=False)
                    
            self._log(f"💾 Saved summary with {len(self._all_matrices)} visualization sheets to {os.path.basename(path)}")
            messagebox.showinfo("Export complete", f"Comparative summary saved to:\n{path}")
        except Exception as exc:
            logger.error("Comparative export error", exc_info=True)
            messagebox.showerror("Export failed", str(exc))
            self._log(f"❌ Export failed: {exc}")

    def _export_figure(self, fmt: str):
        """Export the current figure as PNG or SVG"""
        if self._current_figure is None:
            messagebox.showinfo("No figure", "Run the analysis first to generate a visualization.")
            return

        ext = f".{fmt}"
        filetypes = [(f"{fmt.upper()} files", f"*{ext}"), ("All files", "*.*")]
        
        # Add plot name to filename
        viz_type = self.viz_type.get()
        safe_name = viz_type.replace(" ", "_").replace("-", "")

        path = filedialog.asksaveasfilename(
            title=f"Save figure as {fmt.upper()}",
            defaultextension=ext,
            filetypes=filetypes,
            initialfile=f"Comparative_{safe_name}{ext}",
        )
        if not path:
            return

        try:
            self._current_figure.savefig(path, format=fmt, dpi=300, bbox_inches='tight')
            self._log(f"💾 Saved figure to {os.path.basename(path)}")
            messagebox.showinfo("Export complete", f"Figure saved to:\n{path}")
        except Exception as exc:
            logger.error(f"Figure export error ({fmt})", exc_info=True)
            messagebox.showerror("Export failed", str(exc))
            self._log(f"❌ Export failed: {exc}")
