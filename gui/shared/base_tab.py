"""Base class for all GUI tabs - provides common interface"""
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from abc import ABC, abstractmethod
import logging
import json

from gui.shared.utils import resolve_runtime_config_path


logger = logging.getLogger(__name__)


# Static flag to ensure styles are configured only once
_styles_configured = False


def _setup_global_styles():
    """Configure all ttk styles globally (runs once)"""
    global _styles_configured
    if _styles_configured:
        return
    
    style = ttk.Style()
    
    # Sub-tabs styling for ID Annotation tab
    style.configure('IdAnnotation.TNotebook.Tab', 
                   font=('Arial', 10, 'bold'), 
                   padding=8)
    
    _styles_configured = True

class BaseTab(ABC):
    """
    Abstract base class for all tab implementations.
    
    Provides common interface and utilities for all tabs.
    Each tab should inherit from this class.
    
    Example:
        class MyTab(BaseTab):
            def __init__(self, parent, data_manager):
                super().__init__(parent, data_manager)
                self.setup_ui()
            
            def setup_ui(self):
                # Create UI elements here
                label = ttk.Label(self.frame, text="My Tab")
                label.pack()
    """
    
    def __init__(self, parent, data_manager):
        """
        Initialize base tab.
        
        Args:
            parent: Parent ttk.Notebook widget
            data_manager: Shared DataManager instance for data persistence between tabs
        """
        self.parent = parent
        self.data_manager = data_manager
        self.frame = ttk.Frame(parent)
        self._is_initialized = False
        
        # Store reference to this tab instance on the frame for inter-tab communication
        self.frame.tab_instance = self
        
        # CRITICAL: Make memory_store a shared reference to DataManager's memory_store
        # This ensures all tabs see the same data dictionary
        self.memory_store = self.data_manager.memory_store
        
        logger.debug(f"Initializing {self.__class__.__name__}")
    
    @abstractmethod
    def setup_ui(self):
        """
        Setup the UI for this tab.
        
        Subclasses MUST implement this method.
        Called during __init__ to create all widgets.
        """
        pass
    
    def get_frame(self):
        """
        Return the tab's frame for adding to notebook.
        
        Returns:
            ttk.Frame: The frame containing this tab's UI
        """
        return self.frame
    
    def mark_initialized(self):
        """Mark this tab as fully initialized"""
        self._is_initialized = True
    
    def is_initialized(self):
        """Check if tab is fully initialized"""
        return self._is_initialized
    
    # ========== Common utilities for all tabs ==========
    
    def show_info(self, title, message):
        """Show information message dialog"""
        from tkinter import messagebox
        messagebox.showinfo(title, message)
    
    def show_error(self, title, message):
        """Show error message dialog"""
        from tkinter import messagebox
        messagebox.showerror(title, message)
    
    def show_warning(self, title, message):
        """Show warning message dialog"""
        from tkinter import messagebox
        messagebox.showwarning(title, message)
    
    def ask_yes_no(self, title, message):
        """Ask yes/no question"""
        from tkinter import messagebox
        return messagebox.askyesno(title, message)
    
    def ask_file_open(self, **kwargs):
        """Open file dialog"""
        from tkinter import filedialog
        return filedialog.askopenfilename(**kwargs)
    
    def ask_file_save(self, **kwargs):
        """Save file dialog"""
        from tkinter import filedialog
        return filedialog.asksaveasfilename(**kwargs)
    
    def ask_directory(self, **kwargs):
        """Open directory dialog"""
        from tkinter import filedialog
        return filedialog.askdirectory(**kwargs)
    
    def create_labeled_entry(self, parent, label_text, **entry_kwargs):
        """
        Helper to create label + entry widget pair.
        
        Args:
            parent: Parent widget
            label_text: Text for label
            **entry_kwargs: Keyword arguments for ttk.Entry
        
        Returns:
            tuple: (ttk.Label, ttk.Entry)
        """
        frame = ttk.Frame(parent)
        frame.pack(fill=tk.X, padx=5, pady=5)
        
        label = ttk.Label(frame, text=label_text, width=20)
        label.pack(side=tk.LEFT, padx=5)
        
        entry = ttk.Entry(frame, **entry_kwargs)
        entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        
        return label, entry
    
    def create_labeled_combobox(self, parent, label_text, values=None, **combo_kwargs):
        """
        Helper to create label + combobox widget pair.
        
        Args:
            parent: Parent widget
            label_text: Text for label
            values: List of values for combobox
            **combo_kwargs: Keyword arguments for ttk.Combobox
        
        Returns:
            tuple: (ttk.Label, ttk.Combobox)
        """
        if values is None:
            values = []
        
        frame = ttk.Frame(parent)
        frame.pack(fill=tk.X, padx=5, pady=5)
        
        label = ttk.Label(frame, text=label_text, width=20)
        label.pack(side=tk.LEFT, padx=5)
        
        combo = ttk.Combobox(frame, values=values, **combo_kwargs)
        combo.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)
        
        return label, combo
    
    def create_button(self, parent, text, command, **kwargs):
        """
        Helper to create a button.
        
        Args:
            parent: Parent widget
            text: Button text
            command: Callback function
            **kwargs: Keyword arguments for ttk.Button
        
        Returns:
            ttk.Button: The created button
        """
        button = ttk.Button(parent, text=text, command=command, **kwargs)
        return button
    
    def _create_tooltip(self, widget, text):
        """
        Create a tooltip for a widget.
        
        Args:
            widget: The tkinter widget to attach tooltip to
            text: Tooltip text to display
        """
        if not text or text.strip() == '':
            return
        
        def on_enter(event):
            tooltip = tk.Toplevel()
            tooltip.wm_overrideredirect(True)
            tooltip.wm_geometry(f"+{event.x_root + 10}+{event.y_root + 10}")
            label = tk.Label(
                tooltip,
                text=text,
                background="#ffffe0",
                relief=tk.SOLID,
                borderwidth=1,
                font=("Helvetica", 9),
                wraplength=300
            )
            label.pack()
            
            widget.tooltip = tooltip
        
        def on_leave(event):
            if hasattr(widget, 'tooltip'):
                widget.tooltip.destroy()
                delattr(widget, 'tooltip')
        
        widget.bind("<Enter>", on_enter)
        widget.bind("<Leave>", on_leave)

    def _get_workers_count(self, var: tk.Variable, default: int = 1) -> int:
        """Safely parse a tk.StringVar (or similar) storing worker count.

        Returns a sensible int fallback when the UI control is empty or invalid.
        """
        try:
            if var is None:
                return max(1, default)
            # Some controls may be StringVar or other variable-like
            raw = None
            try:
                raw = var.get()
            except Exception:
                raw = str(var)
            if raw is None:
                return max(1, default)
            raw = str(raw).strip()
            if raw == "":
                # if empty, fall back to a reasonable default (min of 6 or cpu count)
                import multiprocessing
                cpu_count = multiprocessing.cpu_count()
                return max(1, min(6, cpu_count))
            val = int(raw)
            if val < 1:
                return 1
            return val
        except Exception:
            return max(1, default)

   # --- Global settings gather/apply ---
    def _settings_file(self):
        return resolve_runtime_config_path('gui_saved_settings.json')

    def _gather_all_settings(self):
        """Collect all relevant GUI settings into a JSON-serializable dict."""
        settings = {}
        try:
            # Visualization params (will be updated from GUI values)
            self._update_viz_params()
            settings['viz_params'] = {}
            for k, v in self.viz_params.items():
                if v is not None and hasattr(v, '__dict__'):
                    settings['viz_params'][k] = v.__dict__.copy()
                else:
                    settings['viz_params'][k] = None
            
            # Viz groups and mapping
            settings['viz_group_definitions'] = getattr(self, 'viz_group_definitions', {})
            settings['viz_group_mapping'] = getattr(self, 'viz_group_mapping', {})
            
            # Handle color map conversion - only if matplotlib is available
            try:
                from matplotlib import colors as mcolors
                settings['viz_color_map'] = {g: mcolors.to_hex(c) if not isinstance(c, str) else c for g, c in getattr(self, 'viz_color_map', {}).items()}
            except ImportError:
                # Fallback if matplotlib not available
                viz_cmap = getattr(self, 'viz_color_map', {})
                settings['viz_color_map'] = {g: c for g, c in viz_cmap.items()}
            
            settings['viz_preferred_group_order'] = self.viz_preferred_group_order.get() if hasattr(self, 'viz_preferred_group_order') else ''
            
            # Visualization selection checkboxes
            settings['viz_selected'] = {}
            if hasattr(self, 'viz_selected'):
                for k, v in self.viz_selected.items():
                    try:
                        settings['viz_selected'][k] = v.get()
                    except Exception:
                        settings['viz_selected'][k] = False
            
            # Lipid class visualization options
            settings['viz_include_lipid_class'] = {}
            if hasattr(self, 'viz_include_lipid_class'):
                for k, v in self.viz_include_lipid_class.items():
                    try:
                        settings['viz_include_lipid_class'][k] = v.get()
                    except Exception:
                        settings['viz_include_lipid_class'][k] = False
            
            # Gather all GUI control variables for reproducibility
            gui_vars = {}
            
            # Helper function to safely get var value
            def get_var_value(var_name):
                if hasattr(self, var_name):
                    try:
                        var = getattr(self, var_name)
                        if hasattr(var, 'get'):
                            return var.get()
                    except Exception:
                        pass
                return None
            
            # PCA variables
            pca_vars = [
                'pca_components', 'pca_3d', 'pca_interactive', 'pca_scree', 'pca_loadings',
                'pca_loadings_k', 'pca_fig_width', 'pca_fig_height', 'pca_fig_dpi',
                'pca_point_size_2d', 'pca_point_size_3d', 'pca_view_azim', 'pca_view_elev',
                'pca_save_2d', 'pca_save_3d', 'pca_save_excel', 'pca_xlabel_fontsize', 
                'pca_ylabel_fontsize', 'pca_title_fontsize', 'pca_tick_fontsize', 
                'pca_legend_fontsize'
            ]
            for var in pca_vars:
                gui_vars[var] = get_var_value(var)
            
            # Volcano variables
            volcano_vars = [
                'volcano_p_thresh', 'volcano_fc_thresh', 'volcano_skip_fc', 'volcano_annotate',
                'volcano_top_n', 'volcano_fig_width', 'volcano_fig_height', 'volcano_fig_dpi',
                'volcano_point_size_sig', 'volcano_point_size_nonsig', 'volcano_xlabel_fontsize',
                'volcano_ylabel_fontsize', 'volcano_title_fontsize', 'volcano_tick_fontsize',
                'volcano_count_fontsize', 'volcano_total_fontsize', 'volcano_legend_fontsize',
                'volcano_count_background', 'volcano_save_excel'
            ]
            for var in volcano_vars:
                gui_vars[var] = get_var_value(var)
            
            # Boxplot variables
            boxplot_vars = [
                'boxplot_top_n', 'boxplot_no_limit', 'boxplot_annotate', 'boxplot_p_thresh',
                'boxplot_fc_thresh', 'boxplot_skip_fc', 'boxplot_filter_mode', 'boxplot_filter_pairs',
                'boxplot_use_custom_only', 'boxplot_custom_list', 'boxplot_fig_width',
                'boxplot_fig_height', 'boxplot_fig_dpi', 'boxplot_xlabel_fontsize',
                'boxplot_ylabel_fontsize', 'boxplot_title_fontsize', 'boxplot_tick_fontsize',
                'boxplot_rotate_xticks', 'boxplot_xtick_rotation', 'boxplot_save_excel'
            ]
            for var in boxplot_vars:
                gui_vars[var] = get_var_value(var)
            
            # Heatmap variables
            heatmap_vars = [
                'heatmap_max', 'heatmap_show_fc_divider', 'heatmap_combined', 'heatmap_p_thresh',
                'heatmap_fc_thresh', 'heatmap_skip_fc', 'heatmap_filter_mode', 'heatmap_filter_pairs',
                'heatmap_use_custom_only', 'heatmap_auto_scale', 'heatmap_vmin', 'heatmap_vmax',
                'heatmap_combined_mode', 'heatmap_auto_size', 'heatmap_fig_width', 'heatmap_fig_height',
                'heatmap_fig_dpi', 'heatmap_custom_list', 'heatmap_feature_fontsize',
                'heatmap_sample_fontsize', 'heatmap_save_excel', 'heatmap_skip_unlisted', 'heatmap_cluster'
            ]
            for var in heatmap_vars:
                gui_vars[var] = get_var_value(var)
            
            # ROC variables
            roc_vars = [
                'roc_all_pairs', 'roc_max_metabolites', 'roc_min_auc', 'roc_fig_width',
                'roc_fig_height', 'roc_fig_dpi', 'roc_p_thresh', 'roc_fc_thresh',
                'roc_skip_fc', 'roc_filter_mode', 'roc_filter_pairs', 'roc_use_custom_only',
                'roc_custom_list', 'roc_xlabel_fontsize', 'roc_ylabel_fontsize',
                'roc_title_fontsize', 'roc_tick_fontsize', 'roc_legend_fontsize',
                'roc_save_excel', 'roc_skip_unlisted'
            ]
            for var in roc_vars:
                gui_vars[var] = get_var_value(var)
            
            settings['gui_vars'] = gui_vars
            
            # Log detailed summary of what's being saved
            saved_items = []
            saved_items.append(f"  ✓ Visualization parameters: {len([k for k in settings.get('viz_params', {}) if settings['viz_params'][k]])} types")
            if settings.get('viz_group_definitions'):
                saved_items.append(f"  ✓ Group definitions: {len(settings['viz_group_definitions'])} groups ({', '.join(settings['viz_group_definitions'].values())})")
            if settings.get('viz_group_mapping'):
                saved_items.append(f"  ✓ Sample-to-group mapping: {len(settings['viz_group_mapping'])} samples")
            if settings.get('viz_color_map'):
                saved_items.append(f"  ✓ Color map: {len(settings['viz_color_map'])} colors (with custom hex codes)")
            if settings.get('viz_preferred_group_order'):
                saved_items.append(f"  ✓ Preferred group order: {settings['viz_preferred_group_order']}")
            enabled_viz = [k for k, v in settings.get('viz_selected', {}).items() if v]
            if enabled_viz:
                saved_items.append(f"  ✓ Enabled visualizations: {', '.join(enabled_viz)}")
            enabled_lipid = [k for k, v in settings.get('viz_include_lipid_class', {}).items() if v]
            if enabled_lipid:
                saved_items.append(f"  ✓ Lipid class plots: {', '.join(enabled_lipid)}")
            non_none_gui_vars = {k: v for k, v in gui_vars.items() if v is not None}
            if non_none_gui_vars:
                saved_items.append(f"  ✓ GUI variables: {len(non_none_gui_vars)} settings from Visualization, Statistics & Pathway tabs")
            
            logger.info("=" * 60)
            logger.info("SETTINGS BEING SAVED:")
            for item in saved_items:
                logger.info(item)
            logger.info("=" * 60)
            
        except Exception as e:
            logger.error(f"Failed to gather settings: {e}")
            import traceback
            logger.error(traceback.format_exc())
        return settings

    def _apply_all_settings(self, data: dict):
        """Apply loaded settings with detailed logging."""
        loaded_items = []
        try:
            # Apply viz params
            viz_params = data.get('viz_params') or {}
            for k, v in (viz_params.items() if isinstance(viz_params, dict) else []):
                if k in self.viz_params and v is not None:
                    # Restore each attribute to the param object
                    if hasattr(self.viz_params[k], '__dict__'):
                        param_count = 0
                        for attr, val in v.items():
                            try:
                                setattr(self.viz_params[k], attr, val)
                                param_count += 1
                            except Exception as ex:
                                logger.warning(f"Could not set {k}.{attr}: {ex}")
                        if param_count > 0:
                            loaded_items.append(f"  ✓ {k} parameters: {param_count} settings")
            
            # Groups and mapping
            gd = data.get('viz_group_definitions') or {}
            gm = data.get('viz_group_mapping') or {}
            if gd:
                self.viz_group_definitions = gd
                loaded_items.append(f"  ✓ Group definitions: {len(gd)} groups ({', '.join(gd.values())})")
            if gm:
                self.viz_group_mapping = gm
                loaded_items.append(f"  ✓ Sample-to-group mapping: {len(gm)} samples assigned")
            
            # Color map
            cmap = data.get('viz_color_map') or {}
            if cmap:
                try:
                    # Normalize colors to hex strings
                    self.viz_color_map = {g: cmap[g] for g in cmap}
                    loaded_items.append(f"  ✓ Color map: {len(cmap)} colors ({', '.join([f'{g}:{c}' for g, c in cmap.items()])})")
                except Exception:
                    pass
            
            # Preferred order
            if 'viz_preferred_group_order' in data:
                order_str = data.get('viz_preferred_group_order', '')
                if order_str:
                    try:
                        self.viz_preferred_group_order.set(order_str)
                        loaded_items.append(f"  ✓ Preferred group order: {order_str}")
                    except Exception:
                        pass
            
            # Visualization selection checkboxes
            viz_sel = data.get('viz_selected') or {}
            if hasattr(self, 'viz_selected') and viz_sel:
                enabled = []
                for k, v in viz_sel.items():
                    if k in self.viz_selected:
                        try:
                            self.viz_selected[k].set(v)
                            if v:
                                enabled.append(k)
                        except Exception:
                            pass
                if enabled:
                    loaded_items.append(f"  ✓ Enabled visualizations: {', '.join(enabled)}")
            
            # Lipid class visualization options
            lipid_opts = data.get('viz_include_lipid_class') or {}
            if hasattr(self, 'viz_include_lipid_class') and lipid_opts:
                enabled_lipid = []
                for k, v in lipid_opts.items():
                    if k in self.viz_include_lipid_class:
                        try:
                            self.viz_include_lipid_class[k].set(v)
                            if v:
                                enabled_lipid.append(k)
                        except Exception:
                            pass
                if enabled_lipid:
                    loaded_items.append(f"  ✓ Lipid class plots enabled for: {', '.join(enabled_lipid)}")
            
            # Restore all GUI variables
            gui_vars = data.get('gui_vars') or {}
            gui_var_count = 0
            for var_name, value in gui_vars.items():
                if value is not None and hasattr(self, var_name):
                    try:
                        var = getattr(self, var_name)
                        if hasattr(var, 'set'):
                            var.set(value)
                            gui_var_count += 1
                    except Exception as ex:
                        logger.warning(f"Could not restore {var_name}: {ex}")
            if gui_var_count > 0:
                loaded_items.append(f"  ✓ GUI variables: {gui_var_count} settings restored")
            
            # Refresh UI widgets - COMPREHENSIVE UPDATE
            try:
                # Refresh group displays
                if hasattr(self, 'refresh_viz_group_display'):
                    self.refresh_viz_group_display()
                
                # Update color controls if groups were loaded
                if gd or gm:
                    groups = list(gd.values()) if gd else list(set(gm.values()))
                    if groups and hasattr(self, 'update_group_color_controls'):
                        self.update_group_color_controls(groups)
                
                # Update viz params display if needed
                if hasattr(self, '_update_viz_params'):
                    self._update_viz_params()
                
                loaded_items.append(f"  ✓ UI refreshed successfully")
            except Exception as e:
                logger.warning(f"Could not fully refresh UI: {e}")
            
            # Log detailed summary
            if loaded_items:
                logger.info("=" * 60)
                logger.info("SETTINGS LOADED SUCCESSFULLY:")
                for item in loaded_items:
                    logger.info(item)
                logger.info("=" * 60)
                self.log_viz_message(f"✓ Settings loaded successfully: {len(loaded_items)} categories")
            else:
                logger.warning("Settings file loaded but no items were applied")
                self.log_viz_message("⚠ Settings file loaded but no items were applied", tag='WARNING')
            
        except Exception as e:
            logger.error(f"Failed to apply settings: {e}")
            import traceback
            logger.error(traceback.format_exc())
    
    # ========== Tab Switching and Auto-Load Helpers ==========
    
    def switch_to_tab(self, tab_name: str):
        """
        Switch to a tab by name.
        
        Args:
            tab_name: The display name of the tab (e.g., "🔬 ID Annotation")
        
        Returns:
            bool: True if tab was switched, False otherwise
        """
        try:
            # Get GUI instance from notebook
            if hasattr(self.parent, 'gui_instance'):
                gui_instance = self.parent.gui_instance
                if hasattr(gui_instance, 'tab_instances') and tab_name in gui_instance.tab_instances:
                    # Find the tab index in the notebook
                    notebook = self.parent
                    for i in range(len(notebook.tabs())):
                        tab_frame = notebook.tabs()[i]
                        # Try to find by tab name comparison
                        for name in gui_instance.tab_instances:
                            if gui_instance.tab_instances[name].frame == self.parent.nametowidget(tab_frame):
                                if name == tab_name:
                                    notebook.select(i)
                                    logger.info(f"✅ Switched to tab: {tab_name}")
                                    return True
                    
                    # Alternative method: search by tab index in notebook
                    for i, tab_id in enumerate(notebook.tabs()):
                        try:
                            # The tab might have a name attribute or we can check through GUI
                            notebook.select(i)
                            if notebook.tab(i, "text") == tab_name:
                                logger.info(f"✅ Switched to tab: {tab_name}")
                                return True
                        except Exception:
                            pass
        except Exception as e:
            logger.warning(f"Could not switch to tab '{tab_name}': {e}")
        
        return False
    
    def get_tab_by_name(self, tab_name: str):
        """
        Get a tab instance by name.
        
        Args:
            tab_name: The display name of the tab
        
        Returns:
            BaseTab instance or None if not found
        """
        try:
            if hasattr(self.parent, 'gui_instance'):
                gui_instance = self.parent.gui_instance
                if hasattr(gui_instance, 'tab_instances'):
                    return gui_instance.tab_instances.get(tab_name)
        except Exception as e:
            logger.warning(f"Could not get tab '{tab_name}': {e}")
        
        return None
    
    def notify_data_ready(self, tab_name: str, data_key: str = None):
        """
        Notify another tab that new data is ready in memory_store.
        
        Args:
            tab_name: The name of the tab to notify
            data_key: Optional specific key in memory_store (e.g., 'cleaned_metabolites_df')
        """
        try:
            target_tab = self.get_tab_by_name(tab_name)
            if target_tab and hasattr(target_tab, 'on_data_ready'):
                target_tab.on_data_ready(data_key)
                logger.info(f"✅ Notified {tab_name} that data is ready" + (f" ({data_key})" if data_key else ""))
        except Exception as e:
            logger.warning(f"Could not notify {tab_name}: {e}")
    
    def on_data_ready(self, data_key: str = None):
        """
        Called when another tab has finished processing and data is ready.
        Subclasses can override to customize behavior.
        
        Args:
            data_key: Optional key in memory_store that was just loaded
        """
        # Override in subclass to handle data loading
        logger.debug(f"{self.__class__.__name__}.on_data_ready called with key: {data_key}")

