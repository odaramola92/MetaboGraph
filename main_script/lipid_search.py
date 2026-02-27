"""
python lipid_search.py --input cleaned_lipids --feather lipidmap.feather --out matches.xlsx --ppm 10
"""

import argparse
import math
import pandas as pd
import numpy as np
import re
from typing import Optional

DEFAULT_INPUT = "cleaned_lipids"
DEFAULT_FEATHER = "lipidmap.feather"

# Columns we expect (try to be flexible if some are missing)
# Priority order: most specific/reliable first
PREFERRED_COLS = [
    "ABBREVIATION",      # Highest priority: short lipid names like PC(16:0/18:0)
    "NAME",              # Standard lipid names
    "SYSTEMATIC_NAME",   # IUPAC/systematic names
    "CATEGORY",          # Broader categories
    "MAIN_CLASS",        # Class-level matching
    "SUB_CLASS",         # Sub-class matching
    "FORMULA",           # Chemical formula (requires confirmation)
    "EXACT_MASS",        # Mass matching (requires confirmation)
]

CONFIRM_PRIORITY = ["ABBREVIATION", "NAME", "SYSTEMATIC_NAME", "FORMULA", "MAIN_CLASS", "SUB_CLASS"]


def mass_within_ppm(m_input, m_candidate, ppm):
    """Return True if m_candidate is within +/- ppm of m_input."""
    try:
        m_input = float(m_input)
        m_candidate = float(m_candidate)
        ppm = float(ppm)
        if m_input == 0:
            return False
        return abs(m_candidate - m_input) <= (abs(m_input) * ppm / 1_000_000)
    except Exception:
        return False

def norm_str(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    return s if s != "" else None

def strip_suffix_brackets(s):
    if pd.isna(s):
        return ""
    s = str(s).strip()
    return re.sub(r"[\[\(\"'].*$", "", s).strip()

def clean_lipid_shorthand(s):
    """Normalize lipid abbreviation/name shorthand to a stable form.
    Examples:
      PC(O-32:1)+H -> PC(O-32:1)  [preserve parentheses for database matching]
      PE(16:1_18:1)+H -> PE(16:1/18:1)
      Cer(d16:0_16:0) -> Cer(d16:0/16:0)
      PC(14:0) -> PC(14:0)  [preserve simple forms]
    Returns stripped string, preserving case.
    """
    if pd.isna(s):
        return ""
    t = str(s).strip()
    # remove adducts like +H, +Na, etc
    t = re.sub(r"\+.*$", "", t)

    # If the shorthand is written like 'LPA 18:1' or 'FA 20:4', convert to 'LPA(18:1)'
    # only do this when there are no existing parentheses
    if '(' not in t and ')' not in t:
        m = re.match(r"^([A-Za-z0-9+-]+)\s+(.+)$", t)
        if m:
            prefix = m.group(1)
            rest = m.group(2).strip()
            # require a colon in the rest to indicate a fatty-acyl pattern like 18:1
            if ':' in rest:
                t = f"{prefix}({rest})"

    # Handle parenthetical content: replace underscores with slashes
    # but preserve the parentheses structure
    def replace_underscores_in_parens(match):
        prefix = match.group(1)
        content = match.group(2)
        # replace underscore with slash for fatty acyl pairs
        content = content.replace("_", "/")
        return f"{prefix}({content})"

    # Apply underscore->slash replacement inside parentheses
    t = re.sub(r"([A-Za-z0-9]+)\(([^)]+)\)", replace_underscores_in_parens, t)

    # collapse multiple spaces
    t = re.sub(r"\s+", " ", t).strip()
    return t

def display_lipid_shorthand(s):
    """Return a human-friendly display form (e.g. 'PA(20:0)+NH4' -> 'PA 20:0').
    This should not be used for matching; use `clean_lipid_shorthand` for normalization.
    """
    if pd.isna(s):
        return ''
    t = str(s).strip()
    # remove adducts like +H, +Na, etc
    t = re.sub(r"\+.*$", "", t)
    # replace parentheses and slashes with spaces for display
    t = t.replace('(', ' ').replace(')', ' ')
    t = t.replace('/', ' ')
    t = re.sub(r"\s+", " ", t).strip()
    return t

def find_matches_for_row(row, df_feather, ppm):
    """Try to find exact matches for a single input row.
    Returns a tuple (matches_df, match_type) where match_type is a short string like
    'NAME_exact' or 'EXACT_MASS_confirmed_by_NAME'. If no match, returns (empty_df, None).
    """
    feather = df_feather.copy()

    # Support common input aliases: allow users to supply alternate column names
    # Map alias -> canonical
    INPUT_ALIASES = {
        'LipidID': 'ABBREVIATION',
        'Class': 'MAIN_CLASS'
    }

    # Normalize the incoming row so canonical names exist when aliases are used.
    try:
        original_row = row
        row = original_row.copy()
        for alias, canon in INPUT_ALIASES.items():
            if alias in original_row.index and canon not in row.index:
                try:
                    row[canon] = original_row[alias]
                except Exception:
                    # If assignment fails for any reason, skip silently
                    pass
    except Exception:
        # If row isn't a Series-like object or copy fails, fall back to original
        pass

    # Special handling for LipidID: search ABBREVIATION first, then NAME if not found
    if 'LipidID' in original_row.index:
        lipidid_val = original_row['LipidID']
        if not pd.isna(lipidid_val) and (isinstance(lipidid_val, str) and lipidid_val.strip() != ""):
            # First try ABBREVIATION
            if 'ABBREVIATION' in feather.columns:
                abbrev_match = _try_string_match(row, feather, 'ABBREVIATION', lipidid_val, ppm, [])
                if abbrev_match is not None:
                    df_match, match_type, tried_list = abbrev_match
                    return df_match, f"LipidID_matched_ABBREVIATION", tried_list
            
            # If no match in ABBREVIATION, try NAME (Preferred_Name)
            if 'NAME' in feather.columns:
                name_match = _try_string_match(row, feather, 'NAME', lipidid_val, ppm, [])
                if name_match is not None:
                    df_match, match_type, tried_list = name_match
                    return df_match, f"LipidID_matched_NAME", tried_list

    # use top-level helpers (strip_suffix_brackets, clean_lipid_shorthand)

    tried = []
    # Track if ABBREVIATION was present and attempted but produced no match
    abbrev_tried_no_match = False
    # build normalized columns (case-sensitive)
    for c in ["FORMULA", "NAME", "ABBREVIATION", "SYSTEMATIC_NAME", "CATEGORY", "MAIN_CLASS", "SUB_CLASS"]:
        if c in feather.columns:
            if c in ("CATEGORY", "MAIN_CLASS", "SUB_CLASS"):
                feather[f"__{c}_NORM"] = feather[c].apply(lambda s: strip_suffix_brackets(s))
            elif c in ("NAME", "ABBREVIATION"):
                feather[f"__{c}_NORM"] = feather[c].fillna("").apply(lambda s: clean_lipid_shorthand(s))
            elif c in ("SYSTEMATIC_NAME",):
                feather[f"__{c}_NORM"] = feather[c].fillna("").astype(str).str.strip()
            else:
                feather[f"__{c}_NORM"] = feather[c].fillna("").astype(str).str.strip()

    # iterate search columns available in input BY PRIORITY ORDER
    # Try each column in priority order and return first match found
    for col in PREFERRED_COLS:
        if col not in row.index:
            continue

        val = row[col]
        # handle duplicate-column values (Series or ndarray)
        if isinstance(val, pd.Series):
            try:
                nonnull = val.dropna()
                val = nonnull.iloc[0] if not nonnull.empty else (val.iloc[0] if len(val) else None)
            except Exception:
                val = val.iloc[0] if len(val) else None
        elif isinstance(val, np.ndarray):
            try:
                arr = val[~pd.isna(val)]
                val = arr[0] if arr.size else (val[0] if val.size else None)
            except Exception:
                val = val[0] if getattr(val, 'size', 0) else None

        if pd.isna(val) or (isinstance(val, str) and val.strip() == ""):
            continue

        # Try to find match in this column
        match_result = None
        
        # EXACT_MASS handling
        if col == "EXACT_MASS":
            match_result = _try_mass_match(row, feather, val, ppm, tried)
        
        # FORMULA handling  
        elif col == "FORMULA":
            match_result = _try_formula_match(row, feather, val, ppm, tried)
        
        # Generic string columns: exact (case-sensitive) match
        else:
            match_result = _try_string_match(row, feather, col, val, ppm, tried)
        
        # If we found a match, return it immediately (priority order)
        if match_result is not None:
            # If we previously tried ABBREVIATION and it yielded no match,
            # but now NAME (or another lower-priority string) matched,
            # annotate the returned match_type to indicate the fallback.
            try:
                df_match, match_type, tried_list = match_result
            except Exception:
                # match_result may sometimes be a different shape; return as-is
                return match_result

            if abbrev_tried_no_match and col in ("NAME", "SYSTEMATIC_NAME"):
                match_type = f"{match_type}_fallback_from_ABBREVIATION"
            return df_match, match_type, tried_list
        
        # No match in this column, continue to next priority column
        # If we just tried ABBREVIATION and it didn't produce a match, remember that
        if col == "ABBREVIATION":
            # only mark tried-no-match if the value was non-empty (we actually attempted)
            if not pd.isna(val) and not (isinstance(val, str) and val.strip() == ""):
                # if match_result is None, that means ABBREVIATION produced no match
                abbrev_tried_no_match = True

    # No matches found in any column
    return pd.DataFrame(), None, tried

def _try_mass_match(row, feather, val, ppm, tried):
    """Try to match by EXACT_MASS with confirmation"""
    try:
        input_mass = float(val)
    except Exception:
        return None
    if "EXACT_MASS" not in feather.columns:
        return None
    mask = feather["EXACT_MASS"].notna() & feather["EXACT_MASS"].apply(lambda x: mass_within_ppm(input_mass, x, ppm))
    candidates = feather[mask]
    if candidates.empty:
        return None
    # require at least one confirming column
    confirmed = []
    for _, cand in candidates.iterrows():
        for conf in CONFIRM_PRIORITY:
            if conf in row.index and conf in feather.columns:
                v_in = norm_str(row[conf])
                v_cand = cand.get(f"__{conf}_NORM")
                if v_in and v_cand:
                    if conf in ("CATEGORY", "MAIN_CLASS", "SUB_CLASS"):
                        v_in_norm = strip_suffix_brackets(v_in)
                    elif conf in ("NAME", "ABBREVIATION"):
                        v_in_norm = clean_lipid_shorthand(v_in)
                    else:
                        v_in_norm = v_in.strip()
                    v_cand_norm = v_cand.strip() if isinstance(v_cand, str) else str(v_cand).strip()
                    if v_in_norm == v_cand_norm:
                        confirmed.append(cand)
                        break
    if confirmed:
        df = pd.DataFrame(confirmed)
        # find which conf matched for the first confirmed candidate
        first = confirmed[0]
        match_conf = None
        for conf in CONFIRM_PRIORITY:
            if conf in row.index and conf in feather.columns:
                v_in = norm_str(row[conf])
                v_cand = first.get(f"__{conf}_NORM")
                if v_in and v_cand:
                    if conf in ("CATEGORY", "MAIN_CLASS", "SUB_CLASS"):
                        v_in_norm = strip_suffix_brackets(v_in)
                    else:
                        v_in_norm = v_in.strip()
                    v_cand_norm = v_cand.strip() if isinstance(v_cand, str) else str(v_cand).strip()
                    if v_in_norm == v_cand_norm:
                        match_conf = conf
                        break
        match_type = f"EXACT_MASS_confirmed_by_{match_conf if match_conf else 'unknown'}"
        return df, match_type, tried
    return None

def _try_formula_match(row, feather, val, ppm, tried):
    """Try to match by FORMULA with confirmation"""
    if "FORMULA" not in feather.columns:
        return None
    v = norm_str(val)
    if not v:
        return None
    mask = feather["__FORMULA_NORM"] == v.strip()
    candidates = feather[mask]
    if candidates.empty:
        return None
    confirmed = []
    for _, cand in candidates.iterrows():
        name_ok = False
        mass_ok = False
        if "NAME" in row.index and "NAME" in feather.columns:
            v_in = norm_str(row["NAME"])
            v_cand = cand.get("__NAME_NORM")
            if v_in and v_cand and clean_lipid_shorthand(v_in) == (v_cand if isinstance(v_cand, str) else str(v_cand)):
                name_ok = True
        if "EXACT_MASS" in row.index and pd.notna(row["EXACT_MASS"]) and "EXACT_MASS" in feather.columns:
            try:
                mass_ok = mass_within_ppm(float(row["EXACT_MASS"]), cand["EXACT_MASS"], ppm)
            except Exception:
                mass_ok = False
        if name_ok or mass_ok:
            confirmed.append(cand)
    if confirmed:
        df = pd.DataFrame(confirmed)
        first = confirmed[0]
        match_conf = None
        if "NAME" in row.index and "NAME" in feather.columns:
            v_in = norm_str(row["NAME"])
            v_cand = first.get("__NAME_NORM")
            if v_in and v_cand and clean_lipid_shorthand(v_in) == (v_cand if isinstance(v_cand, str) else str(v_cand)):
                match_conf = "NAME"
        if match_conf is None and "EXACT_MASS" in row.index and pd.notna(row["EXACT_MASS"]):
            try:
                if mass_within_ppm(float(row["EXACT_MASS"]), first["EXACT_MASS"], ppm):
                    match_conf = "EXACT_MASS"
            except Exception:
                pass
        match_type = f"FORMULA_confirmed_by_{match_conf if match_conf else 'unknown'}"
        return df, match_type, tried
    return None

def _try_string_match(row, feather, col, val, ppm, tried):
    """Try to match by string column (NAME, ABBREVIATION, etc.)"""
    if col not in feather.columns:
        return None
    v = norm_str(val)
    if not v:
        return None
    key = f"__{col}_NORM" if f"__{col}_NORM" in feather.columns else None
    if key:
        if col in ("CATEGORY", "MAIN_CLASS", "SUB_CLASS"):
            v_norm = strip_suffix_brackets(v)
        elif col in ("NAME", "ABBREVIATION"):
            v_norm = clean_lipid_shorthand(v)
        elif col in ("SYSTEMATIC_NAME",):
            v_norm = v.strip()
        else:
            v_norm = v.strip()
        tried.append((col, v, v_norm))
        mask = feather[key] == v_norm
        candidates = feather[mask]
    else:
        if col in ("CATEGORY", "MAIN_CLASS", "SUB_CLASS"):
            v_norm = strip_suffix_brackets(v)
            tried.append((col, v, v_norm))
            mask = feather[col].astype(str).apply(lambda s: strip_suffix_brackets(s)) == v_norm
        elif col in ("NAME", "ABBREVIATION"):
            v_norm = clean_lipid_shorthand(v)
            tried.append((col, v, v_norm))
            mask = feather[col].astype(str).apply(lambda s: clean_lipid_shorthand(s)) == v_norm
        elif col in ("SYSTEMATIC_NAME",):
            v_norm = v.strip()
            tried.append((col, v, v_norm))
            mask = feather[col].astype(str).str.strip() == v_norm
        else:
            v_norm = v.strip()
            tried.append((col, v, v_norm))
            mask = feather[col].astype(str).str.strip() == v_norm
        candidates = feather[mask]

    if not candidates.empty:
        PRIMARY_KEYS = {"ABBREVIATION", "NAME", "SYSTEMATIC_NAME", "MAIN_CLASS", "SUB_CLASS"}
        if col in PRIMARY_KEYS:
            # require secondary confirmation by EXACT_MASS (±ppm) or normalized FORMULA
            # unless neither is provided in the input row
            def f_norm(x):
                if pd.isna(x):
                    return None
                return str(x).replace(" ", "").upper()
            have_mass = ("EXACT_MASS" in row.index and pd.notna(row["EXACT_MASS"]) and "EXACT_MASS" in feather.columns)
            have_formula = ("FORMULA" in row.index and norm_str(row["FORMULA"]) is not None and "FORMULA" in feather.columns)

            if not have_mass and not have_formula:
                match_type = f"{col}_exact"
                return candidates, match_type, tried

            confirmed = []
            confirmed_by = None
            for _, cand in candidates.iterrows():
                ok = False
                if have_mass:
                    try:
                        if mass_within_ppm(float(row["EXACT_MASS"]), cand["EXACT_MASS"], ppm):
                            ok = True
                            if confirmed_by is None:
                                confirmed_by = "EXACT_MASS"
                    except Exception:
                        pass
                if (not ok) and have_formula:
                    try:
                        if f_norm(row["FORMULA"]) == f_norm(cand.get("FORMULA")):
                            ok = True
                            if confirmed_by is None:
                                confirmed_by = "FORMULA"
                    except Exception:
                        pass
                if ok:
                    confirmed.append(cand)
            if confirmed:
                df = pd.DataFrame(confirmed)
                match_type = f"{col}_exact_confirmed_by_{confirmed_by if confirmed_by else 'unknown'}"
                return df, match_type, tried
            # no secondary confirmation -> no match, continue to next priority column
            return None
        else:
            match_type = f"{col}_exact"
            return candidates, match_type, tried
    
    return None

def run_search(input_df: pd.DataFrame, df_feather: pd.DataFrame, ppm: float = 10.0, progress_callback=None, progress_interval: Optional[int]=None, max_workers: int = 1) -> pd.DataFrame:
    """Run lipid search matching and return the output DataFrame.
    This is reusable from other scripts.
    """
    # Small startup log so callers can see whether matching is parallelized
    try:
        mode = 'parallel' if max_workers and max_workers > 1 else 'sequential'
        print(f"[lipid_search] Starting run_search in {mode} mode (max_workers={max_workers})")
    except Exception:
        pass
    # Apply column mapping to ensure consistent column names
    def map_input_columns(df):
        pref_lower = {p.lower(): p for p in PREFERRED_COLS}
        col_map = {}
        for c in df.columns:
            key = c.strip().lower()
            # Do not map class-related columns to NAME/ABBREVIATION for LipidMaps matching
            if 'class' in key:
                continue
            if key in pref_lower:
                col_map[c] = pref_lower[key]
                continue
            # common variants
            if key == 'name':
                col_map[c] = 'NAME'
                continue
            if key == 'abbreviation':
                col_map[c] = 'ABBREVIATION'
                continue
            compact = key.replace(' ', '')
            if compact in pref_lower:
                col_map[c] = pref_lower[compact]
                continue
            # try fuzzy partial match: prefer if preferred token is subset of key
            for plower, canonical in pref_lower.items():
                if plower in key or key in plower:
                    col_map[c] = canonical
                    break
        if col_map:
            df = df.rename(columns=col_map)
        return df, col_map

    # Apply column mapping
    input_df_mapped, col_map = map_input_columns(input_df.copy())
    
    rows_out = []
    input_cols = list(input_df.columns)  # Keep original column names for output
    mapped_cols = list(input_df_mapped.columns)  # Use mapped names for searching
    feather_cols = list(df_feather.columns)

    total = len(input_df_mapped)
    if progress_interval is None:
        try:
            # Report more frequently for smoother progress: every 1% or at least every 5 items
            progress_interval = max(5, total // 100)
        except Exception:
            progress_interval = 5

    # Prepare cleaned display names: prefer ABBREVIATION then NAME then raw first non-empty column
    display_names = []
    for i, r in input_df_mapped.iterrows():
        dn = None
        for c in ("ABBREVIATION", "NAME"):
            if c in r.index:
                val = r[c]
                # Handle scalar check properly for pandas Series
                if pd.notna(val) if not isinstance(val, pd.Series) else val.notna().any():
                    v = str(val).strip()
                    if v:
                        dn = display_lipid_shorthand(v)
                        break
        if dn is None:
            # fallback to first non-empty original column
            for c in input_df.columns:
                try:
                    v = input_df.iloc[i].get(c)
                    if isinstance(v, str) and v.strip():
                        dn = v.strip()
                        break
                except Exception:
                    continue
        display_names.append(dn or '')

    # If max_workers <= 1, run sequentially (preserves original ordering)
    if max_workers <= 1:
        for i, r in input_df_mapped.iterrows():
            matches, match_type, _ = find_matches_for_row(r, df_feather, ppm)

            if matches is None or (hasattr(matches, 'empty') and matches.empty):
                    out = {"matched": False, "match_type": None}
                    orig_row = input_df.iloc[i]
                    # Attach original input columns
                    for c in input_cols:
                        try:
                            out[f"input_{c}"] = orig_row.get(c) if c in orig_row.index else None
                        except Exception:
                            out[f"input_{c}"] = None
                    for c in feather_cols:
                        out[c] = None
                    rows_out.append(out)

            else:
                try:
                    orig_row = input_df.iloc[i]
                    print(f"[lipid_search DEBUG] Matched input row {i}: display='{display_lipid_shorthand(orig_row.get('ABBREVIATION') or orig_row.get('NAME') or '')}', match_type='{match_type}', hits={len(matches)}")
                except Exception:
                    pass
                for _, m in matches.iterrows():
                    out = {"matched": True, "match_type": match_type}
                    orig_row = input_df.iloc[i]
                    for c in input_cols:
                        try:
                            out[f"input_{c}"] = orig_row.get(c) if c in orig_row.index else None
                        except Exception:
                            out[f"input_{c}"] = None
                    for c in feather_cols:
                        out[c] = m.get(c) if c in m.index else None
                    rows_out.append(out)

            # progress callback per input row (report 1-based index)
            try:
                if progress_callback and ((i + 1) % progress_interval == 0 or (i + 1) == total):
                    name_for_msg = display_names[i] if i < len(display_names) else None
                    msg = f"Processed {i+1}/{total}"
                    if name_for_msg:
                        msg += f": {name_for_msg}"
                    progress_callback(i + 1, total, msg)
            except Exception:
                pass
    else:
        # Parallel execution using ThreadPoolExecutor. Each task is independent and returns list of output rows.
        from concurrent.futures import ThreadPoolExecutor, as_completed

        def _process_index(idx_row):
            i, r = idx_row
            try:
                matches, match_type, _ = find_matches_for_row(r, df_feather, ppm)
            except Exception:
                matches = None
                match_type = None

            local_rows = []
            if matches is None or (hasattr(matches, 'empty') and matches.empty):
                out = {"matched": False, "match_type": None}
                orig_row = input_df.iloc[i]
                for c in input_cols:
                    try:
                        out[f"input_{c}"] = orig_row.get(c) if c in orig_row.index else None
                    except Exception:
                        out[f"input_{c}"] = None
                for c in feather_cols:
                    out[c] = None
                local_rows.append((i, out))
            else:
                try:
                    orig_row = input_df.iloc[i]
                    print(f"[lipid_search DEBUG] Matched (parallel) input row {i}: display='{display_lipid_shorthand(orig_row.get('ABBREVIATION') or orig_row.get('NAME') or '')}', match_type='{match_type}', hits={len(matches)}")
                except Exception:
                    pass
                for _, m in matches.iterrows():
                    out = {"matched": True, "match_type": match_type}
                    orig_row = input_df.iloc[i]
                    for c in input_cols:
                        try:
                            out[f"input_{c}"] = orig_row.get(c) if c in orig_row.index else None
                        except Exception:
                            out[f"input_{c}"] = None
                    for c in feather_cols:
                        out[c] = m.get(c) if c in m.index else None
                    local_rows.append((i, out))
            return local_rows

        indices = list(input_df_mapped.iterrows())
        completed = 0
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_idx = {executor.submit(_process_index, ir): ir[0] for ir in indices}
            for fut in as_completed(future_to_idx):
                idx = future_to_idx[fut]
                try:
                    local_rows = fut.result()
                    # append local rows (keep original input order by storing index)
                    for t in local_rows:
                        rows_out.append(t[1])
                except Exception:
                    pass

                completed += 1
                try:
                    # report progress per completed task (not per output row)
                    if progress_callback and ((completed) % progress_interval == 0 or completed == total):
                        # attempt to get a cleaned display name for message from the input row
                        name_for_msg = display_names[idx] if idx < len(display_names) else None
                        msg = f"Processed {completed}/{total}"
                        if name_for_msg:
                            msg += f": {name_for_msg}"
                        progress_callback(completed, total, msg)
                        
                        # Small sleep to allow GUI event loop to process (prevent freezing)
                        import time
                        time.sleep(0.001)  # 1ms pause every progress_interval items
                except Exception:
                    pass

    df_out = pd.DataFrame(rows_out)
    # Synthesize a stable LipidID column for output:
    # - Prefer the original input column named 'LipidID' if present (as 'input_LipidID')
    # - Otherwise, fallback to ABBREVIATION from the matched feather row
    # - As a last resort, use a cleaned display from input NAME/ABBREVIATION
    try:
        lipidid_series = None
        if 'input_LipidID' in df_out.columns:
            lipidid_series = df_out['input_LipidID']
        elif 'ABBREVIATION' in df_out.columns:
            lipidid_series = df_out['ABBREVIATION']
        else:
            # Build from any available input column that looks like a name
            def _make_display(row):
                for c in ['input_ABBREVIATION', 'input_NAME', 'input_Name']:
                    if c in row and pd.notna(row[c]) and str(row[c]).strip():
                        try:
                            return display_lipid_shorthand(row[c])
                        except Exception:
                            return str(row[c]).strip()
                return ''
            lipidid_series = df_out.apply(_make_display, axis=1)
        df_out.insert(0, 'LipidID', lipidid_series)
    except Exception:
        # If anything goes wrong, continue without failing
        pass

    # Drop MoleculeName from final output if present
    if 'MoleculeName' in df_out.columns:
        df_out = df_out.drop(columns=['MoleculeName'])
    return df_out

def main():
    parser = argparse.ArgumentParser(description="Exact-match lipid search against lipidmap.feather")
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Input Excel file (default: cleaned_lipids.xlsx)")
    parser.add_argument("--feather", default=DEFAULT_FEATHER, help="Feather db to search (default: lipidmap.feather)")
    parser.add_argument("--out", default=None, help="Optional output Excel/CSV file to write matches")
    parser.add_argument("--ppm", type=float, default=10.0, help="PPM tolerance for EXACT_MASS (default 10)")
    parser.add_argument("--workers", "-w", type=int, default=1, help="Number of parallel workers (default 1 = sequential)")
    args = parser.parse_args()

    # load files
    print(f"Loading feather: {args.feather}")
    df_feather = pd.read_feather(args.feather)
    print(f"Feather columns: {list(df_feather.columns)}")
    print(f"Unique categories: {df_feather['CATEGORY'].nunique() if 'CATEGORY' in df_feather.columns else 'N/A'}")

    # Resolve input path: allow users to pass base name without extension (e.g. 'cleaned_lipids')
    from pathlib import Path

    def _find_input_file(inp: str) -> str:
        p = Path(str(inp))
        # If exact path exists, return it
        if p.exists():
            return str(p)
        # Try common extensions
        for ext in ('.xlsx', '.xls', '.csv'):
            cand = p.with_suffix(ext)
            if cand.exists():
                return str(cand)
        # Try adding .xlsx if user passed a basename without suffix
        for name in (f"{inp}.xlsx", f"{inp}.xls", f"{inp}.csv"):
            if Path(name).exists():
                return name
        # Not found
        raise FileNotFoundError(f"Input file not found: '{inp}'. Tried common extensions: .xlsx, .xls, .csv")

    print(f"Loading input: {args.input}")
    try:
        input_path = _find_input_file(args.input)
    except FileNotFoundError as e:
        # Give a helpful error and exit
        print(f"Error: {e}")
        print("Please provide the input file path (e.g. --input cleaned_lipids.xlsx) or place a file named cleaned_lipids.xlsx/csv in the working directory.")
        raise

    # Prefer Excel reading for .xls/.xlsx, otherwise fallback to CSV
    try:
        if str(input_path).lower().endswith(('.xls', '.xlsx')):
            df_input = pd.read_excel(input_path)
        else:
            # attempt CSV first for .csv or other
            df_input = pd.read_csv(input_path)
    except Exception:
        # Last-resort: try read_excel then read_csv to be forgiving
        try:
            df_input = pd.read_excel(input_path)
        except Exception:
            df_input = pd.read_csv(input_path)
    # normalize input column names to canonical expected names (case-insensitive)
    def map_input_columns(df):
        pref_lower = {p.lower(): p for p in PREFERRED_COLS}
        col_map = {}
        for c in df.columns:
            key = c.strip().lower()
            if key in pref_lower:
                col_map[c] = pref_lower[key]
                continue
            # common variants
            if key == 'name':
                col_map[c] = 'NAME'
                continue
            compact = key.replace(' ', '')
            if compact in pref_lower:
                col_map[c] = pref_lower[compact]
                continue
            # try fuzzy partial match: prefer if preferred token is subset of key
            for plower, canonical in pref_lower.items():
                if plower in key or key in plower:
                    col_map[c] = canonical
                    break
        if col_map:
            df = df.rename(columns=col_map)
        return df, col_map

    df_input, col_map = map_input_columns(df_input)
    print(f"Input columns found: {list(df_input.columns)}")
    if col_map:
        print(f"Mapped input columns: {col_map}")
    print(f"Input rows: {len(df_input)}")

    # Run reusable search (pass workers to enable parallel mode)
    df_out = run_search(df_input, df_feather, args.ppm, max_workers=args.workers)

    # Decide output path: use --out if provided, otherwise save next to the input with suffix _matches.xlsx
    from pathlib import Path
    if args.out:
        out_path = args.out
    else:
        inp = Path(args.input)
        out_path = str(inp.with_name(inp.stem + "_matches.xlsx"))

    # Write the file (Excel when extension looks like xls/xlsx, otherwise CSV)
    if out_path.lower().endswith(('.xls', '.xlsx')):
        df_out.to_excel(out_path, index=False)
    else:
        df_out.to_csv(out_path, index=False)
    print(f"Wrote match summary to {out_path}")


if __name__ == '__main__':
    main()
