import numpy as np
import pandas as pd
from scipy import stats
from collections import defaultdict
from pathlib import Path
 
ROWS      = list("ABCDEFGH")
COLS      = list(range(1, 13))
STANDARDS = {r: 10 ** (8 - i) for i, r in enumerate(ROWS)}  # A-1e8 … H-1e1
 
 
# Excel loading
 
def _normalise_col(s):
    return str(s).strip().replace("\u0442", "T").replace("\u0422", "T").upper()
 
 
def load_results(filepath):
    """Read the Results sheet; return a tidy DataFrame with Well, Row, Col, Ct."""
    raw = pd.read_excel(filepath, sheet_name="Results", header=None)
    header_row = next(
        i for i, row in raw.iterrows()
        if any(str(v).strip() == "Well" for v in row)
    )
    df = pd.read_excel(filepath, sheet_name="Results", header=header_row)
    df.columns = df.columns.str.strip()
 
    ct_col = next((c for c in df.columns if _normalise_col(c) == "CT"), None)
    if ct_col is None:
        raise ValueError(f"Ct column not found. Available: {list(df.columns)}")
 
    df = df[df["Well"].notna()].copy()
    df["Row"] = df["Well"].str[0].str.upper()
    df["Col"] = pd.to_numeric(df["Well"].str[1:], errors="coerce")
    df["Ct"]  = pd.to_numeric(df[ct_col], errors="coerce")
    df = df[df["Row"].isin(ROWS) & df["Col"].isin(COLS)]
 
    return df[["Well", "Row", "Col", "Ct"]].reset_index(drop=True)
 
 
def build_plate_grid(df):
    """Pivot tidy DataFrame into an 8×12 plate grid of formatted strings."""
    plate = pd.DataFrame("", index=ROWS, columns=COLS, dtype=object)
    for _, r in df.iterrows():
        plate.at[r["Row"], r["Col"]] = "" if pd.isna(r["Ct"]) else f"{r['Ct']:.2f}"
    return plate
 
 
# Plate display
 
def print_plate(plate):
    col_w = 13
    print("\n96-WELL PLATE - Ct VALUES")
    print("-" * (4 + col_w * 12))
    print("    " + "".join(str(c).center(col_w) for c in COLS))
    print("    " + "-" * (col_w * 12))
    for row_label in ROWS:
        vals = [plate.at[row_label, c] for c in COLS]
        print(f"  {row_label} |" + "".join(
            ("Und" if v == "" else v).center(col_w) for v in vals
        ))
    print()
 
 
# Well selection helpers
 
def parse_cols(raw):
    """'1,2' or '1-4' or '3'  -  sorted list of ints."""
    cols = []
    for part in str(raw).replace(" ", "").split(","):
        if "-" in part:
            a, b = part.split("-", 1)
            cols.extend(range(int(a), int(b) + 1))
        else:
            cols.append(int(part))
    return sorted(set(cols))
 
 
def parse_rows(raw):
    """'A-D' or 'A,B,C' or 'A'  -  sorted list of row letters."""
    rows = []
    for part in str(raw).upper().replace(" ", "").split(","):
        if "-" in part:
            a, b = part.split("-", 1)
            rows.extend(ROWS[ROWS.index(a): ROWS.index(b) + 1])
        else:
            rows.append(part)
    return sorted(set(rows), key=lambda r: ROWS.index(r))
 
 
# Standard curve fitting
 
def fit_curve(df, std_cols, std_rows):
    """
    Fit a standard curve for the given well selection.
 
    Parameters
    ----------
    df        : DataFrame from load_results()
    std_cols  : list of ints  (e.g. [1, 2])
    std_rows  : list of strs  (e.g. ['A','B','C','D','E','F','G','H'])
 
    Returns
    -------
    slope, intercept, r2, efficiency, ct_pairs
    """
    mask   = df["Col"].isin(std_cols) & df["Row"].isin(std_rows) & df["Ct"].notna()
    std_df = df[mask].copy()
    std_df["CopyNumber"] = std_df["Row"].map(STANDARDS)
    ct_pairs = [(r["CopyNumber"], r["Ct"]) for _, r in std_df.iterrows()]
 
    if len(ct_pairs) < 3:
        raise ValueError("Fewer than 3 valid standard points - cannot fit curve.")
 
    copies = np.array([p[0] for p in ct_pairs])
    cts    = np.array([p[1] for p in ct_pairs])
    slope, intercept, r, *_ = stats.linregress(np.log10(copies), cts)
    r2         = r ** 2
    efficiency = (10 ** (-1 / slope) - 1) * 100
 
    return slope, intercept, r2, efficiency, ct_pairs
 
 
def ct_to_copies(ct, slope, intercept):
    return 10 ** ((ct - intercept) / slope)
 
 
# Dilution series back-calculation 
 
def calc_dilution_series(df, unk_cols, unk_rows, slope, intercept,
                         ct_min, ct_max, sample_volume_ml, dilution_factor=10):
    """
    For duplicate columns representing a serial dilution series, average the
    duplicates per row, calculate copies in each reaction, then back-calculate
    the copies/mL in the original undiluted sample.
 
    Parameters
    ----------
    df                : DataFrame from load_results()
    unk_cols          : list of 2 column ints  (the duplicate columns)
    unk_rows          : list of row letters in dilution order (first = most concentrated)
    slope, intercept  : from fit_curve()
    ct_min, ct_max    : Ct range of the standard curve
    sample_volume_ml  : volume of original sample added to the first (most concentrated) row
    dilution_factor   : fold dilution between each row (default 10)
 
    Returns
    -------
    list of dicts, one per row, plus a summary dict
    """
    mask = df["Col"].isin(unk_cols) & df["Row"].isin(unk_rows) & df["Ct"].notna()
    sub  = df[mask].copy()
 
    rows_data = []
    for i, row_label in enumerate(unk_rows):
        row_wells = sub[sub["Row"] == row_label]
        if row_wells.empty:
            continue
 
        cts      = row_wells["Ct"].values
        mean_ct  = float(np.mean(cts))
        sd_ct    = float(np.std(cts, ddof=1)) if len(cts) > 1 else 0.0
        wells    = ", ".join(sorted(row_wells["Well"].values))
 
        # Volume of original sample effectively in this reaction:
        # Row 0: sample_volume_ml
        # Row 1: sample_volume_ml / dilution_factor
        # Row i: sample_volume_ml / (dilution_factor ** i)
        effective_vol_ml = sample_volume_ml / (dilution_factor ** i)
 
        in_range = ct_min <= mean_ct <= ct_max
 
        if in_range:
            copies_in_rxn    = ct_to_copies(mean_ct, slope, intercept)
            copies_per_ml    = copies_in_rxn / effective_vol_ml
        else:
            copies_in_rxn = None
            copies_per_ml = None
 
        rows_data.append({
            "row":             row_label,
            "wells":           wells,
            "n":               len(cts),
            "mean_ct":         mean_ct,
            "sd_ct":           sd_ct,
            "dilution_step":   i,
            "effective_vol_ml": effective_vol_ml,
            "copies_in_rxn":   copies_in_rxn,
            "copies_per_ml":   copies_per_ml,
            "in_range":        in_range,
        })
 
    # Average copies/mL across all in-range rows
    valid_concs = [r["copies_per_ml"] for r in rows_data if r["copies_per_ml"] is not None]
    mean_copies_per_ml = float(np.mean(valid_concs))   if valid_concs else None
    sd_copies_per_ml   = float(np.std(valid_concs, ddof=1)) if len(valid_concs) > 1 else 0.0
 
    summary = {
        "mean_copies_per_ml": mean_copies_per_ml,
        "sd_copies_per_ml":   sd_copies_per_ml,
        "n_valid_rows":        len(valid_concs),
    }
 
    return rows_data, summary
 
 
def print_dilution_series(name, rows_data, summary,
                          sample_volume_ml, dilution_factor):
    print(f"\n  DILUTION SERIES - {name}")
    print(f"  Sample volume in row 1: {sample_volume_ml} mL  |  "
          f"Dilution factor: 1:{dilution_factor:.0f}")
    print("  " + "-" * 72)
    print(f"  {'Row':<5} {'Wells':<12} {'Mean Ct':<10} {'SD':<8} "
          f"{'Copies/rxn':<16} {'Copies/mL orig':<18} {'Note'}")
    print("  " + "-" * 72)
 
    for r in rows_data:
        if r["copies_per_ml"] is not None:
            rxn_str = f"{r['copies_in_rxn']:.3e}"
            cpm_str = f"{r['copies_per_ml']:.3e}"
            note    = ""
        else:
            rxn_str = "-"
            cpm_str = "-"
            note    = "outside std range"
        print(f"  {r['row']:<5} {r['wells']:<12} {r['mean_ct']:<10.2f} "
              f"{r['sd_ct']:<8.3f} {rxn_str:<16} {cpm_str:<18} {note}")
 
    print("  " + "-" * 72)
    if summary["mean_copies_per_ml"] is not None:
        print(f"  Average copies/mL (original sample): "
              f"{summary['mean_copies_per_ml']:.3e}  "
              f"± {summary['sd_copies_per_ml']:.3e}  "
              f"(n={summary['n_valid_rows']} rows)")
    else:
        print("  No in-range rows - cannot calculate original concentration.")
 
 
# Result printing
 
def print_curve_summary(name, slope, intercept, r2, efficiency, ct_pairs):
    print("\n" + "═" * 52)
    print(f"  STANDARD CURVE - {name}")
    print("═" * 52)
    print(f"  Equation  : Ct = {slope:.4f} × log₁₀(copies) + {intercept:.4f}")
    print(f"  R²        : {r2:.5f}  {'✓' if r2 >= 0.98 else '✗  (should be ≥ 0.98)'}")
    print(f"  Efficiency: {efficiency:.1f}%  {'✓' if 90 <= efficiency <= 110 else '✗  (should be 90–110%)'}")
    print(f"  Slope     : {slope:.4f}  {'✓' if -3.6 <= slope <= -3.1 else '✗  (should be -3.6 to -3.1)'}")
 
    grouped = defaultdict(list)
    for copies, ct in ct_pairs:
        grouped[copies].append(ct)
 
    print(f"\n  {'Copy #':<12} {'n':<4} {'Mean Ct':<10} {'SD':<8} {'CV%'}")
    print("  " + "-" * 46)
    for copies in sorted(grouped.keys(), reverse=True):
        cts_list = grouped[copies]
        mean = np.mean(cts_list)
        sd   = np.std(cts_list, ddof=1) if len(cts_list) > 1 else 0.0
        cv   = (sd / mean * 100) if mean else 0.0
        print(f"  {copies:<12.0e} {len(cts_list):<4} {mean:<10.2f} {sd:<8.3f} {cv:.1f}%")
 
 
def print_unknowns(name, unknowns, slope, intercept, ct_min, ct_max):
    print(f"\n  UNKNOWNS - {name}")
    print("  " + "-" * 50)
    print(f"  {'Well':<8} {'Ct':<10} {'Gene Copies':<18} {'Note'}")
    print("  " + "-" * 50)
    for well, ct in unknowns:
        copies = ct_to_copies(ct, slope, intercept)
        note   = "outside std range" if not (ct_min <= ct <= ct_max) else ""
        print(f"  {well:<8} {ct:<10.2f} {copies:<18.3e} {note}")
 
 
# TSV export
 
def export_tsv(filepath, plate, assay_results):
    """Save plate layout + all assay results as a TSV next to the input file."""
    lines = []
 
    # Plate layout
    lines.append("PLATE LAYOUT (Ct values)")
    lines.append("\t" + "\t".join(str(c) for c in COLS))
    for row_label in ROWS:
        vals    = [plate.at[row_label, c] for c in COLS]
        display = ["Und" if v == "" else v for v in vals]
        lines.append(row_label + "\t" + "\t".join(display))
    lines.append("")
 
    for res in assay_results:
        name = res["name"]
 
        # Standard curve
        grouped = defaultdict(list)
        for copies, ct in res["ct_pairs"]:
            grouped[copies].append(ct)
 
        lines.append(f"STANDARD CURVE - {name}")
        lines.append(f"Equation\tCt = {res['slope']:.4f} x log10(copies) + {res['intercept']:.4f}")
        lines.append(f"R2\t{res['r2']:.5f}\t{'PASS' if res['r2'] >= 0.98 else 'FAIL'}")
        lines.append(f"Efficiency (%)\t{res['efficiency']:.2f}\t{'PASS' if 90 <= res['efficiency'] <= 110 else 'FAIL'}")
        lines.append(f"Slope\t{res['slope']:.4f}\t{'PASS' if -3.6 <= res['slope'] <= -3.1 else 'FAIL'}")
        lines.append("")
        lines.append("Copy #\tn\tMean Ct\tSD\tCV%")
        for copies in sorted(grouped.keys(), reverse=True):
            cts_list = grouped[copies]
            mean = np.mean(cts_list)
            sd   = np.std(cts_list, ddof=1) if len(cts_list) > 1 else 0.0
            cv   = (sd / mean * 100) if mean else 0.0
            lines.append(f"{copies:.0e}\t{len(cts_list)}\t{mean:.2f}\t{sd:.3f}\t{cv:.1f}")
        lines.append("")
 
        # Regular unknowns
        if res.get("unknowns"):
            lines.append(f"UNKNOWNS - {name}")
            lines.append("Well\tCt\tGene Copies\tNote")
            for well, ct in res["unknowns"]:
                copies = ct_to_copies(ct, res["slope"], res["intercept"])
                note   = "outside std range" if not (res["ct_min"] <= ct <= res["ct_max"]) else ""
                lines.append(f"{well}\t{ct:.2f}\t{copies:.3e}\t{note}")
            lines.append("")
 
        # Dilution series unknowns
        for ds in res.get("dilution_series", []):
            lines.append(f"DILUTION SERIES - {name}  |  sample: {ds['sample_name']}")
            lines.append(f"Sample volume in row 1 (mL)\t{ds['sample_volume_ml']}")
            lines.append(f"Dilution factor\t1:{ds['dilution_factor']:.0f}")
            lines.append("")
            lines.append("Row\tWells\tMean Ct\tSD\tCopies/rxn\tCopies/mL original\tNote")
            for r in ds["rows_data"]:
                if r["copies_per_ml"] is not None:
                    rxn_str = f"{r['copies_in_rxn']:.3e}"
                    cpm_str = f"{r['copies_per_ml']:.3e}"
                    note    = ""
                else:
                    rxn_str = ""
                    cpm_str = ""
                    note    = "outside std range"
                lines.append(
                    f"{r['row']}\t{r['wells']}\t{r['mean_ct']:.2f}\t{r['sd_ct']:.3f}"
                    f"\t{rxn_str}\t{cpm_str}\t{note}"
                )
            s = ds["summary"]
            if s["mean_copies_per_ml"] is not None:
                lines.append(
                    f"AVERAGE\t\t\t\t\t{s['mean_copies_per_ml']:.3e}"
                    f"\t(SD {s['sd_copies_per_ml']:.3e}, n={s['n_valid_rows']})"
                )
            lines.append("")
 
        lines.append("")
 
    out = Path(filepath).parent / (Path(filepath).stem + ".results.tsv")
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"\n  Results saved {out}\n")
 
 
# Full analysis pipeline
 
def run_analysis(filepath, assays):
    """
    Run a complete qPCR analysis.
 
    Parameters
    ----------
    filepath : str
        Path to the ABI Excel export (.xlsx).
 
    assays : list of dicts, each with keys:
        name       – label for this assay (str)
        std_cols   – standard columns as string or list, e.g. '1,2' or [1,2]
        std_rows   – standard rows   as string or list, e.g. 'A-H'
        unk_cols   – (optional) plain unknown columns
        unk_rows   – (optional) plain unknown rows (default 'A-H')
        dilutions  – (optional) list of dilution series dicts, each with:
            sample_name      – label for this sample
            cols             – the two duplicate columns, e.g. '7,8'
            rows             – rows in dilution order, e.g. 'A-H'
            sample_volume_ml – mL of original sample added to first row
            dilution_factor  – fold dilution per row (default 10)
 
    Example
    -------
    run_analysis("data.xlsx", [
        {
            "name":     "ITR",
            "std_cols": "1,2",
            "std_rows": "A-H",
            "unk_cols": "5,6",      # plain unknowns (no back-calculation)
            "unk_rows": "A-H",
            "dilutions": [
                {
                    "sample_name":      "SampleA",
                    "cols":             "7,8",
                    "rows":             "A-H",
                    "sample_volume_ml": 0.01,
                    "dilution_factor":  10,
                },
            ],
        },
    ])
    """
    print(f"\nLoading {filepath} ...")
    df    = load_results(filepath)
    plate = build_plate_grid(df)
    print_plate(plate)
 
    assay_results = []
 
    for assay in assays:
        name     = assay["name"]
        std_cols = parse_cols(assay["std_cols"]) if isinstance(assay["std_cols"], str) else assay["std_cols"]
        std_rows = parse_rows(assay["std_rows"]) if isinstance(assay["std_rows"], str) else assay["std_rows"]
 
        try:
            slope, intercept, r2, efficiency, ct_pairs = fit_curve(df, std_cols, std_rows)
        except ValueError as e:
            print(f"\n  ✗ {name}: {e}")
            continue
 
        ct_min = min(ct for _, ct in ct_pairs)
        ct_max = max(ct for _, ct in ct_pairs)
 
        print_curve_summary(name, slope, intercept, r2, efficiency, ct_pairs)
 
        # Plain unknowns 
        unknowns = []
        if assay.get("unk_cols"):
            unk_cols = parse_cols(assay["unk_cols"]) if isinstance(assay["unk_cols"], str) else assay["unk_cols"]
            unk_rows = parse_rows(assay.get("unk_rows", "A-H")) if isinstance(assay.get("unk_rows", "A-H"), str) else assay.get("unk_rows")
            unk_mask = df["Col"].isin(unk_cols) & df["Row"].isin(unk_rows) & df["Ct"].notna()
            unknowns = [(r["Well"], r["Ct"]) for _, r in df[unk_mask].iterrows()]
            if unknowns:
                print_unknowns(name, unknowns, slope, intercept, ct_min, ct_max)
 
        # Dilution series
        dilution_results = []
        for ds_def in assay.get("dilutions", []):
            ds_cols    = parse_cols(ds_def["cols"]) if isinstance(ds_def["cols"], str) else ds_def["cols"]
            ds_rows    = parse_rows(ds_def["rows"]) if isinstance(ds_def["rows"], str) else ds_def["rows"]
            vol_ml     = float(ds_def["sample_volume_ml"])
            dil_factor = float(ds_def.get("dilution_factor", 10))
            sname      = ds_def.get("sample_name", "Sample")
 
            rows_data, summary = calc_dilution_series(
                df, ds_cols, ds_rows, slope, intercept,
                ct_min, ct_max, vol_ml, dil_factor
            )
            print_dilution_series(f"{name} / {sname}", rows_data, summary,
                                  vol_ml, dil_factor)
            dilution_results.append({
                "sample_name":      sname,
                "cols":             ds_cols,
                "rows":             ds_rows,
                "sample_volume_ml": vol_ml,
                "dilution_factor":  dil_factor,
                "rows_data":        rows_data,
                "summary":          summary,
            })
 
        assay_results.append({
            "name":            name,
            "slope":           slope,
            "intercept":       intercept,
            "r2":              r2,
            "efficiency":      efficiency,
            "ct_pairs":        ct_pairs,
            "ct_min":          ct_min,
            "ct_max":          ct_max,
            "unknowns":        unknowns,
            "dilution_series": dilution_results,
        })
 
    if assay_results:
        export_tsv(filepath, plate, assay_results)
