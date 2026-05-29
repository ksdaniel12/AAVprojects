from qpcr_functions import run_analysis


FILE_PATH = r"C:\Users\ksdan\OneDrive\Desktop\AAV20260526KD_data.xlsx"   # Windows

# ── Define your assays ─────────────────────────────────────────────────────
#
# Each assay needs:
#   name      – a label shown in the output
#   std_cols  – the two duplicate standard columns  (e.g. "1,2")
#   std_rows  – which rows have standards           (e.g. "A-H")
#   unk_cols  – columns containing your unknowns    (e.g. "5,6")
#   unk_rows  – which rows to include as unknowns   (e.g. "B-H")
#
# Add or remove dicts to match your plate layout.

ASSAYS = [
    {
        "name":     "Assay1",
        "std_cols": "3,4",     # columns 1 and 2 are the duplicate standard curve
        "std_rows": "A-H",     # rows A through H
        "unk_cols": "7,8",     # columns 5 and 6 are the unknowns
        "unk_rows": "A-H",     # rows A through H 
    },
    # Add a second assay by uncommenting and filling in below:
    # {
    #     "name":     "Assay2",
    #     "std_cols": "3,4",
    #     "std_rows": "A-H",
    #     "unk_cols": "7,8",
    #     "unk_rows": "B-H",
    # },
]

# ── Run ────────────────────────────────────────────────────────────────────

run_analysis(FILE_PATH, ASSAYS)
