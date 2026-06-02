from qpcr_functions import run_analysis

# Path to file
FILE_PATH = r"C:\Users\ksdan\OneDrive\Desktop\results\AAV20260526KD_data.xlsx"


# Define assays
#
# Each assay needs:
#   name      – a label shown in the output
#   std_cols  – the two duplicate standard columns       e.g. "1,2"
#   std_rows  – which rows contain standards             e.g. "A-H"
#
# Then add whichever unknown types apply:
#
# OPTION A — plain unknowns (just gives copies per reaction, no back-calculation):
#   unk_cols  – columns containing unknowns              e.g. "5,6"
#   unk_rows  – which rows to include                    e.g. "A-H"
#
# OPTION B — dilution series (back-calculates copies/mL in original sample):
#   dilutions – list of samples, each with:
#     sample_name      – label for this sample
#     cols             – the two duplicate columns       e.g. "7,8"
#     rows             – rows in dilution order          e.g. "A-H"
#     sample_volume_ml – mL of original sample in the first (most concentrated) row
#     dilution_factor  – fold dilution between each row  (default: 10)
#

ASSAYS = [
    {
        "name":     "Primer1",
        "std_cols": "1,2",      # columns are the duplicate standard curve
        "std_rows": "A-H",

        # Plain unknowns — just reports copies per reaction
        "unk_cols": "5,6",
        "unk_rows": "A-H",

        # Dilution series — back-calculates copies/mL in the original sample
        "dilutions": [
            {
                "sample_name":      "SampleA",
                "cols":             "5,6",   # duplicate columns for this sample
                "rows":             "A-H",   # rows A-H are the 1:10 serial dilution
                "sample_volume_ml": 0.01,    # mL of original sample added to row A
                "dilution_factor":  10,      # each row is 10x more dilute than the last
            },
            # Add more samples by copying the block above:
            # {
            #     "sample_name":      "SampleB",
            #     "cols":             "9,10",
            #     "rows":             "A-H",
            #     "sample_volume_ml": 0.01,
            #     "dilution_factor":  10,
            # },
        ],
    },

    # Add a second assay (new primers, etc) by uncommenting below
     {
         "name":     "Primer2",
         "std_cols": "3,4",
         "std_rows": "A-H",
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
]

# Run

run_analysis(FILE_PATH, ASSAYS)

