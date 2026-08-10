#!/usr/bin/env python3

import json
import os

# ============================================================
# Mass points to create
# (MHc, MA)
# ============================================================
MASS_POINTS = [
    (130, 18),
    (130, 25),
    (130, 36),
    (130, 50),
    (130, 70),
    (130, 110),
    (130, 125)
    # 원하는 mass point 계속 추가
]

# ============================================================
# Output directory
# ============================================================
OUTPUT_DIR = (
    "/data6/Users/eunsu/SKNanoAnalyzer/"
    "data/Run3_v15_Run2_v15/2024/Sample/ForSNU"
)

for mhc, ma in MASS_POINTS:

    sample_name = f"TTToHcToWAToBB-MHc{mhc}_MA{ma}_SingleLepFilter"

    output_file = os.path.join(
        OUTPUT_DIR,
        f"{sample_name}.json"
    )

    sample = {
        "name": sample_name,
        "isMC": 1,
        "PD": (
            f"TTToHcToWA_AToBB_MHc-{mhc}_MA-{ma}"
            "_SingleLepFilter_TuneCP5_cff"
        ),
        "xsec": -1,
        "sumsign": -1,
        "sumW": -1,
        "nmc": -1,
        "path": []
    }

    with open(output_file, "w") as f:
        json.dump(sample, f, indent=4)

    print(f"Created: {output_file}")
