import os
import json

def make_file_paths():
    base = "/gv0/DATA/SKNano/NanoAODv15/2024/WtoLNu-4Jets_Bin-4J_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_MC_2024_WtoLNu-4Jets_Bin-4J_TuneCP5_13p6TeV_madgraphMLM-pythia8/260213_202800"
    paths = []

    ranges = [
        ("0000", 1, 644)
    ]

    for subdir, start, end in ranges:
        for i in range(start, end + 1):
            paths.append(f"{base}/{subdir}/tree_{i}.root")

    return paths

def main():
    out_dir = "data/Run3_v15_Run2_v15/2024/Sample/ForSNU"
    os.makedirs(out_dir, exist_ok=True)

    payload = {
        "name": "WtoLNu_4J",
        "isMC": 1,
        "PD": "WtoLNu-4Jets_Bin-4J_TuneCP5_13p6TeV_madgraphMLM-pythia8",
        "xsec": 417.8,
        "sumsign": 85935724.0,
        "sumW": 622683209078.499,
        "nmc": 85935724.0,
        "path": make_file_paths(),
    }

    out_path = os.path.join(out_dir, "WtoLNu_4J.json")
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=4)

    print(f"Created: {out_path}")
    print(f"Number of ROOT files: {len(payload['path'])}")

if __name__ == "__main__":
    main()
