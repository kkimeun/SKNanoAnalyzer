#!/usr/bin/env python3
import os
import json

era = "2024"

dirs = [
    "/gv0/Users/eunsu/SKNano/2024/store/mc/RunIII2024Summer24NanoAODv15/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/2560000",
    "/gv0/Users/eunsu/SKNano/2024/store/mc/RunIII2024Summer24NanoAODv15/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/2810000",
]

paths = []

for d in dirs:
    print("Reading:", d, flush=True)
    for f in os.listdir(d):
        if f.endswith(".root"):
            paths.append(os.path.join(d, f))

paths = sorted(paths)

out = os.path.join(
    os.environ["SKNANO_DATA"],
    era,
    "Sample",
    "ForSNU",
    "ttZtoQQ.json"
)

data = {
    "name": "ttZtoQQ",
    "isMC": 1,
    "PD": "TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8",
    "xsec": 0.6603,
    "sumsign": 35753.0,
    "sumW": 115925.54177975655,
    "nmc": 90817.0,
    "path": paths
}

with open(out, "w") as f:
    json.dump(data, f, indent=4)

print("DONE", flush=True)
print("file:", out, flush=True)
print("nfiles:", len(paths), flush=True)
