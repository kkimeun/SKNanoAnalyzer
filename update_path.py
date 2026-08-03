#!/usr/bin/env python3
import json
import glob
import os
import re

INPUT_JSON = "data/Run3_v15_Run2_v15/2024/Sample/ForSNU/TTToHcToWAToBB-MHc70_MA15_SingleLepFilter.json"

NANO_DIR = "/gv0/Users/eunsu/TTToHcToWA_AToBB_SampleProduction/TTToHcToWA_AToBB_MHc-70_MA-15_SingleLepFilter_TuneCP5_cff/chain_RunIII2024Summer24/20260727_233803"

def extract_number(path):
    name = os.path.basename(path)
    m = re.search(r"_(\d+)\.root$", name)
    return int(m.group(1)) if m else 999999999

root_files = glob.glob(os.path.join(NANO_DIR, "NanoAODv15_*.root"))
root_files = sorted(root_files, key=extract_number)

with open(INPUT_JSON, "r") as f:
    data = json.load(f)

data["path"].extend(root_files)

with open(INPUT_JSON, "w") as f:
    json.dump(data, f, indent=4)

print(f"Updated {INPUT_JSON}")
print(f"Number of root files: {len(root_files)}")
print("First file:", root_files[0] if root_files else "None")
print("Last file:", root_files[-1] if root_files else "None")
