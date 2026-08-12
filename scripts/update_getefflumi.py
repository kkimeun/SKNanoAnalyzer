#!/usr/bin/env python3

import os
import sys
import json
import glob
import shutil
import argparse
import tempfile

import ROOT


# ============================================================
# Default paths
# ============================================================

SKNANO_HOME = os.environ.get("SKNANO_HOME")

if not SKNANO_HOME:
    raise RuntimeError("SKNANO_HOME is not defined.")

ERA = "2024"

DEFAULT_GETEFFLUMI_DIR = (
    f"/data9/Users/eunsu/SKNanoOutput/GetEffLumi/{ERA}"
)

DEFAULT_SAMPLE_DIR = os.path.join(
    SKNANO_HOME,
    "data",
    "Run3_v15_Run2_v15",
    ERA,
    "Sample",
)

DEFAULT_COMMON_JSON = os.path.join(
    DEFAULT_SAMPLE_DIR,
    "CommonSampleInfo.json",
)

DEFAULT_FORSNU_DIR = os.path.join(
    DEFAULT_SAMPLE_DIR,
    "ForSNU",
)


# ============================================================
# ROOT helpers
# ============================================================

def find_root_object(directory, target_name):
    """
    Recursively search a ROOT TFile/TDirectory for an object
    whose name exactly matches target_name.

    This works whether GetEffLumi stores histograms:
        sumW
    directly at file level, or inside:
        <sample>/sumW
    """

    keys = directory.GetListOfKeys()

    if not keys:
        return None

    for key in keys:
        name = key.GetName()
        obj = key.ReadObj()

        if name == target_name:
            return obj

        if obj.InheritsFrom("TDirectory"):
            found = find_root_object(obj, target_name)

            if found:
                return found

    return None


def get_hist_value(hist, name):
    """
    Get the accumulated value stored in bin 1.
    """
    if hist is None:
        raise RuntimeError(f"Histogram '{name}' was not found.")

    if not hist.InheritsFrom("TH1"):
        raise RuntimeError(
            f"Object '{name}' exists but is not a TH1 histogram."
        )

    return float(hist.GetBinContent(1))


def read_getefflumi(root_path):
    """
    Read:
        NEvents
        sumSign
        sumW

    from a GetEffLumi ROOT output.

    nmc is primarily taken from NEvents bin content.
    A consistency check with sumW.GetEntries() is also printed.
    """

    if not os.path.isfile(root_path):
        raise FileNotFoundError(
            f"GetEffLumi ROOT file does not exist:\n"
            f"  {root_path}"
        )

    f = ROOT.TFile.Open(root_path, "READ")

    if not f or f.IsZombie():
        raise RuntimeError(
            f"Could not open ROOT file:\n"
            f"  {root_path}"
        )

    try:
        h_nevents = find_root_object(f, "NEvents")
        h_sumsign = find_root_object(f, "sumSign")
        h_sumw = find_root_object(f, "sumW")

        if h_sumsign is None:
            raise RuntimeError(
                f"'sumSign' histogram was not found in {root_path}"
            )

        if h_sumw is None:
            raise RuntimeError(
                f"'sumW' histogram was not found in {root_path}"
            )

        sumsign = get_hist_value(h_sumsign, "sumSign")
        sumw = get_hist_value(h_sumw, "sumW")

        # ----------------------------------------------------
        # nmc
        #
        # Prefer the explicit NEvents histogram.
        # If unavailable, reproduce the existing
        # sampleManager.py convention:
        #
        #   nmc = h_sumW.GetEntries()
        # ----------------------------------------------------
        if h_nevents is not None:
            nmc = get_hist_value(h_nevents, "NEvents")
            nmc_source = "NEvents bin content"
        else:
            nmc = float(h_sumw.GetEntries())
            nmc_source = "sumW.GetEntries() fallback"

        sumw_entries = float(h_sumw.GetEntries())
        sumsign_entries = float(h_sumsign.GetEntries())

        return {
            "nmc": nmc,
            "sumsign": sumsign,
            "sumW": sumw,
            "nmc_source": nmc_source,
            "sumW_entries": sumw_entries,
            "sumSign_entries": sumsign_entries,
        }

    finally:
        f.Close()


# ============================================================
# JSON helpers
# ============================================================

def load_json(path):
    with open(path, "r") as f:
        return json.load(f)


def write_json_atomic(path, data):
    """
    Atomic JSON update:
    write a temporary file first, then replace original.
    """

    directory = os.path.dirname(path)

    fd, tmp_path = tempfile.mkstemp(
        prefix=".tmp_update_getefflumi_",
        suffix=".json",
        dir=directory,
        text=True,
    )

    try:
        with os.fdopen(fd, "w") as f:
            json.dump(
                data,
                f,
                indent=4,
                ensure_ascii=False,
            )
            f.write("\n")

        os.replace(tmp_path, path)

    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def backup_file(path):
    """
    Create <file>.bak before modifying it.
    Existing .bak is overwritten so that it always represents
    the state immediately before the current update.
    """

    backup_path = path + ".bak"
    shutil.copy2(path, backup_path)

    return backup_path


# ============================================================
# Find corresponding ForSNU JSON
# ============================================================

def find_forsnu_json(sample_name, forsnu_dir):
    """
    Search recursively under ForSNU.

    Expected example:
      ForSNU/TTToHcToWAToBB-MHc130_MA/
        TTToHcToWAToBB-MHc130_MA18_SingleLepFilter.json
    """

    exact = os.path.join(
        forsnu_dir,
        "**",
        sample_name + ".json",
    )

    matches = glob.glob(
        exact,
        recursive=True,
    )

    matches = [
        x for x in matches
        if os.path.isfile(x)
    ]

    if len(matches) == 0:
        raise FileNotFoundError(
            f"No ForSNU JSON found for sample:\n"
            f"  {sample_name}\n"
            f"searched under:\n"
            f"  {forsnu_dir}"
        )

    if len(matches) > 1:
        raise RuntimeError(
            f"Multiple ForSNU JSON files found for '{sample_name}':\n"
            + "\n".join(f"  {x}" for x in matches)
        )

    return matches[0]


# ============================================================
# Validation
# ============================================================

def validate_sample_name(sample_name):
    """
    Safety guard:
    this script is intentionally restricted to
    TTToHcToWAToBB signal samples.
    """

    if not sample_name.startswith("TTToHcToWAToBB-"):
        raise RuntimeError(
            "Safety check failed.\n"
            "This script only updates samples beginning with:\n"
            "  TTToHcToWAToBB-\n"
            f"Received:\n"
            f"  {sample_name}"
        )


def validate_forsnu_json(sample_name, data, path):
    json_name = data.get("name")

    if json_name != sample_name:
        raise RuntimeError(
            f"ForSNU JSON name mismatch:\n"
            f"  requested : {sample_name}\n"
            f"  JSON name : {json_name}\n"
            f"  file      : {path}"
        )


# ============================================================
# One sample update
# ============================================================

def update_one_sample(
    sample_name,
    getefflumi_dir,
    common_json_path,
    forsnu_dir,
    dry_run=False,
    make_backup=True,
):
    validate_sample_name(sample_name)

    root_path = os.path.join(
        getefflumi_dir,
        sample_name + ".root",
    )

    print()
    print("=" * 80)
    print("Sample:", sample_name)
    print("=" * 80)

    # --------------------------------------------------------
    # 1. Read GetEffLumi ROOT
    # --------------------------------------------------------
    values = read_getefflumi(root_path)

    nmc = values["nmc"]
    sumsign = values["sumsign"]
    sumw = values["sumW"]

    print("GetEffLumi ROOT:")
    print(" ", root_path)

    print()
    print("Read values:")
    print(f"  nmc     = {nmc}")
    print(f"  sumsign = {sumsign}")
    print(f"  sumW    = {sumw}")

    print()
    print("Consistency information:")
    print(f"  nmc source          = {values['nmc_source']}")
    print(f"  sumW entries        = {values['sumW_entries']}")
    print(f"  sumSign entries     = {values['sumSign_entries']}")

    if abs(nmc - values["sumW_entries"]) > 0.5:
        print()
        print(
            "[WARNING] NEvents and sumW.GetEntries() differ:"
        )
        print(f"  NEvents      = {nmc}")
        print(f"  sumW entries = {values['sumW_entries']}")

    # --------------------------------------------------------
    # 2. Load CommonSampleInfo
    # --------------------------------------------------------
    common = load_json(common_json_path)

    if sample_name not in common:
        raise KeyError(
            f"Sample does not exist in CommonSampleInfo.json:\n"
            f"  {sample_name}"
        )

    common_entry = common[sample_name]

    # --------------------------------------------------------
    # 3. Find/load matching ForSNU JSON
    # --------------------------------------------------------
    forsnu_json_path = find_forsnu_json(
        sample_name,
        forsnu_dir,
    )

    forsnu = load_json(forsnu_json_path)

    validate_forsnu_json(
        sample_name,
        forsnu,
        forsnu_json_path,
    )

    # --------------------------------------------------------
    # 4. Show old -> new values
    # --------------------------------------------------------
    print()
    print("CommonSampleInfo.json:")
    print(f"  nmc     : {common_entry.get('nmc')} -> {nmc}")
    print(f"  sumsign : {common_entry.get('sumsign')} -> {sumsign}")
    print(f"  sumW    : {common_entry.get('sumW')} -> {sumw}")

    print()
    print("ForSNU JSON:")
    print(" ", forsnu_json_path)
    print(f"  nmc     : {forsnu.get('nmc')} -> {nmc}")
    print(f"  sumsign : {forsnu.get('sumsign')} -> {sumsign}")
    print(f"  sumW    : {forsnu.get('sumW')} -> {sumw}")

    # --------------------------------------------------------
    # 5. Modify ONLY these three fields.
    #
    # xsec, filterEff, PD, path, etc. are untouched.
    # --------------------------------------------------------
    common_entry["nmc"] = nmc
    common_entry["sumsign"] = sumsign
    common_entry["sumW"] = sumw

    forsnu["nmc"] = nmc
    forsnu["sumsign"] = sumsign
    forsnu["sumW"] = sumw

    # --------------------------------------------------------
    # 6. Dry run
    # --------------------------------------------------------
    if dry_run:
        print()
        print("[DRY RUN] No JSON files were modified.")
        return

    # --------------------------------------------------------
    # 7. Backups
    # --------------------------------------------------------
    if make_backup:
        common_backup = backup_file(common_json_path)
        forsnu_backup = backup_file(forsnu_json_path)

        print()
        print("Backups:")
        print(" ", common_backup)
        print(" ", forsnu_backup)

    # --------------------------------------------------------
    # 8. Write both JSONs atomically
    # --------------------------------------------------------
    write_json_atomic(
        common_json_path,
        common,
    )

    write_json_atomic(
        forsnu_json_path,
        forsnu,
    )

    print()
    print("[UPDATED]")
    print(" ", common_json_path)
    print(" ", forsnu_json_path)


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Update nmc, sumsign and sumW for TTToHcToWAToBB "
            "samples using GetEffLumi ROOT outputs."
        )
    )

    parser.add_argument(
        "samples",
        nargs="+",
        help=(
            "Sample name(s), e.g. "
            "TTToHcToWAToBB-MHc130_MA18_SingleLepFilter"
        ),
    )

    parser.add_argument(
        "--getefflumi-dir",
        default=DEFAULT_GETEFFLUMI_DIR,
        help=(
            "Directory containing GetEffLumi ROOT outputs "
            f"(default: {DEFAULT_GETEFFLUMI_DIR})"
        ),
    )

    parser.add_argument(
        "--common-json",
        default=DEFAULT_COMMON_JSON,
        help=(
            "CommonSampleInfo.json path "
            f"(default: {DEFAULT_COMMON_JSON})"
        ),
    )

    parser.add_argument(
        "--forsnu-dir",
        default=DEFAULT_FORSNU_DIR,
        help=(
            "ForSNU directory "
            f"(default: {DEFAULT_FORSNU_DIR})"
        ),
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print changes without modifying JSON files.",
    )

    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not create .bak copies before updating.",
    )

    args = parser.parse_args()

    if not os.path.isfile(args.common_json):
        raise FileNotFoundError(
            f"CommonSampleInfo.json not found:\n"
            f"  {args.common_json}"
        )

    if not os.path.isdir(args.forsnu_dir):
        raise NotADirectoryError(
            f"ForSNU directory not found:\n"
            f"  {args.forsnu_dir}"
        )

    if not os.path.isdir(args.getefflumi_dir):
        raise NotADirectoryError(
            f"GetEffLumi directory not found:\n"
            f"  {args.getefflumi_dir}"
        )

    # --------------------------------------------------------
    # Update only explicitly requested samples.
    # No automatic scanning/modification of unrelated samples.
    # --------------------------------------------------------
    for sample_name in args.samples:
        update_one_sample(
            sample_name=sample_name,
            getefflumi_dir=args.getefflumi_dir,
            common_json_path=args.common_json,
            forsnu_dir=args.forsnu_dir,
            dry_run=args.dry_run,
            make_backup=not args.no_backup,
        )

    print()
    print("=" * 80)

    if args.dry_run:
        print("Finished dry run.")
    else:
        print("Finished updating requested sample(s).")

    print("=" * 80)


if __name__ == "__main__":
    main()
