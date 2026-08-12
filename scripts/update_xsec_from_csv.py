#!/usr/bin/env python3

import os
import re
import csv
import json
import glob
import shutil
import argparse
import tempfile


# ============================================================
# Default paths
# ============================================================

SKNANO_HOME = os.environ.get("SKNANO_HOME")

if not SKNANO_HOME:
    raise RuntimeError("SKNANO_HOME is not defined.")

ERA = "2024"

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
# JSON helpers
# ============================================================

def load_json(path):
    with open(path, "r") as f:
        return json.load(f)


def write_json_atomic(path, data):
    """
    Write JSON safely via temporary file + atomic replacement.
    """

    directory = os.path.dirname(path)

    fd, tmp_path = tempfile.mkstemp(
        prefix=".tmp_update_xsec_",
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
    Save original file as <filename>.bak.
    Existing .bak is overwritten.
    """

    backup_path = path + ".bak"
    shutil.copy2(path, backup_path)
    return backup_path


# ============================================================
# Sample-name parsing
# ============================================================

def parse_mass_point(sample_name):
    """
    Parse:

      TTToHcToWAToBB-MHc130_MA18_SingleLepFilter

    -> MHc = 130
       MA  = 18

    Also works for samples without _SingleLepFilter.
    """

    pattern = r"MHc(\d+)_MA(\d+)"

    match = re.search(pattern, sample_name)

    if not match:
        raise RuntimeError(
            f"Could not parse MHc/MA from sample name:\n"
            f"  {sample_name}"
        )

    mhc = int(match.group(1))
    ma = int(match.group(2))

    return mhc, ma


def validate_sample_name(sample_name):
    """
    Safety guard: only TTToHcToWAToBB samples are allowed.
    """

    if not sample_name.startswith("TTToHcToWAToBB-"):
        raise RuntimeError(
            "Safety check failed.\n"
            "This script only updates samples beginning with:\n"
            "  TTToHcToWAToBB-\n"
            f"Received:\n"
            f"  {sample_name}"
        )


# ============================================================
# CSV reader
# ============================================================

def read_xsec_csv(csv_path):
    """
    Read gridpack cross sections from a CSV.

    Supported format:

        MA,xsec_pb,error_pb
        18,8171,15.29
        25,...

    or simply:

        MA,xsec_pb
        18,8171
        25,...

    Optional more-general format:

        MHc,MA,xsec_pb,error_pb
        130,18,8171,15.29
        ...

    Returns:
        dictionary keyed by either
          (MHc, MA) when MHc column exists
        or
          MA when only MA exists.
    """

    if not os.path.isfile(csv_path):
        raise FileNotFoundError(
            f"xsec CSV does not exist:\n"
            f"  {csv_path}"
        )

    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)

        if not reader.fieldnames:
            raise RuntimeError(
                f"CSV has no header:\n"
                f"  {csv_path}"
            )

        required = {"MA", "xsec_pb"}

        missing = required - set(reader.fieldnames)

        if missing:
            raise RuntimeError(
                f"CSV is missing required column(s): {sorted(missing)}\n"
                f"Found columns: {reader.fieldnames}"
            )

        has_mhc = "MHc" in reader.fieldnames

        result = {}

        for row in reader:

            # Skip blank rows
            if not row.get("MA") or not row.get("xsec_pb"):
                continue

            ma = int(row["MA"])
            xsec = float(row["xsec_pb"])

            if xsec <= 0:
                raise RuntimeError(
                    f"Invalid non-positive xsec in CSV:\n"
                    f"  MA={ma}, xsec={xsec}"
                )

            if has_mhc:
                mhc = int(row["MHc"])
                key = (mhc, ma)
            else:
                key = ma

            if key in result:
                raise RuntimeError(
                    f"Duplicate mass point in CSV:\n"
                    f"  {key}"
                )

            result[key] = xsec

    if not result:
        raise RuntimeError(
            f"No valid xsec entries were read from:\n"
            f"  {csv_path}"
        )

    return result, has_mhc


# ============================================================
# Find ForSNU JSON
# ============================================================

def find_forsnu_json(sample_name, forsnu_dir):

    pattern = os.path.join(
        forsnu_dir,
        "**",
        sample_name + ".json",
    )

    matches = [
        path
        for path in glob.glob(pattern, recursive=True)
        if os.path.isfile(path)
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
            + "\n".join(f"  {path}" for path in matches)
        )

    return matches[0]


# ============================================================
# Numerical formatting
# ============================================================

def format_number(value):
    """
    Make xsec_formula readable while retaining useful precision.

    Examples:
      8171.0       -> "8171"
      0.665814     -> "0.665814"
      8170.912345  -> "8170.912345"
    """

    return f"{float(value):.12g}"


# ============================================================
# One-sample update
# ============================================================

def update_one_sample(
    sample_name,
    xsec_map,
    csv_has_mhc,
    common_json_path,
    forsnu_dir,
    dry_run=False,
):
    validate_sample_name(sample_name)

    mhc, ma = parse_mass_point(sample_name)

    print()
    print("=" * 80)
    print("Sample:", sample_name)
    print("=" * 80)

    # --------------------------------------------------------
    # Load CommonSampleInfo.json
    # --------------------------------------------------------

    common = load_json(common_json_path)

    if sample_name not in common:
        raise KeyError(
            f"Sample does not exist in CommonSampleInfo.json:\n"
            f"  {sample_name}"
        )

    common_entry = common[sample_name]

    # --------------------------------------------------------
    # Read filter efficiency
    # --------------------------------------------------------

    if "filterEff" not in common_entry:
        raise KeyError(
            f"'filterEff' is missing from CommonSampleInfo.json for:\n"
            f"  {sample_name}"
        )

    filter_eff = float(common_entry["filterEff"])

    if not (0.0 < filter_eff <= 1.0):
        raise RuntimeError(
            f"Invalid filterEff for {sample_name}:\n"
            f"  filterEff = {filter_eff}"
        )

    # --------------------------------------------------------
    # Find gridpack xsec in CSV
    # --------------------------------------------------------

    if csv_has_mhc:
        key = (mhc, ma)
    else:
        key = ma

    if key not in xsec_map:
        raise KeyError(
            f"No gridpack xsec found in CSV for:\n"
            f"  MHc={mhc}, MA={ma}"
        )

    gridpack_xsec = float(xsec_map[key])

    # --------------------------------------------------------
    # Calculate effective filtered xsec
    # --------------------------------------------------------

    effective_xsec = gridpack_xsec * filter_eff

    xsec_formula = (
        f"{format_number(gridpack_xsec)}"
        f"*{format_number(filter_eff)}"
    )

    print(f"MHc                   = {mhc}")
    print(f"MA                    = {ma}")
    print(f"Gridpack xsec [pb]    = {gridpack_xsec}")
    print(f"filterEff             = {filter_eff}")
    print(f"xsec_formula          = {xsec_formula}")
    print(f"Effective xsec [pb]   = {effective_xsec}")

    # --------------------------------------------------------
    # Find corresponding ForSNU JSON
    # --------------------------------------------------------

    forsnu_json_path = find_forsnu_json(
        sample_name,
        forsnu_dir,
    )

    forsnu = load_json(forsnu_json_path)

    json_name = forsnu.get("name")

    if json_name != sample_name:
        raise RuntimeError(
            f"ForSNU JSON name mismatch:\n"
            f"  requested : {sample_name}\n"
            f"  JSON name : {json_name}\n"
            f"  file      : {forsnu_json_path}"
        )

    # --------------------------------------------------------
    # Display old -> new
    # --------------------------------------------------------

    print()
    print("CommonSampleInfo.json:")
    print(
        f"  xsec_formula : "
        f"{common_entry.get('xsec_formula', '<missing>')} "
        f"-> {xsec_formula}"
    )
    print(
        f"  xsec         : "
        f"{common_entry.get('xsec')} "
        f"-> {effective_xsec}"
    )

    print()
    print("ForSNU JSON:")
    print(" ", forsnu_json_path)
    print(
        f"  xsec         : "
        f"{forsnu.get('xsec')} "
        f"-> {effective_xsec}"
    )

    # --------------------------------------------------------
    # Dry run
    # --------------------------------------------------------

    if dry_run:
        print()
        print("[DRY RUN] No files modified.")
        return

    # --------------------------------------------------------
    # Modify ONLY xsec and xsec_formula
    #
    # filterEff, sumW, sumsign, nmc, PD, path, etc.
    # are intentionally untouched.
    # --------------------------------------------------------

    common_entry["xsec_formula"] = xsec_formula
    common_entry["xsec"] = effective_xsec

    forsnu["xsec"] = effective_xsec

    # --------------------------------------------------------
    # Backup
    # --------------------------------------------------------

    common_backup = backup_file(common_json_path)
    forsnu_backup = backup_file(forsnu_json_path)

    print()
    print("Backups:")
    print(" ", common_backup)
    print(" ", forsnu_backup)

    # --------------------------------------------------------
    # Write JSON
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
            "Read gridpack cross sections from CSV, multiply by "
            "filterEff from CommonSampleInfo.json, and update "
            "xsec/xsec_formula in CommonSampleInfo and ForSNU JSON."
        )
    )

    parser.add_argument(
        "samples",
        nargs="+",
        help=(
            "Sample name(s), e.g.\n"
            "TTToHcToWAToBB-MHc130_MA18_SingleLepFilter"
        ),
    )

    parser.add_argument(
        "--csv",
        required=True,
        help=(
            "CSV containing gridpack cross sections. "
            "Required columns: MA,xsec_pb "
            "(optional MHc,error_pb allowed)."
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
        help="Show proposed changes without modifying files.",
    )

    args = parser.parse_args()

    # --------------------------------------------------------
    # Validate input paths
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # Read CSV once
    # --------------------------------------------------------

    xsec_map, csv_has_mhc = read_xsec_csv(args.csv)

    print("Loaded xsec CSV:")
    print(" ", args.csv)
    print("Number of mass points:", len(xsec_map))

    # --------------------------------------------------------
    # Explicitly requested samples only
    # --------------------------------------------------------

    for sample_name in args.samples:
        update_one_sample(
            sample_name=sample_name,
            xsec_map=xsec_map,
            csv_has_mhc=csv_has_mhc,
            common_json_path=args.common_json,
            forsnu_dir=args.forsnu_dir,
            dry_run=args.dry_run,
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
