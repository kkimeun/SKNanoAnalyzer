#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


DEFAULT_INPUT_DIR = Path(
    "/data9/Users/eunsu/SKNanoOutput/Tutorial_reco_tt/2024"
)

DEFAULT_OUTPUT_DIR = Path(
    "/data9/Users/eunsu/SKNanoOutput/Tutorial_reco_tt/2024/"
    "single_lepton_filter_pt_comparison"
)

MASS_POINTS = [
    ("MHc70_MA15", "70", "15"),
    ("MHc100_MA60", "100", "60"),
    ("MHc130_MA90", "130", "90"),
    ("MHc160_MA155", "160", "155"),
]

CHANNELS = {
    "Muon": {
        "hist_suffix": "muon_pt0",
        "display": r"$\mu$ channel",
    },
    "Electron": {
        "hist_suffix": "electron_pt0",
        "display": r"$e$ channel",
    },
}

STAGES = {
    "recoBeforeTrigger": "Before trigger",
    "recoAfterTriggerBeforePtCut": "After trigger, before safe-$p_T$ cut",
    "recoAfterTrigger": "After trigger and safe-$p_T$ cut",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare reco lepton pT before/after SingleLepFilter."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Directory containing Tutorial_reco_tt ROOT outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory in which plots and CSV are saved.",
    )
    parser.add_argument(
        "--xmax",
        type=float,
        default=250.0,
        help="Maximum lepton pT shown in the plots.",
    )
    parser.add_argument(
        "--logy",
        action="store_true",
        help="Use logarithmic y axis in the main distribution panels.",
    )
    return parser.parse_args()


def root_file_paths(
    input_dir: Path,
    mass_tag: str,
) -> tuple[Path, Path]:
    unfiltered = input_dir / f"TTToHcToWAToBB-{mass_tag}.root"
    filtered = input_dir / f"TTToHcToWAToBB-{mass_tag}_SingleLepFilter.root"
    return unfiltered, filtered


def find_histogram(
    root_file: uproot.ReadOnlyDirectory,
    channel: str,
    stage: str,
):
    """
    Find a TH1 below, for example:
      Muon/recoBeforeTrigger_muon_pt0/histo_...
    """
    suffix = CHANNELS[channel]["hist_suffix"]
    expected_directory = f"{channel}/{stage}_{suffix}/"

    candidates: list[str] = []

    for key in root_file.keys(recursive=True, cycle=False):
        if key.startswith(expected_directory):
            candidates.append(key)

    if not candidates:
        available = [
            key
            for key in root_file.keys(recursive=True, cycle=False)
            if key.startswith(f"{channel}/")
        ]
        raise KeyError(
            f"No histogram found below '{expected_directory}'.\n"
            f"Available {channel} keys include:\n"
            + "\n".join(f"  {key}" for key in available[:30])
        )

    # There should normally be exactly one histogram in the directory.
    for key in candidates:
        obj = root_file[key]
        if obj.classname.startswith("TH1"):
            return obj, key

    raise TypeError(
        f"Objects found below '{expected_directory}', but none was TH1: "
        f"{candidates}"
    )


def read_histogram(
    file_path: Path,
    channel: str,
    stage: str,
) -> dict:
    if not file_path.is_file():
        raise FileNotFoundError(f"Missing ROOT file: {file_path}")

    with uproot.open(file_path) as root_file:
        hist, hist_key = find_histogram(root_file, channel, stage)

        counts, edges = hist.to_numpy(flow=False)
        counts = np.asarray(counts, dtype=float)
        edges = np.asarray(edges, dtype=float)

        centers = 0.5 * (edges[:-1] + edges[1:])
        widths = np.diff(edges)

        integral = float(np.sum(counts))

        # ROOT fEntries is the number of Fill calls for an ordinary,
        # unit-weight histogram. Fall back to the integral if unavailable.
        try:
            entries = float(hist.member("fEntries"))
        except Exception:
            entries = integral

        if integral > 0:
            mean = float(np.sum(counts * centers) / integral)
            variance = float(
                np.sum(counts * (centers - mean) ** 2) / integral
            )
            stddev = float(np.sqrt(max(variance, 0.0)))

            # Density normalized so that integral over pT is one.
            density = counts / (integral * widths)
        else:
            mean = np.nan
            stddev = np.nan
            density = np.zeros_like(counts)

        return {
            "counts": counts,
            "density": density,
            "edges": edges,
            "centers": centers,
            "widths": widths,
            "integral": integral,
            "entries": entries,
            "mean": mean,
            "stddev": stddev,
            "key": hist_key,
        }


def step_values(values: np.ndarray) -> np.ndarray:
    """Make bin values compatible with plt.stairs."""
    return np.asarray(values, dtype=float)


def shape_distance(
    density_a: np.ndarray,
    density_b: np.ndarray,
    widths: np.ndarray,
) -> float:
    """
    Total variation distance between two normalized histograms.
    0 means identical; 1 means completely non-overlapping.
    """
    return float(0.5 * np.sum(np.abs(density_a - density_b) * widths))


def draw_individual_plot(
    no_filter: dict,
    filtered: dict,
    channel: str,
    stage: str,
    mass_label: str,
    output_path: Path,
    xmax: float,
    logy: bool,
) -> None:
    fig = plt.figure(figsize=(8, 7))
    grid = fig.add_gridspec(
        2,
        1,
        height_ratios=(3.2, 1.0),
        hspace=0.05,
    )

    ax = fig.add_subplot(grid[0])
    ratio_ax = fig.add_subplot(grid[1], sharex=ax)

    edges = no_filter["edges"]

    ax.stairs(
        no_filter["density"],
        edges,
        linewidth=1.8,
        label=(
            "No filter\n"
            f"Entries={no_filter['entries']:.0f}, "
            rf"$\langle p_T\rangle$={no_filter['mean']:.2f} GeV"
        ),
    )

    ax.stairs(
        filtered["density"],
        edges,
        linewidth=1.8,
        label=(
            "SingleLepFilter\n"
            f"Entries={filtered['entries']:.0f}, "
            rf"$\langle p_T\rangle$={filtered['mean']:.2f} GeV"
        ),
    )

    ax.set_ylabel(r"Normalized events / GeV")
    ax.set_title(
        f"{mass_label} · {CHANNELS[channel]['display']}\n"
        f"{STAGES[stage]}"
    )
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    if logy:
        ax.set_yscale("log")
        positive = np.concatenate(
            [
                no_filter["density"][no_filter["density"] > 0],
                filtered["density"][filtered["density"] > 0],
            ]
        )
        if positive.size:
            ax.set_ylim(bottom=max(np.min(positive) * 0.4, 1e-7))

    ratio = np.full_like(filtered["density"], np.nan)
    valid = no_filter["density"] > 0
    ratio[valid] = (
        filtered["density"][valid]
        / no_filter["density"][valid]
    )

    ratio_ax.stairs(
        ratio,
        edges,
        linewidth=1.5,
    )
    ratio_ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ratio_ax.set_ylabel("Filtered /\nNo filter")
    ratio_ax.set_xlabel(r"Reco lepton $p_T$ [GeV]")
    ratio_ax.set_ylim(0.5, 1.5)
    ratio_ax.grid(alpha=0.25)

    ax.set_xlim(0.0, xmax)

    tv_distance = shape_distance(
        no_filter["density"],
        filtered["density"],
        no_filter["widths"],
    )

    raw_ratio = (
        filtered["integral"] / no_filter["integral"]
        if no_filter["integral"] > 0
        else np.nan
    )

    ax.text(
        0.97,
        0.58,
        "\n".join(
            [
                rf"$\Delta\langle p_T\rangle$ = "
                f"{filtered['mean'] - no_filter['mean']:+.2f} GeV",
                f"Shape distance = {tv_distance:.3f}",
                f"Raw integral ratio = {raw_ratio:.3f}",
            ]
        ),
        transform=ax.transAxes,
        horizontalalignment="right",
        verticalalignment="top",
        fontsize=9,
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "alpha": 0.8,
        },
    )

    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def draw_masspoint_summary(
    all_results: dict,
    mass_tag: str,
    mass_label: str,
    output_path: Path,
    xmax: float,
    logy: bool,
) -> None:
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(16, 9),
        sharex=True,
    )

    for row, channel in enumerate(CHANNELS):
        for column, stage in enumerate(STAGES):
            ax = axes[row, column]
            result = all_results[(mass_tag, channel, stage)]

            no_filter = result["no_filter"]
            filtered = result["filtered"]

            ax.stairs(
                no_filter["density"],
                no_filter["edges"],
                linewidth=1.6,
                label="No filter",
            )
            ax.stairs(
                filtered["density"],
                filtered["edges"],
                linewidth=1.6,
                label="SingleLepFilter",
            )

            ax.set_xlim(0.0, xmax)
            ax.grid(alpha=0.25)

            if logy:
                ax.set_yscale("log")

            if row == 0:
                ax.set_title(STAGES[stage])

            if column == 0:
                ax.set_ylabel(
                    f"{CHANNELS[channel]['display']}\n"
                    "Normalized events / GeV"
                )

            if row == 1:
                ax.set_xlabel(r"Reco lepton $p_T$ [GeV]")

            delta_mean = filtered["mean"] - no_filter["mean"]
            tv_distance = shape_distance(
                no_filter["density"],
                filtered["density"],
                no_filter["widths"],
            )

            ax.text(
                0.97,
                0.95,
                (
                    rf"$\Delta\langle p_T\rangle$={delta_mean:+.2f} GeV"
                    "\n"
                    f"shape distance={tv_distance:.3f}"
                ),
                transform=ax.transAxes,
                horizontalalignment="right",
                verticalalignment="top",
                fontsize=8,
            )

    axes[0, 0].legend(fontsize=9)

    fig.suptitle(
        f"{mass_label}: reco lepton $p_T$, "
        "SingleLepFilter comparison",
        fontsize=15,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    individual_dir = output_dir / "individual_24"
    summary_dir = output_dir / "masspoint_summaries"

    individual_dir.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.is_dir():
        print(f"[ERROR] Input directory does not exist: {input_dir}")
        return 1

    results: dict = {}
    rows: list[dict] = []

    for mass_tag, mhc, ma in MASS_POINTS:
        no_filter_path, filtered_path = root_file_paths(
            input_dir,
            mass_tag,
        )

        mass_label = rf"$M_{{H^\pm}}={mhc}$ GeV, $M_A={ma}$ GeV"

        print(f"\n[Mass point] {mass_tag}")
        print(f"  no filter : {no_filter_path}")
        print(f"  filtered  : {filtered_path}")

        for channel in CHANNELS:
            for stage in STAGES:
                try:
                    no_filter = read_histogram(
                        no_filter_path,
                        channel,
                        stage,
                    )
                    filtered = read_histogram(
                        filtered_path,
                        channel,
                        stage,
                    )
                except Exception as exc:
                    print(
                        f"[ERROR] {mass_tag}, {channel}, {stage}: {exc}",
                        file=sys.stderr,
                    )
                    continue

                result = {
                    "no_filter": no_filter,
                    "filtered": filtered,
                }
                results[(mass_tag, channel, stage)] = result

                raw_integral_ratio = (
                    filtered["integral"] / no_filter["integral"]
                    if no_filter["integral"] > 0
                    else np.nan
                )

                tv_distance = shape_distance(
                    no_filter["density"],
                    filtered["density"],
                    no_filter["widths"],
                )

                rows.append(
                    {
                        "mass_point": mass_tag,
                        "MHc_GeV": mhc,
                        "MA_GeV": ma,
                        "channel": channel,
                        "stage": stage,
                        "no_filter_entries": no_filter["entries"],
                        "filtered_entries": filtered["entries"],
                        "no_filter_integral": no_filter["integral"],
                        "filtered_integral": filtered["integral"],
                        "raw_integral_ratio_filtered_over_no_filter":
                            raw_integral_ratio,
                        "no_filter_mean_pt_GeV": no_filter["mean"],
                        "filtered_mean_pt_GeV": filtered["mean"],
                        "delta_mean_pt_GeV":
                            filtered["mean"] - no_filter["mean"],
                        "shape_total_variation_distance": tv_distance,
                        "no_filter_hist_key": no_filter["key"],
                        "filtered_hist_key": filtered["key"],
                    }
                )

                output_name = (
                    f"{mass_tag}_{channel}_{stage}_"
                    "filtered_vs_unfiltered.png"
                )

                draw_individual_plot(
                    no_filter=no_filter,
                    filtered=filtered,
                    channel=channel,
                    stage=stage,
                    mass_label=mass_label,
                    output_path=individual_dir / output_name,
                    xmax=args.xmax,
                    logy=args.logy,
                )

                print(
                    f"  {channel:8s} {stage:35s} "
                    f"mean: {no_filter['mean']:.2f} -> "
                    f"{filtered['mean']:.2f} GeV, "
                    f"shape distance={tv_distance:.4f}"
                )

        expected_keys = [
            (mass_tag, channel, stage)
            for channel in CHANNELS
            for stage in STAGES
        ]

        if all(key in results for key in expected_keys):
            draw_masspoint_summary(
                all_results=results,
                mass_tag=mass_tag,
                mass_label=mass_label,
                output_path=(
                    summary_dir
                    / f"{mass_tag}_all_lepton_pt_comparison.png"
                ),
                xmax=args.xmax,
                logy=args.logy,
            )

    csv_path = output_dir / "lepton_pt_comparison_summary.csv"

    if rows:
        with csv_path.open("w", newline="") as output_file:
            writer = csv.DictWriter(
                output_file,
                fieldnames=list(rows[0].keys()),
            )
            writer.writeheader()
            writer.writerows(rows)

    print("\n[Done]")
    print(f"  Individual plots: {individual_dir}")
    print(f"  Mass-point summaries: {summary_dir}")
    print(f"  Numerical summary: {csv_path}")
    print(f"  Number of individual comparisons made: {len(rows)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
