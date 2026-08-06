#!/usr/bin/env python3

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


# ======================================================================
# Configuration
# ======================================================================

SKNANO_HOME = Path(
    os.environ.get(
        "SKNANO_HOME",
        "/data6/Users/eunsu/SKNanoAnalyzer",
    )
)

JSON_DIR = (
    SKNANO_HOME
    / "data"
    / "Run3_v15_Run2_v15"
    / "2024"
    / "Sample"
    / "ForSNU"
)

OUTPUT_DIR = Path.cwd() / "nanoaod_filter_lepton_pt_comparison"

MASS_POINTS: List[Tuple[int, int]] = [
    (70, 15),
    (100, 60),
    (130, 90),
    (160, 155),
]

CHANNELS = {
    "Muon": {
        "mask": (
            "Muon_pt > 15.f && "
            "abs(Muon_eta) < 2.4f && "
            "Muon_tightId && "
            "Muon_pfRelIso04_all < 0.15f"
        ),
        "selected_pt_name": "SelectedMuonPt",
        "leading_pt_name": "LeadingSelectedMuonPt",
        "pt_branch": "Muon_pt",
        "nbins": 80,
        "xmin": 0.0,
        "xmax": 400.0,
        "x_title": "Leading selected reco muon p_{T} [GeV]",
    },

    "Electron": {
        "mask": (
            "Electron_pt > 15.f && "
            "abs(Electron_eta) < 2.5f && "
            "Electron_cutBased >= 4"
        ),
        "selected_pt_name": "SelectedElectronPt",
        "leading_pt_name": "LeadingSelectedElectronPt",
        "pt_branch": "Electron_pt",
        "nbins": 80,
        "xmin": 0.0,
        "xmax": 400.0,
        "x_title": "Leading selected reco electron p_{T} [GeV]",
    },

}

# ======================================================================
# Required NanoAOD branches for file validation
# ======================================================================

REQUIRED_BRANCHES = {
    "Muon": [
        "nMuon",
        "Muon_pt",
        "Muon_eta",
        "Muon_tightId",
        "Muon_pfRelIso04_all",
    ],
    "Electron": [
        "nElectron",
        "Electron_pt",
        "Electron_eta",
        "Electron_cutBased",
    ],
}

# True:
#   각 histogram을 면적 1로 normalize하여 shape만 비교
#
# False:
#   JSON에 포함된 실제 processed event 수 그대로 비교
#
# 필터 bias를 확인하는 첫 그림은 True가 적절함.
NORMALIZE = True

# overflow까지 normalization에 포함할지 여부
INCLUDE_OVERFLOW_IN_NORMALIZATION = True


# ======================================================================
# JSON and TChain helpers
# ======================================================================

def json_path(mhc: int, ma: int, filtered: bool) -> Path:
    suffix = "_SingleLepFilter" if filtered else ""

    return JSON_DIR / (
        f"TTToHcToWAToBB-MHc{mhc}_MA{ma}{suffix}.json"
    )


def load_root_paths(path: Path) -> List[str]:
    if not path.is_file():
        raise FileNotFoundError(f"JSON 파일을 찾을 수 없음: {path}")

    with path.open("r", encoding="utf-8") as json_file:
        sample_info = json.load(json_file)

    root_paths = sample_info.get("path")

    if not isinstance(root_paths, list) or len(root_paths) == 0:
        raise RuntimeError(
            f"'path' 목록이 없거나 비어 있음: {path}"
        )

    valid_paths: List[str] = []
    missing_paths: List[str] = []

    for root_path in root_paths:
        # root:// 주소는 os.path.isfile로 검사할 수 없으므로 그대로 추가
        if root_path.startswith("root://"):
            valid_paths.append(root_path)
            continue

        if os.path.isfile(root_path):
            valid_paths.append(root_path)
        else:
            missing_paths.append(root_path)

    if missing_paths:
        print(
            f"[WARNING] {path.name}: 존재하지 않는 ROOT 파일 "
            f"{len(missing_paths)}개"
        )

        for missing in missing_paths[:5]:
            print(f"          {missing}")

        if len(missing_paths) > 5:
            print(
                f"          ... and {len(missing_paths) - 5} more"
            )

    if len(valid_paths) == 0:
        raise RuntimeError(
            f"사용 가능한 ROOT 파일이 하나도 없음: {path}"
        )

    return valid_paths

def is_usable_root_file(root_path: str, channel: str) -> bool:
    """
    ROOT 파일을 TChain에 넣기 전에 검사한다.

    검사 항목:
      1. 파일을 정상적으로 열 수 있는지
      2. recovered/zombie 파일이 아닌지
      3. Events tree가 있는지
      4. 해당 채널에 필요한 branch가 모두 있는지
      5. 첫 event를 읽을 수 있는지

    문제가 있으면 False를 반환하여 해당 파일을 건너뛴다.
    """
    root_file = None

    try:
        root_file = ROOT.TFile.Open(root_path, "READ")

        if not root_file or root_file.IsZombie():
            print(f"[SKIP] Cannot open ROOT file: {root_path}")
            return False

        if root_file.TestBit(ROOT.TFile.kRecovered):
            print(f"[SKIP] Recovered/corrupted ROOT file: {root_path}")
            return False

        tree = root_file.Get("Events")

        if not tree or not tree.InheritsFrom("TTree"):
            print(f"[SKIP] Missing Events tree: {root_path}")
            return False

        missing_branches = [
            branch_name
            for branch_name in REQUIRED_BRANCHES[channel]
            if not tree.GetBranch(branch_name)
        ]

        if missing_branches:
            print(
                f"[SKIP] Missing branches for {channel}: {root_path}\n"
                f"       {', '.join(missing_branches)}"
            )
            return False

        # 빈 tree는 꼭 오류는 아니지만 분석에는 기여하지 않으므로 제외
        entries = tree.GetEntries()

        if entries <= 0:
            print(f"[SKIP] Empty Events tree: {root_path}")
            return False

        # 최소한 첫 event의 basket을 읽을 수 있는지 확인
        if tree.GetEntry(0) < 0:
            print(f"[SKIP] Cannot read first event: {root_path}")
            return False

        return True

    except Exception as error:
        print(
            f"[SKIP] Exception while checking ROOT file: {root_path}\n"
            f"       {error}"
        )
        return False

    finally:
        if root_file:
            root_file.Close()


def make_chain(
    mhc: int,
    ma: int,
    filtered: bool,
    channel: str,
) -> ROOT.TChain:
    sample_json = json_path(mhc, ma, filtered)
    root_paths = load_root_paths(sample_json)

    chain = ROOT.TChain("Events")

    added_files = 0
    skipped_files = 0

    for root_path in root_paths:
        # ==============================================================
        # [ADD] 깨진 파일 또는 필요한 branch가 없는 파일 건너뛰기
        # ==============================================================
        if not is_usable_root_file(root_path, channel):
            skipped_files += 1
            continue

        result = chain.Add(root_path)

        if result > 0:
            added_files += 1
        else:
            skipped_files += 1
            print(f"[SKIP] TChain.Add failed: {root_path}")

    if added_files == 0:
        raise RuntimeError(
            f"No usable ROOT files:\n"
            f"  sample  = {sample_json}\n"
            f"  channel = {channel}"
        )

    # TChain이 첫 tree를 실제로 로드하도록 강제
    if chain.LoadTree(0) < 0:
        raise RuntimeError(
            f"Failed to load first valid tree:\n"
            f"  sample  = {sample_json}\n"
            f"  channel = {channel}"
        )

    # chain 수준에서도 branch를 마지막으로 확인
    missing_chain_branches = [
        branch_name
        for branch_name in REQUIRED_BRANCHES[channel]
        if not chain.GetBranch(branch_name)
    ]

    if missing_chain_branches:
        raise RuntimeError(
            f"TChain is missing branches for {channel}:\n"
            f"  {', '.join(missing_chain_branches)}"
        )

    print(
        f"[CHAIN] MHc={mhc}, MA={ma}, "
        f"filtered={filtered}, channel={channel}: "
        f"{added_files} usable files, "
        f"{skipped_files} skipped files, "
        f"{chain.GetEntries()} events"
    )

    return chain

# ======================================================================
# Histogram helpers
# ======================================================================

def make_histogram(
    chain: ROOT.TChain,
    channel: str,
    hist_name: str,
) -> ROOT.TH1D:
    config = CHANNELS[channel]

    dataframe = ROOT.RDataFrame(chain)

    selected_pt_name = str(config["selected_pt_name"])
    leading_pt_name = str(config["leading_pt_name"])

    selected_df = (
        dataframe
        .Define(
            selected_pt_name,
            f"{config['pt_branch']}[{config['mask']}]"
        )
        .Filter(
            f"{selected_pt_name}.size() > 0"
        )
        .Define(
            leading_pt_name,
            f"ROOT::VecOps::Max({selected_pt_name})"
        )
    )

    hist_result = selected_df.Histo1D(
        (
            hist_name,
            "",
            int(config["nbins"]),
            float(config["xmin"]),
            float(config["xmax"]),
        ),
        leading_pt_name,
    )

    histogram = hist_result.GetValue().Clone(hist_name + "_clone")
    histogram.SetDirectory(0)
    histogram.Sumw2()

    selected_entries = int(selected_df.Count().GetValue())

    print(
        f"[HIST] {hist_name}: "
        f"selected events={selected_entries}, "
        f"integral={histogram.Integral():.0f}"
    )

    return histogram


def histogram_integral(histogram: ROOT.TH1) -> float:
    if INCLUDE_OVERFLOW_IN_NORMALIZATION:
        return histogram.Integral(
            0,
            histogram.GetNbinsX() + 1,
        )

    return histogram.Integral(
        1,
        histogram.GetNbinsX(),
    )


def normalize_histogram(histogram: ROOT.TH1) -> None:
    integral = histogram_integral(histogram)

    if integral == 0.0:
        print(
            f"[WARNING] integral이 0이라 normalize하지 않음: "
            f"{histogram.GetName()}"
        )
        return

    histogram.Scale(1.0 / integral)


def make_ratio(
    numerator: ROOT.TH1,
    denominator: ROOT.TH1,
    name: str,
) -> ROOT.TH1:
    ratio = numerator.Clone(name)
    ratio.SetDirectory(0)
    ratio.Divide(denominator)

    return ratio


# ======================================================================
# Plot one mass point
# ======================================================================

def draw_mass_point(
    channel: str,
    mhc: int,
    ma: int,
) -> None:
    
    chain_no_filter = make_chain(
        mhc=mhc,
        ma=ma,
        filtered=False,
        channel=channel,  # [ADD]
    )

    chain_filtered = make_chain(
        mhc=mhc,
        ma=ma,
        filtered=True,
        channel=channel,  # [ADD]
    )

    h_no_filter = make_histogram(
        chain_no_filter,
        channel,
        f"h_{channel}_MHc{mhc}_MA{ma}_NoFilter",
    )

    h_filtered = make_histogram(
        chain_filtered,
        channel,
        f"h_{channel}_MHc{mhc}_MA{ma}_Filtered",
    )

    no_filter_entries = histogram_integral(h_no_filter)
    filtered_entries = histogram_integral(h_filtered)

    if NORMALIZE:
        normalize_histogram(h_no_filter)
        normalize_histogram(h_filtered)

    # Style
    h_no_filter.SetLineColor(ROOT.kBlack)
    h_no_filter.SetMarkerColor(ROOT.kBlack)
    h_no_filter.SetLineWidth(3)

    h_filtered.SetLineColor(ROOT.kRed + 1)
    h_filtered.SetMarkerColor(ROOT.kRed + 1)
    h_filtered.SetLineWidth(3)

    canvas = ROOT.TCanvas(
        f"canvas_{channel}_MHc{mhc}_MA{ma}",
        "",
        850,
        850,
    )

    upper_pad = ROOT.TPad(
        "upper_pad",
        "",
        0.0,
        0.30,
        1.0,
        1.0,
    )

    lower_pad = ROOT.TPad(
        "lower_pad",
        "",
        0.0,
        0.00,
        1.0,
        0.30,
    )

    upper_pad.SetLeftMargin(0.13)
    upper_pad.SetRightMargin(0.04)
    upper_pad.SetTopMargin(0.08)
    upper_pad.SetBottomMargin(0.02)

    lower_pad.SetLeftMargin(0.13)
    lower_pad.SetRightMargin(0.04)
    lower_pad.SetTopMargin(0.04)
    lower_pad.SetBottomMargin(0.34)
    lower_pad.SetGridy()

    upper_pad.Draw()
    lower_pad.Draw()

    # --------------------------------------------------------------
    # Upper panel
    # --------------------------------------------------------------
    upper_pad.cd()

    maximum = max(
        h_no_filter.GetMaximum(),
        h_filtered.GetMaximum(),
    )

    h_no_filter.SetMaximum(
        1.35 * maximum if maximum > 0.0 else 1.0
    )
    h_no_filter.SetMinimum(0.0)

    y_title = (
        "Normalized events"
        if NORMALIZE
        else "Events"
    )

    h_no_filter.SetTitle(
        f"M_{{H^{{#pm}}}}={mhc} GeV, "
        f"M_{{A}}={ma} GeV;"
        f"{CHANNELS[channel]['x_title']};"
        f"{y_title}"
    )

    h_no_filter.GetXaxis().SetLabelSize(0.0)
    h_no_filter.GetXaxis().SetTitleSize(0.0)

    h_no_filter.GetYaxis().SetTitleSize(0.055)
    h_no_filter.GetYaxis().SetTitleOffset(1.10)
    h_no_filter.GetYaxis().SetLabelSize(0.045)

    h_no_filter.Draw("HIST")
    h_filtered.Draw("HIST SAME")

    legend = ROOT.TLegend(
        0.54,
        0.69,
        0.92,
        0.88,
    )
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.038)

    legend.AddEntry(
        h_no_filter,
        f"No filter ({int(no_filter_entries)} reco entries)",
        "l",
    )
    legend.AddEntry(
        h_filtered,
        f"SingleLepFilter ({int(filtered_entries)} reco entries)",
        "l",
    )
    legend.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.043)
    text.DrawLatex(
        0.15,
        0.93,
        f"{channel}: NanoAOD before analysis selection",
    )

    # --------------------------------------------------------------
    # Ratio panel
    # --------------------------------------------------------------
    lower_pad.cd()

    ratio = make_ratio(
        h_filtered,
        h_no_filter,
        f"ratio_{channel}_MHc{mhc}_MA{ma}",
    )

    ratio.SetTitle("")
    ratio.SetLineColor(ROOT.kRed + 1)
    ratio.SetMarkerColor(ROOT.kRed + 1)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(0.65)

    ratio.GetYaxis().SetTitle(
        "Filter / No filter"
    )
    ratio.GetYaxis().SetRangeUser(0.5, 1.5)
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().SetTitleSize(0.10)
    ratio.GetYaxis().SetTitleOffset(0.58)
    ratio.GetYaxis().SetLabelSize(0.085)

    ratio.GetXaxis().SetTitle(
        str(CHANNELS[channel]["x_title"])
    )
    ratio.GetXaxis().SetTitleSize(0.12)
    ratio.GetXaxis().SetTitleOffset(1.05)
    ratio.GetXaxis().SetLabelSize(0.10)

    ratio.Draw("E1")

    unity = ROOT.TLine(
        ratio.GetXaxis().GetXmin(),
        1.0,
        ratio.GetXaxis().GetXmax(),
        1.0,
    )
    unity.SetLineStyle(2)
    unity.SetLineWidth(2)
    unity.Draw("SAME")

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    mode = "shape" if NORMALIZE else "entries"

    output_base = OUTPUT_DIR / (
        f"{channel}_NanoAOD_leadingPt_"
        f"MHc{mhc}_MA{ma}_{mode}"
    )

    canvas.SaveAs(str(output_base) + ".png")
    canvas.SaveAs(str(output_base) + ".pdf")

    print(f"[SAVED] {output_base}.png")
    print(f"[SAVED] {output_base}.pdf")


# ======================================================================
# Four-mass-point summary
# ======================================================================

def draw_summary(channel: str) -> None:
    canvas = ROOT.TCanvas(
        f"summary_{channel}",
        "",
        1500,
        1150,
    )

    canvas.Divide(2, 2)

    # ROOT garbage collection 방지
    keep_alive = []

    for pad_index, (mhc, ma) in enumerate(
        MASS_POINTS,
        start=1,
    ):
        canvas.cd(pad_index)

        ROOT.gPad.SetLeftMargin(0.14)
        ROOT.gPad.SetRightMargin(0.04)
        ROOT.gPad.SetTopMargin(0.10)
        ROOT.gPad.SetBottomMargin(0.13)

        chain_no_filter = make_chain(
            mhc,
            ma,
            False,
            channel,
        )

        chain_filtered = make_chain(
            mhc,
            ma,
            True,
            channel,
        )

        h_no_filter = make_histogram(
            chain_no_filter,
            channel,
            f"h_summary_{channel}_{mhc}_{ma}_no",
        )

        h_filtered = make_histogram(
            chain_filtered,
            channel,
            f"h_summary_{channel}_{mhc}_{ma}_filter",
        )

        if NORMALIZE:
            normalize_histogram(h_no_filter)
            normalize_histogram(h_filtered)

        h_no_filter.SetLineColor(ROOT.kBlack)
        h_no_filter.SetLineWidth(3)

        h_filtered.SetLineColor(ROOT.kRed + 1)
        h_filtered.SetLineWidth(3)

        maximum = max(
            h_no_filter.GetMaximum(),
            h_filtered.GetMaximum(),
        )

        h_no_filter.SetMaximum(
            1.32 * maximum if maximum > 0 else 1.0
        )
        h_no_filter.SetMinimum(0.0)

        y_title = (
            "Normalized events"
            if NORMALIZE
            else "Events"
        )

        h_no_filter.SetTitle(
            f"M_{{H^{{#pm}}}}={mhc} GeV, "
            f"M_{{A}}={ma} GeV;"
            f"{CHANNELS[channel]['x_title']};"
            f"{y_title}"
        )

        h_no_filter.Draw("HIST")
        h_filtered.Draw("HIST SAME")

        legend = ROOT.TLegend(
            0.52,
            0.73,
            0.91,
            0.88,
        )
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.035)
        legend.AddEntry(
            h_no_filter,
            "No generator filter",
            "l",
        )
        legend.AddEntry(
            h_filtered,
            "SingleLepFilter",
            "l",
        )
        legend.Draw()

        keep_alive.extend(
            [
                chain_no_filter,
                chain_filtered,
                h_no_filter,
                h_filtered,
                legend,
            ]
        )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    mode = "shape" if NORMALIZE else "entries"

    output_base = OUTPUT_DIR / (
        f"{channel}_NanoAOD_leadingPt_"
        f"all_masspoints_{mode}"
    )

    canvas.SaveAs(str(output_base) + ".png")
    canvas.SaveAs(str(output_base) + ".pdf")

    print(f"[SAVED] {output_base}.png")
    print(f"[SAVED] {output_base}.pdf")


# ======================================================================
# Main
# ======================================================================

def main() -> int:
    print(f"[INFO] SKNANO_HOME = {SKNANO_HOME}")
    print(f"[INFO] JSON_DIR    = {JSON_DIR}")
    print(f"[INFO] OUTPUT_DIR  = {OUTPUT_DIR}")
    print(f"[INFO] NORMALIZE   = {NORMALIZE}")

    if not JSON_DIR.is_dir():
        print(
            f"[ERROR] JSON directory가 존재하지 않음: {JSON_DIR}",
            file=sys.stderr,
        )
        return 1

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    for channel in CHANNELS:
        for mhc, ma in MASS_POINTS:
            try:
                draw_mass_point(
                    channel=channel,
                    mhc=mhc,
                    ma=ma,
                )
            except Exception as error:
                print(
                    f"[ERROR] {channel}, "
                    f"MHc={mhc}, MA={ma}: {error}",
                    file=sys.stderr,
                )

        try:
            draw_summary(channel)
        except Exception as error:
            print(
                f"[ERROR] {channel} summary: {error}",
                file=sys.stderr,
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
