import os
import ctypes
import math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

base_dir = "/data9/Users/eunsu/SKNanoOutput/Tutorial_reco_tt/2024"

mc_files = [
    "ST_tch_antitop_lep.root",
    "ST_tch_top_lep.root",
    "ST_tW_antitop_Lep.root",
    "ST_tW_antitop_Semilep.root",
    "ST_tW_top_Lep.root",
    "ST_tW_top_Semilep.root",
    "TTLJ_powheg.root",
    "TTLL_powheg.root",
    "WtoLNu_1J.root",
    "WtoLNu_2J.root",
    "WtoLNu_3J.root",
    "WtoLNu_4J.root"
]

def is_hist(obj):
    return (
        obj
        and obj.InheritsFrom("TH1")
        and not obj.InheritsFrom("TH2")
        and not obj.InheritsFrom("TProfile2D")
        and not obj.InheritsFrom("TProfile3D")
    )

def collect_hist_paths(tdir, current_path="", hist_paths=None):
    if hist_paths is None:
        hist_paths = []

    for key in tdir.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()
        full_path = f"{current_path}/{name}" if current_path else name

        if obj.InheritsFrom("TDirectory"):
            collect_hist_paths(obj, full_path, hist_paths)
        elif is_hist(obj):
            hist_paths.append(full_path)

    return hist_paths

def summarize_hist(h):
    err = ctypes.c_double(0.0)
    integral = h.IntegralAndError(1, h.GetNbinsX(), err)

    entries = h.GetEntries()
    eff_entries = h.GetEffectiveEntries()
    mean = h.GetMean()
    rms = h.GetRMS()
    err_val = float(err.value)
    rel_err = err_val / integral if integral > 0 else 0.0

    return {
        "integral": integral,
        "error": err.value,
        "rel_err": rel_err,
        "entries": entries,
        "eff_entries": eff_entries,
        "mean": mean,
        "rms": rms,
    }

def print_summary_for_hist(hist_path):
    print(f"\n===== Histogram: {hist_path} =====")
    print(
        f"{'sample':35s} "
        f"{'Integral':>12s} {'Err':>12s} {'RelErr':>10s} "
        f"{'Entries':>12s} {'EffEntries':>12s} {'Mean':>10s} {'RMS':>10s}"
    )

    for fname in mc_files:
        fpath = os.path.join(base_dir, fname)
        f = ROOT.TFile.Open(fpath)
        if not f or f.IsZombie():
            print(f"{fname:35s}  [FAILED TO OPEN]")
            continue

        h = f.Get(hist_path)
        if not h or not is_hist(h):
            print(f"{fname:35s}  [HIST NOT FOUND]")
            f.Close()
            continue

        s = summarize_hist(h)
        print(
            f"{fname:35s} "
            f"{s['integral']:12.4f} {s['error']:12.4f} {s['rel_err']:10.4f} "
            f"{s['entries']:12.1f} {s['eff_entries']:12.1f} "
            f"{s['mean']:10.3f} {s['rms']:10.3f}"
        )
        f.Close()

def print_all_hist_summaries():
    ref_file = ROOT.TFile.Open(os.path.join(base_dir, mc_files[0]))
    if not ref_file or ref_file.IsZombie():
        raise RuntimeError("Cannot open reference file")

    hist_paths = collect_hist_paths(ref_file)
    ref_file.Close()

    print(f"[INFO] Found {len(hist_paths)} histograms")
    for hist_path in hist_paths:
        print_summary_for_hist(hist_path)

if __name__ == "__main__":
    # 1) 특정 히스토그램만 보고 싶으면 여기 지정
    # 예시:
    # print_summary_for_hist("Central__DiMuCut__had_top_mass__Central__mc")

    # 2) 전부 보고 싶으면 이거 사용
    print_all_hist_summaries()
