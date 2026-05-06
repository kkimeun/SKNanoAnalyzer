#!/usr/bin/env python3
import os
import re
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

INPUT_ROOT = "output.root"
OUTDIR = "eff_png"

def safe_name(name: str) -> str:
    # 파일명으로 쓰기 안전하게 변환
    name = name.replace("/", "_")
    name = name.replace(":", "_")
    name = name.replace(";", "_")
    name = name.replace(" ", "_")
    name = name.replace("#", "_")
    name = re.sub(r"_+", "_", name)
    return name.strip("_")

def main():
    if not os.path.exists(INPUT_ROOT):
        raise FileNotFoundError(f"Input ROOT file not found: {INPUT_ROOT}")

    os.makedirs(OUTDIR, exist_ok=True)

    f = ROOT.TFile.Open(INPUT_ROOT)
    if not f or f.IsZombie():
        raise RuntimeError(f"Failed to open ROOT file: {INPUT_ROOT}")

    keys = f.GetListOfKeys()
    if not keys:
        print("No keys found in ROOT file.")
        return

    n_saved = 0

    for key in keys:
        obj = key.ReadObj()
        name = obj.GetName()

        # eff 히스토그램만 저장
        if "##eff" not in name:
            continue

        if not obj.InheritsFrom("TH1"):
            continue

        canvas = ROOT.TCanvas("c", "c", 900, 700)

        if obj.InheritsFrom("TH2"):
            obj.SetTitle(name)
            obj.GetXaxis().SetTitle("|#eta|")
            obj.GetYaxis().SetTitle("p_{T} [GeV]")
            obj.Draw("COLZ TEXT")
        else:
            obj.SetTitle(name)
            obj.Draw("HIST")

        outname = os.path.join(OUTDIR, safe_name(name) + ".png")
        canvas.SaveAs(outname)
        canvas.Close()

        print(f"Saved: {outname}")
        n_saved += 1

    f.Close()
    print(f"\nDone. Saved {n_saved} efficiency plot(s).")

if __name__ == "__main__":
    main()
