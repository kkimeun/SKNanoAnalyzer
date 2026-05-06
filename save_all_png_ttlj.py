import os
import ROOT

ROOT.gROOT.SetBatch(True)

input_file = "/data9/Users/eunsu/SKNanoOutput/Tutorial_reco_tt/2024/WtoLNu_4J.root"
output_dir = "/data9/Users/eunsu/SKNanoOutput/Tutorial_reco_tt/2024/png_WtoLNu_4J"

os.makedirs(output_dir, exist_ok=True)

f = ROOT.TFile.Open(input_file)
if not f or f.IsZombie():
    raise RuntimeError(f"Cannot open file: {input_file}")

def save_dir(tdir, outdir):
    os.makedirs(outdir, exist_ok=True)

    for key in tdir.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()

        if obj.InheritsFrom("TDirectory"):
            save_dir(obj, os.path.join(outdir, name))
            continue

        if obj.InheritsFrom("TH1") or obj.InheritsFrom("TH2") or obj.InheritsFrom("TProfile"):
            c = ROOT.TCanvas("c", "c", 800, 600)

            if obj.InheritsFrom("TH2"):
                obj.Draw("COLZ")
            else:
                obj.Draw("HIST")

            png_path = os.path.join(outdir, f"{name}.png")
            c.SaveAs(png_path)
            c.Close()

save_dir(f, output_dir)
f.Close()

print(f"Saved PNGs under: {output_dir}")
