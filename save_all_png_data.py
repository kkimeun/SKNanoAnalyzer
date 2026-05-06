import os
import ROOT

ROOT.gROOT.SetBatch(True)

input_dir = "/data9/Users/eunsu/SKNanoOutput/Tutorial_reco_tt/2024"
output_dir = os.path.join(input_dir, "png_data_merged")
os.makedirs(output_dir, exist_ok=True)

# 합칠 파일 목록
eras = ["C", "D", "E", "F", "G", "H", "I"]
samples = ["Muon0", "Muon1"]

input_files = [
    os.path.join(input_dir, f"{sample}_{era}.root")
    for sample in samples
    for era in eras
]

# 열 수 있는 파일만 사용
root_files = []
for path in input_files:
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        print(f"[WARNING] Cannot open file: {path}")
        continue
    print(f"[INFO] Opened: {path}")
    root_files.append(f)

if not root_files:
    raise RuntimeError("No valid ROOT files were opened.")

# path -> summed histogram 저장
merged_hists = {}

def collect_and_merge(tdir, current_path=""):
    for key in tdir.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()

        full_path = f"{current_path}/{name}" if current_path else name

        if obj.InheritsFrom("TDirectory"):
            collect_and_merge(obj, full_path)

        elif (
            obj.InheritsFrom("TH1")
            or obj.InheritsFrom("TH2")
            or obj.InheritsFrom("TProfile")
        ):
            # 파일이 닫혀도 살아있도록 clone
            if full_path not in merged_hists:
                merged_hists[full_path] = obj.Clone(full_path.replace("/", "__"))
                merged_hists[full_path].SetDirectory(0)
            else:
                merged_hists[full_path].Add(obj)

# 모든 파일에서 histogram 수집 및 합치기
for f in root_files:
    collect_and_merge(f)

# 플롯 저장
for full_path, hist in merged_hists.items():
    parts = full_path.split("/")
    hist_name = parts[-1]
    subdirs = parts[:-1]

    outdir = os.path.join(output_dir, *subdirs)
    os.makedirs(outdir, exist_ok=True)

    c = ROOT.TCanvas("c", "c", 800, 600)

    if hist.InheritsFrom("TH2"):
        hist.Draw("COLZ")
    else:
        hist.Draw("HIST")

    png_path = os.path.join(outdir, f"{hist_name}.png")
    c.SaveAs(png_path)
    c.Close()

# 파일 닫기
for f in root_files:
    f.Close()

print(f"\n[INFO] Saved merged PNGs under: {output_dir}")
