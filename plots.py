import os
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)

# =========================
# paths
# =========================
base_dir = "/data9/Users/eunsu/SKNanoOutput/Tutorial_reco_tt/2024"
mc_file_path = os.path.join(base_dir, "TTLJ_powheg.root")
output_dir = os.path.join(base_dir, "compare_data_mc")
os.makedirs(output_dir, exist_ok=True)

eras = ["C", "D", "E", "F", "G", "H", "I"]
data_files = []
for era in eras:
    data_files.append(os.path.join(base_dir, f"Muon0_{era}.root"))
    data_files.append(os.path.join(base_dir, f"Muon1_{era}.root"))

# =========================
# open files
# =========================
mc_file = ROOT.TFile.Open(mc_file_path)
if not mc_file or mc_file.IsZombie():
    raise RuntimeError(f"Cannot open MC file: {mc_file_path}")

opened_data_files = []
for path in data_files:
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        print(f"[WARNING] Cannot open data file: {path}")
        continue
    opened_data_files.append(f)
    print(f"[INFO] Opened data file: {path}")

if not opened_data_files:
    raise RuntimeError("No valid data ROOT files were opened.")

# =========================
# utilities
# =========================
def is_hist(obj):
    return (
        obj.InheritsFrom("TH1")
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

def get_obj_from_path(root_file, path):
    obj = root_file.Get(path)
    if not obj:
        return None
    return obj

def merge_data_hist(path, files):
    merged = None
    for f in files:
        h = get_obj_from_path(f, path)
        if not h:
            continue
        if not is_hist(h):
            continue

        if merged is None:
            merged = h.Clone(path.replace("/", "__") + "_data")
            merged.SetDirectory(0)
        else:
            # binning mismatch 방지용 체크
            if (
                merged.GetNbinsX() != h.GetNbinsX()
                or abs(merged.GetXaxis().GetXmin() - h.GetXaxis().GetXmin()) > 1e-9
                or abs(merged.GetXaxis().GetXmax() - h.GetXaxis().GetXmax()) > 1e-9
            ):
                print(f"[WARNING] Binning mismatch for {path} in file {f.GetName()}, skip")
                continue
            merged.Add(h)

    return merged

def make_ratio_hist(h_data, h_mc, name):
    ratio = h_data.Clone(name)
    ratio.SetDirectory(0)
    ratio.Divide(h_mc)
    return ratio

def make_mc_unc_band(h_mc, name):
    band = h_mc.Clone(name)
    band.SetDirectory(0)
    for i in range(1, band.GetNbinsX() + 1):
        content = h_mc.GetBinContent(i)
        error = h_mc.GetBinError(i)
        if content > 0:
            band.SetBinContent(i, 1.0)
            band.SetBinError(i, error / content)
        else:
            band.SetBinContent(i, 0.0)
            band.SetBinError(i, 0.0)
    return band

def style_mc(hist):
    hist.SetFillColor(ROOT.kRed)
    hist.SetLineColor(ROOT.kRed)
    hist.SetLineWidth(2)

def style_data(hist):
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1.0)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetLineColor(ROOT.kBlack)
    hist.SetLineWidth(2)

def style_unc(hist):
    hist.SetFillColor(ROOT.kGray + 2)
    hist.SetFillStyle(3004)
    hist.SetLineColor(ROOT.kGray + 2)
    hist.SetMarkerSize(0)

def draw_compare_plot(h_data, h_mc, out_path, hist_title=""):
    c = ROOT.TCanvas("c", "c", 800, 800)
    c.cd()

    pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0)
    pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.30)

    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.05)
    pad1.SetTicks()

    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.14)
    pad2.SetRightMargin(0.05)
    pad2.SetTicks()

    pad1.Draw()
    pad2.Draw()

    # -------------------------
    # upper pad
    # -------------------------
    pad1.cd()

    style_mc(h_mc)
    style_data(h_data)

    mc_unc = h_mc.Clone("mc_unc")
    mc_unc.SetDirectory(0)
    style_unc(mc_unc)

    ymax = max(h_mc.GetMaximum(), h_data.GetMaximum()) * 1.4
    h_mc.SetMaximum(ymax)
    h_mc.SetMinimum(0)

    h_mc.SetTitle("")
    h_mc.GetYaxis().SetTitle("Events / bin")
    h_mc.GetYaxis().SetTitleSize(0.06)
    h_mc.GetYaxis().SetLabelSize(0.05)
    h_mc.GetYaxis().SetTitleOffset(1.1)
    h_mc.GetXaxis().SetLabelSize(0)

    h_mc.Draw("HIST")
    mc_unc.Draw("E2 SAME")
    h_data.Draw("E1 SAME")

    # CMS-like text
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(62)
    latex.SetTextSize(0.06)
    latex.DrawLatex(0.16, 0.83, "CMS")

    latex.SetTextFont(52)
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.16, 0.78, "Preliminary")

    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.16, 0.73, "Single Muon")

    latex.SetTextSize(0.045)
    latex.DrawLatex(0.75, 0.93, "110 fb^{-1} (13 TeV)")

#    if hist_title:
#        latex.SetTextSize(0.035)
#        latex.DrawLatex(0.16, 0.68, hist_title)

    leg = ROOT.TLegend(0.60, 0.72, 0.90, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(h_data, "Data", "ep")
    leg.AddEntry(h_mc, "TTLJ powheg", "f")
    leg.AddEntry(mc_unc, "Stat. Unc.", "f")
    leg.Draw()

    # -------------------------
    # lower pad
    # -------------------------
    pad2.cd()

    ratio = make_ratio_hist(h_data, h_mc, "ratio")
    ratio_unc = make_mc_unc_band(h_mc, "ratio_unc")
    style_unc(ratio_unc)

    ratio.SetTitle("")
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(0.9)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetLineColor(ROOT.kBlack)

    ratio.GetYaxis().SetTitle("Data / MC")
    ratio.GetYaxis().SetTitleSize(0.10)
    ratio.GetYaxis().SetLabelSize(0.09)
    ratio.GetYaxis().SetTitleOffset(0.55)
    ratio.GetYaxis().SetNdivisions(505)

    ratio.GetXaxis().SetTitle(h_data.GetXaxis().GetTitle())
    ratio.GetXaxis().SetTitleSize(0.12)
    ratio.GetXaxis().SetLabelSize(0.10)
    ratio.GetXaxis().SetTitleOffset(1.1)

    ratio.SetMinimum(0.5)
    ratio.SetMaximum(1.5)

    ratio.Draw("E1")
    ratio_unc.Draw("E2 SAME")
    ratio.Draw("E1 SAME")

    line = ROOT.TLine(
        ratio.GetXaxis().GetXmin(), 1.0,
        ratio.GetXaxis().GetXmax(), 1.0
    )
    line.SetLineStyle(2)
    line.Draw("SAME")

    c.SaveAs(out_path)
    c.Close()

# =========================
# main loop
# =========================
mc_hist_paths = collect_hist_paths(mc_file)

print(f"[INFO] Found {len(mc_hist_paths)} MC histograms")

n_done = 0
n_skip = 0

for path in mc_hist_paths:
    h_mc = get_obj_from_path(mc_file, path)
    if not h_mc or not is_hist(h_mc):
        n_skip += 1
        continue

    h_mc = h_mc.Clone(path.replace("/", "__") + "_mc")
    h_mc.SetDirectory(0)

    h_data = merge_data_hist(path, opened_data_files)
    if not h_data:
        print(f"[SKIP] No matching data hist: {path}")
        n_skip += 1
        continue

    # empty hist skip
    if h_mc.Integral() <= 0 or h_data.Integral() <= 0:
        print(f"[SKIP] Empty hist: {path}")
        n_skip += 1
        continue

    # output path
    parts = path.split("/")
    hist_name = parts[-1]
    subdirs = parts[:-1]
    out_subdir = os.path.join(output_dir, *subdirs)
    os.makedirs(out_subdir, exist_ok=True)
    out_path = os.path.join(out_subdir, f"{hist_name}.png")

    draw_compare_plot(h_data, h_mc, out_path, hist_title=hist_name)
    print(f"[DONE] {path}")
    n_done += 1

# close files
mc_file.Close()
for f in opened_data_files:
    f.Close()

print(f"\n[INFO] Finished. Done = {n_done}, Skipped = {n_skip}")
print(f"[INFO] Output directory: {output_dir}")
