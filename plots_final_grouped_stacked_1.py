import os
import ROOT
import cmsstyle

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)

# =========================================================
# paths
# =========================================================
base_dir = "/data9/Users/eunsu/SKNanoOutput/Tutorial_reco_tt/2024"
output_dir = os.path.join(base_dir, "compare_data_grouped_mc_stack")
os.makedirs(output_dir, exist_ok=True)

# =========================================================
# data files
# =========================================================
eras = ["C", "D", "E", "F", "G", "H", "I"]
data_paths = []
for era in eras:
    data_paths.append(os.path.join(base_dir, f"Muon0_{era}.root"))
    data_paths.append(os.path.join(base_dir, f"Muon1_{era}.root"))

# =========================================================
# colors
# =========================================================
TT_COLOR = ROOT.TColor.GetColor("#ffbac4")       # flower

# =========================================================
# grouped MC config
# stack order: first added = bottom
# =========================================================
mc_groups = [
    {
        "name": "tt",
        "legend": "t#bar{t}",
        "color": TT_COLOR,
        "files": [
            "TTLJ_powheg.root",
        ],
    }
]

# =========================================================
# CMS style
# =========================================================
cmsstyle.setCMSStyle()
cmsstyle.SetLumi(110.0)
cmsstyle.SetEnergy(13)
cmsstyle.SetExtraText("Preliminary")

# =========================================================
# open files
# =========================================================
opened_data_files = []
for path in data_paths:
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        print(f"[WARNING] Cannot open data file: {path}")
        continue
    opened_data_files.append(f)
    print(f"[INFO] Opened data file: {path}")

if not opened_data_files:
    raise RuntimeError("No valid data ROOT files were opened.")

for group in mc_groups:
    group["handles"] = []
    for fname in group["files"]:
        full_path = os.path.join(base_dir, fname)
        f = ROOT.TFile.Open(full_path)
        if not f or f.IsZombie():
            print(f"[WARNING] Cannot open MC file: {full_path}")
            continue
        group["handles"].append(f)
        print(f"[INFO] Opened MC file: {full_path}")

if not any(len(group["handles"]) > 0 for group in mc_groups):
    raise RuntimeError("No valid MC ROOT files were opened.")

# =========================================================
# utilities
# =========================================================
def is_hist(obj):
    return (
        obj
        and obj.InheritsFrom("TH1")
        and not obj.InheritsFrom("TH2")
        and not obj.InheritsFrom("TProfile2D")
        and not obj.InheritsFrom("TProfile3D")
    )

def clone_hist(h, name):
    out = h.Clone(name)
    out.SetDirectory(0)
    out.SetStats(0)
    return out

def same_binning(h1, h2):
    return (
        h1.GetNbinsX() == h2.GetNbinsX()
        and abs(h1.GetXaxis().GetXmin() - h2.GetXaxis().GetXmin()) < 1e-9
        and abs(h1.GetXaxis().GetXmax() - h2.GetXaxis().GetXmax()) < 1e-9
    )

def get_obj_from_path(root_file, path):
    obj = root_file.Get(path)
    if not obj:
        return None
    return obj

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

def merge_hists_from_files(path, files, tag):
    merged = None
    for i, f in enumerate(files):
        h = get_obj_from_path(f, path)
        if not h or not is_hist(h):
            continue

        h = clone_hist(h, f"{path.replace('/', '__')}_{tag}_{i}")

        if merged is None:
            merged = clone_hist(h, f"{path.replace('/', '__')}_{tag}_merged")
        else:
            if not same_binning(merged, h):
                print(f"[WARNING] Binning mismatch for {path} in {f.GetName()}, skip")
                continue
            merged.Add(h)

    return merged

def merge_data_hist(path, files):
    return merge_hists_from_files(path, files, "data")

def merge_group_hist(path, group):
    if len(group["handles"]) == 0:
        return None
    return merge_hists_from_files(path, group["handles"], group["name"])

def build_group_hists(path, groups):
    out = []
    for group in groups:
        h = merge_group_hist(path, group)
        if not h:
            continue
        if h.Integral() <= 0:
            continue
        out.append((h, group))
    return out

def build_total_mc(group_hists, name):
    total = None
    for h, group in group_hists:
        if total is None:
            total = clone_hist(h, name)
        else:
            if not same_binning(total, h):
                return None
            total.Add(h)
    return total

def make_ratio_hist(h_data, h_mc, name):
    ratio = clone_hist(h_data, name)
    for i in range(1, ratio.GetNbinsX() + 1):
        d = h_data.GetBinContent(i)
        de = h_data.GetBinError(i)
        m = h_mc.GetBinContent(i)

        if m > 0:
            ratio.SetBinContent(i, d / m)
            ratio.SetBinError(i, de / m)
        else:
            ratio.SetBinContent(i, 0.0)
            ratio.SetBinError(i, 0.0)
    return ratio

def make_unc_band(h_mc, name):
    band = clone_hist(h_mc, name)
    for i in range(1, band.GetNbinsX() + 1):
        c = h_mc.GetBinContent(i)
        e = h_mc.GetBinError(i)
        if c > 0:
            band.SetBinContent(i, 1.0)
            band.SetBinError(i, e / c)
        else:
            band.SetBinContent(i, 0.0)
            band.SetBinError(i, 0.0)
    return band

def style_data(hist):
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.SetMarkerSize(1.0)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetLineColor(ROOT.kBlack)
    hist.SetLineWidth(2)

def style_unc(hist):
    hist.SetFillStyle(3345)
    hist.SetFillColor(12)
    hist.SetLineWidth(0)
    hist.SetMarkerSize(0)

def style_mc_component(hist, color):
    hist.SetFillColor(color)
    hist.SetLineColor(ROOT.kBlack)
    hist.SetLineWidth(1)

def draw_stack_with_ratio(h_data, group_hists, h_mc_total, out_path, label=""):
    if not group_hists or not h_mc_total:
        return

    # apply styles to each component
    for h, group in group_hists:
        style_mc_component(h, group["color"])

    stack = ROOT.THStack("stack", "")

    # IMPORTANT:
    # first add tt -> bottom
    # then add single top -> top
    for h, group in group_hists:
        stack.Add(h)

    x_min = h_data.GetXaxis().GetXmin()
    x_max = h_data.GetXaxis().GetXmax()
    x_title = h_data.GetXaxis().GetTitle()
    if not x_title:
        x_title = "plot"

    y_max = max(h_data.GetMaximum(), h_mc_total.GetMaximum()) * 1.8

    c = ROOT.TCanvas("c", "c", 800, 800)
    c.cd()

    pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0)
    pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.30)

    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.08)
    pad1.SetTicks()

    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.14)
    pad2.SetRightMargin(0.05)
    pad2.SetTicks()

    pad1.Draw()
    pad2.Draw()

    # ---------------- upper pad ----------------
    pad1.cd()

    stack.Draw("HIST")
    stack.SetMinimum(0.0)
    stack.SetMaximum(y_max)

    stack.GetYaxis().SetTitle("Events / bin")
    stack.GetYaxis().SetTitleSize(0.06)
    stack.GetYaxis().SetLabelSize(0.05)
    stack.GetYaxis().SetTitleOffset(1.1)

    stack.GetXaxis().SetTitle("")
    stack.GetXaxis().SetLabelSize(0)

    h_mc_unc = clone_hist(h_mc_total, "h_mc_unc")
    style_unc(h_mc_unc)
    h_mc_unc.Draw("E2 SAME")

    style_data(h_data)
    h_data.Draw("E SAME")

    # CMS-like labels
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
    latex.DrawLatex(0.72, 0.93, "110 fb^{-1} (13 TeV)")

    if label:
        latex.SetTextSize(0.035)
        latex.DrawLatex(0.18, 0.68, label)

    leg = ROOT.TLegend(0.45, 0.55, 0.92, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.040)
    leg.SetNColumns(2)

    leg.AddEntry(h_data, "Data", "lpe")

    # legend order: Data / tt / single top / unc
    for h, group in group_hists:
        leg.AddEntry(h, group["legend"], "f")

    leg.AddEntry(h_mc_unc, "Stat. Unc.", "f")
    leg.Draw()

    pad1.RedrawAxis()

    # ---------------- lower pad ----------------
    pad2.cd()

    data_ratio = make_ratio_hist(h_data, h_mc_total, "data_ratio")
    ratio_unc = make_unc_band(h_mc_total, "ratio_unc")

    style_data(data_ratio)
    style_unc(ratio_unc)

    ratio_unc.SetTitle("")
    ratio_unc.GetYaxis().SetTitle("Data / MC")
    ratio_unc.GetYaxis().SetTitleSize(0.10)
    ratio_unc.GetYaxis().SetLabelSize(0.09)
    ratio_unc.GetYaxis().SetTitleOffset(0.55)
    ratio_unc.GetYaxis().SetNdivisions(505)

    ratio_unc.GetXaxis().SetTitle(x_title)
    ratio_unc.GetXaxis().SetTitleSize(0.12)
    ratio_unc.GetXaxis().SetLabelSize(0.10)
    ratio_unc.GetXaxis().SetTitleOffset(1.1)

    ratio_unc.SetMinimum(0.5)
    ratio_unc.SetMaximum(1.5)

    ratio_unc.Draw("E2")
    data_ratio.Draw("E SAME")

    line = ROOT.TLine(x_min, 1.0, x_max, 1.0)
    line.SetLineStyle(2)
    line.Draw("SAME")

    pad2.RedrawAxis()

    c.SaveAs(out_path)
    c.Close()

# =========================================================
# main loop
# =========================================================
reference_file = None
for group in mc_groups:
    if len(group["handles"]) > 0:
        reference_file = group["handles"][0]
        break

if reference_file is None:
    raise RuntimeError("No valid reference MC file found.")

mc_hist_paths = collect_hist_paths(reference_file)
print(f"[INFO] Found {len(mc_hist_paths)} MC histograms in reference file")

n_done = 0
n_skip = 0

for path in mc_hist_paths:
    h_data = merge_data_hist(path, opened_data_files)
    if not h_data:
        print(f"[SKIP] No matching data hist: {path}")
        n_skip += 1
        continue

    group_hists = build_group_hists(path, mc_groups)
    if not group_hists:
        print(f"[SKIP] No MC histograms found: {path}")
        n_skip += 1
        continue

    h_mc_total = build_total_mc(group_hists, path.replace("/", "__") + "_mc_total")
    if not h_mc_total:
        print(f"[SKIP] Failed to build total MC: {path}")
        n_skip += 1
        continue

    if h_data.Integral() <= 0 or h_mc_total.Integral() <= 0:
        print(f"[SKIP] Empty hist: {path}")
        n_skip += 1
        continue

    parts = path.split("/")
    hist_name = parts[-1]
    subdirs = parts[:-1]

    out_subdir = os.path.join(output_dir, *subdirs)
    os.makedirs(out_subdir, exist_ok=True)
    out_path = os.path.join(out_subdir, f"{hist_name}.png")

    draw_stack_with_ratio(
        h_data=h_data,
        group_hists=group_hists,
        h_mc_total=h_mc_total,
        out_path=out_path,
        label=""
    )

    print(f"[DONE] {path}")
    n_done += 1

# =========================================================
# close files
# =========================================================
for f in opened_data_files:
    f.Close()

for group in mc_groups:
    for f in group["handles"]:
        f.Close()

print(f"\n[INFO] Finished. Done = {n_done}, Skipped = {n_skip}")
print(f"[INFO] Output directory: {output_dir}")
