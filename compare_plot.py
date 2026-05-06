import os
import ROOT
import cmsstyle as CMS

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

# =========================================================
# user config
# =========================================================
base_dir = "/data9/Users/eunsu/SKNanoOutput/Tutorial_reco_tt/2024"
mc_file_path = os.path.join(base_dir, "TTLJ_powheg.root")
output_dir = os.path.join(base_dir, "compare_data_mc_cmsstyle")

eras = ["C", "D", "E", "F", "G", "H", "I"]
data_file_paths = []
for era in eras:
    data_file_paths.append(os.path.join(base_dir, f"Muon0_{era}.root"))
    data_file_paths.append(os.path.join(base_dir, f"Muon1_{era}.root"))

lumi_text = "110 fb^{-1} (13 TeV)"
channel_text = "Single Muon"
mc_legend_name = "TTLJ powheg"

os.makedirs(output_dir, exist_ok=True)

# =========================================================
# global style
# =========================================================
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetEndErrorSize(4)
ROOT.gStyle.SetHatchesSpacing(1.0)
ROOT.gStyle.SetHatchesLineWidth(1)

# cmsstyle text setup
CMS.SetExtraText("Preliminary")
CMS.SetLumi(lumi_text)

# =========================================================
# utilities
# =========================================================
def open_root_file(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        return None
    return f

def is_1d_hist(obj):
    if not obj:
        return False
    if not obj.InheritsFrom("TH1"):
        return False
    if obj.InheritsFrom("TH2"):
        return False
    if obj.InheritsFrom("TProfile"):
        return False
    return True

def collect_hist_paths(tdir, current_path="", hist_paths=None):
    if hist_paths is None:
        hist_paths = []

    for key in tdir.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()
        full_path = f"{current_path}/{name}" if current_path else name

        if obj.InheritsFrom("TDirectory"):
            collect_hist_paths(obj, full_path, hist_paths)
        elif is_1d_hist(obj):
            hist_paths.append(full_path)

    return hist_paths

def get_hist(file_handle, path, clone_name=None):
    obj = file_handle.Get(path)
    if not obj or not is_1d_hist(obj):
        return None

    if clone_name is None:
        clone_name = path.replace("/", "__")

    h = obj.Clone(clone_name)
    h.SetDirectory(0)
    h.SetStats(0)
    return h

def same_binning(h1, h2):
    if h1.GetNbinsX() != h2.GetNbinsX():
        return False

    ax1 = h1.GetXaxis()
    ax2 = h2.GetXaxis()

    if abs(ax1.GetXmin() - ax2.GetXmin()) > 1e-9:
        return False
    if abs(ax1.GetXmax() - ax2.GetXmax()) > 1e-9:
        return False

    for i in range(1, h1.GetNbinsX() + 2):
        if abs(ax1.GetBinLowEdge(i) - ax2.GetBinLowEdge(i)) > 1e-9:
            return False

    return True

def merge_data_hist(hist_path, data_files):
    merged = None

    for idx, f in enumerate(data_files):
        h = get_hist(f, hist_path, clone_name=f"data_{idx}__{hist_path.replace('/', '__')}")
        if h is None:
            continue

        if merged is None:
            merged = h.Clone(f"merged__{hist_path.replace('/', '__')}")
            merged.SetDirectory(0)
            merged.SetStats(0)
        else:
            if not same_binning(merged, h):
                print(f"[WARNING] Binning mismatch, skip: {f.GetName()} :: {hist_path}")
                continue
            merged.Add(h)

    return merged

def prettify_hist_name(name):
    return name.replace("__", " ")

def style_mc(hist):
    hist.SetFillColor(ROOT.kRed)
    hist.SetLineColor(ROOT.kRed)
    hist.SetLineWidth(2)
    hist.SetMarkerSize(0)
    hist.SetStats(0)

def style_data(hist):
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1.0)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetLineColor(ROOT.kBlack)
    hist.SetLineWidth(2)
    hist.SetStats(0)

def style_unc_band(hist):
    hist.SetFillColor(ROOT.kGray + 2)
    hist.SetFillStyle(3354)
    hist.SetLineColor(ROOT.kGray + 2)
    hist.SetLineWidth(1)
    hist.SetMarkerSize(0)
    hist.SetStats(0)

def fold_overflow(hist):
    nb = hist.GetNbinsX()

    # underflow -> first bin
    c1 = hist.GetBinContent(1)
    e1 = hist.GetBinError(1)
    cu = hist.GetBinContent(0)
    eu = hist.GetBinError(0)
    hist.SetBinContent(1, c1 + cu)
    hist.SetBinError(1, (e1**2 + eu**2) ** 0.5)

    # overflow -> last bin
    cn = hist.GetBinContent(nb)
    en = hist.GetBinError(nb)
    co = hist.GetBinContent(nb + 1)
    eo = hist.GetBinError(nb + 1)
    hist.SetBinContent(nb, cn + co)
    hist.SetBinError(nb, (en**2 + eo**2) ** 0.5)

def make_ratio_hist(h_data, h_mc):
    ratio = h_data.Clone("ratio")
    ratio.SetDirectory(0)
    ratio.SetStats(0)

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

def make_ratio_unc_band(h_mc):
    band = h_mc.Clone("ratio_unc_band")
    band.SetDirectory(0)
    band.SetStats(0)

    for i in range(1, band.GetNbinsX() + 1):
        m = h_mc.GetBinContent(i)
        me = h_mc.GetBinError(i)

        if m > 0:
            band.SetBinContent(i, 1.0)
            band.SetBinError(i, me / m)
        else:
            band.SetBinContent(i, 0.0)
            band.SetBinError(i, 0.0)

    return band

def set_axis_style_upper(hist):
    hist.SetTitle("")
    hist.GetYaxis().SetTitle("Events / bin")
    hist.GetYaxis().SetTitleFont(42)
    hist.GetYaxis().SetLabelFont(42)
    hist.GetYaxis().SetTitleSize(0.065)
    hist.GetYaxis().SetLabelSize(0.050)
    hist.GetYaxis().SetTitleOffset(0.95)

    hist.GetXaxis().SetTitle("")
    hist.GetXaxis().SetLabelSize(0)
    hist.GetXaxis().SetTitleSize(0)

def set_axis_style_ratio(hist, x_title):
    hist.SetTitle("")
    hist.GetYaxis().SetTitle("Data / MC")
    hist.GetYaxis().SetTitleFont(42)
    hist.GetYaxis().SetLabelFont(42)
    hist.GetYaxis().SetTitleSize(0.11)
    hist.GetYaxis().SetLabelSize(0.09)
    hist.GetYaxis().SetTitleOffset(0.50)
    hist.GetYaxis().SetNdivisions(505)

    hist.GetXaxis().SetTitle(x_title)
    hist.GetXaxis().SetTitleFont(42)
    hist.GetXaxis().SetLabelFont(42)
    hist.GetXaxis().SetTitleSize(0.12)
    hist.GetXaxis().SetLabelSize(0.10)
    hist.GetXaxis().SetTitleOffset(1.10)

def draw_header_text(hist_title, channel_text, lumi_text):
    latex = ROOT.TLatex()
    latex.SetNDC()

    latex.SetTextAlign(11)

    # CMS
    latex.SetTextFont(62)
    latex.SetTextSize(0.070)
    latex.DrawLatex(0.16, 0.88, "CMS")

    # Preliminary
    latex.SetTextFont(52)
    latex.SetTextSize(0.050)
    latex.DrawLatex(0.16, 0.81, "Preliminary")

    # Channel
    latex.SetTextFont(42)
    latex.SetTextSize(0.052)
    latex.DrawLatex(0.16, 0.74, channel_text)

    # Hist title
    if hist_title:
        latex.SetTextFont(42)
        latex.SetTextSize(0.040)
        latex.DrawLatex(0.16, 0.66, hist_title)

    # Lumi
    latex.SetTextAlign(31)
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.95, 0.93, lumi_text)

def draw_compare_plot(h_data, h_mc, out_path, hist_title):
    # safety clones
    h_data = h_data.Clone("h_data_draw")
    h_mc = h_mc.Clone("h_mc_draw")
    h_data.SetDirectory(0)
    h_mc.SetDirectory(0)
    h_data.SetStats(0)
    h_mc.SetStats(0)

    fold_overflow(h_data)
    fold_overflow(h_mc)

    c = ROOT.TCanvas("c", "c", 800, 800)
    c.cd()

    pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.30, 1.0, 1.0)
    pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.00, 1.0, 0.30)

    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.08)
    pad1.SetTicks(1, 1)

    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.14)
    pad2.SetRightMargin(0.05)
    pad2.SetTicks(1, 1)

    pad1.Draw()
    pad2.Draw()

    # ---------------- upper ----------------
    pad1.cd()

    style_mc(h_mc)
    style_data(h_data)

    mc_unc = h_mc.Clone("mc_unc")
    mc_unc.SetDirectory(0)
    style_unc_band(mc_unc)

    set_axis_style_upper(h_mc)

    ymax = max(h_mc.GetMaximum(), h_data.GetMaximum()) * 1.45
    h_mc.SetMaximum(ymax)
    h_mc.SetMinimum(0.0)

    h_mc.Draw("HIST")
    mc_unc.Draw("E2 SAME")
    h_data.Draw("E1 SAME")

    draw_header_text(hist_title, channel_text, lumi_text)

    leg = ROOT.TLegend(0.58, 0.72, 0.90, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.040)
    leg.AddEntry(h_data, "Data", "ep")
    leg.AddEntry(h_mc, mc_legend_name, "f")
    leg.AddEntry(mc_unc, "Stat. Unc.", "f")
    leg.Draw()

    pad1.RedrawAxis()

    # ---------------- lower ----------------
    pad2.cd()

    ratio = make_ratio_hist(h_data, h_mc)
    ratio_unc = make_ratio_unc_band(h_mc)

    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(0.95)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(1)

    style_unc_band(ratio_unc)
    set_axis_style_ratio(ratio, h_data.GetXaxis().GetTitle())

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
    line.SetLineWidth(2)
    line.SetLineColor(ROOT.kBlack)
    line.Draw("SAME")

    pad2.RedrawAxis()

    c.SaveAs(out_path)
    c.Close()

# =========================================================
# main
# =========================================================
print("[INFO] Opening MC file...")
mc_file = open_root_file(mc_file_path)
if mc_file is None:
    raise RuntimeError(f"Cannot open MC ROOT file: {mc_file_path}")

print("[INFO] Opening data files...")
data_files = []
for path in data_file_paths:
    f = open_root_file(path)
    if f is None:
        print(f"[WARNING] Cannot open data file: {path}")
        continue
    print(f"[INFO] Opened: {path}")
    data_files.append(f)

if len(data_files) == 0:
    raise RuntimeError("No valid data ROOT files opened.")

print("[INFO] Collecting MC histogram paths...")
mc_hist_paths = collect_hist_paths(mc_file)
print(f"[INFO] Found {len(mc_hist_paths)} TH1 histograms in MC file.")

n_done = 0
n_skip = 0

for hist_path in mc_hist_paths:
    h_mc = get_hist(mc_file, hist_path, clone_name=f"mc__{hist_path.replace('/', '__')}")
    if h_mc is None:
        n_skip += 1
        continue

    h_data = merge_data_hist(hist_path, data_files)
    if h_data is None:
        print(f"[SKIP] No matching data histogram: {hist_path}")
        n_skip += 1
        continue

    if not same_binning(h_data, h_mc):
        print(f"[SKIP] Data/MC binning mismatch: {hist_path}")
        n_skip += 1
        continue

    if h_data.Integral() <= 0 or h_mc.Integral() <= 0:
        print(f"[SKIP] Empty histogram: {hist_path}")
        n_skip += 1
        continue

    parts = hist_path.split("/")
    hist_name = parts[-1]
    subdirs = parts[:-1]

    out_subdir = os.path.join(output_dir, *subdirs)
    os.makedirs(out_subdir, exist_ok=True)

    out_path = os.path.join(out_subdir, f"{hist_name}.png")
    pretty_title = prettify_hist_name(hist_name)

    draw_compare_plot(
        h_data=h_data,
        h_mc=h_mc,
        out_path=out_path,
        hist_title=pretty_title,
    )

    print(f"[DONE] {hist_path}")
    n_done += 1

mc_file.Close()
for f in data_files:
    f.Close()

print("")
print(f"[INFO] Finished. Done = {n_done}, Skipped = {n_skip}")
print(f"[INFO] Output directory: {output_dir}")
