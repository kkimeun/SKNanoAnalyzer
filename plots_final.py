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
output_dir = os.path.join(base_dir, "compare_data_stacked_mc")
os.makedirs(output_dir, exist_ok=True)

# -------------------------
# data files
# -------------------------
eras = ["C", "D", "E", "F", "G", "H", "I"]
data_paths = []
for era in eras:
    data_paths.append(os.path.join(base_dir, f"Muon0_{era}.root"))
    data_paths.append(os.path.join(base_dir, f"Muon1_{era}.root"))

# -------------------------
# MC samples to stack
# -------------------------
mc_configs = [
    {
        "file": "ST_tch_antitop_lep.root",
        "legend": "ST t-ch antitop lep",
        "color": cmsstyle.p10.kBrown,
    },
    {
        "file": "ST_tch_top_lep.root",
        "legend": "ST t-ch top lep",
        "color": cmsstyle.p10.kOrange,
    },
    {
        "file": "ST_tW_antitop_Lep.root",
        "legend": "ST tW antitop Lep",
        "color": cmsstyle.p10.kYellow,
    },
    {
        "file": "ST_tW_antitop_Semilep.root",
        "legend": "ST tW antitop Semilep",
        "color": cmsstyle.p10.kAsh,
    },
    {
        "file": "ST_tW_top_Lep.root",
        "legend": "ST tW top Lep",
        "color": cmsstyle.p10.kCyan,
    },
    {
        "file": "ST_tW_top_Semilep.root",
        "legend": "ST tW top Semilep",
        "color": cmsstyle.p10.kViolet,
    },
    {
        "file": "TTLJ_powheg.root",
        "legend": "TT LJ",
        "color": cmsstyle.p10.kRed,
    },
    {
        "file": "TTLL_powheg.root",
        "legend": "TT LL",
        "color": cmsstyle.p10.kBlue,
    },
]

# =========================================================
# CMS style setup
# =========================================================
cmsstyle.setCMSStyle()
cmsstyle.SetLumi(110.0, run=None)
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

opened_mc_files = []
for cfg in mc_configs:
    full_path = os.path.join(base_dir, cfg["file"])
    f = ROOT.TFile.Open(full_path)
    if not f or f.IsZombie():
        print(f"[WARNING] Cannot open MC file: {full_path}")
        continue
    cfg["handle"] = f
    cfg["full_path"] = full_path
    opened_mc_files.append(cfg)
    print(f"[INFO] Opened MC file: {full_path}")

if not opened_mc_files:
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

def same_binning(h1, h2):
    return (
        h1.GetNbinsX() == h2.GetNbinsX()
        and abs(h1.GetXaxis().GetXmin() - h2.GetXaxis().GetXmin()) < 1e-9
        and abs(h1.GetXaxis().GetXmax() - h2.GetXaxis().GetXmax()) < 1e-9
    )

def clone_hist(h, name):
    out = h.Clone(name)
    out.SetDirectory(0)
    out.SetStats(0)
    return out

def merge_data_hist(path, files):
    merged = None
    for i, f in enumerate(files):
        h = get_obj_from_path(f, path)
        if not h or not is_hist(h):
            continue

        if merged is None:
            merged = clone_hist(h, path.replace("/", "__") + "_data")
        else:
            if not same_binning(merged, h):
                print(f"[WARNING] Data binning mismatch for {path} in {f.GetName()}, skip")
                continue
            merged.Add(h)

    return merged

def get_mc_hists_for_path(path, mc_cfgs):
    mc_hists = []
    for i, cfg in enumerate(mc_cfgs):
        h = get_obj_from_path(cfg["handle"], path)
        if not h or not is_hist(h):
            continue

        h = clone_hist(h, path.replace("/", "__") + f"_mc_{i}")
        if h.Integral() <= 0:
            continue

        mc_hists.append((h, cfg))
    return mc_hists

def build_total_mc(mc_hists, name):
    total = None
    for h, cfg in mc_hists:
        if total is None:
            total = clone_hist(h, name)
        else:
            if not same_binning(total, h):
                return None
            total.Add(h)
    return total

def make_ratio_hist(h_data, h_mc_total, name):
    ratio = clone_hist(h_data, name)
    for i in range(1, ratio.GetNbinsX() + 1):
        d = h_data.GetBinContent(i)
        de = h_data.GetBinError(i)
        m = h_mc_total.GetBinContent(i)

        if m > 0:
            ratio.SetBinContent(i, d / m)
            ratio.SetBinError(i, de / m)
        else:
            ratio.SetBinContent(i, 0.0)
            ratio.SetBinError(i, 0.0)
    return ratio

def make_mc_unc_band(h_mc_total, name):
    band = clone_hist(h_mc_total, name)
    for i in range(1, band.GetNbinsX() + 1):
        c = h_mc_total.GetBinContent(i)
        e = h_mc_total.GetBinError(i)
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

def draw_stack_with_ratio(h_data, mc_hists, h_mc_total, out_path, label=""):
    if not mc_hists or not h_mc_total:
        return

    # build THStack
    stack_hists = [h for h, cfg in mc_hists]
    stack_colors = [cfg["color"] for h, cfg in mc_hists]
    stack_labels = [cfg["legend"] for h, cfg in mc_hists]

    hs = cmsstyle.buildTHStack(
        stack_hists,
        stack_colors,
        LineColor=-1,
        FillColor=-1
    )

    ymax = max(h_data.GetMaximum(), h_mc_total.GetMaximum()) * 1.8
    xmin = h_data.GetXaxis().GetXmin()
    xmax = h_data.GetXaxis().GetXmax()
    xtitle = h_data.GetXaxis().GetTitle()
    if not xtitle:
        xtitle = "plot"

    c = cmsstyle.cmsDiCanvas(
        "c",
        xmin, xmax,
        0.0, ymax,
        0.5, 1.5,
        xtitle,
        "Events / bin",
        "Data / MC"
    )

    # legend
    plotlegend = cmsstyle.cmsLeg(0.45, 0.50, 0.93, 0.88, textSize=0.035, columns=2)
    cmsstyle.addToLegend(plotlegend, (h_data, "Data", "lpe"))
    for h, cfg in mc_hists:
        cmsstyle.addToLegend(plotlegend, (h, cfg["legend"], "f"))

    h_mc_unc = clone_hist(h_mc_total, "h_mc_unc")
    style_unc(h_mc_unc)
    cmsstyle.addToLegend(plotlegend, (h_mc_unc, "Stat. Unc.", "f"))

    # ---------------- upper pad ----------------
    c.cd(1)

    cmsstyle.cmsObjectDraw(hs, "HIST")
    cmsstyle.GetCmsCanvasHist(ROOT.gPad).GetYaxis().SetMaxDigits(3)

    cmsstyle.cmsObjectDraw(
        h_mc_unc,
        "E2",
        FillStyle=3345,
        LineWidth=0,
        FillColor=12,
        MarkerSize=0
    )

    style_data(h_data)
    cmsstyle.cmsObjectDraw(h_data, "E", MarkerStyle=ROOT.kFullCircle)

    if label:
        extraLabel = ROOT.TLatex(0.18, 0.70, label)
        extraLabel.SetNDC()
        cmsstyle.cmsObjectDraw(extraLabel, TextFont=cmsstyle.additionalInfoFont)

    plotlegend.Draw()
    cmsstyle.UpdatePad()

    # ---------------- lower pad ----------------
    c.cd(2)

    data_ratio = make_ratio_hist(h_data, h_mc_total, "data_ratio")
    pred_ratio = clone_hist(h_mc_total, "pred_ratio")
    pred_ratio.Divide(h_mc_total)

    ratio_unc = make_mc_unc_band(h_mc_total, "ratio_unc")
    style_unc(ratio_unc)

    style_data(data_ratio)

    cmsstyle.cmsObjectDraw(
        ratio_unc,
        "E2",
        FillStyle=3345,
        LineWidth=0,
        FillColor=12,
        MarkerSize=0
    )
    cmsstyle.cmsObjectDraw(data_ratio, "E", MarkerStyle=ROOT.kFullCircle)

    cmsstyle.UpdatePad(c)
    c.SaveAs(out_path)

# =========================================================
# main loop
# =========================================================
reference_mc_file = opened_mc_files[0]["handle"]
mc_hist_paths = collect_hist_paths(reference_mc_file)

print(f"[INFO] Found {len(mc_hist_paths)} MC histograms in reference file")

n_done = 0
n_skip = 0

for path in mc_hist_paths:
    h_data = merge_data_hist(path, opened_data_files)
    if not h_data:
        print(f"[SKIP] No matching data hist: {path}")
        n_skip += 1
        continue

    mc_hists = get_mc_hists_for_path(path, opened_mc_files)
    if not mc_hists:
        print(f"[SKIP] No MC histograms found: {path}")
        n_skip += 1
        continue

    h_mc_total = build_total_mc(mc_hists, path.replace("/", "__") + "_mc_total")
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
        mc_hists=mc_hists,
        h_mc_total=h_mc_total,
        out_path=out_path,
        label=""
    )

    print(f"[DONE] {path}")
    n_done += 1

# =========================================================
# close files
# =========================================================
for cfg in opened_mc_files:
    cfg["handle"].Close()

for f in opened_data_files:
    f.Close()

print(f"\n[INFO] Finished. Done = {n_done}, Skipped = {n_skip}")
print(f"[INFO] Output directory: {output_dir}")
