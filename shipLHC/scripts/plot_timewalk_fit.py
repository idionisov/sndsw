#!/usr/bin/env python3
"""
plot_timewalk_fit.py

Reproduces Figure 4.7 from Andrew Conaboy's PhD dissertation:
Upper panel: 2D dtvqdc histogram, profiled data points with error bars, and the fitted
             5-parameter rational time-walk curve overlaid.
Lower panel: Ratio of data points to the fitted function (y_data / y_fit) vs QDC.
"""

import os
import sys
import json
import argparse
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kRainBow)


def tw_func(qdc, t0, alpha, beta, qdc0, gamma):
    denom = np.maximum(beta * qdc - qdc0, 1e-4)
    return t0 + alpha / denom + gamma * qdc


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot single-channel Time-Walk 2D correlation, fit, and ratio (Thesis Figure 4.7)."
    )
    parser.add_argument(
        "-r", "--runNumber", type=int, default=6640, help="Run number (e.g. 6640)"
    )
    parser.add_argument(
        "-p",
        "--path",
        type=str,
        default="/afs/cern.ch/work/i/idioniso/sndVetoUS-physics2022/",
        help="Base path to physics directory",
    )
    parser.add_argument(
        "-c",
        "--channel",
        type=str,
        default="24005_4",
        help="Fixed channel ID (e.g. 24005_4, 20000_0, 22004_1)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="",
        help="Output PNG/PDF path (default: <path>/plots/run<runNr>/tw/fig4_7_<channel>.png)",
    )
    return parser.parse_args()


def find_root_file(base_path, run_str, fixed_ch):
    candidates = [
        os.path.join(
            base_path,
            f"sndVetoUS-physics2022/rootfiles/run{run_str}/timewalk_{fixed_ch}.root",
        ),
        os.path.join(
            base_path, f"rootfiles/run{run_str}/timewalk_{fixed_ch}.root"
        ),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return None


def main():
    args = parse_args()
    run_str = f"{args.runNumber:06d}"
    base_path = args.path.rstrip("/") + "/"
    fixed_ch = args.channel

    root_file = find_root_file(base_path, run_str, fixed_ch)
    if not root_file:
        print(
            f"Error: Could not find timewalk_{fixed_ch}.root under {base_path}"
        )
        return

    f_in = ROOT.TFile.Open(root_file, "READ")
    histname = f"dtvqdc_{fixed_ch}_uncorrected"
    h2 = f_in.Get(histname)
    if not h2:
        histname += "-5cmFiducialCut"
        h2 = f_in.Get(histname)
    if not h2:
        print(f"Error: Histogram {histname} not found in {root_file}")
        f_in.Close()
        return

    # Load fitted parameters from JSON
    json_path = os.path.join(
        base_path, f"Polyparams/run{run_str}/polyparams9_{fixed_ch}.json"
    )
    if not os.path.exists(json_path):
        json_path = os.path.join(
            base_path,
            f"sndVetoUS-physics2022/Polyparams/run{run_str}/polyparams9_{fixed_ch}.json",
        )

    t0, alpha, beta, qdc0, gamma = 0.0, 2.0, 1.0, 0.1, 0.005
    sigma_sys = 0.350
    chi2_ndf = 1.0
    rms_res = 0.04
    qdc_min_fit, qdc_max_fit = 1.0, 40.0

    if os.path.exists(json_path):
        with open(json_path) as jf:
            d = json.load(jf)
            params = d["params"]
            t0 = params["t_0"]
            alpha = params["alpha"]
            beta = params["beta"]
            qdc0 = params["QDC_0"]
            gamma = params["gamma"]
            sigma_sys = d.get("alpha", 0.350)
            chi2 = d.get("chi2", 1.0)
            ndf = d.get("ndf", 1)
            chi2_ndf = chi2 / max(1, ndf)
            rms_res = d.get("rms_residual", 0.04)
            qdc_min_fit, qdc_max_fit = d.get("limits", [1.0, 40.0])

    # Extract 1D profile
    prof = h2.ProfileX(f"prof_{fixed_ch}", 1, -1, "")
    nbins = prof.GetNbinsX()

    x_pts, y_pts, ey_pts = [], [], []
    for b in range(1, nbins + 1):
        x_val = prof.GetBinCenter(b)
        entries = prof.GetBinEntries(b)
        if (
            x_val < qdc_min_fit
            or x_val > qdc_max_fit
            or entries < 10
            or prof.GetBinError(b) <= 0
        ):
            continue
        x_pts.append(x_val)
        y_pts.append(prof.GetBinContent(b))
        ey_pts.append(np.sqrt(prof.GetBinError(b) ** 2 + sigma_sys**2))

    x_arr = np.array(x_pts)
    y_arr = np.array(y_pts)
    ey_arr = np.array(ey_pts)

    y_fit = tw_func(x_arr, t0, alpha, beta, qdc0, gamma)
    ratio = y_arr / y_fit
    eratio = ey_arr / y_fit

    # Setup Canvas with 2 vertical pads (Upper 70%, Lower 30%)
    c = ROOT.TCanvas("c_fig4_7", f"Figure 4.7 Reproduction - {fixed_ch}", 1000, 900)

    pad1 = ROOT.TPad("pad1", "Upper Pad", 0.0, 0.30, 1.0, 1.0)
    pad2 = ROOT.TPad("pad2", "Lower Pad", 0.0, 0.0, 1.0, 0.30)

    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.12)
    pad1.SetRightMargin(0.12)
    pad1.SetGrid()

    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.30)
    pad2.SetLeftMargin(0.12)
    pad2.SetRightMargin(0.12)
    pad2.SetGrid()

    pad1.Draw()
    pad2.Draw()

    # Draw Upper Pad
    pad1.cd()
    h2.SetTitle(
        f"US Channel {fixed_ch} (Run {run_str});;dt^{{ToF}}_{{SiPM}} [ns]"
    )
    h2.GetXaxis().SetRangeUser(0, max(50.0, qdc_max_fit * 1.1))
    h2.GetXaxis().SetLabelSize(0)
    h2.GetYaxis().SetTitleSize(0.05)
    h2.GetYaxis().SetTitleOffset(1.1)
    h2.GetYaxis().SetRangeUser(
        min(y_arr) - 3.0 * sigma_sys - 0.5, max(y_arr) + 1.5
    )
    h2.Draw("COLZ")

    # Graph of profiled data points
    g_data = ROOT.TGraphErrors(
        len(x_arr),
        x_arr.astype(float),
        y_arr.astype(float),
        np.zeros_like(x_arr),
        ey_arr.astype(float),
    )
    g_data.SetMarkerStyle(20)
    g_data.SetMarkerSize(0.8)
    g_data.SetMarkerColor(ROOT.kBlack)
    g_data.SetLineColor(ROOT.kBlack)
    g_data.Draw("P SAME")

    # TF1 for fitted curve
    fit_tf1 = ROOT.TF1(
        "fit_tw",
        "[0] + [1] / ([2]*x - [3]) + [4]*x",
        qdc_min_fit,
        qdc_max_fit,
    )
    fit_tf1.SetParameters(t0, alpha, beta, qdc0, gamma)
    fit_tf1.SetLineColor(ROOT.kRed + 1)
    fit_tf1.SetLineWidth(3)
    fit_tf1.Draw("SAME")

    # Legend
    leg = ROOT.TLegend(0.48, 0.65, 0.86, 0.88)
    leg.SetBorderSize(1)
    leg.SetFillColorAlpha(ROOT.kWhite, 0.85)
    leg.AddEntry(g_data, "Data (SEM #oplus #sigma_{sys})", "lep")
    leg.AddEntry(
        fit_tf1,
        f"Fit: #chi^{{2}}/#nu = {chi2_ndf:.2f}, #sigma_{{sys}} = {sigma_sys*1000:.0f} ps",
        "l",
    )
    leg.AddEntry(
        "",
        f"t_{{0}}={t0:.2f} ns, #alpha={alpha:.2f}, #beta={beta:.2f}, QDC_{{0}}={qdc0:.2f}",
        "",
    )
    leg.Draw()

    # Draw Lower Pad (Ratio)
    pad2.cd()
    g_ratio = ROOT.TGraphErrors(
        len(x_arr),
        x_arr.astype(float),
        ratio.astype(float),
        np.zeros_like(x_arr),
        eratio.astype(float),
    )
    g_ratio.SetTitle(";QDC_{SiPM} [a.u.];Data / Fit")
    g_ratio.SetMarkerStyle(20)
    g_ratio.SetMarkerSize(0.8)
    g_ratio.SetMarkerColor(ROOT.kBlack)
    g_ratio.SetLineColor(ROOT.kBlack)

    # Set ratio axes
    g_ratio.GetXaxis().SetRangeUser(0, max(50.0, qdc_max_fit * 1.1))
    g_ratio.GetYaxis().SetRangeUser(0.92, 1.08)
    g_ratio.GetXaxis().SetTitleSize(0.11)
    g_ratio.GetXaxis().SetLabelSize(0.09)
    g_ratio.GetXaxis().SetTitleOffset(1.0)
    g_ratio.GetYaxis().SetTitleSize(0.10)
    g_ratio.GetYaxis().SetLabelSize(0.08)
    g_ratio.GetYaxis().SetTitleOffset(0.55)
    g_ratio.GetYaxis().SetNdivisions(505)

    g_ratio.Draw("AP")

    # Reference line at ratio = 1.0
    line1 = ROOT.TLine(0, 1.0, max(50.0, qdc_max_fit * 1.1), 1.0)
    line1.SetLineColor(ROOT.kRed + 1)
    line1.SetLineStyle(2)
    line1.SetLineWidth(2)
    line1.Draw("SAME")

    c.Update()

    if not args.output:
        out_plot_dir = os.path.join(base_path, f"plots/run{run_str}/tw")
        os.makedirs(out_plot_dir, exist_ok=True)
        out_png = os.path.join(out_plot_dir, f"fig4_7_{fixed_ch}.png")
    else:
        out_png = args.output

    c.SaveAs(out_png)
    print(f"Figure 4.7 reproduced and saved to: {out_png}")
    f_in.Close()


if __name__ == "__main__":
    main()
