#!/usr/bin/env python3
"""
plot_timewalk_fit.py

Reproduces Figure 4.7 from Andrew Conaboy's PhD dissertation:
Upper panel: 2D dtvqdc histogram, profiled data points with error bars, and the fitted
             5-parameter rational time-walk curve overlaid.
Lower panel: Ratio of data points to the fitted function (y_data / y_fit) vs QDC.

Saves all objects (TH2F, TGraphErrors, TF1, TLine) and the combined TCanvas into a .root file.
"""

import os
import sys
import json
import argparse
import numpy as np
from scipy.optimize import curve_fit, minimize_scalar
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kRainBow)


def tw_func(qdc, t0, alpha, beta, qdc0, gamma):
    """
    Time-walk parameterisation function (Thesis eq. 4.1.9):
    f(QDC) = t0 + alpha / (beta * QDC - qdc0) + gamma * QDC
    """
    denom = beta * qdc - qdc0
    denom = np.where(np.abs(denom) < 1e-4, 1e-4, denom)
    return t0 + alpha / denom + gamma * qdc


def nll_objective(sigma_sys, x_data, y_data, sigma_stat, params):
    y_pred = tw_func(x_data, *params)
    var_total = sigma_stat**2 + sigma_sys**2
    return 0.5 * np.sum(((y_data - y_pred) ** 2) / var_total + np.log(var_total))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot single-channel Time-Walk 2D correlation, fit, and ratio into a ROOT file (Thesis Figure 4.7)."
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
        help="Output ROOT file path (default: <path>/rootfiles/run<runNr>/timewalk_fit_<channel>.root)",
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

    # 1. Profile along QDC with SEM error
    prof = h2.ProfileX(f"prof_{fixed_ch}", 1, -1, "")
    nbins = prof.GetNbinsX()

    x_pts, y_pts, ey_pts = [], [], []
    for b in range(1, nbins + 1):
        x_val = prof.GetBinCenter(b)
        entries = prof.GetBinEntries(b)
        if x_val < 1.0 or entries < 10 or prof.GetBinError(b) <= 0:
            continue
        x_pts.append(x_val)
        y_pts.append(prof.GetBinContent(b))
        ey_pts.append(prof.GetBinError(b))

    if len(x_pts) < 10:
        print(f"Error: Insufficient points ({len(x_pts)}) for channel {fixed_ch}")
        f_in.Close()
        return

    x_arr = np.array(x_pts)
    y_arr = np.array(y_pts)
    ey_arr = np.array(ey_pts)

    qdc_min_fit = float(x_arr[0])
    qdc_max_fit = float(x_arr[-1])

    # 2. Fit with Thesis Table 4.2 parameter bounds
    high_qdc_mask = x_arr > (0.6 * qdc_max_fit)
    t0_est = (
        float(np.median(y_arr[high_qdc_mask]))
        if np.sum(high_qdc_mask) > 3
        else float(y_arr[-1])
    )

    p0 = [t0_est, -50.0, 5.0, -2.0, 0.005]
    bounds_lower = [
        0.5 * t0_est if t0_est > 0 else 1.5 * t0_est,
        -500.0,
        0.01,
        -10.0,
        0.0,
    ]
    bounds_upper = [
        1.5 * t0_est if t0_est > 0 else 0.5 * t0_est,
        0.0,
        30.0,
        0.0,
        1.0,
    ]

    try:
        popt, _ = curve_fit(
            tw_func,
            x_arr,
            y_arr,
            p0=p0,
            bounds=(bounds_lower, bounds_upper),
            maxfev=10000,
        )
    except Exception as e:
        print(f"Fit error: {e}")
        popt = p0

    # 3. NLL error floor extraction
    res_nll = minimize_scalar(
        nll_objective,
        bounds=(0.001, 1.0),
        method="bounded",
        args=(x_arr, y_arr, ey_arr, popt),
    )
    sigma_sys = float(res_nll.x) if res_nll.success else 0.350

    t0, alpha, beta, qdc0, gamma = popt
    y_fit = tw_func(x_arr, *popt)
    ey_tot = np.sqrt(ey_arr**2 + sigma_sys**2)
    ndf = max(1, len(x_arr) - 5)
    chi2 = float(np.sum(((y_arr - y_fit) ** 2) / (ey_tot**2)))
    chi2_ndf = chi2 / ndf
    rms_res = float(np.sqrt(np.mean((y_arr / y_fit - 1.0) ** 2)))

    ratio = y_arr / y_fit
    eratio = ey_tot / y_fit

    print("=" * 60)
    print(f"Fitted Parameters for Channel {fixed_ch}:")
    print(f"  t0:        {t0:.3f} ns")
    print(f"  alpha:     {alpha:.3f}")
    print(f"  beta:      {beta:.3f}")
    print(f"  QDC0:      {qdc0:.3f}")
    print(f"  gamma:     {gamma:.5f}")
    print(f"  sigma_sys: {sigma_sys*1000:.1f} ps")
    print(f"  chi2/ndf:  {chi2_ndf:.2f}")
    print(f"  RMS res:   {rms_res*100:.2f}%")
    print("=" * 60)

    # 4. Setup ROOT Canvas with 2 pads (Upper 70%, Lower 30%)
    ROOT.TH1.AddDirectory(False)
    c = ROOT.TCanvas(
        f"c_timewalk_fit_{fixed_ch}",
        f"US Time-Walk Fit - Channel {fixed_ch} (Run {run_str})",
        1000,
        900,
    )

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
    h2_clone = h2.Clone(f"h2_{fixed_ch}")
    h2_clone.SetTitle(
        f"US Channel {fixed_ch} (Run {run_str});;dt^{{ToF}}_{{SiPM}} [ns]"
    )
    h2_clone.GetXaxis().SetRangeUser(0, max(50.0, qdc_max_fit * 1.1))
    h2_clone.GetXaxis().SetLabelSize(0)
    h2_clone.GetYaxis().SetTitleSize(0.05)
    h2_clone.GetYaxis().SetTitleOffset(1.1)
    h2_clone.GetYaxis().SetRangeUser(
        min(y_arr) - 3.0 * sigma_sys - 0.5, max(y_arr) + 1.5
    )
    h2_clone.Draw("COLZ")

    # Data points graph
    g_data = ROOT.TGraphErrors(
        len(x_arr),
        x_arr.astype(float),
        y_arr.astype(float),
        np.zeros_like(x_arr),
        ey_tot.astype(float),
    )
    g_data.SetName("g_data_profile")
    g_data.SetMarkerStyle(20)
    g_data.SetMarkerSize(0.8)
    g_data.SetMarkerColor(ROOT.kBlack)
    g_data.SetLineColor(ROOT.kBlack)
    g_data.Draw("P SAME")

    # Fitted TF1 curve
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
    leg = ROOT.TLegend(0.46, 0.65, 0.86, 0.88)
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
    g_ratio.SetName("g_ratio")
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

    # 5. Output ROOT File
    if not args.output:
        out_root_dir = os.path.join(base_path, f"rootfiles/run{run_str}")
        os.makedirs(out_root_dir, exist_ok=True)
        out_root = os.path.join(out_root_dir, f"timewalk_fit_{fixed_ch}.root")
    else:
        out_root = args.output
        if not out_root.endswith(".root"):
            out_root += ".root"
        os.makedirs(os.path.dirname(os.path.abspath(out_root)), exist_ok=True)

    f_out = ROOT.TFile.Open(out_root, "RECREATE")
    c.Write()
    h2_clone.Write("h2_dtvqdc")
    prof.Write("prof_dtvqdc")
    g_data.Write("g_data_profile")
    fit_tf1.Write("fit_tw")
    g_ratio.Write("g_ratio")
    line1.Write("line_unity")
    f_out.Close()
    f_in.Close()

    print(f"ROOT File successfully saved to:\n  {out_root}")
    print("Contains: TCanvas 'c_fig4_7', TH2F, TProfile, TGraphErrors, TF1, TLine.")


if __name__ == "__main__":
    main()
