#!/usr/bin/env python3
"""
plot_timewalk_fit.py

Plots and exports the Time-Walk 2D correlation, fitted curve, dual error representations
(statistical SEM points + shaded total uncertainty band), and data/fit ratio pane (Thesis Fig 4.7).

Supports single-channel export (-c 24005_4) or full detector export (-c all), where every channel
is saved under its own TDirectory inside a single consolidated ROOT file.
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
        description="Plot Time-Walk 2D correlation, fit, and ratio into ROOT file with TDirectories."
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
        help="Fixed channel ID (e.g. 24005_4, 20000_0) or 'all' for all 599 channels",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="",
        help="Output ROOT file path (default: timewalk_fit_<channel>.root or timewalk_fits_all.root)",
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


def process_channel(
    base_path, run_str, fixed_ch, target_dir=None, verbose=True
):
    """Processes a single channel and writes all graphs/canvas into target_dir or returns them."""
    root_file = find_root_file(base_path, run_str, fixed_ch)
    if not root_file:
        if verbose:
            print(
                f"Error: Could not find timewalk_{fixed_ch}.root under {base_path}"
            )
        return False

    f_in = ROOT.TFile.Open(root_file, "READ")
    if not f_in or f_in.IsZombie():
        return False

    histname = f"dtvqdc_{fixed_ch}_uncorrected"
    h2 = f_in.Get(histname)
    if not h2:
        histname += "-5cmFiducialCut"
        h2 = f_in.Get(histname)
    if not h2:
        f_in.Close()
        return False

    # 1. 1D Timing Projections per QDC bin: Extract Modal Value and SEM
    nbins = h2.GetNbinsX()
    x_pts, y_pts, ey_pts = [], [], []

    for b in range(1, nbins + 1):
        x_val = h2.GetXaxis().GetBinCenter(b)
        if x_val < 1.0:
            continue

        py = h2.ProjectionY(f"py_{fixed_ch}_{b}", b, b)
        entries = py.GetEntries()
        if entries < 10:
            continue

        max_bin = py.GetMaximumBin()
        y_val = py.GetBinCenter(max_bin)

        rms_val = py.GetRMS()
        ey_val = rms_val / np.sqrt(entries) if entries > 0 else 0.0

        if ey_val <= 0:
            continue

        x_pts.append(x_val)
        y_pts.append(y_val)
        ey_pts.append(ey_val)

    if len(x_pts) < 10:
        f_in.Close()
        return False

    x_arr = np.array(x_pts)
    y_arr = np.array(y_pts)
    ey_arr = np.array(ey_pts)

    qdc_min_fit = float(x_arr[0])
    qdc_max_fit = float(x_arr[-1])

    # Load fitted parameters from JSON if available
    json_path = os.path.join(
        base_path, f"Polyparams/run{run_str}/polyparams9_{fixed_ch}.json"
    )
    if not os.path.exists(json_path):
        json_path = os.path.join(
            base_path,
            f"sndVetoUS-physics2022/Polyparams/run{run_str}/polyparams9_{fixed_ch}.json",
        )

    popt = None
    sigma_sys = 0.040

    if os.path.exists(json_path):
        try:
            with open(json_path) as jf:
                d = json.load(jf)
                params = d["params"]
                popt = [
                    params["t_0"],
                    params["alpha"],
                    params["beta"],
                    params["QDC_0"],
                    params["gamma"],
                ]
                sigma_sys = d.get("alpha", 0.040)
        except Exception:
            popt = None

    # If parameters not loaded, fit on the fly with Table 4.2 bounds
    if popt is None:
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
        except Exception:
            popt = p0

        res_nll = minimize_scalar(
            nll_objective,
            bounds=(0.001, 1.0),
            method="bounded",
            args=(x_arr, y_arr, ey_arr, popt),
        )
        sigma_sys = float(res_nll.x) if res_nll.success else 0.040

    t0, alpha, beta, qdc0, gamma = [float(v) for v in popt]
    y_fit = tw_func(x_arr, *popt)
    ey_tot = np.sqrt(ey_arr**2 + sigma_sys**2)
    ndf = max(1, len(x_arr) - 5)

    # Chi2 before and after NLL floor
    chi2_stat = float(np.sum(((y_arr - y_fit) ** 2) / (ey_arr**2)))
    chi2_stat_ndf = chi2_stat / ndf

    chi2_tot = float(np.sum(((y_arr - y_fit) ** 2) / (ey_tot**2)))
    chi2_tot_ndf = chi2_tot / ndf
    rms_res = float(np.sqrt(np.mean((y_arr / y_fit - 1.0) ** 2)))

    ratio = y_arr / y_fit
    eratio_stat = ey_arr / y_fit
    eratio_tot = ey_tot / y_fit

    dx_band = np.zeros_like(x_arr)
    if len(x_arr) > 1:
        dx_val = 0.5 * float(np.median(np.diff(x_arr)))
        dx_band.fill(dx_val)

    if verbose:
        print("=" * 60)
        print(f"Fitted Parameters for Channel {fixed_ch}:")
        print(f"  t0:              {t0:.3f} ns")
        print(f"  alpha:           {alpha:.3f}")
        print(f"  beta:            {beta:.3f}")
        print(f"  QDC0:            {qdc0:.3f}")
        print(f"  gamma:           {gamma:.5f}")
        print(f"  sigma_sys:       {sigma_sys*1000:.1f} ps")
        print(f"  chi2_stat / ndf: {chi2_stat_ndf:.2f}  (before error floor)")
        print(f"  chi2_tot / ndf:  {chi2_tot_ndf:.2f}   (after NLL error floor)")
        print(f"  RMS residual:    {rms_res*100:.2f}%")
        print("=" * 60)

    # Setup Canvas with 2 pads (Upper 70%, Lower 30%)
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
    h2_clone = h2.Clone(f"h2_dtvqdc_{fixed_ch}")
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

    # Shaded total error band
    g_data_tot = ROOT.TGraphErrors(
        len(x_arr),
        x_arr.astype(float),
        y_arr.astype(float),
        dx_band.astype(float),
        ey_tot.astype(float),
    )
    g_data_tot.SetName("g_data_total_error")
    g_data_tot.SetFillColorAlpha(ROOT.kAzure + 1, 0.40)
    g_data_tot.SetLineColor(ROOT.kAzure + 1)
    g_data_tot.Draw("E2 SAME")

    # Statistical points (pure SEM)
    g_data_stat = ROOT.TGraphErrors(
        len(x_arr),
        x_arr.astype(float),
        y_arr.astype(float),
        np.zeros_like(x_arr),
        ey_arr.astype(float),
    )
    g_data_stat.SetName("g_data_stat_error")
    g_data_stat.SetMarkerStyle(20)
    g_data_stat.SetMarkerSize(0.7)
    g_data_stat.SetMarkerColor(ROOT.kBlack)
    g_data_stat.SetLineColor(ROOT.kBlack)
    g_data_stat.Draw("P SAME")

    # TF1 fitted curve
    fit_tf1 = ROOT.TF1(
        f"fit_tw_{fixed_ch}",
        "[0] + [1] / ([2]*x - [3]) + [4]*x",
        qdc_min_fit,
        qdc_max_fit,
    )
    fit_tf1.SetParameters(t0, alpha, beta, qdc0, gamma)
    fit_tf1.SetLineColor(ROOT.kRed + 1)
    fit_tf1.SetLineWidth(3)
    fit_tf1.Draw("SAME")

    # Legend
    leg = ROOT.TLegend(0.38, 0.52, 0.88, 0.89)
    leg.SetBorderSize(1)
    leg.SetFillColorAlpha(ROOT.kWhite, 0.88)
    leg.SetTextSize(0.030)
    leg.AddEntry(g_data_stat, "Data (Statistical SEM)", "lep")
    leg.AddEntry(g_data_tot, "Total Uncertainty (SEM #oplus #sigma_{sys})", "f")
    leg.AddEntry(
        fit_tf1,
        f"Fit: #chi^{{2}}_{{stat}}/#nu = {chi2_stat_ndf:.1f}  #rightarrow  #chi^{{2}}_{{tot}}/#nu = {chi2_tot_ndf:.2f}",
        "l",
    )
    leg.AddEntry(
        "",
        f"t_{{0}} = {t0:.2f} ns,  #alpha = {alpha:.2f},  #beta = {beta:.2f}",
        "",
    )
    leg.AddEntry(
        "",
        f"QDC_{{0}} = {qdc0:.2f},  #gamma = {gamma:.5f}",
        "",
    )
    leg.AddEntry(
        "",
        f"#sigma_{{sys}} = {sigma_sys*1000:.1f} ps,  RMS res = {rms_res*100:.2f}%",
        "",
    )
    leg.Draw()

    # Draw Lower Pad (Ratio)
    pad2.cd()
    g_ratio_tot = ROOT.TGraphErrors(
        len(x_arr),
        x_arr.astype(float),
        ratio.astype(float),
        dx_band.astype(float),
        eratio_tot.astype(float),
    )
    g_ratio_tot.SetName("g_ratio_total_error")
    g_ratio_tot.SetTitle(";QDC_{SiPM} [a.u.];Data / Fit")
    g_ratio_tot.SetFillColorAlpha(ROOT.kAzure + 1, 0.40)
    g_ratio_tot.SetLineColor(ROOT.kAzure + 1)

    g_ratio_tot.GetXaxis().SetRangeUser(0, max(50.0, qdc_max_fit * 1.1))
    g_ratio_tot.GetYaxis().SetRangeUser(0.92, 1.08)
    g_ratio_tot.GetXaxis().SetTitleSize(0.11)
    g_ratio_tot.GetXaxis().SetLabelSize(0.09)
    g_ratio_tot.GetXaxis().SetTitleOffset(1.0)
    g_ratio_tot.GetYaxis().SetTitleSize(0.10)
    g_ratio_tot.GetYaxis().SetLabelSize(0.08)
    g_ratio_tot.GetYaxis().SetTitleOffset(0.55)
    g_ratio_tot.GetYaxis().SetNdivisions(505)

    g_ratio_tot.Draw("AE2")

    g_ratio_stat = ROOT.TGraphErrors(
        len(x_arr),
        x_arr.astype(float),
        ratio.astype(float),
        np.zeros_like(x_arr),
        eratio_stat.astype(float),
    )
    g_ratio_stat.SetName("g_ratio_stat_error")
    g_ratio_stat.SetMarkerStyle(20)
    g_ratio_stat.SetMarkerSize(0.7)
    g_ratio_stat.SetMarkerColor(ROOT.kBlack)
    g_ratio_stat.SetLineColor(ROOT.kBlack)
    g_ratio_stat.Draw("P SAME")

    line1 = ROOT.TLine(0, 1.0, max(50.0, qdc_max_fit * 1.1), 1.0)
    line1.SetLineColor(ROOT.kRed + 1)
    line1.SetLineStyle(2)
    line1.SetLineWidth(2)
    line1.Draw("SAME")

    c.Update()

    # Write objects into target_dir if provided
    if target_dir is not None:
        target_dir.cd()
        c.Write(f"c_timewalk_fit_{fixed_ch}")
        h2_clone.Write("h2_dtvqdc")
        g_data_stat.Write("g_data_stat_error")
        g_data_tot.Write("g_data_total_error")
        fit_tf1.Write("fit_tw")
        g_ratio_stat.Write("g_ratio_stat_error")
        g_ratio_tot.Write("g_ratio_total_error")
        line1.Write("line_unity")

    f_in.Close()
    return c, h2_clone, g_data_stat, g_data_tot, fit_tf1, g_ratio_stat, g_ratio_tot, line1


def main():
    args = parse_args()
    run_str = f"{args.runNumber:06d}"
    base_path = args.path.rstrip("/") + "/"
    mode_ch = args.channel.strip().lower()

    left_sipms = [0, 1, 3, 4, 6, 7]
    right_sipms = [8, 9, 11, 12, 14, 15]
    all_large_sipms = left_sipms + right_sipms

    ROOT.TH1.AddDirectory(False)

    if mode_ch in ["all", "all_channels", "*"]:
        # Run across all 599 channels into one ROOT file with per-channel TDirectories
        if not args.output:
            out_root_dir = os.path.join(base_path, f"rootfiles/run{run_str}")
            os.makedirs(out_root_dir, exist_ok=True)
            out_root = os.path.join(out_root_dir, "timewalk_fits_all.root")
        else:
            out_root = args.output
            if not out_root.endswith(".root"):
                out_root += ".root"
            os.makedirs(
                os.path.dirname(os.path.abspath(out_root)), exist_ok=True
            )

        print("=" * 60)
        print("SND@LHC US Time-Walk Canvas Generation for ALL Channels")
        print(f"  Run Number:  {run_str}")
        print(f"  Target File: {out_root}")
        print("=" * 60)

        f_out = ROOT.TFile.Open(out_root, "RECREATE")
        ch_count = 0
        total_ch = 5 * 10 * 12

        for plane in range(5):
            for bar in range(10):
                detID = int(f"2{plane}00{bar}")
                for sipm in all_large_sipms:
                    fixed_ch = f"{detID}_{sipm}"
                    if fixed_ch == "21005_15":
                        continue

                    # Create TDirectory for this channel
                    dir_ch = f_out.mkdir(fixed_ch)
                    res = process_channel(
                        base_path,
                        run_str,
                        fixed_ch,
                        target_dir=dir_ch,
                        verbose=False,
                    )
                    if res:
                        ch_count += 1
                        if ch_count % 50 == 0 or ch_count == 599:
                            print(
                                f"  [{ch_count:3d}/599] Processed channel {fixed_ch} -> TDirectory '{fixed_ch}/'"
                            )

        f_out.Close()
        print("\n" + "=" * 60)
        print(f"Successfully processed and stored {ch_count} channels into:")
        print(f"  {out_root}")
        print("  Each channel has its own TDirectory with TCanvas and graphs.")
        print("=" * 60)

    else:
        # Single channel mode
        fixed_ch = args.channel
        if not args.output:
            out_root_dir = os.path.join(base_path, f"rootfiles/run{run_str}")
            os.makedirs(out_root_dir, exist_ok=True)
            out_root = os.path.join(
                out_root_dir, f"timewalk_fit_{fixed_ch}.root"
            )
        else:
            out_root = args.output
            if not out_root.endswith(".root"):
                out_root += ".root"
            os.makedirs(
                os.path.dirname(os.path.abspath(out_root)), exist_ok=True
            )

        f_out = ROOT.TFile.Open(out_root, "RECREATE")
        res = process_channel(
            base_path, run_str, fixed_ch, target_dir=f_out, verbose=True
        )
        if res:
            f_out.Close()
            print(f"\nROOT File successfully saved to:\n  {out_root}")
            print(
                f"Contains: TCanvas 'c_timewalk_fit_{fixed_ch}', TH2F, TGraphErrors, TF1, TLine."
            )
        else:
            f_out.Close()
            print(f"Failed to process channel {fixed_ch}.")


if __name__ == "__main__":
    main()
