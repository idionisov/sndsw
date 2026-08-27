#!/usr/bin/env python3
"""
extract_twparams.py

Extracts the 5-parameter Time-Walk (TW) calibration curves from merged 2D
correlation histograms (dtvqdc) produced in mode 'tof' for the Upstream (US)
system. Implements the Two-Stage Fit and Negative Log-Likelihood (NLL) systematic
error-floor minimization described in the SND@LHC timing thesis (Sections 4.1.4 & 4.3.4).

Outputs JSON and CSV calibration files formatted for direct use by AnalysisFunctions.py.
"""

import os
import sys
import json
import csv
import argparse
import numpy as np
from scipy.optimize import curve_fit, minimize_scalar
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


def tw_func(qdc, t0, alpha, beta, qdc0, gamma):
    """
    Time-walk parameterisation function (Thesis eq. 4.1.9):
    f(QDC) = t0 + alpha / (beta * QDC - qdc0) + gamma * QDC
    """
    denom = beta * qdc - qdc0
    # Guard against division by zero or negative denominator
    denom = np.where(denom <= 1e-4, 1e-4, denom)
    return t0 + alpha / denom + gamma * qdc


def nll_objective(sigma_sys, x_data, y_data, sigma_stat, params):
    """
    Negative Log-Likelihood objective function (Thesis eq. 4.3.6):
    NLL(sigma_sys) = 0.5 * sum [ (y_i - f(x_i))^2 / (sigma_stat_i^2 + sigma_sys^2)
                                 + ln(sigma_stat_i^2 + sigma_sys^2) ]
    """
    y_pred = tw_func(x_data, *params)
    var_total = sigma_stat**2 + sigma_sys**2
    nll = 0.5 * np.sum(((y_data - y_pred) ** 2) / var_total + np.log(var_total))
    return nll


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract Time-Walk (TW) parameters from mode 'tof' ROOT files."
    )
    parser.add_argument(
        "-r", "--runNumber", type=int, default=6640, help="Run number (e.g. 6640)"
    )
    parser.add_argument(
        "-p",
        "--path",
        type=str,
        default="/afs/cern.ch/work/i/idioniso/sndVetoUS-physics2022/",
        help="Base path to physics data directory containing rootfiles/",
    )
    parser.add_argument(
        "--fiducialCut",
        action="store_true",
        default=False,
        help="Use 5cm fiducial cut histogram (dtvqdc_...-5cmFiducialCut)",
    )
    parser.add_argument(
        "--qdcMin",
        type=float,
        default=1.0,
        help="Lower bound for QDC fit range (default: 1.0 a.u.)",
    )
    parser.add_argument(
        "--minBinEntries",
        type=int,
        default=10,
        help="Minimum entries per QDC bin to include in fit (default: 10)",
    )
    parser.add_argument(
        "--minTotalEntries",
        type=int,
        default=500,
        help="Minimum total histogram entries required (default: 500)",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate validation and summary plots",
    )
    parser.add_argument(
        "--plotDir",
        type=str,
        default="",
        help="Directory to save generated plots (default: <path>/plots/run<runNr>/tw/)",
    )
    return parser.parse_args()


def find_root_file(base_path, run_str, fixed_ch):
    """Search for channel ROOT file in standard and nested directories."""
    candidates = [
        os.path.join(base_path, f"rootfiles/run{run_str}/timewalk_{fixed_ch}.root"),
        os.path.join(
            base_path,
            f"sndVetoUS-physics2022/rootfiles/run{run_str}/timewalk_{fixed_ch}.root",
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
    out_dir = os.path.join(base_path, f"Polyparams/run{run_str}")
    os.makedirs(out_dir, exist_ok=True)

    if args.plotDir:
        plot_dir = args.plotDir
    else:
        plot_dir = os.path.join(base_path, f"plots/run{run_str}/tw")
    if args.plot:
        os.makedirs(plot_dir, exist_ok=True)

    print("=" * 60)
    print("SND@LHC US Time-Walk Parameter Extraction & NLL Minimisation")
    print(f"  Run Number:      {run_str}")
    print(f"  Base Path:       {base_path}")
    print(f"  Output Dir:      {out_dir}")
    print(f"  Fiducial Cut:    {args.fiducialCut}")
    print(f"  Min QDC:         {args.qdcMin} a.u.")
    print("=" * 60)

    # Large SiPM channels
    left_sipms = [0, 1, 3, 4, 6, 7]
    right_sipms = [8, 9, 11, 12, 14, 15]
    all_large_sipms = left_sipms + right_sipms

    summary_twparams = {}
    all_sigma_sys = []
    all_chi2_ndf = []
    all_rms_residuals = []

    n_processed = 0
    n_failed = 0

    # For multi-graph / diagnostic display
    sample_channels = []
    sample_fit_data = {}

    for plane in range(5):
        for bar in range(10):
            detID = int(f"2{plane}00{bar}")
            for sipm in all_large_sipms:
                fixed_ch = f"{detID}_{sipm}"

                if fixed_ch == "21005_15":
                    continue

                root_file = find_root_file(base_path, run_str, fixed_ch)
                if not root_file:
                    continue

                f = ROOT.TFile.Open(root_file, "READ")
                if not f or f.IsZombie():
                    continue

                histname = f"dtvqdc_{fixed_ch}_uncorrected"
                if args.fiducialCut:
                    histname += "-5cmFiducialCut"

                h2 = f.Get(histname)
                if not h2 or h2.IsZombie() or h2.GetEntries() < args.minTotalEntries:
                    # Fallback to plain if fiducial cut is missing
                    if args.fiducialCut:
                        histname = f"dtvqdc_{fixed_ch}_uncorrected"
                        h2 = f.Get(histname)
                    if not h2 or h2.IsZombie() or h2.GetEntries() < args.minTotalEntries:
                        f.Close()
                        continue

                # 1. Profile along QDC with standard error on mean (SEM)
                prof = h2.ProfileX(f"prof_tw_{fixed_ch}", 1, -1, "")
                nbins = prof.GetNbinsX()

                x_pts = []
                y_pts = []
                ey_pts = []

                for b_idx in range(1, nbins + 1):
                    x_val = prof.GetBinCenter(b_idx)
                    entries = prof.GetBinEntries(b_idx)
                    if x_val < args.qdcMin:
                        continue
                    if entries < args.minBinEntries:
                        continue
                    y_val = prof.GetBinContent(b_idx)
                    ey_val = prof.GetBinError(b_idx)
                    if ey_val <= 0:
                        continue

                    x_pts.append(x_val)
                    y_pts.append(y_val)
                    ey_pts.append(ey_val)

                f.Close()

                if len(x_pts) < 10:
                    n_failed += 1
                    continue

                x_arr = np.array(x_pts)
                y_arr = np.array(y_pts)
                ey_arr = np.array(ey_pts)

                qdc_min_fit = float(x_arr[0])
                qdc_max_fit = float(x_arr[-1])

                # 2. Stage 1: Unweighted Least Squares Fit (Shape Extraction)
                # Estimate t0 from high-QDC asymptotic region
                high_qdc_mask = x_arr > (0.6 * qdc_max_fit)
                if np.sum(high_qdc_mask) > 3:
                    t0_est = float(np.median(y_arr[high_qdc_mask]))
                else:
                    t0_est = float(y_arr[-1])

                p0 = [t0_est, 2.0, 1.0, 0.1, 0.005]
                bounds_lower = [t0_est - 20.0, 0.0, 0.01, 0.0, 0.0]
                bounds_upper = [t0_est + 20.0, 50.0, 10.0, 5.0, 0.5]

                try:
                    popt, _ = curve_fit(
                        tw_func,
                        x_arr,
                        y_arr,
                        p0=p0,
                        bounds=(bounds_lower, bounds_upper),
                        maxfev=5000,
                    )
                except Exception:
                    n_failed += 1
                    continue

                # 3. Stage 2: MLE / NLL Minimisation for Systematic Error Floor
                res_nll = minimize_scalar(
                    nll_objective,
                    bounds=(0.001, 0.5),
                    method="bounded",
                    args=(x_arr, y_arr, ey_arr, popt),
                )
                sigma_sys = float(res_nll.x) if res_nll.success else 0.045

                # 4. Compute Metrics: chi2 / ndf and RMS Residual
                y_fitted = tw_func(x_arr, *popt)
                var_tot = ey_arr**2 + sigma_sys**2
                ndf = max(1, len(x_arr) - 5)
                chi2 = float(np.sum(((y_arr - y_fitted) ** 2) / var_tot))
                chi2_ndf = chi2 / ndf

                # RMS residual (Thesis eq. 4.3.4)
                rms_res = float(np.sqrt(np.mean((y_arr / y_fitted - 1.0) ** 2)))

                all_sigma_sys.append(sigma_sys)
                all_chi2_ndf.append(chi2_ndf)
                all_rms_residuals.append(rms_res)

                t0, alpha, beta, qdc0, gamma = [float(v) for v in popt]
                tw_params_list = [t0, alpha, beta, qdc0, gamma]
                summary_twparams[fixed_ch] = tw_params_list

                # Save sample fits for plotting (1 per plane)
                if sipm == 0 and bar == 5:
                    sample_channels.append(fixed_ch)
                    sample_fit_data[fixed_ch] = (
                        x_arr,
                        y_arr,
                        ey_arr,
                        sigma_sys,
                        popt,
                        chi2_ndf,
                    )

                # 5. Write Individual Channel JSON (polyparams9_<fixed_ch>.json)
                ch_json_file = os.path.join(
                    out_dir, f"polyparams9_{fixed_ch}.json"
                )
                ch_data = {
                    "params": {
                        "t_0": t0,
                        "alpha": alpha,
                        "beta": beta,
                        "QDC_0": qdc0,
                        "gamma": gamma,
                    },
                    "limits": [qdc_min_fit, qdc_max_fit],
                    "chi2": chi2,
                    "ndf": ndf,
                    "alpha": sigma_sys,
                    "rms_residual": rms_res,
                }
                with open(ch_json_file, "w") as jf:
                    json.dump(ch_data, jf, indent=2)

                # 6. Write Individual CSV (polyparams9_<fixed_ch>.csv)
                csv_file = os.path.join(out_dir, f"polyparams9_{fixed_ch}.csv")
                with open(csv_file, "w", newline="") as cf:
                    writer = csv.writer(cf)
                    writer.writerow(
                        [
                            f"{t0:.5f}",
                            f"{alpha:.5f}",
                            f"{beta:.5f}",
                            f"{qdc0:.5f}",
                            f"{gamma:.5f}",
                            f"{sigma_sys:.5f}",
                            f"{chi2:.3f}",
                            f"{ndf}",
                        ]
                    )

                n_processed += 1

    # 7. Write Summary JSON files
    summary_file = os.path.join(out_dir, "twparams.json")
    with open(summary_file, "w") as jf:
        json.dump(summary_twparams, jf, indent=2)

    global_summary_file = os.path.join(out_dir, f"run{run_str}_twparams.json")
    with open(global_summary_file, "w") as jf:
        json.dump(summary_twparams, jf, indent=2)

    print("\n" + "=" * 60)
    print("Time-Walk Fit Summary:")
    print(f"  Total Channels Fitted:    {n_processed}")
    print(f"  Channels Failed / Skipped:{n_failed}")
    if all_sigma_sys:
        mean_sig = np.mean(all_sigma_sys) * 1000.0
        std_sig = np.std(all_sigma_sys) * 1000.0
        mean_chi2 = np.mean(all_chi2_ndf)
        mean_rms = np.mean(all_rms_residuals)
        print(f"  Mean Error Floor sigma_sys: {mean_sig:.1f} +/- {std_sig:.1f} ps")
        print(f"  Mean chi2 / ndf:            {mean_chi2:.2f}")
        print(f"  Mean RMS Residual:          {mean_rms:.4f} (<= 1% expected)")
    print(f"  TW Parameter dict saved:  {summary_file}")
    if not all_sigma_sys:
        print("\n[WARNING] No channels were successfully fitted. Skipping plot generation.")
        return

    # 8. Save ROOT Summary File (TH1F, TGraph, TCanvas)
    ROOT.TH1.AddDirectory(False)
    root_summary_file = os.path.join(out_dir, "twparams_summary.root")
    f_root = ROOT.TFile.Open(root_summary_file, "RECREATE")

    h_sigma = ROOT.TH1F(
        "h_sigma_sys",
        f"US Time-Walk Error Floor #sigma_{{sys}} (Run {run_str});#sigma_{{sys}} [ps];Counts",
        60,
        0.0,
        120.0,
    )
    h_chi2 = ROOT.TH1F(
        "h_chi2_ndf",
        f"US Time-Walk Fit #chi^{{2}}/#nu (Run {run_str});#chi^{{2}}/#nu;Counts",
        50,
        0.0,
        3.0,
    )
    h_rms = ROOT.TH1F(
        "h_rms_residual",
        f"US Time-Walk Fit RMS Residuals (Run {run_str});RMS Residual;Counts",
        50,
        0.0,
        0.02,
    )

    h_sigma.SetLineColor(ROOT.kBlue + 1)
    h_sigma.SetFillColorAlpha(ROOT.kBlue + 1, 0.35)
    h_chi2.SetLineColor(ROOT.kOrange + 7)
    h_chi2.SetFillColorAlpha(ROOT.kOrange + 7, 0.35)
    h_rms.SetLineColor(ROOT.kGreen + 2)
    h_rms.SetFillColorAlpha(ROOT.kGreen + 2, 0.35)

    for s_val in all_sigma_sys:
        h_sigma.Fill(s_val * 1000.0)
    for c_val in all_chi2_ndf:
        h_chi2.Fill(c_val)
    for r_val in all_rms_residuals:
        h_rms.Fill(r_val)

    h_sigma.Write()
    h_chi2.Write()
    h_rms.Write()

    # Create 4-Pad Summary TCanvas
    c_summary = ROOT.TCanvas(
        "c_tw_summary", f"US Time-Walk Summary - Run {run_str}", 1400, 1000
    )
    c_summary.Divide(2, 2, 0.01, 0.01)

    c_summary.cd(1)
    ROOT.gPad.SetGrid()
    h_sigma.Draw("HIST")

    c_summary.cd(2)
    ROOT.gPad.SetGrid()
    h_chi2.Draw("HIST")

    c_summary.cd(3)
    ROOT.gPad.SetGrid()
    h_rms.Draw("HIST")

    # Pad 4: Graph of sigma_sys vs Channel Index
    c_summary.cd(4)
    ROOT.gPad.SetGrid()
    g_sigma_ch = ROOT.TGraph()
    g_sigma_ch.SetName("g_sigma_sys_vs_channel")
    g_sigma_ch.SetTitle(
        f"US #sigma_{{sys}} vs Channel Index (Run {run_str});Channel Index;#sigma_{{sys}} [ps]"
    )
    g_sigma_ch.SetMarkerStyle(20)
    g_sigma_ch.SetMarkerSize(0.6)
    g_sigma_ch.SetMarkerColor(ROOT.kAzure + 2)

    pt_i = 0
    for plane in range(5):
        for bar in range(10):
            detID = int(f"2{plane}00{bar}")
            for sipm in all_large_sipms:
                fixed_ch = f"{detID}_{sipm}"
                if fixed_ch not in summary_twparams:
                    continue
                ch_idx = (plane * 10 + bar) * 16 + sipm
                if pt_i < len(all_sigma_sys):
                    ch_sig = all_sigma_sys[pt_i] * 1000.0
                    g_sigma_ch.SetPoint(pt_i, ch_idx, ch_sig)
                    pt_i += 1

    g_sigma_ch.Draw("AP")
    g_sigma_ch.Write()

    c_summary.Update()
    c_summary.Write()

    # Save PNG BEFORE closing TFile so objects remain in memory
    if args.plot:
        plot_png = os.path.join(plot_dir, "twparams_summary.png")
        c_summary.SaveAs(plot_png)
        print(f"  Plot PNG saved:           {plot_png}")

    f_root.Close()
    print(f"  ROOT Summary File saved:  {root_summary_file}")


if __name__ == "__main__":
    main()
