#!/usr/bin/env python3
"""
extract_cscint.py

Extracts the effective speed of light (c_scint / c_SiPM) in the Upstream (US)
scintillators from 2D correlation histograms (dtvxpred) produced in mode 'zeroth'
(or mode 'tw'). Evaluates the fit parameters and the 5-segment systematic
uncertainty (Sections 4.1.3 & 4.3.2).

Supports exporting per-channel 2-pad validation canvases into dedicated TDirectory
subdirectories in the output ROOT summary file (--all / -c all).

Outputs JSON, CSV, and consolidated ROOT files formatted for direct use by AnalysisFunctions.py.
"""

import os
import sys
import json
import csv
import argparse
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gStyle.SetOptStat(0)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract effective speed of light in US scintillators."
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
        "-s",
        "--state",
        type=str,
        default="uncorrected",
        choices=["uncorrected", "corrected"],
        help="State of the histograms to analyze: 'uncorrected' (zeroth) or 'corrected' (tw)",
    )
    parser.add_argument(
        "-c",
        "--channel",
        type=str,
        default="all",
        help="Target channel (e.g. '24005_4') or 'all' for all detector channels (default: all)",
    )
    parser.add_argument(
        "--all",
        dest="all_channels",
        action="store_true",
        default=True,
        help="Process and export all detector channels with per-channel TDirectories (default: True)",
    )
    parser.add_argument(
        "--fitMin",
        type=float,
        default=-65.0,
        help="Lower bound in x_pred [cm] for the linear fit (default: -65.0)",
    )
    parser.add_argument(
        "--fitMax",
        type=float,
        default=-5.0,
        help="Upper bound in x_pred [cm] for the linear fit (default: -5.0)",
    )
    parser.add_argument(
        "--minEntries",
        type=int,
        default=200,
        help="Minimum histogram entries required to perform the fit (default: 200)",
    )
    parser.add_argument(
        "--nSegments",
        type=int,
        default=5,
        help="Number of sub-segments for systematic uncertainty calculation (default: 5)",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate validation and summary PNG plots",
    )
    parser.add_argument(
        "--plotDir",
        type=str,
        default="",
        help="Directory to save generated plots (default: <path>/plots/run<runNr>/)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    run_str = f"{args.runNumber:06d}"
    base_path = args.path.rstrip("/") + "/"
    root_dir = os.path.join(base_path, f"rootfiles/run{run_str}")
    out_dir = os.path.join(base_path, f"cscintvalues/run{run_str}")
    os.makedirs(out_dir, exist_ok=True)

    if args.plotDir:
        plot_dir = args.plotDir
    else:
        plot_dir = os.path.join(base_path, f"plots/run{run_str}/cscint")
    if args.plot:
        os.makedirs(plot_dir, exist_ok=True)

    print("=" * 60)
    print("SND@LHC US Scintillator Effective Signal Speed Extraction")
    print(f"  Run Number:      {run_str}")
    print(f"  State:           {args.state}")
    print(f"  Target Channel:  {args.channel}")
    print(f"  ROOT Files Dir:  {root_dir}")
    print(f"  Output JSON Dir: {out_dir}")
    print(f"  Fit Range [x]:   [{args.fitMin}, {args.fitMax}] cm")
    print(f"  Sub-segments:    {args.nSegments}")
    print("=" * 60)

    if not os.path.exists(root_dir):
        print(f"Error: Directory {root_dir} does not exist!")
        sys.exit(1)

    # Large SiPM channel indices (exclude small SiPMs: 2, 5, 10, 13)
    left_sipms = [0, 1, 3, 4, 6, 7]
    right_sipms = [8, 9, 11, 12, 14, 15]
    all_large_sipms = left_sipms + right_sipms

    # Determine channel list
    if args.channel != "all":
        channels_to_process = [args.channel]
    else:
        channels_to_process = []
        for plane in range(5):
            for bar in range(10):
                detID = int(f"2{plane}00{bar}")
                for sipm in all_large_sipms:
                    fixed_ch = f"{detID}_{sipm}"
                    if fixed_ch != "21005_15":  # Skip known dead channel
                        channels_to_process.append(fixed_ch)

    summary_cscint = {}
    speeds_left = []
    speeds_right = []

    n_processed = 0
    n_failed = 0

    # Output ROOT summary file
    root_summary_file = os.path.join(out_dir, f"cscint_summary_{args.state}.root")
    f_root = ROOT.TFile.Open(root_summary_file, "RECREATE")

    for fixed_ch in channels_to_process:
        parts = fixed_ch.split("_")
        detID = int(parts[0])
        sipm = int(parts[1])
        plane = (detID // 1000) % 10
        bar = detID % 10

        root_file = os.path.join(root_dir, f"timewalk_{fixed_ch}.root")
        if not os.path.exists(root_file):
            continue

        try:
            f = ROOT.TFile.Open(root_file, "READ")
        except Exception:
            n_failed += 1
            continue

        if not f or f.IsZombie():
            n_failed += 1
            continue

        histname = f"dtvxpred_{fixed_ch}_{args.state}"
        h2_src = f.Get(histname)
        if not h2_src or h2_src.IsZombie() or h2_src.GetEntries() < args.minEntries:
            f.Close()
            continue

        # Detach histogram from input file
        h2 = h2_src.Clone(f"h2_{histname}")
        h2.SetDirectory(0)

        # 1. Project profile along X with error on the mean
        prof = h2.ProfileX(f"prof_{fixed_ch}_{args.state}", 1, -1, "s")
        prof.SetDirectory(0)

        # 2. Perform Full-Range Linear Fit
        fit_func = ROOT.TF1(f"fit_{fixed_ch}", "pol1", args.fitMin, args.fitMax)
        prof.Fit(fit_func, "QNS0R")

        slope = fit_func.GetParameter(1)
        slope_err = fit_func.GetParError(1)
        const = fit_func.GetParameter(0)
        const_err = fit_func.GetParError(0)
        chi2 = fit_func.GetChisquare()
        ndf = max(1, fit_func.GetNDF())

        if abs(slope) < 1e-7 or slope_err <= 0:
            f.Close()
            n_failed += 1
            continue

        # Sign convention: +1/m for left (signal -> +x), -1/m for right (signal -> -x)
        side_sign = 1.0 if sipm < 8 else -1.0
        c_sipm = side_sign / slope
        delta_c_fit = abs(slope_err / (slope**2))

        # Discard unphysical fit results
        if c_sipm <= 0 or c_sipm > 50.0:
            f.Close()
            n_failed += 1
            continue

        # 3. Perform 5-Segment Sub-Fits for Systematic Uncertainty
        seg_bounds = np.linspace(args.fitMin, args.fitMax, args.nSegments + 1)
        subrange_speeds = []

        for s_idx in range(args.nSegments):
            seg_func = ROOT.TF1(
                f"seg_{fixed_ch}_{s_idx}",
                "pol1",
                seg_bounds[s_idx],
                seg_bounds[s_idx + 1],
            )
            prof.Fit(seg_func, "QNS0R")
            m_s = seg_func.GetParameter(1)
            if abs(m_s) > 1e-7:
                c_s = side_sign / m_s
                if 5.0 < c_s < 40.0:
                    subrange_speeds.append(float(c_s))

        if len(subrange_speeds) >= 2:
            delta_c_syst = float(np.std(subrange_speeds, ddof=1))
        else:
            delta_c_syst = float(delta_c_fit)

        f.Close()

        # Record statistics
        if sipm < 8:
            speeds_left.append(c_sipm)
        else:
            speeds_right.append(c_sipm)

        summary_cscint[fixed_ch] = [
            float(c_sipm),
            float(delta_c_fit),
            float(delta_c_syst),
        ]

        # 4. Write Individual Channel JSON
        ch_json_file = os.path.join(out_dir, f"cscint_{fixed_ch}.json")
        ch_json_data = {}
        if os.path.exists(ch_json_file):
            try:
                with open(ch_json_file, "r") as jf:
                    ch_json_data = json.load(jf)
            except Exception:
                ch_json_data = {}

        ch_json_data[args.state] = [
            float(c_sipm),
            float(delta_c_fit),
            float(const),
            float(const_err),
            float(chi2),
            int(ndf),
            float(delta_c_syst),
        ]

        with open(ch_json_file, "w") as jf:
            json.dump(ch_json_data, jf, indent=2)

        # 5. Write Subranges JSON
        subrange_json_file = os.path.join(
            out_dir, f"subrange-cscints_{fixed_ch}.json"
        )
        subrange_json_data = {}
        if os.path.exists(subrange_json_file):
            try:
                with open(subrange_json_file, "r") as jf:
                    subrange_json_data = json.load(jf)
            except Exception:
                subrange_json_data = {}

        subrange_json_data[args.state] = subrange_speeds
        with open(subrange_json_file, "w") as jf:
            json.dump(subrange_json_data, jf, indent=2)

        # 6. Save Per-Channel TDirectory in ROOT summary file
        ch_dir = f_root.mkdir(fixed_ch)
        ch_dir.cd()

        # Build 2-pad Canvas for channel fit
        side_str = "Left" if sipm < 8 else "Right"
        c_ch = ROOT.TCanvas(f"c_cscint_fit_{fixed_ch}", f"c_scint Fit {fixed_ch}", 900, 800)
        c_ch.Divide(1, 2)
        
        # Upper Pad: 2D Hist + Fit Line + Profile
        pad1 = c_ch.cd(1)
        pad1.SetPad(0.0, 0.30, 1.0, 1.0)
        pad1.SetBottomMargin(0.02)
        pad1.SetGrid()
        
        h2.SetTitle(f"US Plane {plane}, Bar {bar}, SiPM {sipm} ({side_str});;t_{{0}}^{{DS}} - t_{{SiPM}} [ns]")
        h2.Draw("COLZ")
        
        prof.SetMarkerStyle(20)
        prof.SetMarkerSize(0.6)
        prof.SetMarkerColor(ROOT.kBlack)
        prof.SetLineColor(ROOT.kBlack)
        prof.Draw("P SAME")
        
        fit_func.SetLineColor(ROOT.kRed + 1)
        fit_func.SetLineWidth(2)
        fit_func.Draw("SAME")
        
        leg = ROOT.TLegend(0.48, 0.68, 0.88, 0.88)
        leg.SetBorderSize(1)
        leg.SetFillColor(ROOT.kWhite)
        leg.AddEntry(prof, "Profile Data (#pm SEM)", "ep")
        leg.AddEntry(fit_func, f"Linear Fit: c = {c_sipm:.2f} #pm {delta_c_fit:.2f} #pm {delta_c_syst:.2f} cm/ns", "l")
        leg.AddEntry("", f"#chi^{{2}}/#nu = {chi2/ndf:.2f}, Offset = {const:.2f} #pm {const_err:.2f} ns", "")
        leg.Draw()

        # Lower Pad: Residuals
        pad2 = c_ch.cd(2)
        pad2.SetPad(0.0, 0.0, 1.0, 0.30)
        pad2.SetTopMargin(0.03)
        pad2.SetBottomMargin(0.28)
        pad2.SetGrid()

        g_res = ROOT.TGraphErrors()
        g_res.SetName(f"g_res_{fixed_ch}")
        g_res.SetTitle(f";x_{{predicted}} [cm];Data - Fit [ns]")
        pt_idx = 0
        for b_idx in range(1, prof.GetNbinsX() + 1):
            x_val = prof.GetBinCenter(b_idx)
            if args.fitMin <= x_val <= args.fitMax:
                y_val = prof.GetBinContent(b_idx)
                y_err = prof.GetBinError(b_idx)
                if y_val != 0 or y_err != 0:
                    res_val = y_val - fit_func.Eval(x_val)
                    g_res.SetPoint(pt_idx, x_val, res_val)
                    g_res.SetPointError(pt_idx, 0, y_err)
                    pt_idx += 1

        g_res.SetMarkerStyle(20)
        g_res.SetMarkerSize(0.6)
        g_res.SetMarkerColor(ROOT.kBlack)
        g_res.Draw("AP")
        g_res.GetXaxis().SetLimits(args.fitMin - 5.0, args.fitMax + 5.0)
        g_res.GetXaxis().SetTitleSize(0.10)
        g_res.GetXaxis().SetLabelSize(0.09)
        g_res.GetYaxis().SetTitleSize(0.10)
        g_res.GetYaxis().SetLabelSize(0.09)
        g_res.GetYaxis().SetTitleOffset(0.45)

        line0 = ROOT.TLine(args.fitMin - 5.0, 0, args.fitMax + 5.0, 0)
        line0.SetLineColor(ROOT.kRed + 1)
        line0.SetLineStyle(2)
        line0.Draw("SAME")

        c_ch.Write()
        h2.Write(f"h2_{histname}")
        prof.Write(f"prof_{fixed_ch}")
        fit_func.Write(f"fit_{fixed_ch}")

        if args.plot and args.channel != "all":
            ch_png = os.path.join(plot_dir, f"cscint_fit_{fixed_ch}_{args.state}.png")
            c_ch.SaveAs(ch_png)

        f_root.cd()
        n_processed += 1

    # 7. Write Consolidated JSON
    summary_file = os.path.join(
        out_dir, f"run{run_str}_cscintvalues_{args.state}.json"
    )
    with open(summary_file, "w") as sf:
        json.dump(summary_cschist := summary_cscint, sf, indent=2)

    # Legacy copy for backward compatibility
    legacy_file = os.path.join(out_dir, f"run{run_str}_cscintvalues.json")
    with open(legacy_file, "w") as lf:
        json.dump(summary_cscint, lf, indent=2)

    print("\n" + "=" * 60)
    print("Extraction Summary:")
    print(f"  Total Channels Processed: {n_processed}")
    all_speeds = speeds_left + speeds_right
    if not all_speeds:
        print("\n[WARNING] No channels were successfully processed.")
        print(f"Check that ROOT files in {root_dir} contain histogram 'dtvxpred_<channel>_{args.state}'.")
        f_root.Close()
        return

    mean_all = np.mean(all_speeds)
    std_all = np.std(all_speeds)
    mean_l = np.mean(speeds_left) if speeds_left else 0.0
    mean_r = np.mean(speeds_right) if speeds_right else 0.0
    print(f"  All Channels c_scint:     {mean_all:.2f} +/- {std_all:.2f} cm/ns")
    print(f"  Left-side Channels:       {mean_l:.2f} cm/ns (N={len(speeds_left)})")
    print(f"  Right-side Channels:      {mean_r:.2f} cm/ns (N={len(speeds_right)})")
    print(f"  Summary saved to:         {summary_file}")
    print("=" * 60)

    # 8. Save Global ROOT Histograms, TGraphErrors, and TCanvas at File Root
    f_root.cd()

    h_all = ROOT.TH1F(
        "h_cscint_all",
        f"US c_{{scint}} ({args.state}) - All Channels;c_{{scint}} [cm/ns];Counts",
        60,
        8.0,
        20.0,
    )
    h_left = ROOT.TH1F(
        "h_cscint_left",
        f"US c_{{scint}} ({args.state}) - Left SiPMs;c_{{scint}} [cm/ns];Counts",
        60,
        8.0,
        20.0,
    )
    h_right = ROOT.TH1F(
        "h_cscint_right",
        f"US c_{{scint}} ({args.state}) - Right SiPMs;c_{{scint}} [cm/ns];Counts",
        60,
        8.0,
        20.0,
    )

    h_all.SetLineColor(ROOT.kBlack)
    h_all.SetLineWidth(2)
    h_left.SetLineColor(ROOT.kBlue + 1)
    h_left.SetFillColorAlpha(ROOT.kBlue + 1, 0.35)
    h_left.SetLineWidth(2)
    h_right.SetLineColor(ROOT.kOrange + 7)
    h_right.SetFillColorAlpha(ROOT.kOrange + 7, 0.35)
    h_right.SetLineWidth(2)

    for v in speeds_left:
        h_left.Fill(v)
        h_all.Fill(v)
    for v in speeds_right:
        h_right.Fill(v)
        h_all.Fill(v)

    h_all.Write()
    h_left.Write()
    h_right.Write()

    # Graphs vs Channel Index
    g_left = ROOT.TGraphErrors()
    g_left.SetName("g_cscint_left")
    g_left.SetTitle(f"Left SiPMs c_{{scint}} vs Channel;Channel Index;c_{{scint}} [cm/ns]")
    g_left.SetMarkerStyle(20)
    g_left.SetMarkerSize(0.7)
    g_left.SetMarkerColor(ROOT.kBlue + 1)
    g_left.SetLineColor(ROOT.kBlue + 1)

    g_right = ROOT.TGraphErrors()
    g_right.SetName("g_cscint_right")
    g_right.SetTitle(f"Right SiPMs c_{{scint}} vs Channel;Channel Index;c_{{scint}} [cm/ns]")
    g_right.SetMarkerStyle(21)
    g_right.SetMarkerSize(0.7)
    g_right.SetMarkerColor(ROOT.kOrange + 7)
    g_right.SetLineColor(ROOT.kOrange + 7)

    pt_l = 0
    pt_r = 0
    for plane in range(5):
        for bar in range(10):
            detID = int(f"2{plane}00{bar}")
            for sipm in all_large_sipms:
                fixed_ch = f"{detID}_{sipm}"
                if fixed_ch not in summary_cscint:
                    continue
                ch_idx = (plane * 10 + bar) * 16 + sipm
                c_val, d_fit, d_syst = summary_cscint[fixed_ch]
                d_total = np.sqrt(d_fit**2 + d_syst**2)

                if sipm < 8:
                    g_left.SetPoint(pt_l, ch_idx, c_val)
                    g_left.SetPointError(pt_l, 0, d_total)
                    pt_l += 1
                else:
                    g_right.SetPoint(pt_r, ch_idx, c_val)
                    g_right.SetPointError(pt_r, 0, d_total)
                    pt_r += 1

    g_left.Write()
    g_right.Write()

    # Create Summary TCanvas with 2 Pads (Distribution + Channel Scatter)
    c_summary = ROOT.TCanvas(
        "c_cscint_summary",
        f"US Scintillator Speed of Light Summary - Run {run_str}",
        1200,
        900,
    )
    c_summary.Divide(1, 2, 0.01, 0.01)

    # Pad 1: Overlaid Distributions
    c_summary.cd(1)
    ROOT.gPad.SetGrid()
    max_y = max(h_left.GetMaximum(), h_right.GetMaximum(), h_all.GetMaximum()) * 1.25
    h_all.SetMaximum(max_y)
    h_all.Draw("HIST")
    h_left.Draw("HIST SAME")
    h_right.Draw("HIST SAME")

    leg1 = ROOT.TLegend(0.62, 0.65, 0.88, 0.88)
    leg1.SetBorderSize(1)
    leg1.SetFillColor(ROOT.kWhite)
    leg1.AddEntry(h_all, f"All: {mean_all:.2f} #pm {std_all:.2f} cm/ns", "l")
    leg1.AddEntry(h_left, f"Left: {mean_l:.2f} cm/ns (N={len(speeds_left)})", "f")
    leg1.AddEntry(h_right, f"Right: {mean_r:.2f} cm/ns (N={len(speeds_right)})", "f")
    leg1.Draw()

    # Pad 2: Scatter vs Channel
    c_summary.cd(2)
    ROOT.gPad.SetGrid()
    mg = ROOT.TMultiGraph()
    mg.SetTitle(f"US Scintillator c_{{scint}} as a Function of SiPM Channel (Run {run_str});SiPM Channel Index;c_{{scint}} [cm/ns]")
    mg.Add(g_left, "P")
    mg.Add(g_right, "P")
    mg.Draw("A")
    mg.GetYaxis().SetRangeUser(8.0, 22.0)

    leg2 = ROOT.TLegend(0.72, 0.72, 0.88, 0.88)
    leg2.SetBorderSize(1)
    leg2.SetFillColor(ROOT.kWhite)
    leg2.AddEntry(g_left, "Left SiPMs", "p")
    leg2.AddEntry(g_right, "Right SiPMs", "p")
    leg2.Draw()

    c_summary.Write()
    f_root.Close()
    print(f"  ROOT Summary File saved:  {root_summary_file}")

    # Also export TCanvas to image file if requested
    if args.plot:
        plot_png = os.path.join(plot_dir, f"cscint_distribution_{args.state}.png")
        c_summary.SaveAs(plot_png)
        print(f"  Plot PNG saved:           {plot_png}")


if __name__ == "__main__":
    main()
