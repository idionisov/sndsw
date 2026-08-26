#!/usr/bin/env python3
"""
extract_cscint.py

Extracts the effective speed of light (c_scint / c_SiPM) in the Upstream (US)
scintillators from 2D correlation histograms (dtvxpred) produced in mode 'zeroth'
(or mode 'tw'). Evaluates the fit parameters and the 5-segment systematic
uncertainty as described in the SND@LHC timing thesis (Sections 4.1.3 & 4.3.2).

Outputs JSON and CSV files formatted for direct use by AnalysisFunctions.py.
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
        help="Generate validation and summary plots",
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

    summary_cscint = {}
    speeds_left = []
    speeds_right = []

    n_processed = 0
    n_failed = 0

    for plane in range(5):
        for bar in range(10):
            detID = int(f"2{plane}00{bar}")
            for sipm in all_large_sipms:
                fixed_ch = f"{detID}_{sipm}"

                # Skip known dead channel
                if fixed_ch == "21005_15":
                    continue

                root_file = os.path.join(root_dir, f"timewalk_{fixed_ch}.root")
                if not os.path.exists(root_file):
                    continue

                f = ROOT.TFile.Open(root_file, "READ")
                if not f or f.IsZombie():
                    continue

                histname = f"dtvxpred_{fixed_ch}_{args.state}"
                if not hasattr(f, histname):
                    f.Close()
                    continue

                h2 = f.Get(histname)
                if not h2 or h2.GetEntries() < args.minEntries:
                    f.Close()
                    continue

                # 1. Project profile along X with error on the mean
                prof = h2.ProfileX(f"prof_{fixed_ch}_{args.state}", 1, -1, "s")

                # 2. Perform Full-Range Linear Fit
                fit_func = ROOT.TF1(f"fit_{fixed_ch}", "pol1", args.fitMin, args.fitMax)
                fit_result = prof.Fit(fit_func, "QNS0R")

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

                # 6. Write/Update CSV file for Getcscint_chi2pNDF
                csv_file = os.path.join(out_dir, f"cscint_{fixed_ch}.csv")
                csv_rows = []
                if os.path.exists(csv_file):
                    with open(csv_file, "r") as cf:
                        reader = csv.reader(cf)
                        csv_rows = [row for row in reader]

                new_row = [
                    f"{c_sipm:.4f}",
                    f"{delta_c_fit:.4f}",
                    f"{const:.4f}",
                    f"{const_err:.4f}",
                    f"{chi2:.4f}",
                    f"{ndf}",
                ]

                if args.state == "uncorrected":
                    if len(csv_rows) == 0:
                        csv_rows = [new_row]
                    else:
                        csv_rows[0] = new_row
                else:  # corrected
                    if len(csv_rows) == 0:
                        csv_rows = [new_row, new_row]
                    elif len(csv_rows) == 1:
                        csv_rows.append(new_row)
                    else:
                        csv_rows[1] = new_row

                with open(csv_file, "w", newline="") as cf:
                    writer = csv.writer(cf)
                    writer.writerows(csv_rows)

                n_processed += 1

    # 7. Write Global Summary File
    summary_file = os.path.join(out_dir, f"run{run_str}_cscintvalues.json")
    with open(summary_file, "w") as jf:
        json.dump(summary_cscint, jf, indent=2)

    print("\n" + "=" * 60)
    print("Extraction Summary:")
    print(f"  Total Channels Processed: {n_processed}")
    print(f"  Channels Failed / Skipped: {n_failed}")

    all_speeds = speeds_left + speeds_right
    if all_speeds:
        mean_all = np.mean(all_speeds)
        std_all = np.std(all_speeds)
        mean_l = np.mean(speeds_left) if speeds_left else 0.0
        mean_r = np.mean(speeds_right) if speeds_right else 0.0
        print(f"  All Channels c_scint:     {mean_all:.2f} +/- {std_all:.2f} cm/ns")
        print(f"  Left-side Channels:       {mean_l:.2f} cm/ns (N={len(speeds_left)})")
        print(f"  Right-side Channels:      {mean_r:.2f} cm/ns (N={len(speeds_right)})")
    print(f"  Summary saved to:         {summary_file}")
    print("=" * 60)

    # 8. Optional Plotting
    if args.plot and all_speeds:
        try:
            import matplotlib.pyplot as plt

            # Distribution Histogram
            fig, ax = plt.subplots(figsize=(8, 6))
            bins = np.linspace(8.0, 20.0, 50)
            ax.hist(
                speeds_left,
                bins=bins,
                alpha=0.6,
                label=f"Left SiPMs (mean: {np.mean(speeds_left):.2f})",
                color="tab:blue",
            )
            ax.hist(
                speeds_right,
                bins=bins,
                alpha=0.6,
                label=f"Right SiPMs (mean: {np.mean(speeds_right):.2f})",
                color="tab:orange",
            )
            ax.axvline(
                mean_all,
                color="red",
                linestyle="--",
                label=f"Overall Mean: {mean_all:.2f} cm/ns",
            )
            ax.set_xlabel(r"$c_{\text{scint}}$ [cm/ns]", fontsize=14)
            ax.set_ylabel("Counts", fontsize=14)
            ax.set_title(
                f"US Scintillator Speed of Light ({args.state}) — Run {run_str}",
                fontsize=15,
            )
            ax.legend(fontsize=12)
            ax.grid(True, linestyle="--", alpha=0.6)
            plot_path = os.path.join(
                plot_dir, f"cscint_distribution_{args.state}.png"
            )
            fig.savefig(plot_path, dpi=200, bbox_inches="tight")
            plt.close(fig)
            print(f"Plot saved: {plot_path}")

        except ImportError:
            print("matplotlib not available; skipping PNG plot generation.")


if __name__ == "__main__":
    main()
