#!/usr/bin/env python3
"""
filter_tridents_mc.py
---------------------
High-throughput standalone filter script for extracting and storing Trident events
(mu + Z -> mu + mu+ + mu- + Z) from Monte Carlo rock muon datasets into new ROOT
files with equivalent structure.

Features:
- Pure signal event extraction (no event displays generated).
- High performance ROOT RDataFrame JIT-compiled C++ filtering.
- Selects tridents with pair vertex Z < 360.0 cm (upstream rock boundary).
- Includes both primary-muon and secondary-muon induced trident pairs.
- Preserves identical 'cbmsim' TTree structure with all branches.
- Copies detector geometry ('ShipGeo') if present.
- Supports single files or wildcard patterns with per-file output naming.
"""

import os
import sys
import glob
import time
import argparse
import ROOT

def init_root_jit(rock_z_boundary=360.0, num_threads=6):
    """
    Initializes ROOT multi-threading and declares JIT C++ filter functions.
    """
    if num_threads > 0:
        ROOT.EnableImplicitMT(num_threads)

    sndsw_path = os.environ.get("SNDSW_ROOT", "")
    if sndsw_path:
        ROOT.gInterpreter.ProcessLine(f'#include "{sndsw_path}/analysis/tools/sndSciFiTools.h"')
        ROOT.gInterpreter.ProcessLine(f'#include "{sndsw_path}/analysis/tools/sndTchainGetter.h"')
        ROOT.gInterpreter.ProcessLine(f'#include "{sndsw_path}/analysis/tools/sndGeometryGetter.h"')

    cxx_code = f"""
#ifndef FILTER_TRIDENTS_MC_H
#define FILTER_TRIDENTS_MC_H

#include "sndScifiHit.h"
#include "ShipMCTrack.h"
#include <TClonesArray.h>
#include <atomic>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <vector>
#include <cmath>

namespace TridentFilter {{
    std::atomic<long long> rdf_events{{0}};
    std::atomic<long long> total_rdf_events{{0}};
    std::chrono::steady_clock::time_point start_time;
    std::mutex cout_mutex;
    double g_rock_z_boundary = {rock_z_boundary};

    void set_rock_boundary(double z_boundary) {{
        g_rock_z_boundary = z_boundary;
    }}

    void reset_progress(long long total = 0) {{
        rdf_events = 0;
        total_rdf_events = total;
        start_time = std::chrono::steady_clock::now();
    }}

    bool print_progress() {{
        long long c = ++rdf_events;
        long long total = total_rdf_events.load();

        if (c % 100000 == 0 || (total > 0 && c == total)) {{
            auto now = std::chrono::steady_clock::now();
            auto elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();

            long long hours = elapsed_sec / 3600;
            long long minutes = (elapsed_sec % 3600) / 60;
            long long seconds = elapsed_sec % 60;

            std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << "  Processed " << c;
            if (total > 0) {{
                double pct = (100.0 * c) / total;
                std::cout << " / " << total << " (" << std::fixed << std::setprecision(1) << pct << "%)";
            }}
            std::cout << " | Elapsed: "
                      << std::setfill('0') << std::setw(2) << hours << ":"
                      << std::setfill('0') << std::setw(2) << minutes << ":"
                      << std::setfill('0') << std::setw(2) << seconds
                      << std::endl;
        }}
        return true;
    }}

    struct FilterResult {{
        bool is_signal;
        double weight;
        bool is_secondary_induced;
        double z_vertex;
    }};

    FilterResult check_trident(const TClonesArray& mc_tracks) {{
        std::vector<int> incident_muons;

        struct PairMuon {{
            int index;
            int mother_id;
            int pdg;
            double z_vert;
            double weight;
        }};
        std::vector<PairMuon> pair_muons;

        int n_tracks = mc_tracks.GetEntriesFast();
        for (int i = 0; i < n_tracks; ++i) {{
            auto* tr = static_cast<ShipMCTrack*>(mc_tracks.At(i));
            if (!tr) continue;

            int pdg = tr->GetPdgCode();
            int mother_id = tr->GetMotherId();
            int proc_id = tr->GetProcID();

            if (std::abs(pdg) == 13) {{
                incident_muons.push_back(i);
            }}
            if (std::abs(pdg) == 13 && proc_id == 5) {{
                pair_muons.push_back({{i, mother_id, pdg, tr->GetStartZ(), tr->GetWeight()}});
            }}
        }}

        for (int m_idx : incident_muons) {{
            auto* m_tr = static_cast<ShipMCTrack*>(mc_tracks.At(m_idx));
            bool is_sec = (m_tr && m_tr->GetMotherId() != -1);

            std::vector<PairMuon> daughters;
            for (const auto& d : pair_muons) {{
                if (d.mother_id == m_idx) daughters.push_back(d);
            }}

            bool has_m = false, has_p = false;
            bool zv_ok = true;
            double w = 0.0;
            double z_v = 0.0;

            for (const auto& d : daughters) {{
                if (d.pdg == 13) has_m = true;
                if (d.pdg == -13) has_p = true;
                if (d.z_vert >= g_rock_z_boundary) zv_ok = false;
                w = d.weight;
                z_v = d.z_vert;
            }}

            if (has_m && has_p && zv_ok && daughters.size() >= 2) {{
                return {{true, w, is_sec, z_v}};
            }}
        }}

        return {{false, 0.0, false, 0.0}};
    }}

    bool is_rock_trident(const TClonesArray& mc_tracks) {{
        return check_trident(mc_tracks).is_signal;
    }}

    double get_event_weight(const TClonesArray& mc_tracks) {{
        return check_trident(mc_tracks).weight;
    }}

    bool is_secondary_induced(const TClonesArray& mc_tracks) {{
        return check_trident(mc_tracks).is_secondary_induced;
    }}

    double get_vertex_z(const TClonesArray& mc_tracks) {{
        return check_trident(mc_tracks).z_vertex;
    }}
}}
#endif
    """
    ROOT.gInterpreter.Declare(cxx_code)
    ROOT.TridentFilter.set_rock_boundary(rock_z_boundary)

def extract_file_tag(filename):
    """Extracts tag such as digCPP-200 from filename."""
    base = os.path.basename(filename)
    if "digCPP-" in base:
        parts = base.split("digCPP-")
        if len(parts) > 1:
            tag = "digCPP-" + parts[1].replace(".root", "")
            return tag
    name_no_ext = os.path.splitext(base)[0]
    return name_no_ext

def process_single_file(input_file, output_file):
    """
    Processes a single ROOT file using RDataFrame, filtering trident events and
    saving 'cbmsim' and 'ShipGeo' into output_file.
    """
    t0 = time.time()

    f_test = ROOT.TFile.Open(input_file, "READ")
    if not f_test or f_test.IsZombie():
        print(f"Error: Could not open '{input_file}'")
        return 0, 0.0, 0

    tree = f_test.Get("cbmsim")
    if not tree:
        print(f"Error: 'cbmsim' tree not found in '{input_file}'")
        f_test.Close()
        return 0, 0.0, 0

    total_events = tree.GetEntries()
    geo_obj = f_test.Get("ShipGeo")
    f_test.Close()

    ROOT.TridentFilter.reset_progress(total_events)

    df = ROOT.RDataFrame("cbmsim", input_file)
    df_filtered = (
        df.Filter("TridentFilter::print_progress()")
          .Filter("TridentFilter::is_rock_trident(MCTrack)")
          .Define("weight", "TridentFilter::get_event_weight(MCTrack)")
          .Define("is_sec_induced", "TridentFilter::is_secondary_induced(MCTrack)")
          .Define("z_vertex", "TridentFilter::get_vertex_z(MCTrack)")
    )

    opts = ROOT.RDF.RSnapshotOptions()
    opts.fMode = "RECREATE"
    df_filtered.Snapshot("cbmsim", output_file, "", opts)

    # Re-open input and output files to check counts and copy all metadata objects
    f_in = ROOT.TFile.Open(input_file, "READ")
    metadata_keys = []
    if f_in and not f_in.IsZombie():
        for k in f_in.GetListOfKeys():
            kname = k.GetName()
            kclass = k.GetClassName()
            if kclass != "TTree" and kname != "cbmsim":
                obj = f_in.Get(kname)
                if obj:
                    metadata_keys.append((kname, obj))

    f_out = ROOT.TFile.Open(output_file, "UPDATE")
    out_tree = f_out.Get("cbmsim")
    selected_count = out_tree.GetEntries() if out_tree else 0

    total_weighted_yield = 0.0
    sec_induced_count = 0

    if selected_count > 0:
        for ev in out_tree:
            total_weighted_yield += ev.weight
            if hasattr(ev, "is_sec_induced") and ev.is_sec_induced:
                sec_induced_count += 1

    f_out.cd()
    for kname, obj in metadata_keys:
        obj.Write(kname, ROOT.TObject.kSingleKey | ROOT.TObject.kOverwrite)

    f_out.Close()
    if f_in:
        f_in.Close()

    elapsed = time.time() - t0
    rate = total_events / elapsed if elapsed > 0 else 0
    print(f"  Processed {total_events:,} events in {elapsed:.1f}s ({rate:.0f} ev/s)")
    print(f"  Selected: {selected_count:,} events (Yield: {total_weighted_yield:.6f}) | Secondary-Induced: {sec_induced_count}")

    return total_events, selected_count, total_weighted_yield, sec_induced_count

def main():
    default_input = "/eos/experiment/sndlhc/MonteCarlo/ThreeMuons/sndLHC.Ntuple-TGeant4_boost100LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_filteredAtScoringPlane_digCPP-2*.root"
    default_output = "trimuon_filtered_%s.root"

    parser = argparse.ArgumentParser(
        description="Extract and save rock trident signal events into ROOT files of equivalent structure."
    )
    parser.add_argument(
        "-i", "--input",
        dest="input_pattern",
        default=default_input,
        help="Input ROOT file path or wildcard pattern (default: %(default)s)"
    )
    parser.add_argument(
        "-o", "--output",
        dest="output_pattern",
        default=default_output,
        help="Output ROOT file path or pattern with %%s for tag (default: %(default)s)"
    )
    parser.add_argument(
        "-z", "--rock_boundary",
        dest="rock_boundary",
        type=float,
        default=360.0,
        help="Rock volume Z boundary in cm (default: %(default)s cm)"
    )
    parser.add_argument(
        "-j", "--threads",
        dest="num_threads",
        type=int,
        default=6,
        help="Number of worker threads for RDataFrame (default: %(default)s)"
    )

    args = parser.parse_args()

    # Find matching files
    matched_files = sorted(glob.glob(args.input_pattern))
    if not matched_files:
        if os.path.exists(args.input_pattern):
            matched_files = [args.input_pattern]
        else:
            print(f"Error: No files found matching pattern: {args.input_pattern}")
            sys.exit(1)

    print("=" * 60)
    print("SND@LHC MC Trident Signal Event Extractor")
    print(f"Input Pattern     : {args.input_pattern}")
    print(f"Files Found       : {len(matched_files)}")
    print(f"Rock Z Boundary   : Z < {args.rock_boundary:.1f} cm")
    print(f"Worker Threads    : {args.num_threads}")
    print("=" * 60)

    init_root_jit(rock_z_boundary=args.rock_boundary, num_threads=args.num_threads)

    overall_t0 = time.time()
    grand_total_events = 0
    grand_total_selected = 0
    grand_total_yield = 0.0
    grand_total_sec = 0
    created_files = []

    for idx, input_file in enumerate(matched_files, 1):
        tag = extract_file_tag(input_file)
        if "%s" in args.output_pattern:
            out_file = args.output_pattern % tag
        elif len(matched_files) > 1:
            base_name, ext = os.path.splitext(args.output_pattern)
            out_file = f"{base_name}_{tag}{ext}"
        else:
            out_file = args.output_pattern

        print(f"\n[{idx}/{len(matched_files)}] Filtering: '{input_file}' -> '{out_file}'")
        t_evts, sel_evts, w_yield, sec_evts = process_single_file(input_file, out_file)

        grand_total_events += t_evts
        grand_total_selected += sel_evts
        grand_total_yield += w_yield
        grand_total_sec += sec_evts
        created_files.append(out_file)

    overall_elapsed = time.time() - overall_t0

    print("\n" + "=" * 60)
    print("FILTERING COMPLETE SUMMARY")
    print("=" * 60)
    print(f"Total Input Files Processed : {len(matched_files)}")
    print(f"Total Output Files Created   : {len(created_files)}")
    print(f"Total Events Scanned        : {grand_total_events:,}")
    print(f"Total Signal Events Stored  : {grand_total_selected:,}")
    print(f"Total Weighted Yield        : {grand_total_yield:.6f}")
    print(f"Secondary-Induced Tridents  : {grand_total_sec:,}")
    print(f"Total Processing Time       : {overall_elapsed:.1f} s")
    print("=" * 60)

if __name__ == "__main__":
    main()
