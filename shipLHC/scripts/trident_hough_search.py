import os
import random
import ROOT
import argparse
import csv
import datetime
import time
from array import array
import numpy as np
import pandas as pd
import rootUtils as ut

ROOT.gROOT.SetBatch(True)
for library in ["libBase", "libShipData", "libshipLHC", "libsnd_analysis_tools"]:
    ROOT.gSystem.Load(library)

import SndlhcGeo
import SndlhcMuonReco
import SndlhcTracking

sndsw_path = os.environ['SNDSW_ROOT']
ROOT.gInterpreter.ProcessLine(f'#include "{sndsw_path}/analysis/tools/sndSciFiTools.h"')
ROOT.gInterpreter.ProcessLine(f'#include "{sndsw_path}/analysis/tools/sndGeometryGetter.h"')

histograms = {}
projections = {1: 'xz', 2: 'yz'}

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file',  type=str, required=True, help='Path to input file.')
    parser.add_argument('-csv', '--output-csv', type=str, required=True, help='Output CSV filename')
    parser.add_argument('-o', '--output-root', type=str, help='Output ROOT filename for selected events')
    parser.add_argument('-s', '--start-event', type=int, help='Start event number.', default=0)
    parser.add_argument('-n', '--n-events', type=int, help='Number of events.', default=1000000)
    parser.add_argument('-f', '--fraction', type=float, default=1.0, help='Fraction of events to process')
    parser.add_argument('-par', '--parFile', type=str, default=os.path.join(sndsw_path, "python/TrackingParams.xml"), help='Tracking parameter file')
    parser.add_argument('-zmin', '--z-vtx-min', type=float, default=None, help='Minimum z-vertex intersection')
    parser.add_argument('-zmax', '--z-vtx-max', type=float, default=None, help='Maximum z-vertex intersection')
    return parser.parse_args()

def get_scifi_hit_density(scifi_hits, radius=40, min_check=False, min_hit_density=1000000):
    density = 0
    for hit in scifi_hits:
        channel_id = hit.GetChannelID()
        density += ROOT.snd.analysis_tools.densityScifi(channel_id, scifi_hits, radius, min_hit_density, min_check)
    return density

def get_scifi_total_qdc(scifi_hits):
    total_qdc = 0
    for hit in scifi_hits:
        total_qdc += hit.GetSignal()
    return total_qdc

def initialize_event_display(geo):
    """Setup the canvas and background histograms for event display."""
    ut.bookCanvas(histograms, key='simpleDisplay', title='simple event display', nx=1200, ny=1016, cx=1, cy=2)

    # Determine Z range based on detector floor position
    z_start = 250. if geo.snd_geo.Floor.z != 0 else 60.
    z_end = z_start + 350.
    x_start, y_start = -100., -30.

    histograms.update({
        'xmin': x_start, 'xmax': x_start + 110.,
        'ymin': y_start, 'ymax': y_start + 110.,
        'zmin': z_start, 'zmax': z_end
    })

    for projection, title in [('xz', '; z [cm]; x [cm]'), ('yz', '; z [cm]; y [cm]')]:
        ut.bookHist(histograms, projection, title, 500, histograms['zmin'], histograms['zmax'], 100, histograms[f'{projection[0]}min'], histograms[f'{projection[0]}max'])
        histograms[projection].SetStats(0)

def draw_detector_geometry(geo):
    """Draw geometry boxes for display."""
    mu_filter = geo.snd_geo.MuFilter
    scifi = geo.snd_geo.Scifi
    emulsion = geo.snd_geo.EmulsionDet
    navigator = ROOT.gGeoManager.GetCurrentNavigator()

    nodes = {'volMuFilter_1/volFeBlockEnd_1': ROOT.kGreen-6}

    for i in range(mu_filter.NVetoPlanes):
        nodes[f'volVeto_1/volVetoPlane_{i}_{i}'] = ROOT.kRed
        for j in range(mu_filter.NVetoBars):
            bar_suffix = 'Bar_1' if i < 2 else 'Bar_ver_1'
            nodes[f'volVeto_1/volVetoPlane_{i}_{i}/volVeto{bar_suffix}{i}{j:0>3d}'] = ROOT.kRed
        box_type = "3" if i == 2 else ""
        nodes[f'volVeto_1/subVeto{box_type}Box_{i}'] = ROOT.kGray+1

    for i in range(scifi.nscifi):
        nodes[f'volTarget_1/ScifiVolume{i+1}_{i+1}000000'] = ROOT.kBlue+1
        nodes[f'volTarget_1/volFeTarget{i+1}_1'] = ROOT.kGreen-6
    for i in range(emulsion.wall):
        nodes[f'volTarget_1/volWallborder_{i}'] = ROOT.kGray

    ds_offset = mu_filter.NVetoPlanes + mu_filter.NUpstreamPlanes
    for i in range(mu_filter.NDownstreamPlanes):
        ds_plane_name = f'volMuFilter_1/volMuDownstreamDet_{i}_{i+ds_offset}'
        nodes[ds_plane_name] = ROOT.kBlue+1
        for j in range(mu_filter.NDownstreamBars):
            nodes[f'{ds_plane_name}/volMuDownstreamBar_ver_3{i}{j+mu_filter.NDownstreamBars:0>3d}'] = ROOT.kBlue+1
            if i < 3:
                nodes[f'{ds_plane_name}/volMuDownstreamBar_hor_3{i}{j:0>3d}'] = ROOT.kBlue+1
        nodes[f'volMuFilter_1/subDSBox_{i+ds_offset}'] = ROOT.kGray+1

    for i in range(mu_filter.NUpstreamPlanes):
        nodes[f'volMuFilter_1/subUSBox_{i+mu_filter.NVetoPlanes}'] = ROOT.kGray+1
        nodes[f'volMuFilter_1/volMuUpstreamDet_{i}_{i+mu_filter.NVetoPlanes}'] = ROOT.kBlue+1
        for j in range(mu_filter.NUpstreamBars):
            nodes[f'volMuFilter_1/volMuUpstreamDet_{i}_{i+mu_filter.NVetoPlanes}/volMuUpstreamBar_2{i}00{j}'] = ROOT.kBlue+1
        nodes[f'volMuFilter_1/volFeBlock_{i}'] = ROOT.kGreen-6

    for i in range(ds_offset, ds_offset + mu_filter.NDownstreamPlanes):
        nodes[f'volMuFilter_1/volFeBlock_{i}'] = ROOT.kGreen-6

    pass_nodes = {'Block', 'Wall', 'FeTarget'}
    transverse_nodes = {'UpstreamBar', 'VetoBar', 'hor'}
    axis_map = {'X': 0, 'Y': 1}

    for node_path, color in nodes.items():
        full_path = f'/cave_1/Detector_0/{node_path}'
        for view in ['X', 'Y']:
            if not navigator.CheckPath(full_path):
                continue

            navigator.cd(full_path)
            shape = navigator.GetCurrentNode().GetVolume().GetShape()
            dx, dy, dz = shape.GetDX(), shape.GetDY(), shape.GetDZ()
            ox, oy, oz = shape.GetOrigin()

            corners = {}
            if view == 'X' and (not any(tn in full_path for tn in transverse_nodes) or 'VetoBar_ver' in full_path):
                corners = {'LB': [-dx+ox, oy, -dz+oz], 'LT': [dx+ox, oy, -dz+oz], 'RB': [-dx+ox, oy, dz+oz], 'RT': [dx+ox, oy, dz+oz]}
            elif view == 'Y' and 'ver' not in full_path:
                corners = {'LB': [ox, -dy+oy, -dz+oz], 'LT': [ox, dy+oy, -dz+oz], 'RB': [ox, -dy+oy, dz+oz], 'RT': [ox, dy+oy, dz+oz]}

            if not corners: continue

            polyline = ROOT.TPolyLine()
            projection_axis = axis_map[view]
            for i, corner_key in enumerate(['LB', 'LT', 'RT', 'RB', 'LB']):
                master_coords = array('d', [0, 0, 0])
                navigator.LocalToMaster(array('d', corners[corner_key]), master_coords)
                polyline.SetPoint(i, master_coords[2], master_coords[projection_axis])

            polyline.SetLineColor(color)
            polyline.SetLineWidth(1)
            histograms['simpleDisplay'].cd(projection_axis + 1)

            if any(pn in full_path for pn in pass_nodes):
                polyline.SetFillColorAlpha(color, 0.5)
                polyline.Draw('f&&same')
            polyline.Draw('same')
            histograms[f"{full_path}_{view}"] = polyline

def draw_event_hits_and_tracks(event, geo, hough_lines):
    """Draw hits and hough lines on the canvas."""
    for p_idx in projections:
        histograms['simpleDisplay'].cd(p_idx)
        histograms[projections[p_idx]].Draw('b')

    draw_detector_geometry(geo)

    hit_graphs = {view: {system: ROOT.TGraphErrors() for system in ['Scifi', 'MuFilter']} for view in ['X', 'Y']}
    counts = {view: {system: 0 for system in ['Scifi', 'MuFilter']} for view in ['X', 'Y']}

    scifi_geo = geo.snd_geo.Scifi
    for hit in event.Digi_ScifiHits:
        if not hit.isValid(): continue
        pos_a, pos_b = ROOT.TVector3(), ROOT.TVector3()
        geo.modules['Scifi'].GetSiPMPosition(hit.GetDetectorID(), pos_a, pos_b)
        view = 'X' if hit.isVertical() else 'Y'
        graph, idx = hit_graphs[view]['Scifi'], counts[view]['Scifi']
        graph.SetPoint(idx, pos_a.Z(), pos_a.X() if view == 'X' else pos_a.Y())
        graph.SetPointError(idx, scifi_geo.scifimat_z/2, scifi_geo.channel_width/2)
        counts[view]['Scifi'] += 1

    mu_filter_geo = geo.snd_geo.MuFilter
    for hit in event.Digi_MuFilterHits:
        if not hit.isValid(): continue
        pos_a, pos_b = ROOT.TVector3(), ROOT.TVector3()
        geo.modules['MuFilter'].GetPosition(hit.GetDetectorID(), pos_a, pos_b)
        view = 'X' if hit.isVertical() else 'Y'
        graph, idx = hit_graphs[view]['MuFilter'], counts[view]['MuFilter']
        system = hit.GetSystem()

        if system == 1:
            dx, dy, dz = (mu_filter_geo.Veto3BarX/2 if hasattr(mu_filter_geo, "Veto3BarX") else mu_filter_geo.VetoBarX/2), mu_filter_geo.VetoBarY/2, mu_filter_geo.VetoBarZ/2
        elif system == 2:
            dx, dy, dz = mu_filter_geo.UpstreamBarX/2, mu_filter_geo.UpstreamBarY/2, mu_filter_geo.UpstreamBarZ/2
        else: # DS
            dx, dy, dz = mu_filter_geo.DownstreamBarX_ver/2, mu_filter_geo.DownstreamBarY/2, mu_filter_geo.DownstreamBarZ/2

        graph.SetPoint(idx, (pos_a.Z() + pos_b.Z())/2, pos_a.X() if view == 'X' else pos_a.Y())
        graph.SetPointError(idx, dz, dx if view == 'X' else dy)
        counts[view]['MuFilter'] += 1

    for view, p_idx in [('X', 1), ('Y', 2)]:
        histograms['simpleDisplay'].cd(p_idx)
        for system, color, style in [('Scifi', ROOT.kBlue+2, 20), ('MuFilter', ROOT.kRed+1, 21)]:
            g = hit_graphs[view][system]
            g.SetMarkerStyle(style); g.SetMarkerSize(0.8); g.SetMarkerColor(color)
            g.Draw('sameP')
            histograms[f'hits_{system.lower()[:2]}_{view}'] = g

        proj_name = f'{view}Z'
        for i, (slope, intercept) in enumerate(hough_lines[proj_name]):
            z_min, z_max = histograms['zmin'], histograms['zmax']
            line = ROOT.TLine(z_min, slope * z_min + intercept, z_max, slope * z_max + intercept)
            line.SetLineColor(ROOT.kCyan+2); line.SetLineWidth(2)
            line.Draw("same")
            histograms[f'line_{proj_name}_{i}'] = line

def run_hough_transform(muon_reco_task, event, geo, z_vtx_min=None, z_vtx_max=None):
    """Identify tracks using Hough transform with vertex constraints."""
    hit_collection = {
        "pos": [[], [], []], "d": [[], [], []],
        "vert": [], "system": [], "detectorID": []
    }
    pos_a, pos_b = ROOT.TVector3(), ROOT.TVector3()

    for hit in event.Digi_ScifiHits:
        if not hit.isValid():
            continue

        geo.modules['Scifi'].GetSiPMPosition(hit.GetDetectorID(), pos_a, pos_b)
        for i in range(3):
            hit_collection["pos"][i].append(pos_a[i])

        hit_collection["d"][0].append(muon_reco_task.Scifi_dx)
        hit_collection["d"][1].append(muon_reco_task.Scifi_dy)
        hit_collection["d"][2].append(muon_reco_task.Scifi_dz)
        hit_collection["vert"].append(hit.isVertical())
        hit_collection["system"].append(0)
        hit_collection["detectorID"].append(hit.GetDetectorID())

    if not hit_collection['pos'][0]:
        return 0, {'XZ': [], 'YZ': []}

    for k in ["pos", "d"]:
        hit_collection[k] = np.array(hit_collection[k], dtype=np.float32)
    for k, dt in [("vert", np.bool_), ("system", np.int32), ("detectorID", np.int32)]:
        hit_collection[k] = np.array(hit_collection[k], dtype=dt)

    counts, lines = {'XZ': 0, 'YZ': 0}, {'XZ': [], 'YZ': []}

    for projection_name in ['XZ', 'YZ']:
        is_vertical, axis = (projection_name == 'XZ'), (0 if projection_name == 'XZ' else 1)
        hough_object = muon_reco_task.h_ZX if is_vertical else muon_reco_task.h_ZY
        hits_used = np.zeros(len(hit_collection['pos'][0]), dtype=np.bool_)

        valid_lines_found = 0
        attempts = 0
        max_attempts = 20

        while valid_lines_found < muon_reco_task.max_reco_muons and attempts < max_attempts:
            attempts += 1
            mask = np.logical_and(hit_collection['vert'] == is_vertical, ~hits_used)
            if not np.any(mask): break

            # Perform Hough fit
            fit_result = hough_object.fit_randomize(
                np.dstack([hit_collection['pos'][2][mask], hit_collection['pos'][axis][mask]])[0],
                np.dstack([hit_collection['d'][2][mask], hit_collection['d'][axis][mask]])[0],
                muon_reco_task.n_random, False, False
            )

            if fit_result[0] in [-1, -999]: break

            new_slope, new_intercept = fit_result[0], fit_result[1]
            related_hits = SndlhcMuonReco.hit_finder(
                new_slope, new_intercept,
                np.dstack([hit_collection['pos'][2][mask], hit_collection['pos'][axis][mask]]),
                np.dstack([hit_collection['d'][2][mask], hit_collection['d'][axis][mask]]),
                muon_reco_task.tolerance
            )

            if len(related_hits) == 0: break

            if SndlhcMuonReco.numPlanesHit(hit_collection['system'][mask][related_hits], hit_collection['detectorID'][mask][related_hits]) >= muon_reco_task.min_planes_hit:
                skip_track = False
                conflict_params = None

                for existing_line in lines[projection_name]:
                    ext_m, ext_c = existing_line[0], existing_line[1]
                    if abs(new_slope - ext_m) > 1e-6:
                        z_vertex = (ext_c - new_intercept) / (new_slope - ext_m)
                        if (z_vtx_min is not None and z_vertex < z_vtx_min) or (z_vtx_max is not None and z_vertex > z_vtx_max):
                            skip_track = True
                            conflict_params = (ext_m, ext_c)
                            break

                if skip_track:
                    if conflict_params:
                        # Find closest hit to conflicting track and exclude it
                        global_indices = np.where(mask)[0][related_hits]
                        z_bad, c_bad = hit_collection['pos'][2][global_indices], hit_collection['pos'][axis][global_indices]
                        dist = np.abs(c_bad - (conflict_params[0] * z_bad + conflict_params[1]))
                        hits_used[global_indices[np.argmin(dist)]] = True
                    else:
                        hits_used[np.where(mask)[0][related_hits]] = True
                    continue

                # Valid track found
                counts[projection_name] += 1
                lines[projection_name].append(fit_result)

                # Hit exclusion logic (transverse + Z plane)
                selected_global_idx = np.where(mask)[0][related_hits]
                projection_idx = np.where(hit_collection['vert'] == is_vertical)[0]

                z_sel, c_sel = hit_collection['pos'][2][selected_global_idx], hit_collection['pos'][axis][selected_global_idx]
                z_all, c_all = hit_collection['pos'][2][projection_idx], hit_collection['pos'][axis][projection_idx]

                dz = z_all[:, np.newaxis] - z_sel[np.newaxis, :]
                dc = c_all[:, np.newaxis] - c_sel[np.newaxis, :]
                close_hits_mask = np.any((np.abs(dz) < 1e-3) & (np.abs(dc) < muon_reco_task.tolerance), axis=1)
                hits_used[projection_idx[close_hits_mask]] = True
                valid_lines_found += 1
            else:
                break

    return max(counts.values()), lines

def main():
    args = get_arguments()
    start_time_process = time.time()

    # Open input file and initialize geometry
    input_file = ROOT.TFile.Open(args.input_file if args.input_file.endswith(".root") else args.input_file + ".root")
    input_tree = input_file.Get("rawConv")
    input_tree.GetEvent(0)
    run_number = input_tree.EventHeader.GetRunId()

    geo = SndlhcGeo.GeoInterface(ROOT.snd.analysis_tools.GetGeoPath(run_number))
    initialize_event_display(geo)

    fair_run = ROOT.FairRunAna()
    io_manager = ROOT.FairRootManager.Instance()
    io_manager.SetTreeName("rawConv")
    fair_run.SetSource(ROOT.FairFileSource(input_file))
    fair_run.SetSink(ROOT.FairRootFileSink(ROOT.TMemFile('dummy','CREATE')))

    muon_reco_task = SndlhcMuonReco.MuonReco()
    fair_run.AddTask(muon_reco_task)
    muon_reco_task.SetParFile(args.parFile)
    muon_reco_task.SetHoughSpaceFormat("linearSlopeIntercept")
    muon_reco_task.SetTrackingCase('muon_trident_Sf')
    fair_run.Init()

    input_tree = io_manager.GetInTree()
    if input_tree.GetBranch('Digi_MuFilterHit'):
        input_tree.Digi_MuFilterHits = input_tree.Digi_MuFilterHit

    output_root_file, output_tree = None, None
    if args.output_root:
        output_root_file = ROOT.TFile(args.output_root if args.output_root.endswith(".root") else f"{args.output_root}.root", "RECREATE")
        output_tree = input_tree.CloneTree(0)

    num_events_to_process = min(args.n_events, input_tree.GetEntries() - args.start_event)
    progress_step = 0

    with open(args.output_csv if args.output_csv.endswith(".csv") else f"{args.output_csv}.csv", mode='w', newline='') as csv_file:
        fieldnames = [
            'run', 'event_number', 'n_lines',
            'xz_m1', 'xz_c1', 'xz_m2', 'xz_c2', 'xz_m3', 'xz_c3',
            'yz_m1', 'yz_c1', 'yz_m2', 'yz_c2', 'yz_m3', 'yz_c3',
            'sum_hit_weight_density', 'sum_qdc'
        ]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        for entry_index in range(args.start_event, args.start_event + num_events_to_process):
            input_tree.GetEvent(entry_index)
            event_number = input_tree.EventHeader.GetEventNumber()

            if random.random() >= args.fraction:
                continue

            if ((entry_index - args.start_event) * 100 / num_events_to_process) >= progress_step:
                elapsed = int(time.time() - start_time_process)
                h, rem = divmod(elapsed, 3600)
                m, s = divmod(rem, 60)
                print(f"[{progress_step}%] \t {entry_index-args.start_event:,}/{num_events_to_process:,} \t {h:02d}:{m:02d}:{s:02d}")
                progress_step += 1

            if input_tree.EventHeader.ClassName() == 'SNDLHCEventHeader':
                geo.modules['Scifi'].InitEvent(input_tree.EventHeader)

            n_lines, track_lines = run_hough_transform(
                muon_reco_task, input_tree, geo,
                z_vtx_min=args.z_vtx_min, z_vtx_max=args.z_vtx_max
            )

            if n_lines >= 2:
                print(f"Match: Event {event_number}, Lines {n_lines}")
                if output_tree:
                    output_tree.Fill()
                    draw_event_hits_and_tracks(input_tree, geo, track_lines)
                    canvas = histograms['simpleDisplay']
                    canvas.SetName(f"c_Run{run_number}_{event_number}")
                    canvas.SetTitle(f"Event {event_number} - {n_lines} Lines")
                    output_root_file.cd()
                    canvas.Write()

            if n_lines >= 1:
                def get_line_params(projection, index):
                    if index < len(track_lines[projection]):
                        return track_lines[projection][index][0], track_lines[projection][index][1]
                    return np.nan, np.nan

                writer.writerow({
                    'run': run_number, 'event_number': event_number, 'n_lines': n_lines,
                    'xz_m1': get_line_params('XZ', 0)[0], 'xz_c1': get_line_params('XZ', 0)[1],
                    'xz_m2': get_line_params('XZ', 1)[0], 'xz_c2': get_line_params('XZ', 1)[1],
                    'xz_m3': get_line_params('XZ', 2)[0], 'xz_c3': get_line_params('XZ', 2)[1],
                    'yz_m1': get_line_params('YZ', 0)[0], 'yz_c1': get_line_params('YZ', 0)[1],
                    'yz_m2': get_line_params('YZ', 1)[0], 'yz_c2': get_line_params('YZ', 1)[1],
                    'yz_m3': get_line_params('YZ', 2)[0], 'yz_c3': get_line_params('YZ', 2)[1],
                    'sum_hit_weight_density': get_scifi_hit_density(input_tree.Digi_ScifiHits),
                    'sum_qdc': get_scifi_total_qdc(input_tree.Digi_ScifiHits)
                })
                csv_file.flush()

    if output_root_file:
        output_root_file.cd()
        output_tree.Write()
        output_root_file.Close()
    print(f"Finished in {time.time() - start_time_process:.2f}s.")

if __name__ == "__main__":
    main()
