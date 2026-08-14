#!/usr/bin/env python3
"""
filter_and_display_tridents.py
------------------------------
Generates 2D (XZ, YZ) and 3D event displays for all events in input ROOT files
without filtering, adopting the official SND@LHC 2dEventDisplay.py detector & hit
geometry methodology, and saving the displays directly into output ROOT files.
"""

import os
import sys
import glob
import time
import array
import argparse
import ROOT
import SndlhcGeo

ROOT.gROOT.SetBatch(True)

# Geometry and coordinate constants
Z_3D_MIN = 0.0
Z_MAX = 600.0
ROCK_BOUNDARY = 360.0

def make_array(lst):
    return array.array('d', lst)

def style_detector_volumes(vol):
    """Recursively styles TGeo volumes for standard 3D representation."""
    for i in range(vol.GetNdaughters()):
        dnode = vol.GetNode(i)
        dvol = dnode.GetVolume()
        ROOT.SetOwnership(dvol, False)
        name = dnode.GetName()
        vname = dvol.GetName()

        dvol.SetVisibility(True)
        if 'Wall' in name and 'border' not in name:
            dvol.SetLineColor(ROOT.kGray+2)
            dvol.SetFillColor(ROOT.kGray+1)
            dvol.SetTransparency(30)
        elif 'Scifi' in name or 'ScifiVolume' in vname:
            dvol.SetLineColor(ROOT.kCyan-2)
            dvol.SetFillColor(ROOT.kCyan-6)
            dvol.SetTransparency(20)
        elif 'volFeBlock' in name or 'volFeTarget' in name:
            dvol.SetLineColor(ROOT.kGreen-3)
            dvol.SetFillColor(ROOT.kGreen-7)
            dvol.SetTransparency(20)
        elif 'volMuUpstreamDet' in name or 'volMuDownstreamDet' in name:
            dvol.SetLineColor(ROOT.kBlue-3)
            dvol.SetFillColor(ROOT.kBlue-7)
            dvol.SetTransparency(20)
        elif 'Veto' in name:
            dvol.SetLineColor(ROOT.kOrange+2)
            dvol.SetFillColor(ROOT.kOrange-3)
            dvol.SetTransparency(20)
        else:
            dvol.SetLineColor(ROOT.kGray)
            dvol.SetTransparency(60)

        if dvol.GetNdaughters() > 0:
            style_detector_volumes(dvol)

def find_geofile(geofile_filename):
    """Searches multiple candidate locations for geofile to ensure it is always found."""
    if not geofile_filename:
        geofile_filename = "geofile_full.Ntuple-TGeant4_boost100.0.root"

    script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
    sndsw_root = os.environ.get("SNDSW_ROOT", "")

    candidates = [
        geofile_filename,
        os.path.basename(geofile_filename),
        os.path.join(script_dir, geofile_filename),
        os.path.join(script_dir, os.path.basename(geofile_filename)),
        os.path.join(sndsw_root, geofile_filename),
        os.path.join(sndsw_root, "python", os.path.basename(geofile_filename)),
        os.path.join(sndsw_root, "geofile_full.Ntuple-TGeant4_boost100.0.root"),
        "/afs/cern.ch/user/i/idioniso/snd_master/sndsw/python/geofile_full.Ntuple-TGeant4_boost100.0.root",
        "/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V0_2022.root"
    ]

    for c in candidates:
        if c and os.path.exists(c):
            return os.path.abspath(c)

    return geofile_filename

def load_detector_geometry(geofile_path):
    """
    Loads detector geometry using SndlhcGeo.GeoInterface.
    Extracts 2D bounding boxes and 3D visual structures.
    """
    geofile_path = find_geofile(geofile_path)
    if not os.path.exists(geofile_path):
        print(f"Warning: Geofile not found at '{geofile_path}'.")
        return [], None, None, None, None

    snd_geo = SndlhcGeo.GeoInterface(geofile_path)
    f_geo = ROOT.TFile.Open(geofile_path, "READ")
    if not f_geo or f_geo.IsZombie():
        print(f"Warning: Could not open geometry file '{geofile_path}'.")
        return [], None, None, None, None

    ggeo = ROOT.gGeoManager
    top_vol = ggeo.GetTopVolume()
    ROOT.SetOwnership(top_vol, False)
    style_detector_volumes(top_vol)

    scifi_module = snd_geo.modules.get('Scifi', None)
    nav = ggeo.GetCurrentNavigator()

    geo_elements = []

    def extract_node_bounds(node_path, label, color, alpha=0.25):
        if not nav.CheckPath(node_path):
            return
        nav.cd(node_path)
        node = nav.GetCurrentNode()
        vol = node.GetVolume()
        shape = vol.GetShape()
        dx, dy, dz = shape.GetDX(), shape.GetDY(), shape.GetDZ()
        ox, oy, oz = shape.GetOrigin()[0], shape.GetOrigin()[1], shape.GetOrigin()[2]

        corners = [
            array.array('d', [-dx + ox, -dy + oy, -dz + oz]),
            array.array('d', [ dx + ox, -dy + oy, -dz + oz]),
            array.array('d', [-dx + ox,  dy + oy, -dz + oz]),
            array.array('d', [ dx + ox,  dy + oy, -dz + oz]),
            array.array('d', [-dx + ox, -dy + oy,  dz + oz]),
            array.array('d', [ dx + ox, -dy + oy,  dz + oz]),
            array.array('d', [-dx + ox,  dy + oy,  dz + oz]),
            array.array('d', [ dx + ox,  dy + oy,  dz + oz]),
        ]

        master_pts = []
        for c in corners:
            m = array.array('d', [0.0, 0.0, 0.0])
            nav.LocalToMaster(c, m)
            master_pts.append(m)

        xs = [p[0] for p in master_pts]
        ys = [p[1] for p in master_pts]
        zs = [p[2] for p in master_pts]

        geo_elements.append({
            'name': node_path,
            'label': label,
            'color': color,
            'alpha': alpha,
            'x': (min(xs), max(xs)),
            'y': (min(ys), max(ys)),
            'z': (min(zs), max(zs)),
        })

    # Veto Detector
    for i in range(2):
        extract_node_bounds(f"/cave_1/Detector_0/volVeto_1/subVetoBox_{i}", "Veto Detector", ROOT.kOrange+7, 0.35)
    extract_node_bounds("/cave_1/Detector_0/volVeto_1/subVeto3Box_2", "Veto Detector", ROOT.kOrange+7, 0.35)

    # SciFi Target Emulsion Walls
    for i in range(5):
        extract_node_bounds(f"/cave_1/Detector_0/volTarget_1/volWallborder_{i}", "Target Emulsion Wall", ROOT.kGray+1, 0.30)

    # SciFi Stations
    for i in range(1, 6):
        extract_node_bounds(f"/cave_1/Detector_0/volTarget_1/ScifiVolume{i}_{i}000000", "SciFi Station", ROOT.kAzure-4, 0.40)

    # MuFilter Iron Blocks & Active Detector Planes
    for i in range(5):
        extract_node_bounds(f"/cave_1/Detector_0/volMuFilter_1/volFeBlock_{i}", "MuFilter Iron Block", ROOT.kGreen-6, 0.35)
        extract_node_bounds(f"/cave_1/Detector_0/volMuFilter_1/volMuUpstreamDet_{i}_{i+2}", "MuFilter Active Plane", ROOT.kBlue-4, 0.35)

    for i in range(4):
        extract_node_bounds(f"/cave_1/Detector_0/volMuFilter_1/volFeBlock_{i+5}", "MuFilter Iron Block", ROOT.kGreen-6, 0.35)
        extract_node_bounds(f"/cave_1/Detector_0/volMuFilter_1/volMuDownstreamDet_{i}_{i+7}", "MuFilter Active Plane", ROOT.kBlue-4, 0.35)

    extract_node_bounds("/cave_1/Detector_0/volMuFilter_1/volFeBlockEnd_1", "MuFilter Iron Block", ROOT.kGreen-6, 0.35)

    print(f"Loaded {len(geo_elements)} detector geometry elements from '{geofile_path}'.")
    return geo_elements, ggeo, f_geo, scifi_module, snd_geo

def build_mc_points_map(event):
    """
    Collects true simulated detector steps from ScifiPoint, MuFilterPoint, and EmulsionDetPoint.
    Maps trackID -> list of (z, x, y, pz, px, py) sorted along Z.
    """
    mc_points_map = {}
    for col_name in ['ScifiPoint', 'MuFilterPoint', 'EmulsionDetPoint']:
        if hasattr(event, col_name):
            col = getattr(event, col_name)
            for pt in col:
                trID = pt.GetTrackID()
                if trID not in mc_points_map:
                    mc_points_map[trID] = []
                mc_points_map[trID].append((pt.GetZ(), pt.GetX(), pt.GetY(), pt.GetPz(), pt.GetPx(), pt.GetPy()))
    for trID in mc_points_map:
        mc_points_map[trID].sort(key=lambda p: p[0])
    return mc_points_map

def get_track_trajectory_points(tr_i, track, z_min_plot, z_max_plot, mc_points_map=None):
    """
    Calculates polyline coordinates (x_arr, y_arr, z_arr) for a track.
    Uses true simulated Geant4 step coordinates in the detector to account for multiple
    Coulomb scattering through the rock. Falls back to straight line if no points present.
    """
    z0, x0, y0 = track.GetStartZ(), track.GetStartX(), track.GetStartY()
    pz, px, py = track.GetPz(), track.GetPy(), track.GetPy()

    det_pts = mc_points_map.get(tr_i, []) if (mc_points_map is not None and tr_i is not None and tr_i >= 0) else []
    pts = []

    if det_pts:
        first_pt = det_pts[0]
        if z0 >= z_min_plot:
            pts.append((z0, x0, y0))
        else:
            pz_f, px_f, py_f = first_pt[3], first_pt[4], first_pt[5]
            if pz_f != 0:
                x_enter = first_pt[1] + (px_f / pz_f) * (z_min_plot - first_pt[0])
                y_enter = first_pt[2] + (py_f / pz_f) * (z_min_plot - first_pt[0])
            else:
                x_enter, y_enter = first_pt[1], first_pt[2]
            pts.append((z_min_plot, x_enter, y_enter))

        for p in det_pts:
            if not pts or abs(p[0] - pts[-1][0]) > 0.05:
                pts.append((p[0], p[1], p[2]))

        last_pt = det_pts[-1]
        pz_l, px_l, py_l = last_pt[3], last_pt[4], last_pt[5]
        if last_pt[0] < z_max_plot and pz_l > 0:
            x_exit = last_pt[1] + (px_l / pz_l) * (z_max_plot - last_pt[0])
            y_exit = last_pt[2] + (py_l / pz_l) * (z_max_plot - last_pt[0])
            pts.append((z_max_plot, x_exit, y_exit))
    else:
        if pz != 0:
            x_start = x0 if z0 >= z_min_plot else x0 + (px / pz) * (z_min_plot - z0)
            y_start = y0 if z0 >= z_min_plot else y0 + (py / pz) * (z_min_plot - z0)
            z_s = max(z0, z_min_plot)
            x_end = x0 + (px / pz) * (z_max_plot - z0)
            y_end = y0 + (py / pz) * (z_max_plot - z0)
            pts.append((z_s, x_start, y_start))
            pts.append((z_max_plot, x_end, y_end))
        else:
            pts.append((z0, x0, y0))
            pts.append((z_max_plot, x0, y0))

    z_arr = [p[0] for p in pts]
    x_arr = [p[1] for p in pts]
    y_arr = [p[2] for p in pts]
    return x_arr, y_arr, z_arr

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

def draw_detector_geometry(geo_elements, z_range, x_range, y_range, ggeo=None, file_tag="", evt_num=0, i_event=0):
    """
    Draws the detector geometry volumes on XZ, YZ, and 3D views.
    Silences line/marker rendering on dummy frames to avoid diagonal artifacts.
    """
    z_min_plot, z_max_plot = z_range
    x_min_plot, x_max_plot = x_range
    y_min_plot, y_max_plot = y_range

    # ----------------- XZ Projection -----------------
    c_xz = ROOT.TCanvas(f"XZ_Ev_{evt_num}_F_{file_tag}_Idx_{i_event}", f"XZ Projection (Event #{evt_num})", 950, 650)
    ROOT.SetOwnership(c_xz, False)
    mg_xz = ROOT.TMultiGraph()
    ROOT.SetOwnership(mg_xz, False)
    mg_xz.SetTitle(";Z [cm];X [cm]")

    gr_dummy_xz = ROOT.TGraph(2)
    ROOT.SetOwnership(gr_dummy_xz, False)
    gr_dummy_xz.SetPoint(0, z_min_plot, x_min_plot)
    gr_dummy_xz.SetPoint(1, z_max_plot, x_max_plot)
    gr_dummy_xz.SetLineColor(0)
    gr_dummy_xz.SetLineWidth(0)
    gr_dummy_xz.SetMarkerColor(0)
    gr_dummy_xz.SetMarkerSize(0)
    mg_xz.Add(gr_dummy_xz)
    mg_xz.Draw("AP")
    mg_xz.GetXaxis().SetLimits(z_min_plot, z_max_plot)
    mg_xz.GetYaxis().SetRangeUser(x_min_plot, x_max_plot)
    c_xz.Update()

    boxes_xz = []
    legend_xz = ROOT.TLegend(0.12, 0.70, 0.44, 0.89)
    ROOT.SetOwnership(legend_xz, False)
    legend_xz.SetBorderSize(1)
    legend_xz.SetFillStyle(1001)
    legend_xz.SetFillColorAlpha(ROOT.kWhite, 0.88)
    legend_xz.SetTextFont(42)
    legend_xz.SetTextSize(0.024)

    for elem in geo_elements:
        box = ROOT.TBox(elem['z'][0], elem['x'][0], elem['z'][1], elem['x'][1])
        ROOT.SetOwnership(box, False)
        box.SetFillStyle(1001)
        box.SetFillColorAlpha(elem['color'], elem['alpha'])
        box.SetLineColor(elem['color'])
        box.SetLineWidth(1)
        box.Draw("f same")
        box.Draw("l same")
        boxes_xz.append(box)

    # ----------------- YZ Projection -----------------
    c_yz = ROOT.TCanvas(f"YZ_Ev_{evt_num}_F_{file_tag}_Idx_{i_event}", f"YZ Projection (Event #{evt_num})", 950, 650)
    ROOT.SetOwnership(c_yz, False)
    mg_yz = ROOT.TMultiGraph()
    ROOT.SetOwnership(mg_yz, False)
    mg_yz.SetTitle(";Z [cm];Y [cm]")

    gr_dummy_yz = ROOT.TGraph(2)
    ROOT.SetOwnership(gr_dummy_yz, False)
    gr_dummy_yz.SetPoint(0, z_min_plot, y_min_plot)
    gr_dummy_yz.SetPoint(1, z_max_plot, y_max_plot)
    gr_dummy_yz.SetLineColor(0)
    gr_dummy_yz.SetLineWidth(0)
    gr_dummy_yz.SetMarkerColor(0)
    gr_dummy_yz.SetMarkerSize(0)
    mg_yz.Add(gr_dummy_yz)
    mg_yz.Draw("AP")
    mg_yz.GetXaxis().SetLimits(z_min_plot, z_max_plot)
    mg_yz.GetYaxis().SetRangeUser(y_min_plot, y_max_plot)
    c_yz.Update()

    boxes_yz = []
    legend_yz = ROOT.TLegend(0.12, 0.70, 0.44, 0.89)
    ROOT.SetOwnership(legend_yz, False)
    legend_yz.SetBorderSize(1)
    legend_yz.SetFillStyle(1001)
    legend_yz.SetFillColorAlpha(ROOT.kWhite, 0.88)
    legend_yz.SetTextFont(42)
    legend_yz.SetTextSize(0.024)

    for elem in geo_elements:
        box = ROOT.TBox(elem['z'][0], elem['y'][0], elem['z'][1], elem['y'][1])
        ROOT.SetOwnership(box, False)
        box.SetFillStyle(1001)
        box.SetFillColorAlpha(elem['color'], elem['alpha'])
        box.SetLineColor(elem['color'])
        box.SetLineWidth(1)
        box.Draw("f same")
        box.Draw("l same")
        boxes_yz.append(box)

    # ----------------- 3D View -----------------
    c_3d = ROOT.TCanvas(f"3D_Ev_{evt_num}_F_{file_tag}_Idx_{i_event}", f"3D View (Event #{evt_num})", 850, 850)
    ROOT.SetOwnership(c_3d, False)

    h3 = ROOT.TH3F(f"axis_box_3d_{evt_num}_F_{file_tag}_Idx_{i_event}", ";Z [cm];X [cm];Y [cm]",
                   10, Z_3D_MIN, Z_MAX,
                   10, x_min_plot, x_max_plot,
                   10, y_min_plot, y_max_plot)
    ROOT.SetOwnership(h3, False)
    h3.SetStats(0)
    h3.Draw()

    geo_3d_lines = []
    for elem in geo_elements:
        xmin, xmax = elem['x']
        ymin, ymax = elem['y']
        zmin, zmax = elem['z']
        corners = [
            (zmin, xmin, ymin), (zmax, xmin, ymin), (zmax, xmax, ymin), (zmin, xmax, ymin), (zmin, xmin, ymin),
            (zmin, xmin, ymax), (zmax, xmin, ymax), (zmax, xmax, ymax), (zmin, xmax, ymax), (zmin, xmin, ymax)
        ]
        pl = ROOT.TPolyLine3D(len(corners))
        ROOT.SetOwnership(pl, False)
        for idx, (zc, xc, yc) in enumerate(corners):
            pl.SetPoint(idx, zc, xc, yc)
        pl.SetLineColor(elem['color'])
        pl.SetLineWidth(1)
        pl.Draw("same")
        geo_3d_lines.append(pl)

        for (zc1, xc1, yc1), (zc2, xc2, yc2) in [
            ((zmax, xmin, ymin), (zmax, xmin, ymax)),
            ((zmax, xmax, ymin), (zmax, xmax, ymax)),
            ((zmin, xmax, ymin), (zmin, xmax, ymax))
        ]:
            pl_edge = ROOT.TPolyLine3D(2)
            ROOT.SetOwnership(pl_edge, False)
            pl_edge.SetPoint(0, zc1, xc1, yc1)
            pl_edge.SetPoint(1, zc2, xc2, yc2)
            pl_edge.SetLineColor(elem['color'])
            pl_edge.SetLineWidth(1)
            pl_edge.Draw("same")
            geo_3d_lines.append(pl_edge)

    return {
        'c_xz': c_xz, 'c_yz': c_yz, 'c_3d': c_3d,
        'mg_xz': mg_xz, 'mg_yz': mg_yz, 'h3': h3,
        'legend_xz': legend_xz, 'legend_yz': legend_yz,
        'boxes_xz': boxes_xz, 'boxes_yz': boxes_yz,
        'geo_3d_lines': geo_3d_lines
    }

def draw_scifi_hits(event, scifi_module, geo_ctx, snd_geo=None, marker_color=ROOT.kBlue+2, marker_style=20, marker_size=1.2):
    """
    Draws Digi_ScifiHits and Digi_MuFilterHits adopting the official 2dEventDisplay.py methodology:
    - Queries SiPMPosition using scifi_module.GetSiPMPosition(detID, A, B).
    - Vertical hits -> XZ view with channel width error bars (TGraphErrors) / marker.
    - Horizontal hits -> YZ view with channel width error bars (TGraphErrors) / marker.
    - 3D fiber segments -> 3D view.
    """
    if not scifi_module or not hasattr(event, "Digi_ScifiHits") or len(event.Digi_ScifiHits) == 0:
        return {}

    c_xz = geo_ctx['c_xz']
    c_yz = geo_ctx['c_yz']
    c_3d = geo_ctx['c_3d']
    legend_xz = geo_ctx['legend_xz']
    legend_yz = geo_ctx['legend_yz']

    A = ROOT.TVector3()
    B = ROOT.TVector3()

    gr_xz = ROOT.TGraphErrors()
    ROOT.SetOwnership(gr_xz, False)
    gr_yz = ROOT.TGraphErrors()
    ROOT.SetOwnership(gr_yz, False)

    si = snd_geo.snd_geo.Scifi if (snd_geo and hasattr(snd_geo, "snd_geo")) else None
    sY = si.channel_width if si else 0.025
    sZ = si.scifimat_z if si else 0.135

    poly3d_list = []
    n_x = 0
    n_y = 0

    for hit in event.Digi_ScifiHits:
        detID = hit.GetDetectorID()
        scifi_module.GetSiPMPosition(detID, A, B)

        z_mid = (A.Z() + B.Z()) / 2.0
        x_mid = (A.X() + B.X()) / 2.0
        y_mid = (A.Y() + B.Y()) / 2.0

        if hit.isVertical():
            gr_xz.SetPoint(n_x, z_mid, x_mid)
            gr_xz.SetPointError(n_x, sZ, sY)
            n_x += 1
        else:
            gr_yz.SetPoint(n_y, z_mid, y_mid)
            gr_yz.SetPointError(n_y, sZ, sY)
            n_y += 1

        pl3d = ROOT.TPolyLine3D(2)
        ROOT.SetOwnership(pl3d, False)
        pl3d.SetPoint(0, A.Z(), A.X(), A.Y())
        pl3d.SetPoint(1, B.Z(), B.X(), B.Y())
        pl3d.SetLineColor(marker_color)
        pl3d.SetLineWidth(2)
        poly3d_list.append(pl3d)

    if n_x > 0:
        c_xz.cd()
        gr_xz.SetMarkerStyle(marker_style)
        gr_xz.SetMarkerSize(marker_size)
        gr_xz.SetMarkerColor(marker_color)
        gr_xz.SetLineColor(marker_color)
        gr_xz.Draw("sameP")
        legend_xz.AddEntry(gr_xz, "SciFi Hits (Digi_ScifiHits)", "p")

    if n_y > 0:
        c_yz.cd()
        gr_yz.SetMarkerStyle(marker_style)
        gr_yz.SetMarkerSize(marker_size)
        gr_yz.SetMarkerColor(marker_color)
        gr_yz.SetLineColor(marker_color)
        legend_yz.AddEntry(gr_yz, "SciFi Hits (Digi_ScifiHits)", "p")
        gr_yz.Draw("sameP")

    c_3d.cd()
    for pl in poly3d_list:
        pl.Draw("same")

    return {
        'gr_xz': gr_xz,
        'gr_yz': gr_yz,
        'poly3d_list': poly3d_list
    }

def draw_mctracks(event, primary_idx, primary_track, inducing_idx, inducing_track, my_daughters, geo_ctx, z_range, show_all_mctracks=False, mc_points_map=None):
    """
    Draws MCTracks (Primary Muon, Inducing Secondary, Trident Daughters, and optional background MCTracks)
    on XZ, YZ, and 3D canvases.
    Accurately accounts for Geant4 multiple Coulomb scattering in rock by following true MC detector points.
    """
    c_xz = geo_ctx['c_xz']
    c_yz = geo_ctx['c_yz']
    c_3d = geo_ctx['c_3d']
    legend_xz = geo_ctx['legend_xz']
    legend_yz = geo_ctx['legend_yz']

    z_min_2d_plot, z_max_2d_plot = z_range
    is_secondary_inducing = (inducing_track is not None and primary_track is not None and inducing_track != primary_track)

    if mc_points_map is None:
        mc_points_map = build_mc_points_map(event)

    x_p, y_p, z_p = get_track_trajectory_points(primary_idx, primary_track, z_min_2d_plot, Z_MAX, mc_points_map) if primary_track else ([], [], [])

    x_ind, y_ind, z_ind = None, None, None
    if is_secondary_inducing and inducing_track:
        ind_z_start = max(inducing_track.GetStartZ(), z_min_2d_plot)
        x_ind, y_ind, z_ind = get_track_trajectory_points(inducing_idx, inducing_track, ind_z_start, Z_MAX, mc_points_map)

    d_trajs = []
    for d in my_daughters:
        d_i = d[0]
        d_pdg = d[2]
        d_track = d[3]
        d_z_start = max(d_track.GetStartZ(), z_min_2d_plot)
        x_d, y_d, z_d = get_track_trajectory_points(d_i, d_track, d_z_start, Z_MAX, mc_points_map)
        d_trajs.append((d_pdg, x_d, y_d, z_d))

    other_trajs_2d = []
    other_trajs_3d = []
    if show_all_mctracks and hasattr(event, 'MCTrack'):
        signal_track_indices = set()
        if primary_idx is not None and primary_idx >= 0:
            signal_track_indices.add(primary_idx)
        if inducing_idx is not None and inducing_idx >= 0:
            signal_track_indices.add(inducing_idx)
        for d in my_daughters:
            signal_track_indices.add(d[0])

        for tr_i, tr_obj in enumerate(event.MCTrack):
            if tr_i in signal_track_indices:
                continue
            z0 = tr_obj.GetStartZ()
            pz = tr_obj.GetPz()

            children = [c for c in event.MCTrack if c.GetMotherId() == tr_i]
            if children:
                earliest_child = min(children, key=lambda c: c.GetStartZ())
                z_end = earliest_child.GetStartZ()
            else:
                if pz > 0:
                    z_end = Z_MAX
                elif pz < 0:
                    z_end = max(0.0, z0 - 10.0)
                else:
                    z_end = z0

            if z_end > z_min_2d_plot or z0 > z_min_2d_plot:
                x_oth, y_oth, z_oth = get_track_trajectory_points(tr_i, tr_obj, max(z0, z_min_2d_plot), z_end, mc_points_map)
                if len(z_oth) >= 2:
                    other_trajs_2d.append((x_oth, y_oth, z_oth))
                    x_oth_3d, y_oth_3d, z_oth_3d = get_track_trajectory_points(tr_i, tr_obj, max(z0, Z_3D_MIN), min(z_end, Z_MAX), mc_points_map)
                    if len(z_oth_3d) >= 2:
                        other_trajs_3d.append((x_oth_3d, y_oth_3d, z_oth_3d))

    # --- Draw Other Tracks (behind signal tracks) ---
    gr_other_xz = []
    gr_other_yz = []
    if show_all_mctracks:
        c_xz.cd()
        for x_ends, y_ends, z_ends in other_trajs_2d:
            gr_oth = ROOT.TGraph(len(z_ends), make_array(z_ends), make_array(x_ends))
            ROOT.SetOwnership(gr_oth, False)
            gr_oth.SetLineColorAlpha(ROOT.kGray+1, 0.4)
            gr_oth.SetLineWidth(1)
            gr_oth.Draw("L same")
            gr_other_xz.append(gr_oth)

        c_yz.cd()
        for x_ends, y_ends, z_ends in other_trajs_2d:
            gr_oth = ROOT.TGraph(len(z_ends), make_array(z_ends), make_array(y_ends))
            ROOT.SetOwnership(gr_oth, False)
            gr_oth.SetLineColorAlpha(ROOT.kGray+1, 0.4)
            gr_oth.SetLineWidth(1)
            gr_oth.Draw("L same")
            gr_other_yz.append(gr_oth)

        c_3d.cd()
        for x_ends, y_ends, z_ends in other_trajs_3d:
            pl_oth = ROOT.TPolyLine3D(len(z_ends))
            ROOT.SetOwnership(pl_oth, False)
            for idx_p in range(len(z_ends)):
                pl_oth.SetPoint(idx_p, z_ends[idx_p], x_ends[idx_p], y_ends[idx_p])
            pl_oth.SetLineColorAlpha(ROOT.kGray+1, 0.4)
            pl_oth.SetLineWidth(1)
            pl_oth.Draw("same")

    # --- Draw Primary Track ---
    gr_p_xz, gr_p_yz, pl_p_3d = None, None, None
    if primary_track and len(z_p) >= 2:
        c_xz.cd()
        gr_p_xz = ROOT.TGraph(len(z_p), make_array(z_p), make_array(x_p))
        ROOT.SetOwnership(gr_p_xz, False)
        gr_p_xz.SetLineColor(ROOT.kBlack)
        gr_p_xz.SetLineWidth(2)
        gr_p_xz.Draw("L same")
        legend_xz.AddEntry(gr_p_xz, "Primary #mu^{-}", "l")

        c_yz.cd()
        gr_p_yz = ROOT.TGraph(len(z_p), make_array(z_p), make_array(y_p))
        ROOT.SetOwnership(gr_p_yz, False)
        gr_p_yz.SetLineColor(ROOT.kBlack)
        gr_p_yz.SetLineWidth(2)
        gr_p_yz.Draw("L same")
        legend_yz.AddEntry(gr_p_yz, "Primary #mu^{-}", "l")

        c_3d.cd()
        x_p_3d, y_p_3d, z_p_3d = get_track_trajectory_points(primary_idx, primary_track, Z_3D_MIN, Z_MAX, mc_points_map)
        pl_p_3d = ROOT.TPolyLine3D(len(z_p_3d))
        ROOT.SetOwnership(pl_p_3d, False)
        for idx_p in range(len(z_p_3d)):
            pl_p_3d.SetPoint(idx_p, z_p_3d[idx_p], x_p_3d[idx_p], y_p_3d[idx_p])
        pl_p_3d.SetLineColor(ROOT.kBlack)
        pl_p_3d.SetLineWidth(3)
        pl_p_3d.Draw("same")

    # --- Draw Secondary Inducing Track ---
    gr_ind_xz, gr_ind_yz, pl_ind_3d = None, None, None
    if is_secondary_inducing and x_ind is not None and len(z_ind) >= 2:
        c_xz.cd()
        gr_ind_xz = ROOT.TGraph(len(z_ind), make_array(z_ind), make_array(x_ind))
        ROOT.SetOwnership(gr_ind_xz, False)
        gr_ind_xz.SetLineColor(ROOT.kMagenta+2)
        gr_ind_xz.SetLineWidth(2)
        gr_ind_xz.Draw("L same")
        legend_xz.AddEntry(gr_ind_xz, "Inducing Secondary #mu", "l")

        c_yz.cd()
        gr_ind_yz = ROOT.TGraph(len(z_ind), make_array(z_ind), make_array(y_ind))
        ROOT.SetOwnership(gr_ind_yz, False)
        gr_ind_yz.SetLineColor(ROOT.kMagenta+2)
        gr_ind_yz.SetLineWidth(2)
        gr_ind_yz.Draw("L same")
        legend_yz.AddEntry(gr_ind_yz, "Inducing Secondary #mu", "l")

        c_3d.cd()
        ind_z_start_3d = max(inducing_track.GetStartZ(), Z_3D_MIN)
        x_ind_3d, y_ind_3d, z_ind_3d = get_track_trajectory_points(inducing_idx, inducing_track, ind_z_start_3d, Z_MAX, mc_points_map)
        pl_ind_3d = ROOT.TPolyLine3D(len(z_ind_3d))
        ROOT.SetOwnership(pl_ind_3d, False)
        for idx_p in range(len(z_ind_3d)):
            pl_ind_3d.SetPoint(idx_p, z_ind_3d[idx_p], x_ind_3d[idx_p], y_ind_3d[idx_p])
        pl_ind_3d.SetLineColor(ROOT.kMagenta+2)
        pl_ind_3d.SetLineWidth(3)
        pl_ind_3d.Draw("same")

    # --- Draw Daughter Tracks ---
    gr_ds_xz = []
    gr_ds_yz = []
    pls_3d = []
    for d_pdg, x_d, y_d, z_d in d_trajs:
        if len(z_d) < 2:
            continue
        c_xz.cd()
        gr_x = ROOT.TGraph(len(z_d), make_array(z_d), make_array(x_d))
        ROOT.SetOwnership(gr_x, False)
        gr_x.SetLineWidth(2)
        gr_x.SetLineStyle(2)
        gr_x.SetLineColor(ROOT.kBlue if d_pdg == 13 else ROOT.kRed)
        gr_x.Draw("L same")
        gr_ds_xz.append(gr_x)

        c_yz.cd()
        gr_y = ROOT.TGraph(len(z_d), make_array(z_d), make_array(y_d))
        ROOT.SetOwnership(gr_y, False)
        gr_y.SetLineWidth(2)
        gr_y.SetLineStyle(2)
        gr_y.SetLineColor(ROOT.kBlue if d_pdg == 13 else ROOT.kRed)
        gr_y.Draw("L same")
        gr_ds_yz.append(gr_y)

        c_3d.cd()
        d_z_start = max(z_d[0], Z_3D_MIN)
        pl_3 = ROOT.TPolyLine3D(len(z_d))
        ROOT.SetOwnership(pl_3, False)
        for idx_p in range(len(z_d)):
            pl_3.SetPoint(idx_p, z_d[idx_p], x_d[idx_p], y_d[idx_p])
        pl_3.SetLineColor(ROOT.kBlue if d_pdg == 13 else ROOT.kRed)
        pl_3.SetLineStyle(2)
        pl_3.SetLineWidth(2)
        pl_3.Draw("same")
        pls_3d.append(pl_3)

    if gr_ds_xz:
        legend_xz.AddEntry(gr_ds_xz[0], "Trident #mu / #mu^{+} daughters", "l")
        legend_yz.AddEntry(gr_ds_yz[0], "Trident #mu / #mu^{+} daughters", "l")

    if show_all_mctracks and gr_other_xz:
        legend_xz.AddEntry(gr_other_xz[0], "Other MC Tracks", "l")
        legend_yz.AddEntry(gr_other_yz[0], "Other MC Tracks", "l")

    return {
        'gr_p_xz': gr_p_xz, 'gr_p_yz': gr_p_yz, 'pl_p_3d': pl_p_3d,
        'gr_ind_xz': gr_ind_xz, 'gr_ind_yz': gr_ind_yz, 'pl_ind_3d': pl_ind_3d,
        'gr_ds_xz': gr_ds_xz, 'gr_ds_yz': gr_ds_yz, 'pls_3d': pls_3d,
        'gr_other_xz': gr_other_xz, 'gr_other_yz': gr_other_yz
    }

def create_event_display(event, evt_num, i_event, geo_elements, ggeo, scifi_module, snd_geo,
                         primary_idx, primary_track, inducing_idx, inducing_track, my_daughters, file_tag="",
                         show_all_mctracks=False, draw_hits=True, draw_tracks=True):
    """
    High-level orchestrator that creates canvases, calls draw_detector_geometry,
    draw_scifi_hits, and draw_mctracks, draws legends/banners, and returns the canvases.
    """
    is_secondary_inducing = (inducing_track is not None and primary_track is not None and inducing_track != primary_track)

    mc_points_map = build_mc_points_map(event)

    z_min_2d = 250.0
    x_p, y_p, z_p = get_track_trajectory_points(primary_idx, primary_track, z_min_2d, Z_MAX, mc_points_map) if primary_track else ([], [], [])

    x_ind, y_ind, z_ind = None, None, None
    if is_secondary_inducing and inducing_track:
        ind_z_start = max(inducing_track.GetStartZ(), z_min_2d)
        x_ind, y_ind, z_ind = get_track_trajectory_points(inducing_idx, inducing_track, ind_z_start, Z_MAX, mc_points_map)

    d_trajs = []
    for d in my_daughters:
        d_i = d[0]
        d_pdg = d[2]
        d_track = d[3]
        d_z_start = max(d_track.GetStartZ(), z_min_2d)
        x_d, y_d, z_d = get_track_trajectory_points(d_i, d_track, d_z_start, Z_MAX, mc_points_map)
        d_trajs.append((d_pdg, x_d, y_d, z_d))

    all_x = list(x_p) + [x for d in d_trajs for x in d[1]] + [elem['x'][0] for elem in geo_elements] + [elem['x'][1] for elem in geo_elements]
    all_y = list(y_p) + [y for d in d_trajs for y in d[2]] + [elem['y'][0] for elem in geo_elements] + [elem['y'][1] for elem in geo_elements]
    all_z = list(z_p) + [z for d in d_trajs for z in d[3]] + [elem['z'][0] for elem in geo_elements] + [elem['z'][1] for elem in geo_elements]
    if is_secondary_inducing and x_ind is not None:
        all_x += list(x_ind)
        all_y += list(y_ind)
        all_z += list(z_ind)

    x_min_plot, x_max_plot = min(all_x) - 5, max(all_x) + 5
    y_min_plot, y_max_plot = min(all_y) - 5, max(all_y) + 5
    z_min_2d_plot, z_max_2d_plot = min(all_z) - 5, max(all_z) + 5

    z_vert_val = my_daughters[0][3].GetStartZ() if my_daughters else z_min_2d
    tag_suffix = "_SecMuon" if is_secondary_inducing else ""
    evt_dir_name = f"Event_{evt_num}_Zvert{int(z_vert_val)}cm{tag_suffix}"

    banner_text = f"Event #{evt_num}  (Entry {i_event})"
    if is_secondary_inducing:
        banner_text += "   #color[616]{[Secondary Muon Induced]}"

    # 1. Draw Detector Geometry
    geo_ctx = draw_detector_geometry(
        geo_elements,
        z_range=(z_min_2d_plot, z_max_2d_plot),
        x_range=(x_min_plot, x_max_plot),
        y_range=(y_min_plot, y_max_plot),
        ggeo=ggeo, file_tag=file_tag, evt_num=evt_num, i_event=i_event
    )

    # 2. Draw SciFi Hits (adopting official 2dEventDisplay.py style)
    scifi_ctx = {}
    if draw_hits and scifi_module:
        scifi_ctx = draw_scifi_hits(event, scifi_module, geo_ctx, snd_geo=snd_geo)

    # 3. Draw MCTracks
    mctrack_ctx = {}
    if draw_tracks:
        mctrack_ctx = draw_mctracks(
            event, primary_idx, primary_track, inducing_idx, inducing_track, my_daughters,
            geo_ctx, z_range=(z_min_2d_plot, z_max_2d_plot),
            show_all_mctracks=show_all_mctracks, mc_points_map=mc_points_map
        )

    # 4. Draw Legends and Banners (2dEventDisplay.py style)
    c_xz = geo_ctx['c_xz']
    c_yz = geo_ctx['c_yz']
    c_3d = geo_ctx['c_3d']

    for leg in [geo_ctx['legend_xz'], geo_ctx['legend_yz']]:
        n_entries = leg.GetNRows() if hasattr(leg, "GetNRows") else 3
        leg.SetTextFont(42)
        leg.SetTextSize(0.024)
        leg.SetBorderSize(1)
        leg.SetFillStyle(1001)
        leg.SetFillColorAlpha(ROOT.kWhite, 0.88)
        leg.SetMargin(0.25)

        x1 = 0.13
        x2 = 0.43
        y2 = 0.88
        height = max(0.06, min(0.24, 0.038 * max(1, n_entries) + 0.01))
        y1 = y2 - height
        leg.SetX1NDC(x1)
        leg.SetX2NDC(x2)
        leg.SetY1NDC(y1)
        leg.SetY2NDC(y2)

    banner_event = f"Event #{evt_num}  (Entry {i_event})"
    if is_secondary_inducing:
        banner_event += f"   #color[616]{{[Secondary Muon Induced, Z_{{vtx}}={z_vert_val:.1f} cm]}}"

    c_xz.cd()
    geo_ctx['legend_xz'].Draw()
    banner_xz = ROOT.TLatex()
    ROOT.SetOwnership(banner_xz, False)
    banner_xz.SetNDC()
    banner_xz.SetTextFont(42)
    banner_xz.SetTextSize(0.035)
    banner_xz.DrawLatex(0.12, 0.92, banner_event)

    c_yz.cd()
    geo_ctx['legend_yz'].Draw()
    banner_yz = ROOT.TLatex()
    ROOT.SetOwnership(banner_yz, False)
    banner_yz.SetNDC()
    banner_yz.SetTextFont(42)
    banner_yz.SetTextSize(0.035)
    banner_yz.DrawLatex(0.12, 0.92, banner_event)

    c_3d.cd()
    banner_3d = ROOT.TLatex()
    ROOT.SetOwnership(banner_3d, False)
    banner_3d.SetNDC()
    banner_3d.SetTextFont(42)
    banner_3d.SetTextSize(0.035)
    banner_3d.DrawLatex(0.12, 0.92, banner_event)

    return c_xz, c_yz, c_3d, evt_dir_name, is_secondary_inducing, z_vert_val

def process_single_file(input_file, output_file, geo_elements, ggeo, scifi_module, snd_geo, max_displays=0, show_all_mctracks=False):
    """
    Processes all events in the input file without filtering:
    - Copies 'cbmsim' and all metadata ('ShipGeo', etc.) to output_file.
    - Generates event displays for every event (or up to max_displays).
    - Stores event displays inside the output ROOT file under EventDisplays/ directory.
    """
    t0 = time.time()

    import shutil
    if os.path.abspath(input_file) != os.path.abspath(output_file):
        shutil.copyfile(input_file, output_file)

    f_out = ROOT.TFile.Open(output_file, "UPDATE")
    if not f_out or f_out.IsZombie():
        print(f"Error: Could not open output file '{output_file}'")
        return 0, 0, []

    tree = f_out.Get("cbmsim")
    if not tree:
        print(f"Error: 'cbmsim' TTree not found in '{output_file}'")
        f_out.Close()
        return 0, 0, []

    total_events = tree.GetEntries()
    file_tag = extract_file_tag(input_file)

    displays_dir = f_out.Get("EventDisplays")
    if not displays_dir:
        displays_dir = f_out.mkdir("EventDisplays")
    ROOT.SetOwnership(displays_dir, False)

    events_to_process = total_events if (max_displays <= 0 or max_displays > total_events) else max_displays
    events_plotted = 0
    sec_inducing_events = []

    for i_event in range(events_to_process):
        tree.GetEntry(i_event)
        event = tree

        evt_num = None
        if hasattr(event, "EventHeader") and event.EventHeader:
            try:
                evt_num = event.EventHeader.GetEventNumber()
            except Exception:
                pass
        if evt_num is None and hasattr(event, "MCEventHeader") and event.MCEventHeader:
            try:
                evt_num = event.MCEventHeader.GetEventID()
            except Exception:
                pass
        if evt_num is None or evt_num < 0:
            evt_num = i_event

        primary_idx = -1
        primary_track = None
        inducing_idx = -1
        inducing_track = None
        my_daughters = []

        if hasattr(event, "MCTrack"):
            incident_muons = {}
            pair_daughters = []
            for i, tr in enumerate(event.MCTrack):
                pdg = tr.GetPdgCode()
                mother_id = tr.GetMotherId()
                proc_id = tr.GetProcID()
                if abs(pdg) == 13:
                    incident_muons[i] = tr
                    if mother_id == -1:
                        primary_idx = i
                if abs(pdg) == 13 and proc_id == 5:
                    pair_daughters.append((i, mother_id, pdg, tr))

            for m_idx, m_tr in incident_muons.items():
                ds = [d for d in pair_daughters if d[1] == m_idx]
                has_m = any(d[2] == 13 for d in ds)
                has_p = any(d[2] == -13 for d in ds)
                if has_m and has_p:
                    inducing_idx = m_idx
                    my_daughters = ds
                    break

            if primary_idx != -1:
                primary_track = event.MCTrack[primary_idx]
            elif inducing_idx != -1:
                curr_tr = incident_muons[inducing_idx]
                while curr_tr.GetMotherId() != -1 and curr_tr.GetMotherId() in incident_muons:
                    curr_tr = incident_muons[curr_tr.GetMotherId()]
                primary_track = curr_tr

            if inducing_idx != -1:
                inducing_track = incident_muons[inducing_idx]

        c_xz, c_yz, c_3d, evt_dir_name, is_sec_ind, z_vert_val = create_event_display(
            event, evt_num, i_event, geo_elements, ggeo, scifi_module, snd_geo,
            primary_idx, primary_track, inducing_idx, inducing_track, my_daughters, file_tag=file_tag,
            show_all_mctracks=show_all_mctracks, draw_hits=True, draw_tracks=True
        )

        evt_dir = displays_dir.mkdir(evt_dir_name)
        if not evt_dir:
            evt_dir = displays_dir.Get(evt_dir_name)
        ROOT.SetOwnership(evt_dir, False)
        evt_dir.cd()

        if is_sec_ind:
            sec_inducing_events.append((evt_num, i_event, z_vert_val, primary_idx, inducing_idx, output_file))
            print(f"  [Secondary-Induced Trident] Event {evt_num} (Entry {i_event}): Pair induced by Secondary Muon Track {inducing_idx} at Z={z_vert_val:.1f} cm")

        c_xz.Write("XZ_Projection")
        c_yz.Write("YZ_Projection")
        c_3d.Write("3D_View")
        events_plotted += 1

    f_out.Close()

    elapsed = time.time() - t0
    print(f"  Processed {events_to_process:,} events in {elapsed:.1f}s | Displays generated: {events_plotted:,}")
    return total_events, events_plotted, sec_inducing_events

def main():
    default_input = "trimuon_filtered_*.root"
    default_output = "trimuon_displays_%s.root"
    default_geofile = "python/geofile_full.Ntuple-TGeant4_boost100.0.root"

    parser = argparse.ArgumentParser(description="Generate Event Displays for Events in Input ROOT Files (No Filtering)")
    parser.add_argument("-i", "--input", dest="input_pattern", default=default_input, help="Input ROOT file path or wildcard pattern (default: %(default)s)")
    parser.add_argument("-g", "--geofile", dest="geofile_filename", default=default_geofile, help="ROOT geometry file")
    parser.add_argument("-o", "--output", dest="output_pattern", default=default_output, help="Output ROOT file path or pattern with %%s for tag (default: %(default)s)")
    parser.add_argument("-n", "--max_displays", dest="max_displays", type=int, default=0, help="Max event displays to generate per file (0 = all events, default: %(default)s)")
    parser.add_argument("-a", "--all_mctracks", "--show_all_mctracks", dest="show_all_mctracks", action="store_true", default=False, help="Include all other MCTracks as thin, transparent grey lines")

    args = parser.parse_args()

    # Load geometry
    geo_elements, ggeo, f_geo, scifi_module, snd_geo = load_detector_geometry(args.geofile_filename)

    # Find matching files
    matched_files = sorted(glob.glob(args.input_pattern))
    if not matched_files:
        if os.path.exists(args.input_pattern):
            matched_files = [args.input_pattern]
        else:
            print(f"Error: No files found matching pattern: '{args.input_pattern}'")
            sys.exit(1)

    print("=" * 60)
    print("SND@LHC Event Display Generator (No Filtering)")
    print(f"Input Pattern     : {args.input_pattern}")
    print(f"Files Found       : {len(matched_files)}")
    print(f"Max Displays/File : {'All Events' if args.max_displays == 0 else args.max_displays}")
    print(f"Show All MCTracks : {args.show_all_mctracks}")
    print("=" * 60)

    overall_t0 = time.time()
    grand_total_events = 0
    grand_total_displays = 0
    all_sec_events = []
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

        print(f"\n[{idx}/{len(matched_files)}] Generating Displays: '{input_file}' -> '{out_file}'")
        t_evts, n_disp, sec_list = process_single_file(
            input_file, out_file, geo_elements, ggeo, scifi_module, snd_geo,
            max_displays=args.max_displays, show_all_mctracks=args.show_all_mctracks
        )

        grand_total_events += t_evts
        grand_total_displays += n_disp
        all_sec_events.extend(sec_list)
        created_files.append(out_file)

    overall_elapsed = time.time() - overall_t0

    print("\n" + "=" * 60)
    print("EVENT DISPLAY GENERATION SUMMARY")
    print("=" * 60)
    print(f"Total Input Files Processed : {len(matched_files)}")
    print(f"Total Output Files Created   : {len(created_files)}")
    print(f"Total Events Processed      : {grand_total_events:,}")
    print(f"Total Displays Generated    : {grand_total_displays:,}")
    print(f"Secondary-Induced Events    : {len(all_sec_events):,}")
    print(f"Total Execution Time        : {overall_elapsed:.1f} s")
    print("=" * 60)

    if all_sec_events:
        print("\nList of Secondary-Muon Induced Trident Events:")
        for evt_num, i_evt, zv, p_idx, s_idx, out_f in all_sec_events:
            print(f"  - File: {out_f} | Event #{evt_num} (Entry {i_evt}) | Primary Track: {p_idx}, Inducing Track: {s_idx} | Vertex Z = {zv:.1f} cm")
        print("=" * 60)

if __name__ == "__main__":
    main()
