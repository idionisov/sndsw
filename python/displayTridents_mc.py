import ROOT
import array
import argparse
import os
import numpy as np

DEFAULT_INPUT = "trimuon_events_0.root"
DEFAULT_GEOFILE = "python/geofile_full.Ntuple-TGeant4_boost100.0.root"
DEFAULT_OUTPUT = "trident_trajectories_root.root"
DEFAULT_MAX_EVENTS = 10
Z_3D_MIN = 0.0
Z_MAX = 600.0
ROCK_BOUNDARY = 260.0

ROOT.gROOT.SetBatch(True)

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
    If the track deposited MC points in ScifiPoint / MuFilterPoint / EmulsionDetPoint,
    uses the true simulated step coordinates in the detector to account for multiple
    Coulomb scattering through the rock. Otherwise falls back to straight-line extrapolation.
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

def get_trajectory_endpoints(track, z_start, z_end):
    """
    Calculates the (x, y, z) endpoints for a straight-line trajectory.
    Returns (x_start, x_end), (y_start, y_end), (z_start, z_end)
    """
    x0, y0, z0 = track.GetStartX(), track.GetStartY(), track.GetStartZ()
    px, py, pz = track.GetPx(), track.GetPy(), track.GetPz()

    if pz == 0:
        return (x0, x0), (y0, y0), (z_start, z_end)

    x_start = x0 + (px / pz) * (z_start - z0)
    y_start = y0 + (py / pz) * (z_start - z0)

    x_end = x0 + (px / pz) * (z_end - z0)
    y_end = y0 + (py / pz) * (z_end - z0)

    return (x_start, x_end), (y_start, y_end), (z_start, z_end)

def make_array(lst):
    return array.array('d', lst)

def style_detector_volumes(vol):
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
            dvol.SetLineColor(ROOT.kBlue+2)
            dvol.SetFillColor(ROOT.kAzure+7)
            dvol.SetTransparency(25)
        elif 'Veto' in name or 'subVeto' in vname:
            dvol.SetLineColor(ROOT.kOrange+2)
            dvol.SetFillColor(ROOT.kOrange+7)
            dvol.SetTransparency(30)
        elif 'FeBlock' in name:
            dvol.SetLineColor(ROOT.kGreen+4)
            dvol.SetFillColor(ROOT.kGreen+3)
            dvol.SetTransparency(45)
        elif 'subUSBox' in name or 'subDSBox' in name:
            dvol.SetLineColor(ROOT.kBlue+2)
            dvol.SetFillColor(ROOT.kBlue-4)
            dvol.SetTransparency(30)

        if dvol.GetNdaughters() > 0:
            style_detector_volumes(dvol)

def load_detector_geometry(geofile_path):
    if not os.path.exists(geofile_path):
        alt_path = os.path.join(os.path.dirname(__file__), os.path.basename(geofile_path))
        if os.path.exists(alt_path):
            geofile_path = alt_path
        else:
            print(f"Warning: Geofile '{geofile_path}' not found. Geometry will not be plotted.")
            return [], None, None, None

    snd_geo = None
    scifi_module = None
    try:
        import SndlhcGeo
        snd_geo = SndlhcGeo.GeoInterface(geofile_path)
        geo = snd_geo.sGeo
        f_geo = snd_geo.fgeo
        if hasattr(snd_geo, 'modules') and 'Scifi' in snd_geo.modules:
            scifi_module = snd_geo.modules['Scifi']
    except Exception as e:
        print(f"Warning: Could not initialize SndlhcGeo: {e}")
        f_geo = ROOT.TFile.Open(geofile_path, "READ")
        if not f_geo or f_geo.IsZombie():
            print(f"Error opening geofile '{geofile_path}'")
            return [], None, None, None
        geo = f_geo.Get("FAIRGeom")
        if not geo:
            print(f"Error: FAIRGeom TGeoManager object not found in '{geofile_path}'")
            return [], None, None, None

    ROOT.SetOwnership(f_geo, False)
    ROOT.SetOwnership(geo, False)
    ROOT.gGeoManager = geo
    elements = []

    def process_subsystem(parent_path):
        if not geo.cd(parent_path):
            return
        pnode = geo.GetCurrentNode()
        pvol = pnode.GetVolume()
        for i in range(pvol.GetNdaughters()):
            node = pvol.GetNode(i)
            name = node.GetName()
            subpath = f"{parent_path}/{name}"

            label, color, alpha = None, None, 0.2
            if 'Wall' in name and 'border' not in name:
                label, color, alpha = 'Target Emulsion Wall (Grey)', ROOT.kGray+1, 0.20
            elif 'ScifiVolume' in name:
                label, color, alpha = 'SciFi Station (Blue)', ROOT.kAzure+7, 0.25
            elif 'subVetoBox' in name or 'volVetoPlane' in name:
                label, color, alpha = 'Veto Detector (Orange)', ROOT.kOrange+7, 0.20
            elif 'FeBlock' in name:
                label, color, alpha = 'MuFilter Iron Block (Dark Green)', ROOT.kGreen+3, 0.15
            elif 'subUSBox' in name or 'subDSBox' in name:
                label, color, alpha = 'MuFilter Active Plane (Blue)', ROOT.kBlue-4, 0.25

            if label is not None:
                geo.cd(subpath)
                matrix = geo.GetCurrentMatrix()
                vol = node.GetVolume()
                shape = vol.GetShape()
                try:
                    dx, dy, dz = shape.GetDX(), shape.GetDY(), shape.GetDZ()
                    orig = shape.GetOrigin()
                except AttributeError:
                    continue
                master_corners = []
                for xs in [-1, 1]:
                    for ys in [-1, 1]:
                        for zs in [-1, 1]:
                            local = array.array('d', [orig[0] + xs*dx, orig[1] + ys*dy, orig[2] + zs*dz])
                            master = array.array('d', [0.0, 0.0, 0.0])
                            matrix.LocalToMaster(local, master)
                            master_corners.append(list(master))
                mc = np.array(master_corners)
                elements.append({
                    'name': name,
                    'label': label,
                    'color': color,
                    'alpha': alpha,
                    'z': (float(np.min(mc[:, 2])), float(np.max(mc[:, 2]))),
                    'x': (float(np.min(mc[:, 0])), float(np.max(mc[:, 0]))),
                    'y': (float(np.min(mc[:, 1])), float(np.max(mc[:, 1])))
                })

    process_subsystem('/cave_1/Detector_0/volVeto_1')
    process_subsystem('/cave_1/Detector_0/volTarget_1')
    process_subsystem('/cave_1/Detector_0/volMuFilter_1')

    top_det = geo.GetVolume("Detector")
    if top_det:
        ROOT.SetOwnership(top_det, False)
        style_detector_volumes(top_det)

    print(f"Loaded {len(elements)} detector geometry elements from '{geofile_path}'.")
    return elements, geo, f_geo, scifi_module

def make_3d_box(z_min, z_max, x_min, x_max, y_min, y_max, color, style=1, width=1):
    pts = [
        (z_min, x_min, y_min), (z_max, x_min, y_min), (z_max, x_max, y_min), (z_min, x_max, y_min), (z_min, x_min, y_min),
        (z_min, x_min, y_max), (z_max, x_min, y_max), (z_max, x_max, y_max), (z_min, x_max, y_max), (z_min, x_min, y_max),
        (z_max, x_min, y_max), (z_max, x_min, y_min),
        (z_max, x_max, y_min), (z_max, x_max, y_max),
        (z_min, x_max, y_max), (z_min, x_max, y_min)
    ]
    pl = ROOT.TPolyLine3D(len(pts))
    ROOT.SetOwnership(pl, False)
    for i, (z, x, y) in enumerate(pts):
        pl.SetPoint(i, z, x, y)
    pl.SetLineColor(color)
    pl.SetLineStyle(style)
    pl.SetLineWidth(width)
    return pl

def draw_detector_geometry(geo_elements, z_range, x_range, y_range, ggeo=None, file_tag="", evt_num=0, i_event=0):
    """
    Creates XZ, YZ, and 3D canvases and draws detector geometry volumes and boundaries.
    Returns a dictionary holding the canvases, axes, legends, and graphics objects.
    """
    z_min_plot, z_max_plot = z_range
    x_min_plot, x_max_plot = x_range
    y_min_plot, y_max_plot = y_range

    # ----------------- XZ Projection -----------------
    c_xz = ROOT.TCanvas(f"XZ_Ev_{evt_num}_Idx_{i_event}", f"XZ Projection (Event #{evt_num})", 950, 650)
    ROOT.SetOwnership(c_xz, False)
    mg_xz = ROOT.TMultiGraph()
    ROOT.SetOwnership(mg_xz, False)
    mg_xz.SetTitle(";Z [cm];X [cm]")

    gr_dummy_xz = ROOT.TGraph(2)
    ROOT.SetOwnership(gr_dummy_xz, False)
    gr_dummy_xz.SetPoint(0, z_min_plot, x_min_plot)
    gr_dummy_xz.SetPoint(1, z_max_plot, x_max_plot)
    mg_xz.Add(gr_dummy_xz)
    mg_xz.Draw("A")
    mg_xz.GetXaxis().SetLimits(z_min_plot, z_max_plot)
    mg_xz.GetYaxis().SetRangeUser(x_min_plot, x_max_plot)
    c_xz.Update()

    boxes_xz = []
    drawn_labels_xz = set()
    legend_xz = ROOT.TLegend(0.12, 0.65, 0.52, 0.88)
    ROOT.SetOwnership(legend_xz, False)
    legend_xz.SetBorderSize(1)
    legend_xz.SetFillStyle(1001)

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

        if elem['label'] not in drawn_labels_xz:
            legend_xz.AddEntry(box, elem['label'], "f")
            drawn_labels_xz.add(elem['label'])

    line_xz = ROOT.TLine(ROCK_BOUNDARY, x_min_plot, ROCK_BOUNDARY, x_max_plot)
    ROOT.SetOwnership(line_xz, False)
    line_xz.SetLineColor(ROOT.kGray+2)
    line_xz.SetLineStyle(3)
    line_xz.Draw()

    # ----------------- YZ Projection -----------------
    c_yz = ROOT.TCanvas(f"YZ_Ev_{evt_num}_Idx_{i_event}", f"YZ Projection (Event #{evt_num})", 950, 650)
    ROOT.SetOwnership(c_yz, False)
    mg_yz = ROOT.TMultiGraph()
    ROOT.SetOwnership(mg_yz, False)
    mg_yz.SetTitle(";Z [cm];Y [cm]")

    gr_dummy_yz = ROOT.TGraph(2)
    ROOT.SetOwnership(gr_dummy_yz, False)
    gr_dummy_yz.SetPoint(0, z_min_plot, y_min_plot)
    gr_dummy_yz.SetPoint(1, z_max_plot, y_max_plot)
    mg_yz.Add(gr_dummy_yz)
    mg_yz.Draw("A")
    mg_yz.GetXaxis().SetLimits(z_min_plot, z_max_plot)
    mg_yz.GetYaxis().SetRangeUser(y_min_plot, y_max_plot)
    c_yz.Update()

    boxes_yz = []
    drawn_labels_yz = set()
    legend_yz = ROOT.TLegend(0.12, 0.65, 0.52, 0.88)
    ROOT.SetOwnership(legend_yz, False)
    legend_yz.SetBorderSize(1)
    legend_yz.SetFillStyle(1001)

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

        if elem['label'] not in drawn_labels_yz:
            legend_yz.AddEntry(box, elem['label'], "f")
            drawn_labels_yz.add(elem['label'])

    line_yz = ROOT.TLine(ROCK_BOUNDARY, y_min_plot, ROCK_BOUNDARY, y_max_plot)
    ROOT.SetOwnership(line_yz, False)
    line_yz.SetLineColor(ROOT.kGray+2)
    line_yz.SetLineStyle(3)
    line_yz.Draw()

    # ----------------- 3D View -----------------
    c_3d = ROOT.TCanvas(f"3D_Ev_{evt_num}_Idx_{i_event}", f"3D View (Event #{evt_num})", 850, 850)
    ROOT.SetOwnership(c_3d, False)

    h3 = ROOT.TH3F(f"box_{evt_num}_Idx_{i_event}", ";Z [cm];X [cm];Y [cm]",
                   10, Z_3D_MIN, Z_MAX,
                   10, x_min_plot, x_max_plot,
                   10, y_min_plot, y_max_plot)
    ROOT.SetOwnership(h3, False)
    h3.SetStats(0)
    h3.Draw("AXIS")

    if ggeo and ggeo.GetVolume("Detector"):
        top_det_vol = ggeo.GetVolume("Detector")
        ROOT.SetOwnership(top_det_vol, False)
        top_det_vol.Draw("same")

    geo_3d_lines = []
    for elem in geo_elements:
        box_3d = make_3d_box(
            elem['z'][0], elem['z'][1],
            elem['x'][0], elem['x'][1],
            elem['y'][0], elem['y'][1],
            elem['color'], style=1, width=1
        )
        box_3d.Draw("same")
        geo_3d_lines.append(box_3d)

    plane_3d = ROOT.TPolyLine3D(5)
    ROOT.SetOwnership(plane_3d, False)
    plane_3d.SetPoint(0, ROCK_BOUNDARY, x_min_plot, y_min_plot)
    plane_3d.SetPoint(1, ROCK_BOUNDARY, x_max_plot, y_min_plot)
    plane_3d.SetPoint(2, ROCK_BOUNDARY, x_max_plot, y_max_plot)
    plane_3d.SetPoint(3, ROCK_BOUNDARY, x_min_plot, y_max_plot)
    plane_3d.SetPoint(4, ROCK_BOUNDARY, x_min_plot, y_min_plot)
    plane_3d.SetLineColor(ROOT.kGray+1)
    plane_3d.SetLineStyle(3)
    plane_3d.Draw("same")

    return {
        'c_xz': c_xz, 'c_yz': c_yz, 'c_3d': c_3d,
        'mg_xz': mg_xz, 'mg_yz': mg_yz, 'h3': h3,
        'legend_xz': legend_xz, 'legend_yz': legend_yz,
        'boxes_xz': boxes_xz, 'boxes_yz': boxes_yz,
        'line_xz': line_xz, 'line_yz': line_yz,
        'geo_3d_lines': geo_3d_lines, 'plane_3d': plane_3d
    }

def draw_scifi_hits(event, scifi_module, geo_ctx, marker_color=ROOT.kGreen+2, marker_style=20, marker_size=0.8):
    """
    Draws Digi_ScifiHits from event onto the XZ, YZ, and 3D canvases.
    - Vertical fibers are plotted as centroid markers in XZ view.
    - Horizontal fibers are plotted as centroid markers in YZ view.
    - 3D fiber segments from A to B are plotted on 3D view.
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

    xz_pts = []
    yz_pts = []
    poly3d_list = []

    for hit in event.Digi_ScifiHits:
        detID = hit.GetDetectorID()
        scifi_module.GetSiPMPosition(detID, A, B)

        z_mid = (A.Z() + B.Z()) / 2.0
        x_mid = (A.X() + B.X()) / 2.0
        y_mid = (A.Y() + B.Y()) / 2.0

        if hit.isVertical():
            xz_pts.append((z_mid, x_mid))
        else:
            yz_pts.append((z_mid, y_mid))

        pl3d = ROOT.TPolyLine3D(2)
        ROOT.SetOwnership(pl3d, False)
        pl3d.SetPoint(0, A.Z(), A.X(), A.Y())
        pl3d.SetPoint(1, B.Z(), B.X(), B.Y())
        pl3d.SetLineColor(marker_color)
        pl3d.SetLineWidth(2)
        poly3d_list.append(pl3d)

    gr_xz = None
    if xz_pts:
        gr_xz = ROOT.TGraph(len(xz_pts))
        ROOT.SetOwnership(gr_xz, False)
        for idx, (z, x) in enumerate(xz_pts):
            gr_xz.SetPoint(idx, z, x)
        gr_xz.SetMarkerStyle(marker_style)
        gr_xz.SetMarkerSize(marker_size)
        gr_xz.SetMarkerColor(marker_color)
        c_xz.cd()
        gr_xz.Draw("P same")
        legend_xz.AddEntry(gr_xz, "SciFi Hits (Digi_ScifiHits)", "p")

    gr_yz = None
    if yz_pts:
        gr_yz = ROOT.TGraph(len(yz_pts))
        ROOT.SetOwnership(gr_yz, False)
        for idx, (z, y) in enumerate(yz_pts):
            gr_yz.SetPoint(idx, z, y)
        gr_yz.SetMarkerStyle(marker_style)
        gr_yz.SetMarkerSize(marker_size)
        gr_yz.SetMarkerColor(marker_color)
        c_yz.cd()
        gr_yz.Draw("P same")
        legend_yz.AddEntry(gr_yz, "SciFi Hits (Digi_ScifiHits)", "p")

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

    x_p, y_p, z_p = get_track_trajectory_points(primary_idx, primary_track, z_min_2d_plot, Z_MAX, mc_points_map)

    x_ind, y_ind, z_ind = None, None, None
    if is_secondary_inducing:
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
    if is_secondary_inducing and x_ind is not None:
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

def create_event_display(event, evt_num, i_event, geo_elements, ggeo, scifi_module,
                         primary_idx, primary_track, inducing_idx, inducing_track, my_daughters, file_tag="",
                         show_all_mctracks=False, draw_hits=True, draw_tracks=True):
    """
    High-level orchestrator that creates canvases, calls draw_detector_geometry,
    draw_scifi_hits, and draw_mctracks, draws legends/banners, and returns the canvases.
    """
    is_secondary_inducing = (inducing_track is not None and primary_track is not None and inducing_track != primary_track)

    mc_points_map = build_mc_points_map(event)

    # Calculate preliminary plot bounding box
    z_min_2d = 250.0  # Show detector region starting before Veto
    x_p, y_p, z_p = get_track_trajectory_points(primary_idx, primary_track, z_min_2d, Z_MAX, mc_points_map)

    x_ind, y_ind, z_ind = None, None, None
    if is_secondary_inducing:
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

    # 2. Draw SciFi Hits
    scifi_ctx = {}
    if draw_hits and scifi_module:
        scifi_ctx = draw_scifi_hits(event, scifi_module, geo_ctx)

    # 3. Draw MCTracks
    mctrack_ctx = {}
    if draw_tracks:
        mctrack_ctx = draw_mctracks(
            event, primary_idx, primary_track, inducing_idx, inducing_track, my_daughters,
            geo_ctx, z_range=(z_min_2d_plot, z_max_2d_plot),
            show_all_mctracks=show_all_mctracks, mc_points_map=mc_points_map
        )

    # 4. Draw Legends and Banners
    c_xz = geo_ctx['c_xz']
    c_yz = geo_ctx['c_yz']
    c_3d = geo_ctx['c_3d']

    c_xz.cd()
    geo_ctx['legend_xz'].Draw()
    banner_xz = ROOT.TLatex()
    ROOT.SetOwnership(banner_xz, False)
    banner_xz.SetNDC()
    banner_xz.SetTextFont(42)
    banner_xz.SetTextSize(0.035)
    banner_xz.DrawLatex(0.12, 0.92, banner_text)

    c_yz.cd()
    geo_ctx['legend_yz'].Draw()
    banner_yz = ROOT.TLatex()
    ROOT.SetOwnership(banner_yz, False)
    banner_yz.SetNDC()
    banner_yz.SetTextFont(42)
    banner_yz.SetTextSize(0.035)
    banner_yz.DrawLatex(0.12, 0.92, banner_text)

    c_3d.cd()
    banner_3d = ROOT.TLatex()
    ROOT.SetOwnership(banner_3d, False)
    banner_3d.SetNDC()
    banner_3d.SetTextFont(42)
    banner_3d.SetTextSize(0.035)
    banner_3d.DrawLatex(0.12, 0.92, banner_text)

    return c_xz, c_yz, c_3d, evt_dir_name, is_secondary_inducing, z_vert_val

def main():
    parser = argparse.ArgumentParser(description="Display Trident Muon Trajectories with Detector Geometry")
    parser.add_argument("-i", "--input", dest="input_filename", default=DEFAULT_INPUT, help="Input ROOT event file")
    parser.add_argument("-g", "--geofile", dest="geofile_filename", default=DEFAULT_GEOFILE, help="ROOT geometry file")
    parser.add_argument("-o", "--output", dest="output_filename", default=DEFAULT_OUTPUT, help="Output ROOT file")
    parser.add_argument("-n", "--max_events", dest="max_events", type=int, default=DEFAULT_MAX_EVENTS, help="Max events to plot")
    parser.add_argument("-a", "--all_mctracks", "--show_all_mctracks", dest="show_all_mctracks", action="store_true", default=False, help="Include all other MCTracks as thin, transparent grey lines")
    parser.add_argument("-v", "--view", dest="interactive_view", action="store_true", help="Open 3D WebGL HTML plot in default web browser")
    args = parser.parse_args()

    geo_elements, ggeo, f_geo, scifi_module = load_detector_geometry(args.geofile_filename)

    try:
        f_in = ROOT.TFile.Open(args.input_filename, "READ")
        if not f_in or f_in.IsZombie():
            print(f"Error: Could not open input file {args.input_filename}")
            return
        ROOT.SetOwnership(f_in, False)
    except Exception as e:
        print(f"Error: Could not open input file '{args.input_filename}': {e}")
        return

    tree = f_in.Get("cbmsim")
    if not tree:
        print(f"Error: 'cbmsim' TTree not found in {args.input_filename}")
        return

    f_out = ROOT.TFile.Open(args.output_filename, "RECREATE")
    if not f_out or f_out.IsZombie():
        print(f"Error: Could not create output file {args.output_filename}")
        return
    ROOT.SetOwnership(f_out, False)

    events_plotted = 0

    for i_event, event in enumerate(tree):
        if events_plotted >= args.max_events:
            break

        # Extract Event Number from EventHeader (or MCEventHeader / fallback to index)
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

        inducing_idx = -1
        my_daughters = []

        for m_idx, m_tr in incident_muons.items():
            ds = [d for d in pair_daughters if d[1] == m_idx]
            has_m = any(d[2] == 13 for d in ds)
            has_p = any(d[2] == -13 for d in ds)
            zv_ok = all(d[3].GetStartZ() < ROCK_BOUNDARY for d in ds) if ds else False
            if has_m and has_p and zv_ok:
                inducing_idx = m_idx
                my_daughters = ds
                break

        if (primary_idx != -1 or inducing_idx != -1) and len(my_daughters) >= 2:
            inducing_track = incident_muons[inducing_idx] if inducing_idx != -1 else None
            if primary_idx != -1:
                primary_track = event.MCTrack[primary_idx]
            else:
                curr_tr = inducing_track
                while curr_tr.GetMotherId() != -1 and curr_tr.GetMotherId() in incident_muons:
                    curr_tr = incident_muons[curr_tr.GetMotherId()]
                primary_track = curr_tr

            c_xz, c_yz, c_3d, evt_dir_name, is_sec_ind, z_vert_val = create_event_display(
                event, evt_num, i_event, geo_elements, ggeo, scifi_module,
                primary_idx, primary_track, inducing_idx, inducing_track, my_daughters, file_tag="",
                show_all_mctracks=args.show_all_mctracks, draw_hits=True, draw_tracks=True
            )

            c_xz.Write(f"XZ_Projection_Ev_{evt_num}")
            c_yz.Write(f"YZ_Projection_Ev_{evt_num}")
            c_3d.Write(f"3D_View_Ev_{evt_num}")
            events_plotted += 1
            print(f"Saved Event {events_plotted} (TTree Index: {i_event}) to ROOT file.")

    f_out.Close()
    f_in.Close()
    print(f"\nDone! All canvases saved in '{args.output_filename}'.")

if __name__ == "__main__":
    main()
