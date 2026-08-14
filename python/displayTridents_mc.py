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
ROCK_BOUNDARY = 360.0

ROOT.gROOT.SetBatch(True)

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
            return [], None, None

    f_geo = ROOT.TFile.Open(geofile_path, "READ")
    if not f_geo or f_geo.IsZombie():
        print(f"Error opening geofile '{geofile_path}'")
        return [], None, None

    ROOT.SetOwnership(f_geo, False)
    geo = f_geo.Get("FAIRGeom")
    if not geo:
        print(f"Error: FAIRGeom TGeoManager object not found in '{geofile_path}'")
        return [], None, None

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

    # Apply solid fill colors and transparency to Detector TGeoVolume
    top_det = geo.GetVolume("Detector")
    if top_det:
        ROOT.SetOwnership(top_det, False)
        style_detector_volumes(top_det)

    print(f"Loaded {len(elements)} detector geometry elements from '{geofile_path}'.")
    return elements, geo, f_geo

def make_3d_box(z_min, z_max, x_min, x_max, y_min, y_max, color, style=1, width=1):
    """
    Creates a TPolyLine3D representing the 12 wireframe edges of a 3D bounding box.
    """
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

def main():
    parser = argparse.ArgumentParser(description="Display Trident Muon Trajectories with Detector Geometry")
    parser.add_argument("-i", "--input", dest="input_filename", default=DEFAULT_INPUT, help="Input ROOT event file")
    parser.add_argument("-g", "--geofile", dest="geofile_filename", default=DEFAULT_GEOFILE, help="ROOT geometry file")
    parser.add_argument("-o", "--output", dest="output_filename", default=DEFAULT_OUTPUT, help="Output ROOT file")
    parser.add_argument("-n", "--max_events", dest="max_events", type=int, default=DEFAULT_MAX_EVENTS, help="Max events to plot")
    parser.add_argument("-v", "--view", dest="interactive_view", action="store_true", help="Open 3D WebGL HTML plot in default web browser")
    args = parser.parse_args()

    geo_elements, ggeo, f_geo = load_detector_geometry(args.geofile_filename)

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
            inducing_track = incident_muons[inducing_idx]
            if primary_idx != -1:
                primary_track = event.MCTrack[primary_idx]
            else:
                curr_tr = inducing_track
                while curr_tr.GetMotherId() != -1 and curr_tr.GetMotherId() in incident_muons:
                    curr_tr = incident_muons[curr_tr.GetMotherId()]
                primary_track = curr_tr

            z_vert_val = my_daughters[0][3].GetStartZ() if my_daughters else Z_MAX
            is_secondary_inducing = (inducing_idx != primary_idx and primary_idx != -1)

            z_min_2d = primary_track.GetStartZ()
            x_p, y_p, z_p = get_trajectory_endpoints(primary_track, z_min_2d, Z_MAX)

            x_ind, y_ind, z_ind = None, None, None
            if is_secondary_inducing:
                ind_z_start = inducing_track.GetStartZ()
                x_ind, y_ind, z_ind = get_trajectory_endpoints(inducing_track, ind_z_start, Z_MAX)

            d_trajs = []
            for d in my_daughters:
                d_pdg = d[2]
                d_track = d[3]
                d_z_start = d_track.GetStartZ()
                x_d, y_d, z_d = get_trajectory_endpoints(d_track, d_z_start, Z_MAX)
                d_trajs.append((d_pdg, x_d, y_d, z_d))

            all_x = list(x_p) + [x for d in d_trajs for x in d[1]] + [elem['x'][0] for elem in geo_elements] + [elem['x'][1] for elem in geo_elements]
            all_y = list(y_p) + [y for d in d_trajs for y in d[2]] + [elem['y'][0] for elem in geo_elements] + [elem['y'][1] for elem in geo_elements]
            all_z = list(z_p) + [z for d in d_trajs for z in d[3]] + [elem['z'][0] for elem in geo_elements] + [elem['z'][1] for elem in geo_elements]
            if is_secondary_inducing:
                all_x += list(x_ind)
                all_y += list(y_ind)
                all_z += list(z_ind)

            x_min_plot, x_max_plot = min(all_x) - 5, max(all_x) + 5
            y_min_plot, y_max_plot = min(all_y) - 5, max(all_y) + 5
            z_min_2d_plot, z_max_2d_plot = min(all_z) - 5, max(all_z) + 5

            banner_text = f"Event #{evt_num}  (Entry {i_event})"
            if is_secondary_inducing:
                banner_text += "   #color[616]{[Secondary Muon Induced]}"

            # ----------------- XZ Projection -----------------
            c_xz = ROOT.TCanvas(f"XZ_Ev_{evt_num}_Idx_{i_event}", f"XZ Projection (Event #{evt_num})", 950, 650)
            ROOT.SetOwnership(c_xz, False)
            mg_xz = ROOT.TMultiGraph()
            ROOT.SetOwnership(mg_xz, False)
            mg_xz.SetTitle(";Z [cm];X [cm]")

            gr_p_xz = ROOT.TGraph(2, make_array(z_p), make_array(x_p))
            ROOT.SetOwnership(gr_p_xz, False)
            gr_p_xz.SetLineColor(ROOT.kBlack)
            gr_p_xz.SetLineWidth(2)
            mg_xz.Add(gr_p_xz)

            gr_ind_xz = None
            if is_secondary_inducing:
                gr_ind_xz = ROOT.TGraph(2, make_array(z_ind), make_array(x_ind))
                ROOT.SetOwnership(gr_ind_xz, False)
                gr_ind_xz.SetLineColor(ROOT.kMagenta+2)
                gr_ind_xz.SetLineWidth(2)
                mg_xz.Add(gr_ind_xz)

            gr_ds_xz = []
            for d_pdg, x_d, y_d, z_d in d_trajs:
                gr = ROOT.TGraph(2, make_array(z_d), make_array(x_d))
                ROOT.SetOwnership(gr, False)
                gr.SetLineWidth(2)
                gr.SetLineStyle(2)
                gr.SetLineColor(ROOT.kBlue if d_pdg == 13 else ROOT.kRed)
                mg_xz.Add(gr)
                gr_ds_xz.append(gr)

            mg_xz.Draw("A")
            mg_xz.GetXaxis().SetLimits(z_min_2d_plot, z_max_2d_plot)
            mg_xz.GetYaxis().SetRangeUser(x_min_plot, x_max_plot)
            c_xz.Update()

            boxes_xz = []
            drawn_labels_xz = set()
            legend_xz = ROOT.TLegend(0.12, 0.65, 0.52, 0.88)
            ROOT.SetOwnership(legend_xz, False)
            legend_xz.SetBorderSize(1)
            legend_xz.SetFillStyle(1001)

            legend_xz.AddEntry(gr_p_xz, "Primary #mu^{-}", "l")
            if is_secondary_inducing and gr_ind_xz:
                legend_xz.AddEntry(gr_ind_xz, "Inducing Secondary #mu", "l")
            if gr_ds_xz:
                legend_xz.AddEntry(gr_ds_xz[0], "Trident #mu / #mu^{+} daughters", "l")

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

            gr_p_xz.Draw("L same")
            if is_secondary_inducing and gr_ind_xz:
                gr_ind_xz.Draw("L same")
            for gr in gr_ds_xz:
                gr.Draw("L same")

            line_xz = ROOT.TLine(ROCK_BOUNDARY, x_min_plot, ROCK_BOUNDARY, x_max_plot)
            ROOT.SetOwnership(line_xz, False)
            line_xz.SetLineColor(ROOT.kGray+2)
            line_xz.SetLineStyle(3)
            line_xz.Draw()

            legend_xz.Draw()

            banner_xz = ROOT.TLatex()
            ROOT.SetOwnership(banner_xz, False)
            banner_xz.SetNDC()
            banner_xz.SetTextFont(42)
            banner_xz.SetTextSize(0.035)
            banner_xz.DrawLatex(0.12, 0.92, banner_text)

            c_xz.Write("XZ_Projection")

            # ----------------- YZ Projection -----------------
            c_yz = ROOT.TCanvas(f"YZ_Ev_{evt_num}_Idx_{i_event}", f"YZ Projection (Event #{evt_num})", 950, 650)
            ROOT.SetOwnership(c_yz, False)
            mg_yz = ROOT.TMultiGraph()
            ROOT.SetOwnership(mg_yz, False)
            mg_yz.SetTitle(";Z [cm];Y [cm]")

            gr_p_yz = ROOT.TGraph(2, make_array(z_p), make_array(y_p))
            ROOT.SetOwnership(gr_p_yz, False)
            gr_p_yz.SetLineColor(ROOT.kBlack)
            gr_p_yz.SetLineWidth(2)
            mg_yz.Add(gr_p_yz)

            gr_ind_yz = None
            if is_secondary_inducing:
                gr_ind_yz = ROOT.TGraph(2, make_array(z_ind), make_array(y_ind))
                ROOT.SetOwnership(gr_ind_yz, False)
                gr_ind_yz.SetLineColor(ROOT.kMagenta+2)
                gr_ind_yz.SetLineWidth(2)
                mg_yz.Add(gr_ind_yz)

            gr_ds_yz = []
            for d_pdg, x_d, y_d, z_d in d_trajs:
                gr = ROOT.TGraph(2, make_array(z_d), make_array(y_d))
                ROOT.SetOwnership(gr, False)
                gr.SetLineWidth(2)
                gr.SetLineStyle(2)
                gr.SetLineColor(ROOT.kBlue if d_pdg == 13 else ROOT.kRed)
                mg_yz.Add(gr)
                gr_ds_yz.append(gr)

            mg_yz.Draw("A")
            mg_yz.GetXaxis().SetLimits(z_min_2d_plot, z_max_2d_plot)
            mg_yz.GetYaxis().SetRangeUser(y_min_plot, y_max_plot)
            c_yz.Update()

            boxes_yz = []
            drawn_labels_yz = set()
            legend_yz = ROOT.TLegend(0.12, 0.65, 0.52, 0.88)
            ROOT.SetOwnership(legend_yz, False)
            legend_yz.SetBorderSize(1)
            legend_yz.SetFillStyle(1001)

            legend_yz.AddEntry(gr_p_yz, "Primary #mu^{-}", "l")
            if is_secondary_inducing and gr_ind_yz:
                legend_yz.AddEntry(gr_ind_yz, "Inducing Secondary #mu", "l")
            if gr_ds_yz:
                legend_yz.AddEntry(gr_ds_yz[0], "Trident #mu / #mu^{+} daughters", "l")

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

            gr_p_yz.Draw("L same")
            if is_secondary_inducing and gr_ind_yz:
                gr_ind_yz.Draw("L same")
            for gr in gr_ds_yz:
                gr.Draw("L same")

            line_yz = ROOT.TLine(ROCK_BOUNDARY, y_min_plot, ROCK_BOUNDARY, y_max_plot)
            ROOT.SetOwnership(line_yz, False)
            line_yz.SetLineColor(ROOT.kGray+2)
            line_yz.SetLineStyle(3)
            line_yz.Draw()

            legend_yz.Draw()

            banner_yz = ROOT.TLatex()
            ROOT.SetOwnership(banner_yz, False)
            banner_yz.SetNDC()
            banner_yz.SetTextFont(42)
            banner_yz.SetTextSize(0.035)
            banner_yz.DrawLatex(0.12, 0.92, banner_text)

            c_yz.Write("YZ_Projection")

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

            p_3d_start_z = max(primary_track.GetStartZ(), Z_3D_MIN)
            x_p_3d, y_p_3d, z_p_3d = get_trajectory_endpoints(primary_track, p_3d_start_z, Z_MAX)
            pl_p = ROOT.TPolyLine3D(2)
            ROOT.SetOwnership(pl_p, False)
            pl_p.SetPoint(0, z_p_3d[0], x_p_3d[0], y_p_3d[0])
            pl_p.SetPoint(1, z_p_3d[1], x_p_3d[1], y_p_3d[1])
            pl_p.SetLineColor(ROOT.kBlack)
            pl_p.SetLineWidth(3)
            pl_p.Draw("same")

            if is_secondary_inducing:
                ind_z_start_3d = max(inducing_track.GetStartZ(), Z_3D_MIN)
                x_ind_3d, y_ind_3d, z_ind_3d = get_trajectory_endpoints(inducing_track, ind_z_start_3d, Z_MAX)
                pl_ind = ROOT.TPolyLine3D(2)
                ROOT.SetOwnership(pl_ind, False)
                pl_ind.SetPoint(0, z_ind_3d[0], x_ind_3d[0], y_ind_3d[0])
                pl_ind.SetPoint(1, z_ind_3d[1], x_ind_3d[1], y_ind_3d[1])
                pl_ind.SetLineColor(ROOT.kMagenta+2)
                pl_ind.SetLineWidth(3)
                pl_ind.Draw("same")

            pls = []
            for d_pdg, x_d, y_d, z_d in d_trajs:
                d_z_start = max(z_d[0], Z_3D_MIN)
                pl = ROOT.TPolyLine3D(2)
                ROOT.SetOwnership(pl, False)
                pl.SetPoint(0, d_z_start, x_d[0], y_d[0])
                pl.SetPoint(1, z_d[1], x_d[1], y_d[1])
                pl.SetLineColor(ROOT.kBlue if d_pdg == 13 else ROOT.kRed)
                pl.SetLineStyle(2)
                pl.SetLineWidth(2)
                pl.Draw("same")
                pls.append(pl)

            plane = ROOT.TPolyLine3D(5)
            ROOT.SetOwnership(plane, False)
            plane.SetPoint(0, ROCK_BOUNDARY, x_min_plot, y_min_plot)
            plane.SetPoint(1, ROCK_BOUNDARY, x_max_plot, y_min_plot)
            plane.SetPoint(2, ROCK_BOUNDARY, x_max_plot, y_max_plot)
            plane.SetPoint(3, ROCK_BOUNDARY, x_min_plot, y_max_plot)
            plane.SetPoint(4, ROCK_BOUNDARY, x_min_plot, y_min_plot)
            plane.SetLineColor(ROOT.kGray+1)
            plane.SetLineStyle(3)
            plane.Draw("same")

            banner_3d = ROOT.TLatex()
            ROOT.SetOwnership(banner_3d, False)
            banner_3d.SetNDC()
            banner_3d.SetTextFont(42)
            banner_3d.SetTextSize(0.035)
            banner_3d.DrawLatex(0.12, 0.92, f"#bf{{SND@LHC Trident}}  Event #{evt_num}  (Entry {i_event})")

            c_3d.Write("3D_View")
            events_plotted += 1
            print(f"Saved Event {events_plotted} (TTree Index: {i_event}) to ROOT file.")

    f_out.Close()
    f_in.Close()
    print(f"\nDone! All canvases saved in '{args.output_filename}'.")

if __name__ == "__main__":
    main()
