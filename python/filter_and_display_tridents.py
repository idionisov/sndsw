import os
import sys
import array
import argparse
import warnings
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)

try:
    ROOT.EnableImplicitMT(6)
except Exception:
    pass

sndsw_path = os.environ.get("SNDSW_ROOT", "/afs/cern.ch/user/i/idioniso/snd_master/sndsw")

ROOT.gInterpreter.ProcessLine('#include "sndScifiHit.h"')
ROOT.gInterpreter.ProcessLine('#include "ShipMCTrack.h"')

for tool_hdr in ["sndSciFiTools.h", "sndTchainGetter.h", "sndGeometryGetter.h"]:
    hdr_path = os.path.join(sndsw_path, "analysis", "tools", tool_hdr)
    if os.path.exists(hdr_path):
        ROOT.gInterpreter.ProcessLine(f'#include "{hdr_path}"')

ROOT.gInterpreter.Declare("""
#include "sndScifiHit.h"
#include "ShipMCTrack.h"
#include <atomic>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <vector>
#include <cmath>

std::atomic<long long> rdf_events{0};
std::atomic<long long> total_rdf_events{0};
std::chrono::steady_clock::time_point start_time;
std::mutex cout_mutex;

void reset_progress(long long total = 0) {
    rdf_events = 0;
    total_rdf_events = total;
    start_time = std::chrono::steady_clock::now();
}

bool print_progress() {
    long long c = ++rdf_events;
    long long total = total_rdf_events.load();

    if (c % 100000 == 0 || (total > 0 && c == total)) {
        auto now = std::chrono::steady_clock::now();
        auto elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();

        long long hours = elapsed_sec / 3600;
        long long minutes = (elapsed_sec % 3600) / 60;
        long long seconds = elapsed_sec % 60;

        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "Processed " << c;
        if (total > 0) {
            double pct = (100.0 * c) / total;
            std::cout << " / " << total << " (" << std::fixed << std::setprecision(1) << pct << "%)";
        }
        std::cout << " | Elapsed: "
                  << std::setfill('0') << std::setw(2) << hours << ":"
                  << std::setfill('0') << std::setw(2) << minutes << ":"
                  << std::setfill('0') << std::setw(2) << seconds
                  << std::endl;
    }
    return true;
}

const double ROCK_Z_BOUNDARY = 360.0;

struct TrimuonResult {
    bool is_rock_trimuon;
    double weight;
};

TrimuonResult check_rock_trimuon(const TClonesArray& mc_tracks) {
    std::vector<int> incident_muons;

    struct PairProducedMuon {
        int index;
        int mother_id;
        int pdg;
        double z_vert;
        double weight;
    };
    std::vector<PairProducedMuon> pair_produced_muons;

    int n_tracks = mc_tracks.GetEntriesFast();
    for (int i = 0; i < n_tracks; ++i) {
        auto* tr = static_cast<ShipMCTrack*>(mc_tracks.At(i));
        if (!tr) continue;

        int pdg = tr->GetPdgCode();
        int mother_id = tr->GetMotherId();
        int proc_id = tr->GetProcID();

        if (std::abs(pdg) == 13) {
            incident_muons.push_back(i);
        }
        if (std::abs(pdg) == 13 && proc_id == 5) {
            pair_produced_muons.push_back({i, mother_id, pdg, tr->GetStartZ(), tr->GetWeight()});
        }
    }

    for (int mu_idx : incident_muons) {
        bool has_mu_minus = false;
        bool has_mu_plus = false;
        bool all_zv_ok = true;
        double first_weight = 1.0;
        bool found_daughter = false;

        for (const auto& d : pair_produced_muons) {
            if (d.mother_id == mu_idx) {
                if (!found_daughter) {
                    first_weight = d.weight;
                    found_daughter = true;
                }
                if (d.pdg == 13) has_mu_minus = true;
                if (d.pdg == -13) has_mu_plus = true;
                if (d.z_vert >= ROCK_Z_BOUNDARY) {
                    all_zv_ok = false;
                }
            }
        }

        if (has_mu_minus && has_mu_plus && all_zv_ok) {
            return {true, first_weight};
        }
    }

    return {false, 1.0};
}

bool is_rock_trimuon_event(const TClonesArray& mc_tracks) {
    return check_rock_trimuon(mc_tracks).is_rock_trimuon;
}

double get_trimuon_weight(const TClonesArray& mc_tracks) {
    return check_rock_trimuon(mc_tracks).weight;
}
""")

class DataSet:
    def __init__(self, name, path_pattern, is_mc=True, tree_name="cbmsim"):
        self.name = name
        self.is_mc = is_mc
        self.tree_name = tree_name
        self.file_vec = ROOT.std.vector('std::string')()

        patterns = path_pattern if isinstance(path_pattern, (list, tuple)) else [path_pattern]
        for pat in patterns:
            self.add_pattern(pat)

    def add_pattern(self, pattern):
        if "*" in pattern or "?" in pattern:
            cmd = f"ls -1 {pattern} 2>/dev/null"
            raw_output = ROOT.gSystem.GetFromPipe(cmd).Data()
            matched_files = [f.strip() for f in raw_output.split("\n") if f.strip().endswith(".root")]
        else:
            matched_files = [pattern] if os.path.exists(pattern) else []

        if not matched_files:
            warnings.warn(f"DataSet '{self.name}': No files found matching pattern '{pattern}'.", UserWarning)

        for f in matched_files:
            self.file_vec.push_back(f)

    @property
    def n_files(self):
        return self.file_vec.size()

    def to_rdf(self):
        if self.file_vec.empty():
            raise RuntimeError(f"DataSet '{self.name}' has 0 matching files.")
        return ROOT.RDataFrame(self.tree_name, self.file_vec)

# Geometry calculation and plotting constants
Z_3D_MIN = 0.0
Z_MAX = 600.0
ROCK_BOUNDARY = 360.0

def get_trajectory_endpoints(track, z_start, z_end):
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

    top_det = geo.GetVolume("Detector")
    if top_det:
        ROOT.SetOwnership(top_det, False)
        style_detector_volumes(top_det)

    print(f"Loaded {len(elements)} detector geometry elements from '{geofile_path}'.")
    return elements, geo, f_geo

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

def get_output_filename(input_path, output_arg):
    base_name = os.path.basename(input_path)
    name_no_ext = os.path.splitext(base_name)[0]

    match = re.search(r'(digCPP-\d+|\d+)$', name_no_ext)
    tag = match.group(1) if match else name_no_ext

    if '%s' in output_arg:
        return output_arg % tag

    out_dir = os.path.dirname(output_arg)
    out_base = os.path.basename(output_arg)
    out_stem, out_ext = os.path.splitext(out_base)
    if not out_ext:
        out_ext = '.root'

    new_base = f"{out_stem}_{tag}{out_ext}"
    return os.path.join(out_dir, new_base) if out_dir else new_base

def process_single_file(input_file, output_file, geo_elements, ggeo, f_geo, max_displays=0, file_tag=""):
    print(f"\n------------------------------------------")
    print(f"Processing Input File : '{input_file}'")
    print(f"Output Destination    : '{output_file}'")
    print(f"------------------------------------------")

    is_prefiltered = False
    try:
        f_chk = ROOT.TFile.Open(input_file, "READ")
        if f_chk and not f_chk.IsZombie():
            t_chk = f_chk.Get("cbmsim")
            if t_chk and t_chk.GetBranch("event_weight"):
                is_prefiltered = True
            f_chk.Close()
    except Exception:
        pass

    if is_prefiltered:
        print(f"Input file is already pre-filtered.")
        if input_file != output_file:
            print(f"Copying '{input_file}' -> '{output_file}'...")
            os.system(f"cp {input_file} {output_file}")
        unweighted_selected = -1
        total_weighted_yield = 0.0
    else:
        rdf = ROOT.RDataFrame("cbmsim", input_file)
        total_events = rdf.Count().GetValue()

        col_names = [str(col) for col in rdf.GetColumnNames()]
        if "event_weight" in col_names:
            rdf_filtered = rdf.Filter("is_rock_trimuon_event(MCTrack)")
        else:
            rdf_filtered = rdf.Filter("is_rock_trimuon_event(MCTrack)").Define(
                "event_weight", "get_trimuon_weight(MCTrack)"
            )

        count_res = rdf_filtered.Count()
        sum_res = rdf_filtered.Sum("event_weight")

        snapshot_opts = ROOT.RDF.RSnapshotOptions()
        snapshot_opts.fMode = "RECREATE"
        snapshot_opts.fLazy = True
        snapshot_res = rdf_filtered.Snapshot("cbmsim", output_file, "", snapshot_opts)

        unweighted_selected = count_res.GetValue()
        total_weighted_yield = sum_res.GetValue()
        snapshot_res.GetValue()

        print(f"  Processed {total_events:,} events | Selected: {unweighted_selected} (Yield: {total_weighted_yield:.4f})")

        f_in = ROOT.TFile.Open(input_file, "READ")
        if f_in and not f_in.IsZombie():
            f_out_meta = ROOT.TFile.Open(output_file, "UPDATE")
            if f_out_meta and not f_out_meta.IsZombie():
                for key in f_in.GetListOfKeys():
                    obj_name = key.GetName()
                    obj_class = key.GetClassName()
                    if obj_class != "TTree" and not f_out_meta.Get(obj_name):
                        obj = f_in.Get(obj_name)
                        if obj:
                            f_out_meta.cd()
                            obj.Write(obj_name, ROOT.TObject.kSingleKey)
                f_out_meta.Close()
            f_in.Close()

    if unweighted_selected == 0:
        print(f"  No trident events selected in '{input_file}'.")
        return 0, 0.0, 0, []

    f_out = ROOT.TFile.Open(output_file, "UPDATE")
    if not f_out or f_out.IsZombie():
        print(f"Error: Could not open output file '{output_file}' for writing displays.")
        return unweighted_selected, total_weighted_yield, 0, []
    ROOT.SetOwnership(f_out, False)

    tree = f_out.Get("cbmsim")
    if not tree:
        print(f"Error: 'cbmsim' TTree not found in '{output_file}'.")
        return unweighted_selected, total_weighted_yield, 0, []

    displays_dir = f_out.mkdir("EventDisplays")
    if not displays_dir:
        displays_dir = f_out.Get("EventDisplays")
    ROOT.SetOwnership(displays_dir, False)

    max_to_plot = max_displays if max_displays > 0 else tree.GetEntries()
    events_plotted = 0
    sec_inducing_events = []

    for i_event, event in enumerate(tree):
        if events_plotted >= max_to_plot:
            break

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

            # Inside a single file, EventHeader event number is 100% unique!
            z_vert_val = my_daughters[0][3].GetStartZ() if my_daughters else z_min_2d
            tag_suffix = "_SecMuon" if is_secondary_inducing else ""
            evt_dir_name = f"Event_{evt_num}_Zvert{int(z_vert_val)}cm{tag_suffix}"
            evt_dir = displays_dir.mkdir(evt_dir_name)
            if not evt_dir:
                evt_dir = displays_dir.Get(evt_dir_name)
            ROOT.SetOwnership(evt_dir, False)
            evt_dir.cd()

            if is_secondary_inducing:
                sec_inducing_events.append((evt_num, i_event, z_vert_val, primary_idx, inducing_idx, output_file))
                print(f"  [Secondary-Induced Trident] Event {evt_num} (Entry {i_event}): Pair induced by Secondary Muon Track {inducing_idx} at Z={z_vert_val:.1f} cm")

            banner_text = f"Event #{evt_num}  (Entry {i_event})"
            if is_secondary_inducing:
                banner_text += "   #color[616]{[Secondary Muon Induced]}"

            # ----------------- XZ Projection -----------------
            c_xz = ROOT.TCanvas(f"XZ_Ev_{evt_num}_F_{file_tag}_Idx_{i_event}", f"XZ Projection (Event #{evt_num})", 950, 650)
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
            c_yz = ROOT.TCanvas(f"YZ_Ev_{evt_num}_F_{file_tag}_Idx_{i_event}", f"YZ Projection (Event #{evt_num})", 950, 650)
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
            c_3d = ROOT.TCanvas(f"3D_Ev_{evt_num}_F_{file_tag}_Idx_{i_event}", f"3D View (Event #{evt_num})", 850, 850)
            ROOT.SetOwnership(c_3d, False)

            h3 = ROOT.TH3F(f"axis_box_3d_{evt_num}_F_{file_tag}_Idx_{i_event}", ";Z [cm];X [cm];Y [cm]",
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
            banner_3d.DrawLatex(0.12, 0.92, banner_text)

            c_3d.Write("3D_View")

            events_plotted += 1

    f_out.Close()
    return unweighted_selected, total_weighted_yield, events_plotted, sec_inducing_events

def main():
    default_input = "/eos/experiment/sndlhc/MonteCarlo/ThreeMuons/sndLHC.Ntuple-TGeant4_boost100LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_filteredAtScoringPlane_digCPP-2*.root"
    default_output = "trimuon_events_selected.root"
    default_geofile = "python/geofile_full.Ntuple-TGeant4_boost100.0.root"

    parser = argparse.ArgumentParser(description="Filter Trident Events and Save TTree + Displays per Input File")
    parser.add_argument("-i", "--input", dest="input_pattern", default=default_input, help="Input ROOT file path or wildcard pattern")
    parser.add_argument("-g", "--geofile", dest="geofile_filename", default=default_geofile, help="ROOT geometry file")
    parser.add_argument("-o", "--output", dest="output_filename", default=default_output, help="Output ROOT file base or pattern (e.g. trimuon_events_selected_%%s.root)")
    parser.add_argument("-n", "--max_displays", dest="max_displays", type=int, default=0, help="Max event displays to generate per file (0 for all)")
    args = parser.parse_args()

    ds = DataSet(name="mc_3mu_rock", path_pattern=args.input_pattern, is_mc=True)
    if ds.n_files == 0:
        print(f"Error: No matching input files found for pattern: {args.input_pattern}")
        sys.exit(1)

    file_list = [str(f) for f in ds.file_vec]

    print(f"==========================================")
    print(f"Filter & Display Dataset: {ds.name}")
    print(f"Input Files Found      : {len(file_list)}")
    print(f"Processing Mode        : 1 Output File per Input File")
    print(f"==========================================")

    geo_elements, ggeo, f_geo = load_detector_geometry(args.geofile_filename)

    total_unweighted = 0
    total_weighted = 0.0
    total_displays = 0
    all_sec_inducing = []
    created_output_files = []

    for idx, input_file in enumerate(file_list, 1):
        if len(file_list) == 1 and not ("%s" in args.output_filename or len(file_list) > 1):
            out_file = args.output_filename
        else:
            out_file = get_output_filename(input_file, args.output_filename)

        file_tag = os.path.splitext(os.path.basename(input_file))[0]
        print(f"\n[{idx}/{len(file_list)}]")
        sel_count, sel_yield, disp_count, sec_list = process_single_file(
            input_file, out_file, geo_elements, ggeo, f_geo, args.max_displays, file_tag
        )

        if sel_count > 0:
            total_unweighted += sel_count
        if sel_yield > 0:
            total_weighted += sel_yield
        total_displays += disp_count
        all_sec_inducing.extend(sec_list)
        created_output_files.append(out_file)

    print(f"\n==========================================")
    print(f"SUCCESS: Processing Complete Across All Files!")
    print(f"Total Output Files Created : {len(created_output_files)}")
    print(f"Total Selected Events      : {total_unweighted}")
    if total_weighted > 0:
        print(f"Total Weighted Yield       : {total_weighted:.4f}")
    print(f"Total Displays Generated   : {total_displays}")
    print(f"Secondary-Induced Events   : {len(all_sec_inducing)}")
    print(f"==========================================")
    if all_sec_inducing:
        print("\nList of Secondary-Muon Induced Trident Events:")
        for ev_num, entry_idx, zv, p_idx, ind_idx, out_f in all_sec_inducing:
            print(f"  - File: {out_f} | Event #{ev_num} (Entry {entry_idx}) | Primary Track: {p_idx}, Inducing Track: {ind_idx} | Vertex Z = {zv:.1f} cm")
        print("==========================================")

if __name__ == "__main__":
    main()
