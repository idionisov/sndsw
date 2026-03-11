import os
import ROOT
import argparse
import pandas as pd
import csv
import time
import datetime
import numpy as np

# Load SND@LHC specific libraries
for lib in ["libBase", "libShipData", "libshipLHC"]:
    ROOT.gSystem.Load(lib)

import SndlhcGeo
import SndlhcMuonReco
import SndlhcTracking

sndsw_path = os.environ['SNDSW_ROOT']
ROOT.gInterpreter.ProcessLine(f'#include "{sndsw_path}/analysis/tools/sndSciFiTools.h"')

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file',  type=str, required=True, help='Path to input file.')
    parser.add_argument('-o', '--output-file', type=str, required=True, help='Output CSV filename')
    parser.add_argument('-or', '--output-root', type=str, help='Output ROOT filename for selected events')
    parser.add_argument('-s', '--start-event', type=int, help='Start event number.', default=0)
    parser.add_argument('-n', '--n-events', type=int, help='Number of events.', default=1000000)
    parser.add_argument('-f', '--fraction', type=float, default=1.0, help='Fraction of events to process')
    return parser.parse_args()

def get_geofile(year: int) -> str:
    if year == 2022:
        return "/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V4_2022.root"
    elif year == 2023:
        return "/eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V3_2023.root"
    elif year == 2024:
        return "/eos/experiment/sndlhc/convertedData/physics/2024/geofile_sndlhc_TI18_V12_2024.root"
    elif year == 2025:
        return "/eos/experiment/sndlhc/convertedData/physics/2024/geofile_sndlhc_TI18_V8_2025.root"
    return None

def veto_is_activated(mf_hits: ROOT.TClonesArray) -> bool:
    for mf_hit in mf_hits:
        if mf_hit.GetSystem() == 1:
            return True
    return False

def get_sum_hit_weight_density(sf_hits: ROOT.TClonesArray, radius: float = 40, min_check: bool = False, min_hit_density: int = 1000000): 
    density = 0
    for sf in sf_hits:
        id_ = sf.GetChannelID()
        density += ROOT.snd.analysis_tools.densityScifi(id_, sf_hits, radius, min_hit_density, min_check)
    return density

def get_sf_qdc(sf_hits: ROOT.TClonesArray) -> float:
    qdc = 0 
    for sf_hit in sf_hits: qdc += sf_hit.GetSignal()
    return qdc

def getHoughLines(ht_task, event, geo):
    """ Returns (n_tracks, max_dangle) """
    hit_collection = {"pos": [[], [], []], "d": [[], [], []], "vert": [], "system": [], "detectorID": []}
    a, b = ROOT.TVector3(), ROOT.TVector3()
    for sfHit in event.Digi_ScifiHits:
        if not sfHit.isValid(): continue
        geo.modules['Scifi'].GetSiPMPosition(sfHit.GetDetectorID(), a, b)
        hit_collection["pos"][0].append(a.X()); hit_collection["pos"][1].append(a.Y()); hit_collection["pos"][2].append(a.Z())
        hit_collection["d"][0].append(ht_task.Scifi_dx); hit_collection["d"][1].append(ht_task.Scifi_dy); hit_collection["d"][2].append(ht_task.Scifi_dz)
        hit_collection["vert"].append(sfHit.isVertical()); hit_collection["system"].append(0); hit_collection["detectorID"].append(sfHit.GetDetectorID())
    if not hit_collection['pos'][0]: return 0, 0.0
    for k in ["pos", "d"]: hit_collection[k] = np.array(hit_collection[k], dtype=np.float32)
    hit_collection["vert"] = np.array(hit_collection["vert"], dtype=np.bool_)
    hit_collection["system"] = np.array(hit_collection["system"], dtype=np.int32)
    hit_collection["detectorID"] = np.array(hit_collection["detectorID"], dtype=np.int32)
    counts, slopes = {'XZ': 0, 'YZ': 0}, {'XZ': [], 'YZ': []}
    for proj in ['XZ', 'YZ']:
        is_v, ax = (proj == 'XZ'), (0 if proj == 'XZ' else 1)
        h_obj = ht_task.h_ZX if is_v else ht_task.h_ZY
        proj_used = np.zeros(len(hit_collection['pos'][0]), dtype=np.bool_)
        for i_muon in range(ht_task.max_reco_muons):
            m = np.logical_and(hit_collection['vert'] == is_v, ~proj_used)
            if not np.any(m): break
            res = h_obj.fit_randomize(np.dstack([hit_collection['pos'][2][m], hit_collection['pos'][ax][m]])[0],
                                      np.dstack([hit_collection['d'][2][m], hit_collection['d'][ax][m]])[0],
                                      ht_task.n_random, False, False)
            if res[0] in [-1, -999]: break
            hits_rel = SndlhcMuonReco.hit_finder(res[0], res[1], np.dstack([hit_collection['pos'][2][m], hit_collection['pos'][ax][m]]),
                                               np.dstack([hit_collection['d'][2][m], hit_collection['d'][ax][m]]), ht_task.tolerance)
            if len(hits_rel) == 0: break
            if SndlhcMuonReco.numPlanesHit(hit_collection['system'][m][hits_rel], hit_collection['detectorID'][m][hits_rel]) >= ht_task.min_planes_hit:
                counts[proj] += 1; slopes[proj].append(res[0]); proj_used[np.where(m)[0][hits_rel]] = True
            else: break
    max_dangle_XZ = np.ptp(np.arctan(slopes['XZ'])) if counts['XZ'] >= 2 else 0.0
    max_dangle_YZ = np.ptp(np.arctan(slopes['YZ'])) if counts['YZ'] >= 2 else 0.0
    if counts['XZ'] > counts['YZ']: return counts['XZ'], max_dangle_XZ
    elif counts['YZ'] > counts['XZ']: return counts['YZ'], max_dangle_YZ
    else: return counts['XZ'], max(max_dangle_XZ, max_dangle_YZ)

def main():
    args = get_args()
    t0 = time.time()
    INPUT_FILE = args.input_file if args.input_file.endswith(".root") else args.input_file + ".root"
    OUT_CSV = args.output_file if args.output_file.endswith(".csv") else args.output_file + ".csv"
    
    f_in = ROOT.TFile.Open(INPUT_FILE)
    tree_in = f_in.Get("rawConv")
    tree_in.GetEvent(0)
    
    # 1. Initialize Run and Manager
    run = ROOT.FairRunAna()
    ioman = ROOT.FairRootManager.Instance(); ioman.SetTreeName("rawConv")
    run.SetSource(ROOT.FairFileSource(f_in))
    sink = ROOT.FairRootFileSink(ROOT.TMemFile('dummy','CREATE')); run.SetSink(sink)
    xrdb = ROOT.FairRuntimeDb.instance()
    xrdb.getContainer("FairBaseParSet").setStatic(); xrdb.getContainer("FairGeoParSet").setStatic()

    # 2. Geometry Setup (Follow 2dEventDisplay pattern)
    year = datetime.datetime.fromtimestamp(tree_in.EventHeader.GetUTCtimestamp()).year
    geo = SndlhcGeo.GeoInterface(get_geofile(year))
    lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
    # Explicitly replace existing TObject detectors with correctly-typed ones
    for m_name in ['Scifi', 'MuFilter']:
        obj = lsOfGlobals.FindObject(m_name)
        if obj: lsOfGlobals.Remove(obj)
        lsOfGlobals.Add(geo.modules[m_name])

    # 3. Add Tasks
    HT_tasks = {'muon_reco_task_Sf': SndlhcMuonReco.MuonReco(), 'muon_reco_task_DS': SndlhcMuonReco.MuonReco(), 'muon_reco_task_nuInt': SndlhcMuonReco.MuonReco()}
    for task in HT_tasks.values():
        run.AddTask(task)
        task.SetParFile(f"/afs/cern.ch/user/i/idioniso/snd_master/sndsw/trackingParams.xml")
        task.SetHoughSpaceFormat("linearSlopeIntercept"); task.ForceGenfitTrackFormat()

    HT_tasks['muon_reco_task_Sf'].SetTrackingCase('passing_mu_Sf')
    HT_tasks['muon_reco_task_DS'].SetTrackingCase('passing_mu_DS')
    HT_tasks['muon_reco_task_nuInt'].SetTrackingCase('nu_interaction_products')
    run.AddTask(SndlhcTracking.Tracking())
    
    # 4. Initialize Run
    run.Init()

    # Update tree reference
    tree_in = ioman.GetInTree()
    if tree_in.GetBranch('Digi_MuFilterHit'): tree_in.Digi_MuFilterHits = tree_in.Digi_MuFilterHit

    # 5. Whitelist Setup
    df = pd.read_csv("tridents.csv")
    WHITELIST = set(df.query(f"run == {tree_in.EventHeader.GetRunId()}")["event_number"])
    del df

    # 6. Output ROOT Setup
    f_out_root, tree_out = None, None
    if args.output_root:
        out_root_path = args.output_root if args.output_root.endswith(".root") else args.output_root + ".root"
        f_out_root = ROOT.TFile(out_root_path, "RECREATE")
        tree_out = tree_in.CloneTree(0)

    n_break = min(args.n_events, tree_in.GetEntries() - args.start_event)
    pr = 0

    with open(OUT_CSV, mode='a', newline='') as csv_file:
        fieldnames = ['run', 'event_number', 'activated_veto', 'n_tracks', 'max_dangle', 'sum_hit_weight_density', 'sum_qdc']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        if not os.path.isfile(OUT_CSV) or os.stat(OUT_CSV).st_size == 0: writer.writeheader()

        for i_entry in range(args.start_event, args.start_event + n_break):
            if i_entry >= tree_in.GetEntries(): break

            tree_in.GetEvent(i_entry)
            event_number = tree_in.EventHeader.GetEventNumber()
            if event_number not in WHITELIST: continue

            if ((i_entry - args.start_event) * 100 / n_break) >= pr:
                h, rem = divmod(int(time.time() - t0), 3600); m, s = divmod(rem, 60)
                print(f"[{int(pr)}%] \t {i_entry-args.start_event:,}/{n_break:,} \t {h:02d}:{m:02d}:{s:02d}"); pr += 1

            if tree_in.EventHeader.ClassName() == 'SNDLHCEventHeader':
                geo.modules['Scifi'].InitEvent(tree_in.EventHeader)

            n_tracks, max_dangle = getHoughLines(HT_tasks['muon_reco_task_Sf'], tree_in, geo)
            
            if n_tracks >= 3:
                if tree_out: tree_out.Fill()
                print(f"Match: Event {event_number}, Lines {n_tracks}, DAngle {max_dangle:.4f}")

            writer.writerow({
                'run': tree_in.EventHeader.GetRunId(), 'event_number': event_number,
                'activated_veto': veto_is_activated(tree_in.Digi_MuFilterHits),
                'n_tracks': n_tracks, 'max_dangle': max_dangle,
                'sum_hit_weight_density': get_sum_hit_weight_density(tree_in.Digi_ScifiHits),
                'sum_qdc': get_sf_qdc(tree_in.Digi_ScifiHits)
            })

    if f_out_root:
        f_out_root.cd(); tree_out.Write(); f_out_root.Close()
        print(f"Selected events saved to ROOT file.")
    print(f"Finished in {time.time() - t0:.2f}s.")

if __name__ == "__main__":
    main()
