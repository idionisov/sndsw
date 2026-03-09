import os
import random
import ROOT
import SndlhcGeo
import SndlhcMuonReco
import SndlhcTracking
import argparse
import pandas as pd
import csv
import time
import datetime

sndsw_path = os.environ['SNDSW_ROOT']
ROOT.gInterpreter.ProcessLine(f'#include "{sndsw_path}/analysis/tools/sndSciFiTools.h"')

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file',  type=str, help='Path to input file.')
    parser.add_argument('-o', '--output-file', type=str, help='Output CSV filename')
    parser.add_argument('-s', '--start-event', type=int, help='Start event number.', default=0)
    parser.add_argument('-n', '--n-events', type=int, help='Number of events.', default=1000000)
    parser.add_argument('-f', '--fraction', type=float, default=1.0, help='Fraction of events to process (0.0 to 1.0)')
    return parser.parse_args()

def get_geofile(year: int) -> str:
    if year==2022:
        return "/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V4_2022.root"
    elif year==2023:
        return "/eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V3_2023.root"
    elif year==2024:
        return "/eos/experiment/sndlhc/convertedData/physics/2024/geofile_sndlhc_TI18_V12_2024.root"
    elif year==2025:
        return "/eos/experiment/sndlhc/convertedData/physics/2024/geofile_sndlhc_TI18_V8_2025.root"
    else:
        return None


def veto_is_activated(mf_hits: ROOT.TClonesArray) -> bool:
    for mf_hit in mf_hits:
        if mf_hit.GetSystem() == 1:
            return True
    return False

def getNhitsPerPlane(sf_hits: ROOT.TClonesArray):
    nV = {plane: 0 for plane in range(1, 6)}
    nH = {plane: 0 for plane in range(1, 6)}

    for sf_hit in sf_hits:
        station = sf_hit.GetStation()
        if sf_hit.isVertical():
            nV[station] += 1
        else:
            nH[station] += 1

    return nV, nH


def get_sum_hit_weight_density(
    sf_hits: ROOT.TClonesArray,
    radius: float = 40,
    min_check: bool = False,
    min_hit_density: int = 1000000
): 
    density = 0
    for sf in sf_hits:
        id_ = sf.GetChannelID()
        density += ROOT.snd.analysis_tools.densityScifi(id_, sf_hits, radius, min_hit_density, min_check)
    return density


def get_sf_qdc(sf_hits: ROOT.TClonesArray) -> float:
    qdc = 0 
    for sf_hit in sf_hits:
        qdc += sf_hit.GetSignal()
    return qdc
    
    


def main():
    args = get_args()
    t0 = time.time()  # Used later as start_time


    INPUT_FILE = args.input_file
    if not INPUT_FILE.endswith(".root"):
        INPUT_FILE += ".root"

    OUT_CSV = args.output_file
    if not OUT_CSV.endswith(".csv"):
        OUT_CSV += ".csv"

    # 1. Open the file manually first to detect the Year and Run ID
    f = ROOT.TFile.Open(INPUT_FILE)
    if not f or f.IsZombie():
        print(f"Error: Could not open {INPUT_FILE}")
        return

    tree = f.Get("rawConv")
    tree.GetEvent(0)
    daq_run = tree.EventHeader.GetRunId()

    df = pd.read_csv("tridents.csv")
    WHITELIST = list(df.query(f"run == {daq_run}")["event_number"])
    print(WHITELIST)
    del df

    # Correctly detect year to get the right geofile
    utc_timestamp = tree.EventHeader.GetUTCtimestamp()
    year = datetime.datetime.fromtimestamp(utc_timestamp).year
    GEOFILE = get_geofile(year)

    # 2. Initialize Geometry BEFORE run.Init()
    geo = SndlhcGeo.GeoInterface(GEOFILE)
    lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
    lsOfGlobals.Add(geo.modules['Scifi'])
    lsOfGlobals.Add(geo.modules['MuFilter'])

    # 3. Setup FairRunAna
    run = ROOT.FairRunAna()
    ioman = ROOT.FairRootManager.Instance()
    ioman.SetTreeName("rawConv")

    # Reuse the open file
    source = ROOT.FairFileSource(f)
    run.SetSource(source)

    outFile = ROOT.TMemFile('dummy','CREATE')
    sink = ROOT.FairRootFileSink(outFile)
    run.SetSink(sink)

    xrdb = ROOT.FairRuntimeDb.instance()
    xrdb.getContainer("FairBaseParSet").setStatic()
    xrdb.getContainer("FairGeoParSet").setStatic()

    # 4. Add Tracking Tasks (Match 2dEventDisplay setup)
    HT_tasks = {'muon_reco_task_Sf': SndlhcMuonReco.MuonReco(),
                'muon_reco_task_DS': SndlhcMuonReco.MuonReco(),
                'muon_reco_task_nuInt': SndlhcMuonReco.MuonReco()}
    for task in HT_tasks.values():
        run.AddTask(task)
        task.SetParFile(f"/afs/cern.ch/user/i/idioniso/snd_master/sndsw/trackingParams.xml")
        task.SetHoughSpaceFormat("linearSlopeIntercept")
        task.ForceGenfitTrackFormat()

    HT_tasks['muon_reco_task_Sf'].SetTrackingCase('passing_mu_Sf')
    HT_tasks['muon_reco_task_DS'].SetTrackingCase('passing_mu_DS')
    HT_tasks['muon_reco_task_nuInt'].SetTrackingCase('nu_interaction_products')

    # Add simple tracking task for genfit/material effects initialization consistency
    import SndlhcTracking
    trackTask = SndlhcTracking.Tracking()
    trackTask.SetName('simpleTracking')
    run.AddTask(trackTask)

    # 5. Now Initialize the Run
    run.Init()

    # 6. Initialize Alignment and Tree References
    OT = sink.GetOutTree()
    tree = ioman.GetInTree()
    tree.GetEvent(0)

    # Backward compatibility for early converted events
    if tree.GetBranch('Digi_MuFilterHit'):
        tree.Digi_MuFilterHits = tree.Digi_MuFilterHit

    # Prepare CSV and Counters
    pr = 0
    I = 0
    n_break = args.n_events if args.n_events else tree.GetEntries()
    start_ev = args.start_event if args.start_event else 0

    # Get the tree once
    tree = ioman.GetInTree()

    with open(OUT_CSV, mode='a', newline='') as csv_file:
        fieldnames = ['run', 'event_number', 'activated_veto', 'n_tracks', 'sum_hit_weight_density', 'sum_qdc']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        if not os.path.isfile(OUT_CSV):
            writer.writeheader()

        # Loop over specific range
        start_event = args.start_event
        for i_entry in range(start_ev, start_ev + n_break):

            if i_entry >= tree.GetEntries():
                break

            rc = source.GetInTree().GetEvent(i_entry)
            event_number = tree.EventHeader.GetEventNumber()
            if event_number not in WHITELIST:
                continue
    #                if random.random() > args.fraction:
    #                    continue

            # Progress bar logic
            current_progress = (i_entry - start_ev)
            if (current_progress * 100 / n_break) >= pr:
                elapsed_time = time.time() - t0

                # Use standard divmod to get hours, minutes, and seconds
                h, rem = divmod(int(elapsed_time), 3600)
                m, s = divmod(rem, 60)

                # Correct percentage calculation
                percent = int((current_progress * 100) / n_break)

                print(f"[{percent}%] \t {current_progress:,}/{n_break:,} \t {h:02d}:{m:02d}:{s:02d}")
                pr += 1

            activated_veto = veto_is_activated(tree.Digi_MuFilterHits)
    #            if not if not activated_veto:
    #                continue

            nV, nH = getNhitsPerPlane(tree.Digi_ScifiHits)
            nVpl = sum(1 for v in nV.values() if v != 0)
            nHpl = sum(1 for v in nH.values() if v != 0)

            # Logic check: Typically trident search needs 3 vertical AND 3 horizontal 
            # or a similar strict geometry. Adjust if your logic requires AND vs OR.
            if not ((nVpl >= 3 and nHpl >= 3) or (nVpl >= 3 and nHpl >= 3)):
                continue


            sum_density = get_sum_hit_weight_density(tree.Digi_ScifiHits)
            qdc = get_sf_qdc(tree.Digi_ScifiHits)

            daq_run = tree.EventHeader.GetRunId()

            # Process tracking
            # Re-fetch event for Hough tracking as done in 2dEventDisplay
            rc = source.GetInTree().GetEvent(i_entry)

            # Re-apply alignment per event
            if tree.EventHeader.ClassName() == 'SNDLHCEventHeader':
                geo.modules['Scifi'].InitEvent(tree.EventHeader)
                geo.modules['MuFilter'].InitEvent(tree.EventHeader)

            # Initialize track collection for this event
            OT.Reco_MuonTracks = ROOT.TObjArray(10)

            # Execute SciFi reco task
            sf_task = HT_tasks['muon_reco_task_Sf']
            sf_task.kalman_tracks.Delete()
            sf_task.Exec(0)

            for ht_trk in sf_task.kalman_tracks:
                OT.Reco_MuonTracks.Add(ht_trk)

            n_tracks = OT.Reco_MuonTracks.GetEntries()
            
            # Write results
            writer.writerow({
                'run': daq_run,
                'event_number': event_number,
                'activated_veto': activated_veto,
                'n_tracks': n_tracks,
                'sum_hit_weight_density': sum_density,
                'sum_qdc': qdc
            })

            if n_tracks > 2:
                print(f"Match: Run {daq_run}, Event {event_number}, Tracks {n_tracks}, density: {sum_density}, qdc: {qdc}")
            
            I += 1

    end_time = time.time()
    duration = end_time - t0  # Changed from start_time to t0
    print(f"\nExecution finished in {duration:.2f} seconds.")
    print(f"Results saved to {OUT_CSV}")


if __name__ == "__main__":

    main()

