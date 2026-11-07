#!/usr/bin/env python
import ROOT,os,sys,subprocess,atexit,time,csv
from XRootD import client
from XRootD.client.flags import DirListFlags, OpenFlags, MkDirFlags, QueryCode
import Monitor, SndlhcMuonReco, TimeWalk

def pyExit():
    os._exit(0)
atexit.register(pyExit)

def _lookup_path_from_csv(run_number, csv_name):
    """Retrieve run mapping from analysis/tools/<csv_name>."""
    sndsw_root = os.environ.get("SNDSW_ROOT", "")
    if not sndsw_root:
        sndsw_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    csv_file = os.path.join(sndsw_root, "analysis/tools", csv_name)
    if not os.path.exists(csv_file):
        return None
    try:
        with open(csv_file, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                min_run = int(row["min_run_number"])
                max_run = int(row["max_run_number"])
                if min_run <= run_number <= max_run:
                    path = row["path"]
                    if "root://" in path:
                        path = "/" + path.split("//", 2)[-1].lstrip("/")
                    return path
    except Exception as e:
        print(f"Warning: could not parse {csv_file}: {e}")
    return None

from argparse import ArgumentParser
import argparse

parser = ArgumentParser(description="SND@LHC Muon System TimeWalk and Time Alignment runner")

# Monitor.py compatibility settings
parser.add_argument("-M", "--online", dest="online", help="online mode", default=False, action="store_true")
parser.add_argument("-A", "--auto", dest="auto", help="run in auto mode online monitoring", default=False, action="store_true")
parser.add_argument("--Nupdate", dest="Nupdate", help="frequency of updating online plots", default=100, type=int)
parser.add_argument("--Nlast", dest="Nlast", help="last N events to analyze on file", default=10, type=int)
parser.add_argument("--Cosmics", dest="cosmics", help="use default data stream if no beam", action="store_true", default=False)
parser.add_argument("--sudo", dest="sudo", help="update files on EOS", default=False, action="store_true")
parser.add_argument("-cpp", "--convRawCPP", action="store_true", dest="FairTask_convRaw", help="convert raw data using ConvRawData FairTask", default=False)
parser.add_argument("--withCalibration", action="store_true", dest="makeCalibration", help="make QDC and TDC calibration, not taking from raw data", default=False)
parser.add_argument("--customEventChain", dest="customEventChain", action="store_true")
parser.add_argument("-c", "--command", dest="command", help="command", default="")
parser.add_argument("--goodEvents", dest="goodEvents", action="store_true", default=False)
parser.add_argument("--withTrack", dest="withTrack", action="store_true", default=False)
parser.add_argument("--nTracks", dest="nTracks", default=0, type=int)
parser.add_argument("--save", dest="save", action="store_true", default=False)
parser.add_argument("--postScale", dest="postScale", help="post scale events, 1..10..100", default=-1, type=int)

# Run & IO settings
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int, default=-1)
parser.add_argument("-y", "--year", dest="year", help="data taking year (e.g. 2022, 2023, 2024, 2025) to override automatic detection", default=None, type=str)
parser.add_argument("-p", "--path", dest="path", help="path to converted data", default="./", type=str)
parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC / custom data", type=str, default=None)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geometry file", required=False, default=False)
parser.add_argument("-P", "--partition", dest="partition", help="partition of data", type=int, default=-1)
parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events", default=-1, type=int)
parser.add_argument("-s", "--nStart", dest="nStart", help="first event", default=0, type=int)
parser.add_argument("-b", "--heartBeat", dest="heartBeat", help="heart beat", default=10000, type=int)
parser.add_argument("-d", "--debug", dest="debug", help="debug mode", type=int, default=0)
parser.add_argument("--batch", dest="batch", help="batch mode", default=False, action="store_true")
parser.add_argument("--server", dest="server", help="xrootd server", default=os.environ.get("EOSSHIP", ""))

# Timing directories & parameters
parser.add_argument("-o", "--outputDir", dest="outputDir", help="base output directory for histograms and results", default="./Timing", type=str)
parser.add_argument("--outFile", dest="outFile", help="explicit output ROOT file path (overrides default naming)", default=None, type=str)
parser.add_argument("--paramDir", dest="paramDir", help="directory holding alignment and timewalk calibration json files", default=None, type=str)
parser.add_argument("--histFormatFile", dest="histFormatFile", help="path to custom histogram formatting JSON file", default=None, type=str)
parser.add_argument("-D", "--datalocation", dest="datalocation", help="data location environment (physics, commissioning, H8)", type=str, default="physics")

parser.add_argument("--Task", dest="Task", help="Task name (TimeWalk)", default="TimeWalk")
parser.add_argument("--referencesystem", dest="referencesystem", help="Reference system: 0 (both, default), 1 (SciFi), 3 (DS)", default="0")
parser.add_argument("-t", "--trackType", dest="trackType", help="Track type (DS, Scifi, or ScifiDS; defaults according to --referencesystem)", default=None)
parser.add_argument("--nStations", dest="nStations", help="Number of DS stations used in track fit", type=int, default=3)
parser.add_argument("--TWCorrectionRun", dest="TWCorrectionRun", help="Run to take TW correction parameters from", type=str, default="005408")
parser.add_argument("--AlignmentRun", dest="AlignmentRun", help="Alignment run number", type=int, default=None)
parser.add_argument("--state", dest="state", help="State (uncorrected / corrected)", type=str, default="uncorrected")
parser.add_argument("--CorrectionType", dest="CorrectionType", help="Type of correction parameterization", default=5, type=int)

# Flags & selections
parser.add_argument("--XT", dest="XT", help="Cross-talk removal / fiducial selection", action="store_true")
parser.add_argument("--OneHitPerSystem", dest="OneHitPerSystem", action="store_true")
parser.add_argument("--SlopesCut", dest="SlopesCut", action="store_true")
parser.add_argument("--nSiPMsCut", dest="nSiPMsCut", action="store_true")
parser.add_argument("--simulation", dest="simulation", action="store_true")
parser.add_argument("--simMode", dest="simMode", type=str, default=None)
parser.add_argument("--simTest", dest="simTest", action="store_true")
parser.add_argument("--simLumiFile", dest="simLumiFile", help="path to simulation luminosities json", default=None, type=str)
parser.add_argument("--numuStudy", dest="numuStudy", action="store_true")
parser.add_argument("--nueStudy", dest="nueStudy", action="store_true")
parser.add_argument("--interactive", dest="interactive", action="store_true", default=False)

# Effective speed of light fitting options
parser.add_argument("--minEntries", dest="minEntries", help="Minimum entries in 2D histogram for c_scint fit", type=int, default=1000)
parser.add_argument("--fitMin", dest="fitMin", help="Lower bound in x_pred [cm] for c_scint fit", type=float, default=-55.0)
parser.add_argument("--fitMax", dest="fitMax", help="Upper bound in x_pred [cm] for c_scint fit", type=float, default=-5.0)
parser.add_argument("--nSegments", dest="nSegments", help="Number of sub-segments for systematic uncertainty", type=int, default=5)
parser.add_argument("--minSliceEntries", dest="minSliceEntries", help="Minimum entries per x-slice for MPV fit", type=int, default=15)
parser.add_argument("--noCscintFit", dest="noCscintFit", help="Disable c_scint speed fitting in zeroth mode", action="store_true")
parser.add_argument("--fitSmallSiPMs", dest="fitSmallSiPMs", help="Enable c_scint fitting for small SiPM channels", action="store_true")
parser.add_argument("--qdcMin", dest="qdcMin", help="Lower bound for QDC fit range (default: 1.0)", type=float, default=1.0)
parser.add_argument("--minBinEntries", dest="minBinEntries", help="Minimum entries per QDC bin for time-walk fit", type=int, default=10)
parser.add_argument("--minTotalEntries", dest="minTotalEntries", help="Minimum total entries in dtvqdc histogram for time-walk fit", type=int, default=500)
parser.add_argument("--saveSlices", dest="saveSlices", help="Save individual 1D slice histograms and Landau fits in channel subdirectories", action="store_true", default=False)
parser.add_argument("--speedEvents", dest="speedEvents", help="Number of events to process in Step 1 for c_scint (default: 200000, 0=use all)", type=int, default=200000)
parser.add_argument("--cscintFile", dest="cscintFile", help="Path to precomputed c_scint JSON to skip Step 1", type=str, default=None)
parser.add_argument("--noFastFilter", dest="noFastFilter", help="Disable fast IP1/US-hit pre-filtering before tracking", action="store_true", default=False)
parser.add_argument("--fillOnly", dest="fillOnly", help="Only accumulate 2D histograms and save to raw_histograms.root (skip fitting & JSON)", action="store_true", default=False)
parser.add_argument("--fitOnly", dest="fitOnly", help="Path to merged histograms ROOT file to run fitting and produce final outputs", type=str, default=None)
parser.add_argument("--speedOnly", dest="speedOnly", help="Only extract c_scint constants in Step 1, export JSON, and exit", action="store_true", default=False)

options = parser.parse_args()
options.fastTrackingFilter = not options.noFastFilter

if options.trackType is None:
    ref_arg = str(options.referencesystem).lower()
    if ref_arg in ('1', 'sf', 'scifi', 'sf_only'):
        options.trackType = 'Scifi'
    elif ref_arg in ('3', '2', 'ds', 'ds_only'):
        options.trackType = 'DS'
    else:  # 0, 'both', default
        options.trackType = 'ScifiDS'

if not options.interactive:
    ROOT.gROOT.SetBatch(True)

# Normalize work/output directories
if not options.outputDir:
    options.outputDir = "./Timing"
if not hasattr(options, 'saveTo') or not options.saveTo:
    options.saveTo = os.path.join(options.outputDir, f"run{str(options.runNumber).zfill(6)}") + "/"
# if no geofile given, use defaults according to run number

# Year defaults mapping
YEAR_DEFAULTS = {
    "2022": ("/eos/experiment/sndlhc/convertedData/physics/2022/", "geofile_sndlhc_TI18_V4_2022.root"),
    "2023": ("/eos/experiment/sndlhc/convertedData/physics/2023/", "geofile_sndlhc_TI18_V3_2023.root"),
    "2024": ("/eos/experiment/sndlhc/convertedData/physics/2024/run_241/", "geofile_sndlhc_TI18_V12_2024.root"),
    "2025": ("/eos/experiment/sndlhc/convertedData/physics/2025/run_251/", "geofile_sndlhc_TI18_V8_2025.root"),
}

if options.year:
    y_str = str(options.year)
    if y_str in YEAR_DEFAULTS:
        def_path, def_geo = YEAR_DEFAULTS[y_str]
        if options.path == "./" or not options.path:
            options.path = def_path
        if not options.geoFile:
            options.geoFile = def_geo
    else:
        print(f"Warning: year {options.year} not recognized in defaults.")

if options.runNumber > 0:
    csv_data_path = _lookup_path_from_csv(options.runNumber, "data_paths.csv")
    csv_geo_path = _lookup_path_from_csv(options.runNumber, "geo_paths.csv")

    if options.path == "./" or not options.path:
        if csv_data_path:
            options.path = csv_data_path

    if not options.geoFile and csv_geo_path:
        options.geoFile = os.path.basename(csv_geo_path)

    # Cross-check
    if csv_data_path and csv_geo_path:
        for y in ("2022", "2023", "2024", "2025"):
            expected_prefix = f"/eos/experiment/sndlhc/convertedData/physics/{y}"
            if expected_prefix in csv_data_path and expected_prefix not in csv_geo_path:
                print(f"Warning: Year mismatch between data_paths ({csv_data_path}) and geo_paths ({csv_geo_path})!")

if options.runNumber < 0 and not options.geoFile:
    print('No run number given and no geoFile. Do not know what to do. Exit.')
    exit()

if not options.geoFile:
    for y in ("2022", "2023", "2024", "2025"):
        if f"/eos/experiment/sndlhc/convertedData/physics/{y}" in options.path:
            options.geoFile = YEAR_DEFAULTS[y][1]
            break

# works only for runs on EOS
# if not options.server.find('eos')<0 and not options.simulation:
#     if options.path.find('2023')!=-1:
#         rawDataPath='/eos/experiment/sndlhc/raw_data/physics/2023/'
#     elif options.path.find('2022')!=-1:
#         rawDataPath='/eos/experiment/sndlhc/raw_data/physics/2022/'
#     else:
#         rawDataPath='/eos/experiment/sndlhc/raw_data/commissioning/TI18/data/'
#     options.rawDataPath=rawDataPath
#     print(f'rawDataPath: {rawDataPath}')
#     runDir=rawDataPath+'run_'+str(options.runNumber).zfill(6)
#     jname = "run_timestamps.json"
#     dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+runDir,shell=True) )
#     if jname in dirlist:
#        with client.File() as f:
#           f.open(options.server+runDir+"/run_timestamps.json")
#           status, jsonStr = f.read()
#        exec("date = "+jsonStr.decode())
#        options.startTime = date['start_time'].replace('Z','')
#        if 'stop_time' in date:
#            options.startTime += " - "+ date['stop_time'].replace('Z','')

# prepare tasks:
if options.fitOnly:
    print("\n" + "=" * 70)
    print(f"[FIT-ONLY] Running in fit-only mode on merged histograms:")
    print(f"           Input file: {options.fitOnly}")
    print("=" * 70)
    # Initialize minimal Monitor with empty FairTasks (loads geometry without tracking)
    M = Monitor.Monitoring(options, [])
    options.mode = 'calibration'
    tw_task = TimeWalk.TimeWalk(options, M)
    if options.cscintFile:
        tw_task.LoadCscintFile(options.cscintFile)
    if not tw_task.LoadHistogramsFromRoot(options.fitOnly):
        print(f"[ERROR] Failed to load histograms from {options.fitOnly}. Exiting.")
        exit(1)

    tw_task.WriteOutHistograms()
    print("[SUCCESS] Fit-only completed successfully.")
    exit(0)

FairTasks = []
houghTransform = False # under construction, not yet tested
if houghTransform:
    muon_reco_task = SndlhcMuonReco.MuonReco()
    muon_reco_task.SetName("houghTransform")
    FairTasks.append(muon_reco_task)
else:
    import SndlhcTracking
    trackTask = SndlhcTracking.Tracking()
    trackTask.SetName('simpleTracking')
    FairTasks.append(trackTask)

M = Monitor.Monitoring(options,FairTasks)
monitorTasks = {}

if options.nEvents < 0 or options.nEvents > M.GetEntries():   options.nEvents = M.GetEntries()

options.mode = 'calibration'
tw_task = TimeWalk.TimeWalk(options, M)
monitorTasks['TimeWalk'] = tw_task

start, nEvents = options.nStart, options.nEvents
deciles = [i / 10 for i in range(11)]

# Optimize tree branch reading: activate only needed detector branches
if hasattr(M, 'eventTree') and M.eventTree:
    M.eventTree.SetBranchStatus('*', 0)
    for b in ('EventHeader*', 'Digi_MuFilterHits*', 'Digi_ScifiHits*'):
        M.eventTree.SetBranchStatus(b, 1)

# Check if precomputed c_scint file is provided
skip_step1 = False
if options.cscintFile:
    if tw_task.LoadCscintFile(options.cscintFile):
        skip_step1 = True

if not skip_step1:
    n_speed = min(nEvents, options.speedEvents) if options.speedEvents > 0 else nEvents
    print("\n" + "=" * 70)
    print(f"[STEP 1/2] Processing {n_speed} events for Signal Speed (c_scint)...")
    print("=" * 70)
    tw_task.step = 1
    for n in range(start, start + n_speed):
        event = M.GetEvent(n)
        if ((n - start) / n_speed) in deciles:
            progstr = f'Step 1 Progress: {n}/{start+n_speed}, {int(100*(n-start)/n_speed)}%'
            print(progstr)
        tw_task.ExecuteEvent(M.eventTree)

    print("\n" + "=" * 70)
    print("Step 1 complete. Extracting c_scint constants for photon ToF correction...")
    print("=" * 70)
    tw_task.FitAndSaveCscint(root_file=None)

    cscint_out = os.path.join(options.outputDir, f"run{str(options.runNumber).zfill(6)}_cscintvalues.json")
    os.makedirs(options.outputDir, exist_ok=True)
    import json
    with open(cscint_out, 'w') as jf:
        json.dump(tw_task.cscintvalues_summary, jf, indent=2)
    print(f"[INFO] Exported c_scint constants to {cscint_out}")

    if options.speedOnly:
        print("[INFO] --speedOnly requested: Step 1 completed successfully, exiting.")
        os._exit(0)
else:
    print("\n" + "=" * 70)
    print("Skipping Step 1 (loaded precomputed c_scint constants).")
    print("=" * 70)

print("\n" + "=" * 70)
print(f"[STEP 2/2] Processing {nEvents} events for Scintillator ToF & Time-Walk...")
print("=" * 70)
tw_task.step = 2
for n in range(start, start + nEvents):
    event = M.GetEvent(n)
    if ((n - start) / nEvents) in deciles:
        progstr = f'Step 2 Progress: {n}/{start+nEvents}, {int(100*(n-start)/nEvents)}%'
        print(progstr)
    tw_task.ExecuteEvent(M.eventTree)

if options.fillOnly:
    print("\n" + "=" * 70)
    print("Saving raw histograms (fillOnly mode, skipping fits & JSON)...")
    print("=" * 70)
    tw_task.SaveRawHistograms()
elif not options.debug:
    print("\n" + "=" * 70)
    print("Fitting Time-Walk parameters & saving consolidated outputs...")
    print("=" * 70)
    tw_task.WriteOutHistograms()

# Clean up presenter file created by Monitor.py if not used
if hasattr(M, 'presenterFile') and M.presenterFile.IsOpen():
    pres_path = M.presenterFile.GetName()
    M.presenterFile.Close()
    if os.path.exists(pres_path) and (not os.path.isabs(pres_path) or pres_path.startswith('run') or 'run' + str(options.runNumber).zfill(6) + '_' in pres_path):
        try:
            os.remove(pres_path)
        except Exception:
            pass
