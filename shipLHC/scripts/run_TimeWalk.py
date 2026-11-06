#!/usr/bin/env python
import ROOT,os,sys,subprocess,atexit,time,csv
from XRootD import client
from XRootD.client.flags import DirListFlags, OpenFlags, MkDirFlags, QueryCode
import Monitor, SndlhcMuonReco, TimeWalk

def pyExit():
    print("Make suicide until solution found for freezing")
    os.system('kill '+str(os.getpid()))
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

# Task & Mode configuration
parser.add_argument("--Task", dest="Task", help="Task name (TimeWalk)", default="TimeWalk")
parser.add_argument("-m", "--mode", dest="mode", help="Analysis mode (zeroth, tof, tw, res, systemalignment, tds0-studies, etc.)", type=str, default="zeroth")
parser.add_argument("--referencesystem", dest="referencesystem", help="Reference system: 3 for DS, 1 for SciFi", type=int, default=3)
parser.add_argument("-t", "--trackType", dest="trackType", help="Track type (DS or Scifi)", default="DS")
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

options = parser.parse_args()

if not options.interactive:
    ROOT.gROOT.SetBatch(True)

# Normalize work/output directories
if not options.outputDir:
    options.outputDir = "./Timing"
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

if options.Task=='TimeWalk':
    if not options.mode:
        print('=='*20+f'\nNo mode given for time walk task. Give mode as zeroth, tof or tw.\n'+'=='*20)
        pyExit()
    monitorTasks['TimeWalk'] = TimeWalk.TimeWalk(options, M)

start, nEvents=options.nStart, options.nEvents
deciles=[i/10 for i in range(11)]
for n in range(start,start+nEvents):
    event = M.GetEvent(n)

    if ((n-start)/nEvents) in deciles:
        progstr=f'Progress: {n}/{start+nEvents}, {int(100*(n-start)/nEvents)}%'
        print('='*len(progstr)+'\n')
        print(progstr, '\n')

    for m in monitorTasks:
        monitorTasks[m].ExecuteEvent(M.eventTree)
if 'TimeWalk' in monitorTasks:
    if not options.debug: monitorTasks['TimeWalk'].WriteOutHistograms()
