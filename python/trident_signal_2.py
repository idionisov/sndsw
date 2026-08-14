import os
import warnings
import ROOT

ROOT.EnableImplicitMT(6)

sndsw_path = os.environ["SNDSW_ROOT"]

ROOT.gInterpreter.ProcessLine('#include "sndScifiHit.h"')
ROOT.gInterpreter.ProcessLine('#include "ShipMCTrack.h"')
ROOT.gInterpreter.ProcessLine(f'#include "{sndsw_path}/analysis/tools/sndSciFiTools.h"')
ROOT.gInterpreter.ProcessLine(f'#include "{sndsw_path}/analysis/tools/sndTchainGetter.h"')
ROOT.gInterpreter.ProcessLine(f'#include "{sndsw_path}/analysis/tools/sndGeometryGetter.h"')

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

    // Print every 500,000 events
    if (c % 500000 == 0 || (total > 0 && c == total)) {
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
    std::vector<int> primary_muons;

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

        if (std::abs(pdg) == 13 && mother_id == -1) {
            primary_muons.push_back(i);
        } else if (std::abs(pdg) == 13 && proc_id == 5) {
            pair_produced_muons.push_back({i, mother_id, pdg, tr->GetStartZ(), tr->GetWeight()});
        }
    }

    for (int p_mu_idx : primary_muons) {
        bool has_mu_minus = false;
        bool has_mu_plus = false;
        bool all_zv_ok = true;
        double first_weight = 1.0;
        bool found_daughter = false;

        for (const auto& d : pair_produced_muons) {
            if (d.mother_id == p_mu_idx) {
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
    def __init__(self, name, path_pattern, is_mc=False, tree_name=None, metadata=None):
        self.name = name
        self.is_mc = is_mc
        self.metadata = metadata or {}
        self.tree_name = tree_name or ("cbmsim" if is_mc else "rawConv")
        
        # C++ vector to store all resolved file paths for RDataFrame
        self.file_vec = ROOT.std.vector('std::string')()
        
        patterns = path_pattern if isinstance(path_pattern, (list, tuple)) else [path_pattern]
        for pat in patterns:
            self.add_pattern(pat)

    def add_pattern(self, pattern):
        """Uses ROOT's C++ system pipe to expand EOS wildcards natively without Python's glob."""
        if "*" in pattern or "?" in pattern:
            cmd = f"ls -1 {pattern} 2>/dev/null"
            raw_output = ROOT.gSystem.GetFromPipe(cmd).Data()
            matched_files = [f.strip() for f in raw_output.split("\n") if f.strip().endswith(".root")]
        else:
            matched_files = [pattern]

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


# Target pattern containing wildcards
input_path = "/eos/experiment/sndlhc/MonteCarlo/ThreeMuons/sndLHC.Ntuple-TGeant4_boost100LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_filteredAtScoringPlane_digCPP-2*.root"
output_filename = "trimuon_events_0.root"

TRI = DataSet(
    name="mc_3mu_rock",
    path_pattern=input_path,
    is_mc=True,
    tree_name="cbmsim",
    metadata={"process": "ThreeMuons_boost100LHC_rock"},
)

print(f"==========================================")
print(f"DataSet    : {TRI.name}")
print(f"Files Found: {TRI.n_files} (resolved natively via ROOT gSystem)")
print(f"==========================================")

print(f"Building RDataFrame graph for dataset '{TRI.name}'...")

rdf = TRI.to_rdf()

# Get total event count directly from RDataFrame to avoid TChain deadlocks
total_events = rdf.Count().GetValue()

rdf_progress = rdf.Filter("print_progress()")
rdf_filtered = rdf_progress.Filter("is_rock_trimuon_event(MCTrack)").Define(
    "event_weight", "get_trimuon_weight(MCTrack)"
)

count_res = rdf_filtered.Count()
sum_res = rdf_filtered.Sum("event_weight")

snapshot_opts = ROOT.RDF.RSnapshotOptions()
snapshot_opts.fMode = "RECREATE"
snapshot_opts.fLazy = True  # Ensures Snapshot execution is delayed until GetValue()
snapshot_res = rdf_filtered.Snapshot("cbmsim", output_filename, "", snapshot_opts)

print(f"\n---> Starting dataset '{TRI.name}' ({total_events:,} events across {TRI.n_files} files)")
print("Processing events...")

ROOT.reset_progress(total_events)

# Single execution trigger for all lazy actions
unweighted_selected = count_res.GetValue()
total_weighted_yield = sum_res.GetValue()
snapshot_res.GetValue()

print(f"\nTotal unweighted selected events: {unweighted_selected}")
print(f"Total weighted physical yield: {total_weighted_yield}")

# Copy metadata/keys from the first file into the output ROOT file
if TRI.n_files > 0:
    first_file_path = str(TRI.file_vec[0])
    f_in = ROOT.TFile.Open(first_file_path, "READ")
    if f_in and not f_in.IsZombie():
        f_out = ROOT.TFile.Open(output_filename, "UPDATE")
        if f_out and not f_out.IsZombie():
            for key in f_in.GetListOfKeys():
                obj_name = key.GetName()
                obj_class = key.GetClassName()
                if obj_class != "TTree" and not f_out.Get(obj_name):
                    obj = f_in.Get(obj_name)
                    if obj:
                        f_out.cd()
                        obj.Write(obj_name, ROOT.TObject.kSingleKey)
            f_out.Close()
        f_in.Close()

print(f"Successfully saved selected trimuon events to '{output_filename}' maintaining full file structure.")
