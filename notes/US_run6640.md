# Upstream Scintillator (US) Timing Calibration: Run 6640

This note documents the complete operational pipeline, new configuration/batch files, execution commands, and physics results for the timing calibration of the **Upstream Scintillator (US)** system on **Run 6640** (Physics 2022 dataset).

---

## 1. Directory Structure & Newly Added Files

To automate the full calibration on HTCondor and standardise histogram formatting, two directories were added to the workspace:

```text
/afs/cern.ch/work/i/idioniso/sndVetoUS/
  ├── htcondor/                   # HTCondor batch submission & parallel merging
  │     ├── 0_timewalk.sub        # HTCondor job description file
  │     ├── run_timewalk.sh       # Worker execution wrapper (loads CVMFS + runs TimeWalk.py)
  │     ├── generate_args_timewalk.py # Parameterized batch chunk generator
  │     ├── merge_files.sh        # 16-worker parallel merger for 600 split channel ROOT files
  │     ├── args_timewalk_zeroth.txt # 50 jobs for mode 'zeroth' (10M events)
  │     ├── args_timewalk_tof.txt    # 50 jobs for mode 'tof' (10M events)
  │     └── args_timewalk_tw.txt     # 50 jobs for mode 'tw' (10M events)
  │
  ├── config/                     # Calibration configuration files
  │     └── TWhistogramformatting.json # Binning & ranges for dtvxpred, dtvqdc, sidetime, etc.
  │
  └── sndsw/shipLHC/scripts/      # Calibration analysis & plotting scripts
        ├── extract_cscint.py     # Speed of light extraction (raw & corrected)
        ├── extract_twparams.py   # 599-channel time-walk & NLL minimization pipeline
        └── plot_timewalk_fit.py  # 2-pad validation ROOT canvas generator (-c all)
```

---

## 2. Calibration Configuration (`config/TWhistogramformatting.json`)

Controls the 2D histogram axis binning and ranges loaded by `TimeWalk.py`:

* **`dtvxpred` (Time vs. Position)**:
  * `uncorrected`: $x \in [-100, 10]\text{ cm}$ (110 bins), $dt \in [-10, 30]\text{ ns}$ (400 bins)
  * `corrected`: $x \in [-100, 10]\text{ cm}$ (110 bins), $dt \in [-5, 5]\text{ ns}$ (200 bins)
* **`dtvqdc` (Time vs. Amplitude)**:
  * `uncorrected`: $\text{QDC} \in [0, 200]\text{ a.u.}$ (200 bins), $dt \in [-10, 30]\text{ ns}$ (400 bins)
  * `corrected`: $\text{QDC} \in [0, 200]\text{ a.u.}$ (200 bins), $dt \in [-5, 5]\text{ ns}$ (200 bins)

---

## 3. Chronological Execution Pipeline

```text
[Stage 1: Mode 'zeroth'] ──► [merge_files.sh] ──► [extract_cscint.py -s uncorrected]
                                                        │
                                                        ▼
                                           ⟨c_scint⟩ = 14.32 ± 0.62 cm/ns
                                                        │
[Stage 2: Mode 'tof']    ──► [merge_files.sh] ──► [extract_twparams.py]
                                                        │
                                                        ▼
                                           σ_sys = 40.4 ± 17.6 ps, χ²/ν = 1.19
                                           twparams.json (599 channels)
                                                        │
                                           ┌────────────┴────────────┐
                                           ▼                         ▼
                                [plot_timewalk_fit.py]     [Stage 4: Mode 'tw']
                                (timewalk_fits_all.root)   (condor_submit 0_timewalk.sub)
                                                                     │
                                                                     ▼
                                                           [merge_files.sh]
                                                                     │
                                                                     ▼
                                                           [extract_cscint.py -s corrected]
                                                                     │
                                                                     ▼
                                                           ⟨c_scint^corr⟩ ≈ 15.5 cm/ns
```

---

## 4. Step-by-Step Executed Stages & Results

### 4.1 Stage 1: Uncorrected Scintillator Speed Extraction (Mode `zeroth`)
* **Job Generation & Submission**:
  ```bash
  cd /afs/cern.ch/work/i/idioniso/sndVetoUS/htcondor
  python3 generate_args_timewalk.py -r 6640 -m zeroth -n 10000000 -c 200000 -o args_timewalk_zeroth.txt
  condor_submit 0_timewalk.sub
  ```
* **Merging**:
  ```bash
  ./merge_files.sh 6640 zeroth
  ```
* **Extraction Command**:
  ```bash
  python3 /afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/extract_cscint.py \
    -r 6640 \
    -p /afs/cern.ch/work/i/idioniso/sndVetoUS-physics2022/ \
    -s uncorrected \
    --plot
  ```
* **Physical Results**:
  * Total Fitted Channels: **599 / 599** (0 failed)
  * Mean Scintillator Speed: **$\langle c_{\text{scint}} \rangle = 14.32 \pm 0.62\text{ cm/ns}$**
    * Left-side SiPMs: $14.24\text{ cm/ns}$
    * Right-side SiPMs: $14.40\text{ cm/ns}$
  * Output Files:
    * `cscintvalues/run006640/run006640_cscintvalues_uncorrected.json`
    * `cscintvalues/run006640/cscint_summary_uncorrected.root`

---

### 4.2 Stage 2: Photon Propagation Correction (Mode `tof`)
* **Job Submission & Merging**:
  ```bash
  python3 generate_args_timewalk.py -r 6640 -m tof -n 10000000 -c 200000 -o args_timewalk_tof.txt
  condor_submit 0_timewalk.sub
  ./merge_files.sh 6640 tof
  ```
* **Output**:
  Generates `dtvqdc_{fixed_ch}_uncorrected` 2D histograms (photon ToF subtracted along the bar).

---

### 4.3 Stage 3: Time-Walk Parameter Extraction & NLL Optimization
* **Extraction Command**:
  ```bash
  python3 /afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/extract_twparams.py \
    -r 6640 \
    -p /afs/cern.ch/work/i/idioniso/sndVetoUS-physics2022/
  ```
* **Physical Results**:
  * Total Channels Fitted: **599 / 599** (0 failed)
  * Mean Error Floor: **$\langle \sigma_{\text{sys}} \rangle = 40.4 \pm 17.6\text{ ps}$** (Matches thesis target $\approx 40\text{ ps}$)
  * Mean $\chi^2/\nu$: **$1.19$** (Matches target $\approx 1.0$)
  * Mean RMS Residual: **$0.81\%$** (Matches target $< 1.0\%$)
  * Output Files:
    * `Polyparams/run006640/twparams.json` (Global dictionary loaded by `AnalysisFunctions.py`)
    * `Polyparams/run006640/polyparams9_{fixed_ch}.json` / `.csv`
    * `Polyparams/run006640/twparams_summary.root`

---

### 4.4 Diagnostic Summary Plots (`twparams_summary.root`)

The summary file `twparams_summary.root` provides global detector Quality Assurance (QA) metrics:

* **`h_sigma_sys`** (`TH1F`): Error floor distribution ($\langle \sigma_{\text{sys}} \rangle = 40.4 \pm 17.6\text{ ps}$). *Importance:* Confirms consistent, sub-50 ps electronics stability across the full detector.
* **`h_chi2_ndf`** (`TH1F`): Reduced chi-square distribution ($\langle \chi^2/\nu \rangle = 1.19$). *Importance:* Verifies that the rational function accurately models time-walk response without under- or over-fitting.
* **`h_rms_residual`** (`TH1F`): Fractional RMS residuals ($\langle \text{RMS} \rangle = 0.81\% < 1.0\%$). *Importance:* Demonstrates sub-percent calibration accuracy across all channels.
* **`g_sigma_sys_vs_channel`** (`TGraph`): $\sigma_{\text{sys}}$ vs Channel Index (0 to 598). *Importance:* Validates uniform detector response across all planes and bars, immediately flagging any anomalous SiPMs.
* **`c_tw_summary`** (`TCanvas`): 4-pad canvas displaying all four QA metrics together.

---

### 4.5 Stage 4: Publication-Quality Validation ROOT Canvases
* **Single-Channel Validation**:
  ```bash
  python3 /afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/plot_timewalk_fit.py \
    -r 6640 \
    -p /afs/cern.ch/work/i/idioniso/sndVetoUS-physics2022/ \
    -c 24005_4
  ```
* **Full Detector Export (599 Channels)**:
  ```bash
  python3 /afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/plot_timewalk_fit.py \
    -r 6640 \
    -p /afs/cern.ch/work/i/idioniso/sndVetoUS-physics2022/ \
    -c all
  ```
* **Output**:
  `/afs/cern.ch/work/i/idioniso/sndVetoUS-physics2022/rootfiles/run006640/timewalk_fits_all.root` (Contains 599 per-channel `TDirectory` folders with 2-pad dual-error canvases).

---

### 4.6 Stage 5: Time-Walk Corrected Scintillator Speed (Mode `tw` — Active Stage)
* **Job Submission**:
  ```bash
  cd /afs/cern.ch/work/i/idioniso/sndVetoUS/htcondor
  condor_submit 0_timewalk.sub
  ```
* **Post-Processing (After Jobs Complete)**:
  ```bash
  ./merge_files.sh 6640 tw
  python3 /afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/extract_cscint.py \
    -r 6640 \
    -p /afs/cern.ch/work/i/idioniso/sndVetoUS-physics2022/ \
    -s corrected \
    --plot
  ```
* **Physical Target**: $\langle c_{\text{scint}}^{\text{corr}} \rangle \approx \mathbf{15.5\text{ cm/ns}}$ (True bulk speed of light in polystyrene).

---

## 5. Upcoming Calibration Stages

1. **Fine Alignment ($d_{\text{SiPM}}$)**:
   * Slices $[t_{\text{mode}} - 1\sigma, t_{\text{mode}} + 1\sigma]$ on $dt_{\text{SiPM}}^{\text{TW, ToF}}$ to compute truncated mean $d_{\text{SiPM}}$.
   * Saves to `Alignmentparams/run006640/alignmentparameterDS.json`.
2. **System Alignment & Multi-SiPM Averaging (Mode `systemalignment`)**:
   * Evaluates bar-side covariance matrices and single-bar time resolution ($\sigma_{\text{bar}} \approx \mathbf{245\text{ ps}}$).
