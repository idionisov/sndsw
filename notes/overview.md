# Upstream Scintillator (US) Timing Calibration & Alignment Guide

This document provides a comprehensive overview of the timing calibration, time-walk (TW) correction, and channel alignment procedure for the Upstream Scintillator (US) system in **SND@LHC**, as developed in [dissertation_conaboy_andrew.pdf](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/dissertation_conaboy_andrew.pdf) (Chapter 4) and implemented in the [`sndsw`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw) framework.

---

## 1. Overview of the Procedure Workflow

The calibration transforms raw SiPM TDC timestamps into synchronized physics times via a 4-step sequence:

```text
[Raw SiPM Hits & Tracks]
        │
        ├──► Step 0: Reference Event Time Definition (t₀ᴰˢ or t₀ˢᶠ)
        │
        ├──► Step 1: Effective Speed of Light (c_SiPM) ──► Photon ToF Correction
        │            (Removes position-dependent propagation delay along bar)
        │
        ├──► Step 2: Time-Walk (TW) Correction (Two-Stage Fit + NLL Minimisation)
        │            (Removes amplitude-dependent discriminator delay)
        │
        ├──► Step 3: Channel-by-Channel Fine Alignment (d_SiPM)
        │            (Absorbs cable lengths, ASIC delays, and clock skews to global t₀)
        │
        └──► Step 4: Multi-SiPM Averaging (Side, Bar, Plane, US System)
                     (Combines 6 large SiPMs/side; accounts for PCB covariance)
```

---

## 2. Step 0: Reference Event Time ($t_0$)

To study single-channel timing responses without event-to-event LHC bunch jitter or global start time offsets, an event-wise reference time $t_0$ is defined:

### 2.1 Downstream Reference Time ($t_0^{DS}$) — Equation 4.1.1
Used for single-muon events (`referencesystem = 3`):
$$t_0^{DS} = \frac{1}{N_{DSH}} \sum_{i=1}^{N_{DSH}} t_{DSH, i}$$
where each horizontal DS bar average is:
$$t_{DSH, i} = \frac{1}{2} (t_{i, \text{left}} + t_{i, \text{right}})$$

> **Geometric Cancellation:** Because scintillation light travels to both ends with effective speed $c$, the average $t_{DSH} = \frac{1}{2} [ (t_{hit} + x/c) + (t_{hit} + (L-x)/c) ] = t_{hit} + \frac{L}{2c}$ is completely **independent of the horizontal hit position $x$**. This produces a clean reference time with intrinsic resolution $\sigma(t_0^{DS}) \approx 263\text{ ps}$.

### 2.2 SciFi Reference Time ($t_0^{SF}$) — Equation 4.1.2
Used when no muon reaches the DS (e.g., zero-muon neutrino interactions, `referencesystem = 1`):
$$t_0^{SF} = \frac{1}{N} \sum_{i=1}^{N} t_{SiPM, i}^{SF}$$

For every US SiPM channel, raw relative times are defined as:
$$dt_{SiPM} = t_0^{DS} - t_{SiPM}$$

---

## 3. Step 1: Effective Speed of Light Determination ($c_{SiPM}$)

### 3.1 Physics Motivation
Scintillation photons emitted at track position $x$ travel along the bar before reaching the SiPM. Across an $82.5\text{ cm}$ bar, propagation time varies by up to $\sim 6\text{ ns}$. If this is not corrected first, it washes out the amplitude-dependent time-walk effect.

### 3.2 Methodology (Section 4.1.3)
1. Extrapolate reconstructed DS tracks to the $z$-position of the US plane to obtain $(x_{track}, y_{track})$.
2. Correlate raw relative time $dt_{SiPM} = t_0^{DS} - t_{SiPM}$ against $x_{track}$ in a 2D histogram.
3. For each $x$-slice, extract the Most Probable Value (MPV) of $dt_{SiPM}$.
4. Fit a linear relationship (Equation 4.1.4):
   $$dt_{SiPM} = \pm \frac{x}{c_{SiPM}} + \text{const}$$
   *(Positive slope for left-side SiPMs, negative slope for right-side SiPMs due to coordinate system orientation).*
5. Extract $c_{SiPM} \approx 14\text{--}16\text{ cm/ns}$.
6. **Systematic Uncertainty:** The fit range $[A - 5\text{ cm}, L_\mu]$ is divided into 5 equal segments; the variance across sub-measurements gives the systematic uncertainty on $c_{SiPM}$ (Equation 4.1.6).

### 3.3 Photon Time-of-Flight (ToF) Correction (Equations 4.1.7 & 4.1.8)
Correct the time to the reference bar center ($x_{ref} = \frac{A_x + B_x}{2}$):
$$t_\gamma = \frac{x_{track} - x_{ref}}{c_{SiPM}}$$
$$t_{SiPM}^{ToF} = \begin{cases} t_{SiPM} + t_\gamma & \text{for left SiPMs} \\ t_{SiPM} - t_\gamma & \text{for right SiPMs} \end{cases}$$
$$dt_{SiPM}^{ToF} = t_0^{DS} - t_{SiPM}^{ToF}$$

---

## 4. Step 2: Time-Walk (TW) Correction & Two-Stage Fitting

### 4.1 Physics Motivation
Due to the leading-edge discriminator in the TOFPET2 ASIC, smaller signals (low QDC) cross the threshold later than larger signals (high QDC).

### 4.2 Time-Walk Parameterisation (Equation 4.1.9)
$$\Delta t_{TW}(QDC) = t_0 + \frac{\alpha}{\beta \cdot QDC - QDC_0} + \gamma \cdot QDC$$

* Parameter bounds established from Table 4.2:
  * $t_0 \in [0.5, 1.5] \times t_0^{\text{estimate}}$
  * $\alpha \in [-500, 0]$ (negative amplitude bend)
  * $\beta \in [0, 30]$
  * $QDC_0 \in [-10, 0]$
  * $\gamma \in [0, 1]$

### 4.3 Why Two-Stage Fitting? (Section 4.3.4)
Because the QDC follows a convolved Landau-Gaussian distribution, high-density bins near the MPV have huge statistics and tiny standard errors ($\sigma_{stat} = \text{SEM} \sim \text{few ps}$). A direct $\chi^2$ fit results in unphysically massive $\chi^2/\text{ndf} \sim 100\text{--}2500$ because systematic variations dominate.

1. **Stage 1 (Shape Extraction):**
   * Slices the 2D correlation into 1D timing projections (`ProjectionY`) for each QDC bin, taking the **modal value (peak bin)** as $y_i$ and the Standard Error on the Mean ($\text{SEM} = \frac{\text{RMS}}{\sqrt{N}}$) as statistical weight.
   * Fits the function unweighted (minimising direct residuals, ignoring statistical error bars) over $[0, QDC_{max}]$.
   * Establishes the underlying physical parameters $(\alpha, \beta, QDC_0, \gamma, t_0)$.
2. **Stage 2 (MLE for Systematic Error Floor):**
   * Construct the Negative Log-Likelihood (NLL, Equation 4.3.6) assuming total error $\sigma^2 = \sigma_{stat}^2 + \sigma_{sys}^2$:
     $$\text{NLL}(\theta) = \frac{1}{2} \sum_{i=1}^{N} \left[ \frac{(y_i - f(x_i; \theta))^2}{\sigma_{stat, i}^2 + \sigma_{sys}^2} + \ln(\sigma_{stat, i}^2 + \sigma_{sys}^2) \right]$$
   * Numerically minimise NLL using SciPy to find the optimal error floor $\sigma_{sys}$ ($\approx 40.4\text{ ps}$).
   * Re-evaluate $\chi^2/\text{ndf}$ with total errors $\sigma_i = \sqrt{\sigma_{stat, i}^2 + \sigma_{sys}^2}$, yielding proper $\chi^2/\text{ndf} \approx 1.19$.

### 4.4 Applying the TW Correction (Equation 4.1.11 / [`AnalysisFunctions.py`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/python/AnalysisFunctions.py#L53))
$$t_{SiPM}^{TW} = t_{SiPM} + \left( \frac{\alpha}{\beta \cdot QDC - QDC_0} + \gamma \cdot QDC \right)$$

---

## 5. Step 3: Channel-by-Channel Fine Alignment ($d_{SiPM}$)

### 5.1 Physics Motivation (Section 4.1.7)
Even after ToF and TW corrections, individual SiPM channels have static time offsets due to differences in readout cable lengths, TOFPET ASIC internal propagation delays, SiPM transit times, and DAQ clock distribution offsets.

### 5.2 Alignment Parameter Extraction ($d_{SiPM}$)
1. Form the 1D projection of fully corrected residual: $dt_{SiPM}^{TW, ToF} = t_0^{DS} - t_{SiPM}^{TW, ToF}$.
2. Compute the **truncated mean** within $\pm 1\sigma$ around the modal value.
3. The alignment constant $d_{SiPM}$ is this truncated mean.
4. Correct the event-wise time:
   $$dt_{SiPM}^{aligned} = t_0^{DS} - t_{SiPM}^{TW, ToF} - d_{SiPM} \approx 0\text{ ns}$$

---

## 6. Scripts & Processing Utilities Reference

The following table summarizes all custom scripts developed in this repository to automate the US timing calibration pipeline:

| Script / File | Location | Primary Purpose & Features |
| :--- | :--- | :--- |
| **`generate_args_timewalk.py`** | [`htcondor/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/htcondor/) | Generates HTCondor job argument text files (`args_timewalk_<mode>.txt`) splitting $N$ events into chunked ranges (`run, start, nevents, mode`). |
| **`run_timewalk.sh`** | [`htcondor/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/htcondor/) | HTCondor executable wrapper; sets up the CVMFS environment, executes `run_TimeWalk.py`, and transfers chunk outputs. |
| **`0_timewalk.sub`** | [`htcondor/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/htcondor/) | HTCondor submission configuration for batch job queuing. |
| **`merge_files.sh`** | [`htcondor/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/htcondor/) | Multi-threaded merging script (16 parallel workers) that merges 600 channel split ROOT files into consolidated `timewalk_<ch>.root` files. |
| **`extract_cscint.py`** | [`shipLHC/scripts/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/) | Extracts effective speed of light ($c_{\text{scint}}$) from $dt$ vs $x_{\text{pred}}$ correlations; fits linear slope $m = \pm 1/c_{\text{scint}}$, computes 5-segment systematic errors, and outputs JSON/ROOT summaries. Supports `-s uncorrected` and `-s corrected`. |
| **`extract_twparams.py`** | [`shipLHC/scripts/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/) | Complete 599-channel time-walk calibration pipeline; projects 1D modal timing slices, runs Stage 1 unweighted least squares fit, runs Stage 2 SciPy NLL error floor extraction, and exports `twparams.json` and per-channel JSON/CSV files. |
| **`plot_timewalk_fit.py`** | [`shipLHC/scripts/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/) | Dissertation Fig. 4.7 reproducer and ROOT diagnostic generator; draws upper 2D histogram + statistical points + shaded total error envelope ($\text{SEM} \oplus \sigma_{\text{sys}}$) + fit curve, lower Data/Fit ratio pane, displays $\chi^2_{\text{stat}}/\nu \rightarrow \chi^2_{\text{tot}}/\nu$, and supports `-c all` export into per-channel `TDirectory` folders. |
| **`AnalysisFunctions.py`** | [`python/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/python/) | Core framework library updated with run-aware `MakeTWCorrectionDict` and `MuFilterCorrectedTime` for event-by-event time corrections. |

---

## 7. Concrete Operational Pipeline & Executed Workflow (Run 6640)

Below is the chronological execution pipeline for **Run 6640** (Physics 2022 dataset):

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

### Step-by-Step Executed Commands

#### 1. Uncorrected Scintillator Speed Extraction (Mode `zeroth`)
* **Job Submission**:
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
* **Result**:
  * Total Fitted Channels: **599 / 599**
  * Mean Scintillator Speed: **$\langle c_{\text{scint}} \rangle = 14.32 \pm 0.62\text{ cm/ns}$** (Left: $14.24\text{ cm/ns}$, Right: $14.40\text{ cm/ns}$)
  * Files: `cscintvalues/run006640/run006640_cscintvalues_uncorrected.json`, `cscint_summary_uncorrected.root`.

---

#### 2. Photon Propagation Correction (Mode `tof`)
* **Job Submission & Merging**:
  ```bash
  python3 generate_args_timewalk.py -r 6640 -m tof -n 10000000 -c 200000 -o args_timewalk_tof.txt
  condor_submit 0_timewalk.sub
  ./merge_files.sh 6640 tof
  ```
* **Output**:
  Generates `dtvqdc_{fixed_ch}_uncorrected` 2D histograms with photon ToF subtracted.

---

#### 3. Time-Walk Calibration & NLL Optimization
* **Extraction Command**:
  ```bash
  python3 /afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/extract_twparams.py \
    -r 6640 \
    -p /afs/cern.ch/work/i/idioniso/sndVetoUS-physics2022/
  ```
* **Result**:
  * Total Channels Fitted: **599 / 599**
  * Mean Systematic Error Floor: **$\langle \sigma_{\text{sys}} \rangle = 40.4 \pm 17.6\text{ ps}$** (Target: $\approx 40\text{ ps}$, [thesis §4.1.4](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/dissertation_text.txt#L5090))
  * Mean $\chi^2/\nu$: **$1.19$** (Target: $\approx 1.0$)
  * Mean RMS Residual: **$0.81\%$** (Target: $< 1.0\%$)
  * Files: `Polyparams/run006640/twparams.json`, `polyparams9_{fixed_ch}.json`, `twparams_summary.root`.

---

#### 4. Validation ROOT Canvases
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

#### 5. Time-Walk Corrected Scintillator Speed (Mode `tw` — Active Stage)
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
* **Physical Target**: $\langle c_{\text{scint}}^{\text{corr}} \rangle \approx \mathbf{15.5\text{ cm/ns}}$ ([thesis Fig. 4.26](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/dissertation_text.txt#L6591)).
