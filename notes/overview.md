# Upstream Scintillator (US) Timing Calibration & Alignment Guide

This document provides a comprehensive reference of the general timing calibration, time-walk (TW) correction, and channel alignment procedure for the Upstream Scintillator (US) system in **SND@LHC**, based on [dissertation_conaboy_andrew.pdf](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/dissertation_conaboy_andrew.pdf) (Chapter 4) and implemented in the [`sndsw`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw) framework.

---

## 1. Overview of the Calibration Workflow

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

> **Note on $t_0$:** The parameter $t_0$ is discarded when applying the TW correction because absolute offsets are calibrated in the fine alignment step.

---

## 5. Step 3: Channel-by-Channel Fine Alignment ($d_{SiPM}$)

### 5.1 Physics Motivation (Section 4.1.7)
Even after ToF and TW corrections, individual SiPM channels have static time offsets due to differences in readout cable lengths, TOFPET ASIC internal propagation delays, SiPM transit times, and DAQ clock distribution offsets.

### 5.2 Alignment Parameter Extraction ($d_{SiPM}$)
1. Form the 1D projection of fully corrected residual: $dt_{SiPM}^{TW, ToF} = t_0^{DS} - t_{SiPM}^{TW, ToF}$.
2. Select $\pm 1\sigma$ truncated window around the mode: $[t_{\text{mode}} - 1\sigma, t_{\text{mode}} + 1\sigma]$ (removes low-QDC negative noise tails).
3. Compute the **truncated mean** as $d_{SiPM}$ with uncertainty $\text{SEM} = \frac{\sigma_{\text{trunc}}}{\sqrt{N_{\text{trunc}}}}$.
4. Correct the event-wise time:
   $$dt_{SiPM}^{aligned} = t_0^{DS} - t_{SiPM}^{TW, ToF} - d_{SiPM} \equiv 0.00\text{ ns}$$

---

## 6. Step 4: Multi-SiPM Averaging & System Time Resolution

### 6.1 Channel Selection & Bar-Side Average
* Each US bar end has 6 large ($6 \times 6\text{ mm}^2$) and 2 small ($1 \times 1\text{ mm}^2$) SiPMs. Small SiPMs are excluded from timing averaging due to lower light collection.
* **Bar-Side Average** (Equation 4.1.13):
  $$dt_{side}^{ToF} = \frac{1}{6} \sum_{i \in side} dt_{\text{SiPM}, i}^{\text{aligned}}$$
  * Single SiPM resolution: $\sigma \approx 313\text{ ps}$
  * Bar-side resolution: $\sigma \approx 261\text{ ps}$ *(limited by PCB cross-talk / covariance)*

### 6.2 Full Bar Average (Equation 4.3.11)
$$dt_{bar}^{ToF} = \frac{1}{2} \left( dt_{left}^{ToF} + dt_{right}^{ToF} \right)$$
* Full bar resolution: $\sigma \approx \mathbf{245\text{ ps}}$.
* **Position Cancellation:** When averaging left and right sides without track extrapolation ($x$ unknown), signal ToF along the bar cancels out automatically:
  $$t_{left} \propto \frac{x}{c}, \quad t_{right} \propto \frac{L-x}{c} \implies \frac{t_{left} + t_{right}}{2} \propto \frac{L}{2c}$$

---

## 7. Scripts & Processing Utilities Reference

| Script / File | Location | Primary Purpose & Features |
| :--- | :--- | :--- |
| **`generate_args_timewalk.py`** | [`htcondor/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/htcondor/) | Generates HTCondor job argument text files (`args_timewalk_<mode>.txt`) splitting $N$ events into chunked ranges (`run, start, nevents, mode`). |
| **`run_timewalk.sh`** | [`htcondor/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/htcondor/) | HTCondor executable wrapper; sets up CVMFS environment, runs `run_TimeWalk.py`, and transfers chunk outputs. |
| **`0_timewalk.sub`** | [`htcondor/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/htcondor/) | HTCondor submission configuration for batch job queuing. |
| **`merge_files.sh`** | [`htcondor/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/htcondor/) | Multi-threaded merging script (16 parallel workers) that merges 600 channel split ROOT files into consolidated `timewalk_<ch>.root` files. |
| **`extract_cscint.py`** | [`shipLHC/scripts/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/) | Extracts effective speed of light ($c_{\text{scint}}$) from $dt$ vs $x_{\text{pred}}$ correlations; fits linear slope $m = \pm 1/c_{\text{scint}}$, computes 5-segment systematic errors, and outputs JSON/ROOT summaries. Supports `-s uncorrected` and `-s corrected`. |
| **`extract_twparams.py`** | [`shipLHC/scripts/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/) | Complete 599-channel time-walk calibration pipeline; projects 1D modal timing slices, runs Stage 1 unweighted least squares fit, runs Stage 2 SciPy NLL error floor extraction, and exports `twparams.json` and per-channel JSON/CSV files. |
| **`plot_timewalk_fit.py`** | [`shipLHC/scripts/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/shipLHC/scripts/) | Dissertation Fig. 4.7 reproducer and ROOT diagnostic generator; draws upper 2D histogram + statistical points + shaded total error envelope ($\text{SEM} \oplus \sigma_{\text{sys}}$) + fit curve, lower Data/Fit ratio pane, displays $\chi^2_{\text{stat}}/\nu \rightarrow \chi^2_{\text{tot}}/\nu$, and supports `-c all` export into per-channel `TDirectory` folders. |
| **`AnalysisFunctions.py`** | [`python/`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/python/) | Core framework library updated with run-aware `MakeTWCorrectionDict` and `MuFilterCorrectedTime` for event-by-event time corrections. |

---

## 8. Related Run Notes & Subsystems

* **[US Run 6640 Calibration Log](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/notes/US_run6640.md)**: Concrete execution steps, commands, and physics results for Run 6640.
* **[Veto Timing Calibration Plan](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/notes/veto.md)**: Adaptation plan for the Veto system (pre-2024 horizontal planes and post-2024 vertical plane).
