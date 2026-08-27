# Veto Timing Calibration & Alignment Guide

This document outlines the adaptation of the 4-step timing calibration pipeline to the **Veto System** in **SND@LHC**, based on the methodologies in [dissertation_conaboy_andrew.pdf](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/dissertation_conaboy_andrew.pdf) (Chapter 4).

---

## 1. Veto Detector Geometry & Geometry Upgrades

The Veto system is located upstream of the Target / SciFi tracker to tag incoming charged particles (predominantly punch-through muons from the LHC collision point).

```text
Pre-2024 Layout (Planes 0 & 1):
  ◄── LEFT (x < 0) ───[ Horizontal Bar: 1x6x42 cm³ ]─── RIGHT (x > 0) ──►
      8 Large SiPMs (0..7)                              8 Large SiPMs (8..15)

Post-2024 Upgrade (Plane 2 Added):
      ▲ TOP (y > 0) ──► 8 Large SiPMs (0..7)
      │
      │   [ Vertical Bar: 1x6x46 cm³ ]
      │   (Single-ended readout on TOP only)
      ▼
```

### Geometry Comparison: US vs. Veto

| Parameter | Upstream (US) | Veto (Planes 0 & 1: Pre-2024) | Veto (Plane 2: Added 2024) |
| :--- | :---: | :---: | :---: |
| **System ID ($s$)** | `2` | `1` | `1` |
| **Plane Indices ($p$)** | `0, 1, 2, 3, 4` (5 planes) | `0, 1` (2 planes) | `2` (1 plane) |
| **Bars per Plane ($b$)** | 10 horizontal bars | 7 horizontal bars | 7 **vertical** bars |
| **Bar Dimensions** | $1 \times 6 \times 82.5\text{ cm}^3$ | $1 \times 6 \times 42\text{ cm}^3$ | $1 \times 6 \times 46\text{ cm}^3$ |
| **Bar Orientation** | Horizontal ($x$-axis) | Horizontal ($x$-axis) | **Vertical ($y$-axis)** |
| **Readout Ends** | Both Ends (Left & Right) | Both Ends (Left & Right) | **Top End Only** |
| **SiPMs per Bar** | 16 (6 large + 2 small / side) | 16 (**8 large** / side, all identical) | **8 large** (Top only) |
| **Total Channels** | 599 channels | 224 channels | +56 channels (280 total) |

---

## 2. Adapting the 4 Calibration Steps for Veto

```text
[Veto Calibration Pipeline]
   │
   ├──► Step 1: Speed of Light (c_scint) & Photon ToF
   │      • Planes 0 & 1 (Horizontal): Linear fit dt vs x_track (42 cm range)
   │      • Plane 2 (Vertical, 2024+): Linear fit dt vs y_track (46 cm range, top readout only)
   │
   ├──► Step 2: Time-Walk Calibration (TW)
   │      • 1D modal slice + Two-stage rational fit + SciPy NLL error floor
   │      • Applied identically to all 224 (or 280) Veto channels
   │
   ├──► Step 3: Channel-by-Channel Fine Alignment (d_SiPM)
   │      • Truncated mean in ±1σ window around mode → centers all channels at 0.00 ns
   │
   └──► Step 4: Multi-SiPM Averaging & Veto Time Resolution
          • Planes 0 & 1: 8 SiPMs/side → Left/Right average → Full bar average
          • Plane 2: 8 Top SiPMs average
```

---

### 2.1 Step 1: Effective Speed of Light ($c_{\text{scint}}^{\text{Veto}}$) & Photon ToF

* **Horizontal Planes 0 & 1 ($p \in \{0, 1\}$)**:
  * Extrapolate DS tracks to $z_{\text{Veto}, p}$ to obtain $(x_{\text{track}}, y_{\text{track}})$.
  * Correlate $dt_{\text{SiPM}} = t_0^{\text{DS}} - t_{\text{SiPM}}$ against $x_{\text{track}}$ over $[-21, +21]\text{ cm}$.
  * Linear fit slope $m \implies c_{\text{scint}} = \pm 1/m$.
  * Photon ToF to bar center ($x_{\text{ref}} = 0\text{ cm}$):
    $$t_\gamma = \frac{x_{\text{track}} - x_{\text{ref}}}{c_{\text{scint}}}, \quad t_{\text{SiPM}}^{\text{ToF}} = \begin{cases} t_{\text{SiPM}} + t_\gamma & \text{for left (0..7)} \\ t_{\text{SiPM}} - t_\gamma & \text{for right (8..15)} \end{cases}$$

* **Vertical Plane 2 ($p = 2$, Post-2024 Upgrade)**:
  * Propagation delay occurs along the **$y$-axis**, not the $x$-axis.
  * Correlate $dt_{\text{SiPM}}$ against $y_{\text{track}}$ over $[-23, +23]\text{ cm}$.
  * Since readout is **only at the top** ($y = +y_{\text{top}}$), light always travels upward towards $+y$:
    $$t_\gamma = \frac{y_{\text{top}} - y_{\text{track}}}{c_{\text{scint}}}, \quad t_{\text{SiPM}}^{\text{ToF}} = t_{\text{SiPM}} + t_\gamma$$

---

### 2.2 Step 2: Time-Walk (TW) Correction

* The physical silicon response and TOFPET2 ASIC discriminator parameterisation are identical:
  $$\Delta t_{\text{TW}}(\text{QDC}) = t_0 + \frac{\alpha}{\beta \cdot \text{QDC} - \text{QDC}_0} + \gamma \cdot \text{QDC}$$
* Slices 1D timing projections along $\text{QDC}$, runs Stage 1 unweighted shape fit with Table 4.2 bounds, and Stage 2 NLL minimization for $\sigma_{\text{sys}}$.
* Saves parameters to `Polyparams/runXXXXXX/twparams_veto.json`.

---

### 2.3 Step 3: Channel-by-Channel Fine Alignment ($d_{\text{SiPM}}$)

* For each Veto channel, project the fully ToF- and TW-corrected residual:
  $$dt_{\text{SiPM}}^{\text{TW, ToF}} = t_0^{\text{DS}} - t_{\text{SiPM}}^{\text{TW, ToF}}$$
* Select the $\pm 1\sigma$ truncated window around the mode: $[t_{\text{mode}} - 1\sigma, t_{\text{mode}} + 1\sigma]$.
* Compute truncated mean $d_{\text{SiPM}}$ and save to `Alignmentparams/runXXXXXX/alignmentparams_veto.json`.
* Event-wise subtraction centers all Veto channels at $0.00\text{ ns}$:
  $$dt_{\text{SiPM}}^{\text{aligned}} = t_0^{\text{DS}} - t_{\text{SiPM}}^{\text{TW, ToF}} - d_{\text{SiPM}} \equiv 0.00\text{ ns}$$

---

### 2.4 Step 4: Multi-SiPM Averaging & Veto Time Resolution

* In Veto, **all 8 SiPMs per side are large ($6 \times 6\text{ mm}^2$)**, so all 8 participate in the timing average:
  $$dt_{\text{side}} = \frac{1}{8} \sum_{i=0}^7 dt_{\text{SiPM}, i}^{\text{aligned}}$$
* Expected Resolutions:
  * **Single SiPM**: $\sigma_t \approx \mathbf{320\text{ ps}}$
  * **Bar-Side Average** (8 SiPMs): $\sigma_{\text{side}} \approx \mathbf{230\text{ ps}}$ (improved over US due to 8 active channels)
  * **Full Bar Average** (Left + Right): $\sigma_{\text{bar}} \approx \mathbf{200\text{ ps}}$

---

## 3. Required Codebase Adaptations in `sndsw`

1. **[`AnalysisFunctions.py`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/python/AnalysisFunctions.py)**:
   * Update system dictionaries for Veto (`s = 1`):
     ```python
     self.systemAndPlanes = {1: 2, 2: 5, 3: 7}  # {1: 3} for post-2024 data
     self.systemAndBars = {1: 7, 2: 10, 3: 60}
     self.systemAndSiPMs = {
         1: list(range(16)),  # 0..15 for planes 0,1; 0..7 for plane 2
         2: [0, 1, 3, 4, 6, 7, 8, 9, 11, 12, 14, 15],  # US large SiPMs
     }
```
   * Ensure `correct_ToF(...)` computes $y$-propagation for Veto Plane 2.
2. **Extraction Scripts (`extract_cscint.py`, `extract_twparams.py`, `extract_alignment.py`)**:
   * Add `--system {veto, us, ds}` command-line arguments to select the target sub-detector and geometry bounds.
