# Downstream (DS) Muon Filter Timing Calibration Guide

This document outlines the architecture, mathematical formulations, and adaptation plan to perform full timing calibration and alignment on the **Downstream (DS) Muon Filter** in **SND@LHC**, based on [dissertation_conaboy_andrew.pdf](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/dissertation_conaboy_andrew.pdf) (Chapter 4 and Appendix A.6).

---

## 1. DS Detector Architecture & Channel Mapping

The Downstream Muon Filter consists of **4 detector stations** interleaved between iron absorber blocks. 
* Stations 1–3 contain both a **Horizontal (DSH)** and a **Vertical (DSV)** plane.
* Station 4 contains only a **Vertical (DSV)** plane.

```text
DS Station Layout:
  • Station 1: Plane 0 (DSH 1: Horizontal)  +  Plane 1 (DSV 1: Vertical)
  • Station 2: Plane 2 (DSH 2: Horizontal)  +  Plane 3 (DSV 2: Vertical)
  • Station 3: Plane 4 (DSH 3: Horizontal)  +  Plane 5 (DSV 3: Vertical)
  • Station 4: Plane 6 (DSV 4: Vertical)

Readout Topologies:
  • DSH (Horizontal): ◄── LEFT (SiPM 0) ───[ 1x1x82.5 cm³ ]─── RIGHT (SiPM 1) ──►  (Double-ended)
  • DSV (Vertical):   ▲ TOP (SiPM 0) ──► [ 1x1x62.5 cm³ ]                          (Top single-ended)
```

### Channel Count Breakdown

| Subsystem Component | Planes | Bars per Plane | Bar Dimensions | Readout Ends | SiPMs/Bar | Channel Count |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **DSH (Horizontal)** | 3 (`p = 0, 2, 4`) | 60 | $1 \times 1 \times 82.5\text{ cm}^3$ | Both (Left & Right) | 2 (1 per side) | $3 \times 60 \times 2 = \mathbf{360}$ |
| **DSV (Vertical)** | 4 (`p = 1, 3, 5, 6`) | 60 | $1 \times 1 \times 62.5\text{ cm}^3$ | Top End Only | 1 (Top) | $4 \times 60 \times 1 = \mathbf{240}$ |
| **Total DS System** | **7 planes** | — | — | — | — | **$\mathbf{600\text{ channels}}$** |

---

## 2. Self-Consistent Reference Clock for DS Alignment

In the US calibration, $t_0^{\text{DS}}$ was the average of the 3 horizontal DS planes. To align the DS system itself without circular bias:

1. **Leave-One-Out Reference for Horizontal Planes ($p_i \in \{0, 2, 4\}$)**:
   When calibrating DSH plane $p_i$, build the reference clock using the **other two DSH planes** ($j \neq i$):
   $$t_{0, -i}^{\text{DS}} = \frac{1}{2} \sum_{j \neq i} t_{\text{DSH}, j}, \quad \text{where } t_{\text{DSH}, j} = \frac{t_{j, \text{left}} + t_{j, \text{right}}}{2}$$
2. **Reference Clock for Vertical Planes ($p \in \{1, 3, 5, 6\}$)**:
   Use the full 3-plane average:
   $$t_0^{\text{DS}} = \frac{1}{3} \left( t_{\text{DSH}, 1} + t_{\text{DSH}, 2} + t_{\text{DSH}, 3} \right)$$

---

## 3. Adapting the 4 Calibration Steps for DS

```text
[DS Calibration Pipeline]
   │
   ├──► Step 1: Speed of Light (c_scint) & Photon Propagation ToF
   │      • DSH (Planes 0, 2, 4): Linear fit dt vs x_track (82.5 cm range, left/right readout)
   │      • DSV (Planes 1, 3, 5, 6): Linear fit dt vs y_track (62.5 cm range, top readout)
   │
   ├──► Step 2: Time-Walk (TW) Calibration
   │      • 1D modal slice + Two-stage rational fit + SciPy NLL error floor
   │      • Extracted for all 600 single-SiPM channels
   │
   ├──► Step 3: Channel-by-Channel Fine Alignment (d_SiPM)
   │      • Truncated mean in ±1σ window around mode → centers all channels at 0.00 ns
   │
   └──► Step 4: Full Bar Averaging & DS Resolution
          • DSH: (t_left + t_right)/2 → Automatic position cancellation (σ_bar ≈ 263 ps)
          • Combined DS Reference: σ(t₀ᴰˢ) = 263 ps / √3 ≈ 150 ps
```

---

### 3.1 Step 1: Effective Speed of Light ($c_{\text{scint}}^{\text{DS}}$) & Photon ToF

* **Horizontal Planes (DSH, $p \in \{0, 2, 4\}$)**:
  * Correlate $dt_{\text{SiPM}} = t_{0, -i}^{\text{DS}} - t_{\text{SiPM}}$ against $x_{\text{track}}$ over $[-40, +40]\text{ cm}$.
  * Slope gives $c_{\text{scint}} = \pm 1/m$.
  * Photon ToF to bar center ($x_{\text{ref}} = 0\text{ cm}$):

  $$t_\gamma = \frac{x_{\text{track}} - x_{\text{ref}}}{c_{\text{scint}}}, \quad t_{\text{SiPM}}^{\text{ToF}} = \begin{cases} t_{\text{SiPM}} + t_\gamma & \text{for Left (SiPM 0)} \\ t_{\text{SiPM}} - t_\gamma & \text{for Right (SiPM 1)} \end{cases}$$

* **Vertical Planes (DSV, $p \in \{1, 3, 5, 6\}$)**:
  * Propagation delay occurs along the **$y$-axis** over $[-30, +30]\text{ cm}$.
  * Since readout is **only on top** ($y = +y_{\text{top}}$), light travels upward towards $+y$:
    $$t_\gamma = \frac{y_{\text{top}} - y_{\text{track}}}{c_{\text{scint}}}, \quad t_{\text{SiPM}}^{\text{ToF}} = t_{\text{SiPM}} + t_\gamma$$

---

### 3.2 Step 2: Time-Walk (TW) Correction

* Fitted with the 5-parameter rational function:
  $$\Delta t_{\text{TW}}(\text{QDC}) = t_0 + \frac{\alpha}{\beta \cdot \text{QDC} - \text{QDC}_0} + \gamma \cdot \text{QDC}$$
* Slices 1D timing projections along $\text{QDC}$, fits unweighted shape, and runs SciPy NLL error floor extraction.
* Output: `Polyparams/runXXXXXX/twparams_ds.json` (600 channels).

---

### 3.3 Step 3: Channel-by-Channel Fine Alignment ($d_{\text{SiPM}}$)

* For each DS channel, project the residual:
  $$dt_{\text{SiPM}}^{\text{TW, ToF}} = t_0^{\text{DS}} - t_{\text{SiPM}}^{\text{TW, ToF}}$$
* Select the $\pm 1\sigma$ truncated window around the mode: $[t_{\text{mode}} - 1\sigma, t_{\text{mode}} + 1\sigma]$.
* Alignment parameter $d_{\text{SiPM}} = \text{Truncated Mean}$.
* Event-wise subtraction centers all 600 DS channels at $0.00\text{ ns}$:
  $$dt_{\text{SiPM}}^{\text{aligned}} = t_0^{\text{DS}} - t_{\text{SiPM}}^{\text{TW, ToF}} - d_{\text{SiPM}} \equiv 0.00\text{ ns}$$

---

### 3.4 Step 4: Full Bar Averaging & DS Reference Resolution

* **DSH Horizontal Bar Average** (Equation 4.1.1):
  $$t_{\text{DSH}} = \frac{1}{2} \left( t_{\text{left}}^{\text{aligned}} + t_{\text{right}}^{\text{aligned}} \right)$$
  * **Automatic Position Cancellation:** Scintillation light propagation along the bar cancels out exactly:
    $$\frac{(t_{\text{hit}} + x/c) + (t_{\text{hit}} + (L-x)/c)}{2} = t_{\text{hit}} + \frac{L}{2c}$$
  * Single DSH Bar Resolution: $\sigma(t_{\text{DSH}}) \approx \mathbf{263\text{ ps}}$.
* **3-Plane DS Global Reference Time** ($t_0^{\text{DS}}$):
  $$t_0^{\text{DS}} = \frac{1}{3} \sum_{i=1}^3 t_{\text{DSH}, i} \implies \sigma(t_0^{\text{DS}}) \approx \frac{263\text{ ps}}{\sqrt{3}} \approx \mathbf{150\text{ ps}}$$

---

## 4. Codebase Adaptations in `sndsw`

1. **[`AnalysisFunctions.py`](file:///afs/cern.ch/work/i/idioniso/sndVetoUS/sndsw/python/AnalysisFunctions.py)**:
   * Define DS geometry dictionary (`s = 3`):
     ```python
     self.systemAndPlanes = {1: 2, 2: 5, 3: 7}
     self.systemAndBars = {1: 7, 2: 10, 3: 60}
     self.systemAndSiPMs = {
         1: list(range(16)),
         2: [0, 1, 3, 4, 6, 7, 8, 9, 11, 12, 14, 15],
         3: [0, 1],  # 0,1 for DSH; 0 for DSV
     }
```
   * Add `GetDSH_leave_one_out(excluded_plane)` for unbiased DSH calibration.
2. **Extraction Scripts**:
   * Add `--system ds` flag to handle 60 bars $\times$ 7 planes and select $x$ (DSH) vs $y$ (DSV) track correlations.
