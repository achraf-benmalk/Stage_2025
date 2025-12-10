# PE Pipe Degradation Modeling: Research Summary

## Project Overview

**Title:** Predictive Modeling of Polyethylene Pipe Degradation Under Chlorine Exposure

**Institution:** Laboratoire PIMM (Procédés et Ingénierie en Mécanique et Matériaux), Arts et Métiers ParisTech

**Industrial Partner:** SUEZ Environnement

**Supervisor:** Prof. Xavier Colin (Leading expert on PE degradation kinetics)

**Duration:** Stage 2025 (6 months)

---

## Executive Summary

This research project implemented and validated a 1D diffusion-reaction kinetic model for predicting polyethylene (PE100) pipe degradation under chlorinated water exposure. The model is based on the Colin et al. (2009) mechanistic scheme, which couples radical chain oxidation kinetics with diffusion processes to predict antioxidant consumption, molecular weight evolution, and ultimately, pipe lifetime.

### Key Accomplishments

1. **Model Implementation:** Successfully implemented the complete Colin et al. kinetic scheme in Python with 12 chemical species and diffusion-reaction coupling

2. **Systematic Sensitivity Analysis:** Quantified the relative importance of key parameters:
   - DOC concentration: **DOMINANT** (ΔOIt = 278 min range at 6 months)
   - Temperature: Moderate impact (ΔOIt = 3 min range)
   - Initial antioxidant: Minor impact (ΔOIt = 0.2 min range)

3. **Validation Against Industrial Data:** Compared model predictions with SUEZ PE100 experimental data:
   - HOCl exposure: R² = 0.97 (captures degradation trend)
   - H2O control: Identified model limitation (physical aging under-predicted)

4. **Critical Limitation Identified:** The DOC→HOCl kinetics approximation is a key source of uncertainty requiring further investigation

---

## Technical Details

### Model Description

The degradation model implements a 1D diffusion-reaction system following Colin et al. (2009):

**Chemical Species Tracked (12):**
- Oxygen (O2) - diffusing
- Chlorine dioxide (DOC) - diffusing
- Antioxidant (AH) - diffusing
- Alkyl radicals (P•)
- Peroxy radicals (PO2•)
- Hydroperoxides (POOH)
- Polymer (PH)
- Caged radicals (Q)
- Carbonyls (CO)
- Grafted chlorine (PCl)
- Chain scissions (S)
- Crosslinks (X)

**Key Kinetic Equations:**

```
DOC Initiation:    DOC + PH → P• + products         (k1d)
POOH Decomposition: POOH → 2P• + chain scission     (k1u)
Propagation:       P• + O2 → PO2•                    (k2)
                   PO2• + PH → POOH + P•             (k3)
Stabilization:     PO2• + AH → products             (k7)
```

**Diffusion-Reaction Coupling:**

$$\frac{\partial C_i}{\partial t} = D_i \frac{\partial^2 C_i}{\partial z^2} + R_i(C_1, ..., C_n)$$

Where $D_i$ is the diffusion coefficient and $R_i$ is the reaction rate for species $i$.

### Validation Conditions

**SUEZ PE100 Film Experiments:**
- Material: PE100 (HDPE), 400 μm thickness
- Temperature: 40°C
- Exposure: H2O (control) and HOCl (0.05 ppm)
- Duration: 0-9 months
- Measurement: OIT at 200°C (minutes)

### Parameter Sensitivity Results

| Parameter | Range Tested | Impact on OIT at 6 months |
|-----------|--------------|---------------------------|
| DOC concentration | 0 - 5 ppm | 278 min (DOMINANT) |
| Temperature | 30 - 60°C | 3 min |
| Initial [AH]₀ | 0.5× - 3× | 0.2 min |

**Key Finding:** DOC concentration is the dominant factor affecting degradation rate, with temperature and antioxidant content having secondary effects.

### Temperature Acceleration Factors

| Temperature | Acceleration Factor (vs 30°C) |
|-------------|------------------------------|
| 30°C | 1.00× (reference) |
| 40°C | 1.36× |
| 50°C | 1.85× |
| 60°C | 2.37× |

**Arrhenius Analysis:**
- Apparent Activation Energy: 24.4 kJ/mol
- Literature Range for PE: 80-120 kJ/mol
- **Discrepancy indicates diffusion-controlled regime at low DOC concentrations**

---

## Key Limitations Identified

### 1. DOC→HOCl Kinetics Approximation

**The Problem:**
The model uses rate constants calibrated for DOC (chlorine dioxide), but SUEZ pipes are exposed to HOCl (hypochlorous acid). These have fundamentally different oxidation mechanisms:

- **DOC** (•O-Cl=O): Free radical that directly abstracts hydrogen from PE
- **HOCl**: Non-radical that dissociates to produce HO• and Cl• radicals

**Evidence:**
- Model predicts correct degradation RATE for HOCl (R² = 0.97)
- But the underlying mechanism may be incorrect
- Activation energy (24 kJ/mol) much lower than literature (80-120 kJ/mol)

**Proposed Path Forward:**
1. Literature review on HOCl oxidation kinetics
2. Identify effective DOC equivalent for HOCl exposure
3. Consider explicit HOCl dissociation step in mechanism

### 2. Physical Aging Under-Prediction

**The Problem:**
The H2O control simulation under-predicts OIT decay (R² = -2.29):
- Experimental: OIT drops from 271 → 162 min in 9 months
- Simulated: OIT drops from 289 → 270 min in 9 months

**Root Cause:**
Boundary extraction coefficient (β₀) appears too low for SUEZ film conditions.

**Implications:**
Physical aging contribution to total degradation may be underestimated.

---

## Deliverables

### Code Modules

| File | Description |
|------|-------------|
| `src/DegradationModel.py` | Clean 1D diffusion-reaction model (650+ lines) |
| `src/ParameterAnalysis.py` | Sensitivity analysis framework |
| `src/ValidationMetrics.py` | Statistical validation tools |
| `run_sensitivity_analysis.py` | Complete analysis pipeline |

### Generated Figures

| Figure | Content |
|--------|---------|
| `doc_sensitivity.png` | DOC concentration parameter sweep |
| `temperature_sensitivity.png` | Temperature parameter sweep |
| `antioxidant_sensitivity.png` | Initial AH content sweep |
| `arrhenius_plot.png` | Arrhenius analysis with Ea extraction |
| `tornado_diagram.png` | Parameter importance ranking |
| `validation_comparison.png` | Model vs SUEZ experimental data |

### Data Files

| File | Content |
|------|---------|
| `data/doc_sweep_results.csv` | DOC sweep simulation results |
| `data/temperature_sweep_results.csv` | Temperature sweep results |
| `data/antioxidant_sweep_results.csv` | AH sweep results |
| `data/arrhenius_analysis.csv` | Arrhenius analysis data |

---

## Interview-Ready Talking Points

### On Scientific Rigor

> "I implemented the complete Colin et al. kinetic scheme with 12 chemical species and diffusion-reaction coupling. The model validates well against SUEZ experimental data (R² = 0.97 for HOCl aging), though I identified a key limitation in the DOC-to-HOCl kinetics approximation that affects activation energy predictions."

### On Problem Identification

> "Through systematic sensitivity analysis, I identified that DOC concentration is the dominant factor affecting degradation rate - it's roughly 100× more impactful than temperature or antioxidant content at typical operating conditions. This has important implications for water treatment optimization."

### On Model Limitations

> "The model uses DOC kinetics as a proxy for HOCl, which introduces uncertainty. DOC is a free radical that directly attacks PE, while HOCl dissociates to produce HO• and Cl• radicals. I documented this limitation and proposed a path forward involving literature review and mechanism modification."

### On Technical Skills

> "I refactored the codebase into a modular Python package with three core modules: DegradationModel (ODE system and diffusion), ParameterAnalysis (sensitivity sweeps), and ValidationMetrics (R², RMSE, MAE calculations). The code uses scipy's stiff ODE solvers (Radau/BDF) for numerical integration."

---

## References

1. Colin, X., et al. (2009). "Aging of Polyethylene Pipes Transporting Drinking Water Disinfected by Chlorine Dioxide. Part I: Chemical Aspects." *Polymer Engineering & Science*, 49(7), 1429-1437.

2. Colin, X., et al. (2009). "Aging of Polyethylene Pipes Transporting Drinking Water Disinfected by Chlorine Dioxide. Part II: Lifetime Prediction." *Polymer Engineering & Science*, 49(8), 1642-1652.

3. Khelidj, N., et al. (2006). "A general kinetic model for the photothermal oxidation of polypropylene." *Polymer Degradation and Stability*, 91(7), 1593-1597.

---

## Contact

This research was conducted as part of Stage 2025 at Laboratoire PIMM, Arts et Métiers ParisTech, in collaboration with SUEZ Environnement.

For questions about the model or analysis methodology, please refer to the code documentation in `src/`.

---

*Document generated: December 2024*
