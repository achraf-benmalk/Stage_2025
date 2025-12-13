# Multi-Physics Modeling of Polyethylene Pipe Degradation Under Chlorinated Water Exposure

---

## **Stage 2025 Research Project | Laboratoire PIMM, ENSAM | Industrial Partner: SUEZ**

---

<table>
<tr>
<td width="33%" valign="top">

## 1. Introduction & Objectives

### Problem Statement
Polyethylene (PE) pipes in drinking water networks undergo accelerated degradation when exposed to chlorine-based disinfectants. Predicting service lifetime requires understanding the coupled chemical-physical processes.

### Research Objectives
1. Implement and validate Colin et al. (2009) diffusion-reaction model
2. Quantify parameter sensitivity for engineering applications
3. Characterize HOCl vs DOC kinetics gap
4. Develop uncertainty quantification framework
5. Establish chemo-mechanical coupling (Stage 3)

### Model System
- **Material:** PE100 (high-density polyethylene)
- **Geometry:** 400 μm film / 4mm pipe wall
- **Exposure:** Chlorinated water (HOCl/ClO₂)
- **Temperature:** 25-60°C

---

## 2. Methodology

### 2.1 Kinetic Framework
12-species diffusion-reaction model:
- O₂, DOC/HOCl, AH (antioxidant)
- P•, PO₂•, POOH (radicals/hydroperoxides)
- PH, Q, CO, PCl, S, X (polymer species)

### 2.2 Key Reactions
```
Initiation:   PH + DOC → P• + HClO₂ (k₁ᵈ)
Propagation:  P• + O₂ → PO₂• (k₂)
              PO₂• + PH → POOH + P• (k₃)
Stabilization: PO₂• + AH → POOH + A• (k₇)
```

### 2.3 OIT Calculation
```
OIT(t) = ti₀ × [AH(t)/AH₀]
```
Where ti₀ = 291 min (initial OIT)

### 2.4 Analysis Methods
- Parameter sweeps (DOC, T, AH₀)
- Monte Carlo uncertainty (N=500)
- Arrhenius lifetime fitting
- Chemo-mechanical FEM coupling

</td>
<td width="34%" valign="top">

## 3. Key Results

### 3.1 Parameter Sensitivity

| Parameter | Sensitivity | Impact |
|-----------|-------------|--------|
| **DOC conc.** | 278 min/ppm | **DOMINANT** |
| Temperature | 3.0 min/°C | Secondary |
| AH₀ loading | Linear | Proportional |

**Finding:** Disinfectant concentration is the primary control lever for pipe lifetime.

---

### 3.2 HOCl Kinetics Gap

| DOC (ppm) | Ratio | RMSE | Assessment |
|-----------|-------|------|------------|
| 0.0125 | 4× | 52.0 | Literature: POOR |
| **0.039** | **1.28×** | **14.8** | **Best fit** |
| 0.050 | 1× | 17.3 | Original |

**Root Cause:** DOC→HOCl is a **mechanistic difference**, not simple concentration scaling.

---

### 3.3 Uncertainty Quantification

Monte Carlo analysis (N=500 samples):

| Month | Mean OIT | CV (%) | 90% CI |
|-------|----------|--------|--------|
| 1 | 130.8 min | 10% | [110, 153] |
| 3 | 13.2 min | 41% | [6, 23] |
| 6 | 0.3 min | 92% | [0, 0.7] |

---

### 3.4 Lifetime Correlation

**Engineering Equation:**
```
t_life = 4.1×10⁻⁵ × exp(25,700/RT) × DOC⁻⁰·⁴⁸
```

| Parameter | Value | Unit |
|-----------|-------|------|
| Eₐ (apparent) | 25.7 ± 1.1 | kJ/mol |
| DOC exponent | -0.48 ± 0.02 | — |
| R² | 0.984 | — |

---

### 3.5 Chemo-Mechanical Coupling

**Failure Criterion Comparison:**

| Criterion | Lifetime | Ratio |
|-----------|----------|-------|
| OIT < 10 min | 3.4 months | 1× |
| ε_b < 50% | 24.0 months | **7×** |

**Implication:** OIT criterion provides 7× conservative safety margin.

</td>
<td width="33%" valign="top">

## 4. Validation Results

### SUEZ Experimental Comparison
**Conditions:** PE100 film, 40°C, 9 months

| Condition | R² | RMSE (min) |
|-----------|-----|------------|
| H₂O Control | -3.84 | 30.5 |
| HOCl (orig) | 0.96 | 17.2 |
| **HOCl (cal)** | **0.97** | **14.8** |

---

## 5. Stage 3: FEM Integration

### Coupled Analysis Framework
1. **Chemistry:** 1D diffusion-reaction (OIT profile)
2. **Material:** Mw evolution → property degradation
3. **Structural:** 2D FEM stress analysis

### Material Property Evolution
- Molecular weight: Mw(t) from Saito equation
- Yield stress: σ_y = f(Mw, crystallinity)
- Elongation: ductile-brittle transition

### Structural Analysis
- Internal pressure: 6 bar (standard)
- Hoop stress distribution
- Safety factor evolution

---

## 6. Conclusions

### Scientific Contributions
1. First systematic sensitivity analysis of Colin model
2. DOC→HOCl requires new kinetics, not scaling
3. Monte Carlo framework for PE degradation
4. Chemo-mechanical coupling for lifetime prediction

### Engineering Recommendations
1. Use calibrated DOC = 0.039 ppm for HOCl = 0.05 ppm
2. OIT criterion provides 7× conservative margin
3. Include ±20 min uncertainty beyond Month 3

### Future Work
- HOCl kinetics calibration
- Multi-temperature validation
- Field data comparison

---

## 7. References

1. Colin, X., et al. (2009). *Polym. Eng. Sci.*, 49(7-8)
2. Fayolle, B., et al. (2007). *Polym. Degrad. Stab.*, 92
3. Mikdam, A., et al. (2017). *Polym. Degrad. Stab.*, 146

---

**Acknowledgments:** Prof. Xavier Colin (PIMM), SUEZ Environnement

</td>
</tr>
</table>

---

## Visual Summary

### Figure References
| Figure | Description | Location |
|--------|-------------|----------|
| Sensitivity Analysis | Parameter importance tornado diagram | `figures/tornado_diagram.png` |
| HOCl Gap Analysis | DOC calibration comparison | `figures/hocl_kinetics_gap_analysis.png` |
| Uncertainty Bands | Monte Carlo 90% confidence intervals | `figures/mc_uncertainty_bands.png` |
| Lifetime Correlation | Arrhenius fitting validation | `figures/lifetime_correlation.png` |
| Chemo-Mechanical | Stage 3 coupling results | `figures/chemo_mechanical_coupling.png` |
| FEM Analysis | Structural stress distribution | `figures/structural_stress_comparison.png` |
| Validation | SUEZ experimental comparison | `figures/validation_summary.png` |

---

*Stage 2025 Research Project | December 2024*
