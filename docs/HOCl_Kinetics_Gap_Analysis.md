# HOCl Kinetics Gap Analysis

## Executive Summary

This document presents a comprehensive analysis of the DOC→HOCl kinetics gap identified in the Colin et al. (2009) PE degradation model. The key finding is that **the discrepancy is fundamentally a mechanistic difference, not a concentration scaling issue**.

## 1. Background

### Problem Statement
The Colin et al. kinetic model was developed for chlorine dioxide (DOC/ClO₂) exposure, but SUEZ PE100 pipes are exposed to hypochlorous acid (HOCl). Literature suggests DOC is approximately 4× more aggressive than HOCl for stabilizer consumption (Mikdam et al. 2017).

### Initial Hypothesis
If DOC is 4× more aggressive, then:
```
DOC_effective = HOCl_nominal / 4 = 0.05 ppm / 4 = 0.0125 ppm
```

## 2. Analysis Results

### 2.1 DOC Concentration Sweep

| DOC (ppm) | Aggression Ratio | R² | RMSE (min) | Assessment |
|-----------|------------------|-----|------------|------------|
| 0.0050 | 10.0× | -0.30 | 99.3 | Way too slow |
| 0.0100 | 5.0× | 0.48 | 62.9 | Too slow |
| **0.0125** | **4.0×** | **0.64** | **52.0** | **Literature expected** |
| 0.0200 | 2.5× | 0.87 | 31.3 | Better |
| 0.0300 | 1.67× | 0.96 | 18.0 | Good |
| **0.0390** | **1.28×** | **0.97** | **14.8** | **Best empirical fit** |
| 0.0500 | 1.0× | 0.96 | 17.3 | Original |
| 0.0750 | 0.67× | 0.91 | 26.3 | Too fast |

**Key Finding:** The literature-expected value (DOC=0.0125 ppm) gives WORSE fit than the original (RMSE=52 vs 17 min).

### 2.2 Error Pattern Analysis

| Time Point | Original (DOC=0.05) | Literature (DOC=0.0125) | Optimal (DOC=0.039) |
|------------|---------------------|-------------------------|---------------------|
| Month 1 | -3.8 min | +88.6 min | +15.0 min |
| Month 4 | -20.4 min | +27.4 min | -17.4 min |
| Month 9 | -9.1 min | -7.0 min | -9.1 min |

**Critical Observation:** The error pattern REVERSES between early and late times:
- Low DOC (literature): Too slow early, better late
- High DOC (original): Good early, too fast late
- No DOC value correctly captures the degradation SHAPE

## 3. Root Cause Analysis

### 3.1 DOC (ClO₂) Mechanism
- **Direct H-abstraction:** ClO₂ + P-H → P• + HClO₂
- First-order kinetics: d[AH]/dt ∝ k₁ₐ[DOC][AH]
- Rapid initiation → fast exponential decay curve

### 3.2 HOCl Mechanism
- **Dissociation-mediated:** HOCl ⇌ H⁺ + OCl⁻ (pH dependent, pKa = 7.5)
- Multiple reactive species: HOCl, OCl⁻, possibly HO• from photolysis
- Different rate law: possibly fractional order or saturation kinetics
- Slower initiation but sustained propagation

### 3.3 Consequence
The Colin et al. rate constants (k1d, k7, k8d) embed DOC-specific reaction mechanisms. Simply scaling concentration preserves the wrong curve shape.

## 4. Scientific Conclusions

1. **The DOC→HOCl gap is NOT fixable by concentration scaling alone**

2. **Minimum achievable RMSE is ~15 minutes** at DOC=0.039 ppm
   - Even at optimal DOC, errors range from +15 to -17 minutes
   - Systematic over-prediction at late times (model exhausts AH too completely)

3. **Literature 4× ratio doesn't apply to degradation kinetics**
   - The 4× ratio from Colin et al./Mikdam refers to stabilizer consumption RATE
   - Degradation MECHANISM is fundamentally different

4. **For accurate HOCl modeling, would need:**
   - New rate constants calibrated specifically for HOCl
   - Or explicit HOCl dissociation step with new intermediates
   - This requires new experimental data (beyond scope of current work)

## 5. Engineering Recommendations

### For Design Screening (Conservative)
Use empirical calibration: **DOC = 0.039 ppm for HOCl = 0.05 ppm**
- Provides R² = 0.97 (acceptable correlation)
- RMSE = 14.8 min (reasonable for screening)
- Known to over-predict late-time degradation (conservative for safety)

### Uncertainty Statement
Model predictions beyond Month 3 should include ±20 min uncertainty band.
Late-time predictions (6-9 months) may over-predict degradation by up to 100%.

### For Detailed Design
Additional experimental validation required at:
- Multiple HOCl concentrations (0.01, 0.05, 0.1 ppm)
- Multiple temperatures (30°C, 40°C, 50°C)
- Different PE grades if relevant

## 6. Generated Files

| File | Description |
|------|-------------|
| `figures/hocl_kinetics_gap_analysis.png` | Four-panel analysis figure |
| `figures/hocl_kinetics_gap_analysis.pdf` | Vector format for publication |
| `src/ModifiedKinetics.py` | HOCl calibration module |
| `src/analyze_hocl_kinetics_gap.py` | Analysis script |

## 7. References

1. Colin, X., et al. (2009). "Aging of Polyethylene Pipes Transporting Drinking Water Disinfected by Chlorine Dioxide." *Polymer Engineering & Science*, 49(7-8).

2. Mikdam, A., et al. (2017). "A kinetic model for predicting the oxidative degradation of additive free polyethylene in bleach desinfected water." *Polymer Degradation and Stability*, 146, 78-94.

---
*Document generated: December 2024*
*Stage 2025 Research Project - Section 3.1*
