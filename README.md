# PE Pipe Degradation Model

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**1D Diffusion-Reaction Kinetic Model for Polyethylene Pipe Degradation Under Chlorine Exposure**

This repository contains the implementation and analysis of a mechanistic model for predicting PE pipe degradation, based on the Colin et al. (2009) kinetic framework.

## Project Structure

```
Stage_2025/
├── src/                          # Core Python modules
│   ├── __init__.py
│   ├── DegradationModel.py       # Main diffusion-reaction model
│   ├── ParameterAnalysis.py      # Sensitivity analysis framework
│   └── ValidationMetrics.py      # Statistical validation tools
├── figures/                      # Generated figures (PNG/PDF)
├── data/                         # Analysis results (CSV)
├── notebooks/                    # Jupyter notebooks (if any)
├── run_sensitivity_analysis.py   # Main analysis script
├── RESEARCH_SUMMARY.md           # Key findings for CV/interviews
└── README.md                     # This file
```

## Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/achraf-benmalk/Stage_2025.git
cd Stage_2025

# Install dependencies
pip install numpy scipy pandas matplotlib
```

### Run Analysis

```bash
# Run complete sensitivity analysis
python run_sensitivity_analysis.py
```

This generates:
- Parameter sweep figures in `figures/`
- CSV data files in `data/`
- Console output with validation metrics and key findings

### Basic Usage

```python
from src import create_suez_film_model, SimulationParams

# Create model for SUEZ PE100 film (400 μm)
model = create_suez_film_model()

# Run simulation
params = SimulationParams(
    T_celsius=40.0,
    DOC_ppm=0.05,
    t_end_years=0.75
)
result = model.simulate(params)

# Get OIT profile
times, oit = model.calculate_average_oit(result)
print(f"Final OIT: {oit[-1]:.1f} min")
```

## Model Description

### Physical System

The model describes degradation of polyethylene pipes exposed to chlorinated drinking water. Key processes include:

1. **Diffusion** of O2, DOC (chlorine dioxide), and antioxidant (AH) through the polymer
2. **Chemical reactions** including:
   - DOC-initiated radical chain oxidation
   - Hydroperoxide decomposition (uni- and bimolecular)
   - Antioxidant consumption
   - Chain scissions and crosslinking

### Chemical Species (12)

| Species | Description | Mobile |
|---------|-------------|--------|
| O2 | Dissolved oxygen | Yes |
| DOC | Chlorine dioxide | Yes |
| AH | Phenolic antioxidant | Yes |
| P• | Alkyl radicals | No |
| PO2• | Peroxy radicals | No |
| POOH | Hydroperoxides | No |
| PH | Polymer (methylenes) | No |
| Q | Caged radical pair | No |
| CO | Carbonyls | No |
| PCl | Grafted chlorine | No |
| S | Chain scissions | No |
| X | Crosslinks | No |

### Key Equations

**OIT Calculation:**
```
OIT(t) = ti0 × [AH](t) / [AH]0
```

**Molecular Weight Evolution (Saito):**
```
1/Mw = 1/Mw0 + S/2 - 2X
```

**Embrittlement Criterion:**
```
Mw < 70 kg/mol → Material failure
```

## Analysis Results

### Parameter Sensitivity (at 6 months, 40°C)

| Parameter | Range | Impact |
|-----------|-------|--------|
| DOC concentration | 0-5 ppm | 278 min OIT range |
| Temperature | 30-60°C | 3 min OIT range |
| Initial [AH]₀ | 0.5-3× | 0.2 min OIT range |

**Finding:** DOC concentration dominates degradation kinetics.

### Temperature Acceleration

| Temperature | Factor vs 30°C |
|-------------|----------------|
| 30°C | 1.0× |
| 40°C | 1.4× |
| 50°C | 1.9× |
| 60°C | 2.4× |

### Validation

| Condition | R² | Quality |
|-----------|-----|---------|
| HOCl (0.05 ppm) | 0.97 | Excellent |
| H2O (control) | -2.29 | Under-predicts decay |

## Limitations

1. **DOC-HOCl Kinetics:** Model uses DOC rate constants as proxy for HOCl. Different oxidation mechanisms may affect accuracy.

2. **Physical Aging:** Pure water degradation under-predicted. Boundary extraction coefficient may need adjustment.

3. **Temperature Range:** Validated only at 40°C. Extrapolation to other temperatures requires caution.

## API Reference

### DegradationModel

```python
class DegradationModel:
    def __init__(self, L, nz, simulation_mode, material):
        """
        Initialize model.

        Args:
            L: Thickness (m)
            nz: Grid points
            simulation_mode: 'pipe' or 'film'
            material: MaterialProperties object
        """

    def simulate(self, params) -> dict:
        """Run simulation with given parameters."""

    def calculate_oit(self, AH_profile, AH0_ref) -> np.ndarray:
        """Calculate OIT from antioxidant profile."""
```

### SimulationParams

```python
@dataclass
class SimulationParams:
    T_celsius: float      # Temperature (°C)
    DOC_ppm: float        # DOC concentration (ppm)
    t_end_years: float    # Duration (years)
    AH0_mult: float = 1.0 # Antioxidant multiplier
```

### ParameterAnalysis

```python
class ParameterAnalysis:
    def sweep_doc_concentration(self, doc_values, ...) -> SweepResult
    def sweep_temperature(self, T_values, ...) -> SweepResult
    def sweep_antioxidant(self, AH_multipliers, ...) -> SweepResult
    def arrhenius_analysis(self, T_sweep_result, ...) -> dict
```

## References

1. Colin, X., et al. (2009). *Polymer Engineering & Science*, 49(7), 1429-1437.
2. Colin, X., et al. (2009). *Polymer Engineering & Science*, 49(8), 1642-1652.

## License

MIT License - see LICENSE file for details.

## Acknowledgments

- **Supervisor:** Prof. Xavier Colin, PIMM Laboratory
- **Industrial Partner:** SUEZ Environnement
- **Institution:** Arts et Métiers ParisTech

---

*Stage 2025 Research Project*
