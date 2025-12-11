"""
Validation Metrics Module

This module provides:
- Statistical validation metrics (R², RMSE, MAE, MAPE)
- SUEZ experimental data constants
- Comparison utilities for model validation
- Summary table generation

Author: Stage 2025 Research Project
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


# =============================================================================
# SUEZ Experimental Data (PE100 Film, 400μm, 40°C)
# =============================================================================

@dataclass
class SUEZExperimentalData:
    """
    SUEZ experimental data for PE100 film degradation.

    Data from SUEZ aging platform:
    - Material: PE100 film (400 μm thickness)
    - Temperature: 40°C
    - Conditions: H2O (control) and HOCl (0.05 ppm)
    - Measurement: OIT at 200°C (minutes)
    """

    # Time points (months)
    times_h2o: np.ndarray = None
    times_hocl: np.ndarray = None

    # OIT values (minutes)
    oit_h2o: np.ndarray = None
    oit_hocl: np.ndarray = None

    # Initial OIT (reference)
    ti0_oit: float = 291.07  # From measured initial value

    # Carbonyl index (if available)
    carbonyl_hocl: np.ndarray = None

    def __post_init__(self):
        # H2O control experiment (physical aging only)
        # Raw values minus 20 min "time under azote" correction from original code.ipynb
        self.times_h2o = np.array([0, 1, 2, 3, 4, 9])
        self.oit_h2o = np.array([271.07, 267.95, 258.70, 253.00, 248.03, 229.30])

        # HOCl exposure experiment (0.05 ppm)
        # Raw values minus 20 min "time under azote" correction from original code.ipynb
        self.times_hocl = np.array([0, 1, 2, 3, 4, 6, 9])
        self.oit_hocl = np.array([271.07, 126.81, 61.44, 32.29, 23.82, 19.06, 9.07])

        # Initial OIT before any aging (measured value)
        self.ti0_oit = 291.07

    @staticmethod
    def get_instance():
        """Get a singleton instance of SUEZ data."""
        return SUEZExperimentalData()


# =============================================================================
# Validation Metrics Functions
# =============================================================================

def calculate_r2(observed: np.ndarray, predicted: np.ndarray) -> float:
    """
    Calculate coefficient of determination (R²).

    R² = 1 - SS_res / SS_tot

    Args:
        observed: Observed (experimental) values
        predicted: Predicted (simulated) values

    Returns:
        R² value (can be negative for poor fits)
    """
    valid = ~(np.isnan(observed) | np.isnan(predicted))
    if np.sum(valid) < 2:
        return np.nan

    obs = observed[valid]
    pred = predicted[valid]

    ss_res = np.sum((obs - pred)**2)
    ss_tot = np.sum((obs - np.mean(obs))**2)

    if ss_tot < 1e-12:
        return 0.0

    return 1.0 - ss_res / ss_tot


def calculate_rmse(observed: np.ndarray, predicted: np.ndarray) -> float:
    """
    Calculate Root Mean Square Error.

    RMSE = sqrt(mean((obs - pred)²))

    Args:
        observed: Observed values
        predicted: Predicted values

    Returns:
        RMSE value
    """
    valid = ~(np.isnan(observed) | np.isnan(predicted))
    if np.sum(valid) < 1:
        return np.nan

    obs = observed[valid]
    pred = predicted[valid]

    return np.sqrt(np.mean((obs - pred)**2))


def calculate_mae(observed: np.ndarray, predicted: np.ndarray) -> float:
    """
    Calculate Mean Absolute Error.

    MAE = mean(|obs - pred|)

    Args:
        observed: Observed values
        predicted: Predicted values

    Returns:
        MAE value
    """
    valid = ~(np.isnan(observed) | np.isnan(predicted))
    if np.sum(valid) < 1:
        return np.nan

    obs = observed[valid]
    pred = predicted[valid]

    return np.mean(np.abs(obs - pred))


def calculate_mape(observed: np.ndarray, predicted: np.ndarray) -> float:
    """
    Calculate Mean Absolute Percentage Error.

    MAPE = mean(|obs - pred| / |obs|) * 100%

    Args:
        observed: Observed values
        predicted: Predicted values

    Returns:
        MAPE value (%)
    """
    valid = ~(np.isnan(observed) | np.isnan(predicted)) & (np.abs(observed) > 1e-12)
    if np.sum(valid) < 1:
        return np.nan

    obs = observed[valid]
    pred = predicted[valid]

    return np.mean(np.abs((obs - pred) / obs)) * 100


def calculate_all_metrics(
    observed: np.ndarray,
    predicted: np.ndarray
) -> Dict[str, float]:
    """
    Calculate all validation metrics.

    Args:
        observed: Observed values
        predicted: Predicted values

    Returns:
        Dictionary with R², RMSE, MAE, MAPE
    """
    return {
        'R2': calculate_r2(observed, predicted),
        'RMSE': calculate_rmse(observed, predicted),
        'MAE': calculate_mae(observed, predicted),
        'MAPE': calculate_mape(observed, predicted)
    }


# =============================================================================
# Validation Comparison Class
# =============================================================================

class ValidationComparison:
    """
    Comprehensive validation comparison for model results.

    Compares simulation outputs against experimental data
    and generates summary tables/reports.
    """

    def __init__(self, suez_data: Optional[SUEZExperimentalData] = None):
        """
        Initialize validation comparison.

        Args:
            suez_data: SUEZ experimental data (uses default if None)
        """
        self.suez_data = suez_data or SUEZExperimentalData()
        self.comparison_results = {}

    def compare_h2o(
        self,
        sim_times_months: np.ndarray,
        sim_oit_avg: np.ndarray
    ) -> Dict[str, any]:
        """
        Compare simulation with H2O (control) experimental data.

        Args:
            sim_times_months: Simulation time array (months)
            sim_oit_avg: Simulated average OIT (minutes)

        Returns:
            Dictionary with metrics and comparison data
        """
        exp_times = self.suez_data.times_h2o
        exp_oit = self.suez_data.oit_h2o

        # Interpolate simulation to experimental times
        sim_oit_at_exp_times = np.interp(exp_times, sim_times_months, sim_oit_avg)

        # Calculate metrics
        metrics = calculate_all_metrics(exp_oit, sim_oit_at_exp_times)

        result = {
            'condition': 'H2O',
            'exp_times': exp_times,
            'exp_oit': exp_oit,
            'sim_oit': sim_oit_at_exp_times,
            'residuals': exp_oit - sim_oit_at_exp_times,
            'metrics': metrics
        }

        self.comparison_results['H2O'] = result
        return result

    def compare_hocl(
        self,
        sim_times_months: np.ndarray,
        sim_oit_avg: np.ndarray
    ) -> Dict[str, any]:
        """
        Compare simulation with HOCl experimental data.

        Args:
            sim_times_months: Simulation time array (months)
            sim_oit_avg: Simulated average OIT (minutes)

        Returns:
            Dictionary with metrics and comparison data
        """
        exp_times = self.suez_data.times_hocl
        exp_oit = self.suez_data.oit_hocl

        # Interpolate simulation to experimental times
        sim_oit_at_exp_times = np.interp(exp_times, sim_times_months, sim_oit_avg)

        # Calculate metrics
        metrics = calculate_all_metrics(exp_oit, sim_oit_at_exp_times)

        result = {
            'condition': 'HOCl',
            'exp_times': exp_times,
            'exp_oit': exp_oit,
            'sim_oit': sim_oit_at_exp_times,
            'residuals': exp_oit - sim_oit_at_exp_times,
            'metrics': metrics
        }

        self.comparison_results['HOCl'] = result
        return result

    def generate_comparison_table(self) -> pd.DataFrame:
        """
        Generate a summary comparison table.

        Returns:
            DataFrame with point-by-point comparison
        """
        rows = []

        for condition, result in self.comparison_results.items():
            for i in range(len(result['exp_times'])):
                rows.append({
                    'Condition': condition,
                    'Time (months)': result['exp_times'][i],
                    'OIT Exp (min)': result['exp_oit'][i],
                    'OIT Sim (min)': result['sim_oit'][i],
                    'Residual': result['residuals'][i],
                    'Error (%)': result['residuals'][i] / result['exp_oit'][i] * 100
                        if result['exp_oit'][i] > 1e-6 else np.nan
                })

        return pd.DataFrame(rows)

    def generate_metrics_summary(self) -> pd.DataFrame:
        """
        Generate summary table of validation metrics.

        Returns:
            DataFrame with metrics for each condition
        """
        rows = []

        for condition, result in self.comparison_results.items():
            m = result['metrics']
            rows.append({
                'Condition': condition,
                'R²': m['R2'],
                'RMSE (min)': m['RMSE'],
                'MAE (min)': m['MAE'],
                'MAPE (%)': m['MAPE'],
                'N points': len(result['exp_times'])
            })

        return pd.DataFrame(rows)

    def print_summary(self):
        """Print a formatted validation summary."""
        print("\n" + "="*70)
        print("VALIDATION SUMMARY")
        print("="*70)

        for condition, result in self.comparison_results.items():
            m = result['metrics']
            print(f"\n{condition}:")
            print(f"  R²    = {m['R2']:.4f}")
            print(f"  RMSE  = {m['RMSE']:.2f} min")
            print(f"  MAE   = {m['MAE']:.2f} min")
            print(f"  MAPE  = {m['MAPE']:.1f}%")

            # Quality assessment
            if m['R2'] >= 0.95:
                quality = "EXCELLENT"
            elif m['R2'] >= 0.90:
                quality = "GOOD"
            elif m['R2'] >= 0.75:
                quality = "ACCEPTABLE"
            elif m['R2'] >= 0.50:
                quality = "POOR"
            else:
                quality = "INADEQUATE"

            print(f"  Quality: {quality}")

        print("\n" + "="*70)


# =============================================================================
# Additional Analysis Functions
# =============================================================================

def analyze_model_limitations(
    h2o_comparison: Dict,
    hocl_comparison: Dict
) -> str:
    """
    Generate analysis text identifying model limitations.

    Args:
        h2o_comparison: H2O validation result
        hocl_comparison: HOCl validation result

    Returns:
        Analysis text string
    """
    h2o_r2 = h2o_comparison['metrics']['R2']
    hocl_r2 = hocl_comparison['metrics']['R2']

    analysis = []
    analysis.append("="*70)
    analysis.append("MODEL LIMITATIONS ANALYSIS")
    analysis.append("="*70)

    # H2O analysis
    analysis.append(f"\n1. H2O (Control) Validation: R² = {h2o_r2:.4f}")
    if h2o_r2 >= 0.90:
        analysis.append("   -> Physical aging (AH loss by diffusion/evaporation) is well captured")
        analysis.append("   -> Diffusion coefficients and boundary conditions are appropriate")
    else:
        analysis.append("   -> Physical aging mechanism needs review")

    # HOCl analysis
    analysis.append(f"\n2. HOCl Validation: R² = {hocl_r2:.4f}")
    if hocl_r2 < 0.70:
        analysis.append("   -> SIGNIFICANT DISCREPANCY IDENTIFIED")
        analysis.append("\n   ROOT CAUSE: DOC→HOCl kinetics approximation")
        analysis.append("   - Model uses DOC rate constants (k1d, k7, k8d)")
        analysis.append("   - HOCl has different oxidation mechanism:")
        analysis.append("     * DOC is a radical (•O-Cl=O) - direct H-abstraction")
        analysis.append("     * HOCl dissociates: HOCl ⇌ HO• + Cl•")
        analysis.append("     * HO• is more reactive than DOC")
        analysis.append("   - Observed: Model predicts TOO FAST degradation")
        analysis.append("   - This suggests effective DOC equivalent is LOWER than assumed")

    # Residual analysis
    analysis.append("\n3. Residual Analysis:")
    h2o_residuals = h2o_comparison['residuals']
    hocl_residuals = hocl_comparison['residuals']

    if np.all(h2o_residuals >= 0):
        analysis.append("   H2O: Model consistently UNDER-predicts OIT")
    elif np.all(h2o_residuals <= 0):
        analysis.append("   H2O: Model consistently OVER-predicts OIT")
    else:
        analysis.append("   H2O: Residuals change sign - possible temporal bias")

    if np.all(hocl_residuals >= 0):
        analysis.append("   HOCl: Model consistently UNDER-predicts OIT (predicts faster degradation)")
    elif np.all(hocl_residuals <= 0):
        analysis.append("   HOCl: Model consistently OVER-predicts OIT (predicts slower degradation)")
    else:
        # Check if residuals change sign over time
        early_residuals = hocl_residuals[:2] if len(hocl_residuals) > 2 else hocl_residuals
        late_residuals = hocl_residuals[-2:] if len(hocl_residuals) > 2 else hocl_residuals

        if np.mean(early_residuals) > 0 and np.mean(late_residuals) < 0:
            analysis.append("   HOCl: Early match THEN increasing divergence")
            analysis.append("        -> Suggests initiation rate is acceptable but propagation differs")

    # Recommendations
    analysis.append("\n4. RECOMMENDATIONS:")
    if hocl_r2 < 0.70:
        analysis.append("   a) Literature review on HOCl oxidation kinetics")
        analysis.append("   b) Determine effective DOC equivalent for HOCl")
        analysis.append("   c) Consider explicit HOCl dissociation step")
        analysis.append("   d) Validate at multiple HOCl concentrations")

    return "\n".join(analysis)


def quantify_model_gap(
    hocl_comparison: Dict,
    threshold_months: List[float] = [1, 2, 3, 4, 6]
) -> pd.DataFrame:
    """
    Quantify the gap between model and experiment at specific time points.

    Args:
        hocl_comparison: HOCl validation result
        threshold_months: Time points to analyze

    Returns:
        DataFrame with gap analysis
    """
    rows = []

    exp_times = hocl_comparison['exp_times']
    exp_oit = hocl_comparison['exp_oit']
    sim_oit = hocl_comparison['sim_oit']

    for t in threshold_months:
        idx = np.argmin(np.abs(exp_times - t))
        if abs(exp_times[idx] - t) < 0.5:  # Within 0.5 months
            exp_val = exp_oit[idx]
            sim_val = sim_oit[idx]
            gap = exp_val - sim_val
            gap_pct = gap / exp_val * 100 if exp_val > 1e-6 else np.nan

            # Correction factor needed
            correction = exp_val / sim_val if sim_val > 1e-6 else np.nan

            rows.append({
                'Time (months)': exp_times[idx],
                'OIT Exp (min)': exp_val,
                'OIT Sim (min)': sim_val,
                'Gap (min)': gap,
                'Gap (%)': gap_pct,
                'Correction Factor': correction
            })

    return pd.DataFrame(rows)
