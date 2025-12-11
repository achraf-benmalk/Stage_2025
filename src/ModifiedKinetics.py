"""
Modified Kinetics Module for HOCl Exposure

This module addresses the DOC→HOCl kinetics gap identified in model validation.

Literature Findings:
-------------------
1. ClO2 (DOC) is MORE aggressive than HOCl (hypochlorite):
   - ClO2 stabilizer consumption rate is 4× greater than chlorinated water
   - ClO2 diffuses faster as a dissolved gas
   - Source: Colin et al., Mikdam et al. (2017)

2. Therefore, to model HOCl exposure using DOC kinetics:
   - The effective DOC concentration should be LOWER than nominal HOCl
   - Conversion factor: DOC_eff = HOCl_nominal / aggression_ratio

3. Optimization approach:
   - Back-calculate effective DOC that minimizes RMSE vs SUEZ data
   - Verify the ratio is physically reasonable (~4× based on literature)

Author: Stage 2025 Research Project
"""

import numpy as np
from scipy.optimize import minimize_scalar, minimize
from typing import Tuple, Dict, Optional
import sys
import os

# Add parent to path for imports
sys.path.insert(0, os.path.dirname(__file__))

from DegradationModel import DegradationModel, SimulationParams, create_suez_film_model
from ValidationMetrics import SUEZExperimentalData, calculate_r2, calculate_rmse, calculate_mae


class HOClKineticsCalibrator:
    """
    Calibrate effective DOC concentration for HOCl exposure.

    The base model uses DOC kinetics. To model HOCl, we need to find
    an effective DOC concentration that matches experimental data.
    """

    def __init__(self, model: Optional[DegradationModel] = None):
        """
        Initialize calibrator.

        Args:
            model: DegradationModel instance (creates SUEZ film model if None)
        """
        self.model = model or create_suez_film_model()
        self.suez_data = SUEZExperimentalData()
        self.calibration_results = {}

    def _run_simulation(self, doc_ppm: float, verbose: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """Run simulation and return (times_months, oit_avg)."""
        params = SimulationParams(
            T_celsius=40.0,
            DOC_ppm=doc_ppm,
            t_end_years=0.75,
            n_timepoints=100
        )
        result = self.model.simulate(params, verbose=verbose)

        if not result['success']:
            return None, None

        times_m, oit_avg = self.model.calculate_average_oit(result)
        return times_m, oit_avg

    def _objective_rmse(self, doc_ppm: float) -> float:
        """Objective function: RMSE vs SUEZ HOCl data."""
        times_m, oit_avg = self._run_simulation(doc_ppm)

        if times_m is None:
            return 1e10

        # Interpolate at experimental times
        sim_oit = np.interp(self.suez_data.times_hocl, times_m, oit_avg)

        # Calculate RMSE
        rmse = np.sqrt(np.mean((self.suez_data.oit_hocl - sim_oit)**2))
        return rmse

    def _objective_weighted(self, doc_ppm: float, late_weight: float = 2.0) -> float:
        """
        Weighted objective: emphasize late-time fit.

        The model diverges at late times, so we weight those points more.
        """
        times_m, oit_avg = self._run_simulation(doc_ppm)

        if times_m is None:
            return 1e10

        # Interpolate at experimental times
        sim_oit = np.interp(self.suez_data.times_hocl, times_m, oit_avg)

        # Weight: early (0-2 months) = 1, late (>2 months) = late_weight
        weights = np.where(self.suez_data.times_hocl <= 2, 1.0, late_weight)

        # Weighted RMSE
        errors_sq = (self.suez_data.oit_hocl - sim_oit)**2
        wmse = np.sum(weights * errors_sq) / np.sum(weights)

        return np.sqrt(wmse)

    def calibrate_effective_doc(
        self,
        doc_range: Tuple[float, float] = (0.001, 0.1),
        method: str = 'rmse',
        verbose: bool = True
    ) -> Dict:
        """
        Find optimal effective DOC concentration.

        Args:
            doc_range: (min, max) DOC values to search
            method: 'rmse' or 'weighted'
            verbose: Print progress

        Returns:
            Dictionary with calibration results
        """
        if verbose:
            print("="*70)
            print("HOCl KINETICS CALIBRATION")
            print("="*70)
            print(f"Searching DOC range: {doc_range[0]:.4f} - {doc_range[1]:.4f} ppm")
            print(f"Objective: {method}")
            print()

        # Select objective function
        if method == 'rmse':
            objective = self._objective_rmse
        else:
            objective = self._objective_weighted

        # Grid search first to find approximate minimum
        doc_values = np.linspace(doc_range[0], doc_range[1], 20)
        obj_values = []

        if verbose:
            print("Grid search...")

        for doc in doc_values:
            obj = objective(doc)
            obj_values.append(obj)
            if verbose:
                print(f"  DOC={doc:.4f} ppm: objective={obj:.2f}")

        # Find best from grid
        best_idx = np.argmin(obj_values)
        doc_init = doc_values[best_idx]

        if verbose:
            print(f"\nGrid search best: DOC={doc_init:.4f} ppm, obj={obj_values[best_idx]:.2f}")

        # Refine with bounded optimization
        if verbose:
            print("\nRefining with bounded optimization...")

        result = minimize_scalar(
            objective,
            bounds=doc_range,
            method='bounded',
            options={'xatol': 1e-5}
        )

        doc_optimal = result.x
        obj_optimal = result.fun

        # Run final simulation with optimal DOC
        times_m, oit_avg = self._run_simulation(doc_optimal, verbose=False)
        sim_oit = np.interp(self.suez_data.times_hocl, times_m, oit_avg)

        # Calculate metrics
        r2 = calculate_r2(self.suez_data.oit_hocl, sim_oit)
        rmse = calculate_rmse(self.suez_data.oit_hocl, sim_oit)
        mae = calculate_mae(self.suez_data.oit_hocl, sim_oit)

        # Calculate aggression ratio (DOC effective / HOCl nominal)
        hocl_nominal = 0.05  # ppm (SUEZ experimental condition)
        aggression_ratio = hocl_nominal / doc_optimal

        results = {
            'doc_optimal': doc_optimal,
            'hocl_nominal': hocl_nominal,
            'aggression_ratio': aggression_ratio,
            'r2': r2,
            'rmse': rmse,
            'mae': mae,
            'sim_times': times_m,
            'sim_oit': oit_avg,
            'sim_oit_at_exp': sim_oit,
            'exp_times': self.suez_data.times_hocl,
            'exp_oit': self.suez_data.oit_hocl
        }

        self.calibration_results = results

        if verbose:
            print(f"\n{'='*70}")
            print("CALIBRATION RESULTS")
            print("="*70)
            print(f"Optimal effective DOC: {doc_optimal:.5f} ppm")
            print(f"HOCl nominal: {hocl_nominal:.2f} ppm")
            print(f"Aggression ratio (HOCl/DOC_eff): {aggression_ratio:.2f}×")
            print(f"  (Literature suggests ~4× based on stabilizer consumption rates)")
            print()
            print(f"Validation metrics:")
            print(f"  R² = {r2:.4f}")
            print(f"  RMSE = {rmse:.2f} min")
            print(f"  MAE = {mae:.2f} min")
            print()
            print("Point-by-point comparison:")
            print(f"{'Month':<8} {'Exp OIT':<12} {'Sim OIT':<12} {'Error':<10} {'Error %':<10}")
            print("-"*52)
            for i, t in enumerate(self.suez_data.times_hocl):
                exp = self.suez_data.oit_hocl[i]
                sim = sim_oit[i]
                err = sim - exp
                err_pct = err / exp * 100 if exp > 0.1 else np.nan
                print(f"{t:<8.0f} {exp:<12.2f} {sim:<12.2f} {err:<10.2f} {err_pct:>9.1f}%")

        return results

    def compare_original_vs_calibrated(self, verbose: bool = True) -> Dict:
        """
        Compare original DOC=0.05 vs calibrated effective DOC.

        Returns:
            Dictionary with comparison results
        """
        if not self.calibration_results:
            raise ValueError("Must run calibrate_effective_doc first")

        # Original simulation (DOC=0.05 ppm as-is)
        times_orig, oit_orig = self._run_simulation(0.05)
        sim_orig = np.interp(self.suez_data.times_hocl, times_orig, oit_orig)
        r2_orig = calculate_r2(self.suez_data.oit_hocl, sim_orig)
        rmse_orig = calculate_rmse(self.suez_data.oit_hocl, sim_orig)

        # Calibrated results
        r2_calib = self.calibration_results['r2']
        rmse_calib = self.calibration_results['rmse']
        sim_calib = self.calibration_results['sim_oit_at_exp']

        comparison = {
            'original': {
                'doc_ppm': 0.05,
                'r2': r2_orig,
                'rmse': rmse_orig,
                'sim_oit': sim_orig
            },
            'calibrated': {
                'doc_ppm': self.calibration_results['doc_optimal'],
                'r2': r2_calib,
                'rmse': rmse_calib,
                'sim_oit': sim_calib
            },
            'improvement': {
                'r2_delta': r2_calib - r2_orig,
                'rmse_delta': rmse_orig - rmse_calib,  # Positive = improvement
                'rmse_pct_improvement': (rmse_orig - rmse_calib) / rmse_orig * 100
            }
        }

        if verbose:
            print("\n" + "="*70)
            print("COMPARISON: ORIGINAL vs CALIBRATED")
            print("="*70)
            print(f"\n{'Metric':<20} {'Original':<20} {'Calibrated':<20} {'Improvement':<15}")
            print("-"*75)
            print(f"{'DOC (ppm)':<20} {0.05:<20.4f} {self.calibration_results['doc_optimal']:<20.5f}")
            print(f"{'R²':<20} {r2_orig:<20.4f} {r2_calib:<20.4f} {r2_calib - r2_orig:>+14.4f}")
            print(f"{'RMSE (min)':<20} {rmse_orig:<20.2f} {rmse_calib:<20.2f} {rmse_orig - rmse_calib:>+14.2f}")

            print("\n" + "-"*75)
            print("Point-by-point errors:")
            print(f"{'Month':<8} {'Exp':<10} {'Orig':<10} {'Calib':<10} {'Orig Err':<12} {'Calib Err':<12}")
            print("-"*62)
            for i, t in enumerate(self.suez_data.times_hocl):
                exp = self.suez_data.oit_hocl[i]
                orig = sim_orig[i]
                calib = sim_calib[i]
                print(f"{t:<8.0f} {exp:<10.1f} {orig:<10.1f} {calib:<10.1f} {orig-exp:>+11.1f} {calib-exp:>+11.1f}")

        return comparison


def create_hocl_conversion_factor() -> float:
    """
    Calculate the HOCl→DOC effective conversion factor.

    Based on literature:
    - ClO2 is ~4× more aggressive than HOCl for stabilizer consumption
    - Therefore: DOC_effective = HOCl_nominal / 4

    Returns:
        Conversion factor (multiply HOCl ppm by this to get effective DOC)
    """
    # Literature value: DOC is 4× more aggressive
    aggression_ratio = 4.0
    return 1.0 / aggression_ratio


class ModifiedDegradationModel:
    """
    Wrapper around DegradationModel with HOCl kinetics correction.

    Automatically converts HOCl concentration to effective DOC.
    """

    def __init__(
        self,
        base_model: Optional[DegradationModel] = None,
        conversion_factor: Optional[float] = None
    ):
        """
        Initialize modified model.

        Args:
            base_model: Base DegradationModel (creates SUEZ model if None)
            conversion_factor: HOCl→DOC conversion (uses calibrated if None)
        """
        self.base_model = base_model or create_suez_film_model()

        if conversion_factor is None:
            # Use literature-based default
            self.conversion_factor = create_hocl_conversion_factor()
        else:
            self.conversion_factor = conversion_factor

    def simulate_hocl(
        self,
        T_celsius: float,
        HOCl_ppm: float,
        t_end_years: float,
        **kwargs
    ) -> dict:
        """
        Run simulation with HOCl exposure (automatically converts to effective DOC).

        Args:
            T_celsius: Temperature
            HOCl_ppm: HOCl concentration in water
            t_end_years: Simulation duration
            **kwargs: Additional SimulationParams arguments

        Returns:
            Simulation result dictionary
        """
        # Convert HOCl to effective DOC
        doc_effective = HOCl_ppm * self.conversion_factor

        params = SimulationParams(
            T_celsius=T_celsius,
            DOC_ppm=doc_effective,
            t_end_years=t_end_years,
            **kwargs
        )

        result = self.base_model.simulate(params, verbose=kwargs.get('verbose', True))

        # Add HOCl info to result
        result['HOCl_ppm'] = HOCl_ppm
        result['DOC_effective'] = doc_effective
        result['conversion_factor'] = self.conversion_factor

        return result

    def calculate_average_oit(self, result: dict) -> Tuple[np.ndarray, np.ndarray]:
        """Convenience method to get OIT from result."""
        return self.base_model.calculate_average_oit(result)


# =============================================================================
# Main execution for testing
# =============================================================================

if __name__ == '__main__':
    print("HOCl Kinetics Calibration\n")

    # Create calibrator
    calibrator = HOClKineticsCalibrator()

    # Run calibration
    results = calibrator.calibrate_effective_doc(
        doc_range=(0.005, 0.05),
        method='rmse',
        verbose=True
    )

    # Compare original vs calibrated
    comparison = calibrator.compare_original_vs_calibrated(verbose=True)

    print("\n" + "="*70)
    print("CONCLUSION")
    print("="*70)
    print(f"\nThe optimal effective DOC is {results['doc_optimal']:.5f} ppm")
    print(f"This corresponds to an aggression ratio of {results['aggression_ratio']:.2f}×")
    print(f"Literature suggests ~4× (DOC is 4× more aggressive than HOCl)")

    if abs(results['aggression_ratio'] - 4.0) < 2.0:
        print("\n✓ Calibrated ratio is physically reasonable!")
    else:
        print(f"\n⚠ Calibrated ratio differs from literature expectation")
