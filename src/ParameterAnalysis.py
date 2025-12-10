"""
Parameter Sensitivity Analysis Framework

This module provides tools for systematic parameter analysis of the PE degradation model:
- Parameter sweeps (DOC concentration, temperature, antioxidant content)
- Sensitivity quantification
- Arrhenius analysis for activation energy extraction
- Results visualization and export

Author: Stage 2025 Research Project
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Callable
from dataclasses import dataclass
import time
import warnings

from DegradationModel import DegradationModel, SimulationParams, create_suez_film_model


@dataclass
class SweepResult:
    """Results from a parameter sweep."""
    parameter_name: str
    parameter_values: np.ndarray
    times_months: np.ndarray
    oit_avg_matrix: np.ndarray      # (n_params, n_times)
    oit_surface_matrix: np.ndarray  # (n_params, n_times)
    co_surface_matrix: np.ndarray   # (n_params, n_times)
    success_flags: List[bool]
    metadata: dict


class ParameterAnalysis:
    """
    Framework for systematic parameter sensitivity analysis.

    Provides methods for:
    - Single parameter sweeps
    - Multi-parameter grid searches
    - Arrhenius analysis
    - Sensitivity metrics calculation
    """

    def __init__(self, model: DegradationModel):
        """
        Initialize analysis framework.

        Args:
            model: Configured DegradationModel instance
        """
        self.model = model
        self.sweep_results: Dict[str, SweepResult] = {}

    def sweep_doc_concentration(
        self,
        doc_values: np.ndarray,
        T_celsius: float = 40.0,
        t_end_years: float = 0.75,
        n_timepoints: int = 50,
        AH0_mult: float = 1.0,
        verbose: bool = True
    ) -> SweepResult:
        """
        Sweep DOC/HOCl concentration values.

        Args:
            doc_values: Array of DOC concentrations (ppm)
            T_celsius: Temperature (°C)
            t_end_years: Simulation duration (years)
            n_timepoints: Number of output time points
            AH0_mult: Antioxidant multiplier
            verbose: Print progress

        Returns:
            SweepResult with all simulation outputs
        """
        return self._run_parameter_sweep(
            parameter_name='DOC_ppm',
            parameter_values=doc_values,
            base_params={
                'T_celsius': T_celsius,
                't_end_years': t_end_years,
                'AH0_mult': AH0_mult,
                'n_timepoints': n_timepoints
            },
            verbose=verbose
        )

    def sweep_temperature(
        self,
        T_values: np.ndarray,
        DOC_ppm: float = 0.05,
        t_end_years: float = 0.75,
        n_timepoints: int = 50,
        AH0_mult: float = 1.0,
        verbose: bool = True
    ) -> SweepResult:
        """
        Sweep temperature values.

        Args:
            T_values: Array of temperatures (°C)
            DOC_ppm: DOC concentration
            t_end_years: Simulation duration
            n_timepoints: Number of output points
            AH0_mult: Antioxidant multiplier
            verbose: Print progress

        Returns:
            SweepResult with all simulation outputs
        """
        return self._run_parameter_sweep(
            parameter_name='T_celsius',
            parameter_values=T_values,
            base_params={
                'DOC_ppm': DOC_ppm,
                't_end_years': t_end_years,
                'AH0_mult': AH0_mult,
                'n_timepoints': n_timepoints
            },
            verbose=verbose
        )

    def sweep_antioxidant(
        self,
        AH_multipliers: np.ndarray,
        T_celsius: float = 40.0,
        DOC_ppm: float = 0.05,
        t_end_years: float = 0.75,
        n_timepoints: int = 50,
        verbose: bool = True
    ) -> SweepResult:
        """
        Sweep initial antioxidant content.

        Args:
            AH_multipliers: Array of AH0 multipliers (1.0 = baseline)
            T_celsius: Temperature
            DOC_ppm: DOC concentration
            t_end_years: Simulation duration
            n_timepoints: Number of output points
            verbose: Print progress

        Returns:
            SweepResult with all simulation outputs
        """
        return self._run_parameter_sweep(
            parameter_name='AH0_mult',
            parameter_values=AH_multipliers,
            base_params={
                'T_celsius': T_celsius,
                'DOC_ppm': DOC_ppm,
                't_end_years': t_end_years,
                'n_timepoints': n_timepoints
            },
            verbose=verbose
        )

    def _run_parameter_sweep(
        self,
        parameter_name: str,
        parameter_values: np.ndarray,
        base_params: dict,
        verbose: bool = True
    ) -> SweepResult:
        """
        Internal method to run a parameter sweep.

        Args:
            parameter_name: Name of parameter being swept
            parameter_values: Array of values to test
            base_params: Dictionary of fixed parameters
            verbose: Print progress

        Returns:
            SweepResult object
        """
        n_values = len(parameter_values)

        if verbose:
            print(f"\n{'='*60}")
            print(f"Parameter Sweep: {parameter_name}")
            print(f"  Values: {parameter_values}")
            print(f"  Base conditions: {base_params}")
            print(f"{'='*60}")

        # Storage for results
        results_list = []
        success_flags = []
        times_months = None

        start_total = time.time()

        for i, val in enumerate(parameter_values):
            # Build params for this run
            run_params = base_params.copy()
            run_params[parameter_name] = val

            params = SimulationParams(**run_params)

            if verbose:
                print(f"\n[{i+1}/{n_values}] {parameter_name} = {val}")

            result = self.model.simulate(params, verbose=verbose)
            results_list.append(result)
            success_flags.append(result['success'])

            if times_months is None and result['success']:
                times_months = result['t_months']

        total_time = time.time() - start_total
        if verbose:
            print(f"\n{'='*60}")
            print(f"Sweep completed in {total_time:.1f}s")
            print(f"  Success rate: {sum(success_flags)}/{n_values}")
            print(f"{'='*60}")

        # Process results into matrices
        n_times = len(times_months) if times_months is not None else 0
        oit_avg_matrix = np.full((n_values, n_times), np.nan)
        oit_surface_matrix = np.full((n_values, n_times), np.nan)
        co_surface_matrix = np.full((n_values, n_times), np.nan)

        for i, result in enumerate(results_list):
            if result['success']:
                _, oit_avg = self.model.calculate_average_oit(result)
                _, oit_profiles = self.model.calculate_oit_profile_evolution(result)
                _, co_surface = self.model.get_surface_concentration(result, 'CO')

                oit_avg_matrix[i, :] = oit_avg
                oit_surface_matrix[i, :] = oit_profiles[0, :]
                co_surface_matrix[i, :] = co_surface

        sweep_result = SweepResult(
            parameter_name=parameter_name,
            parameter_values=parameter_values,
            times_months=times_months,
            oit_avg_matrix=oit_avg_matrix,
            oit_surface_matrix=oit_surface_matrix,
            co_surface_matrix=co_surface_matrix,
            success_flags=success_flags,
            metadata={
                'base_params': base_params,
                'model_config': {
                    'L_mm': self.model.L * 1000,
                    'mode': self.model.simulation_mode,
                    'ti0_oit': self.model.material.ti0_oit
                },
                'computation_time_s': total_time
            }
        )

        # Store result
        self.sweep_results[parameter_name] = sweep_result

        return sweep_result

    # =========================================================================
    # Sensitivity Analysis Methods
    # =========================================================================

    def calculate_sensitivity_metrics(
        self,
        sweep_result: SweepResult,
        reference_idx: Optional[int] = None,
        target_oit: float = 50.0
    ) -> pd.DataFrame:
        """
        Calculate sensitivity metrics from a parameter sweep.

        Metrics computed:
        - Time to reach target OIT for each parameter value
        - Relative change vs reference
        - OIT at fixed time points

        Args:
            sweep_result: SweepResult from a sweep
            reference_idx: Index of reference case (default: middle)
            target_oit: Target OIT value for time calculation

        Returns:
            DataFrame with sensitivity metrics
        """
        if reference_idx is None:
            reference_idx = len(sweep_result.parameter_values) // 2

        results = []

        for i, param_val in enumerate(sweep_result.parameter_values):
            oit_profile = sweep_result.oit_avg_matrix[i, :]
            times = sweep_result.times_months

            # Time to target OIT
            try:
                t_target = self._find_time_to_threshold(times, oit_profile, target_oit)
            except:
                t_target = np.nan

            # OIT at fixed times (1, 3, 6, 9 months)
            oit_1m = np.interp(1.0, times, oit_profile) if times[-1] >= 1.0 else np.nan
            oit_3m = np.interp(3.0, times, oit_profile) if times[-1] >= 3.0 else np.nan
            oit_6m = np.interp(6.0, times, oit_profile) if times[-1] >= 6.0 else np.nan
            oit_9m = np.interp(9.0, times, oit_profile) if times[-1] >= 9.0 else np.nan

            results.append({
                sweep_result.parameter_name: param_val,
                'time_to_oit_50': t_target,
                'oit_1month': oit_1m,
                'oit_3month': oit_3m,
                'oit_6month': oit_6m,
                'oit_9month': oit_9m,
                'success': sweep_result.success_flags[i]
            })

        df = pd.DataFrame(results)

        # Calculate relative changes vs reference
        ref_row = df.iloc[reference_idx]
        for col in ['time_to_oit_50', 'oit_1month', 'oit_3month', 'oit_6month', 'oit_9month']:
            if ref_row[col] > 0:
                df[f'{col}_rel'] = (df[col] - ref_row[col]) / ref_row[col] * 100

        return df

    def calculate_acceleration_factor(
        self,
        sweep_result: SweepResult,
        reference_idx: int = 0,
        target_oit: float = 50.0
    ) -> pd.DataFrame:
        """
        Calculate acceleration factors between parameter values.

        Acceleration factor = t_ref / t_test

        Args:
            sweep_result: SweepResult from a sweep
            reference_idx: Index of reference (baseline) case
            target_oit: Target OIT for time calculation

        Returns:
            DataFrame with acceleration factors
        """
        times = sweep_result.times_months
        results = []

        # Reference time to target
        oit_ref = sweep_result.oit_avg_matrix[reference_idx, :]
        t_ref = self._find_time_to_threshold(times, oit_ref, target_oit)

        for i, param_val in enumerate(sweep_result.parameter_values):
            oit_profile = sweep_result.oit_avg_matrix[i, :]

            try:
                t_val = self._find_time_to_threshold(times, oit_profile, target_oit)
                acc_factor = t_ref / t_val if t_val > 0 else np.inf
            except:
                t_val = np.nan
                acc_factor = np.nan

            results.append({
                sweep_result.parameter_name: param_val,
                'time_to_target': t_val,
                'acceleration_factor': acc_factor
            })

        return pd.DataFrame(results)

    def _find_time_to_threshold(
        self,
        times: np.ndarray,
        values: np.ndarray,
        threshold: float
    ) -> float:
        """Find time when values first cross below threshold."""
        # Find crossing point (values decreasing)
        below = values < threshold
        if not np.any(below):
            return np.inf
        first_below_idx = np.argmax(below)

        if first_below_idx == 0:
            return times[0]

        # Linear interpolation
        t1, t2 = times[first_below_idx - 1], times[first_below_idx]
        v1, v2 = values[first_below_idx - 1], values[first_below_idx]

        if abs(v2 - v1) < 1e-12:
            return t1

        return t1 + (threshold - v1) * (t2 - t1) / (v2 - v1)

    # =========================================================================
    # Arrhenius Analysis
    # =========================================================================

    def arrhenius_analysis(
        self,
        T_sweep_result: SweepResult,
        target_oit: float = 50.0,
        time_metrics: List[float] = [1.0, 3.0, 6.0]
    ) -> dict:
        """
        Perform Arrhenius analysis on temperature sweep results.

        Calculates apparent activation energy from:
        1. Time to reach target OIT vs 1/T
        2. OIT values at fixed times vs 1/T

        Args:
            T_sweep_result: SweepResult from temperature sweep
            target_oit: Target OIT value
            time_metrics: Times (months) for OIT extraction

        Returns:
            Dictionary with Arrhenius parameters and fit quality
        """
        R = 8.314  # J/(mol·K)

        T_values = T_sweep_result.parameter_values
        T_K = T_values + 273.15
        inv_T = 1.0 / T_K

        results = {
            'T_celsius': T_values,
            'T_kelvin': T_K,
            '1/T': inv_T
        }

        # Method 1: Time to target OIT
        times = T_sweep_result.times_months
        t_to_target = []

        for i in range(len(T_values)):
            oit_profile = T_sweep_result.oit_avg_matrix[i, :]
            try:
                t = self._find_time_to_threshold(times, oit_profile, target_oit)
                t_to_target.append(t)
            except:
                t_to_target.append(np.nan)

        t_to_target = np.array(t_to_target)
        results['time_to_oit_target'] = t_to_target

        # Fit ln(1/t) vs 1/T
        valid = ~np.isnan(t_to_target) & (t_to_target > 0)
        if np.sum(valid) >= 2:
            ln_rate = np.log(1.0 / t_to_target[valid])
            slope, intercept = np.polyfit(inv_T[valid], ln_rate, 1)
            Ea_time = -slope * R / 1000  # kJ/mol
            results['Ea_from_time_kJ_mol'] = Ea_time

            # Calculate R²
            fitted = slope * inv_T[valid] + intercept
            ss_res = np.sum((ln_rate - fitted)**2)
            ss_tot = np.sum((ln_rate - np.mean(ln_rate))**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            results['Ea_time_R2'] = r2
        else:
            results['Ea_from_time_kJ_mol'] = np.nan
            results['Ea_time_R2'] = np.nan

        # Method 2: OIT at fixed times
        results['oit_at_times'] = {}
        for t_month in time_metrics:
            oit_values = []
            for i in range(len(T_values)):
                oit_profile = T_sweep_result.oit_avg_matrix[i, :]
                oit_val = np.interp(t_month, times, oit_profile)
                oit_values.append(oit_val)
            results['oit_at_times'][f'{t_month}m'] = np.array(oit_values)

        return results

    # =========================================================================
    # Tornado Diagram Data
    # =========================================================================

    def calculate_tornado_data(
        self,
        parameter_sweeps: Dict[str, SweepResult],
        baseline_values: Dict[str, float],
        metric_time: float = 6.0
    ) -> pd.DataFrame:
        """
        Calculate data for tornado diagram showing parameter importance.

        Args:
            parameter_sweeps: Dict of SweepResults keyed by parameter name
            baseline_values: Dict of baseline values for each parameter
            metric_time: Time (months) at which to evaluate OIT

        Returns:
            DataFrame sorted by impact magnitude
        """
        tornado_data = []

        for param_name, sweep_result in parameter_sweeps.items():
            if param_name not in baseline_values:
                continue

            baseline_val = baseline_values[param_name]
            param_values = sweep_result.parameter_values
            times = sweep_result.times_months

            # Find baseline index
            baseline_idx = np.argmin(np.abs(param_values - baseline_val))
            baseline_oit = np.interp(metric_time, times,
                                     sweep_result.oit_avg_matrix[baseline_idx, :])

            # Find min and max impacts
            min_oit = np.inf
            max_oit = -np.inf
            min_param = baseline_val
            max_param = baseline_val

            for i, pval in enumerate(param_values):
                oit = np.interp(metric_time, times, sweep_result.oit_avg_matrix[i, :])
                if oit < min_oit:
                    min_oit = oit
                    min_param = pval
                if oit > max_oit:
                    max_oit = oit
                    max_param = pval

            impact_range = max_oit - min_oit
            impact_pct = impact_range / baseline_oit * 100 if baseline_oit > 0 else 0

            tornado_data.append({
                'parameter': param_name,
                'baseline_value': baseline_val,
                'baseline_oit': baseline_oit,
                'min_oit': min_oit,
                'max_oit': max_oit,
                'min_param_value': min_param,
                'max_param_value': max_param,
                'oit_range': impact_range,
                'impact_pct': impact_pct
            })

        df = pd.DataFrame(tornado_data)
        return df.sort_values('oit_range', ascending=True)

    # =========================================================================
    # Export Methods
    # =========================================================================

    def export_sweep_results(
        self,
        sweep_result: SweepResult,
        filepath: str,
        format: str = 'csv'
    ):
        """
        Export sweep results to file.

        Args:
            sweep_result: SweepResult to export
            filepath: Output file path
            format: 'csv' or 'json'
        """
        # Create DataFrame with all results
        data = {
            sweep_result.parameter_name: [],
            'time_months': [],
            'oit_avg': [],
            'oit_surface': [],
            'co_surface': []
        }

        for i, param_val in enumerate(sweep_result.parameter_values):
            for j, t in enumerate(sweep_result.times_months):
                data[sweep_result.parameter_name].append(param_val)
                data['time_months'].append(t)
                data['oit_avg'].append(sweep_result.oit_avg_matrix[i, j])
                data['oit_surface'].append(sweep_result.oit_surface_matrix[i, j])
                data['co_surface'].append(sweep_result.co_surface_matrix[i, j])

        df = pd.DataFrame(data)

        if format == 'csv':
            df.to_csv(filepath, index=False)
        elif format == 'json':
            df.to_json(filepath, orient='records')

        print(f"Results exported to: {filepath}")


# =============================================================================
# Convenience functions
# =============================================================================

def run_standard_sensitivity_study(
    model: Optional[DegradationModel] = None,
    doc_values: np.ndarray = np.array([0, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0]),
    T_values: np.ndarray = np.array([30, 40, 50, 60]),
    AH_values: np.ndarray = np.array([0.5, 1.0, 1.5, 2.0, 3.0]),
    base_T: float = 40.0,
    base_DOC: float = 0.05,
    base_AH: float = 1.0,
    t_end_years: float = 0.75,
    verbose: bool = True
) -> Dict[str, SweepResult]:
    """
    Run a standard sensitivity study with all three parameter sweeps.

    Args:
        model: Model to use (creates SUEZ film model if None)
        doc_values: DOC concentrations to sweep
        T_values: Temperatures to sweep
        AH_values: AH multipliers to sweep
        base_T: Baseline temperature
        base_DOC: Baseline DOC
        base_AH: Baseline AH multiplier
        t_end_years: Simulation duration
        verbose: Print progress

    Returns:
        Dictionary of SweepResults keyed by parameter name
    """
    if model is None:
        model = create_suez_film_model()

    analyzer = ParameterAnalysis(model)

    results = {}

    # DOC sweep
    if verbose:
        print("\n" + "="*70)
        print("SENSITIVITY STUDY: DOC Concentration")
        print("="*70)

    results['DOC_ppm'] = analyzer.sweep_doc_concentration(
        doc_values=doc_values,
        T_celsius=base_T,
        t_end_years=t_end_years,
        AH0_mult=base_AH,
        verbose=verbose
    )

    # Temperature sweep
    if verbose:
        print("\n" + "="*70)
        print("SENSITIVITY STUDY: Temperature")
        print("="*70)

    results['T_celsius'] = analyzer.sweep_temperature(
        T_values=T_values,
        DOC_ppm=base_DOC,
        t_end_years=t_end_years,
        AH0_mult=base_AH,
        verbose=verbose
    )

    # Antioxidant sweep
    if verbose:
        print("\n" + "="*70)
        print("SENSITIVITY STUDY: Initial Antioxidant")
        print("="*70)

    results['AH0_mult'] = analyzer.sweep_antioxidant(
        AH_multipliers=AH_values,
        T_celsius=base_T,
        DOC_ppm=base_DOC,
        t_end_years=t_end_years,
        verbose=verbose
    )

    return results
