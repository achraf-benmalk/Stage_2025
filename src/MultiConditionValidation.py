"""
Multi-Condition Validation Module

Comprehensive validation of the PE degradation model against SUEZ experimental
data and systematic assessment of model performance across conditions.

Key outputs:
- Validation metrics (R², RMSE, MAE, MAPE) for each condition
- Comparison matrices and heatmaps
- Identification of model limitations

Author: Stage 2025 Research Project
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, Tuple, List, Optional
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from DegradationModel import DegradationModel, SimulationParams, create_suez_film_model
from ValidationMetrics import (SUEZExperimentalData, calculate_r2, calculate_rmse,
                                calculate_mae, calculate_mape)


@dataclass
class ValidationResult:
    """Results from a single validation comparison."""
    condition_name: str
    r2: float
    rmse: float
    mae: float
    mape: float
    n_points: int
    times_exp: np.ndarray
    oit_exp: np.ndarray
    oit_sim: np.ndarray


class MultiConditionValidator:
    """
    Comprehensive model validation across multiple conditions.
    """

    def __init__(self):
        """Initialize validator."""
        self.model = create_suez_film_model()
        self.suez = SUEZExperimentalData()
        self.results = {}

    def validate_h2o_control(self, verbose: bool = True) -> ValidationResult:
        """
        Validate against H2O control experiment (no DOC).

        This tests the physical aging component of the model.
        """
        if verbose:
            print("\n" + "="*60)
            print("H2O CONTROL VALIDATION (Physical Aging Only)")
            print("="*60)

        # Run simulation with no DOC
        params = SimulationParams(
            T_celsius=40.0,
            DOC_ppm=0.0,
            t_end_years=1.0,
            n_timepoints=100
        )

        result = self.model.simulate(params, verbose=False)
        t_months, oit_sim = self.model.calculate_average_oit(result)

        # Interpolate to experimental times
        oit_pred = np.interp(self.suez.times_h2o, t_months, oit_sim)

        # Calculate metrics
        r2 = calculate_r2(self.suez.oit_h2o, oit_pred)
        rmse = calculate_rmse(self.suez.oit_h2o, oit_pred)
        mae = calculate_mae(self.suez.oit_h2o, oit_pred)
        mape = calculate_mape(self.suez.oit_h2o, oit_pred)

        vr = ValidationResult(
            condition_name='H2O_Control',
            r2=r2, rmse=rmse, mae=mae, mape=mape,
            n_points=len(self.suez.times_h2o),
            times_exp=self.suez.times_h2o,
            oit_exp=self.suez.oit_h2o,
            oit_sim=oit_pred
        )

        self.results['H2O_Control'] = vr

        if verbose:
            self._print_validation_result(vr)

        return vr

    def validate_hocl_exposure(self, doc_ppm: float = 0.05,
                               verbose: bool = True) -> ValidationResult:
        """
        Validate against HOCl exposure experiment.

        Args:
            doc_ppm: DOC concentration to use (simulating HOCl)
        """
        if verbose:
            print("\n" + "="*60)
            print(f"HOCl EXPOSURE VALIDATION (DOC={doc_ppm} ppm)")
            print("="*60)

        # Run simulation
        params = SimulationParams(
            T_celsius=40.0,
            DOC_ppm=doc_ppm,
            t_end_years=1.0,
            n_timepoints=100
        )

        result = self.model.simulate(params, verbose=False)
        t_months, oit_sim = self.model.calculate_average_oit(result)

        # Interpolate to experimental times
        oit_pred = np.interp(self.suez.times_hocl, t_months, oit_sim)

        # Calculate metrics
        r2 = calculate_r2(self.suez.oit_hocl, oit_pred)
        rmse = calculate_rmse(self.suez.oit_hocl, oit_pred)
        mae = calculate_mae(self.suez.oit_hocl, oit_pred)
        mape = calculate_mape(self.suez.oit_hocl, oit_pred)

        name = f'HOCl_{doc_ppm}'
        vr = ValidationResult(
            condition_name=name,
            r2=r2, rmse=rmse, mae=mae, mape=mape,
            n_points=len(self.suez.times_hocl),
            times_exp=self.suez.times_hocl,
            oit_exp=self.suez.oit_hocl,
            oit_sim=oit_pred
        )

        self.results[name] = vr

        if verbose:
            self._print_validation_result(vr)

        return vr

    def validate_calibrated_hocl(self, verbose: bool = True) -> ValidationResult:
        """
        Validate using calibrated effective DOC (from HOCl kinetics analysis).

        Optimal DOC = 0.039 ppm from calibration study.
        """
        return self.validate_hocl_exposure(doc_ppm=0.039, verbose=verbose)

    def validate_temperature_sensitivity(self,
                                         temperatures: List[float] = None,
                                         doc_ppm: float = 0.05,
                                         verbose: bool = True) -> Dict[str, ValidationResult]:
        """
        Run validation at multiple temperatures.

        Note: No experimental data for other temperatures, so this
        generates predictions for cross-validation.
        """
        if temperatures is None:
            temperatures = [30, 40, 50, 60]

        if verbose:
            print("\n" + "="*60)
            print("TEMPERATURE SENSITIVITY ANALYSIS")
            print("="*60)

        results = {}

        for T in temperatures:
            params = SimulationParams(
                T_celsius=T,
                DOC_ppm=doc_ppm,
                t_end_years=1.0,
                n_timepoints=100
            )

            result = self.model.simulate(params, verbose=False)
            t_months, oit_sim = self.model.calculate_average_oit(result)

            # Find time to reach OIT=10 min
            idx = np.argmax(oit_sim < 10)
            if idx > 0:
                lifetime = t_months[idx]
            else:
                lifetime = t_months[-1]

            name = f'T{T}C_DOC{doc_ppm}'
            results[name] = {
                't_months': t_months,
                'oit_sim': oit_sim,
                'lifetime': lifetime,
                'T': T,
                'DOC': doc_ppm
            }

            if verbose:
                print(f"  T={T}°C: Lifetime to OIT<10 = {lifetime:.2f} months")

        return results

    def _print_validation_result(self, vr: ValidationResult):
        """Print validation result summary."""
        print(f"\nCondition: {vr.condition_name}")
        print(f"  N points: {vr.n_points}")
        print(f"  R²:   {vr.r2:.4f}")
        print(f"  RMSE: {vr.rmse:.2f} min")
        print(f"  MAE:  {vr.mae:.2f} min")
        print(f"  MAPE: {vr.mape:.1f}%")

        print("\n  Point-by-point comparison:")
        print(f"  {'Time':<8} {'Exp OIT':<10} {'Sim OIT':<10} {'Error':<10} {'Error %':<10}")
        print("  " + "-"*48)
        for i in range(vr.n_points):
            t = vr.times_exp[i]
            exp = vr.oit_exp[i]
            sim = vr.oit_sim[i]
            err = sim - exp
            err_pct = 100 * err / exp if exp != 0 else 0
            print(f"  {t:<8.1f} {exp:<10.2f} {sim:<10.2f} {err:>+10.2f} {err_pct:>+10.1f}%")

    def run_comprehensive_validation(self, verbose: bool = True) -> Dict:
        """
        Run all validation tests.
        """
        if verbose:
            print("="*70)
            print("COMPREHENSIVE MODEL VALIDATION")
            print("="*70)

        # H2O control
        self.validate_h2o_control(verbose=verbose)

        # Original DOC
        self.validate_hocl_exposure(doc_ppm=0.05, verbose=verbose)

        # Calibrated DOC
        self.validate_calibrated_hocl(verbose=verbose)

        # Temperature sensitivity
        temp_results = self.validate_temperature_sensitivity(verbose=verbose)

        return {
            'validation_results': self.results,
            'temperature_sensitivity': temp_results
        }

    def plot_validation_summary(self, save_path: Optional[str] = None) -> plt.Figure:
        """
        Create validation summary figure.
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Panel A: H2O control
        ax = axes[0, 0]
        if 'H2O_Control' in self.results:
            vr = self.results['H2O_Control']
            ax.plot(vr.times_exp, vr.oit_exp, 'ko-', markersize=10,
                    linewidth=2, label='Experimental')
            ax.plot(vr.times_exp, vr.oit_sim, 'b^--', markersize=8,
                    linewidth=2, label=f'Model (R²={vr.r2:.3f})')
            ax.set_xlabel('Time (months)', fontsize=12)
            ax.set_ylabel('OIT (min)', fontsize=12)
            ax.set_title('A) H₂O Control (Physical Aging)', fontsize=12, fontweight='bold')
            ax.legend(fontsize=10)
            ax.grid(True, alpha=0.3)
            ax.set_xlim(-0.5, 10)

        # Panel B: HOCl comparison
        ax = axes[0, 1]

        # Experimental
        ax.plot(self.suez.times_hocl, self.suez.oit_hocl, 'ko-', markersize=10,
                linewidth=2, label='Experimental (HOCl 0.05 ppm)', zorder=10)

        # Original DOC
        if 'HOCl_0.05' in self.results:
            vr = self.results['HOCl_0.05']
            ax.plot(vr.times_exp, vr.oit_sim, 'r--', linewidth=2,
                    label=f'Original DOC=0.05 (R²={vr.r2:.3f})')

        # Calibrated DOC
        if 'HOCl_0.039' in self.results:
            vr = self.results['HOCl_0.039']
            ax.plot(vr.times_exp, vr.oit_sim, 'g-', linewidth=2,
                    label=f'Calibrated DOC=0.039 (R²={vr.r2:.3f})')

        ax.set_xlabel('Time (months)', fontsize=12)
        ax.set_ylabel('OIT (min)', fontsize=12)
        ax.set_title('B) HOCl Exposure Validation', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10, loc='upper right')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-0.5, 10)
        ax.set_ylim(-10, 310)

        # Panel C: Metrics comparison
        ax = axes[1, 0]

        conditions = list(self.results.keys())
        metrics = ['R²', 'RMSE', 'MAE']
        x = np.arange(len(conditions))
        width = 0.25

        r2_vals = [self.results[c].r2 for c in conditions]
        rmse_vals = [self.results[c].rmse / 100 for c in conditions]  # Normalize
        mae_vals = [self.results[c].mae / 100 for c in conditions]

        ax.bar(x - width, r2_vals, width, label='R²', color='blue', alpha=0.7)
        ax.bar(x, rmse_vals, width, label='RMSE/100', color='red', alpha=0.7)
        ax.bar(x + width, mae_vals, width, label='MAE/100', color='green', alpha=0.7)

        ax.set_xlabel('Condition', fontsize=12)
        ax.set_ylabel('Metric Value', fontsize=12)
        ax.set_title('C) Validation Metrics Comparison', fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels([c.replace('_', '\n') for c in conditions], fontsize=9)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, axis='y')

        # Panel D: Residual analysis
        ax = axes[1, 1]

        colors = ['blue', 'red', 'green']
        for i, (name, vr) in enumerate(self.results.items()):
            errors = vr.oit_sim - vr.oit_exp
            rel_errors = 100 * errors / vr.oit_exp

            ax.scatter(vr.times_exp, rel_errors, c=colors[i % 3],
                       s=80, alpha=0.7, label=name, edgecolors='k')

        ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
        ax.axhline(y=20, color='gray', linestyle='--', alpha=0.5)
        ax.axhline(y=-20, color='gray', linestyle='--', alpha=0.5)

        ax.set_xlabel('Time (months)', fontsize=12)
        ax.set_ylabel('Relative Error (%)', fontsize=12)
        ax.set_title('D) Residual Analysis', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9, loc='lower left')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-0.5, 10)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")

        return fig


# =============================================================================
# Main execution
# =============================================================================

if __name__ == '__main__':
    print("MULTI-CONDITION VALIDATION ANALYSIS")
    print()

    validator = MultiConditionValidator()

    # Run comprehensive validation
    results = validator.run_comprehensive_validation()

    # Generate summary figure
    fig = validator.plot_validation_summary(
        save_path='/home/user/Stage_2025/figures/validation_summary.png'
    )
    plt.savefig('/home/user/Stage_2025/figures/validation_summary.pdf',
                dpi=300, bbox_inches='tight')
    plt.close()

    # Summary table
    print()
    print("="*70)
    print("VALIDATION SUMMARY TABLE")
    print("="*70)
    print()
    print(f"{'Condition':<20} {'R²':<10} {'RMSE (min)':<12} {'MAE (min)':<12} {'MAPE (%)':<10}")
    print("-"*64)

    for name, vr in validator.results.items():
        print(f"{name:<20} {vr.r2:<10.4f} {vr.rmse:<12.2f} {vr.mae:<12.2f} {vr.mape:<10.1f}")

    print()
    print("="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
