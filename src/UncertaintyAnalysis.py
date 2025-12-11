"""
Monte Carlo Uncertainty Analysis Module

This module implements Monte Carlo sampling for uncertainty quantification
of the PE degradation model predictions.

Key Features:
- Latin Hypercube Sampling (LHS) for efficient parameter space coverage
- Configurable parameter distributions (normal, uniform, lognormal)
- Parallel execution support
- Statistical output (percentiles, confidence intervals)
- Publication-quality uncertainty band figures

Author: Stage 2025 Research Project
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, uniform, lognorm, qmc
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Callable
import time
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from DegradationModel import DegradationModel, SimulationParams, create_suez_film_model


@dataclass
class ParameterDistribution:
    """Define a parameter distribution for MC sampling."""
    name: str
    nominal: float
    distribution: str  # 'normal', 'uniform', 'lognormal'
    cv: float = 0.1  # Coefficient of variation (for normal/lognormal)
    bounds: Tuple[float, float] = None  # For uniform distribution

    def sample(self, n: int, rng: np.random.Generator = None) -> np.ndarray:
        """Generate n samples from the distribution."""
        if rng is None:
            rng = np.random.default_rng()

        if self.distribution == 'normal':
            sigma = self.nominal * self.cv
            samples = rng.normal(self.nominal, sigma, n)
            # Ensure positive values
            samples = np.clip(samples, self.nominal * 0.1, self.nominal * 3.0)

        elif self.distribution == 'uniform':
            if self.bounds:
                low, high = self.bounds
            else:
                low = self.nominal * (1 - self.cv * 2)
                high = self.nominal * (1 + self.cv * 2)
            samples = rng.uniform(low, high, n)

        elif self.distribution == 'lognormal':
            # Lognormal with specified CV
            sigma_ln = np.sqrt(np.log(1 + self.cv**2))
            mu_ln = np.log(self.nominal) - sigma_ln**2 / 2
            samples = rng.lognormal(mu_ln, sigma_ln, n)

        else:
            raise ValueError(f"Unknown distribution: {self.distribution}")

        return samples


@dataclass
class MonteCarloConfig:
    """Configuration for Monte Carlo analysis."""
    n_samples: int = 500
    seed: int = 42
    use_lhs: bool = True  # Latin Hypercube Sampling
    parameters: List[ParameterDistribution] = field(default_factory=list)

    # Simulation settings
    T_celsius: float = 40.0
    DOC_ppm: float = 0.05
    t_end_years: float = 0.75
    n_timepoints: int = 50


def create_default_parameter_distributions() -> List[ParameterDistribution]:
    """
    Create default parameter distributions based on Colin et al. uncertainty estimates.

    Parameters are represented as SCALE FACTORS centered around 1.0.
    A scale factor of 1.0 means nominal value, 0.5 means half, 2.0 means double.

    Conservative ranges to ensure numerical stability.

    Returns:
        List of ParameterDistribution objects (scale factors)
    """
    return [
        # Most sensitive kinetic parameters (conservative uncertainty)
        ParameterDistribution('k1d', 1.0, 'uniform', bounds=(0.7, 1.4)),   # DOC initiation ±30%
        ParameterDistribution('k3', 1.0, 'uniform', bounds=(0.8, 1.25)),   # Propagation ±20%
        ParameterDistribution('k7', 1.0, 'uniform', bounds=(0.7, 1.4)),    # AH stabilization ±30%

        # Diffusion coefficients (bounded for stability)
        ParameterDistribution('D_O2', 1.0, 'uniform', bounds=(0.8, 1.25)),
        ParameterDistribution('D_DOC', 1.0, 'uniform', bounds=(0.7, 1.4)),

        # Initial concentrations (well-characterized)
        ParameterDistribution('AH0', 1.0, 'uniform', bounds=(0.85, 1.15)),
    ]


class MonteCarloAnalysis:
    """
    Monte Carlo uncertainty analysis for PE degradation model.

    Samples parameter space and propagates uncertainty through simulations.
    """

    def __init__(self, config: Optional[MonteCarloConfig] = None):
        """
        Initialize Monte Carlo analysis.

        Args:
            config: MonteCarloConfig object (uses defaults if None)
        """
        self.config = config or MonteCarloConfig()
        if not self.config.parameters:
            self.config.parameters = create_default_parameter_distributions()

        self.rng = np.random.default_rng(self.config.seed)
        self.samples = None
        self.results = None
        self.base_model = create_suez_film_model()

    def generate_samples(self) -> np.ndarray:
        """
        Generate parameter samples using LHS or random sampling.

        Returns:
            Array of shape (n_samples, n_parameters)
        """
        n = self.config.n_samples
        n_params = len(self.config.parameters)

        if self.config.use_lhs:
            # Latin Hypercube Sampling for better coverage
            sampler = qmc.LatinHypercube(d=n_params, seed=self.config.seed)
            unit_samples = sampler.random(n)

            # Transform to parameter distributions
            samples = np.zeros((n, n_params))
            for i, param in enumerate(self.config.parameters):
                if param.distribution == 'normal':
                    sigma = param.nominal * param.cv
                    samples[:, i] = norm.ppf(unit_samples[:, i], loc=param.nominal, scale=sigma)
                    samples[:, i] = np.clip(samples[:, i], param.nominal * 0.1, param.nominal * 3.0)

                elif param.distribution == 'uniform':
                    if param.bounds:
                        low, high = param.bounds
                    else:
                        low = param.nominal * (1 - param.cv * 2)
                        high = param.nominal * (1 + param.cv * 2)
                    samples[:, i] = uniform.ppf(unit_samples[:, i], loc=low, scale=high-low)

                elif param.distribution == 'lognormal':
                    sigma_ln = np.sqrt(np.log(1 + param.cv**2))
                    mu_ln = np.log(param.nominal) - sigma_ln**2 / 2
                    samples[:, i] = lognorm.ppf(unit_samples[:, i], s=sigma_ln, scale=np.exp(mu_ln))

        else:
            # Simple random sampling
            samples = np.zeros((n, n_params))
            for i, param in enumerate(self.config.parameters):
                samples[:, i] = param.sample(n, self.rng)

        self.samples = samples
        return samples

    def _create_modified_model(self, param_values: np.ndarray) -> DegradationModel:
        """
        Create a model with modified parameters.

        The model uses:
        - k_coeffs: dict of (Ea, A0) tuples for rate constants
        - D_coeffs: dict of (Ea, D0) tuples for diffusion coefficients
        - AH0_conc, POOH0_conc: initial concentrations

        We vary the pre-exponential factors (A0) to represent parameter uncertainty.

        Args:
            param_values: Array of parameter values (scale factors)

        Returns:
            Modified DegradationModel
        """
        model = create_suez_film_model()

        # Map parameter names to model structure
        for i, param in enumerate(self.config.parameters):
            scale_factor = param_values[i]

            # Kinetic rate constants - modify pre-exponential factor A0
            if param.name in model.k_coeffs:
                Ea, A0 = model.k_coeffs[param.name]
                model.k_coeffs[param.name] = (Ea, A0 * scale_factor)

            # Diffusion coefficients - modify pre-exponential D0
            elif param.name == 'D_O2':
                Ea, D0 = model.D_coeffs['O2']
                model.D_coeffs['O2'] = (Ea, D0 * scale_factor)
            elif param.name == 'D_DOC':
                Ea, D0 = model.D_coeffs['DOC']
                model.D_coeffs['DOC'] = (Ea, D0 * scale_factor)
            elif param.name == 'D_AH':
                Ea, D0 = model.D_coeffs['AH']
                model.D_coeffs['AH'] = (Ea, D0 * scale_factor)

            # Initial concentrations - directly modify
            elif param.name == 'AH0':
                model.AH0_conc = model.AH0_conc * scale_factor
            elif param.name == 'POOH0':
                model.POOH0_conc = model.POOH0_conc * scale_factor

        return model

    def run_simulation(self, sample_idx: int) -> Tuple[np.ndarray, np.ndarray, bool]:
        """
        Run a single simulation with sampled parameters.

        Uses LSODA solver with loose tolerances for speed.

        Args:
            sample_idx: Index into samples array

        Returns:
            (times_months, oit_avg, success)
        """
        param_values = self.samples[sample_idx]
        model = self._create_modified_model(param_values)

        # Use BDF solver (good for stiff problems) with moderate tolerances
        params = SimulationParams(
            T_celsius=self.config.T_celsius,
            DOC_ppm=self.config.DOC_ppm,
            t_end_years=self.config.t_end_years,
            n_timepoints=30,  # Fewer output points
            method='BDF',     # Backward differentiation formula (stiff)
            rtol=1e-5,        # Moderate tolerance
            atol=1e-8
        )

        try:
            result = model.simulate(params, verbose=False)
            if result['success']:
                times_m, oit_avg = model.calculate_average_oit(result)
                return times_m, oit_avg, True
        except Exception as e:
            pass

        return None, None, False

    def run_analysis(self, verbose: bool = True) -> Dict:
        """
        Run the complete Monte Carlo analysis.

        Args:
            verbose: Print progress updates

        Returns:
            Dictionary with results
        """
        if verbose:
            print("="*70)
            print("MONTE CARLO UNCERTAINTY ANALYSIS")
            print("="*70)
            print(f"N samples: {self.config.n_samples}")
            print(f"Sampling method: {'Latin Hypercube' if self.config.use_lhs else 'Random'}")
            print(f"Parameters varied: {len(self.config.parameters)}")
            print()

        # Generate samples if not done
        if self.samples is None:
            if verbose:
                print("Generating parameter samples...")
            self.generate_samples()

        # Run simulations
        if verbose:
            print(f"Running {self.config.n_samples} simulations...")

        start_time = time.time()

        # Store results on common time grid
        t_common = np.linspace(0, self.config.t_end_years * 12, self.config.n_timepoints)
        oit_matrix = np.full((self.config.n_samples, len(t_common)), np.nan)

        n_success = 0
        for i in range(self.config.n_samples):
            times_m, oit_avg, success = self.run_simulation(i)

            if success and times_m is not None:
                # Interpolate to common grid
                oit_matrix[i, :] = np.interp(t_common, times_m, oit_avg,
                                             left=oit_avg[0], right=oit_avg[-1])
                n_success += 1

            if verbose and (i + 1) % 50 == 0:
                elapsed = time.time() - start_time
                rate = (i + 1) / elapsed
                eta = (self.config.n_samples - i - 1) / rate
                print(f"  Progress: {i+1}/{self.config.n_samples} ({n_success} successful) - ETA: {eta:.0f}s")

        elapsed_total = time.time() - start_time

        if verbose:
            print(f"\nCompleted in {elapsed_total:.1f}s")
            print(f"Successful simulations: {n_success}/{self.config.n_samples} ({100*n_success/self.config.n_samples:.1f}%)")

        # Calculate statistics
        if verbose:
            print("\nCalculating statistics...")

        # Mask failed simulations
        valid_mask = ~np.isnan(oit_matrix[:, 0])
        oit_valid = oit_matrix[valid_mask, :]

        percentiles = [5, 10, 25, 50, 75, 90, 95]
        stats = {}

        stats['times'] = t_common
        stats['mean'] = np.nanmean(oit_valid, axis=0)
        stats['std'] = np.nanstd(oit_valid, axis=0)
        stats['cv'] = stats['std'] / np.maximum(stats['mean'], 1e-10)

        for p in percentiles:
            stats[f'p{p}'] = np.nanpercentile(oit_valid, p, axis=0)

        # Store results
        self.results = {
            'config': self.config,
            'samples': self.samples,
            'oit_matrix': oit_matrix,
            'valid_mask': valid_mask,
            'n_success': n_success,
            'stats': stats,
            'elapsed_time': elapsed_total
        }

        if verbose:
            self._print_summary()

        return self.results

    def _print_summary(self):
        """Print summary statistics."""
        stats = self.results['stats']

        print()
        print("="*70)
        print("UNCERTAINTY ANALYSIS RESULTS")
        print("="*70)

        # Key time points
        time_points = [0, 1, 2, 3, 6, 9]

        print(f"\n{'Month':<8} {'Mean OIT':<12} {'Std':<10} {'CV (%)':<10} {'90% CI':<20}")
        print("-"*60)

        for month in time_points:
            idx = np.argmin(np.abs(stats['times'] - month))
            mean = stats['mean'][idx]
            std = stats['std'][idx]
            cv = stats['cv'][idx] * 100
            p5 = stats['p5'][idx]
            p95 = stats['p95'][idx]

            print(f"{month:<8} {mean:<12.2f} {std:<10.2f} {cv:<10.1f} [{p5:.1f}, {p95:.1f}]")

        print()
        print(f"Total simulations: {self.config.n_samples}")
        print(f"Successful: {self.results['n_success']}")
        print(f"Parameters varied: {len(self.config.parameters)}")

    def plot_uncertainty_bands(self, save_path: Optional[str] = None) -> plt.Figure:
        """
        Create publication-quality uncertainty band plot.

        Args:
            save_path: Path to save figure (optional)

        Returns:
            matplotlib Figure
        """
        if self.results is None:
            raise ValueError("Must run analysis first")

        stats = self.results['stats']

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Panel A: Full uncertainty bands
        ax = axes[0]

        t = stats['times']

        # 90% CI band
        ax.fill_between(t, stats['p5'], stats['p95'],
                        alpha=0.2, color='blue', label='90% CI')
        # 50% CI band
        ax.fill_between(t, stats['p25'], stats['p75'],
                        alpha=0.4, color='blue', label='50% CI')
        # Median
        ax.plot(t, stats['p50'], 'b-', linewidth=2, label='Median')
        # Mean
        ax.plot(t, stats['mean'], 'b--', linewidth=1.5, label='Mean')

        ax.set_xlabel('Time (months)', fontsize=12)
        ax.set_ylabel('OIT (minutes)', fontsize=12)
        ax.set_title('A) Monte Carlo Uncertainty Bands', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10, loc='upper right')
        ax.set_xlim(0, 9)
        ax.set_ylim(0, 320)
        ax.grid(True, alpha=0.3)

        # Annotation
        ax.annotate(f'N = {self.results["n_success"]} simulations',
                    xy=(0.05, 0.05), xycoords='axes fraction',
                    fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Panel B: CV evolution
        ax = axes[1]

        ax.plot(t, stats['cv'] * 100, 'r-', linewidth=2)
        ax.set_xlabel('Time (months)', fontsize=12)
        ax.set_ylabel('Coefficient of Variation (%)', fontsize=12)
        ax.set_title('B) Prediction Uncertainty Over Time', fontsize=12, fontweight='bold')
        ax.set_xlim(0, 9)
        ax.set_ylim(0, max(stats['cv'] * 100) * 1.1)
        ax.grid(True, alpha=0.3)

        # Add horizontal reference lines
        ax.axhline(y=10, color='green', linestyle='--', alpha=0.7, label='10% CV')
        ax.axhline(y=20, color='orange', linestyle='--', alpha=0.7, label='20% CV')
        ax.axhline(y=50, color='red', linestyle='--', alpha=0.7, label='50% CV')
        ax.legend(fontsize=10, loc='upper left')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")

        return fig

    def plot_parameter_sensitivity(self, save_path: Optional[str] = None) -> plt.Figure:
        """
        Create parameter sensitivity plot showing correlation with OIT.

        Args:
            save_path: Path to save figure

        Returns:
            matplotlib Figure
        """
        if self.results is None:
            raise ValueError("Must run analysis first")

        # Get OIT at month 6
        t_target = 6.0
        t_idx = np.argmin(np.abs(self.results['stats']['times'] - t_target))
        oit_at_6m = self.results['oit_matrix'][:, t_idx]

        valid = ~np.isnan(oit_at_6m)
        oit_valid = oit_at_6m[valid]
        samples_valid = self.samples[valid, :]

        # Calculate correlations
        correlations = []
        for i, param in enumerate(self.config.parameters):
            corr = np.corrcoef(samples_valid[:, i], oit_valid)[0, 1]
            correlations.append((param.name, corr))

        # Sort by absolute correlation
        correlations.sort(key=lambda x: abs(x[1]), reverse=True)

        # Plot
        fig, ax = plt.subplots(figsize=(10, 6))

        names = [c[0] for c in correlations]
        values = [c[1] for c in correlations]
        colors = ['red' if v < 0 else 'blue' for v in values]

        y_pos = np.arange(len(names))
        ax.barh(y_pos, values, color=colors, alpha=0.7)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(names)
        ax.set_xlabel('Correlation with OIT at 6 months', fontsize=12)
        ax.set_title('Parameter Sensitivity (Correlation Analysis)', fontsize=12, fontweight='bold')
        ax.axvline(x=0, color='black', linewidth=0.5)
        ax.set_xlim(-1, 1)
        ax.grid(True, alpha=0.3, axis='x')

        # Add value labels
        for i, (name, val) in enumerate(correlations):
            x_pos = val + 0.02 if val >= 0 else val - 0.02
            ha = 'left' if val >= 0 else 'right'
            ax.text(x_pos, i, f'{val:.2f}', va='center', ha=ha, fontsize=9)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")

        return fig


def run_quick_mc_analysis(n_samples: int = 100, verbose: bool = True) -> Dict:
    """
    Run a quick Monte Carlo analysis with reduced samples.

    Args:
        n_samples: Number of samples (default 100)
        verbose: Print progress

    Returns:
        Analysis results dictionary
    """
    config = MonteCarloConfig(
        n_samples=n_samples,
        seed=42,
        use_lhs=True,
        T_celsius=40.0,
        DOC_ppm=0.05,
        t_end_years=0.75,
        n_timepoints=50
    )

    mc = MonteCarloAnalysis(config)
    results = mc.run_analysis(verbose=verbose)

    return results, mc


# =============================================================================
# Main execution
# =============================================================================

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Monte Carlo Uncertainty Analysis')
    parser.add_argument('-n', '--samples', type=int, default=500,
                        help='Number of MC samples')
    parser.add_argument('-q', '--quick', action='store_true',
                        help='Quick run with 100 samples')
    args = parser.parse_args()

    n_samples = 100 if args.quick else args.samples

    print(f"Running Monte Carlo analysis with {n_samples} samples...")

    config = MonteCarloConfig(
        n_samples=n_samples,
        seed=42,
        use_lhs=True,
        T_celsius=40.0,
        DOC_ppm=0.05,
        t_end_years=0.75
    )

    mc = MonteCarloAnalysis(config)
    results = mc.run_analysis(verbose=True)

    # Generate figures
    print("\nGenerating figures...")

    fig_dir = '/home/user/Stage_2025/figures'

    mc.plot_uncertainty_bands(save_path=f'{fig_dir}/mc_uncertainty_bands.png')
    mc.plot_parameter_sensitivity(save_path=f'{fig_dir}/mc_parameter_sensitivity.png')

    # Save PDF versions
    mc.plot_uncertainty_bands(save_path=f'{fig_dir}/mc_uncertainty_bands.pdf')
    mc.plot_parameter_sensitivity(save_path=f'{fig_dir}/mc_parameter_sensitivity.pdf')

    plt.close('all')

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
