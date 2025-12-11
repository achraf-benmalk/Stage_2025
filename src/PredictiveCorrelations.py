"""
Predictive Correlations Module

Fits empirical lifetime equations from simulation data for engineering predictions.
Implements Arrhenius-type temperature dependence and power-law concentration scaling.

Key outputs:
- Lifetime prediction equations: t_lifetime = f(DOC, T)
- Extrapolation capability to new conditions
- Confidence intervals on predictions

Author: Stage 2025 Research Project
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import t as t_dist
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, List
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from DegradationModel import DegradationModel, SimulationParams, create_suez_film_model


@dataclass
class LifetimeCorrelation:
    """Results from lifetime correlation fitting."""
    # Model form: t_life = A * exp(Ea/(R*T)) * C^(-n)
    A: float            # Pre-exponential factor
    Ea: float           # Activation energy (J/mol)
    n: float            # Concentration exponent

    # Statistics
    r2: float           # Coefficient of determination
    rmse: float         # Root mean square error
    n_points: int       # Number of data points

    # Standard errors
    A_se: float = 0.0
    Ea_se: float = 0.0
    n_se: float = 0.0


class PredictiveCorrelations:
    """
    Fit and use lifetime prediction correlations.

    Uses simulation data to derive empirical equations relating
    lifetime to DOC concentration and temperature.
    """

    R_GAS = 8.314  # J/(mol·K)

    def __init__(self, oit_threshold: float = 10.0):
        """
        Initialize correlations module.

        Args:
            oit_threshold: OIT value (min) defining end-of-life
        """
        self.oit_threshold = oit_threshold
        self.model = create_suez_film_model()
        self.data = None
        self.correlation = None

    def generate_lifetime_data(self,
                               doc_values: List[float] = None,
                               temp_values: List[float] = None,
                               verbose: bool = True) -> Dict:
        """
        Generate lifetime data from simulations.

        Args:
            doc_values: DOC concentrations (ppm)
            temp_values: Temperatures (°C)
            verbose: Print progress

        Returns:
            Dictionary with arrays of DOC, T, lifetime data
        """
        if doc_values is None:
            doc_values = [0.01, 0.02, 0.03, 0.05, 0.075, 0.1]
        if temp_values is None:
            temp_values = [30, 40, 50, 60]

        if verbose:
            print("="*70)
            print("GENERATING LIFETIME DATA")
            print("="*70)
            print(f"DOC values: {doc_values}")
            print(f"Temperatures: {temp_values}")
            print(f"OIT threshold: {self.oit_threshold} min")
            print()

        results = {
            'DOC': [],
            'T_celsius': [],
            'T_kelvin': [],
            'lifetime_months': [],
            'lifetime_years': [],
            'final_oit': []
        }

        total = len(doc_values) * len(temp_values)
        count = 0

        for doc in doc_values:
            for T in temp_values:
                count += 1
                if verbose:
                    print(f"  [{count}/{total}] DOC={doc:.3f} ppm, T={T}°C... ", end='')

                lifetime = self._estimate_lifetime(doc, T)

                if lifetime is not None:
                    results['DOC'].append(doc)
                    results['T_celsius'].append(T)
                    results['T_kelvin'].append(T + 273.15)
                    results['lifetime_months'].append(lifetime)
                    results['lifetime_years'].append(lifetime / 12.0)

                    if verbose:
                        print(f"lifetime = {lifetime:.2f} months")
                else:
                    if verbose:
                        print("FAILED")

        # Convert to arrays
        for key in results:
            results[key] = np.array(results[key])

        self.data = results

        if verbose:
            print()
            print(f"Generated {len(results['DOC'])} data points")

        return results

    def _estimate_lifetime(self, doc_ppm: float, T_celsius: float,
                           max_years: float = 10.0) -> Optional[float]:
        """
        Estimate lifetime for given conditions via binary search.

        Args:
            doc_ppm: DOC concentration
            T_celsius: Temperature
            max_years: Maximum simulation time

        Returns:
            Lifetime in months, or None if failed
        """
        # First run to max time to check if OIT threshold is reached
        params = SimulationParams(
            T_celsius=T_celsius,
            DOC_ppm=doc_ppm,
            t_end_years=max_years,
            n_timepoints=100,
            method='BDF',
            rtol=1e-5,
            atol=1e-8
        )

        try:
            result = self.model.simulate(params, verbose=False)
            if not result['success']:
                return None

            times_m, oit_avg = self.model.calculate_average_oit(result)

            # Check if threshold reached
            if oit_avg[-1] > self.oit_threshold:
                # OIT never drops below threshold
                return max_years * 12  # Return max time

            # Interpolate to find threshold crossing
            idx = np.argmax(oit_avg < self.oit_threshold)
            if idx == 0:
                return 0.0

            # Linear interpolation
            t1, t2 = times_m[idx-1], times_m[idx]
            o1, o2 = oit_avg[idx-1], oit_avg[idx]

            if o1 != o2:
                t_life = t1 + (t2 - t1) * (o1 - self.oit_threshold) / (o1 - o2)
            else:
                t_life = t1

            return t_life

        except Exception as e:
            return None

    def fit_correlation(self, verbose: bool = True) -> LifetimeCorrelation:
        """
        Fit lifetime correlation of form:
            t_life = A * exp(Ea/(R*T)) * C^(-n)

        Or in log form:
            ln(t_life) = ln(A) + Ea/(R*T) - n*ln(C)

        Args:
            verbose: Print fitting results

        Returns:
            LifetimeCorrelation object with fitted parameters
        """
        if self.data is None:
            raise ValueError("Must generate data first using generate_lifetime_data()")

        if verbose:
            print()
            print("="*70)
            print("FITTING LIFETIME CORRELATION")
            print("="*70)

        # Extract data
        T_K = self.data['T_kelvin']
        C = self.data['DOC']
        t_life = self.data['lifetime_months']

        # Remove any zero or negative lifetimes
        valid = t_life > 0
        T_K = T_K[valid]
        C = C[valid]
        t_life = t_life[valid]

        # Log transform
        ln_t = np.log(t_life)
        inv_T = 1.0 / T_K
        ln_C = np.log(C)

        # Linear regression: ln(t) = a0 + a1/T + a2*ln(C)
        # Design matrix
        X = np.column_stack([np.ones_like(T_K), inv_T, ln_C])

        # Least squares
        coeffs, residuals, rank, s = np.linalg.lstsq(X, ln_t, rcond=None)

        ln_A = coeffs[0]
        Ea_over_R = coeffs[1]
        n = -coeffs[2]

        A = np.exp(ln_A)
        Ea = Ea_over_R * self.R_GAS

        # Predictions and statistics
        ln_t_pred = X @ coeffs
        t_pred = np.exp(ln_t_pred)

        ss_res = np.sum((ln_t - ln_t_pred)**2)
        ss_tot = np.sum((ln_t - np.mean(ln_t))**2)
        r2 = 1 - ss_res / ss_tot

        rmse = np.sqrt(np.mean((t_life - t_pred)**2))

        # Standard errors (approximate)
        n_pts = len(t_life)
        n_params = 3
        mse = ss_res / (n_pts - n_params)

        try:
            cov = mse * np.linalg.inv(X.T @ X)
            se = np.sqrt(np.diag(cov))
            A_se = A * se[0]  # Delta method
            Ea_se = se[1] * self.R_GAS
            n_se = se[2]
        except:
            A_se = Ea_se = n_se = np.nan

        self.correlation = LifetimeCorrelation(
            A=A, Ea=Ea, n=n,
            r2=r2, rmse=rmse, n_points=n_pts,
            A_se=A_se, Ea_se=Ea_se, n_se=n_se
        )

        if verbose:
            print()
            print("Fitted model: t_life = A × exp(Ea/(R×T)) × C^(-n)")
            print()
            print(f"  A  = {A:.4e} ± {A_se:.2e} months")
            print(f"  Ea = {Ea/1000:.2f} ± {Ea_se/1000:.2f} kJ/mol")
            print(f"  n  = {n:.3f} ± {n_se:.3f}")
            print()
            print(f"  R² = {r2:.4f}")
            print(f"  RMSE = {rmse:.2f} months")
            print(f"  N = {n_pts} data points")
            print()

            # Engineering form
            print("Engineering equation:")
            print(f"  t_life (months) = {A:.4e} × exp({Ea/1000:.1f}×1000/(8.314×T_K)) × DOC^(-{n:.2f})")

        return self.correlation

    def predict_lifetime(self, doc_ppm: float, T_celsius: float) -> float:
        """
        Predict lifetime using fitted correlation.

        Args:
            doc_ppm: DOC concentration (ppm)
            T_celsius: Temperature (°C)

        Returns:
            Predicted lifetime in months
        """
        if self.correlation is None:
            raise ValueError("Must fit correlation first")

        T_K = T_celsius + 273.15
        c = self.correlation

        t_life = c.A * np.exp(c.Ea / (self.R_GAS * T_K)) * doc_ppm**(-c.n)

        return t_life

    def plot_correlation_fit(self, save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot correlation fit quality.

        Args:
            save_path: Path to save figure

        Returns:
            matplotlib Figure
        """
        if self.data is None or self.correlation is None:
            raise ValueError("Must generate data and fit correlation first")

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Panel A: Parity plot
        ax = axes[0, 0]
        t_actual = self.data['lifetime_months']
        t_pred = np.array([self.predict_lifetime(d, t)
                          for d, t in zip(self.data['DOC'], self.data['T_celsius'])])

        ax.scatter(t_actual, t_pred, c=self.data['T_celsius'], cmap='coolwarm',
                   s=80, edgecolors='k', alpha=0.7)

        lims = [0, max(t_actual.max(), t_pred.max()) * 1.1]
        ax.plot(lims, lims, 'k--', label='Perfect fit')
        ax.set_xlabel('Simulated Lifetime (months)', fontsize=12)
        ax.set_ylabel('Predicted Lifetime (months)', fontsize=12)
        ax.set_title(f'A) Parity Plot (R² = {self.correlation.r2:.4f})',
                     fontsize=12, fontweight='bold')
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.grid(True, alpha=0.3)

        cbar = plt.colorbar(ax.collections[0], ax=ax)
        cbar.set_label('Temperature (°C)')

        # Panel B: Lifetime vs DOC (log-log)
        ax = axes[0, 1]

        temps = np.unique(self.data['T_celsius'])
        colors = plt.cm.coolwarm(np.linspace(0, 1, len(temps)))

        for T, color in zip(temps, colors):
            mask = self.data['T_celsius'] == T
            ax.scatter(self.data['DOC'][mask], t_actual[mask],
                       c=[color], s=80, edgecolors='k', label=f'{T}°C data')

            # Prediction line
            doc_range = np.linspace(0.01, 0.12, 50)
            t_line = [self.predict_lifetime(d, T) for d in doc_range]
            ax.plot(doc_range, t_line, c=color, linestyle='--', alpha=0.7)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('DOC Concentration (ppm)', fontsize=12)
        ax.set_ylabel('Lifetime (months)', fontsize=12)
        ax.set_title('B) Lifetime vs DOC Concentration', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9, loc='upper right')
        ax.grid(True, alpha=0.3, which='both')

        # Panel C: Arrhenius plot
        ax = axes[1, 0]

        # Group by DOC
        docs = np.unique(self.data['DOC'])
        colors = plt.cm.viridis(np.linspace(0, 1, len(docs)))

        for doc, color in zip(docs, colors):
            mask = self.data['DOC'] == doc
            inv_T = 1000 / self.data['T_kelvin'][mask]  # 1000/T for readability
            ln_t = np.log(t_actual[mask])

            ax.scatter(inv_T, ln_t, c=[color], s=80, edgecolors='k',
                       label=f'{doc:.3f} ppm')

        ax.set_xlabel('1000/T (K⁻¹)', fontsize=12)
        ax.set_ylabel('ln(Lifetime)', fontsize=12)
        ax.set_title(f'C) Arrhenius Plot (Ea = {self.correlation.Ea/1000:.1f} kJ/mol)',
                     fontsize=12, fontweight='bold')
        ax.legend(fontsize=9, loc='upper left', title='DOC')
        ax.grid(True, alpha=0.3)

        # Panel D: Residuals
        ax = axes[1, 1]

        residuals = t_actual - t_pred
        rel_err = 100 * residuals / t_actual

        ax.scatter(t_pred, rel_err, c=self.data['T_celsius'], cmap='coolwarm',
                   s=80, edgecolors='k', alpha=0.7)
        ax.axhline(y=0, color='k', linestyle='--')
        ax.axhline(y=10, color='gray', linestyle=':', alpha=0.7)
        ax.axhline(y=-10, color='gray', linestyle=':', alpha=0.7)

        ax.set_xlabel('Predicted Lifetime (months)', fontsize=12)
        ax.set_ylabel('Relative Error (%)', fontsize=12)
        ax.set_title('D) Prediction Residuals', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)

        cbar = plt.colorbar(ax.collections[0], ax=ax)
        cbar.set_label('Temperature (°C)')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")

        return fig


# =============================================================================
# Main execution
# =============================================================================

if __name__ == '__main__':
    print("PREDICTIVE CORRELATIONS - LIFETIME EQUATIONS")
    print()

    # Create correlations object
    pc = PredictiveCorrelations(oit_threshold=10.0)

    # Generate lifetime data
    data = pc.generate_lifetime_data(
        doc_values=[0.01, 0.02, 0.03, 0.05, 0.075, 0.1],
        temp_values=[30, 40, 50, 60]
    )

    # Fit correlation
    corr = pc.fit_correlation()

    # Generate plot
    fig = pc.plot_correlation_fit(
        save_path='/home/user/Stage_2025/figures/lifetime_correlation.png'
    )
    plt.savefig('/home/user/Stage_2025/figures/lifetime_correlation.pdf',
                dpi=300, bbox_inches='tight')
    plt.close()

    # Test predictions
    print()
    print("="*70)
    print("TEST PREDICTIONS")
    print("="*70)
    print()
    print("Predicting lifetime for various conditions:")
    print(f"{'DOC (ppm)':<12} {'T (°C)':<10} {'Lifetime':<15} {'Unit':<10}")
    print("-"*50)

    test_conditions = [
        (0.05, 40),   # SUEZ condition
        (0.05, 25),   # Room temperature
        (0.01, 40),   # Low DOC
        (0.1, 40),    # High DOC
        (0.05, 60),   # High temperature
    ]

    for doc, T in test_conditions:
        t_months = pc.predict_lifetime(doc, T)
        if t_months < 12:
            print(f"{doc:<12.3f} {T:<10.0f} {t_months:<15.1f} {'months':<10}")
        else:
            print(f"{doc:<12.3f} {T:<10.0f} {t_months/12:<15.1f} {'years':<10}")

    print()
    print("="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
