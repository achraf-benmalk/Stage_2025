"""
Stage 3: Chemo-Mechanical Coupling Framework

Connects chemical degradation (AH depletion, chain scission) to mechanical
property evolution for lifetime prediction under combined aging.

Key components:
1. Molecular weight (Mw) evolution from chain scission
2. Mechanical property correlations (yield stress, elongation at break)
3. Critical property thresholds for failure criteria

Reference:
- Colin et al. (2009) PE degradation kinetics
- Fayolle et al. (2007) Mw-mechanical property relationships

Author: Stage 2025 Research Project
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, List
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from DegradationModel import DegradationModel, SimulationParams, create_suez_film_model


@dataclass
class MechanicalProperties:
    """Mechanical properties of PE."""
    # Virgin material properties (PE100 typical)
    sigma_y0: float = 24.0     # Initial yield stress (MPa)
    epsilon_b0: float = 600.0  # Initial elongation at break (%)
    E0: float = 900.0          # Initial Young's modulus (MPa)

    # Critical thresholds for failure
    sigma_y_crit: float = 15.0   # Critical yield stress (MPa)
    epsilon_b_crit: float = 50.0 # Critical elongation (%) - ductile-brittle transition

    # Mw dependence parameters (Fayolle et al. 2007)
    Mw0: float = 150.0           # Initial Mw (kg/mol)
    Mw_crit: float = 40.0        # Critical Mw for embrittlement (kg/mol)
    Mc: float = 5.0              # Entanglement Mw (kg/mol)


class ChemoMechanicalModel:
    """
    Chemo-mechanical coupling model for PE degradation.

    Tracks:
    - Chain scission density from oxidation products
    - Molecular weight evolution
    - Mechanical property degradation
    """

    def __init__(self, mech_props: MechanicalProperties = None):
        """
        Initialize chemo-mechanical model.

        Args:
            mech_props: MechanicalProperties object
        """
        self.props = mech_props or MechanicalProperties()
        self.degradation_model = create_suez_film_model()

    def calculate_chain_scission(self, result: Dict) -> np.ndarray:
        """
        Calculate chain scission density from simulation results.

        Chain scissions come from:
        1. β-scission of PO• radicals (primary source)
        2. Carbonyl formation (CO species)

        Following Colin et al.:
        s(t) ≈ k_s × ∫[CO(t)] dt

        where k_s is scission yield per carbonyl.

        Args:
            result: Simulation result dictionary

        Returns:
            Chain scission density array (mol/L) at each time
        """
        # Get CO concentration (carbonyl)
        idx = self.degradation_model.idx
        y = result['y']  # Shape: (n_species, nz, n_times)

        # Average CO across thickness
        CO = y[idx['CO'], :, :]  # (nz, n_times)
        CO_avg = np.mean(CO, axis=0)  # (n_times,)

        # Scission yield: approximately 1 scission per 2 carbonyls formed
        # (some carbonyls come from non-scission pathways)
        scission_yield = 0.5

        # Chain scission density
        s = scission_yield * CO_avg

        return s

    def calculate_mw_evolution(self, result: Dict) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate molecular weight evolution.

        For random chain scission:
        1/Mn(t) - 1/Mn(0) = s(t)

        And approximately:
        Mw(t) = Mw(0) / (1 + s(t) × Mw(0))

        Args:
            result: Simulation result dictionary

        Returns:
            Tuple of (times in months, Mw array in kg/mol)
        """
        s = self.calculate_chain_scission(result)

        # Convert Mw0 from kg/mol to mol/kg for scission calculation
        Mw0 = self.props.Mw0  # kg/mol

        # For random scission in linear polymer:
        # Mw(t) = Mw0 / (1 + s × Mw0 / ρ)
        # where ρ is density in mol/L of repeat units

        # PE repeat unit: CH2-CH2, M = 28 g/mol
        rho_PE = 950  # g/L
        c_repeat = rho_PE / 28  # mol/L of repeat units ≈ 34 mol/L

        # Relative scission: fraction of chains broken
        s_relative = s / c_repeat

        # Mw evolution
        Mw = Mw0 / (1 + s_relative * Mw0)

        # Get time array
        t_months = result['t'] / (365.25 * 24 * 3600 / 12)  # Convert seconds to months

        return t_months, Mw

    def calculate_mechanical_properties(self,
                                        Mw: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Calculate mechanical properties from Mw.

        Based on Fayolle et al. (2007) correlations for PE:
        - Elongation at break shows sharp drop at Mw_c
        - Yield stress changes moderately

        Args:
            Mw: Molecular weight array (kg/mol)

        Returns:
            Dict with sigma_y, epsilon_b, E arrays
        """
        Mw0 = self.props.Mw0
        Mw_c = self.props.Mw_crit
        Mc = self.props.Mc

        # Elongation at break: ductile-brittle transition
        # For Mw > Mw_c: ductile behavior
        # For Mw < Mw_c: brittle behavior

        # Simplified model (Fayolle et al.):
        # ε_b = ε_b0 × (1 - exp(-α(Mw - Mw_c)/Mc))  for Mw > Mw_c
        # ε_b = ε_b_crit × (Mw/Mw_c)                 for Mw < Mw_c

        epsilon_b = np.zeros_like(Mw)
        mask_ductile = Mw >= Mw_c
        mask_brittle = Mw < Mw_c

        # Ductile regime
        alpha = 0.5  # Shape parameter
        epsilon_b[mask_ductile] = self.props.epsilon_b0 * (
            1 - np.exp(-alpha * (Mw[mask_ductile] - Mw_c) / Mc)
        )

        # Brittle regime
        epsilon_b[mask_brittle] = self.props.epsilon_b_crit * (
            Mw[mask_brittle] / Mw_c
        )

        # Yield stress: moderate Mw dependence
        # σ_y = σ_y0 × (1 - β(1 - Mw/Mw0))
        beta = 0.3
        sigma_y = self.props.sigma_y0 * (1 - beta * (1 - Mw / Mw0))
        sigma_y = np.maximum(sigma_y, self.props.sigma_y_crit * 0.5)

        # Young's modulus: slight increase with crystallinity changes
        # (oxidation can increase crystallinity slightly)
        E = self.props.E0 * (1 + 0.1 * (1 - Mw / Mw0))

        return {
            'sigma_y': sigma_y,
            'epsilon_b': epsilon_b,
            'E': E
        }

    def estimate_mechanical_lifetime(self,
                                     result: Dict,
                                     criterion: str = 'epsilon_b') -> float:
        """
        Estimate mechanical lifetime based on failure criterion.

        Args:
            result: Simulation result
            criterion: 'epsilon_b' or 'sigma_y'

        Returns:
            Lifetime in months
        """
        t_months, Mw = self.calculate_mw_evolution(result)
        mech = self.calculate_mechanical_properties(Mw)

        if criterion == 'epsilon_b':
            # Ductile-brittle transition
            threshold = self.props.epsilon_b_crit
            prop = mech['epsilon_b']
        elif criterion == 'sigma_y':
            threshold = self.props.sigma_y_crit
            prop = mech['sigma_y']
        else:
            raise ValueError(f"Unknown criterion: {criterion}")

        # Find threshold crossing
        idx = np.argmax(prop < threshold)
        if idx == 0 and prop[0] >= threshold:
            return t_months[-1]  # Never reaches threshold

        # Linear interpolation
        if idx > 0:
            t1, t2 = t_months[idx-1], t_months[idx]
            p1, p2 = prop[idx-1], prop[idx]
            if p1 != p2:
                t_life = t1 + (t2 - t1) * (p1 - threshold) / (p1 - p2)
            else:
                t_life = t1
        else:
            t_life = 0.0

        return t_life

    def run_coupled_simulation(self,
                               T_celsius: float = 40.0,
                               DOC_ppm: float = 0.05,
                               t_end_years: float = 2.0,
                               verbose: bool = True) -> Dict:
        """
        Run coupled chemo-mechanical simulation.

        Args:
            T_celsius: Temperature
            DOC_ppm: DOC concentration
            t_end_years: Simulation duration
            verbose: Print progress

        Returns:
            Dictionary with all results
        """
        if verbose:
            print("="*70)
            print("CHEMO-MECHANICAL COUPLED SIMULATION")
            print("="*70)
            print(f"Temperature: {T_celsius}°C")
            print(f"DOC: {DOC_ppm} ppm")
            print(f"Duration: {t_end_years} years")
            print()

        # Run degradation simulation
        params = SimulationParams(
            T_celsius=T_celsius,
            DOC_ppm=DOC_ppm,
            t_end_years=t_end_years,
            n_timepoints=100
        )

        result = self.degradation_model.simulate(params, verbose=verbose)

        if not result['success']:
            return {'success': False}

        # Calculate OIT profile
        t_months, oit_avg = self.degradation_model.calculate_average_oit(result)

        # Calculate chain scission and Mw
        s = self.calculate_chain_scission(result)
        t_m, Mw = self.calculate_mw_evolution(result)

        # Calculate mechanical properties
        mech = self.calculate_mechanical_properties(Mw)

        # Estimate lifetimes
        oit_lifetime = None
        for i, oit in enumerate(oit_avg):
            if oit < 10.0:  # OIT threshold
                oit_lifetime = t_months[i]
                break

        mech_lifetime = self.estimate_mechanical_lifetime(result, 'epsilon_b')

        results = {
            'success': True,
            't_months': t_months,
            'oit_avg': oit_avg,
            'chain_scission': s,
            'Mw': Mw,
            'sigma_y': mech['sigma_y'],
            'epsilon_b': mech['epsilon_b'],
            'E': mech['E'],
            'oit_lifetime': oit_lifetime,
            'mech_lifetime': mech_lifetime,
            'params': params
        }

        if verbose:
            print()
            print("RESULTS:")
            print(f"  OIT lifetime (AH depletion): {oit_lifetime:.1f} months" if oit_lifetime else "  OIT lifetime: > simulation time")
            print(f"  Mechanical lifetime (ε_b < 50%): {mech_lifetime:.1f} months" if mech_lifetime else "  Mechanical lifetime: > simulation time")
            print(f"  Final Mw: {Mw[-1]:.1f} kg/mol (initial: {self.props.Mw0} kg/mol)")
            print(f"  Final ε_b: {mech['epsilon_b'][-1]:.1f}% (initial: {self.props.epsilon_b0}%)")

        return results

    def plot_coupled_results(self,
                             results: Dict,
                             save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot coupled simulation results.

        Args:
            results: Results from run_coupled_simulation
            save_path: Path to save figure

        Returns:
            matplotlib Figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        t = results['t_months']

        # Panel A: OIT evolution
        ax = axes[0, 0]
        ax.plot(t, results['oit_avg'], 'b-', linewidth=2)
        ax.axhline(y=10, color='r', linestyle='--', label='OIT threshold (10 min)')
        ax.set_xlabel('Time (months)', fontsize=12)
        ax.set_ylabel('Average OIT (min)', fontsize=12)
        ax.set_title('A) Antioxidant Depletion (OIT)', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, t[-1])

        # Panel B: Molecular weight
        ax = axes[0, 1]
        ax.plot(t, results['Mw'], 'g-', linewidth=2)
        ax.axhline(y=self.props.Mw_crit, color='r', linestyle='--',
                   label=f'Critical Mw ({self.props.Mw_crit} kg/mol)')
        ax.set_xlabel('Time (months)', fontsize=12)
        ax.set_ylabel('Mw (kg/mol)', fontsize=12)
        ax.set_title('B) Molecular Weight Evolution', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, t[-1])

        # Panel C: Elongation at break
        ax = axes[1, 0]
        ax.plot(t, results['epsilon_b'], 'r-', linewidth=2)
        ax.axhline(y=self.props.epsilon_b_crit, color='k', linestyle='--',
                   label=f'Ductile-Brittle ({self.props.epsilon_b_crit}%)')
        ax.set_xlabel('Time (months)', fontsize=12)
        ax.set_ylabel('Elongation at Break (%)', fontsize=12)
        ax.set_title('C) Mechanical Property Evolution', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, t[-1])

        # Panel D: Yield stress
        ax = axes[1, 1]
        ax.plot(t, results['sigma_y'], 'm-', linewidth=2, label='Yield Stress')
        ax.axhline(y=self.props.sigma_y_crit, color='r', linestyle='--',
                   label=f'Critical σ_y ({self.props.sigma_y_crit} MPa)')
        ax.set_xlabel('Time (months)', fontsize=12)
        ax.set_ylabel('Yield Stress (MPa)', fontsize=12)
        ax.set_title('D) Yield Stress Evolution', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, t[-1])

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")

        return fig


def run_stage3_analysis():
    """Run complete Stage 3 coupling analysis."""
    print("="*70)
    print("STAGE 3: CHEMO-MECHANICAL COUPLING ANALYSIS")
    print("="*70)
    print()

    cm = ChemoMechanicalModel()

    # Run coupled simulation at SUEZ conditions
    results = cm.run_coupled_simulation(
        T_celsius=40.0,
        DOC_ppm=0.05,
        t_end_years=2.0
    )

    if results['success']:
        # Generate plot
        fig = cm.plot_coupled_results(
            results,
            save_path='/home/user/Stage_2025/figures/chemo_mechanical_coupling.png'
        )
        plt.savefig('/home/user/Stage_2025/figures/chemo_mechanical_coupling.pdf',
                    dpi=300, bbox_inches='tight')
        plt.close()

        # Compare OIT vs mechanical lifetime
        print()
        print("="*70)
        print("LIFETIME COMPARISON")
        print("="*70)
        print()
        print("Failure Criteria Comparison:")
        print(f"  OIT-based (AH depletion < 10 min): {results['oit_lifetime']:.1f} months")
        print(f"  Mechanical (ε_b < 50%):            {results['mech_lifetime']:.1f} months")
        print()

        if results['oit_lifetime'] and results['mech_lifetime']:
            ratio = results['mech_lifetime'] / results['oit_lifetime']
            print(f"  Ratio (mechanical/OIT): {ratio:.2f}")
            print()
            if ratio > 1:
                print("  → OIT criterion is MORE CONSERVATIVE (predicts shorter lifetime)")
            else:
                print("  → Mechanical criterion is MORE CONSERVATIVE")

    return results


# =============================================================================
# Main execution
# =============================================================================

if __name__ == '__main__':
    results = run_stage3_analysis()

    print()
    print("="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
