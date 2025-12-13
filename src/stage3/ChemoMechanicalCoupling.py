"""
Chemo-Mechanical Coupling Module

Complete multi-physics coupling chain:
Chemistry (Kinetic Model) → Materials (Mw → Properties) → Structure (FEM)

Integrates:
- PE degradation kinetic model (Colin et al. 2009)
- Fayolle material property correlations
- Axisymmetric FEM structural analysis

Author: Stage 2025 Research Project
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, List
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
sys.path.insert(0, os.path.dirname(__file__))

from DegradedMaterialProperties import DegradedMaterial, PropertyField, MaterialConstants
from StructuralAnalysis import PipeGeometry, PipeStructuralModel, LoadCase


@dataclass
class CouplingConfig:
    """Configuration for chemo-mechanical coupling."""
    # Kinetic model parameters
    T_celsius: float = 40.0
    DOC_ppm: float = 0.05
    thickness_m: float = 0.0045

    # Time points for analysis
    times_months: List[float] = None

    # Structural parameters
    P_operating: float = 1.0e6  # Pa (10 bar)

    def __post_init__(self):
        if self.times_months is None:
            self.times_months = [0, 1, 2, 3, 6, 9, 12]


class MultiPhysicsCoupling:
    """
    Complete multi-physics coupling for PE pipe degradation.

    Flow:
    1. Run kinetic model → Mw(z, t) profiles
    2. Apply Fayolle correlations → E(z,t), σ_y(z,t)
    3. Run FEM at each time → stress, safety factor
    4. Predict failure time and mode
    """

    def __init__(self, config: CouplingConfig = None):
        """
        Initialize coupling framework.

        Args:
            config: CouplingConfig instance
        """
        self.config = config or CouplingConfig()

        # Initialize components
        self.geometry = PipeGeometry(thickness=self.config.thickness_m)
        self.material = DegradedMaterial()
        self.structural = PipeStructuralModel(self.geometry)

        # Results storage
        self.results = None

    def generate_mw_profiles_from_kinetics(self) -> Tuple[np.ndarray, List[np.ndarray]]:
        """
        Generate Mw profiles from kinetic model simulation.

        Uses simplified analytical approximation or loads from saved data.

        Returns:
            Tuple of (z_coords, list of Mw profiles at each time)
        """
        L = self.config.thickness_m
        nz = 51
        z = np.linspace(0, L, nz)

        # Try to load from kinetic model if available
        try:
            from DegradationModel import DegradationModel, SimulationParams, create_suez_film_model
            from ChemoMechanicalCoupling import ChemoMechanicalModel

            print("Loading Mw profiles from kinetic model...")

            cm = ChemoMechanicalModel()
            Mw_profiles = []

            for t_months in self.config.times_months:
                t_years = t_months / 12.0

                if t_years > 0:
                    params = SimulationParams(
                        T_celsius=self.config.T_celsius,
                        DOC_ppm=self.config.DOC_ppm,
                        t_end_years=t_years + 0.01,
                        n_timepoints=50
                    )
                    result = cm.degradation_model.simulate(params, verbose=False)

                    if result['success']:
                        # Get chain scission and Mw
                        t_m, Mw = cm.calculate_mw_evolution(result)
                        # Interpolate to target time
                        Mw_at_t = Mw[-1] * 1000  # Convert kg/mol to g/mol

                        # Create profile (simplified - uniform Mw assumption)
                        # In reality, would extract spatial profile
                        Mw_profile = np.full(nz, Mw_at_t)
                    else:
                        Mw_profile = np.full(nz, 150000)
                else:
                    Mw_profile = np.full(nz, 150000)

                Mw_profiles.append(Mw_profile)

            return z, Mw_profiles

        except ImportError:
            pass

        # Fallback: Use analytical approximation
        print("Using analytical Mw profile approximation...")
        return self._analytical_mw_profiles(z)

    def _analytical_mw_profiles(self, z: np.ndarray) -> Tuple[np.ndarray, List[np.ndarray]]:
        """
        Generate Mw profiles using analytical approximation.

        Based on diffusion-reaction with first-order kinetics:
        Mw(z,t) = Mw0 * exp(-k_eff * t * exp(-z/δ))

        where δ is the penetration depth.

        Args:
            z: Spatial coordinates

        Returns:
            Tuple of (z, list of Mw profiles)
        """
        Mw0 = 150000  # g/mol

        # Effective degradation rate (fitted to kinetic model)
        k_eff = 0.15  # per month at surface

        # Penetration depth (DOC diffusion length)
        delta = 0.0015  # m (1.5 mm)

        Mw_profiles = []

        for t in self.config.times_months:
            # Degradation factor decays from inner surface
            degradation = np.exp(-k_eff * t * np.exp(-z / delta))

            # Mw profile
            Mw = Mw0 * degradation

            # Apply minimum Mw (can't go below ~10 kg/mol)
            Mw = np.maximum(Mw, 10000)

            Mw_profiles.append(Mw)

        return z, Mw_profiles

    def run_coupled_analysis(self, verbose: bool = True) -> Dict:
        """
        Run complete multi-physics coupling analysis.

        Args:
            verbose: Print progress

        Returns:
            Dictionary with all results
        """
        if verbose:
            print("="*70)
            print("MULTI-PHYSICS CHEMO-MECHANICAL COUPLING")
            print("="*70)
            print(f"Temperature: {self.config.T_celsius}°C")
            print(f"DOC: {self.config.DOC_ppm} ppm")
            print(f"Pressure: {self.config.P_operating/1e6:.1f} MPa")
            print(f"Times: {self.config.times_months} months")
            print()

        # Step 1: Generate Mw profiles from chemistry
        if verbose:
            print("Step 1: Generating Mw profiles from kinetic model...")
        z, Mw_profiles = self.generate_mw_profiles_from_kinetics()

        # Step 2: Convert to material properties
        if verbose:
            print("Step 2: Computing material properties (Fayolle correlations)...")

        property_fields = []
        for Mw in Mw_profiles:
            pf = PropertyField(z, Mw, self.material)
            property_fields.append(pf)

        # Step 3: Run structural analysis at each time
        if verbose:
            print("Step 3: Running FEM structural analysis...")

        structural_results = []
        failure_pressures = []

        for i, (t, Mw) in enumerate(zip(self.config.times_months, Mw_profiles)):
            # Analyze at operating pressure
            res = self.structural.analyze_degraded(Mw, z, self.config.P_operating)
            structural_results.append(res)

            # Find failure pressure
            P_fail = self.structural.find_failure_pressure(Mw, z)
            failure_pressures.append(P_fail)

            if verbose:
                print(f"  t={t:3.0f} months: SF_min={res['min_safety_factor']:.2f}, "
                      f"P_fail={P_fail/1e6:.2f} MPa")

        # Compile results
        self.results = {
            'config': self.config,
            'z': z,
            'times': np.array(self.config.times_months),
            'Mw_profiles': Mw_profiles,
            'property_fields': property_fields,
            'structural_results': structural_results,
            'failure_pressures': np.array(failure_pressures),
            'safety_factors': np.array([r['min_safety_factor'] for r in structural_results]),
            'max_stress_vm': np.array([r['sigma_vm'].max() for r in structural_results])
        }

        # Determine failure time
        sf_array = self.results['safety_factors']
        failed_mask = sf_array < 1.0

        if np.any(failed_mask):
            failure_idx = np.argmax(failed_mask)
            self.results['failure_time'] = self.config.times_months[failure_idx]
            self.results['has_failed'] = True
        else:
            self.results['failure_time'] = None
            self.results['has_failed'] = False

        if verbose:
            print()
            print("="*70)
            print("ANALYSIS COMPLETE")
            print("="*70)
            print(f"\nVirgn failure pressure: {failure_pressures[0]/1e6:.2f} MPa")
            print(f"Final failure pressure: {failure_pressures[-1]/1e6:.2f} MPa")
            reduction = 100 * (1 - failure_pressures[-1] / failure_pressures[0])
            print(f"Strength reduction: {reduction:.1f}%")

            if self.results['has_failed']:
                print(f"\nWARNING: Pipe fails at t = {self.results['failure_time']} months!")
            else:
                print(f"\nPipe survives {self.config.times_months[-1]} months (SF > 1)")

        return self.results

    def plot_coupling_results(self, save_path: Optional[str] = None) -> plt.Figure:
        """
        Create comprehensive multi-panel figure showing coupling results.

        Args:
            save_path: Path to save figure

        Returns:
            matplotlib Figure
        """
        if self.results is None:
            raise ValueError("Must run analysis first")

        fig, axes = plt.subplots(2, 3, figsize=(15, 10))

        times = self.results['times']
        z_mm = self.results['z'] * 1000

        # Color map for time evolution
        colors = plt.cm.viridis(np.linspace(0, 1, len(times)))

        # Panel A: Mw profiles through thickness
        ax = axes[0, 0]
        for i, (t, Mw) in enumerate(zip(times, self.results['Mw_profiles'])):
            ax.plot(z_mm, Mw / 1000, color=colors[i], linewidth=2,
                    label=f't = {t:.0f} mo')
        ax.axhline(y=40, color='r', linestyle='--', alpha=0.7, label='Mw_critical')
        ax.set_xlabel('Distance from inner surface (mm)', fontsize=11)
        ax.set_ylabel('Molecular Weight (kg/mol)', fontsize=11)
        ax.set_title('A) Mw Evolution (Chemistry)', fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='lower right')
        ax.grid(True, alpha=0.3)

        # Panel B: E profiles
        ax = axes[0, 1]
        for i, (t, pf) in enumerate(zip(times, self.results['property_fields'])):
            ax.plot(z_mm, pf.E / 1e6, color=colors[i], linewidth=2,
                    label=f't = {t:.0f} mo')
        ax.set_xlabel('Distance from inner surface (mm)', fontsize=11)
        ax.set_ylabel('Young\'s Modulus (MPa)', fontsize=11)
        ax.set_title('B) Modulus Evolution (Materials)', fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='lower right')
        ax.grid(True, alpha=0.3)

        # Panel C: σ_y profiles
        ax = axes[0, 2]
        for i, (t, pf) in enumerate(zip(times, self.results['property_fields'])):
            ax.plot(z_mm, pf.sigma_y / 1e6, color=colors[i], linewidth=2,
                    label=f't = {t:.0f} mo')
        ax.axhline(y=15, color='r', linestyle='--', alpha=0.7, label='σ_y_critical')
        ax.set_xlabel('Distance from inner surface (mm)', fontsize=11)
        ax.set_ylabel('Yield Stress (MPa)', fontsize=11)
        ax.set_title('C) Yield Stress Evolution (Materials)', fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='lower right')
        ax.grid(True, alpha=0.3)

        # Panel D: Von Mises stress profiles
        ax = axes[1, 0]
        for i, (t, res) in enumerate(zip(times, self.results['structural_results'])):
            z_elem = res['z_elem'] * 1000
            ax.plot(z_elem, res['sigma_vm'] / 1e6, color=colors[i], linewidth=2,
                    label=f't = {t:.0f} mo')
        ax.set_xlabel('Distance from inner surface (mm)', fontsize=11)
        ax.set_ylabel('Von Mises Stress (MPa)', fontsize=11)
        ax.set_title('D) Stress Distribution (Structure)', fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.3)

        # Panel E: Safety factor evolution
        ax = axes[1, 1]
        ax.plot(times, self.results['safety_factors'], 'bo-', linewidth=2, markersize=8)
        ax.axhline(y=1.0, color='r', linestyle='--', linewidth=2, label='Failure (SF=1)')
        ax.fill_between(times, 0, 1, alpha=0.2, color='red', label='Failure zone')
        ax.set_xlabel('Time (months)', fontsize=11)
        ax.set_ylabel('Minimum Safety Factor', fontsize=11)
        ax.set_title('E) Safety Factor Evolution', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, max(times))
        ax.set_ylim(0, max(self.results['safety_factors']) * 1.1)

        # Panel F: Failure pressure evolution
        ax = axes[1, 2]
        P0 = self.results['failure_pressures'][0]
        P_norm = self.results['failure_pressures'] / P0 * 100

        ax.plot(times, P_norm, 'go-', linewidth=2, markersize=8)
        ax.axhline(y=100, color='b', linestyle=':', alpha=0.7, label='Virgin')
        ax.set_xlabel('Time (months)', fontsize=11)
        ax.set_ylabel('Failure Pressure (% of virgin)', fontsize=11)
        ax.set_title('F) Strength Retention', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, max(times))
        ax.set_ylim(0, 110)

        # Add annotation
        reduction = 100 - P_norm[-1]
        ax.annotate(f'{reduction:.1f}% loss',
                    xy=(times[-1], P_norm[-1]),
                    xytext=(times[-1] * 0.6, P_norm[-1] + 10),
                    fontsize=11, color='red',
                    arrowprops=dict(arrowstyle='->', color='red'))

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")

        return fig

    def generate_summary_table(self) -> str:
        """
        Generate markdown summary table.

        Returns:
            Markdown-formatted table
        """
        if self.results is None:
            return "No results available"

        lines = [
            "## Multi-Physics Coupling Results",
            "",
            f"**Conditions:** T = {self.config.T_celsius}°C, DOC = {self.config.DOC_ppm} ppm, P = {self.config.P_operating/1e6:.1f} MPa",
            "",
            "| Time (months) | Mw_inner (kg/mol) | E_inner (MPa) | σ_y_inner (MPa) | SF_min | P_fail (MPa) |",
            "|---------------|-------------------|---------------|-----------------|--------|--------------|"
        ]

        for i, t in enumerate(self.results['times']):
            Mw = self.results['Mw_profiles'][i][0] / 1000
            pf = self.results['property_fields'][i]
            E = pf.E[0] / 1e6
            sigma_y = pf.sigma_y[0] / 1e6
            SF = self.results['safety_factors'][i]
            P_fail = self.results['failure_pressures'][i] / 1e6

            lines.append(f"| {t:13.0f} | {Mw:17.1f} | {E:13.1f} | {sigma_y:15.1f} | {SF:6.2f} | {P_fail:12.2f} |")

        lines.extend([
            "",
            f"**Strength Reduction:** {100*(1 - self.results['failure_pressures'][-1]/self.results['failure_pressures'][0]):.1f}%",
            f"**Failure Status:** {'FAILED' if self.results['has_failed'] else 'SAFE'}"
        ])

        return "\n".join(lines)


# =============================================================================
# Main execution
# =============================================================================

if __name__ == '__main__':
    # Run coupling analysis
    config = CouplingConfig(
        T_celsius=40.0,
        DOC_ppm=0.05,
        times_months=[0, 1, 2, 3, 6, 9, 12],
        P_operating=1.0e6
    )

    coupling = MultiPhysicsCoupling(config)
    results = coupling.run_coupled_analysis(verbose=True)

    # Generate figure
    fig = coupling.plot_coupling_results(
        save_path='/home/user/Stage_2025/figures/chemo_mechanical_fem_coupling.png'
    )
    plt.savefig('/home/user/Stage_2025/figures/chemo_mechanical_fem_coupling.pdf',
                dpi=300, bbox_inches='tight')
    plt.close()

    # Print summary table
    print("\n")
    print(coupling.generate_summary_table())
