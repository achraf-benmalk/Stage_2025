"""
Degraded Material Properties Module

Maps molecular weight evolution to mechanical property degradation
using Fayolle et al. (2007) correlations for polyethylene.

Key relationships:
- E(Mw): Elastic modulus depends on entanglement density
- σ_y(Mw): Yield stress decreases with chain scission
- ε_b(Mw): Elongation at break with ductile-brittle transition
- Fracture toughness degradation

Author: Stage 2025 Research Project
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, Callable
import matplotlib.pyplot as plt


@dataclass
class MaterialConstants:
    """Physical constants for PE100 material."""
    # Molecular weights (g/mol)
    Mw0: float = 150_000           # Initial weight-average Mw
    Mn0: float = 15_000            # Initial number-average Mn
    Mc: float = 3_800              # Critical entanglement Mw (PE)
    Mw_critical: float = 40_000   # Ductile-brittle transition Mw

    # Virgin mechanical properties
    E0: float = 1000e6             # Young's modulus (Pa) = 1 GPa
    sigma_y0: float = 25e6         # Yield stress (Pa) = 25 MPa
    epsilon_b0: float = 600.0      # Elongation at break (%)
    nu: float = 0.42               # Poisson's ratio (nearly incompressible)

    # Density
    rho: float = 950.0             # kg/m³

    # Failure thresholds
    sigma_y_critical: float = 15e6  # Critical yield stress (Pa)
    epsilon_b_critical: float = 50.0  # Brittle threshold (%)


class DegradedMaterial:
    """
    Maps molecular weight to mechanical properties.

    Based on Fayolle et al. (2007) correlations for PE:
    - E ~ (Mw/Mc)^0.5 for entangled polymers
    - σ_y has weak Mw dependence
    - ε_b shows sharp ductile-brittle transition

    Attributes:
        constants: MaterialConstants dataclass
    """

    def __init__(self, constants: MaterialConstants = None):
        """
        Initialize degraded material model.

        Args:
            constants: Material constants (uses PE100 defaults if None)
        """
        self.constants = constants or MaterialConstants()
        self._validate_constants()

    def _validate_constants(self):
        """Validate material constants are physically reasonable."""
        c = self.constants
        assert c.Mw0 > c.Mw_critical, "Initial Mw must exceed critical Mw"
        assert c.Mw_critical > c.Mc, "Critical Mw must exceed entanglement Mw"
        assert c.E0 > 0, "Modulus must be positive"
        assert 0 < c.nu < 0.5, "Poisson ratio must be in (0, 0.5)"

    def elastic_modulus(self, Mw: np.ndarray) -> np.ndarray:
        """
        Calculate Young's modulus from molecular weight.

        For entangled polymers:
        E = E0 * (Mw/Mw0)^α where α ≈ 0.5

        Below entanglement threshold, modulus drops sharply.

        Args:
            Mw: Molecular weight (g/mol), scalar or array

        Returns:
            Young's modulus (Pa)
        """
        Mw = np.atleast_1d(Mw)
        c = self.constants

        # Entanglement factor
        entanglement_ratio = np.maximum(Mw / c.Mc, 1.0)
        entanglement_ratio_0 = c.Mw0 / c.Mc

        # E scales with entanglement density
        # E ~ ρ_e ~ Mw/Mc for Mw > Mc
        alpha = 0.5
        E = c.E0 * (entanglement_ratio / entanglement_ratio_0) ** alpha

        # Below critical Mw, additional degradation
        mask_degraded = Mw < c.Mw_critical
        if np.any(mask_degraded):
            degradation_factor = (Mw[mask_degraded] / c.Mw_critical) ** 0.3
            E[mask_degraded] *= degradation_factor

        return np.squeeze(E)

    def yield_stress(self, Mw: np.ndarray) -> np.ndarray:
        """
        Calculate yield stress from molecular weight.

        Yield stress has relatively weak Mw dependence for ductile regime,
        but drops significantly in brittle regime.

        σ_y = σ_y0 * (1 - β*(1 - Mw/Mw0))  for Mw > Mw_critical
        σ_y = σ_y_crit * (Mw/Mw_critical)^0.5  for Mw < Mw_critical

        Args:
            Mw: Molecular weight (g/mol)

        Returns:
            Yield stress (Pa)
        """
        Mw = np.atleast_1d(Mw)
        c = self.constants

        sigma_y = np.zeros_like(Mw, dtype=float)

        # Ductile regime
        mask_ductile = Mw >= c.Mw_critical
        beta = 0.25  # Sensitivity coefficient
        sigma_y[mask_ductile] = c.sigma_y0 * (
            1 - beta * (1 - Mw[mask_ductile] / c.Mw0)
        )

        # Brittle regime - accelerated degradation
        mask_brittle = Mw < c.Mw_critical
        sigma_y[mask_brittle] = c.sigma_y_critical * (
            Mw[mask_brittle] / c.Mw_critical
        ) ** 0.5

        return np.squeeze(sigma_y)

    def elongation_at_break(self, Mw: np.ndarray) -> np.ndarray:
        """
        Calculate elongation at break from molecular weight.

        Shows characteristic ductile-brittle transition:
        - Ductile (Mw > Mw_c): High elongation (100-600%)
        - Brittle (Mw < Mw_c): Low elongation (<50%)

        Args:
            Mw: Molecular weight (g/mol)

        Returns:
            Elongation at break (%)
        """
        Mw = np.atleast_1d(Mw)
        c = self.constants

        epsilon_b = np.zeros_like(Mw, dtype=float)

        # Ductile regime - gradual decrease
        mask_ductile = Mw >= c.Mw_critical
        # Asymptotic approach to epsilon_b0
        gamma = 0.3
        epsilon_b[mask_ductile] = c.epsilon_b0 * (
            1 - np.exp(-gamma * (Mw[mask_ductile] - c.Mw_critical) / c.Mc)
        )

        # Brittle regime - drops to critical
        mask_brittle = Mw < c.Mw_critical
        epsilon_b[mask_brittle] = c.epsilon_b_critical * (
            Mw[mask_brittle] / c.Mw_critical
        )

        return np.squeeze(epsilon_b)

    def poisson_ratio(self, Mw: np.ndarray) -> np.ndarray:
        """
        Calculate Poisson's ratio (weak Mw dependence for PE).

        Args:
            Mw: Molecular weight (g/mol)

        Returns:
            Poisson's ratio (dimensionless)
        """
        Mw = np.atleast_1d(Mw)
        # Nearly constant for PE
        nu = np.full_like(Mw, self.constants.nu, dtype=float)
        return np.squeeze(nu)

    def shear_modulus(self, Mw: np.ndarray) -> np.ndarray:
        """
        Calculate shear modulus G = E / (2*(1+nu)).

        Args:
            Mw: Molecular weight (g/mol)

        Returns:
            Shear modulus (Pa)
        """
        E = self.elastic_modulus(Mw)
        nu = self.poisson_ratio(Mw)
        return E / (2 * (1 + nu))

    def bulk_modulus(self, Mw: np.ndarray) -> np.ndarray:
        """
        Calculate bulk modulus K = E / (3*(1-2*nu)).

        Args:
            Mw: Molecular weight (g/mol)

        Returns:
            Bulk modulus (Pa)
        """
        E = self.elastic_modulus(Mw)
        nu = self.poisson_ratio(Mw)
        return E / (3 * (1 - 2 * nu))

    def lame_parameters(self, Mw: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate Lamé parameters (λ, μ) for elasticity.

        λ = E*ν / ((1+ν)*(1-2ν))
        μ = E / (2*(1+ν)) = G

        Args:
            Mw: Molecular weight (g/mol)

        Returns:
            Tuple of (lambda, mu) in Pa
        """
        E = self.elastic_modulus(Mw)
        nu = self.poisson_ratio(Mw)

        lam = E * nu / ((1 + nu) * (1 - 2 * nu))
        mu = E / (2 * (1 + nu))

        return lam, mu

    def is_brittle(self, Mw: np.ndarray) -> np.ndarray:
        """
        Determine if material is in brittle state.

        Args:
            Mw: Molecular weight (g/mol)

        Returns:
            Boolean array (True = brittle)
        """
        Mw = np.atleast_1d(Mw)
        return Mw < self.constants.Mw_critical

    def damage_index(self, Mw: np.ndarray) -> np.ndarray:
        """
        Calculate normalized damage index (0 = virgin, 1 = fully degraded).

        D = 1 - Mw/Mw0

        Args:
            Mw: Molecular weight (g/mol)

        Returns:
            Damage index (0-1)
        """
        Mw = np.atleast_1d(Mw)
        D = 1 - Mw / self.constants.Mw0
        return np.clip(np.squeeze(D), 0, 1)

    def get_all_properties(self, Mw: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Get all mechanical properties at given Mw.

        Args:
            Mw: Molecular weight (g/mol)

        Returns:
            Dictionary with all properties
        """
        lam, mu = self.lame_parameters(Mw)

        return {
            'E': self.elastic_modulus(Mw),
            'nu': self.poisson_ratio(Mw),
            'G': self.shear_modulus(Mw),
            'K': self.bulk_modulus(Mw),
            'lambda': lam,
            'mu': mu,
            'sigma_y': self.yield_stress(Mw),
            'epsilon_b': self.elongation_at_break(Mw),
            'is_brittle': self.is_brittle(Mw),
            'damage': self.damage_index(Mw)
        }

    def plot_property_evolution(self, save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot mechanical properties vs molecular weight.

        Args:
            save_path: Path to save figure

        Returns:
            matplotlib Figure
        """
        Mw = np.linspace(10000, 160000, 200)

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        c = self.constants

        # Panel A: Elastic modulus
        ax = axes[0, 0]
        E = self.elastic_modulus(Mw) / 1e6  # Convert to MPa
        ax.plot(Mw / 1000, E, 'b-', linewidth=2)
        ax.axvline(x=c.Mw_critical / 1000, color='r', linestyle='--',
                   alpha=0.7, label='Ductile-Brittle')
        ax.axvline(x=c.Mw0 / 1000, color='g', linestyle='--',
                   alpha=0.7, label='Virgin Mw')
        ax.set_xlabel('Molecular Weight (kg/mol)', fontsize=12)
        ax.set_ylabel('Young\'s Modulus E (MPa)', fontsize=12)
        ax.set_title('A) Elastic Modulus vs Mw', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Panel B: Yield stress
        ax = axes[0, 1]
        sigma_y = self.yield_stress(Mw) / 1e6  # Convert to MPa
        ax.plot(Mw / 1000, sigma_y, 'r-', linewidth=2)
        ax.axvline(x=c.Mw_critical / 1000, color='r', linestyle='--', alpha=0.7)
        ax.axhline(y=c.sigma_y_critical / 1e6, color='k', linestyle=':',
                   alpha=0.7, label='Critical σ_y')
        ax.set_xlabel('Molecular Weight (kg/mol)', fontsize=12)
        ax.set_ylabel('Yield Stress σ_y (MPa)', fontsize=12)
        ax.set_title('B) Yield Stress vs Mw', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Panel C: Elongation at break
        ax = axes[1, 0]
        eps_b = self.elongation_at_break(Mw)
        ax.plot(Mw / 1000, eps_b, 'g-', linewidth=2)
        ax.axvline(x=c.Mw_critical / 1000, color='r', linestyle='--',
                   alpha=0.7, label='Ductile-Brittle')
        ax.axhline(y=c.epsilon_b_critical, color='k', linestyle=':',
                   alpha=0.7, label='Brittle threshold')
        ax.set_xlabel('Molecular Weight (kg/mol)', fontsize=12)
        ax.set_ylabel('Elongation at Break (%)', fontsize=12)
        ax.set_title('C) Elongation at Break vs Mw', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 650)

        # Panel D: Damage index
        ax = axes[1, 1]
        D = self.damage_index(Mw)
        is_brittle = self.is_brittle(Mw)

        ax.fill_between(Mw / 1000, 0, 1, where=is_brittle,
                        alpha=0.3, color='red', label='Brittle zone')
        ax.plot(Mw / 1000, D, 'k-', linewidth=2, label='Damage index')
        ax.set_xlabel('Molecular Weight (kg/mol)', fontsize=12)
        ax.set_ylabel('Damage Index D', fontsize=12)
        ax.set_title('D) Damage Index vs Mw', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1.1)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")

        return fig


class PropertyField:
    """
    Represents spatially-varying material properties through thickness.

    Maps z-coordinate to Mw to mechanical properties.
    """

    def __init__(self,
                 z_coords: np.ndarray,
                 Mw_profile: np.ndarray,
                 material: DegradedMaterial = None):
        """
        Initialize property field from Mw profile.

        Args:
            z_coords: Spatial coordinates (m)
            Mw_profile: Molecular weight at each z (g/mol)
            material: DegradedMaterial instance
        """
        self.z = np.atleast_1d(z_coords)
        self.Mw = np.atleast_1d(Mw_profile)
        self.material = material or DegradedMaterial()

        assert len(self.z) == len(self.Mw), "z and Mw must have same length"

        # Pre-compute properties
        self._compute_properties()

    def _compute_properties(self):
        """Compute all properties on the grid."""
        self.E = self.material.elastic_modulus(self.Mw)
        self.nu = self.material.poisson_ratio(self.Mw)
        self.sigma_y = self.material.yield_stress(self.Mw)
        self.epsilon_b = self.material.elongation_at_break(self.Mw)
        self.lam, self.mu = self.material.lame_parameters(self.Mw)
        self.is_brittle = self.material.is_brittle(self.Mw)
        self.damage = self.material.damage_index(self.Mw)

    def interpolate(self, z_new: np.ndarray) -> 'PropertyField':
        """
        Interpolate properties to new z coordinates.

        Args:
            z_new: New spatial coordinates

        Returns:
            New PropertyField at interpolated locations
        """
        Mw_new = np.interp(z_new, self.z, self.Mw)
        return PropertyField(z_new, Mw_new, self.material)

    def get_property(self, name: str) -> np.ndarray:
        """
        Get named property array.

        Args:
            name: Property name ('E', 'nu', 'sigma_y', etc.)

        Returns:
            Property array
        """
        return getattr(self, name)

    def average_property(self, name: str) -> float:
        """
        Get thickness-averaged property.

        Args:
            name: Property name

        Returns:
            Average value
        """
        prop = self.get_property(name)
        return np.trapz(prop, self.z) / (self.z[-1] - self.z[0])

    def minimum_property(self, name: str) -> Tuple[float, float]:
        """
        Get minimum property value and location.

        Args:
            name: Property name

        Returns:
            Tuple of (min_value, z_location)
        """
        prop = self.get_property(name)
        idx = np.argmin(prop)
        return prop[idx], self.z[idx]


# =============================================================================
# Main execution - demonstration
# =============================================================================

if __name__ == '__main__':
    print("="*70)
    print("DEGRADED MATERIAL PROPERTIES - DEMONSTRATION")
    print("="*70)

    # Create material model
    material = DegradedMaterial()

    # Test at various Mw values
    test_Mw = [150000, 100000, 60000, 40000, 20000]

    print("\nProperty Evolution with Degradation:")
    print(f"{'Mw (kg/mol)':<15} {'E (MPa)':<12} {'σ_y (MPa)':<12} {'ε_b (%)':<12} {'Brittle?':<10}")
    print("-"*60)

    for Mw in test_Mw:
        E = material.elastic_modulus(Mw) / 1e6
        sigma_y = material.yield_stress(Mw) / 1e6
        eps_b = material.elongation_at_break(Mw)
        brittle = material.is_brittle(Mw)

        print(f"{Mw/1000:<15.0f} {E:<12.1f} {sigma_y:<12.1f} {eps_b:<12.1f} {str(brittle):<10}")

    # Generate property evolution plot
    fig = material.plot_property_evolution(
        save_path='/home/user/Stage_2025/figures/material_property_evolution.png'
    )
    plt.savefig('/home/user/Stage_2025/figures/material_property_evolution.pdf',
                dpi=300, bbox_inches='tight')
    plt.close()

    print("\nFigure saved: material_property_evolution.png/pdf")
    print("="*70)
