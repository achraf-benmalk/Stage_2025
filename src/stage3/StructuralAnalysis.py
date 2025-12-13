"""
Structural Analysis Module

FEM-based structural analysis for PE pipe with degraded material properties.
Implements axisymmetric elasticity with spatially-varying properties.

Key features:
- 1D through-thickness FEM (axisymmetric assumption)
- Spatially-varying E(z), σ_y(z) from degradation
- Internal pressure and external load cases
- Failure analysis (von Mises vs yield)

Author: Stage 2025 Research Project
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, Tuple, Optional, List
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from DegradedMaterialProperties import DegradedMaterial, PropertyField, MaterialConstants


@dataclass
class PipeGeometry:
    """Pipe geometry parameters."""
    R_inner: float = 0.050      # Inner radius (m) = 50 mm
    thickness: float = 0.0045   # Wall thickness (m) = 4.5 mm
    length: float = 0.100       # Pipe length for analysis (m)

    @property
    def R_outer(self) -> float:
        return self.R_inner + self.thickness

    @property
    def R_mean(self) -> float:
        return self.R_inner + self.thickness / 2


@dataclass
class LoadCase:
    """Loading conditions."""
    P_internal: float = 0.0     # Internal pressure (Pa)
    P_external: float = 0.0     # External pressure (Pa)
    F_axial: float = 0.0        # Axial force (N)
    T_change: float = 0.0       # Temperature change (K)


class AxisymmetricFEM:
    """
    1D Finite Element solver for axisymmetric thick-walled cylinder.

    Solves the equilibrium equations for an axisymmetric body with
    spatially-varying material properties E(r), ν(r).

    The displacement field u(r) satisfies:
    d/dr[r σ_rr] - σ_θθ = 0

    With constitutive relations:
    σ_rr = (λ + 2μ) ε_rr + λ ε_θθ
    σ_θθ = λ ε_rr + (λ + 2μ) ε_θθ
    """

    def __init__(self,
                 geometry: PipeGeometry,
                 n_elements: int = 50):
        """
        Initialize FEM solver.

        Args:
            geometry: PipeGeometry instance
            n_elements: Number of elements through thickness
        """
        self.geometry = geometry
        self.n_elements = n_elements
        self.n_nodes = n_elements + 1

        # Create mesh
        self._create_mesh()

        # Initialize material (uniform initially)
        self.material = DegradedMaterial()
        self.property_field = None

    def _create_mesh(self):
        """Create 1D mesh through thickness."""
        g = self.geometry

        # Nodal coordinates (radial)
        self.r_nodes = np.linspace(g.R_inner, g.R_outer, self.n_nodes)

        # Element connectivity
        self.elements = np.column_stack([
            np.arange(self.n_elements),
            np.arange(1, self.n_nodes)
        ])

        # Element centers and lengths
        self.r_elem = (self.r_nodes[:-1] + self.r_nodes[1:]) / 2
        self.h_elem = np.diff(self.r_nodes)

        # Coordinate through thickness (z = 0 at inner, z = L at outer)
        self.z_nodes = self.r_nodes - g.R_inner

    def set_material_properties(self, property_field: PropertyField):
        """
        Set spatially-varying material properties.

        Args:
            property_field: PropertyField instance with E(z), ν(z), etc.
        """
        # Interpolate to element centers
        self.property_field = property_field.interpolate(self.z_nodes)

        # Get Lamé parameters at element centers
        z_elem = self.r_elem - self.geometry.R_inner
        elem_props = property_field.interpolate(z_elem)
        self.lam_elem = elem_props.lam
        self.mu_elem = elem_props.mu
        self.E_elem = elem_props.E
        self.sigma_y_elem = elem_props.sigma_y

    def set_uniform_material(self, Mw: float = 150000):
        """
        Set uniform material properties.

        Args:
            Mw: Molecular weight (g/mol)
        """
        Mw_profile = np.full(self.n_nodes, Mw)
        pf = PropertyField(self.z_nodes, Mw_profile, self.material)
        self.set_material_properties(pf)

    def _element_stiffness(self, e: int) -> np.ndarray:
        """
        Compute element stiffness matrix for axisymmetric elasticity.

        For linear element with nodes at r1, r2:
        K_e is 2x2 matrix coupling radial displacements

        Args:
            e: Element index

        Returns:
            2x2 element stiffness matrix
        """
        r1, r2 = self.r_nodes[e], self.r_nodes[e + 1]
        h = r2 - r1
        r_c = (r1 + r2) / 2  # Element center

        # Material at element center
        lam = self.lam_elem[e]
        mu = self.mu_elem[e]
        C11 = lam + 2 * mu  # σ_rr / ε_rr coefficient

        # Axisymmetric stiffness matrix (simplified 1D)
        # Based on strain-displacement: ε_rr = du/dr, ε_θθ = u/r

        # Contribution from radial strain (du/dr)
        K_rad = (C11 * r_c / h) * np.array([
            [1, -1],
            [-1, 1]
        ])

        # Contribution from hoop strain (u/r) - adds stiffness
        # Using lumped approximation
        K_hoop = ((lam + 2 * mu) * h / (2 * r_c)) * np.array([
            [1, 0],
            [0, 1]
        ])

        # Cross-coupling term
        K_cross = (lam / 2) * np.array([
            [-1, 1],
            [-1, 1]
        ])

        return K_rad + K_hoop + K_cross

    def _assemble_stiffness(self) -> sparse.csr_matrix:
        """
        Assemble global stiffness matrix.

        Returns:
            Sparse stiffness matrix
        """
        # Initialize sparse matrix
        K = sparse.lil_matrix((self.n_nodes, self.n_nodes))

        # Assembly loop
        for e in range(self.n_elements):
            K_e = self._element_stiffness(e)
            nodes = self.elements[e]

            for i in range(2):
                for j in range(2):
                    K[nodes[i], nodes[j]] += K_e[i, j]

        return K.tocsr()

    def _apply_boundary_conditions(self,
                                   K: sparse.csr_matrix,
                                   f: np.ndarray,
                                   load: LoadCase) -> Tuple[sparse.csr_matrix, np.ndarray]:
        """
        Apply boundary conditions.

        For pressure loading:
        - Traction at inner surface: σ_rr(r=Ri) = -P_internal
        - Traction at outer surface: σ_rr(r=Ro) = -P_external

        Args:
            K: Stiffness matrix
            f: Force vector
            load: LoadCase instance

        Returns:
            Modified (K, f)
        """
        # Convert to lil for modification
        K = K.tolil()

        # Pressure loads on boundaries
        # For axisymmetric: F = 2π r P (force per unit length)
        g = self.geometry

        # Inner surface pressure (inward positive = tension negative)
        if load.P_internal != 0:
            f[0] += -load.P_internal * 2 * np.pi * g.R_inner

        # Outer surface pressure
        if load.P_external != 0:
            f[-1] += load.P_external * 2 * np.pi * g.R_outer

        return K.tocsr(), f

    def solve(self, load: LoadCase) -> Dict:
        """
        Solve the elasticity problem.

        Args:
            load: LoadCase instance

        Returns:
            Dictionary with solution fields
        """
        # Ensure material is set
        if self.property_field is None:
            self.set_uniform_material()

        # Assemble stiffness
        K = self._assemble_stiffness()

        # Initialize force vector
        f = np.zeros(self.n_nodes)

        # Apply boundary conditions
        K, f = self._apply_boundary_conditions(K, f, load)

        # For well-posed problem, fix one node (outer surface)
        # or use penalty method
        penalty = 1e20
        K = K.tolil()
        K[-1, -1] += penalty
        f[-1] += 0 * penalty  # Zero displacement at outer
        K = K.tocsr()

        # Solve
        u = spsolve(K, f)

        # Post-process
        results = self._post_process(u)
        results['u'] = u
        results['load'] = load

        return results

    def _post_process(self, u: np.ndarray) -> Dict:
        """
        Compute stresses and strains from displacement.

        Args:
            u: Nodal displacement vector

        Returns:
            Dictionary with stress/strain fields
        """
        n_elem = self.n_elements

        # Initialize arrays at element centers
        eps_rr = np.zeros(n_elem)
        eps_tt = np.zeros(n_elem)
        sigma_rr = np.zeros(n_elem)
        sigma_tt = np.zeros(n_elem)

        for e in range(n_elem):
            r1, r2 = self.r_nodes[e], self.r_nodes[e + 1]
            u1, u2 = u[e], u[e + 1]
            h = r2 - r1
            r_c = self.r_elem[e]

            # Strains
            eps_rr[e] = (u2 - u1) / h  # du/dr
            eps_tt[e] = (u1 + u2) / (2 * r_c)  # u/r (averaged)

            # Stresses using Lamé
            lam = self.lam_elem[e]
            mu = self.mu_elem[e]

            sigma_rr[e] = (lam + 2 * mu) * eps_rr[e] + lam * eps_tt[e]
            sigma_tt[e] = lam * eps_rr[e] + (lam + 2 * mu) * eps_tt[e]

        # Von Mises stress (plane strain assumption: σ_zz = ν(σ_rr + σ_tt))
        # For axisymmetric: σ_vm = sqrt(σ_rr² - σ_rr*σ_tt + σ_tt²)
        sigma_vm = np.sqrt(sigma_rr**2 - sigma_rr * sigma_tt + sigma_tt**2)

        # Safety factor vs yield
        safety_factor = self.sigma_y_elem / np.maximum(sigma_vm, 1e-10)

        return {
            'r_elem': self.r_elem,
            'z_elem': self.r_elem - self.geometry.R_inner,
            'eps_rr': eps_rr,
            'eps_tt': eps_tt,
            'sigma_rr': sigma_rr,
            'sigma_tt': sigma_tt,
            'sigma_vm': sigma_vm,
            'sigma_y': self.sigma_y_elem,
            'safety_factor': safety_factor,
            'min_safety_factor': np.min(safety_factor),
            'failure_location': self.r_elem[np.argmin(safety_factor)]
        }


class PipeStructuralModel:
    """
    Complete structural model for pipe with degradation.

    Combines:
    - Kinetic model output (Mw profiles)
    - Material property correlations
    - FEM structural analysis
    - Failure prediction
    """

    def __init__(self,
                 geometry: PipeGeometry = None,
                 n_elements: int = 50):
        """
        Initialize structural model.

        Args:
            geometry: PipeGeometry (defaults to PE100 pipe)
            n_elements: FEM mesh density
        """
        self.geometry = geometry or PipeGeometry()
        self.fem = AxisymmetricFEM(self.geometry, n_elements)
        self.material = DegradedMaterial()

    def analyze_virgin(self, P_internal: float) -> Dict:
        """
        Analyze virgin (undegraded) pipe.

        Args:
            P_internal: Internal pressure (Pa)

        Returns:
            Analysis results
        """
        self.fem.set_uniform_material(Mw=150000)
        load = LoadCase(P_internal=P_internal)
        return self.fem.solve(load)

    def analyze_degraded(self,
                         Mw_profile: np.ndarray,
                         z_coords: np.ndarray,
                         P_internal: float) -> Dict:
        """
        Analyze degraded pipe with Mw gradient.

        Args:
            Mw_profile: Molecular weight through thickness
            z_coords: z coordinates (m)
            P_internal: Internal pressure (Pa)

        Returns:
            Analysis results
        """
        pf = PropertyField(z_coords, Mw_profile, self.material)
        self.fem.set_material_properties(pf)
        load = LoadCase(P_internal=P_internal)
        return self.fem.solve(load)

    def find_failure_pressure(self,
                              Mw_profile: np.ndarray = None,
                              z_coords: np.ndarray = None,
                              safety_target: float = 1.0) -> float:
        """
        Find pressure at which safety factor reaches target.

        Uses bisection to find P where min(SF) = safety_target.

        Args:
            Mw_profile: Mw through thickness (None = virgin)
            z_coords: z coordinates
            safety_target: Target safety factor

        Returns:
            Failure pressure (Pa)
        """
        # Set material
        if Mw_profile is not None and z_coords is not None:
            pf = PropertyField(z_coords, Mw_profile, self.material)
            self.fem.set_material_properties(pf)
        else:
            self.fem.set_uniform_material(Mw=150000)

        # Bisection search
        P_low, P_high = 1e4, 1e8  # 0.1 bar to 1000 bar

        for _ in range(50):
            P_mid = (P_low + P_high) / 2
            result = self.fem.solve(LoadCase(P_internal=P_mid))

            if result['min_safety_factor'] > safety_target:
                P_low = P_mid
            else:
                P_high = P_mid

            if abs(P_high - P_low) < 1e3:
                break

        return P_mid

    def time_evolution_analysis(self,
                                times_months: List[float],
                                Mw_profiles: List[np.ndarray],
                                z_coords: np.ndarray,
                                P_internal: float = 1e6) -> Dict:
        """
        Analyze structural response over degradation time.

        Args:
            times_months: List of time points
            Mw_profiles: List of Mw profiles at each time
            z_coords: z coordinates
            P_internal: Operating pressure (Pa)

        Returns:
            Dictionary with time-evolution results
        """
        results = {
            'times': times_months,
            'min_safety_factor': [],
            'failure_pressure': [],
            'max_stress_vm': [],
            'failure_location': [],
            'is_failed': []
        }

        for t, Mw in zip(times_months, Mw_profiles):
            # Set degraded properties
            pf = PropertyField(z_coords, Mw, self.material)
            self.fem.set_material_properties(pf)

            # Analyze at operating pressure
            load = LoadCase(P_internal=P_internal)
            res = self.fem.solve(load)

            results['min_safety_factor'].append(res['min_safety_factor'])
            results['max_stress_vm'].append(np.max(res['sigma_vm']))
            results['failure_location'].append(res['failure_location'])
            results['is_failed'].append(res['min_safety_factor'] < 1.0)

            # Find failure pressure
            P_fail = self.find_failure_pressure(Mw, z_coords)
            results['failure_pressure'].append(P_fail)

        # Convert to arrays
        for key in results:
            if key != 'times':
                results[key] = np.array(results[key])

        return results

    def plot_stress_comparison(self,
                               results_virgin: Dict,
                               results_degraded: Dict,
                               save_path: Optional[str] = None) -> plt.Figure:
        """
        Compare stress distributions: virgin vs degraded.

        Args:
            results_virgin: FEM results for virgin pipe
            results_degraded: FEM results for degraded pipe
            save_path: Path to save figure

        Returns:
            matplotlib Figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        z_v = results_virgin['z_elem'] * 1000  # Convert to mm
        z_d = results_degraded['z_elem'] * 1000

        # Panel A: Radial stress
        ax = axes[0, 0]
        ax.plot(z_v, results_virgin['sigma_rr'] / 1e6, 'b-',
                linewidth=2, label='Virgin')
        ax.plot(z_d, results_degraded['sigma_rr'] / 1e6, 'r--',
                linewidth=2, label='Degraded')
        ax.set_xlabel('Distance from inner surface (mm)', fontsize=12)
        ax.set_ylabel('Radial Stress σ_rr (MPa)', fontsize=12)
        ax.set_title('A) Radial Stress Distribution', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Panel B: Hoop stress
        ax = axes[0, 1]
        ax.plot(z_v, results_virgin['sigma_tt'] / 1e6, 'b-',
                linewidth=2, label='Virgin')
        ax.plot(z_d, results_degraded['sigma_tt'] / 1e6, 'r--',
                linewidth=2, label='Degraded')
        ax.set_xlabel('Distance from inner surface (mm)', fontsize=12)
        ax.set_ylabel('Hoop Stress σ_θθ (MPa)', fontsize=12)
        ax.set_title('B) Hoop Stress Distribution', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Panel C: Von Mises stress vs yield
        ax = axes[1, 0]
        ax.plot(z_v, results_virgin['sigma_vm'] / 1e6, 'b-',
                linewidth=2, label='Virgin σ_vm')
        ax.plot(z_d, results_degraded['sigma_vm'] / 1e6, 'r-',
                linewidth=2, label='Degraded σ_vm')
        ax.plot(z_d, results_degraded['sigma_y'] / 1e6, 'r--',
                linewidth=2, label='Degraded σ_y')
        ax.axhline(y=25, color='b', linestyle=':', alpha=0.7, label='Virgin σ_y')

        ax.set_xlabel('Distance from inner surface (mm)', fontsize=12)
        ax.set_ylabel('Stress (MPa)', fontsize=12)
        ax.set_title('C) Von Mises Stress vs Yield Strength', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Panel D: Safety factor
        ax = axes[1, 1]
        ax.plot(z_v, results_virgin['safety_factor'], 'b-',
                linewidth=2, label='Virgin')
        ax.plot(z_d, results_degraded['safety_factor'], 'r-',
                linewidth=2, label='Degraded')
        ax.axhline(y=1.0, color='k', linestyle='--', label='Failure (SF=1)')

        ax.set_xlabel('Distance from inner surface (mm)', fontsize=12)
        ax.set_ylabel('Safety Factor', fontsize=12)
        ax.set_title('D) Safety Factor Distribution', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, max(results_virgin['safety_factor'].max(),
                          results_degraded['safety_factor'].max()) * 1.1)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")

        return fig


# =============================================================================
# Main execution - demonstration
# =============================================================================

if __name__ == '__main__':
    print("="*70)
    print("STRUCTURAL ANALYSIS - DEMONSTRATION")
    print("="*70)

    # Create model
    geometry = PipeGeometry(R_inner=0.050, thickness=0.0045)
    model = PipeStructuralModel(geometry, n_elements=50)

    print(f"\nPipe Geometry:")
    print(f"  Inner radius: {geometry.R_inner*1000:.1f} mm")
    print(f"  Thickness: {geometry.thickness*1000:.2f} mm")
    print(f"  Outer radius: {geometry.R_outer*1000:.1f} mm")

    # Operating pressure: 10 bar = 1 MPa
    P_op = 1.0e6  # Pa

    print(f"\nOperating pressure: {P_op/1e6:.1f} MPa")

    # Analyze virgin pipe
    print("\n--- VIRGIN PIPE ANALYSIS ---")
    results_virgin = model.analyze_virgin(P_op)
    print(f"  Max von Mises stress: {results_virgin['sigma_vm'].max()/1e6:.2f} MPa")
    print(f"  Min safety factor: {results_virgin['min_safety_factor']:.2f}")

    P_fail_virgin = model.find_failure_pressure()
    print(f"  Failure pressure: {P_fail_virgin/1e6:.1f} MPa")

    # Create degraded Mw profile (simulating 6 months exposure)
    z = np.linspace(0, geometry.thickness, 51)
    # Exponential decay from inner surface
    Mw_degraded = 150000 * (0.3 + 0.7 * (1 - np.exp(-z / 0.001)))

    print("\n--- DEGRADED PIPE ANALYSIS (6 months) ---")
    results_degraded = model.analyze_degraded(Mw_degraded, z, P_op)
    print(f"  Mw at inner surface: {Mw_degraded[0]/1000:.0f} kg/mol")
    print(f"  Mw at outer surface: {Mw_degraded[-1]/1000:.0f} kg/mol")
    print(f"  Max von Mises stress: {results_degraded['sigma_vm'].max()/1e6:.2f} MPa")
    print(f"  Min safety factor: {results_degraded['min_safety_factor']:.2f}")

    P_fail_degraded = model.find_failure_pressure(Mw_degraded, z)
    print(f"  Failure pressure: {P_fail_degraded/1e6:.1f} MPa")

    # Strength reduction
    reduction = 100 * (1 - P_fail_degraded / P_fail_virgin)
    print(f"\n  STRENGTH REDUCTION: {reduction:.1f}%")

    # Generate comparison plot
    fig = model.plot_stress_comparison(
        results_virgin, results_degraded,
        save_path='/home/user/Stage_2025/figures/structural_stress_comparison.png'
    )
    plt.savefig('/home/user/Stage_2025/figures/structural_stress_comparison.pdf',
                dpi=300, bbox_inches='tight')
    plt.close()

    print("\nFigure saved: structural_stress_comparison.png/pdf")
    print("="*70)
