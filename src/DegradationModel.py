"""
PE Pipe Degradation Model - Implementation of Colin et al. (2009)

This module implements a 1D diffusion-reaction kinetic model for polyethylene (PE)
pipe degradation under chlorine dioxide (DOC) exposure, following the mechanistic
scheme of Colin et al. (2009).

The model predicts:
- Antioxidant concentration profiles [AH](z,t)
- Oxidation Induction Time (OIT) profiles
- Carbonyl buildup [CO](z,t)
- Molecular weight evolution Mw(z,t)
- Chain scissions and crosslinks

Key References:
[1] Colin, X., et al. (2009). Polymer Engineering & Science, 49(7), 1429-1437. (Part I)
[2] Colin, X., et al. (2009). Polymer Engineering & Science, 49(8), 1642-1652. (Part II)

Author: Stage 2025 Research Project
Supervisor: Xavier Colin (PIMM Laboratory)
Industrial Partner: SUEZ
"""

import numpy as np
from scipy.integrate import solve_ivp
import time
from dataclasses import dataclass
from typing import Dict, Optional, Tuple, List


@dataclass
class MaterialProperties:
    """Material properties for PE100 pipe/film."""
    Mw0: float = 150.0          # Initial molecular weight (kg/mol)
    dens0: float = 0.95         # Total PE density (kg/L)
    densa: float = 0.85         # Amorphous phase density (kg/L)
    Xc: float = 0.45            # Crystallinity fraction
    ti0_oit: float = 165.0      # Initial OIT at 190°C (min) - default from Colin
    carbon_black: float = 0.02  # Carbon black content (mass fraction)


@dataclass
class SimulationParams:
    """Parameters for a specific simulation run."""
    T_celsius: float            # Temperature (°C)
    DOC_ppm: float              # DOC/HOCl concentration in water (ppm)
    t_end_years: float          # Simulation duration (years)
    AH0_mult: float = 1.0       # Multiplier for initial AH concentration
    O2_sat_mult: float = 1.0    # Multiplier for O2 saturation
    method: str = 'Radau'       # ODE solver method
    rtol: float = 1e-6          # Relative tolerance
    atol: float = 1e-9          # Absolute tolerance
    n_timepoints: int = 100     # Number of output time points


class DegradationModel:
    """
    1D Diffusion-Reaction model for PE degradation under DOC exposure.

    Implements the Colin et al. (2009) mechanistic scheme with:
    - 12 chemical species tracked spatially
    - Diffusion of O2, DOC, and antioxidant (AH)
    - Full radical chain oxidation kinetics
    - Boundary conditions for pipe/film geometry

    Attributes:
        L (float): Sample thickness (m)
        nz (int): Number of spatial grid points
        simulation_mode (str): 'pipe' or 'film'
        material (MaterialProperties): Material properties
    """

    # Physical constants
    R_GAS = 8.314  # J/(mol·K)

    # Species names and indices
    SPECIES = ['O2', 'DOC', 'AH', 'P', 'PO2', 'POOH', 'PH', 'Q', 'CO', 'PCl', 'S', 'X']

    def __init__(
        self,
        L: float = 4.5e-3,
        nz: int = 100,
        simulation_mode: str = 'pipe',
        material: Optional[MaterialProperties] = None
    ):
        """
        Initialize the degradation model.

        Args:
            L: Sample thickness (m). Default 4.5mm for pipe, use 0.4mm for SUEZ films
            nz: Number of spatial grid points
            simulation_mode: 'pipe' (one-sided DOC exposure) or 'film' (two-sided)
            material: Material properties (uses default PE100 if None)
        """
        if simulation_mode not in ['pipe', 'film']:
            raise ValueError("simulation_mode must be 'pipe' or 'film'")

        self.L = L
        self.nz = nz
        self.z = np.linspace(0, L, nz)
        self.dz = L / (nz - 1)
        self.simulation_mode = simulation_mode

        # Material properties
        self.material = material or MaterialProperties()
        self._setup_material_derived_params()

        # Species indexing
        self.n_species = len(self.SPECIES)
        self.idx = {name: i for i, name in enumerate(self.SPECIES)}

        # Initialize kinetic parameters (will be updated for temperature)
        self._setup_kinetic_coefficients()

        # Current state
        self.current_T_K = None
        self.k = {}      # Rate constants at current T
        self.D = {}      # Diffusion coefficients at current T
        self.beta = {}   # Boundary exchange coefficients at current T

    def _setup_material_derived_params(self):
        """Calculate derived material parameters."""
        m = self.material
        self.Tam = 1.0 - m.Xc  # Amorphous fraction
        self.Tav = (m.dens0 / m.densa) * self.Tam if m.densa > 1e-9 else self.Tam

    def _setup_kinetic_coefficients(self):
        """
        Setup Arrhenius parameters for rate constants and diffusion coefficients.

        Format: (Ea [J/mol], A [pre-exponential factor])
        Rate constants: L/(mol·s) or 1/s
        Diffusion coefficients: m²/s
        """
        # Rate constants (Ea, A) from Colin et al. (2009)
        self.k_coeffs = {
            # Initiation
            'k1u': (140e3, 8.0e12),    # POOH unimolecular decomposition
            'k1b': (105e3, 2.8e9),     # POOH bimolecular decomposition
            'k1d': (0, 2.7e-5),        # DOC + PH -> P• (DOC-PE reaction)

            # Propagation
            'k2': (0, 1.0e8),          # P• + O2 -> PO2•
            'k3': (73e3, 1.5e10),      # PO2• + PH -> POOH + P•

            # Termination
            'k4': (0, 8.0e11),         # P• + P• -> inactive
            'k4d': (21.1e3, 6.6e9),    # P• + DOC -> P-Cl
            'k5': (5.9e3, 1.5e12),     # P• + PO2• -> products
            'k60': (80e3, 4.9e19),     # PO2• + PO2• -> Q + O2
            'k61': (0, 2e6),           # Q -> POOP (crosslink)
            'k62': (5e3, 1.2e6),       # Q -> CO + P-OH
            'k63': (17.4e3, 4.8e9),    # Q -> 2P• (cage escape)

            # Stabilization
            'k7': (49.9e3, 1.3e9),     # PO2• + AH -> products
            'k8d': (0, 5.0e-2),        # DOC + AH -> products
        }

        # Diffusion coefficients (Ea, D0)
        self.D_coeffs = {
            'O2': (35e3, 4.3e-5),      # Oxygen diffusion
            'DOC': (0, 2.0e-11),       # DOC diffusion (low activation)
            'AH': (115.7e3, 9.1e4),    # Antioxidant diffusion
        }

        # Boundary exchange coefficients (Ea, beta0)
        self.beta_coeffs = {
            'beta0': (0, 1.9e-9),      # Extraction at water interface (z=0)
            'betaL': (0, 1.0e-10),     # Evaporation at air interface (z=L)
        }

        # For film mode, both surfaces are in water
        if self.simulation_mode == 'film':
            self.beta_coeffs['betaL'] = self.beta_coeffs['beta0']

        # Yield parameters
        self.gamma = {
            'y1s': 1.0,    # Chain scission yield from alkoxy radicals
            'y1co': 0.61,  # Carbonyl yield
            'y4': 0.5,     # Crosslink yield from P• + P•
            'y5': 0.0,     # Crosslink yield from P• + PO2•
        }

        # DOC solubility parameters (from Colin Appendix A)
        self.S0_DOC = 2.6e-9   # mol/(L·Pa) - DOC solubility pre-exponential
        self.ES_DOC = -2690    # K - DOC solubility activation (Hs/R)
        self.pd0_DOC = 5.7e4   # Pa/ppm - Partial pressure conversion
        self.Ep_DOC = 26440    # J/mol - Vaporization energy

        # Stoichiometric parameters
        self.n_AH = 4           # Antioxidant functionality (Irganox 1010)
        self.PH0_conc = 60.0    # Initial PE concentration (mol/L)
        self.POOH0_conc = 1e-2  # Initial POOH concentration (mol/L)
        self.O2_sat_conc = 3.8e-4   # O2 saturation concentration (mol/L)
        self.AH0_conc = 1.8e-3  # Base antioxidant concentration (mol/L)

    def _arrhenius(self, Ea: float, A: float, T_K: float) -> float:
        """Calculate Arrhenius rate constant."""
        if Ea > 0:
            return A * np.exp(-Ea / (self.R_GAS * T_K))
        return A

    def _update_params_for_temp(self, T_K: float):
        """Update all rate constants and diffusion coefficients for temperature."""
        if T_K == self.current_T_K:
            return

        self.current_T_K = T_K

        # Rate constants
        for name, (Ea, A) in self.k_coeffs.items():
            self.k[name] = self._arrhenius(Ea, A, T_K)

        # Diffusion coefficients
        for name, (Ea, D0) in self.D_coeffs.items():
            self.D[name] = self._arrhenius(Ea, D0, T_K)

        # Boundary coefficients
        for name, (Ea, b0) in self.beta_coeffs.items():
            self.beta[name] = self._arrhenius(Ea, b0, T_K)

        # DOC solubility at T
        self.Sd_DOC_T = self.S0_DOC * np.exp(-self.ES_DOC / T_K)
        self.pd_DOC_T = self.pd0_DOC * np.exp(-self.Ep_DOC / (self.R_GAS * T_K))

    def _calc_DOC_boundary_conc(self, DOC_ppm: float) -> float:
        """
        Calculate DOC equilibrium concentration in PE from water concentration.

        Uses Henry's law and DOC solubility model from Colin Appendix A.
        """
        if DOC_ppm <= 0:
            return 0.0

        # DOC concentration in amorphous phase (mol/L_amorphous)
        doc_conc_amorphous = self.Sd_DOC_T * self.pd_DOC_T * DOC_ppm

        # Convert to total PE volume basis (mol/L_PE_total)
        return doc_conc_amorphous * self.Tam

    def _system_equations(self, t: float, y: np.ndarray) -> np.ndarray:
        """
        ODE system for diffusion-reaction equations.

        Implements the full kinetic scheme from Colin et al. (2009) Table 1
        with diffusion-reaction coupling for O2, DOC, and AH.
        """
        nz = self.nz
        dz = self.dz
        idx = self.idx
        k = self.k
        D = self.D
        beta = self.beta
        gamma = self.gamma
        n_AH = self.n_AH

        # Reshape to (n_species, nz)
        C = y.reshape((self.n_species, nz))
        dCdt = np.zeros_like(C)

        # Extract species concentrations
        O2 = C[idx['O2']]
        DOC = C[idx['DOC']]
        AH = C[idx['AH']]
        P = C[idx['P']]
        PO2 = C[idx['PO2']]
        POOH = C[idx['POOH']]
        PH = C[idx['PH']]
        Q = C[idx['Q']]

        # === Interior points (diffusion + reaction) ===
        for i in range(1, nz - 1):
            # Diffusion terms (Fick's second law)
            diff_O2 = D['O2'] * (O2[i+1] - 2*O2[i] + O2[i-1]) / dz**2
            diff_DOC = D['DOC'] * (DOC[i+1] - 2*DOC[i] + DOC[i-1]) / dz**2
            diff_AH = D['AH'] * (AH[i+1] - 2*AH[i] + AH[i-1]) / dz**2

            # O2: Eq. A-1 from Part II
            dCdt[idx['O2'], i] = diff_O2 - k['k2']*O2[i]*P[i] + k['k60']*PO2[i]**2

            # DOC: Eq. A-2
            dCdt[idx['DOC'], i] = diff_DOC \
                - k['k1d']*DOC[i]*PH[i] \
                - k['k4d']*P[i]*DOC[i] \
                - n_AH*k['k8d']*DOC[i]*AH[i]

            # AH: Eq. A-3
            dCdt[idx['AH'], i] = diff_AH \
                - n_AH*k['k8d']*DOC[i]*AH[i] \
                - n_AH*k['k7']*PO2[i]*AH[i]

            # P•: Eq. A-5 / B1
            dCdt[idx['P'], i] = \
                k['k1d']*DOC[i]*PH[i] \
                + 2*k['k1u']*POOH[i] \
                + k['k1b']*POOH[i]**2 \
                - k['k2']*O2[i]*P[i] \
                + k['k3']*PH[i]*PO2[i] \
                - 2*k['k4']*P[i]**2 \
                - k['k4d']*P[i]*DOC[i] \
                - k['k5']*P[i]*PO2[i] \
                + 2*k['k63']*Q[i]

            # PO2•: Eq. A-4 / B2
            dCdt[idx['PO2'], i] = \
                k['k1b']*POOH[i]**2 \
                + k['k2']*O2[i]*P[i] \
                - k['k3']*PH[i]*PO2[i] \
                - k['k5']*P[i]*PO2[i] \
                - 2*k['k60']*PO2[i]**2 \
                - n_AH*k['k7']*PO2[i]*AH[i]

            # POOH: Eq. B3
            dCdt[idx['POOH'], i] = \
                - k['k1u']*POOH[i] \
                - 2*k['k1b']*POOH[i]**2 \
                + k['k3']*PH[i]*PO2[i] \
                + (1 - gamma['y5'])*k['k5']*P[i]*PO2[i]

            # PH: Eq. A-8 / B5
            dCdt[idx['PH'], i] = \
                - k['k1d']*DOC[i]*PH[i] \
                - (2 + gamma['y1s'])*k['k1u']*POOH[i] \
                - (1 + gamma['y1s'])*k['k1b']*POOH[i]**2 \
                - k['k3']*PH[i]*PO2[i] \
                + 2*gamma['y4']*k['k4']*P[i]**2 \
                + (3*gamma['y5'] - 1)*k['k5']*P[i]*PO2[i] \
                + 2*k['k61']*Q[i] \
                - 2*(1 + gamma['y1s'])*k['k63']*Q[i]

            # Q (caged radical pair): Eq. B4
            dCdt[idx['Q'], i] = \
                k['k60']*PO2[i]**2 \
                - (k['k61'] + k['k62'] + k['k63'])*Q[i]

            # CO (carbonyl): Eq. B6
            dCdt[idx['CO'], i] = \
                gamma['y1co']*k['k1u']*POOH[i] \
                + gamma['y1co']*k['k1b']*POOH[i]**2 \
                + k['k62']*Q[i] \
                + 2*gamma['y1co']*k['k63']*Q[i]

            # PCl (grafted chlorine): Eq. B7
            dCdt[idx['PCl'], i] = k['k4d']*P[i]*DOC[i]

            # S (chain scissions): Eq. B8
            dCdt[idx['S'], i] = \
                gamma['y1s']*k['k1u']*POOH[i] \
                + gamma['y1s']*k['k1b']*POOH[i]**2 \
                + 2*gamma['y1s']*k['k63']*Q[i]

            # X (crosslinks): Eq. B9
            dCdt[idx['X'], i] = \
                gamma['y4']*k['k4']*P[i]**2 \
                + gamma['y5']*k['k5']*P[i]*PO2[i] \
                + k['k61']*Q[i]

        # === Boundary conditions ===
        self._apply_boundary_z0(dCdt, C)
        self._apply_boundary_zL(dCdt, C)

        return dCdt.flatten()

    def _apply_boundary_z0(self, dCdt: np.ndarray, C: np.ndarray):
        """Apply boundary conditions at z=0 (inner surface, water interface)."""
        idx = self.idx
        k = self.k
        beta = self.beta
        gamma = self.gamma
        n_AH = self.n_AH

        # Extract local concentrations
        O2, DOC, AH = C[idx['O2'], 0], C[idx['DOC'], 0], C[idx['AH'], 0]
        P, PO2, POOH = C[idx['P'], 0], C[idx['PO2'], 0], C[idx['POOH'], 0]
        PH, Q = C[idx['PH'], 0], C[idx['Q'], 0]

        # Dirichlet for O2 and DOC (fixed concentration)
        dCdt[idx['O2'], 0] = 0
        dCdt[idx['DOC'], 0] = 0

        # AH: extraction + reactions
        dCdt[idx['AH'], 0] = -beta['beta0']*AH \
            - n_AH*k['k8d']*DOC*AH \
            - n_AH*k['k7']*PO2*AH

        # Immobile species: reactions only
        dCdt[idx['P'], 0] = k['k1d']*DOC*PH + 2*k['k1u']*POOH + k['k1b']*POOH**2 \
            - k['k2']*O2*P + k['k3']*PH*PO2 - 2*k['k4']*P**2 \
            - k['k4d']*P*DOC - k['k5']*P*PO2 + 2*k['k63']*Q

        dCdt[idx['PO2'], 0] = k['k1b']*POOH**2 + k['k2']*O2*P \
            - k['k3']*PH*PO2 - k['k5']*P*PO2 - 2*k['k60']*PO2**2 \
            - n_AH*k['k7']*PO2*AH

        dCdt[idx['POOH'], 0] = -k['k1u']*POOH - 2*k['k1b']*POOH**2 \
            + k['k3']*PH*PO2 + (1 - gamma['y5'])*k['k5']*P*PO2

        dCdt[idx['Q'], 0] = k['k60']*PO2**2 - (k['k61'] + k['k62'] + k['k63'])*Q

        dCdt[idx['PH'], 0] = -k['k1d']*DOC*PH - (2 + gamma['y1s'])*k['k1u']*POOH \
            - (1 + gamma['y1s'])*k['k1b']*POOH**2 - k['k3']*PH*PO2 \
            + 2*gamma['y4']*k['k4']*P**2 + (3*gamma['y5'] - 1)*k['k5']*P*PO2 \
            + 2*k['k61']*Q - 2*(1 + gamma['y1s'])*k['k63']*Q

        dCdt[idx['CO'], 0] = gamma['y1co']*k['k1u']*POOH \
            + gamma['y1co']*k['k1b']*POOH**2 + k['k62']*Q + 2*gamma['y1co']*k['k63']*Q

        dCdt[idx['PCl'], 0] = k['k4d']*P*DOC

        dCdt[idx['S'], 0] = gamma['y1s']*k['k1u']*POOH \
            + gamma['y1s']*k['k1b']*POOH**2 + 2*gamma['y1s']*k['k63']*Q

        dCdt[idx['X'], 0] = gamma['y4']*k['k4']*P**2 \
            + gamma['y5']*k['k5']*P*PO2 + k['k61']*Q

    def _apply_boundary_zL(self, dCdt: np.ndarray, C: np.ndarray):
        """Apply boundary conditions at z=L (outer surface)."""
        idx = self.idx
        k = self.k
        D = self.D
        beta = self.beta
        gamma = self.gamma
        n_AH = self.n_AH
        dz = self.dz

        # Extract local concentrations
        O2, DOC, AH = C[idx['O2'], -1], C[idx['DOC'], -1], C[idx['AH'], -1]
        P, PO2, POOH = C[idx['P'], -1], C[idx['PO2'], -1], C[idx['POOH'], -1]
        PH, Q = C[idx['PH'], -1], C[idx['Q'], -1]

        # O2: Dirichlet (fixed concentration)
        dCdt[idx['O2'], -1] = 0

        if self.simulation_mode == 'film':
            # Film: DOC from both sides (Dirichlet)
            dCdt[idx['DOC'], -1] = 0
            # AH: extraction (same as z=0)
            dCdt[idx['AH'], -1] = -beta['betaL']*AH \
                - n_AH*k['k8d']*DOC*AH - n_AH*k['k7']*PO2*AH
        else:
            # Pipe: no DOC from air side (Neumann)
            DOC_prev = C[idx['DOC'], -2]
            dCdt[idx['DOC'], -1] = D['DOC']*2.0*(DOC_prev - DOC)/dz**2 \
                - k['k1d']*DOC*PH - k['k4d']*P*DOC - n_AH*k['k8d']*DOC*AH
            # AH: evaporation
            dCdt[idx['AH'], -1] = -beta['betaL']*AH \
                - n_AH*k['k8d']*DOC*AH - n_AH*k['k7']*PO2*AH

        # Immobile species: reactions only
        dCdt[idx['P'], -1] = k['k1d']*DOC*PH + 2*k['k1u']*POOH + k['k1b']*POOH**2 \
            - k['k2']*O2*P + k['k3']*PH*PO2 - 2*k['k4']*P**2 \
            - k['k4d']*P*DOC - k['k5']*P*PO2 + 2*k['k63']*Q

        dCdt[idx['PO2'], -1] = k['k1b']*POOH**2 + k['k2']*O2*P \
            - k['k3']*PH*PO2 - k['k5']*P*PO2 - 2*k['k60']*PO2**2 \
            - n_AH*k['k7']*PO2*AH

        dCdt[idx['POOH'], -1] = -k['k1u']*POOH - 2*k['k1b']*POOH**2 \
            + k['k3']*PH*PO2 + (1 - gamma['y5'])*k['k5']*P*PO2

        dCdt[idx['Q'], -1] = k['k60']*PO2**2 - (k['k61'] + k['k62'] + k['k63'])*Q

        dCdt[idx['PH'], -1] = -k['k1d']*DOC*PH - (2 + gamma['y1s'])*k['k1u']*POOH \
            - (1 + gamma['y1s'])*k['k1b']*POOH**2 - k['k3']*PH*PO2 \
            + 2*gamma['y4']*k['k4']*P**2 + (3*gamma['y5'] - 1)*k['k5']*P*PO2 \
            + 2*k['k61']*Q - 2*(1 + gamma['y1s'])*k['k63']*Q

        dCdt[idx['CO'], -1] = gamma['y1co']*k['k1u']*POOH \
            + gamma['y1co']*k['k1b']*POOH**2 + k['k62']*Q + 2*gamma['y1co']*k['k63']*Q

        dCdt[idx['PCl'], -1] = k['k4d']*P*DOC

        dCdt[idx['S'], -1] = gamma['y1s']*k['k1u']*POOH \
            + gamma['y1s']*k['k1b']*POOH**2 + 2*gamma['y1s']*k['k63']*Q

        dCdt[idx['X'], -1] = gamma['y4']*k['k4']*P**2 \
            + gamma['y5']*k['k5']*P*PO2 + k['k61']*Q

    def simulate(self, params: SimulationParams, verbose: bool = True) -> dict:
        """
        Run a degradation simulation.

        Args:
            params: SimulationParams object with simulation conditions
            verbose: Print progress information

        Returns:
            Dictionary with:
                - 'success': bool
                - 't': time array (seconds)
                - 't_years': time array (years)
                - 'y': solution array (n_species, nz, n_times)
                - 'params': copy of input parameters
                - 'AH0_used': actual initial AH concentration
                - 'DOC_boundary': DOC boundary concentration
        """
        T_K = params.T_celsius + 273.15
        self._update_params_for_temp(T_K)

        # Calculate boundary/initial concentrations
        O2_boundary = self.O2_sat_conc * params.O2_sat_mult
        AH0_actual = self.AH0_conc * params.AH0_mult
        DOC_boundary = self._calc_DOC_boundary_conc(params.DOC_ppm)

        # Initialize concentration array
        y0 = np.zeros(self.n_species * self.nz)
        C0 = y0.reshape((self.n_species, self.nz))

        # Set initial conditions
        C0[self.idx['O2'], :] = O2_boundary
        C0[self.idx['DOC'], :] = 0.0  # DOC diffuses in from boundary
        C0[self.idx['AH'], :] = AH0_actual
        C0[self.idx['POOH'], :] = self.POOH0_conc
        C0[self.idx['PH'], :] = self.PH0_conc

        # Set boundary values (Dirichlet conditions)
        C0[self.idx['O2'], 0] = O2_boundary
        C0[self.idx['O2'], -1] = O2_boundary
        C0[self.idx['DOC'], 0] = DOC_boundary

        if self.simulation_mode == 'film':
            C0[self.idx['DOC'], -1] = DOC_boundary

        # Time span
        t_end_sec = params.t_end_years * 365.25 * 24 * 3600
        t_eval = np.linspace(0, t_end_sec, max(2, params.n_timepoints))

        if verbose:
            print(f"Running simulation: T={params.T_celsius}°C, DOC={params.DOC_ppm} ppm, "
                  f"t={params.t_end_years:.3f} years, mode='{self.simulation_mode}'")

        start_time = time.time()

        sol = solve_ivp(
            self._system_equations,
            (0, t_end_sec),
            C0.flatten(),
            method=params.method,
            t_eval=t_eval,
            rtol=params.rtol,
            atol=params.atol
        )

        elapsed = time.time() - start_time

        if verbose:
            status = "SUCCESS" if sol.success else f"FAILED: {sol.message}"
            print(f"  Simulation {status} ({elapsed:.2f}s)")

        # Reshape solution
        y_reshaped = sol.y.reshape((self.n_species, self.nz, -1))

        return {
            'success': sol.success,
            't': sol.t,
            't_years': sol.t / (365.25 * 24 * 3600),
            't_months': sol.t / (365.25 * 24 * 3600) * 12,
            'y': y_reshaped,
            'params': params,
            'AH0_used': AH0_actual,
            'DOC_boundary': DOC_boundary,
            'message': sol.message if not sol.success else None
        }

    # =========================================================================
    # Post-processing methods
    # =========================================================================

    def calculate_oit(
        self,
        AH_profile: np.ndarray,
        AH0_reference: float
    ) -> np.ndarray:
        """
        Calculate Oxidation Induction Time from antioxidant profile.

        Uses Eq. (1) from Colin Part II:
            OIT(t) = ti0 * [AH](t) / [AH]0

        Args:
            AH_profile: Antioxidant concentration array
            AH0_reference: Reference initial AH concentration

        Returns:
            OIT array in minutes
        """
        if AH0_reference <= 1e-12:
            return np.zeros_like(AH_profile)

        ratio = np.maximum(0, AH_profile) / AH0_reference
        return self.material.ti0_oit * ratio

    def calculate_oit_profile_evolution(self, result: dict) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate OIT profiles at all time points.

        Args:
            result: Simulation result dictionary

        Returns:
            (times_months, oit_profiles) where oit_profiles is (nz, n_times)
        """
        AH_data = result['y'][self.idx['AH'], :, :]  # (nz, n_times)
        AH0 = result['AH0_used']

        oit_profiles = self.calculate_oit(AH_data, AH0)

        return result['t_months'], oit_profiles

    def calculate_average_oit(self, result: dict) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate average OIT across thickness at all time points.

        Args:
            result: Simulation result dictionary

        Returns:
            (times_months, oit_avg) arrays
        """
        _, oit_profiles = self.calculate_oit_profile_evolution(result)
        oit_avg = np.mean(oit_profiles, axis=0)
        return result['t_months'], oit_avg

    def calculate_molar_mass(self, result: dict) -> np.ndarray:
        """
        Calculate molecular weight evolution using Saito equation.

        Mw(t) = 1 / (1/Mw0 + s/2 - 2x)

        Args:
            result: Simulation result dictionary

        Returns:
            Mw array (nz, n_times) in kg/mol
        """
        S = result['y'][self.idx['S'], :, :]   # Chain scissions
        X = result['y'][self.idx['X'], :, :]   # Crosslinks

        Mw0 = self.material.Mw0
        Mw_inv = 1/Mw0 + S/2 - 2*X

        # Protect against negative or very large values
        Mw_inv = np.maximum(Mw_inv, 1e-6)
        return 1.0 / Mw_inv

    def get_surface_concentration(
        self,
        result: dict,
        species: str,
        surface: str = 'inner'
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get concentration of a species at a surface over time.

        Args:
            result: Simulation result dictionary
            species: Species name (e.g., 'CO', 'AH')
            surface: 'inner' (z=0) or 'outer' (z=L)

        Returns:
            (times_months, concentration) arrays
        """
        idx = 0 if surface == 'inner' else -1
        conc = result['y'][self.idx[species], idx, :]
        return result['t_months'], conc

    def interpolate_at_times(
        self,
        result: dict,
        target_months: np.ndarray,
        quantity: str = 'oit_avg'
    ) -> np.ndarray:
        """
        Interpolate simulation results at specific time points.

        Args:
            result: Simulation result dictionary
            target_months: Array of target times in months
            quantity: 'oit_avg', 'oit_surface', 'co_surface'

        Returns:
            Interpolated values at target times
        """
        sim_months = result['t_months']

        if quantity == 'oit_avg':
            _, values = self.calculate_average_oit(result)
        elif quantity == 'oit_surface':
            _, oit_profiles = self.calculate_oit_profile_evolution(result)
            values = oit_profiles[0, :]  # Inner surface
        elif quantity == 'co_surface':
            _, values = self.get_surface_concentration(result, 'CO', 'inner')
        else:
            raise ValueError(f"Unknown quantity: {quantity}")

        return np.interp(target_months, sim_months, values)


# =============================================================================
# Convenience functions
# =============================================================================

def create_suez_film_model(ti0_oit: float = 291.07) -> DegradationModel:
    """
    Create a model configured for SUEZ PE100 film experiments.

    SUEZ experiments use:
    - 400 μm thick films
    - PE100 material
    - Initial OIT ~291 min at 200°C (measured value with N2 correction)

    Args:
        ti0_oit: Initial OIT in minutes (default from SUEZ PE100)

    Returns:
        Configured DegradationModel instance
    """
    material = MaterialProperties(
        Mw0=150.0,
        dens0=0.95,
        densa=0.85,
        Xc=0.45,
        ti0_oit=ti0_oit,
        carbon_black=0.018  # 1.8% from report
    )

    return DegradationModel(
        L=0.4e-3,  # 400 μm
        nz=50,
        simulation_mode='film',
        material=material
    )


def create_colin_pipe_model() -> DegradationModel:
    """
    Create a model configured for Colin et al. (2009) pipe experiments.

    Returns:
        Configured DegradationModel instance
    """
    material = MaterialProperties(
        Mw0=150.0,
        dens0=0.95,
        densa=0.85,
        Xc=0.45,
        ti0_oit=165.0,  # From Colin paper
        carbon_black=0.025
    )

    return DegradationModel(
        L=4.5e-3,  # 4.5 mm pipe wall
        nz=100,
        simulation_mode='pipe',
        material=material
    )
