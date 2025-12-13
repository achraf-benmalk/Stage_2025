"""
Stage 3: Chemo-Mechanical Coupling Framework

This module implements multi-physics coupling between:
- Chemical degradation (kinetic model → Mw profiles)
- Mechanical property evolution (Fayolle correlations)
- Structural analysis (FEM-based failure prediction)

Components:
- DegradedMaterialProperties: Mw → mechanical property correlations
- StructuralAnalysis: FEM solver for pipe geometry
- ChemoMechanicalCoupling: Complete multi-physics chain
"""

from .DegradedMaterialProperties import DegradedMaterial, PropertyField
from .StructuralAnalysis import PipeStructuralModel, AxisymmetricFEM
from .ChemoMechanicalCoupling import MultiPhysicsCoupling

__all__ = [
    'DegradedMaterial',
    'PropertyField',
    'PipeStructuralModel',
    'AxisymmetricFEM',
    'MultiPhysicsCoupling'
]
