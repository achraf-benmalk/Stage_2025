"""
PE Pipe Degradation Analysis Package

This package provides tools for modeling and analyzing polyethylene pipe
degradation under chlorine exposure, following the Colin et al. (2009) framework.

Modules:
    - DegradationModel: 1D diffusion-reaction model
    - ParameterAnalysis: Sensitivity analysis framework
    - ValidationMetrics: Statistical validation tools
"""

from .DegradationModel import (
    DegradationModel,
    MaterialProperties,
    SimulationParams,
    create_suez_film_model,
    create_colin_pipe_model
)

from .ParameterAnalysis import (
    ParameterAnalysis,
    SweepResult,
    run_standard_sensitivity_study
)

from .ValidationMetrics import (
    SUEZExperimentalData,
    ValidationComparison,
    calculate_r2,
    calculate_rmse,
    calculate_mae,
    calculate_mape,
    calculate_all_metrics,
    analyze_model_limitations,
    quantify_model_gap
)

__version__ = '1.0.0'
__author__ = 'Stage 2025 Research Project'
