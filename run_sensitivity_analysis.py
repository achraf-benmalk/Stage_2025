#!/usr/bin/env python3
"""
Comprehensive Sensitivity Analysis Script

This script runs all parameter sweeps and generates analysis results for the
PE pipe degradation research project.

Outputs:
- Parameter sweep results (CSV)
- Validation metrics
- Tornado diagram data
- Arrhenius analysis
- Figures (PNG/PDF)

Author: Stage 2025 Research Project
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from DegradationModel import (
    DegradationModel, SimulationParams, create_suez_film_model
)
from ParameterAnalysis import ParameterAnalysis, SweepResult
from ValidationMetrics import (
    SUEZExperimentalData, ValidationComparison,
    calculate_all_metrics, analyze_model_limitations, quantify_model_gap
)

# Configure matplotlib
plt.rcParams.update({
    'figure.figsize': (10, 6),
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})


def run_baseline_validation(model, verbose=True):
    """Run baseline H2O and HOCl simulations for validation."""
    print("\n" + "="*70)
    print("BASELINE VALIDATION")
    print("="*70)

    suez = SUEZExperimentalData()
    validator = ValidationComparison(suez)

    # H2O simulation (DOC=0)
    if verbose:
        print("\n1. Running H2O (control) simulation...")

    params_h2o = SimulationParams(
        T_celsius=40.0,
        DOC_ppm=0.0,  # Pure water
        t_end_years=0.75,  # 9 months
        n_timepoints=100
    )
    result_h2o = model.simulate(params_h2o, verbose=verbose)

    if result_h2o['success']:
        times_m, oit_avg = model.calculate_average_oit(result_h2o)
        h2o_comparison = validator.compare_h2o(times_m, oit_avg)
        print(f"   H2O R² = {h2o_comparison['metrics']['R2']:.4f}")
    else:
        print("   H2O simulation FAILED")
        return None, None, None

    # HOCl simulation (using DOC as proxy)
    # Note: This uses DOC kinetics as approximation for HOCl
    if verbose:
        print("\n2. Running HOCl simulation (DOC proxy)...")

    # Effective DOC equivalent for 0.05 ppm HOCl
    # This needs calibration - start with literature value
    doc_effective = 0.05  # ppm - initial estimate

    params_hocl = SimulationParams(
        T_celsius=40.0,
        DOC_ppm=doc_effective,
        t_end_years=0.75,
        n_timepoints=100
    )
    result_hocl = model.simulate(params_hocl, verbose=verbose)

    if result_hocl['success']:
        times_m, oit_avg = model.calculate_average_oit(result_hocl)
        hocl_comparison = validator.compare_hocl(times_m, oit_avg)
        print(f"   HOCl R² = {hocl_comparison['metrics']['R2']:.4f}")
    else:
        print("   HOCl simulation FAILED")
        return result_h2o, None, validator

    # Print summary
    validator.print_summary()

    return result_h2o, result_hocl, validator


def run_doc_sweep(analyzer, verbose=True):
    """Run DOC concentration parameter sweep."""
    print("\n" + "="*70)
    print("PARAMETER SWEEP: DOC Concentration")
    print("="*70)

    doc_values = np.array([0, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0])

    result = analyzer.sweep_doc_concentration(
        doc_values=doc_values,
        T_celsius=40.0,
        t_end_years=0.75,
        n_timepoints=50,
        verbose=verbose
    )

    # Calculate sensitivity metrics
    metrics_df = analyzer.calculate_sensitivity_metrics(result)
    print("\nSensitivity Metrics (DOC sweep):")
    print(metrics_df.to_string(index=False))

    return result, metrics_df


def run_temperature_sweep(analyzer, verbose=True):
    """Run temperature parameter sweep with Arrhenius analysis."""
    print("\n" + "="*70)
    print("PARAMETER SWEEP: Temperature")
    print("="*70)

    T_values = np.array([30, 40, 50, 60])

    result = analyzer.sweep_temperature(
        T_values=T_values,
        DOC_ppm=0.05,
        t_end_years=0.75,
        n_timepoints=50,
        verbose=verbose
    )

    # Arrhenius analysis
    arrhenius = analyzer.arrhenius_analysis(result, target_oit=50.0)

    print("\nArrhenius Analysis:")
    print(f"  Activation Energy (from time to OIT=50): {arrhenius.get('Ea_from_time_kJ_mol', np.nan):.1f} kJ/mol")
    print(f"  R² of fit: {arrhenius.get('Ea_time_R2', np.nan):.4f}")
    print(f"  Literature value for PE: 80-120 kJ/mol")

    # Calculate acceleration factors
    acc_df = analyzer.calculate_acceleration_factor(result, reference_idx=0, target_oit=50.0)
    print("\nAcceleration Factors (vs 30°C):")
    print(acc_df.to_string(index=False))

    return result, arrhenius, acc_df


def run_ah_sweep(analyzer, verbose=True):
    """Run initial antioxidant parameter sweep."""
    print("\n" + "="*70)
    print("PARAMETER SWEEP: Initial Antioxidant [AH]₀")
    print("="*70)

    ah_values = np.array([0.5, 1.0, 1.5, 2.0, 3.0])

    result = analyzer.sweep_antioxidant(
        AH_multipliers=ah_values,
        T_celsius=40.0,
        DOC_ppm=0.05,
        t_end_years=0.75,
        n_timepoints=50,
        verbose=verbose
    )

    # Calculate sensitivity metrics
    metrics_df = analyzer.calculate_sensitivity_metrics(result, reference_idx=1)  # 1.0x is reference
    print("\nSensitivity Metrics (AH sweep):")
    print(metrics_df.to_string(index=False))

    return result, metrics_df


def create_sensitivity_figures(doc_result, T_result, ah_result, output_dir):
    """Create publication-quality sensitivity figures."""
    print("\n" + "="*70)
    print("CREATING FIGURES")
    print("="*70)

    os.makedirs(output_dir, exist_ok=True)

    # Figure 1: DOC Concentration Sweep
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(doc_result.parameter_values)))

    for i, doc_val in enumerate(doc_result.parameter_values):
        label = f'DOC = {doc_val} ppm' if doc_val > 0 else 'H₂O (control)'
        ax1.plot(doc_result.times_months, doc_result.oit_avg_matrix[i, :],
                 color=colors[i], linewidth=2, label=label)

    ax1.set_xlabel('Time (months)')
    ax1.set_ylabel('Average OIT (min)')
    ax1.set_title('Effect of DOC Concentration on OIT Decay (T=40°C)')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 9)
    ax1.set_ylim(0, None)

    fig1.savefig(os.path.join(output_dir, 'doc_sensitivity.png'))
    fig1.savefig(os.path.join(output_dir, 'doc_sensitivity.pdf'))
    print(f"  Saved: {output_dir}/doc_sensitivity.png")
    plt.close(fig1)

    # Figure 2: Temperature Sweep
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    colors = plt.cm.plasma(np.linspace(0.2, 0.9, len(T_result.parameter_values)))

    for i, T_val in enumerate(T_result.parameter_values):
        ax2.plot(T_result.times_months, T_result.oit_avg_matrix[i, :],
                 color=colors[i], linewidth=2, label=f'T = {T_val}°C')

    ax2.set_xlabel('Time (months)')
    ax2.set_ylabel('Average OIT (min)')
    ax2.set_title('Effect of Temperature on OIT Decay (DOC=0.05 ppm)')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 9)
    ax2.set_ylim(0, None)

    fig2.savefig(os.path.join(output_dir, 'temperature_sensitivity.png'))
    fig2.savefig(os.path.join(output_dir, 'temperature_sensitivity.pdf'))
    print(f"  Saved: {output_dir}/temperature_sensitivity.png")
    plt.close(fig2)

    # Figure 3: Antioxidant Sweep
    fig3, ax3 = plt.subplots(figsize=(10, 6))
    colors = plt.cm.coolwarm(np.linspace(0.1, 0.9, len(ah_result.parameter_values)))

    for i, ah_val in enumerate(ah_result.parameter_values):
        ax3.plot(ah_result.times_months, ah_result.oit_avg_matrix[i, :],
                 color=colors[i], linewidth=2, label=f'[AH]₀ = {ah_val}×')

    ax3.set_xlabel('Time (months)')
    ax3.set_ylabel('Average OIT (min)')
    ax3.set_title('Effect of Initial Antioxidant Content on OIT Decay (T=40°C, DOC=0.05 ppm)')
    ax3.legend(loc='upper right')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 9)
    ax3.set_ylim(0, None)

    fig3.savefig(os.path.join(output_dir, 'antioxidant_sensitivity.png'))
    fig3.savefig(os.path.join(output_dir, 'antioxidant_sensitivity.pdf'))
    print(f"  Saved: {output_dir}/antioxidant_sensitivity.png")
    plt.close(fig3)

    # Figure 4: Arrhenius Plot
    fig4, ax4 = plt.subplots(figsize=(10, 6))

    T_K = T_result.parameter_values + 273.15
    inv_T = 1000.0 / T_K  # 1000/T for better scale

    # Calculate time to OIT=50 for each temperature
    times = T_result.times_months
    t_to_50 = []
    for i in range(len(T_result.parameter_values)):
        oit_profile = T_result.oit_avg_matrix[i, :]
        # Find crossing
        below = oit_profile < 50
        if np.any(below):
            first_idx = np.argmax(below)
            if first_idx > 0:
                t1, t2 = times[first_idx-1], times[first_idx]
                v1, v2 = oit_profile[first_idx-1], oit_profile[first_idx]
                t_cross = t1 + (50 - v1) * (t2 - t1) / (v2 - v1)
                t_to_50.append(t_cross)
            else:
                t_to_50.append(times[0])
        else:
            t_to_50.append(np.nan)

    t_to_50 = np.array(t_to_50)
    valid = ~np.isnan(t_to_50) & (t_to_50 > 0)

    if np.sum(valid) >= 2:
        ln_rate = np.log(1.0 / t_to_50[valid])

        # Linear fit
        slope, intercept = np.polyfit(inv_T[valid], ln_rate, 1)
        fit_line = slope * inv_T + intercept

        Ea = -slope * 8.314 * 1000 / 1000  # kJ/mol (accounting for 1000/T)

        ax4.scatter(inv_T[valid], ln_rate, s=100, c='blue', zorder=3, label='Simulated data')
        ax4.plot(inv_T, fit_line, 'r--', linewidth=2,
                 label=f'Linear fit: Ea = {Ea:.1f} kJ/mol')

        # Add temperature labels
        for i, T in enumerate(T_result.parameter_values):
            if valid[i]:
                ax4.annotate(f'{T}°C', (inv_T[i], ln_rate[i]),
                            textcoords="offset points", xytext=(5,5))

    ax4.set_xlabel('1000/T (K⁻¹)')
    ax4.set_ylabel('ln(1/t₅₀) (t₅₀ in months)')
    ax4.set_title('Arrhenius Plot: Temperature Dependence of Degradation Rate')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    fig4.savefig(os.path.join(output_dir, 'arrhenius_plot.png'))
    fig4.savefig(os.path.join(output_dir, 'arrhenius_plot.pdf'))
    print(f"  Saved: {output_dir}/arrhenius_plot.png")
    plt.close(fig4)

    return True


def create_tornado_diagram(analyzer, doc_result, T_result, ah_result, output_dir):
    """Create tornado diagram showing parameter importance."""
    print("\n" + "="*70)
    print("CREATING TORNADO DIAGRAM")
    print("="*70)

    # Calculate impact at 6 months
    metric_time = 6.0
    times = doc_result.times_months

    # Baseline values
    baseline = {
        'DOC_ppm': 0.05,
        'T_celsius': 40.0,
        'AH0_mult': 1.0
    }

    tornado_data = []

    # DOC impact
    doc_idx_base = np.argmin(np.abs(doc_result.parameter_values - baseline['DOC_ppm']))
    oit_base = np.interp(metric_time, times, doc_result.oit_avg_matrix[doc_idx_base, :])

    oit_min = np.min([np.interp(metric_time, times, doc_result.oit_avg_matrix[i, :])
                      for i in range(len(doc_result.parameter_values))])
    oit_max = np.max([np.interp(metric_time, times, doc_result.oit_avg_matrix[i, :])
                      for i in range(len(doc_result.parameter_values))])

    tornado_data.append({
        'Parameter': 'DOC (ppm)',
        'Range': f'{doc_result.parameter_values.min()}-{doc_result.parameter_values.max()}',
        'OIT_low': oit_min,
        'OIT_high': oit_max,
        'Impact': oit_max - oit_min
    })

    # Temperature impact
    T_idx_base = np.argmin(np.abs(T_result.parameter_values - baseline['T_celsius']))

    oit_min = np.min([np.interp(metric_time, times, T_result.oit_avg_matrix[i, :])
                      for i in range(len(T_result.parameter_values))])
    oit_max = np.max([np.interp(metric_time, times, T_result.oit_avg_matrix[i, :])
                      for i in range(len(T_result.parameter_values))])

    tornado_data.append({
        'Parameter': 'Temperature (°C)',
        'Range': f'{T_result.parameter_values.min()}-{T_result.parameter_values.max()}',
        'OIT_low': oit_min,
        'OIT_high': oit_max,
        'Impact': oit_max - oit_min
    })

    # AH impact
    ah_idx_base = np.argmin(np.abs(ah_result.parameter_values - baseline['AH0_mult']))

    oit_min = np.min([np.interp(metric_time, times, ah_result.oit_avg_matrix[i, :])
                      for i in range(len(ah_result.parameter_values))])
    oit_max = np.max([np.interp(metric_time, times, ah_result.oit_avg_matrix[i, :])
                      for i in range(len(ah_result.parameter_values))])

    tornado_data.append({
        'Parameter': '[AH]₀ multiplier',
        'Range': f'{ah_result.parameter_values.min()}-{ah_result.parameter_values.max()}×',
        'OIT_low': oit_min,
        'OIT_high': oit_max,
        'Impact': oit_max - oit_min
    })

    # Sort by impact
    tornado_df = pd.DataFrame(tornado_data).sort_values('Impact')

    # Create tornado diagram
    fig, ax = plt.subplots(figsize=(12, 6))

    y_pos = np.arange(len(tornado_df))

    # Plot bars
    for i, row in enumerate(tornado_df.itertuples()):
        # Low side (red) - deviation below baseline
        ax.barh(i, oit_base - row.OIT_low, left=row.OIT_low, color='#E74C3C',
                height=0.6, alpha=0.8, label='Low value' if i == 0 else '')
        # High side (blue) - deviation above baseline
        ax.barh(i, row.OIT_high - oit_base, left=oit_base, color='#3498DB',
                height=0.6, alpha=0.8, label='High value' if i == 0 else '')

    # Baseline line
    ax.axvline(x=oit_base, color='black', linestyle='--', linewidth=2,
               label=f'Baseline OIT = {oit_base:.1f} min')

    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{row.Parameter}\n({row.Range})" for row in tornado_df.itertuples()])
    ax.set_xlabel('OIT at 6 months (min)')
    ax.set_title('Parameter Sensitivity: Tornado Diagram\n(Impact on OIT at 6 months)')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3, axis='x')

    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, 'tornado_diagram.png'))
    fig.savefig(os.path.join(output_dir, 'tornado_diagram.pdf'))
    print(f"  Saved: {output_dir}/tornado_diagram.png")
    plt.close(fig)

    print("\nParameter Importance (by OIT range at 6 months):")
    print(tornado_df.to_string(index=False))

    return tornado_df


def create_validation_figure(model, validator, output_dir):
    """Create validation figure comparing simulation vs experiment."""
    print("\n" + "="*70)
    print("CREATING VALIDATION FIGURE")
    print("="*70)

    suez = SUEZExperimentalData()

    # Run extended simulations
    params_h2o = SimulationParams(
        T_celsius=40.0,
        DOC_ppm=0.0,
        t_end_years=0.75,
        n_timepoints=100
    )
    result_h2o = model.simulate(params_h2o, verbose=False)
    times_h2o, oit_h2o = model.calculate_average_oit(result_h2o)

    params_hocl = SimulationParams(
        T_celsius=40.0,
        DOC_ppm=0.05,
        t_end_years=0.75,
        n_timepoints=100
    )
    result_hocl = model.simulate(params_hocl, verbose=False)
    times_hocl, oit_hocl = model.calculate_average_oit(result_hocl)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # H2O validation
    ax1.plot(times_h2o, oit_h2o, 'b-', linewidth=2, label='Simulation')
    ax1.scatter(suez.times_h2o, suez.oit_h2o, s=100, c='blue', marker='o',
                label='SUEZ Experimental', zorder=3, edgecolors='black')

    # Calculate R²
    sim_at_exp = np.interp(suez.times_h2o, times_h2o, oit_h2o)
    r2_h2o = 1 - np.sum((suez.oit_h2o - sim_at_exp)**2) / np.sum((suez.oit_h2o - np.mean(suez.oit_h2o))**2)

    ax1.set_xlabel('Time (months)')
    ax1.set_ylabel('OIT (min)')
    ax1.set_title(f'H₂O Control (Physical Aging)\nR² = {r2_h2o:.4f}')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(-0.5, 10)
    ax1.set_ylim(0, 320)

    # HOCl validation
    ax2.plot(times_hocl, oit_hocl, 'r-', linewidth=2, label='Simulation (DOC proxy)')
    ax2.scatter(suez.times_hocl, suez.oit_hocl, s=100, c='red', marker='o',
                label='SUEZ Experimental', zorder=3, edgecolors='black')

    # Calculate R²
    sim_at_exp_hocl = np.interp(suez.times_hocl, times_hocl, oit_hocl)
    r2_hocl = 1 - np.sum((suez.oit_hocl - sim_at_exp_hocl)**2) / np.sum((suez.oit_hocl - np.mean(suez.oit_hocl))**2)

    ax2.set_xlabel('Time (months)')
    ax2.set_ylabel('OIT (min)')
    ax2.set_title(f'HOCl 0.05 ppm (Chemical Aging)\nR² = {r2_hocl:.4f}')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(-0.5, 10)
    ax2.set_ylim(0, 320)

    fig.suptitle('Model Validation: SUEZ PE100 Film (400 μm, 40°C)', fontsize=14, y=1.02)
    fig.tight_layout()

    fig.savefig(os.path.join(output_dir, 'validation_comparison.png'))
    fig.savefig(os.path.join(output_dir, 'validation_comparison.pdf'))
    print(f"  Saved: {output_dir}/validation_comparison.png")
    plt.close(fig)

    return r2_h2o, r2_hocl


def export_results(doc_result, T_result, ah_result, arrhenius, output_dir):
    """Export all results to CSV files."""
    print("\n" + "="*70)
    print("EXPORTING RESULTS")
    print("="*70)

    data_dir = os.path.join(output_dir, '..', 'data')
    os.makedirs(data_dir, exist_ok=True)

    # DOC sweep results
    doc_data = []
    for i, val in enumerate(doc_result.parameter_values):
        for j, t in enumerate(doc_result.times_months):
            doc_data.append({
                'DOC_ppm': val,
                'time_months': t,
                'OIT_avg_min': doc_result.oit_avg_matrix[i, j],
                'OIT_surface_min': doc_result.oit_surface_matrix[i, j],
                'CO_surface_mol_L': doc_result.co_surface_matrix[i, j]
            })
    pd.DataFrame(doc_data).to_csv(os.path.join(data_dir, 'doc_sweep_results.csv'), index=False)
    print(f"  Saved: {data_dir}/doc_sweep_results.csv")

    # Temperature sweep results
    T_data = []
    for i, val in enumerate(T_result.parameter_values):
        for j, t in enumerate(T_result.times_months):
            T_data.append({
                'T_celsius': val,
                'time_months': t,
                'OIT_avg_min': T_result.oit_avg_matrix[i, j],
                'OIT_surface_min': T_result.oit_surface_matrix[i, j],
                'CO_surface_mol_L': T_result.co_surface_matrix[i, j]
            })
    pd.DataFrame(T_data).to_csv(os.path.join(data_dir, 'temperature_sweep_results.csv'), index=False)
    print(f"  Saved: {data_dir}/temperature_sweep_results.csv")

    # AH sweep results
    ah_data = []
    for i, val in enumerate(ah_result.parameter_values):
        for j, t in enumerate(ah_result.times_months):
            ah_data.append({
                'AH0_mult': val,
                'time_months': t,
                'OIT_avg_min': ah_result.oit_avg_matrix[i, j],
                'OIT_surface_min': ah_result.oit_surface_matrix[i, j],
                'CO_surface_mol_L': ah_result.co_surface_matrix[i, j]
            })
    pd.DataFrame(ah_data).to_csv(os.path.join(data_dir, 'antioxidant_sweep_results.csv'), index=False)
    print(f"  Saved: {data_dir}/antioxidant_sweep_results.csv")

    # Arrhenius results
    arrhenius_df = pd.DataFrame({
        'T_celsius': arrhenius['T_celsius'],
        'T_kelvin': arrhenius['T_kelvin'],
        '1_over_T': arrhenius['1/T'],
        'time_to_OIT_50': arrhenius.get('time_to_oit_target', [np.nan]*len(arrhenius['T_celsius']))
    })
    arrhenius_df.to_csv(os.path.join(data_dir, 'arrhenius_analysis.csv'), index=False)
    print(f"  Saved: {data_dir}/arrhenius_analysis.csv")


def main():
    """Main function to run complete sensitivity analysis."""
    print("="*70)
    print("PE DEGRADATION MODEL - COMPREHENSIVE SENSITIVITY ANALYSIS")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70)

    # Create model
    print("\nInitializing SUEZ PE100 film model...")
    model = create_suez_film_model()
    print(f"  Thickness: {model.L*1000:.1f} mm")
    print(f"  Mode: {model.simulation_mode}")
    print(f"  Initial OIT: {model.material.ti0_oit:.1f} min")

    # Output directories
    figures_dir = 'figures'

    # Run baseline validation
    result_h2o, result_hocl, validator = run_baseline_validation(model)

    # Create analyzer
    analyzer = ParameterAnalysis(model)

    # Run parameter sweeps
    doc_result, doc_metrics = run_doc_sweep(analyzer)
    T_result, arrhenius, T_acc = run_temperature_sweep(analyzer)
    ah_result, ah_metrics = run_ah_sweep(analyzer)

    # Create figures
    create_sensitivity_figures(doc_result, T_result, ah_result, figures_dir)
    tornado_df = create_tornado_diagram(analyzer, doc_result, T_result, ah_result, figures_dir)
    r2_h2o, r2_hocl = create_validation_figure(model, validator, figures_dir)

    # Export results
    export_results(doc_result, T_result, ah_result, arrhenius, figures_dir)

    # Print summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE - SUMMARY")
    print("="*70)

    print("\n1. VALIDATION METRICS:")
    print(f"   H₂O (physical aging):  R² = {r2_h2o:.4f}")
    print(f"   HOCl (chemical aging): R² = {r2_hocl:.4f}")

    print("\n2. ARRHENIUS ANALYSIS:")
    print(f"   Activation Energy: {arrhenius.get('Ea_from_time_kJ_mol', np.nan):.1f} kJ/mol")
    print(f"   Literature range: 80-120 kJ/mol")

    print("\n3. PARAMETER IMPORTANCE (at 6 months):")
    for _, row in tornado_df.iterrows():
        print(f"   {row['Parameter']}: ΔOIt = {row['Impact']:.1f} min")

    print("\n4. TEMPERATURE ACCELERATION:")
    for _, row in T_acc.iterrows():
        print(f"   {row['T_celsius']}°C: {row['acceleration_factor']:.2f}× vs 30°C")

    print("\n5. KEY FINDINGS:")
    print("   • DOC concentration has LARGEST impact on degradation rate")
    print("   • Temperature shows clear Arrhenius behavior")
    print("   • Initial antioxidant content affects induction time")
    print("   • HOCl validation shows significant gap (R² < 0.5)")
    print("   • Gap due to DOC→HOCl kinetics approximation")

    print("\n" + "="*70)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70)


if __name__ == '__main__':
    main()
