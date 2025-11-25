## This module contains functions for plotting simulation results, including concentration profiles and OIT validation.

import numpy as np
import matplotlib.pyplot as plt
## This function plots concentration profiles vs. depth for specified species at different times.
def plot_profiles_evolution(model_object, sol, species_to_plot, 
                            title_suffix="", 
                            time_points_years=None, 
                            selected_time_indices=None):
    """
    Plots concentration profiles vs. depth for specified species at different times.
    """
    if not sol.success:
        print("Cannot plot results from failed simulation.")
        return None, None

    # Use attributes from the passed model_object
    n_species_model = model_object.n_species
    nz_model = model_object.nz
    z_coords_m = model_object.z
    species_idx_map = model_object.idx

    times_sec = sol.t
    y_data = np.asarray(sol.y) 
    C_all_times = y_data.reshape((n_species_model, nz_model, len(times_sec)))
    z_mm = z_coords_m * 1000

    if not species_to_plot:
        print("No species specified for plotting.")
        return None, None

    # Select time indices to plot
    if selected_time_indices is not None:
        indices_to_plot = selected_time_indices
    elif time_points_years is not None:
        times_years = times_sec / (365.25 * 24 * 3600)
        indices_to_plot = []
        for target_year in time_points_years:
            idx = np.argmin(np.abs(times_years - target_year))
            indices_to_plot.append(idx)
    else:
        # Default: 5 evenly spaced time points
        num_points = min(5, len(times_sec))
        indices_to_plot = np.linspace(0, len(times_sec)-1, num_points, dtype=int)

    # Setup subplot grid
    num_species_to_plot = len(species_to_plot)
    cols_per_row = 3 

    if num_species_to_plot <= cols_per_row:
        num_rows_plot = 1
        num_cols_plot = num_species_to_plot
    else:
        num_rows_plot = (num_species_to_plot + cols_per_row - 1) // cols_per_row
        num_cols_plot = cols_per_row

    fig, axes_grid = plt.subplots(num_rows_plot, num_cols_plot, 
                                 figsize=(6 * num_cols_plot, 5 * num_rows_plot),
                                 sharex=True, squeeze=False)
    axes_list = axes_grid.flatten()

    # Plot each species
    colors = plt.cm.viridis(np.linspace(0, 1, len(indices_to_plot)))
    
    for i, species_name in enumerate(species_to_plot):
        if i >= len(axes_list):
            break
        
        ax = axes_list[i]
        species_idx = species_idx_map[species_name]
        
        # Plot concentration profiles at different times
        for j, time_idx in enumerate(indices_to_plot):
            time_years = times_sec[time_idx] / (365.25 * 24 * 3600)
            concentration = C_all_times[species_idx, :, time_idx]
            ax.plot(z_mm, concentration, color=colors[j], 
                   label=f't = {time_years:.2f} yr', linewidth=2)
        
        ax.set_title(species_name, fontweight='bold')
        ax.set_ylabel('Concentration')
        ax.grid(True, alpha=0.3)
        
        # Add legend on first plot of each row
        if i % num_cols_plot == 0:
            ax.legend(loc='upper right', fontsize='x-small')
        
        # Add x-label only on bottom row
        current_row = i // num_cols_plot
        if current_row == num_rows_plot - 1:
            ax.set_xlabel('Depth (mm)')
        else:
            ax.tick_params(labelbottom=False)

    # Hide unused subplots
    for j in range(num_species_to_plot, len(axes_list)):
        fig.delaxes(axes_list[j])

    # Set title and layout
    fig_title = f"Concentration Profiles Evolution{' - ' + title_suffix if title_suffix else ''}"
    fig.suptitle(fig_title, fontsize=12)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    return fig, axes_grid
## This function plots OIT validation comparing experimental and simulated data.
def plot_OIT_H20_validation(model_object, experimental_times_months, experimental_oit, simulated_oit, 
                           temp_celsius, pipe_thickness_m, save_filename=None):
    """
    Plot OIT validation comparing experimental and simulated data
    """
    # Use attributes from the passed model_object
    ti0_oit = model_object.ti0_oit
    
    plt.figure(figsize=(10, 6))
    
    # Plot experimental data (filter NaN)
    valid_exp_indices = ~np.isnan(experimental_oit)
    if np.any(valid_exp_indices):
        plt.plot(experimental_times_months[valid_exp_indices], experimental_oit[valid_exp_indices], 
                 marker='o', linestyle='-', color='blue', label='OIT Expérimental SUEZ (Eau Pure)')
    
    # Plot simulated data
    valid_sim_indices = ~np.isnan(simulated_oit)
    if np.any(valid_sim_indices):
        plt.plot(experimental_times_months[valid_sim_indices], simulated_oit[valid_sim_indices], 
                 marker='s', linestyle='--', color='red', label=f'OIT Simulé (Eau Pure, DOC=0ppm)')

    plt.xlabel("Temps de Vieillissement (mois)")
    plt.ylabel("OIT (min)")
    plt.title(f"Validation Modèle: OIT vs. Temps en Eau Pure (PE100 Film {pipe_thickness_m*1000:.1f}mm, {temp_celsius}°C)")
    plt.legend()
    plt.grid(True, linestyle=':')
    
    # Adjust Y axis limits dynamically
    all_valid_oit_values = []
    if np.any(valid_exp_indices): all_valid_oit_values.extend(experimental_oit[valid_exp_indices])
    if np.any(valid_sim_indices): all_valid_oit_values.extend(simulated_oit[valid_sim_indices])
    if not all_valid_oit_values: all_valid_oit_values.append(ti0_oit)
    
    plt.ylim(bottom=0, top=max(all_valid_oit_values) * 1.15 if all_valid_oit_values else ti0_oit * 1.15)
    plt.xlim(left=-0.1, right=9.1)
    plt.xticks(np.arange(0, 10, 1))
    
    plt.tight_layout()
    plt.show()
    
    if save_filename:
        plt.savefig(save_filename, dpi=300, bbox_inches='tight')
## This function calculates OIT for plotting, adapted from the model class.
def calculate_oit_for_plot(AH_profile_at_time_t, AH0_initial_for_material, ti0_oit_material):
    """
    Helper function to calculate OIT.
    (Copied and adapted from your model class, removing 'model_object')
    """
    if AH0_initial_for_material <= 1e-12:
        return np.zeros_like(AH_profile_at_time_t)
    ah_ratio = np.maximum(0, AH_profile_at_time_t) / AH0_initial_for_material
    oit = ti0_oit_material * ah_ratio
    return np.maximum(0, oit)
## This function plots the OIT and Carbonyl validation for HOCl aging, comparing experimental data with model simulations.
def plot_suez_hocl_validation( 
                                 experimental_times_months,
                                 experimental_oit_hocl,
                                 experimental_carbonyl_hocl,
                                 solution_hocl_sim, # L'objet 'sol' de la simulation HOCl
                                 # Paramètres du modèle/matériau nécessaires pour le traitement/légendes
                                 model_n_species,
                                 model_nz,
                                 model_idx_AH,
                                 model_idx_CO,
                                 model_L_meters,
                                 model_ti0_oit, # OIT initial du matériau pour l'échelle de l'axe Y
                                 carbonyl_exp_label="Indice Carbonyle Exp. (u.a.)",
                                 title_extra=""):
    """
    Trace les OIT et les Carbonyles (expérimentaux et simulés) pour le vieillissement HOCl
    en confrontant les données SUEZ avec les simulations du modèle.
    Cette fonction est autonome et destinée à être dans un module.
    """
    if not (solution_hocl_sim and hasattr(solution_hocl_sim, 'success') and solution_hocl_sim.success and solution_hocl_sim.y.size > 0):
        print("plot_suez_hocl_validation: Simulation HOCl non fournie, non réussie, ou avec des résultats vides. Impossible de tracer.")
        return None, None, None

    # --- Extraction des données simulées ---
    sim_params = getattr(solution_hocl_sim, 'sim_params', {})
    sim_times_sec = solution_hocl_sim.t
    sim_times_months = sim_times_sec / (365.25 * 24 * 3600 / 12) 

    C_all_times_sim = solution_hocl_sim.y.reshape((model_n_species, model_nz, len(sim_times_sec)))
    
    AH_profiles_sim = C_all_times_sim[model_idx_AH, :, :]
    AH0_actual_conc_used_sim = sim_params.get('AH0_conc_used', 0) 
    if AH0_actual_conc_used_sim == 0:
        print("plot_suez_hocl_validation: Attention: AH0_conc_used non trouvé dans sim_params ou est zéro. OIT simulé sera incorrect.")
        # Pourrait être un fallback ici si vous avez une valeur de base pour le multiplicateur AH0
        # AH0_actual_conc_used_sim = fallback_AH0_base * sim_params.get('AH0_mult',1.0)


    simulated_oit_profiles_at_each_t_eval = np.array([
        calculate_oit_for_plot(AH_profiles_sim[:, t_idx], AH0_actual_conc_used_sim, model_ti0_oit) 
        for t_idx in range(len(sim_times_sec))
    ])
    simulated_oit_avg_vs_time = np.mean(simulated_oit_profiles_at_each_t_eval, axis=1)

    simulated_co_surface_vs_time = C_all_times_sim[model_idx_CO, 0, :] 

    # --- Création du Graphique ---
    fig, ax1 = plt.subplots(figsize=(12, 7))

    color_oit_exp = 'green' 
    color_oit_sim = 'purple' 
    ax1.set_xlabel("Temps de Vieillissement (mois)")
    ax1.set_ylabel("OIT (min)", color=color_oit_exp)
    
    valid_exp_oit_indices = ~np.isnan(experimental_oit_hocl)
    if np.any(valid_exp_oit_indices):
        ax1.plot(experimental_times_months[valid_exp_oit_indices], experimental_oit_hocl[valid_exp_oit_indices], 
                 marker='o', linestyle='-', color=color_oit_exp, label='OIT Expérimental SUEZ (HOCl)')
    
    valid_sim_oit_indices = ~np.isnan(simulated_oit_avg_vs_time)
    if np.any(valid_sim_oit_indices):
        ax1.plot(sim_times_months[valid_sim_oit_indices], simulated_oit_avg_vs_time[valid_sim_oit_indices], 
                 marker='s', linestyle='--', color=color_oit_sim, label=f'OIT Simulé (HOCl via DOCeff={sim_params.get("DOC_ppm", "N/A")}ppm)')
    ax1.tick_params(axis='y', labelcolor=color_oit_exp)
    
    all_valid_oit = []
    if np.any(valid_exp_oit_indices): all_valid_oit.extend(experimental_oit_hocl[valid_exp_oit_indices])
    if np.any(valid_sim_oit_indices): all_valid_oit.extend(simulated_oit_avg_vs_time[valid_sim_oit_indices])
    if not all_valid_oit: all_valid_oit.append(model_ti0_oit) 
    
    min_oit_val = min(all_valid_oit) if all_valid_oit else 0
    max_oit_val = max(all_valid_oit) if all_valid_oit else model_ti0_oit
    ax1.set_ylim(bottom=max(0, min_oit_val * 0.9 - 10), top=max_oit_val * 1.1 + 10)

    ax2 = ax1.twinx()
    color_co_exp = 'darkorange' 
    color_co_sim = 'brown'      
    ax2.set_ylabel(carbonyl_exp_label, color=color_co_exp)

    valid_exp_co_indices = ~np.isnan(experimental_carbonyl_hocl)
    if np.any(valid_exp_co_indices):
        ax2.plot(experimental_times_months[valid_exp_co_indices], experimental_carbonyl_hocl[valid_exp_co_indices], 
                 marker='^', linestyle='-', color=color_co_exp, label=carbonyl_exp_label)

    valid_sim_co_indices = ~np.isnan(simulated_co_surface_vs_time)
    if np.any(valid_sim_co_indices):
        ax2.plot(sim_times_months[valid_sim_co_indices], simulated_co_surface_vs_time[valid_sim_co_indices], 
                 marker='x', linestyle=':', color=color_co_sim, label='[CO] Simulé surface (mol/L)')
    ax2.tick_params(axis='y', labelcolor=color_co_exp)

    all_valid_co = []
    if np.any(valid_exp_co_indices): all_valid_co.extend(experimental_carbonyl_hocl[valid_exp_co_indices])
    if not all_valid_co and np.any(valid_sim_co_indices): 
        all_valid_co.extend(simulated_co_surface_vs_time[valid_sim_co_indices]) 
        ax2.set_ylabel("[CO] Simulé surface (mol/L)", color=color_co_sim) 
        ax2.tick_params(axis='y', labelcolor=color_co_sim)

    if not all_valid_co: all_valid_co.append(0.01) 
    
    min_co_val = 0 
    max_co_val = max(all_valid_co) if all_valid_co else 0.01
    ax2.set_ylim(bottom=0, top=max_co_val * 1.2 if max_co_val > 1e-9 else 0.01)
    
    T_sim_val = sim_params.get('T_celsius', 'N/A') 
    DOC_eff_sim_val = sim_params.get('DOC_ppm', 'N/A') 
    
    main_title = (f"Validation Modèle: OIT & Carbonyles vs. Temps en HOCl\n"
                  f"{title_extra} (Film {model_L_meters*1000:.1f}mm, {T_sim_val}°C, DOCeff={DOC_eff_sim_val}ppm)")
    ax1.set_title(main_title) 
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.18))

    all_times_months_plot = []
    all_times_months_plot.extend(experimental_times_months)
    all_times_months_plot.extend(sim_times_months) 
    if not all_times_months_plot: all_times_months_plot.append(9.0)
    
    max_x_plot_val_final = max(all_times_months_plot) if all_times_months_plot else 9.0
    ax1.set_xlim(left=-0.1, right=max_x_plot_val_final * 1.05 if max_x_plot_val_final > 0 else 9.1)
    ax1.set_xticks(np.arange(0, int(max_x_plot_val_final) + 2, 1) if max_x_plot_val_final > 0 else np.arange(0,10,1))
    
    ax1.grid(True, linestyle=':') 
    fig.tight_layout(rect=[0, 0, 1, 0.90]) 
    
    return fig, ax1, ax2

# Dans votre fichier plot_functions.py
## Cette fonction générique pour tracer la validation OIT & Carbonyles pour SUEZ.
def plot_suez_validation_combined(condition_label, experimental_times_months, experimental_oit, experimental_carbonyl, solution_sim, model, doc_eff_ppm, title_extra=""):
    """
    Fonction générique pour tracer la validation OIT & Carbonyles pour SUEZ.
    """
    # ... (Tout le code de votre ancienne fonction 'plot_suez_hocl_validation' va ici)
    # ... (Assurez-vous de remplacer les références directes à 'HOCl' par des variables)
    
    # Exemple de modifications à faire à l'intérieur de la fonction :
    
    # 1. Calculer les OITs et Carbonyles simulés à partir de 'solution_sim' et 'model'
    #    (Ce code ne change pas)
    times_sim_years = solution_sim.t
    times_sim_months = times_sim_years * 12
    oit_sim = model.calculate_oit_from_solution(solution_sim)
    carbonyl_sim = model.calculate_surface_concentration(solution_sim, 'CO')
    
    # 2. Créer la figure et les axes (ce code ne change pas)
    fig, ax1 = plt.subplots(figsize=(14, 7))
    ax2 = ax1.twinx()

    # 3. Filtrer les NaN pour les données expérimentales (important pour les données H2O)
    valid_indices_carbonyl = ~np.isnan(experimental_carbonyl)
    
    # 4. Tracer les données
    # OIT (axe gauche)
    ax1.plot(experimental_times_months, experimental_oit, 'o-', color='green', label=f'OIT Expérimental SUEZ ({condition_label})')
    ax1.plot(times_sim_months, oit_sim, 's--', color='purple', label=f'OIT Simulé ({condition_label} via DOCeff={doc_eff_ppm}ppm)')
    
    # Carbonyles (axe droit)
    ax2.plot(experimental_times_months[valid_indices_carbonyl], experimental_carbonyl[valid_indices_carbonyl], '^-', color='darkorange', label='Indice Carbonyle Exp. (Hauteur, u.a.)')
    ax2.plot(times_sim_months, carbonyl_sim, 'x:', color='darkred', label='[CO] Simulé surface (mol/L)')
    
    # 5. Mettre en forme le graphique
    #    Utilisez 'condition_label' pour rendre les titres et labels dynamiques
    ax1.set_title(f"Validation Modèle: OIT & Carbonyles vs. Temps en {condition_label}\n{title_extra}, {model.L*1000:.1f}mm, {model.T_C:.1f}°C, DOCeff={doc_eff_ppm}ppm")
    # ... (le reste des labels et de la mise en forme) ...
    
    return fig, ax1, ax2
## Cette fonction trace la validation OIT & Carbonyles pour SUEZ, combinant les données expérimentales et simulées.
def plot_suez_combined_validation( 
                                 experimental_times_months, # Les temps pour lesquels on a des données exp. OIT ET CO
                                 experimental_oit,
                                 experimental_carbonyl,
                                 solution_sim, # Objet solution de la simulation (peut couvrir plus de temps)
                                 
                                 model_n_species, model_nz, model_idx_AH, model_idx_CO,
                                 model_L_meters, model_ti0_oit,
                                 
                                 nom_fichier_figure="validation_plot.png",
                                 plot_title_prefix="Validation Modèle: OIT & Carbonyles vs. Temps",
                                 oit_exp_label='OIT Expérimental SUEZ',
                                 oit_sim_label='OIT Simulé',
                                 carbonyl_exp_label="Indice Carbonyle Exp. (u.a.)",
                                 co_sim_label='[CO] Simulé surface (mol/L)',
                                 color_oit_exp='blue', color_oit_sim='red',
                                 color_co_exp='green', color_co_sim='purple',
                                 title_extra=""):
    """
    Trace les OIT et les Carbonyles (expérimentaux et simulés) pour une condition donnée.
    Sauvegarde la figure.
    """
    if not (solution_sim and hasattr(solution_sim, 'success') and solution_sim.success and solution_sim.y.size > 0):
        print(f"{plot_title_prefix}: Simulation non fournie ou non réussie. Impossible de tracer.")
        return None, None, None

    sim_params = getattr(solution_sim, 'sim_params', {})
    sim_times_sec = solution_sim.t
    sim_times_months_simulation = sim_times_sec / (365.25 * 24 * 3600 / 12) 

    C_all_times_sim = solution_sim.y.reshape((model_n_species, model_nz, len(sim_times_sec)))
    
    AH_profiles_sim = C_all_times_sim[model_idx_AH, :, :]
    AH0_actual_conc_used_sim = sim_params.get('AH0_conc_used', 0)
    if AH0_actual_conc_used_sim == 0:
        print(f"{plot_title_prefix}: Attention: AH0_conc_used non trouvé. OIT simulé sera incorrect.")

    simulated_oit_profiles_at_each_t_eval = np.array([
        calculate_oit_for_plot(AH_profiles_sim[:, t_idx], AH0_actual_conc_used_sim, model_ti0_oit) 
        for t_idx in range(len(sim_times_sec))
    ])
    simulated_oit_avg_vs_time = np.mean(simulated_oit_profiles_at_each_t_eval, axis=1)
    simulated_co_surface_vs_time = C_all_times_sim[model_idx_CO, 0, :] 

    fig, ax1 = plt.subplots(figsize=(12, 7))

    ax1.set_xlabel("Temps de Vieillissement (mois)")
    ax1.set_ylabel("OIT (min)", color=color_oit_exp)
    
    valid_exp_oit_indices = ~np.isnan(experimental_oit)
    if np.any(valid_exp_oit_indices):
        ax1.plot(experimental_times_months[valid_exp_oit_indices], experimental_oit[valid_exp_oit_indices], 
                 marker='o', linestyle='-', color=color_oit_exp, label=oit_exp_label)
    
    valid_sim_oit_indices = ~np.isnan(simulated_oit_avg_vs_time)
    if np.any(valid_sim_oit_indices):
        ax1.plot(sim_times_months_simulation[valid_sim_oit_indices], simulated_oit_avg_vs_time[valid_sim_oit_indices], 
                 marker='s', linestyle='--', color=color_oit_sim, label=oit_sim_label)
    ax1.tick_params(axis='y', labelcolor=color_oit_exp)
    
    all_valid_oit = []
    if np.any(valid_exp_oit_indices): all_valid_oit.extend(experimental_oit[valid_exp_oit_indices])
    if np.any(valid_sim_oit_indices): all_valid_oit.extend(simulated_oit_avg_vs_time[valid_sim_oit_indices])
    if not all_valid_oit: all_valid_oit.append(model_ti0_oit) 
    
    min_oit_val = min(all_valid_oit) if all_valid_oit else 0
    max_oit_val = max(all_valid_oit) if all_valid_oit else model_ti0_oit
    ax1.set_ylim(bottom=max(0, min_oit_val * 0.9 - 10), top=max_oit_val * 1.15 + 10) # Un peu plus de marge

    ax2 = ax1.twinx()
    ax2.set_ylabel(carbonyl_exp_label, color=color_co_exp)

    valid_exp_co_indices = ~np.isnan(experimental_carbonyl)
    if np.any(valid_exp_co_indices):
        ax2.plot(experimental_times_months[valid_exp_co_indices], experimental_carbonyl[valid_exp_co_indices], 
                 marker='^', linestyle='-', color=color_co_exp, label=carbonyl_exp_label)

    valid_sim_co_indices = ~np.isnan(simulated_co_surface_vs_time)
    if np.any(valid_sim_co_indices):
        ax2.plot(sim_times_months_simulation[valid_sim_co_indices], simulated_co_surface_vs_time[valid_sim_co_indices], 
                 marker='x', linestyle=':', color=color_co_sim, label=co_sim_label)
    ax2.tick_params(axis='y', labelcolor=color_co_exp)

    all_valid_co = []
    if np.any(valid_exp_co_indices): all_valid_co.extend(experimental_carbonyl[valid_exp_co_indices])
    if not all_valid_co and np.any(valid_sim_co_indices): 
        all_valid_co.extend(simulated_co_surface_vs_time[valid_sim_co_indices])
        # Si seulement simulé est tracé, ajuster le label de l'axe
        if not np.any(valid_exp_co_indices):
            ax2.set_ylabel(co_sim_label, color=color_co_sim)
            ax2.tick_params(axis='y', labelcolor=color_co_sim)
            
    if not all_valid_co: all_valid_co.append(0.01) 
    
    max_co_val = max(all_valid_co) if all_valid_co else 0.01
    ax2.set_ylim(bottom=0, top=max_co_val * 1.2 if max_co_val > 1e-9 else 0.01)
    
    T_sim_val = sim_params.get('T_celsius', 'N/A') 
    DOC_sim_val = sim_params.get('DOC_ppm', 'N/A') # Sera 0 pour H2O
    
    title = (f"{plot_title_prefix}\n"
             f"{title_extra} (Film {model_L_meters*1000:.1f}mm, {T_sim_val}°C, DOCsim={DOC_sim_val}ppm)")
    ax1.set_title(title) 
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.22)) # Ajusté bbox

    max_x_plot_val_final = max(np.concatenate((experimental_times_months, sim_times_months_simulation))) if len(experimental_times_months)>0 or len(sim_times_months_simulation)>0 else 9.0
    ax1.set_xlim(left=-0.1, right=max_x_plot_val_final * 1.05 if max_x_plot_val_final > 0 else 9.1)
    ax1.set_xticks(np.arange(0, int(max_x_plot_val_final) + 2, 1) if max_x_plot_val_final > 0 else np.arange(0,10,1))
    
    ax1.grid(True, linestyle=':') 
    fig.tight_layout(rect=[0, 0, 1, 0.90]) 
    
    try:
        fig.savefig(nom_fichier_figure, dpi=300, bbox_inches='tight')
        print(f"\nGraphique sauvegardé sous : {nom_fichier_figure}")
    except Exception as e:
        print(f"Erreur lors de la sauvegarde du graphique '{nom_fichier_figure}': {e}")
    
    return fig, ax1, ax2