# Able to run multiple retrievals at once

import numpy as np
from POSEIDON.core import create_star, create_planet, load_data, wl_grid_constant_R, define_model, read_opacities, make_atmosphere, set_priors
from POSEIDON.constants import R_Sun, R_J
from POSEIDON.utility import read_retrieved_spectrum, plot_collection
from POSEIDON.visuals import plot_spectra_retrieved, plot_data
from POSEIDON.corner import generate_cornerplot, generate_overplot
from POSEIDON.retrieval import run_retrieval

# --- Define stellar properties --- #
R_s = 0.687*R_Sun       # Stellar radius [m]
T_s = 4250.0            # Stellar effective temperature [K]
err_T_s = 70            # Stellar effective temperature error [K]
Met_s = 0.2             # Stellar metallicity [log10(Fe/H_star / Fe/H_solar)]
log_g_s = 4.60          # Stellar log surface gravity [log10(cm/s^20)]

star = create_star(R_s, T_s, log_g_s, Met_s, T_eff_error = err_T_s)


# --- Define planetary properties --- #
planet_name = 'TOI-1130b'
R_p = 0.327*R_J         # Planetary radius [m]
# g_p = 14.4825           # Surface gravity [m/s^2]
log_g_p = 3.1609        # Log surface gravity [log10(g / cm/s^2)]  
T_eq = 825.0            # Equilibrium temperature [K]

planet = create_planet(planet_name, R_p, log_g = log_g_p, T_eq = T_eq)


# --- Model wavelength grid --- #
wl_min = 0.6            # Min and max wavelengths
wl_max = 5.2           
R = 10000               # Spectral resolution of grid

wl = wl_grid_constant_R(wl_min, wl_max, R)


# --- Specify data location and instruments --- #
data_dir = 'home/mlau/POSEIDON/POSEIDON/reference_data/observations/TOI-1130b'
datasets = ['TOI-1130b_Order2.dat', 'TOI-1130b_Order1.dat', 'TOI-1130b_NRS1.dat', 'TOI-1130b_NRS2.dat']
instruments = ['JWST_NIRISS_SOSS_Ord2', 'JWST_NIRISS_SOSS_Ord1', 'JWST_NIRSpec_G395H_NRS1', 'JWST_NIRSpec_G395H_NRS2']

# Load dataset and instrument transmission functions
data = load_data(data_dir, datasets, instruments, wl, 
                 offset_1_datasets = ['TOI-1130b_Order2.dat'], offset_2_datasets = ['TOI-1130b_Order1.dat'], offset_3_datasets = ['TOI-1130b_NRS2.dat'])

# Plot data
fig_data = plot_data(data, planet_name)


# --- Define model --- #
model_name = 'water_only'

bulk_species = ['H2', 'He']
param_species = ['H2O']
# stellar_contam = 'one_spot'

# Create model - comment/uncomment clouds and stellar contamination 
model = define_model(model_name, bulk_species, param_species, 
                     PT_profile = 'isotherm', X_profile = 'isochem',
                     # cloud_model = 'MacMad17', cloud_type = 'deck', cloud_dim = 1, 
                     offsets_applied = 'three_datasets', 
                     # stellar_contam = stellar_contam
                    )


# --- Read opacity data --- #
opacity_treatment = 'opacity_sampling'

# Define fine temp grid, [K], and fine pressure grid, [log10(P/bar)]
Tf_min = 300
Tf_max = 1500
Tf_step = 10

T_fine = np.arange(Tf_min, Tf_max + Tf_step, Tf_step)

log_Pf_min = -6.0
log_Pf_max = 2.0
log_Pf_step = 0.2

log_P_fine = np.arange(log_Pf_min, log_Pf_max + log_Pf_step, log_Pf_step)

# Pre-interpolate opacities
opac = read_opacities(model, wl, opacity_treatment, T_fine, log_P_fine)


# --- Specify fixed atmospheric settings for retrieval --- #
P_min = np.log10(1.0e-7)
P_max = np.log10(100)
N_layers = 100

P = np.logspace(P_max, P_min, N_layers)
P_ref = 10


# --- Set priors for retrieval --- #
print('Free parameters: ' + str(model['param_names']))

# Prior types - uniform or gaussian
prior_types = {}
prior_types['T'] = 'uniform'
prior_types['R_p_ref'] = 'uniform'
prior_types['log_X'] = 'uniform'
# prior_types['log_P_cloud'] = 'uniform'
prior_types['delta_rel_1'] = 'uniform'
prior_types['delta_rel_2'] = 'uniform'
prior_types['delta_rel_3'] = 'uniform'

# Prior ranges - [Min, Max], or [Mean, Stdev]
prior_ranges = {}
prior_ranges['T'] = [300, 1500]
prior_ranges['R_p_ref'] = [0.05*R_p, 5*R_p]
prior_ranges['log_X'] = [-12, -1]
# prior_ranges['log_P_cloud'] = [-8, 1]
prior_ranges['delta_rel_1'] = [-1000, 1000]
prior_ranges['delta_rel_2'] = [-1000, 1000]
prior_ranges['delta_rel_3'] = [-1000, 1000]

priors = set_priors(planet, star, model, data, prior_types, prior_ranges)


# --- Run retrieval --- #
run_retrieval(planet, star, model, opac, data, priors, wl, P, P_ref, R = R, 
              spectrum_type = 'transmission', sampling_algorithm = 'MultiNest', N_live = 1000, resume = False, verbose = True)

# Generate corner plot
fig_corner = generate_cornerplot(planet, model)

# ----------------------------------------------------------------------------------------

# --- Define second model --- # --> comment/uncomment run_retrieval lines depending on how many you want to do
model_name_2 = 'add_stuff'

bulk_species = ['H2', 'He']
param_species_2 = ['H2O', 'CO2', 'SO2', 'CH4']

model_2 = define_model(model_name_2, bulk_species, param_species_2, 
                       PT_profile = 'isotherm', X_profile = 'isochem', 
                       # cloud_model = 'MacMad17', cloud_type = 'deck', cloud_dim = 1,
                       offsets_applied = 'three_datasets', 
                       # stellar_contam = stellar_contam
                      )

# Read opacity data for new model
opac_2 = read_opacities(model_2, wl, opacity_treatment, T_fine, log_P_fine)

# Set new priors
print('Free parameters: ' + str(model_2['param_names']))

prior_types_2 = {}
prior_types_2['T'] = 'uniform'
prior_types_2['R_p_ref'] = 'uniform'
prior_types_2['log_X'] = 'uniform'
# prior_types_2['log_P_cloud'] = 'uniform'
prior_types_2['delta_rel_1'] = 'uniform'
prior_types_2['delta_rel_2'] = 'uniform'
prior_types_2['delta_rel_3'] = 'uniform'

prior_ranges_2 = {}
prior_ranges_2['T'] = [300, 1500]
prior_ranges_2['R_p_ref'] = [0.05*R_p, 5*R_p]
prior_ranges_2['log_X'] = [-12, -1]
# prior_ranges['log_P_cloud'] = [-8, 1]
prior_ranges_2['delta_rel_1'] = [-1000, 1000]
prior_ranges_2['delta_rel_2'] = [-1000, 1000]
prior_ranges_2['delta_rel_3'] = [-1000, 1000]

priors_2 = set_priors(planet, star, model_2, data, prior_types_2, prior_ranges_2)

# Run second retrieval
run_retrieval(planet, star, model_2, opac_2, data, priors_2, wl, P, P_ref, R = R, 
              spectrum_type = 'transmission', sampling_algorithm = 'MultiNest', N_live = 1000, resume = False, verbose = True)

# Generate corner plot
fig_corner = generate_cornerplot(planet, model_2)

# ----------------------------------------------------------------------------------------

# --- Define third model --- # --> comment/uncomment run_retrieval lines depending on how many you want to do
model_name_3 = 'add_morestuff'

bulk_species = ['H2', 'He']
param_species_3 = ['H2O', 'CO2', 'SO2', 'CH4' ,'CO', 'NH3']

model_3 = define_model(model_name_3, bulk_species, param_species_3, 
                       PT_profile = 'isotherm', X_profile = 'isochem', 
                       # cloud_model = 'MacMad17', cloud_type = 'deck', cloud_dim = 1,
                       offsets_applied = 'three_datasets', 
                       # stellar_contam = stellar_contam
                      )

# Read opacity data for new model
opac_3 = read_opacities(model_3, wl, opacity_treatment, T_fine, log_P_fine)

# Set new priors
print('Free parameters: ' + str(model_3['param_names']))

prior_types_3 = {}
prior_types_3['T'] = 'uniform'
prior_types_3['R_p_ref'] = 'uniform'
prior_types_3['log_X'] = 'uniform'
# prior_types_3['log_P_cloud'] = 'uniform'
prior_types_3['delta_rel_1'] = 'uniform'
prior_types_3['delta_rel_2'] = 'uniform'
prior_types_3['delta_rel_3'] = 'uniform'

prior_ranges_3 = {}
prior_ranges_3['T'] = [300, 1500]
prior_ranges_3['R_p_ref'] = [0.05*R_p, 5*R_p]
prior_ranges_3['log_X'] = [-12, -1]
# prior_ranges['log_P_cloud'] = [-8, 1]
prior_ranges_3'delta_rel_1'] = [-1000, 1000]
prior_ranges_3['delta_rel_2'] = [-1000, 1000]
prior_ranges_3['delta_rel_3'] = [-1000, 1000]

priors_3 = set_priors(planet, star, model_3, data, prior_types_3, prior_ranges_3)

# Run second retrieval
run_retrieval(planet, star, model_3, opac_3, data, priors_3, wl, P, P_ref, R = R, 
              spectrum_type = 'transmission', sampling_algorithm = 'MultiNest', N_live = 1000, resume = False, verbose = True)

# Generate corner plot
fig_corner = generate_cornerplot(planet, model_3)

# ----------------------------------------------------------------------------------------

# --- Define fourth model --- # --> comment/uncomment run_retrieval lines depending on how many you want to do
model_name_4 = 'barat_molecules'

bulk_species = ['H2', 'He']
param_species_4 = ['H2O', 'CO2', 'SO2', 'CH4', 'CO', 'NH3', 'H2S', 'HCN', 'CS2']

model_4 = define_model(model_name_4, bulk_species, param_species_4, 
                       PT_profile = 'isotherm', X_profile = 'isochem', 
                       # cloud_model = 'MacMad17', cloud_type = 'deck', cloud_dim = 1,
                       offsets_applied = 'three_datasets', 
                       # stellar_contam = stellar_contam
                      )

# Read opacity data for new model
opac_4 = read_opacities(model_4, wl, opacity_treatment, T_fine, log_P_fine)

# Set new priors
print('Free parameters: ' + str(model_4['param_names']))

prior_types_4 = {}
prior_types_4['T'] = 'uniform'
prior_types_4['R_p_ref'] = 'uniform'
prior_types_4['log_X'] = 'uniform'
# prior_types_4['log_P_cloud'] = 'uniform'
prior_types_4['delta_rel_1'] = 'uniform'
prior_types_4['delta_rel_2'] = 'uniform'
prior_types_4['delta_rel_3'] = 'uniform'

prior_ranges_4 = {}
prior_ranges_4['T'] = [300, 1500]
prior_ranges_4['R_p_ref'] = [0.05*R_p, 5*R_p]
prior_ranges_4['log_X'] = [-12, -1]
# prior_ranges['log_P_cloud'] = [-8, 1]
prior_ranges_4'delta_rel_1'] = [-1000, 1000]
prior_ranges_4['delta_rel_2'] = [-1000, 1000]
prior_ranges_4['delta_rel_3'] = [-1000, 1000]

priors_4 = set_priors(planet, star, model_4, data, prior_types_4, prior_ranges_4)

# Run second retrieval
run_retrieval(planet, star, model_4, opac_4, data, priors_4, wl, P, P_ref, R = R, 
              spectrum_type = 'transmission', sampling_algorithm = 'MultiNest', N_live = 1000, resume = False, verbose = True)

# Generate corner plot
fig_corner = generate_cornerplot(planet, model_4)
