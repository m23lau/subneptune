import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

from POSEIDON.core import create_star, create_planet
from POSEIDON.constants import R_Sun, R_J

from POSEIDON.utility import read_retrieved_spectrum, plot_collection
from POSEIDON.visuals import plot_spectra_retrieved
from POSEIDON.corner import generate_cornerplot
from POSEIDON.core import load_data, wl_grid_constant_R
from POSEIDON.visuals import plot_data


#***** Define stellar properties *****#

R_s = 0.687*R_Sun     # Stellar radius (m)
T_s = 4250.0          # Stellar effective temperature (K)
err_T_s = 70
Met_s = 0.2         # Stellar metallicity [log10(Fe/H_star / Fe/H_solar)]
log_g_s = 4.60        # Stellar log surface gravity (log10(cm/s^2) by convention)

# Create the stellar object
star = create_star(R_s, T_s, log_g_s, Met_s, T_eff_error = err_T_s)

#***** Define planet properties *****#

planet_name = 'TOI-1130b'  # Planet name used for plots, output files etc.

R_p = 0.327*R_J     # Planetary radius (m)
# g_p = 14.4825
log_g_p = 3.1609          # Gravitational field of planet (m/s^2)
T_eq = 825.0       # Equilibrium temperature (K)

# Create the planet object
planet = create_planet(planet_name, R_p, log_g = log_g_p, T_eq = T_eq)

from POSEIDON.core import load_data, wl_grid_constant_R
from POSEIDON.visuals import plot_data

#***** Model wavelength grid *****#

wl_min = 0.6      # Minimum wavelength (um)
wl_max = 4.7      # Maximum wavelength (um)
R = 10000          # Spectral resolution of grid

# We need to provide a model wavelength grid to initialise instrument properties
wl = wl_grid_constant_R(wl_min, wl_max, R)


#***** Specify data location and instruments  *****#
data_dir = '/home/peng/POSEIDON/POSEIDON/reference_data/observations/TOI-1130b'       # Change this to where your data is stored
datasets = ['TOI-1130b_Order2.dat', 'TOI-1130b_Order1.dat', 'TOI-1130b_NRS1.dat', 'TOI-1130b_NRS2_cut.dat']                       # Found in reference_data/observations
instruments =['JWST_NIRISS_SOSS_Ord2', 'JWST_NIRISS_SOSS_Ord1', 'JWST_NIRSpec_G395H_NRS1', 'JWST_NIRSpec_G395H_NRS2']                   # Instruments corresponding to the data

# Load dataset, pre-load instrument PSF and transmission function
#data = load_data(data_dir, datasets, instruments, wl)
data = load_data(data_dir, datasets, instruments, wl, offset_1_datasets =['TOI-1130b_Order2.dat'], offset_2_datasets =['TOI-1130b_Order1.dat'], offset_3_datasets = ['TOI-1130b_NRS2_cut.dat'])

# Plot our data
fig_data = plot_data(data, planet_name)

from POSEIDON.core import define_model

#***** Define model *****#

model_name = 'full_wateronly_cloud'  # Model name used for plots, output files etc.

bulk_species = ['H2', 'He']     # H2 + He comprises the bulk atmosphere
param_species =['H2O']
#param_species = [ 'H2O', 'CO2', 'CO', 'CH4','SO2','NH3', 'C2H2', 'C2H4', 'H2S', 'HCN','CS2','OCS' ]         # The only trace gas is H2O
#stellar_contam = 'one_spot'

# Create the model object
model = define_model(model_name, bulk_species, param_species, 
                             PT_profile = 'isotherm', 
                             X_profile ='isochem',
                             cloud_model = 'MacMad17',
                             cloud_type = 'deck',
                             cloud_dim = 1,
                             offsets_applied = 'single_dataset',
                             #stellar_contam = stellar_contam
                             )

from POSEIDON.core import read_opacities
import numpy as np

#***** Read opacity data *****#

opacity_treatment = 'opacity_sampling'

# Define fine temperature grid (K)
T_fine_min = 300     # Same as prior range for T
T_fine_max = 1500    # Same as prior range for T
T_fine_step = 10     # 10 K steps are a good tradeoff between accuracy and RAM

T_fine = np.arange(T_fine_min, (T_fine_max + T_fine_step), T_fine_step)

# Define fine pressure grid (log10(P/bar))
log_P_fine_min = -6.0   # 1 ubar is the lowest pressure in the opacity database
log_P_fine_max = 100    # 100 bar is the highest pressure in the opacity database
log_P_fine_step = 0.2   # 0.2 dex steps are a good tradeoff between accuracy and RAM

log_P_fine = np.arange(log_P_fine_min, (log_P_fine_max + log_P_fine_step),
                               log_P_fine_step)

# Pre-interpolate the opacities
opac = read_opacities(model, wl, opacity_treatment, T_fine, log_P_fine)

from POSEIDON.retrieval import run_retrieval

#***** Specify fixed atmospheric settings for retrieval *****#

# Atmospheric pressure grid
P_min = 1.0e-7    # 0.1 ubar
P_max = 100       # 100 bar
N_layers = 100    # 100 layers

# Let's space the layers uniformly in log-pressure
P = np.logspace(np.log10(P_max), np.log10(P_min), N_layers)

# Specify the reference pressure
P_ref = 10   # Retrieved R_p_ref parameter will be the radius at 10 bar

from POSEIDON.core import make_atmosphere


# Check the free parameters defining this model
print("Free parameters: " + str(model['param_names']))

from POSEIDON.core import set_priors

#***** Set priors for retrieval *****#

# Initialise prior type dictionary
prior_types = {}

# Specify whether priors are linear, Gaussian, etc.
prior_types['T'] = 'uniform'
prior_types['R_p_ref'] = 'uniform'
prior_types['log_H2O'] = 'uniform'
#prior_types['log_CO'] = 'uniform'
#prior_types['log_CO2'] = 'uniform'
#prior_types['log_CH4'] = 'uniform'
#prior_types['log_SO2'] = 'uniform'
#prior_types['log_NH3'] = 'uniform'
#prior_types['log_C2H2'] = 'uniform'
#prior_types['log_C2H4'] = 'uniform'
#prior_types['log_H2S'] = 'uniform'
#prior_types['log_HCN'] = 'uniform'
#prior_types['log_CS2'] = 'uniform'
#prior_types['log_OCS'] = 'uniform'
prior_types['log_P_cloud'] = 'uniform'
prior_types['delta_rel'] = 'uniform'

#prior_types['T_het'] = 'uniform'
#prior_types['f_het'] = 'uniform'

#prior_types['T_phot'] = 'gaussian'
#prior_types['f_spot'] = 'uniform'
#prior_types['f_fac'] = 'uniform'
#prior_types['T_spot'] = 'uniform'
#prior_types['T_fac'] = 'uniform'

# Initialise prior range dictionary
prior_ranges = {}


# Specify prior ranges for each free parameter
prior_ranges['T'] = [300, 1500]
prior_ranges['R_p_ref'] = [0.05*R_p,5*R_p]
prior_ranges['log_H2O'] = [-12, -0.000000005]
#prior_ranges['log_CO'] = [-12, -0.005]
#prior_ranges['log_CO2'] = [-12, -0.005]
#prior_ranges['log_CH4'] = [-12, -0.005]
#prior_ranges['log_SO2'] = [-12, -0.005]
#prior_ranges['log_NH3'] = [-12, -0.005]
#prior_ranges['log_C2H2'] = [-12, -0.005]
#prior_ranges['log_C2H4'] = [-12, -0.005]
#prior_ranges['log_H2S'] = [-12, -0.005]
#prior_ranges['log_HCN'] = [-12, -0.005]
#prior_ranges['log_CS2'] = [-12, -0.005]
#prior_ranges['log_OCS'] = [-12, -0.005]
prior_ranges['log_P_cloud'] = [-8, 1]
prior_ranges['delta_rel'] =[-250,250]
#prior_ranges['f_het'] = [0.0, 0.5]
#prior_ranges['T_het'] = [2000, 1.2*T_s]
#prior_ranges['log_g_het'] = [3.0, 5.0]
#prior_ranges['T_phot'] = [3500, err_T_s]
#prior_ranges['f_spot'] = [0.0, 0.5]
#prior_ranges['f_fac'] = [0.0, 0.5]
#prior_ranges['T_spot'] = [3500, T_s+3*err_T_s]
#prior_ranges['T_fac'] = [3500, 1.2*T_s]
#prior_ranges['log_g_spot'] = [3.0, 5.0]
#prior_ranges['log_g_fac'] = [3.0, 5.0]

# Create prior object for retrieval
priors = set_priors(planet, star, model, data, prior_types, prior_ranges)


from POSEIDON.retrieval import run_retrieval

#***** Run atmospheric retrieval *****#

run_retrieval(planet, star, model, opac, data, priors, wl, P, P_ref, R = R,
                      spectrum_type = 'transmission', sampling_algorithm = 'MultiNest',
                                    N_live = 1000, resume= False, verbose = True,
                                    )


# Read retrieved spectrum confidence regions
wl, spec_low2, spec_low1, spec_median, \
spec_high1, spec_high2 = read_retrieved_spectrum(planet_name, model_name)

# Create composite spectra objects for plotting
spectra_median = plot_collection(spec_median, wl, collection = [])
spectra_low1 = plot_collection(spec_low1, wl, collection = [])
spectra_low2 = plot_collection(spec_low2, wl, collection = [])
spectra_high1 = plot_collection(spec_high1, wl, collection = [])
spectra_high2 = plot_collection(spec_high2, wl, collection = [])

# Produce figure
fig_spec = plot_spectra_retrieved(spectra_median, spectra_low2, spectra_low1,
                                  spectra_high1, spectra_high2, planet_name,
                                  data, R_to_bin = 100,
                                  data_labels = ['SOSS Order 2', 'SOSS Order 1', 'NIRSpec NRS1', 'NIRSpec NRS2'],
                                  data_colour_list = ['gold', 'lime', 'cyan', 'steelblue'],
                                 add_retrieved_offsets = True, model=model)

#***** Make corner plot *****#

fig_corner = generate_cornerplot(planet, model)




#from POSEIDON.contributions import pressure_contribution, plot_pressure_contribution
#Contribution, norm, \
#        spectrum_contribution_list_names = pressure_contribution(planet, star, model ,atmosphere_contribution_transmission, opac, wl,spectrum_type = 'transmission',bulk_species = False, cloud_contribution = False, total_pressure_contribution = True )

#np.savetxt("/home/gressier/jwst/toi270/toi270b_contribution.txt", np.array(Contribution))


#plot_pressure_contribution(wl, P, planet, Contribution,
#                                   spectrum_contribution_list_names, R = 100,
#                                                              show_log_plot = False, save_fig = True)




