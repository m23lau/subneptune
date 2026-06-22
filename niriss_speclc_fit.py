import juliet
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pickle
import corner
# import transitspectroscopy as ts
from tqdm import tqdm
import csv

import seaborn as sns
sns.set_style('ticks')

# Open data file - CHANGE
path = '/media/peng/KINGSTON/'
with open(path + '_jw03385004001_toi1130_niriss_lcs.pkl', 'rb') as file:
    data = pickle.load(file)

# Define planetary parameters - CHANGE
per = 4.074554             # Period (Days)
t0 = 2460547.1522035785    # Mid-Transit Time (BJD) for Order 1
# t0 = 2460547.1522049229    # for Order 2

u1 = -0.2394599063         # Limb-darkening coeffs. for Order 1
u2 = 1.0756814820 

# u1 = -0.0476971007         # for Order 2
# u2 = 1.0468093463

a = 13.77                  # Semimajor Axis (a/R*)
b = -0.518                 # Impact Parameter
ecc = 0.052162             # Eccentricity
omega = 141.11             # Argument of Periastron (Degrees)

# Define planet and instrument name - CHANGE
star = 'TOI-1130'
pn = 'TOI-1130b'
itm = 'Order1'

t = data['times']
wl = [i for i in data['order'+itm[-1]]['spectral light curves']]   # Reversed order for some reason but it's ok

# Each data['order#']['spectral light curves'][wl[0]]['flux'] is an array that contains the corresponding spectroscopic light curve in ppm
# Now store all spectroscopics in a dict, where each key is a wavelength whose value is another dict. In the nested dict, keys will be flux and flux error 
# --> spec_lcs = {wavelength 1: {flux: [...], err: [...]}, wavelength 2: {flux: [...], err: [...]}, ..., wavelength n: {flux: [...], err: [...]}}
spec_lcs = {}
for i in range(len(wl)):
    wavelength = wl[i]
    slc = data['order'+itm[-1]]['spectral light curves'][wavelength]['flux']
    slcerr = data['order'+itm[-1]]['spectral light curves'][wavelength]['errors']

    # Normalize each light curve
    slc_med = np.nanmedian(slc[:100])
    spec_lcs[wavelength] = {}
    spec_lcs[wavelength]['f'] = slc / slc_med
    spec_lcs[wavelength]['ferr'] = slcerr / slc_med

# Bin light curves since spec_lcs has spectroscopic light curves in instrument's native resolution
bin_slc = {}
bin_size = 10
for i in range(0, len(wl) - bin_size, bin_size):
    wl_mean = np.round(np.mean(wl[i: i+bin_size]), 7)  # This is just to name the bin
    
    bin_slc[wl_mean] = {}
    bin_slc[wl_mean]['f'] = 0
    bin_slc[wl_mean]['ferr'] = 0

    # Add every bin_size spectroscopic light curves together, then normalize flux to 1
    for j in wl[i: i+bin_size]:
        bin_slc[wl_mean]['f'] += spec_lcs[j]['f']
        bin_slc[wl_mean]['ferr'] += spec_lcs[j]['ferr']

    fbin_med = np.nanmedian(bin_slc[wl_mean]['f'][:100])
    bin_slc[wl_mean]['f'] = bin_slc[wl_mean]['f'] / fbin_med
    bin_slc[wl_mean]['ferr'] = bin_slc[wl_mean]['ferr'] / fbin_med

# Define standardized regressors for linear and/or GP fits
def standardized_regressors(x):
    return (x - np.mean(x)) / np.sqrt(np.var(x))

regressors = np.zeros([len(t), 1]) 
regressors[:, 0] = standardized_regressors(t)


# Open csv to write data - CHANGE 
with open(pn+'_'+itm+'.csv', 'w', newline = '') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['wavelength', 'transit depth', 'transit depth error', 'log evidence'])
    
    for wvlngth in bin_slc:
        # Define data points 
        t = t
        f_s = bin_slc[wvlngth]['f']
        ferr_s = bin_slc[wvlngth]['ferr']
    
        # Bin data (visual purposes only)
        tbin, fsbin, ferrsbin = juliet.bin_data(t, f_s, 25)

        # COMMENT/UNCOMMENT FOR ORDER 1 OR 2
        # Order 1
        params =  ['P_p1', 't0_p1','a_p1','b_p1', 'u1_'+itm, 'u2_'+itm ,'ecc_p1','omega_p1', 'p_p1',
                   'mdilution_'+itm, 'mflux_'+itm, 'sigma_w_'+itm, 
                   'GP_sigma_'+itm, 'GP_rho_'+itm]
    
        dists = ['fixed', 'fixed', 'fixed', 'fixed', 'fixed', 'fixed', 'fixed', 'fixed', 'uniform',
                 'fixed', 'normal', 'loguniform',
                'loguniform', 'loguniform']
    
        hyperps = [per, t0, a, b, u1, u2, ecc, omega, [0., 0.2],
                   1.0, [0., 0.1], [10., 1000.], 
                  [1e-8, 1e2], [1e-8, 1e3]]
        
        # # Order 2
        # params =  ['P_p1', 't0_p1','a_p1','b_p1', 'u1_'+itm, 'u2_'+itm ,'ecc_p1','omega_p1', 'p_p1',
        #            'mdilution_'+itm, 'mflux_'+itm, 'sigma_w_'+itm, 
        #            'theta0_'+itm, 'GP_sigma_'+itm, 'GP_rho_'+itm]
    
        # dists = ['fixed', 'fixed', 'fixed', 'fixed', 'fixed', 'fixed', 'fixed', 'fixed', 'uniform',
        #          'fixed', 'normal', 'loguniform',
        #          'uniform', 'loguniform', 'loguniform']
    
        # hyperps = [per, t0, a, b, u1, u2, ecc, omega, [0., 0.2],
        #            1.0, [0., 0.1], [10., 1000.], 
        #           [-10, 10], [1e-6, 1e6], [1e-3, 1e3]]
        
        priors = juliet.generate_priors(params, dists, hyperps)
        
        # Put data into dictionaries
        times, fluxes, fluxes_err, linear_regressors, GP_regressors = {}, {}, {}, {}, {}
        times[itm], fluxes[itm], fluxes_err[itm] = t, f_s, ferr_s
        linear_regressors[itm] = regressors
        GP_regressors[itm] = regressors
        
        # Load dataset and run sampler - COMMENT/UNCOMMENT FOR ORDERS 1 AND 2
        # Order 1
        dataset = juliet.load(priors = priors, t_lc = times, y_lc = fluxes, yerr_lc = fluxes_err, 
                  GP_regressors_lc = GP_regressors,
                  ld_laws = 'quadratic', out_folder = pn+'_'+str(round(wvlngth, 7))+'_fit')
        
        # # Order 2
        # dataset = juliet.load(priors = priors, t_lc = times, y_lc = fluxes, yerr_lc = fluxes_err, 
        #           linear_regressors_lc = linear_regressors, GP_regressors_lc = GP_regressors,
        #           ld_laws = 'quadratic', out_folder = pn+'_'+str(round(wvlngth, 7))+'_fit')
        
        # results = dataset.fit(sampler = 'dynamic_dynesty', nthreads = 3, verbose = True, n_live_points = 1000)
    
        
        # Take jitter values from posterior and multiply into errors
        sigma = np.median(results.posteriors['posterior_samples']['sigma_w_'+itm] * 1e-6)
        mflux = np.median(results.posteriors['posterior_samples']['mflux_'+itm])
        total_errors = np.sqrt(sigma**2 + ferr_s**2)
    
        # Plot the lightcurve fit and residuals relative to mid-transit
        norm_t = (t - t0) * 24
        norm_tbin = (tbin - t0) * 24
        
        # Scatter data
        plt.errorbar(norm_t, f_s, yerr = total_errors, fmt = '.', color = 'lightgrey', ms = 5, elinewidth = 1, label = 'Data', alpha = 0.5, zorder = 1)
        plt.errorbar(norm_tbin, fsbin, yerr = ferrsbin, fmt = '.', color = 'k', ms = 10, elinewidth = 1, label = 'Binned Data', zorder = 5)
        
        # Get best fit and +/- 1 sigma band from posterior, then plot
        transit_model, transit_up68, transit_low68 = results.lc.evaluate(itm, return_err = True)
        
        plt.plot(norm_t, transit_model, color = 'r', lw = 2, label = 'Fit', zorder = 10)
        plt.fill_between(norm_t, transit_up68, transit_low68, color = 'y', alpha = 0.5, zorder = 8)
        
        plt.title(pn+ ' '+str(round(wvlngth, 7))+ ' $\mu$m Light Curve + Fit')
        plt.xlabel('Time after mid-transit [hours]')
        plt.ylabel('Normalized Flux')
        plt.legend()
        plt.savefig(pn+'_'+str(round(wvlngth, 7))+'_fit/'+str(round(wvlngth, 7))+'_lc_fit.png')
        plt.close()
        
        # Residuals plot
        resid = f_s - transit_model
        resid_binx, resid_biny, _ = juliet.bin_data(norm_t, resid, 20)
        
        plt.scatter(norm_t, resid, alpha = 0.1, zorder = 1)
        plt.scatter(resid_binx, resid_biny, alpha = 0.5, color = 'b', zorder = 5)
        plt.axhline(y = 0, ls = '--', color = 'k', alpha = 0.8, zorder = 10)
        
        plt.title('Residuals')
        plt.xlabel('Time after mid-transit [hours]')
        plt.ylabel('Difference in Normalized Flux')
        plt.savefig(pn+'_'+str(round(wvlngth, 7))+'_fit/'+str(round(wvlngth, 7))+'_residuals.png')
        plt.close()
    
    
        # Corner plot
        list_post = []
        names = []
        for f in results.posteriors['posterior_samples']:
            if f not in ['unnamed', 'loglike']:
                names.append(f)
                list_post.append(results.posteriors['posterior_samples'][f])
        tr = np.array(list_post).T
        corner.corner(tr, labels = names, show_titles=True, quantiles = [0.16, 0.50, 0.84])
        plt.savefig(pn+'_'+str(round(wvlngth, 7))+'_fit/'+str(round(wvlngth, 7))+'_corner.png')
        plt.close()
    
        # Get transit depths
        tds = (results.posteriors['posterior_samples']['p_p1']**2) * 100
    
        # Write wavelength, median transit depth, transit depth error and log evidence into csv file
        writer.writerow([wvlngth, np.nanmedian(tds), np.sqrt(np.var(tds)), results.posteriors['lnZ']])

