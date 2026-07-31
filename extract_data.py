# Call these functions in fitting notebooks/files to automatically extract data from transitspectroscopy/eureka outputs

import numpy as np
import matplotlib.pyplot as plt
import pickle
import h5py

# Transitspectroscopy - NIRSpec

# with open('file.pkl', 'rb') as file:
#     data = pickle.load(file)

# data is a dict with 7 entries: 
# 1. metadata: Information about instrument and detector, ex. name: NIRSpec, detector: NRS1
# 2. & 3. tso & tso_err: Not sure, just a bunch of zeroes as far as I can see
# 4. traces: Contains times in random order, x, ycorrected, ysmoothed
# 5 spectra: Contains times in increasing order, original + original_err, corrected + corrected_err, wavelength map (a bunch of nans) and wavelengths that the detector measures
# --> IMPORTANT: data['spectra']['corrected'] and corrected_err contain each spectroscopic lc for the corresponding wavelength in ppm (more detail in next section)
# 6 & 7. whitelight & whitelight_err: Relative flux for white light curve and its errorbars

def tspec_nirspec(path, bin_size):
    """ Extract data from a transitspectroscopy output file for NIRSpec
    Args:
        path (str): Path to file
        bin_size (int): Resolution to bin spectroscopic light curves down to 
    Returns:
        tuple: ndarray, ndarray, ndarray, ndarray, dict, dict
               times, white light flux, wl err, wavelengths, spectroscopic lcs and binned slcs
    """
    with open(path, 'rb') as file:
        data = pickle.load(file)

    spectra = data['spectra']
    t = spectra['times']

    # White light curves
    f_wl = data['whitelight']
    ferr_wl = data['whitelight_err']
    
    # Normalize white light curves
    f_med = np.nanmedian(f_wl[:100])
    f_wl, ferr_wl = f_wl / f_med, ferr_wl / f_med
    
    # Spectroscopic light curves
    wl = spectra['wavelengths']
    slc, slcerr = spectra['corrected'], spectra['corrected_err']

    spec_lcs = {}
    for i in range(len(wl)):
        wavelength = wl[i]
    
        # Normalize each light curve 
        slc_med = np.nanmedian(slc[:, i][:100])
        spec_lcs[wavelength] = {}
        spec_lcs[wavelength]['f'] = slc[:, i] / slc_med
        spec_lcs[wavelength]['ferr'] = slcerr[:, i] / slc_med

    # Binned spectroscopic light curves
    bin_slc = {}
    for i in range(0, len(wl) - bin_size, bin_size):
        wl_mean = np.round(np.mean(wl[i: i+bin_size]), 7)
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

    return t, f_wl, ferr_wl, wl, spec_lcs, bin_slc

# -----------------------------------------------------------------------

# Transitspectroscopy - NIRISS

# with open('file.pkl', 'rb') as file:
#     data = pickle.load(file)
    
# data is a dict with 3 entries:
# 1. order1: Contains spectral light curves and white light for 0.9 to 2.8 micron
# --> each one of these keys contains more dicts, for example, data['order1']['spectral light curves'] shows each {wavelength: {flux: [...], errors: [...]}} or just flux and errors for wl
# 2. order2: same as order 1 but for 0.6 to 1.4 micron
# 3. times: Just an array with times

def tspec_niriss(path, bin_size, order):
    """ Extract data from a transitspectroscopy output file for NIRISS
    Args:
        path (str): Path to file
        bin_size (int): Resolution to bin spectroscopic light curves down to 
        order (int): NIRISS SOSS order number
    Returns:
        tuple: ndarray, ndarray, ndarray, ndarray, dict, dict
               times, white light flux, wl err, wavelengths, spectroscopic lcs and binned slcs
    """
    with open(path, 'rb') as file:
        data = pickle.load(file)

    t = data['times']

    # White light curves
    f_wl = data['order'+str(order)]['white light']['flux']
    ferr_wl = data['order'+str(order)]['white light']['errors']

    # Normalize white light curves
    f_med = np.nanmedian(f_wl[:100])
    f_wl, ferr_wl = f_wl / f_med, ferr_wl / f_med
    
    # Spectroscopic light curves
    wl = [i for i in data['order'+str(order)]['spectral light curves']]

    spec_lcs = {}
    for i in range(len(wl)):
        wavelength = wl[i]
        slc = data['order'+str(order)]['spectral light curves'][wavelength]['flux']
        slcerr = data['order'+str(order)]['spectral light curves'][wavelength]['errors']
    
        # Normalize each light curve
        slc_med = np.nanmedian(slc[:100])
        spec_lcs[wavelength] = {}
        spec_lcs[wavelength]['f'] = slc / slc_med
        spec_lcs[wavelength]['ferr'] = slcerr / slc_med

    # Binned spectroscopic light curves
    bin_slc = {}
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
    
    return t, f_wl, ferr_wl, wl, spec_lcs, bin_slc

# -----------------------------------------------------------------------

# Eureka - S4 is the same for all instruments so the output .h5 file is the same

# data = h5py.File('file.h5', 'r')

# Calling data.keys() will list all keys from the file, but the important ones are
# wavelength: Wavelengths from S3
# time: Times of observation
# data: An mxn array containing spectroscopic light curves, with m wavelengths and n fluxes that correspond to time
# err: A similar array as flux but for slc errors
# flux_white: White light fluxes, just a 1D array with same length as time
# err_white: White light flux errors

# Since you can bin the data within Eureka's S4 ecf, no need to create bin_slc here

def eureka_dat(path):
    """ Extract data from a Eureka's Stage 4 output file
    Args:
        path (str): Path to file
    Returns:
        tuple: ndarray, ndarray, ndarray, ndarray, dict, dict
               times, white light flux, wl err, wavelengths and spectroscopic lcs
    """
    data = h5py.File(path)
    t = np.array(data['time'])

    # White light curves
    f_wl = np.array(data['flux_white'])         
    ferr_wl = np.array(data['err_white'])

     # Normalize white light curves
    f_med = np.nanmedian(f_wl[:100])
    f_wl, ferr_wl = f_wl / f_med, ferr_wl / f_med
    
    # Spectroscopic light curves
    wl = np.array(data['wavelength'])
    slc = np.array(data['data'])
    slcerr = np.array(data['err'])

    spec_lcs = {}
    for i in range(len(wl)):
        wavelength = wl[i]

        # Normalize each light curve
        slc_med = np.nanmedian(slc[i][:100])
        spec_lcs[wavelength] = {}
        spec_lcs[wavelength]['f'] = slc[i] / slc_med
        spec_lcs[wavelength]['ferr'] = slcerr[i] / slc_med

    return t, f_wl, ferr_wl, wl, spec_lcs