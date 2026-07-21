# Call these functions in fitting notebooks/files to automatically extract data from transitspectroscopy/eureka outputs

import numpy as np
import matplotlib.pyplot as plt
import pickle
import h5py


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



def eureka_dat(path, bin_size):
    """ Extract data from a Eureka's Stage 4 output file
    Args:
        path (str): Path to file
        bin_size (int): Resolution to bin spectroscopic light curves down to 
    Returns:
        tuple: ndarray, ndarray, ndarray, ndarray, dict, dict
               times, white light flux, wl err, wavelengths, spectroscopic lcs and binned slcs
    """
    data = h5py.File(path)
    t = np.array(data['time'])

    # White light curves
    f_wl = np.array(data['flux_white'])         
    ferr_wl = np.array(data['err_white'])

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