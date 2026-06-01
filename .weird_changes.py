# I've documented weird changes to Eureka! that I made here when I was playing around with it
# I probably should do this in any file that's not a .py file but too bad

# 1.
# pkg_resources was deprecated after Python 3.12, so I changed a few lines in exotid_ld's ld_grids.py and modelgrid.py
# in Eureka/src/eureka/S5_lightcurve_fitting (lines 45 and 243) as follows:
#
# FROM:     path = pkg_resources.resource_filename('my_package', 'resource.dat')
# TO:       ref = importlib_resources.files('my_package') / 'resource.dat'
#           with importlib_resources.asfile(ref) as path:
#           # and then do something

# -----------------------------------------------------------------------------------------------------------------------

# 2.
# Needed to add these lines in the run_eureka.py files before importing eureka to connect to the server
#
# import os
# os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
# os.environ["CRDS_PATH"] = "/PATH/TO/crds_cache"
# os.environ["CRDS_MODE"] = "remote"

# -----------------------------------------------------------------------------------------------------------------------

# 3. OPTIONAL
# Small problem with simulated NIRISS demo data where I get a KeywordError - missing INTSTART and INTEND in header
# might need to change a few lines (line 103) in Eureka/src/eureka/S1_detector_processing. This isn't an issue with
# downloaded data sets from the notebook.
#
# FROM:       meta.intstart = hdulist[0].header['INTSTART']-1
#             meta.intend = hdulist[0].header['INTEND']
# TO:         meta.intstart = hdulist[0].header.get('INTSTART', 1) -  1
#             meta.intend = hdulist[0].header.get('INTEND', 1) - 1

# -----------------------------------------------------------------------------------------------------------------------

# 4.
# Likely an issue with something being out of date again, in Eureka/src/eureka/S3_data_reduction, needed to change
# background.py (line 39) to not include "skip_bg" since S1_detector_processing/s1_meta.py does not give the S1metaclass
# object an attribute called S1metaclass.skip_bg
#
# FROM:    if meta.bg_deg is None or meta.skip_bg:
# TO:      if meta.bg_deg is None:

# -----------------------------------------------------------------------------------------------------------------------
# ok i think i'm low key going to give up trying to get eureka to work but i guess i'll at least keep this documentation
# if i want to do things in the future with it, just not now since i don't think i have enough time to get it running
# i'll upload this final version
