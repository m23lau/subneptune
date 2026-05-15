# I'll document weird changes that I made here that may or may not be visible
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
# might need to change a few lines (line 103) in Eureka/src/eureka/S1_detector_processing
#
# FROM:       meta.intstart = hdulist[0].header['INTSTART']-1
#             meta.intend = hdulist[0].header['INTEND']
# TO:         meta.intstart = hdulist[0].header.get('INTSTART', 1) -  1
#             meta.intend = hdulist[0].header.get('INTEND', 1) - 1