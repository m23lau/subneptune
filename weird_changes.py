# I'll document weird changes that I made here that may or may not be visible

# 1.
# pkg_resources was deprecated after Python 3.12, so I changed a few lines in exotid_ld's ld_grids.py and modelgrid.py in
# Eureka/src/eureka/S5_lightcurve_fitting as follows:

# FROM:     path = pkg_resources.resource_filename('my_package', 'resource.dat')
# TO:       ref = importlib_resources.files('my_package') / 'resource.dat'
#           with importlib_resources.asfile(ref) as path:
#           # and then do something

# -----------------------------------------------------------------------------------------------------------------------

# 2.
# Needed to add these lines in the run_eureka.py files before importing eureka to connect to the server

# import os
# os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
# os.environ["CRDS_PATH"] = "/PATH/TO/crds_cache"
# os.environ["CRDS_MODE"] = "remote"

# -----------------------------------------------------------------------------------------------------------------------
