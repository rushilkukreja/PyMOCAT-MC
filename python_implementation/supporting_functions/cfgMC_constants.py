"""
Configuration constants for MOCAT-MC simulation

Sets up conversion units and global variables used throughout the simulation.
"""

import numpy as np
import sys
import os

# Add supporting functions to path
sys.path.append(os.path.dirname(__file__))

from getgravc import getgravc


def get_cfg_mc_constants():
    """
    Get MOCAT-MC configuration constants
    
    Returns:
        cfgMC: dictionary containing all constants
    """
    
    cfgMC = {}
    
    # Set conversion units
    cfgMC['DAY2MIN'] = 60 * 24
    cfgMC['DAY2SEC'] = cfgMC['DAY2MIN'] * 60
    cfgMC['YEAR2DAY'] = 365.2425  # days(years(1)) 
    cfgMC['YEAR2MIN'] = cfgMC['YEAR2DAY'] * cfgMC['DAY2MIN']
    cfgMC['rad'] = np.pi / 180
    
    # GLOBAL VARIABLES
    whichconst = 84
    tumin, mu_const, radiusearthkm, xke, j2, j3, j4, j3oj2 = getgravc(whichconst)
    
    cfgMC['tumin'] = tumin
    cfgMC['mu_const'] = mu_const
    cfgMC['radiusearthkm'] = radiusearthkm
    cfgMC['j2'] = j2
    cfgMC['omega_earth'] = 2 * np.pi / cfgMC['DAY2SEC']
    
    return cfgMC


# For backward compatibility, create module-level constants
cfgMC = get_cfg_mc_constants()

# Export individual constants for easy access
DAY2MIN = cfgMC['DAY2MIN']
DAY2SEC = cfgMC['DAY2SEC'] 
YEAR2DAY = cfgMC['YEAR2DAY']
YEAR2MIN = cfgMC['YEAR2MIN']
rad = cfgMC['rad']
tumin = cfgMC['tumin']
mu_const = cfgMC['mu_const']
radiusearthkm = cfgMC['radiusearthkm']
j2 = cfgMC['j2']
omega_earth = cfgMC['omega_earth']