"""
1D Plasma Transport Module

Transp_1d contains:
    Diffusion only
    Ambipolar
    Drfit-Diffusion
    Momentum Solver

    Continuity Eq. dn/dt = -dF/dx * dt + S * det
    Input: depends on transport mode
    Output: dF/dx for continuity equation.
"""

from Constants import KB_EV, EON_MASS, UNIT_CHARGE

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


class React_1d(object):
    """Define the base tranport module/object."""
    
    def __init__(self, pla):
        """Import geometry information."""
        nx = pla.geom.nx
        self.se = np.zeros(nx)  # initial eon flux
        self.si = np.zeros(nx)  # initial ion flux
        
    def __str__(self):
        """Print Transport Module."""
        return f'label = {self.dfluxe}'