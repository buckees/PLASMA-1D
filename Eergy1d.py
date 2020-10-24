"""
1D Plasma Electron Energy Module

Eergy_1d contains:
    Electron energy equation
    d(3/2nekTe) = -dQ/dx + Power_in(ext.) - Power_loss(react)
    Input: ne, Te from Plasma1d, E_ext from field solver
    Output: Te
"""

from Constants import KB_EV, EON_MASS, UNIT_CHARGE

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


class Transp_1d(object):
    """Define the eon energy module/object."""
    
    def __init__(self, pla):
        """Import Plasma1d information."""
        nx = pla.geom.nx
        self.qfluxe = np.zeros(nx)  # initial eon flux
        self.qdfluxe = np.zeros(nx)  # initial eon flux       
        
    def __str__(self):
        """Print eon energy module."""
        return f'label = {self.qdfluxe}'
