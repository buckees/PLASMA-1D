"""
1D Plasma Electron Energy Module

Eergy_1d contains:
    Electron energy equation
    d(3/2nekTe)/dt = -dQ/dx + Power_in(ext.) - Power_loss(react)
    Input: ne, Te from Plasma1d, E_ext from field solver
    Output: Te
"""

from Constants import KB_EV, EON_MASS, UNIT_CHARGE

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


class Eergy_1d(object):
    """Define the eon energy module/object."""
    
    def __init__(self, pla):
        """Import Plasma1d information."""
        nx = pla.geom.nx
        self.Qe = np.zeros(nx)  # initial eon flux
        self.dQe = np.zeros(nx)  # initial eon flux
        self.Te = deepcopy(pla.Te)
        # eon energy = 3/2 * ne * kTe
        self.ergy_e = 1.5*KB_EV*np.multiply(pla.ne, pla.Te)
        
        
    def __str__(self):
        """Print eon energy module."""
        return f'label = {self.qdfluxe}'
    
    def calc_th_cond_coeff(self, pla):
        """
        Calc thermal conduction coefficient.

        pla: Plasma_1d object
             calc uses pla.Te,i and pla.coll_em
        heat_cond_e: W/m/K, heat conductivity for eon
        """
        # calc thermal conductivity for eon
        self.th_cond_e = np.ones(pla.ne)*1.0
    
    def calc_th_flux(self, pla, txp):
        """
        Calc eon thermal flux, Qe
        
        Qe = 5/2kTe * fluxe - ke * dTe/dx
        dQe = 5/2kTe * dfluxe - ke * d2Te/dx2
        """
        # calc convection term
        self.Qe = 2.5*KB_EV*np.multiply(self.Te, txp.fluxe)
        self.dQe = 2.5*KB_EV*np.multiply(self.Te, txp.dfluxe)
        # calc conduction term
        self.dTe = pla.goem.cnt_diff(self.Te)
        self.d2Te = pla.goem.cnt_diff_2nd(self.Te)
        self.Qe -= np.multiply(self.th_cond_e, self.dTe)
        self.dQe -= np.multiply(self.th_cond_e, self.d2Te)
        
    def calc_Te(self, delt, pla, pwr):
        """Calc Te"""
        self.ergy_e += (-self.dQe + pwr.input)*delt
        self.Te = np.divide(self.ergy_e, self.Te)/1.5/KB_EV

