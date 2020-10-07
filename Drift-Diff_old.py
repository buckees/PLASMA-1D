"""
Define Physics Module - Drift-Diffusion.

Input: plasma variables including E-field
Calc diffusion coeff: D = kT/m/v_coll
Calc mobility: Mu = e/m/v_coll
Calc total flux: Flux = -D*dn/dx + q*Mu*E-field
Output: plasma flux
"""

import numpy as np


class Drift_Diff(object):
    """Define Drift-Diffusion Physics."""

    def __init__(self, CCP_1d):
        self.fluxe = CCP_1d.fluxe
        self.fluxi = CCP_1d.fluxi
        # variables below are not outputs
        self.De = CCP_1d.De
        self.Di = CCP_1d.Di

    def __str__(self):
        """Print Drift-Diffusion Approximation."""
        return f'label = {self.fluxe}'

    def calc_transp(self, CCP1d):
        """
        Calc transport coeff.

        Calc diffusion coeff
        Calc mobility
        output: diffusion coeff and mobility
        """
        # self.De =
        # self.Di =
        return self.De, self.Di

    def calc_flux(self, CCP_1d):
        """
        Calc plasma flux.

        The drift-diffusion approximation:
        """
        self.Da = CCP_1d.Di*(1 + np.divide(CCP_1d.Te, CCP_1d.Ti))
        self.Ea = (CCP_1d.Di - CCP_1d.De)/(CCP_1d.Mui + CCP_1d.Mue)
        dnidx = CCP_1d.geom.cnt_diff(CCP_1d.ni)
        self.Ea *= np.divide(dnidx, CCP_1d.ni,
                             out=np.zeros_like(dnidx), where=CCP_1d.ni != 0)
        self.Ea[1] = self.Ea[2]
        self.Ea[-2] = self.Ea[-3]
        return self.Da, self.Ea
