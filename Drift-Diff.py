"""
Define Physics Module - Drift-Diffusion.

Input: plasma variables including E-field
Calc diffusion coeff: D = kT/m/v_coll
Calc mobility: Mu = e/m/v_coll
Calc total flux: Flux = -D*dn/dx + q*Mu*E-field
Output: plasma flux
"""

import numpy as np
import matplotlib.pyplot as plt


class DD_Base():
    """Define the base class for Drift-Diffusion and Ambipolar."""

    def __init__(self, geom):
        """Import geometry information."""
        self.geom = geom
        nx = self.geom.nx
        self.fluxe = np.zeros(nx)  # initial eon flux
        self.fluxi = np.zeros(nx)  # initial ion flux

    def __str__(self):
        """Print Drift-Diffusion Approximation."""
        return f'label = {self.fluxe}'

    def init_transp(self, De=5e-1, Di=5e-3, Mue=1.0, Mui=1e-4):
        """
        Initiate diffusion coefficient and mobility.

        initial De=5e-1, Di=5e-3, in m^2/s
        initial Mue=1.0, Mui=1e-4 in (m/s)*(m/V)
        """
        nx = self.geom.nx
        self.De = np.ones(nx)*De  # initial eon diff coeff
        self.Di = np.ones_like(self.De)*Di  # initial ion diff coeff
        self.Mue = np.ones(nx)*Mue  # initial eon mobility
        self.Mui = np.ones_like(self.Mue)*Mui  # initial ion mobility

    def plot_transp(self):
        """Plot potential, E-field."""
        x = self.geom.x
        fig, axes = plt.subplots(1, 2, figsize=(8, 4),
                                 constrained_layout=True)
        # plot potential
        ax = axes[0]
        ax.plot(x, self.De, 'b-')
        ax.plot(x, self.Di, 'r-')
        ax.legend(['e Diffusion Coeff', 'Ion Diffusion Coeff'])
        # plot E-field
        ax = axes[1]
        ax.plot(x, self.Mue, 'b-')
        ax.plot(x, self.Mui, 'r-')
        ax.legend(['e Mobility', 'Ion Mobility'])
        # show fig
        plt.show(fig)

    def calc_diff(self, CCP_1d):
        """Calc diffusion term in flux."""
        dnedx = self.geom.cnt_diff(CCP_1d.ne)
        self.fluxe = self.De*dnedx
        dnidx = self.geom.cnt_diff(CCP_1d.ni)
        self.fluxi = self.Di*dnidx


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
