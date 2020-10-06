"""
Define CCP_1d.

CCP_1d plasma model.
CCP_1d contains:
    density, flux, potential, E-field
"""

from Mesh import Mesh_1d

import numpy as np
import matplotlib.pyplot as plt


class CCP_1d(object):
    """Define 1d CCP."""

    def __init__(self, geom):
        """
        CCP_1d is defined as a container.

        CCP_1d as a basket containing:
            geometry
            physics.
        """
        self.geom = geom

    def __str__(self):
        """Print 1d mesh information."""
        res = 'CCP_1d:'
        return res

    def init_plasma(self, ne=1e17, nAr=3.3e20, Te=1, Ti=0.1, Se=1e20):
        """
        Initiate plasma attributes.

        At 1 atm, nAr = 0.025e27 m^-3.
        At 1 Torr, nAr = 3.3e22 m^-3.
        At 1 mTorr, nAr = 3.3e19 m^-3.
        """
        nx = self.geom.nx
        self.ne = np.ones(nx)*ne  # initial uniform ne on 1d mesh
        self.ni = np.ones_like(self.ne)*ne  # initial ni to neutralize ne
        self.nAr = np.ones_like(self.ne)*nAr  # initial Ar density
        self.fluxe = np.zeros_like(self.ne)  # initial eon flux
        self.fluxi = np.zeros_like(self.ne)  # initial ion flux
        self.Te = np.ones_like(self.ne)*Te  # initial eon temperature
        self.Ti = np.ones_like(self.ne)*Ti  # initial ion temperature
        self.Se = np.ones_like(self.ne)*Se  # initial e source term
        self.Si = np.ones_like(self.ne)*Se  # initial ion source term

    def plot_plasma(self):
        """Plot 1d mesh in X."""
        fig, axes = plt.subplots(1, 2, figsize=(8, 4),
                                 constrained_layout=True)
        # plot densities
        ax = axes[0]
        ax.plot(self.x, self.ne, 'b-')
        ax.plot(self.x, self.ni, 'r-')
        # add legend
        ax.legend(['E', 'Ion'])
        # plot fluxes
        ax = axes[1]
        ax.plot(self.x, self.fluxe, 'b-')
        ax.plot(self.x, self.fluxi, 'r-')
        # add legend
        ax.legend(['flux E', 'flux Ion'])
        # show fig
        plt.show(fig)

    def add_bndy(self):
        """Add boundaries."""
        pass

    def add_mat(self):
        """Add materials."""
        pass

    def cnt_diff(self, y):
        """
        Caculate dy/dx using central differencing.

        Using x[i+1] - x[i-1] instead of delx makes cnt_diff() compatible
        with non-uniform 1d mesh.
        """
        dy = np.zeros_like(self.x)
        # Although dy[0] and dy[-1] are signed here,
        # they are eventually specified in boundary conditions
        dy[0] = (y[1] - y[0])/(self.x[1] - self.x[0])
        dy[-1] = (y[-1] - y[-2])/(self.x[-1] - self.x[-2])
        for i in range(1, self.nx-1):
            dy[i] = (y[i+1] - y[i-1])/(self.x[i+1] - self.x[i-1])
        return dy


if __name__ == '__main__':
    """Test Mesh."""
    mesh1d = Mesh_1d('CCP_1d', 10e-2, nx=101)
    print(mesh1d)
