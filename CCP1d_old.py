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

    def init_plasma(self, ne=1e17, nn=3.3e20, Te=1, Ti=0.1, Se=1e20):
        """
        Initiate plasma attributes.

        At 1 atm, number density = 0.025e27 m^-3.
        At 1 Torr, number density = 3.3e22 m^-3.
        At 1 mTorr, number density = 3.3e19 m^-3.
        """
        nx = self.geom.nx
        self.ne = np.ones(nx)*ne  # initial uniform ne on 1d mesh
        self.ni = np.ones_like(self.ne)*ne  # initial ni to neutralize ne
        self.nn = np.ones_like(self.ne)*nn  # initial neutral density
        self.fluxe = np.zeros_like(self.ne)  # initial eon flux
        self.fluxi = np.zeros_like(self.ne)  # initial ion flux
        self.Te = np.ones_like(self.ne)*Te  # initial eon temperature
        self.Ti = np.ones_like(self.ne)*Ti  # initial ion temperature
        self.Se = np.ones_like(self.ne)*Se  # initial e source term
        self.Si = np.ones_like(self.ne)*Se  # initial ion source term

    def plot_plasma(self):
        """
        Plot plasma variables vs. position x.

        density, flux, temperature
        """
        x = self.geom.x
        fig, axes = plt.subplots(1, 3, figsize=(12, 3),
                                 constrained_layout=True)
        # plot densities
        ax = axes[0]
        ax.plot(x, self.ne, 'b-')
        ax.plot(x, self.ni, 'r-')
        ax.legend(['E', 'Ion'])
        ax.set_xlabel('Position (m)')
        ax.set_ylabel('Density (m^-3)')
        # plot fluxes
        ax = axes[1]
        ax.plot(x, self.fluxe, 'b-')
        ax.plot(x, self.fluxi, 'r-')
        ax.legend(['flux E', 'flux Ion'])
        ax.set_xlabel('Position (m)')
        ax.set_ylabel('Flux (m^-2s^-1)')
        # plot temperature
        ax = axes[2]
        ax.plot(x, self.Te, 'b-')
        ax.plot(x, self.Ti, 'r-')
        ax.legend(['Te', 'Ti'])
        ax.set_xlabel('Position (m)')
        ax.set_ylabel('Temperature (eV)')
        # show fig
        plt.show(fig)

    def init_pot(self, phi=0.0):
        """Initiate potential attributes."""
        nx = self.geom.nx
        self.pot = np.ones(nx)*phi  # initial uniform potential
        self.ef = np.zeros_like(self.pot)  # initial uniform E-field
        self.ef_ambi = np.zeros_like(self.pot)  # initial ambipolar E-field

    def plot_pot(self):
        """Plot potential, E-field."""
        x = self.geom.x
        fig, axes = plt.subplots(1, 2, figsize=(8, 4),
                                 constrained_layout=True)
        # plot potential
        ax = axes[0]
        ax.plot(x, self.pot, 'y-')
        ax.legend(['Potential'])
        # plot E-field
        ax = axes[1]
        ax.plot(x, self.ef, 'g-')
        ax.legend(['E-field'])
        # show fig
        plt.show(fig)

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


if __name__ == '__main__':
    """Test CCP_1d."""
    mesh1d = Mesh_1d('CCP_1d', 10e-2, nx=101)
    print(mesh1d)
    ccp1d = CCP_1d(mesh1d)
    ccp1d.init_plasma()
    ccp1d.plot_plasma()
    ccp1d.init_pot()
    ccp1d.plot_pot()
    ccp1d.init_transp()
    ccp1d.plot_transp()