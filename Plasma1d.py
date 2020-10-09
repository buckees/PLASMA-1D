"""
Define Plasma_1d.

Plasma_1d plasma model.
Plasma_1d contains:
    density, flux, potential, E-field
"""

from Mesh import Mesh_1d

import numpy as np
import matplotlib.pyplot as plt


class Plasma_1d(object):
    """Define 1d Plasma."""

    def __init__(self, geom):
        """
        Plasma_1d is defined as a container.

        Plasma_1d as a basket containing:
            geometry
            physics.
        """
        self.geom = geom

    def __str__(self):
        """Print 1d mesh information."""
        res = 'Plasma_1d:'
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
        self.bndy_plasma()
        self.limit_plasma()

    def bndy_plasma(self):
        """Impose b.c. on the plasma."""
        self.ne[0], self.ne[-1] = 0.0, 0.0
        self.ni[0], self.ni[-1] = 0.0, 0.0
        self.nn[0], self.nn[-1] = 0.0, 0.0
        self.Te[0], self.Te[-1] = 0.0, 0.0
        self.Ti[0], self.Ti[-1] = 0.0, 0.0
        self.Se[0], self.Se[-1] = 0.0, 0.0
        self.Si[0], self.Si[-1] = 0.0, 0.0

    def limit_plasma(self, n_min=1e11, n_max=1e22, T_min=0.001, T_max=100.0):
        """Limit variables in the plasma."""
        self.ne = np.clip(self.ne, n_min, n_max)
        self.ni = np.clip(self.ni, n_min, n_max)
        self.nn = np.clip(self.nn, n_min, n_max)
        self.Te = np.clip(self.Te, T_min, T_max)
        self.Ti = np.clip(self.Ti, T_min, T_max)

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

    def den_evolve(self, delt):
        """
        Evolve the density in Plasma by solving the continuity equation.

        dn/dt = -dFlux/dx + Se
        dn(t + dt) = dn(t) - dFlux/dx*dt + Se*dt
        """
        dfluxe = self.geom.cnt_diff(self.fluxe)
        self.ne -= dfluxe*delt
        dfluxi = self.geom.cnt_diff(self.fluxi)
        self.ni -= dfluxi*delt


if __name__ == '__main__':
    """Test Plasma_1d."""
    from Drift_Diff import Ambipolar, Diffusion
    mesh1d = Mesh_1d('Plasma_1d', 10e-2, nx=51)
    print(mesh1d)
    Plasma1d = Plasma_1d(mesh1d)
    Plasma1d.init_plasma()
    Plasma1d.plot_plasma()
    Plasma1d.init_pot()
    # Plasma1d.plot_pot()
    ambi = Ambipolar(mesh1d)
    # diff = Diffusion(mesh1d)
    ne_ave, ni_ave = [], []
    time = []
    dt = 1e-6
    niter = 10000
    for itn in range(niter):
        Plasma1d.fluxe, Plasma1d.fluxi = ambi.calc_flux(Plasma1d)
        # Plasma1d.fluxe, Plasma1d.fluxi = diff.calc_flux(Plasma1d)
        Plasma1d.den_evolve(dt)
        Plasma1d.bndy_plasma()
        Plasma1d.limit_plasma()
        ne_ave.append(np.mean(Plasma1d.ne))
        ni_ave.append(np.mean(Plasma1d.ni))
        time.append(dt*(niter+1))
        if not (itn+1) % (niter/10):
            Plasma1d.plot_plasma()
