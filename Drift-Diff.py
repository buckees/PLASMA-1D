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
import copy

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

    def calc_transp(self, CCP1d):
        """
        Calc transport coeff.

        Calc diffusion coeff
        Calc mobility
        output: diffusion coeff and mobility
        """
        self.init_transp()
        # self.De =
        # self.Di =
        pass

    def calc_diff(self, CCP_1d):
        """Calc diffusion term in flux."""
        dnedx = self.geom.cnt_diff(CCP_1d.ne)
        self.fluxe = self.De*dnedx
        dnidx = self.geom.cnt_diff(CCP_1d.ni)
        self.fluxi = self.Di*dnidx
        # self.bndy_flux()

    def bndy_flux(self):
        """Impose b.c. to flux."""
        self.fluxe[1], self.fluxe[-2] = self.fluxe[2], self.fluxe[-3]
        self.fluxi[1], self.fluxi[-2] = self.fluxi[2], self.fluxi[-3]


class Ambipolar(DD_Base):
    """
    Define Ambipolar Physics.

    Input: Geometry, Te/Ti, ne/ni
    Output: Flux, E-field
    """

    def calc_ambi(self, CCP_1d):
        """
        Calc ambipolar diffusion coefficient.

        The ambipolar diffusion assumptions:
            1. steady state, dne/dt = 0. it cannot be used to
            describe plasma decay.
            2. ni is calculated from continuity equation.
            3. plasma is charge neutral, ne = ni
            4. Ionization Se is needed to balance diffusion loss.
        Da = (De*Mui + Di*Mue)/(Mue + Mui)
        Da = Di(1 + Te/Ti).
        Ea = (Di - De)/(Mui + Mue)*dn/dx/n
        Orginal Ambipolar Coeff Da = (De*Mui + Di*Mue)/(Mue + Mui)
        self.Da = (CCP_1d.De*CCP_1d.Mui + CCP_1d.Di*CCP_1d.Mue) / \
                  (CCP_1d.Mue + CCP_1d.Mui)
        Assume Te >> Ti, Ambipolar Coeff can be simplified as
        Da = Di(1 + Te/Ti).
        """
        self.calc_transp(CCP_1d)
        self.Da = self.Di*(1 + CCP_1d.Te / CCP_1d.Ti)
        self.Ea = (self.Di - self.De)/(self.Mui + self.Mue)
        dnidx = self.geom.cnt_diff(CCP_1d.ni)
        self.Ea *= np.divide(dnidx, CCP_1d.ni,
                             out=np.zeros_like(dnidx), where=CCP_1d.ni != 0)
        self.bndy_ambi()

    def bndy_ambi(self):
        """Impose b.c. to Ea."""
        self.Ea[0], self.Ea[-1] = 0.0, 0.0
        self.Ea[1], self.Ea[-2] = self.Ea[2], self.Ea[-3]
        self.Da[0], self.Da[-1] = 0.0, 0.0

    def calc_flux(self, CCP_1d):
        """Calc Ambipolar flux."""
        self.calc_ambi(CCP_1d)
        self.De, self.Di = self.Da, self.Da
        self.calc_diff(CCP_1d)
        self.fluxe = copy.deepcopy(self.fluxi)
        return self.fluxe, self.fluxi

    def calc_ne(self, CCP_1d):
        """Calc ne = sum(ni), ensure charge neutrality."""
        return copy.deepcopy(CCP_1d.ni)


class Drift_Diff(DD_Base):
    """Define Drift-Diffusion Physics."""

    def calc_flux(self, CCP_1d):
        """
        Calc plasma flux.

        The drift-diffusion approximation:
        """
        pass


if __name__ == '__main__':
    """Test the Ambipolar."""
    from Mesh import Mesh_1d
    from CCP1d import CCP_1d
    mesh1d = Mesh_1d('CCP_1d', 10e-2, nx=101)
    print(mesh1d)
    ccp1d = CCP_1d(mesh1d)
    ccp1d.init_plasma()
    # ccp1d.plot_plasma()
    ccp1d.init_pot()
    # ccp1d.plot_pot()
    ambi = Ambipolar(mesh1d)
    # ambi.plot_transp()
    ccp1d.fluxe, ccp1d.fluxi = ambi.calc_flux(ccp1d)
    print(ambi.__dict__)
