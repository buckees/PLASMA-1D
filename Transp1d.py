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


class Transp_1d(object):
    """Define the base tranport module/object."""
    
    def __init__(self, pla):
        """Import geometry information."""
        nx = pla.geom.nx
        self.fluxe = np.zeros(nx)  # initial eon flux
        self.fluxi = np.zeros(nx)  # initial ion flux
        self.dfluxe = np.zeros(nx)  # initial eon flux
        self.dfluxi = np.zeros(nx)  # initial ion flux
        
    def __str__(self):
        """Print Transport Module."""
        return f'label = {self.dfluxe}'

    # De/Di, Mue/Mui, initial values are not corret
    def calc_transp_coeff(self, pla):
        """
        Initiate diffusion coefficient and mobility.

        pla: Plasma_1d object
             calc uses pla.Te,i and pla.coll_em
        De,i: m^2/s, D = k*T/(m*coll_m)
        Mue,i: m^2/(V*s), Mu = q/(m*coll_m)
        """
        # calc diff coeff: D = k*T/(m*coll_m)
        self.De = np.divide(KB_EV*pla.Te, EON_MASS*pla.coll_em)  
        self.Di = np.divide(KB_EV*pla.Ti, pla.Mi*pla.coll_im)  
        # calc mobility: Mu = q/(m*coll_m)
        self.Mue = UNIT_CHARGE/EON_MASS/pla.coll_em
        self.Mui = UNIT_CHARGE/pla.Mi/pla.coll_im
        
    # def bndy_transp(self):
    #     """
    #     Impose b.c. on transport coeff.

    #     Extension b.c.
    #     """
    #     # self.De[0], self.De[-1] = 0.0, 0.0
    #     # self.Di[0], self.Di[-1] = 0.0, 0.0
    #     # self.Mue[0], self.Mue[-1] = 0.0, 0.0
    #     # self.Mui[0], self.Mui[-1] = 0.0, 0.0
    #     self.De[0], self.De[-1] = self.De[1], self.De[-2]
    #     self.Di[0], self.Di[-1] = self.Di[1], self.Di[-2]
    #     self.Mue[0], self.Mue[-1] = self.Mue[1], self.Mue[-2]
    #     self.Mui[0], self.Mui[-1] = self.Mui[1], self.Mui[-2]

    # def limit_transp(self, D_min=1e-6, D_max=1e3, M_min=1e-7, M_max=1e3):
    #     """Limit variables in the plasma."""
    #     self.De = np.clip(self.De, D_min, D_max)
    #     self.Di = np.clip(self.Di, D_min, D_max)
    #     self.Mue = np.clip(self.Mue, M_min, M_max)
    #     self.Mui = np.clip(self.Mui, M_min, M_max)

    def plot_transp_coeff(self, pla):
        """
        Plot transp coeff.
        
        pla: Plasma_1d object
            use pla.geom.x for plot
        """
        x = pla.geom.x
        fig, axes = plt.subplots(1, 2, figsize=(8, 4),
                                 constrained_layout=True)
        # plot potential
        ax = axes[0]
        ax.plot(x, self.De, 'bo-')
        ax.plot(x, self.Di, 'ro-')
        ax.legend(['e Diffusion Coeff', 'Ion Diffusion Coeff'])
        # plot E-field
        ax = axes[1]
        ax.plot(x, self.Mue, 'bo-')
        ax.plot(x, self.Mui, 'ro-')
        ax.legend(['e Mobility', 'Ion Mobility'])
        # show fig
        plt.show(fig)
    
    def plot_flux(self, pla):
        """
        Plot flux and dflux.
        
        pla: Plasma_1d object
            use pla.geom.x for plot
        """
        x = pla.geom.x
        fig, axes = plt.subplots(1, 2, figsize=(8, 4),
                                 constrained_layout=True)
        # plot potential
        ax = axes[0]
        ax.plot(x, self.fluxe, 'bo-')
        ax.plot(x, self.fluxi, 'ro-')
        ax.legend(['e flux', 'Ion flux'])
        # plot E-field
        ax = axes[1]
        ax.plot(x, self.dfluxe, 'bo-')
        ax.plot(x, self.dfluxi, 'ro-')
        ax.legend(['e dflux', 'Ion dflux'])
        # show fig
        plt.show(fig)

    # def calc_diff(self, Plasma_1d):
    #     """Calc diffusion term in flux."""
    #     dnedx = self.geom.cnt_diff(Plasma_1d.ne)
    #     self.fluxe = -self.De*dnedx
    #     dnidx = self.geom.cnt_diff(Plasma_1d.ni)
    #     self.fluxi = -self.Di*dnidx
    #     # self.bndy_flux()

    # def bndy_flux(self):
    #     """Impose b.c. to flux."""
    #     self.fluxe[1], self.fluxe[-2] = 2*self.fluxe[2], 2*self.fluxe[-3]
    #     self.fluxi[1], self.fluxi[-2] = 2*self.fluxi[2], 2*self.fluxi[-3]

class Diff_1d(Transp_1d):
    """
    Calc the dflux for Diffusion Only Module.
    
    dn/dt = -D * d2n/dx2 + Se
    D: m^2/s, diffusion coefficient is calc from Tranps_1d
    Output: D * d2n/dx2
    """
    def calc_diff(self, pla):
        """Calc diffusion term: D * d2n/dx2 and diffusion flux D * dn/dx. """
        # Calc 
        self.fluxe = -self.De * pla.geom.cnt_diff(pla.ne)
        self.fluxi = -self.Di * pla.geom.cnt_diff(pla.ni)
        
        self.dfluxe = -self.De * pla.geom.cnt_diff_2nd(pla.ne)
        self.dfluxi = -self.Di * pla.geom.cnt_diff_2nd(pla.ni)
    
# class Ambipolar(DD_Base):
#     """
#     Define Ambipolar Physics.

#     Input: Geometry, Te/Ti, ne/ni
#     Output: Flux, E-field
#     """

#     def calc_ambi(self, Plasma_1d):
#         """
#         Calc ambipolar diffusion coefficient.

#         The ambipolar diffusion assumptions:
#             1. steady state, dne/dt = 0. it cannot be used to
#             describe plasma decay.
#             2. ni is calculated from continuity equation.
#             3. plasma is charge neutral, ne = ni
#             4. Ionization Se is needed to balance diffusion loss.
#         Da = (De*Mui + Di*Mue)/(Mue + Mui)
#         Da = Di(1 + Te/Ti).
#         Ea = (Di - De)/(Mui + Mue)*dn/dx/n
#         Orginal Ambipolar Coeff Da = (De*Mui + Di*Mue)/(Mue + Mui)
#         self.Da = (Plasma_1d.De*Plasma_1d.Mui + Plasma_1d.Di*Plasma_1d.Mue) / \
#                   (Plasma_1d.Mue + Plasma_1d.Mui)
#         Assume Te >> Ti, Ambipolar Coeff can be simplified as
#         Da = Di(1 + Te/Ti).
#         """
#         self.calc_transp(Plasma_1d)
#         self.Da = self.Di*(1 + Plasma_1d.Te / Plasma_1d.Ti)
#         self.Ea = (self.Di - self.De)/(self.Mui + self.Mue)
#         dnidx = self.geom.cnt_diff(Plasma_1d.ni)
#         self.Ea *= np.divide(dnidx, Plasma_1d.ni,
#                              out=np.zeros_like(dnidx), where=Plasma_1d.ni != 0)
#         self.bndy_ambi()

#     def bndy_ambi(self):
#         """
#         Impose b.c. to Ea.

#         Extension b.c.
#         """
#         self.Ea[0], self.Ea[-1] = self.Ea[1], self.Ea[-2]
#         # self.Ea[1], self.Ea[-2] = self.Ea[2], self.Ea[-3]
#         self.Da[0], self.Da[-1] = self.Da[1], self.Da[-2]

#     def calc_flux(self, Plasma_1d):
#         """Calc Ambipolar flux."""
#         self.calc_ambi(Plasma_1d)
#         self.De, self.Di = self.Da, self.Da
#         self.calc_diff(Plasma_1d)
#         self.fluxe = copy.deepcopy(self.fluxi)
#         return self.fluxe, self.fluxi

#     def calc_ne(self, Plasma_1d):
#         """Calc ne = sum(ni), ensure charge neutrality."""
#         return copy.deepcopy(Plasma_1d.ni)


# class Diffusion(DD_Base):
#     """
#     Diffusion only.

#     No interactions between electrons and ions.
#     They diffuse independently.
#     """

#     def calc_flux(self, Plasma_1d):
#         """Calc diffusion coeff."""
#         self.calc_transp(Plasma_1d)
#         self.calc_diff(Plasma_1d)
#         return self.fluxe, self.fluxi


# class Drift_Diff(DD_Base):
#     pass
#     """Define Drift-Diffusion Physics."""

#     def calc_flux(self, Plasma_1d):
#         """
#         Calc plasma flux.

#         The drift-diffusion approximation:
#         """
#         pass


if __name__ == '__main__':
    """Test the tranp coeff calc."""
    from Mesh import Mesh_1d
    from Plasma1d import Plasma_1d
    mesh1d = Mesh_1d('Plasma_1d', 10e-2, nx=11)
    print(mesh1d)
    plasma1d = Plasma_1d(mesh1d)
    plasma1d.init_plasma()
    # Plasma1d.plot_plasma()
    txp1d = Diff_1d(plasma1d)
    txp1d.calc_transp_coeff(plasma1d)
    txp1d.plot_transp_coeff(plasma1d)
    txp1d.calc_diff(plasma1d)
    txp1d.plot_flux(plasma1d)
    
