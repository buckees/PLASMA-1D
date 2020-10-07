"""
Define Physics Module Ambipolar.

Calc Ambipolar diffusion coeff
Calc Ambipolar E-field
"""

import numpy as np


class Ambipolar(object):
    """Define Ambipolar Physics."""

    def __init__(self, CCP_1d):
        self.Da = CCP_1d.Di
        self.Ea = CCP_1d.ef

    def __str__(self):
        """Print Ambipolar Physics."""
        return f'label = {self.Da}'

    def calc_ambi(self, CCP_1d):
        """
        Calc ambipolar diffusion coefficient.

        The ambipolar diffusion assumptions:
            1. steady state, dne/dt = 0. it cannot be used to
            describe plasma decay.
            2. ni is calculated from continuity equation.
            3. plasma is charge neutral, ne = ni
            4. Ionization Se is needed to balance diffusion loss.
        Da = Di + De*Mui/Mue
        Da = Di(1 + Te/Ti).
        Es = (Di - De)/(Mui + Mue)*dn/dx/n
        """
        self.Da = CCP_1d.Di*(1 + np.divide(CCP_1d.Te, CCP_1d.Ti))
        self.Ea = (CCP_1d.Di - CCP_1d.De)/(CCP_1d.Mui + CCP_1d.Mue)
        dnidx = CCP_1d.geom.cnt_diff(CCP_1d.ni)
        self.Ea *= np.divide(dnidx, CCP_1d.ni,
                             out=np.zeros_like(dnidx), where=CCP_1d.ni != 0)
        self.Ea[1] = self.Ea[2]
        self.Ea[-2] = self.Ea[-3]
        return self.Da, self.Ea
