"""
Define Mesh.

Mesh --> 1d or 2d.
Mesh contains:
    x position
"""

from Geometry import Geom_1d

import numpy as np


class Mesh_1d(Geom_1d):
    """Define 1d Mesh."""

    def __init__(self, label, width, nx=11, res=0.1):
        """Add option to choose nx or res for mesh."""
        super().__init__(label, width)
        self.nx = nx
        self.delx = self.width/(self.nx-1)
        self.x = np.linspace(0.0, self.width, self.nx)

    def __str__(self):
        """Print 1d mesh information."""
        res = 'Mesh_1d:'
        res += '\n ' + super().__str__()
        res += f'\nnx = {self.nx}'
        res += f'\ndelx = {self.delx} m'
        return res

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
    geom1d = Geom_1d('A', 10e-2)
    mesh1d = Mesh_1d(geom1d.label, geom1d.width,
                     nx=101)
    # print(geom1d)
    print(mesh1d)
