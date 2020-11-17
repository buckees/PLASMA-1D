"""
Geometry module 1D, constructing the 1D geometry.

1D geometry is made of intervals.
Assign materials to the shapes, such as
'Metal', 'Quartz', etc.

1D geometry is defined separately from 2D geometry,
but they share the same strucuture.
"""

from Domain import Domain_1d


class Geom_1d(Domain_1d):
    """Define 1d Geometry."""

    def __init__(self, label, width):
        super().__init__(label, width)

    def add_intv(self):
        """Add intervals."""
        pass

    def add_bndy(self):
        """Add boundaries."""
        pass

    def add_mat(self):
        """Add materials."""
        pass

if __name__ == '__main__':
    """Test Geometry."""
    geom1d = Geom_1d('A', 10e-2)
    print(geom1d)
