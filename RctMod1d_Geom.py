"""
Geometry module 1D, constructing the 1D geometry.

1D geometry is made of intervals.
Assign materials to the shapes, such as
'Metal', 'Quartz', etc.

1D geometry is defined separately from 2D geometry,
but they share the same strucuture.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patch
import numpy as np

class Interval():
    """Init the Interval."""
    
    def __init__(self, label):
        """
        Init the Shape.
        
        label: str, var, label of Interval.
        """
        self.label = label

    def __str__(self):
        """Print Shape info."""
        return f'label = {self.label}'

class Domain(Interval):
    """Define the Domain."""
    
    def __init__(self, bndy=(0.0, 1.0)):
        """
        Init the Domain.
        
        bndy: unit in m, (2, ) tuple, defines the domain
        label: Domain label is fixed to 'Plasma'
        """
        self.bndy = np.asarray(bndy)
        super().__init__(label='Plasma')

    def __str__(self):
        """Print Domain info."""
        res = 'Domain:'
        res += f'\nlabel = {self.label}'
        res += f'\nleft and right boundaries = {self.bndy} m'
        return res

class Segment(Interval):
    """Segment is a basic Interval."""
    
    def __init__(self, label, lr):
        """
        Init the Segment.
        
        lr: unit in m, (2, ) tuple, left and right boundaries of the segment
        label: str, var, label of Interval.
        type: str, var, type of Interval
        """
        super().__init__(label)
        self.lr = np.asarray(lr)
        self.type = 'Segment'

    def __str__(self):
        """Print Rectangle info."""
        res = 'Rectangle:'
        res += f'\nleft and right bndy = {self.lr} m'
        return res

    def __contains__(self, posn):
        """
        Determind if a position is inside the Interval.
        
        posn: unit in m, (2, ) array, position as input
        boundaries are not consindered as "Inside"
        """
        return all(self.bl <= posn) and all(posn <= self.ur)

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
