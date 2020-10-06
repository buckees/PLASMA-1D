"""
Define Geometry.

Geometry --> 1d or 2d.
Geometry contains:
    intervals
    materials
    boundaries
"""

from Domain import Domain_1d


class Geom_1d(Domain_1d):
    """Define 1d Geometry."""

    def __init__(self, label, width):
        super().__init__(label, width)

    def add_bndy(self):
        """Add boundaries."""
        pass
