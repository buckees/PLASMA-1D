"""
Define Domain.

Domain --> 1d or 2d.
Domain has boundaries.
"""


class Domain(object):
    """Define Domain label."""

    def __init__(self, label):
        self.label = label

    def __str__(self):
        """Print Domain label."""
        return f'label = {self.label}'


class Domain_1d(Domain):
    """Define 1D Domain."""

    def __init__(self, label, width):
        super(Domain_1d, self).__init__(label)
        self.width = width

    def __str__(self):
        """Print 1d domain information."""
        res = 'Domain_1d:'
        res += '\n' + super().__str__()
        res += f'\nwidth = {self.width} m'
        return res


if __name__ == '__main__':
    """Test Domain."""
    domain1d = Domain_1d('A', 10e-2)
    print(domain1d)
