"""
Classes and functions for element-level operations
"""
import numpy as np
import astropy.units as u
import plasmapy

import fiasco

__all__ = ['Element']


class Element(fiasco.IonCollection):
    """
    Object containing all ions for a particular element.

    The `Element` object provides a way to logically group together ions of the same
    element. This provides an easy way to compute element-level derived quantities such
    as the ionization fraction as a function of temperature.

    Parameters
    ----------
    element_name : `str`, `int`
        Symbol, atomic number, or full name of the element
    temperature : `~astropy.units.Quantity`

    See Also
    --------
    fiasco.Ion : All the same keyword arguments can also be passed here.
    """

    @u.quantity_input
    def __init__(self, element_name, temperature: u.K, **kwargs):
        if type(element_name) is str:
            element_name = element_name.capitalize()
        self.atomic_symbol = plasmapy.atomic.atomic_symbol(element_name)
        self.atomic_number = plasmapy.atomic.atomic_number(element_name)
        self.element_name = plasmapy.atomic.element_name(element_name)

        ion_list = []
        # NOTE: check only the first ion to see if it is in the database as
        # checking all of them is very slow
        for i in range(self.atomic_number + 1):
            ion = fiasco.Ion(f'{self.atomic_symbol} {i+1}', temperature, **kwargs)
            ion_list.append(ion)

        super().__init__(*ion_list)

    @property
    def abundance(self):
        return self[0].abundance

    def _rate_matrix(self):
        rate_matrix = np.zeros(self.temperature.shape+(self.atomic_number+1, self.atomic_number+1))
        rate_unit = self[0].ionization_rate().unit
        rate_matrix = rate_matrix * rate_unit
        for i in range(1, self.atomic_number):
            rate_matrix[:, i, i] = -(self[i].ionization_rate() + self[i].recombination_rate())
            rate_matrix[:, i, i-1] = self[i-1].ionization_rate()
            rate_matrix[:, i, i+1] = self[i+1].recombination_rate()
        rate_matrix[:, 0, 0] = -(self[0].ionization_rate() + self[0].recombination_rate())
        rate_matrix[:, 0, 1] = self[1].recombination_rate()
        rate_matrix[:, -1, -1] = -(self[-1].ionization_rate() + self[-1].recombination_rate())
        rate_matrix[:, -1, -2] = self[-2].ionization_rate()

        return rate_matrix

    def equilibrium_ionization(self, **kwargs):
        """
        Calculate the ionization equilibrium for all ions of the element.

        Calculate the population fractions for every ion of this element as a function of
        temperature, assuming ionization equilibrium.
        """
        rate_matrix = kwargs.get('rate_matrix', None)
        if rate_matrix is None:
            rate_matrix = self._rate_matrix()
        # Solve system of equations using singular value decomposition
        _, _, V = np.linalg.svd(rate_matrix.value)
        # Select columns of V with smallest eigenvalues (returned in descending order)
        # NOTE: must take the absolute value as the SVD solution is only accurate up
        # to the sign. We require that the solutions must be positive.
        ioneq = np.fabs(V[:, -1, :])
        ioneq /= ioneq.sum(axis=1)[:, np.newaxis]

        return u.Quantity(ioneq)

    def __getitem__(self, value):
        if type(value) is str:
            el, ion = value.split()
            if '+' in ion:
                value = int(ion.strip('+'))
            else:
                value = int(ion) - 1
        return super().__getitem__(value)

    def __repr__(self):
        ion_list = '\n'.join([i.ion_name for i in self._ion_list])
        return f"""Element
-------
{self.atomic_symbol} ({self.atomic_number}) -- {self.element_name}

Available Ions
--------------
{ion_list}"""
