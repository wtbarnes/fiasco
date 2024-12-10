"""
Classes and functions for element-level operations
"""
import astropy.units as u
import numpy as np
import plasmapy

from functools import cached_property

import fiasco

from fiasco.util import parse_ion_name

__all__ = ['Element']


class Element(fiasco.IonCollection):
    """
    Collection of all ions for a particular element.

    The `Element` object provides a way to logically group together ions of the same
    element. This makes it easy to compute element-level derived quantities such
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
        if isinstance(element_name, str):
            element_name = element_name.capitalize()
        Z = plasmapy.particles.atomic_number(element_name)
        ion_list = []
        for i in range(Z + 1):
            ion = fiasco.Ion((Z, i+1), temperature, **kwargs)
            ion_list.append(ion)

        super().__init__(*ion_list)

    @property
    def atomic_symbol(self):
        """The standard atomic symbol for the element, e.g. "H" for hydrogen."""
        return self[0].atomic_symbol

    @property
    def atomic_number(self):
        """The atomic number of the element, :math:`Z`."""
        return self[0].atomic_number

    @property
    def element_name(self):
        """The full name of the element, e.g. "hydrogen"."""
        return self[0].element_name

    @property
    def abundance(self):
        "Elemental abundance relative to H."
        return self[0].abundance

    @abundance.setter
    def abundance(self, abundance):
        for _ion in self:
            _ion.abundance = abundance

    @cached_property
    def _rate_matrix(self):
        rate_matrix = np.zeros(self.temperature.shape+(self.atomic_number+1, self.atomic_number+1))
        rate_unit = self[0].ionization_rate.unit
        rate_matrix = rate_matrix * rate_unit
        for i in range(1, self.atomic_number):
            rate_matrix[:, i, i] = -(self[i].ionization_rate + self[i].recombination_rate)
            rate_matrix[:, i, i-1] = self[i-1].ionization_rate
            rate_matrix[:, i, i+1] = self[i+1].recombination_rate
        rate_matrix[:, 0, 0] = -(self[0].ionization_rate + self[0].recombination_rate)
        rate_matrix[:, 0, 1] = self[1].recombination_rate
        rate_matrix[:, -1, -1] = -(self[-1].ionization_rate + self[-1].recombination_rate)
        rate_matrix[:, -1, -2] = self[-2].ionization_rate

        return rate_matrix

    @cached_property
    def equilibrium_ionization(self):
        """
        The ionization fraction, in equilibrium, for all ions of the element.

        The population fractions for every ion of this element, assuming
        ionization equilibrium, calculated as a function of temperature.
        This returns a matrix with dimensions ``(n,Z+1)``, where ``n``
        corresponds to the temperature dimension and ``Z+1`` corresponds to
        the number of ionization stages of the element.

        .. note:: The result here is not simply the tabulated values included in
                  the CHIANTI database. Rather, these values are calculated as a
                  function of temperature on the fly using the ionization and
                  recombination rates associated with each ion in the element.

        Examples
        --------
        >>> temperature = 10**np.arange(3.9, 6.5, 0.01) * u.K
        >>> carbon = Element('C', temperature)
        >>> carbon_ionization = carbon.equilibrium_ionization
        >>> carbon_ionization[:, 4].max()  # max population fraction of C V as a function of temperature
        <Quantity 0.99776769>

        See Also
        --------
        fiasco.Ion.ionization_rate
        fiasco.Ion.recombination_rate
        """
        # Solve system of equations using singular value decomposition
        _, _, V = np.linalg.svd(self._rate_matrix.value)
        # Select columns of V with smallest eigenvalues (returned in descending order)
        # NOTE: must take the absolute value as the SVD solution is only accurate up
        # to the sign. We require that the solutions must be positive.
        ionization_fraction = np.fabs(V[:, -1, :])
        ionization_fraction /= ionization_fraction.sum(axis=1)[:, np.newaxis]

        return u.Quantity(ionization_fraction)

    def __getitem__(self, value):
        if isinstance(value, (str, tuple)):  # NOQA: UP038
            _, value = parse_ion_name(value)
            value -= 1
        return super().__getitem__(value)

    def __repr__(self):
        ion_name_list = '\n'.join([i.ion_name for i in self._ion_list])
        return f"""Element
-------
{self.atomic_symbol} ({self.atomic_number}) -- {self.element_name}
Temperature range: [{self.temperature[0].to(u.MK):.3f}, {self.temperature[-1].to(u.MK):.3f}]

Available Ions
--------------
{ion_name_list}"""
