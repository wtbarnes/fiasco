"""
Classes and functions for element-level operations
"""
import numpy as np
import h5py
import astropy.units as u
import plasmapy

import fiasco

__all__ = ['Element']


class Element(fiasco.IonCollection):
    """
    Object containing all ions for a particular element.

    The Element object provides a way to logically group together ions of the same
    element. This provides easy ways to compute element-level derived quantities such 
    as the population fractions.

    Parameters
    ----------
    element_name : `str`, `int`
        Symbol, atomic number, or full name of the element
    temperature : `~astropy.units.Quantity`
    hdf5_path : `str`, optional
        Path to HDF5 database; defaults to that listed in `~fiasco.defaults`

    Other Parameters
    ----------------
    ion_kwargs : `dict`
        Possible keyword arguments for individual ions

    Examples
    --------
    """

    @u.quantity_input
    def __init__(self, element_name, temperature: u.K, hdf5_path=None, **kwargs):
        self.temperature = temperature
        if type(element_name) is str:
            element_name = element_name.capitalize()
        self.atomic_symbol = plasmapy.atomic.atomic_symbol(element_name)
        self.atomic_number = plasmapy.atomic.atomic_number(element_name)
        self.element_name = plasmapy.atomic.element_name(element_name)
        if hdf5_path is None:
            self.hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
        else:
            self.hdf5_dbase_root = hdf5_path
        ion_kwargs = kwargs.get('ion_kwargs', {})
        ion_kwargs['hdf5_path'] = self.hdf5_dbase_root

        chianti_ions = fiasco.DataIndexer(self.hdf5_dbase_root, self.atomic_symbol.lower()).fields
        ion_list = []
        for i in range(self.atomic_number + 1):
            ion = f'{self.atomic_symbol.lower()}_{i+1}'
            if ion in chianti_ions:
                ion_list.append(fiasco.Ion(f'{self.atomic_symbol} {i+1}',
                                temperature, **ion_kwargs))
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

    def equilibrium_ionization(self, rate_matrix=None):
        """
        Calculate the ionization equilibrium for all ions of the element.

        Calculate the population fractions for every ion of this element as a function of
        temperature, assuming ionization equilibrium.

        Parameters
        ----------
        rate_matrix : `~astropy.units.Quantity`, optional
            Precomputed matrix of total ionization and recombination rates
        """
        if rate_matrix is None:
            rate_matrix = self._rate_matrix()
        # Solve system of equations using singular value decomposition
        _, _, V = np.linalg.svd(rate_matrix.value)
        # Select columns of V with smallest eigenvalues (returned in descending order)
        ioneq = np.fabs(V[:, -1, :])
        ioneq /= np.sum(ioneq, axis=1)[:, np.newaxis]

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
