"""
Element class for logically grouping ions
"""
import numpy as np
import h5py
import astropy.units as u
import plasmapy

import fiasco


class Element(object):
    """
    Object containing all ions for a particular element.

    The Element object provides a way to logically group together ions of the same
    element. This provides easy ways to compute element-level derived quantities. 

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
        self._ion_kwargs = kwargs.get('ion_kwargs', {})
        self._ion_kwargs['hdf5_path'] = self.hdf5_dbase_root

    @property
    def ions(self):
        """
        All ions available in CHIANTI for this element
        """
        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            ions = sorted([i.split('_') for i in hf[self.atomic_symbol.lower()].keys() 
                           if '{}_'.format(self.atomic_symbol.lower()) in i], key=lambda x: int(x[1]))
        return ['_'.join(i) for i in ions]

    def _rate_matrix(self):
        rate_matrix = np.zeros(self.temperature.shape + (self.atomic_number+1, self.atomic_number+1))
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

    def __add__(self, value):
        return fiasco.IonCollection(self, value)

    def __radd__(self, value):
        return fiasco.IonCollection(value, self)

    def __getitem__(self, value):
        if type(value) is int:
            value = self.ions[value]
        return fiasco.Ion(value, self.temperature, **self._ion_kwargs)

    def __contains__(self, value):
        return value in self.ions

    def __repr__(self):
        ion_list = '\n'.join([f"{i.split('_')[0].capitalize()} {i.split('_')[1]}"
                              for i in self.ions])
        return f"""Element
-------
{self.atomic_symbol} ({self.atomic_number}) -- {self.element_name}

Available Ions
--------------
{ion_list}"""
