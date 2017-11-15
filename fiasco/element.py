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
    Logical grouping of ion objects according to their element

    Notes
    -----

    Examples
    --------
    """

    @u.quantity_input
    def __init__(self, element_name, temperature: u.K, hdf5_path=None, **kwargs):
        self.temperature = temperature
        if type(element_name) is str:
            element_name = element_name.capitalize()
        self.atomic_symbol = plasmapy.atomic.element_symbol(element_name)
        self.atomic_number = plasmapy.atomic.atomic_number(element_name)
        self.element_name = plasmapy.atomic.element_name(element_name)
        if hdf5_path is None:
            self.hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
        else:
            self.hdf5_dbase_root = hdf5_path

    @property
    def ions(self):
        """
        All ions available in CHIANTI for this element
        """
        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            ions = sorted([i.split('_') for i in hf[self.atomic_symbol.lower()].keys() 
                           if '{}_'.format(self.atomic_symbol.lower()) in i], key=lambda x: int(x[1]))
        return ['_'.join(i) for i in ions]

    def ionization_equilibrium(self):
        """
        Calculate the ionization equilibrium for all ions of the element.

        Brief explanation and equations about how these equations are solved.
        """
        # Make matrix of ionization and recombination rates
        a_matrix = np.zeros(self.temperature.shape + (self.atomic_number+1, self.atomic_number+1))
        for i in range(1, self.atomic_number):
            a_matrix[:, i, i] = -(self[i].ionization_rate() + self[i].recombination_rate()).value
            a_matrix[:, i, i-1] = self[i-1].ionization_rate().value
            a_matrix[:, i, i+1] = self[i+1].recombination_rate().value
        a_matrix[:, 0, 0] = -(self[0].ionization_rate() + self[0].recombination_rate()).value
        a_matrix[:, 0, 1] = self[1].recombination_rate().value
        a_matrix[:, -1, -1] = -(self[-1].ionization_rate() + self[-1].recombination_rate()).value
        a_matrix[:, -1, -2] = self[-2].ionization_rate().value
        
        # Solve system of equations using SVD and normalize
        _, _, V = np.linalg.svd(a_matrix)
        # Select columns of V with smallest eigenvalues (returned in descending order)
        ioneq = np.fabs(V[:, -1, :])
        ioneq /= np.sum(ioneq, axis=1)[:, np.newaxis]

        return ioneq

    def __add__(self, value):
        return fiasco.IonCollection(self, value)

    def __radd__(self, value):
        return fiasco.IonCollection(value, self)

    def __getitem__(self, value):
        if type(value) is int:
            value = self.ions[value]
        return fiasco.Ion(value, self.temperature, hdf5_path=self.hdf5_dbase_root)

    def __contains__(self, value):
        return value in self.ions

    def __repr__(self):
        ion_list = ['{} {}'.format(i.split('_')[0].capitalize(), i.split('_')[1]) for i in self.ions]
        return '''Element
-------
{} ({}) -- {}

Available Ions
--------------
{}'''.format(self.atomic_symbol, self.atomic_number, self.element_name, '\n'.join(ion_list))