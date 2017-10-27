"""
Element class for logically grouping ions
"""
import h5py
import periodictable
import astropy.units as u

import fiasco
from .ion import Ion


class Element(object):
    """
    Logical grouping of ion objects according to their element
    """

    @u.quantity_input
    def __init__(self, element_name, temperature: u.K, hdf5_path=None, **kwargs):
        self.temperature = temperature
        if type(element_name) is str:
            self.element_name = element_name.lower()
        elif type(element_name) is int:
            self.element_name = periodictable.elements[element_name].symbol.lower()
        else:
            raise ValueError('Unrecognized input type {}'.format(type(element_name)))
        if hdf5_path is None:
            self.hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
        else:
            self.hdf5_dbase_root = hdf5_path

    @property
    def ions(self):
        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            ions = sorted([i.split('_') for i in hf[self.element_name].keys() if '{}_'.format(self.element_name) in i],
                          key=lambda x: int(x[1]))
        return ['_'.join(i) for i in ions]

    def __getitem__(self, x):
        if type(x) is int:
            x = self.ions[x]
        return Ion(x, self.temperature, hdf5_path=self.hdf5_dbase_root)

    def __contains__(self, x):
        return x in self.ions

    def __repr__(self):
        ion_list = ['{} {}'.format(i.split('_')[0].capitalize(), i.split('_')[1]) for i in self.ions]
        return '''Element
-------
{}

Available Ions
--------------
{}'''.format(self.element_name.capitalize(), '\n'.join(ion_list))