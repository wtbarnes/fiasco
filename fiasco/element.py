"""
Element class for logically grouping ions
"""
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
        self.atomic_symbol = plasmapy.atomic.atomic_symbol(element_name)
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