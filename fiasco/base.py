"""
Base class for access to CHIANTI ion data
"""
import os

import numpy as np
import h5py
import astropy.units as u
from astropy.table import QTable
import plasmapy.atomic

import fiasco
from .io.factory import all_subclasses
from .io.generic import GenericParser
from .util import download_dbase, build_hdf5_dbase


class DataIndexer(object):
    
    def __init__(self, hdf5_path, top_level_path):
        self.top_level_path = top_level_path
        self.hdf5_dbase_root = hdf5_path

    @property
    def fields(self):
        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            fields = [k for k in hf[self.top_level_path]]
        return fields

    def as_table(self):
        qt = QTable()
        for field in self.fields:
            qt[field] = self[field]
        return qt
    
    def __contains__(self, key):
        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            key_in_grp = key in hf[self.top_level_path]
        return key_in_grp
    
    def __getitem__(self, key):
        if type(key) is int:
            raise NotImplementedError('Iteration not supported.')
        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            grp = hf[self.top_level_path]
            if key not in grp:
                raise IndexError('{} not found in {} filetype'.format(key, self.top_level_path))
            ds = grp[key]
            if isinstance(ds, h5py.Group):
                data = DataIndexer(self.hdf5_dbase_root, '/'.join([self.top_level_path, key]))
            else:
                if ds.attrs['unit'] == 'SKIP':
                    data = np.array(ds, dtype=ds.dtype)
                else:
                    data = u.Quantity(np.array(ds), ds.attrs['unit'], dtype=ds.dtype)
                if '|S' in data.dtype.str:
                    data = data.astype(str)
        return data
    
    def __repr__(self):
        def ufilter(x):
            return 'unit' not in x.attrs or x.attrs['unit'] == 'SKIP' or x.attrs['unit'] == ''
        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            grp = hf[self.top_level_path]
            var_names = [(key, '', '{}'.format(grp[key].attrs['description'])) if ufilter(grp[key]) 
                         else (key, '({})'.format(grp[key].attrs['unit']), '{}'.format(grp[key].attrs['description'])) 
                         for key in grp]
            footer = grp.attrs['footer']
            
        name_strs = '\n'.join(['{} {} -- {}'.format(*v) for v in var_names])
        return '''{top_level_path}

Fields
------
{vars_and_units}

Footer
------
{footer}'''.format(top_level_path=self.top_level_path, vars_and_units=name_strs, footer=footer)


class IonBase(object):
    """
    Base class for accessing CHIANTI data attached to a particular ion
    """
    
    def __init__(self, ion_name, hdf5_path=None, **kwargs):
        self._ion_name = ion_name
        self.atomic_symbol = ion_name.split('_')[0].capitalize()
        self.atomic_number = plasmapy.atomic.atomic_number(self.atomic_symbol)
        self.element_name = plasmapy.atomic.element_name(self.atomic_symbol)
        self.ionization_stage = int(ion_name.split('_')[-1])
        self.charge_state = self.ionization_stage - 1
        self.ion_name = '{} {}'.format(self.atomic_symbol, self.ionization_stage)
        if hdf5_path is None:
            self.hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
        else:
            self.hdf5_dbase_root = hdf5_path
        download_dbase(fiasco.defaults['ascii_dbase_root'], ask_before=kwargs.get('ask_before', True))
        build_hdf5_dbase(fiasco.defaults['ascii_dbase_root'], self.hdf5_dbase_root,
                         ask_before=kwargs.get('ask_before', True))
       
    @property
    def abundance(self):
        return DataIndexer(self.hdf5_dbase_root, '/'.join([self.atomic_symbol.lower(), 'abundance']))
        

def add_property(cls, filetype):
    """
    Dynamically add filetype properties to base data access class
    """
    def property_template(self):
        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            if '/'.join([self.atomic_symbol.lower(), self._ion_name, filetype]) not in hf:
                return None
        return DataIndexer(self.hdf5_dbase_root, '/'.join([self.atomic_symbol.lower(), self._ion_name, filetype]))

    property_template.__doc__ = 'Data in {} type file'.format(filetype)
    property_template.__name__ = '_'.join(filetype.split('/'))
    setattr(cls, property_template.__name__, property(property_template))
    
# Collect the filetypes and add the methods
all_ext = [cls.filetype for cls in all_subclasses(GenericParser) 
           if hasattr(cls, 'filetype') and cls.filetype not in ['abund']]
for filetype in all_ext:
    add_property(IonBase, filetype)
    add_property(IonBase, '/'.join(['dielectronic', filetype]))
