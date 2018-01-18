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
from .io.generic import GenericParser, GenericIonParser
from .util import download_dbase, build_hdf5_dbase


class DataIndexer(object):
    """
    Data access layer for each distinct CHIANTI dataset

    Acts as an interface layer between `Ion` and the CHIANTI stored in the
    HDF5 database. All data that the user interacts with passes through this layer. 
    """
    
    def __init__(self, hdf5_path, top_level_path):
        self.top_level_path = top_level_path
        self.hdf5_dbase_root = hdf5_path

    @property
    def version(self):
        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            if 'chianti_version' in hf[self.top_level_path].attrs:
                version = hf[self.top_level_path].attrs['chianti_version']
            else:
                version = None
        return version
    
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
                if ds.attrs['unit'] == 'SKIP' or ds.dtype == 'object':
                    data = np.array(ds, dtype=ds.dtype)
                else:
                    data = u.Quantity(np.array(ds), ds.attrs['unit'], dtype=ds.dtype)
                if '|S' in data.dtype.str:
                    data = data.astype(str)
        return data
    
    def __repr__(self):

        def ufilter(x):
            return ('' if 'unit' not in x.attrs or x.attrs['unit'] == 'SKIP' or x.attrs['unit'] == ''
                    else '({})'.format(x.attrs['unit']))

        def dfilter(x):
            return '' if 'description' not in x.attrs else '{}'.format(x.attrs['description'])

        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            grp = hf[self.top_level_path]
            var_names = [(key, ufilter(grp[key]), dfilter(grp[key])) for key in grp]
            footer = '' if 'footer' not in grp.attrs else grp.attrs['footer']
            
        name_strs = '\n'.join(['{} {} -- {}'.format(*v) for v in var_names])
        return '''{top_level_path} {version}

Fields
------
{vars_and_units}

Footer
------
{footer}'''.format(top_level_path=self.top_level_path, vars_and_units=name_strs, footer=footer,
                   version='' if self.version is None else '-- v{}'.format(self.version))


class IonBase(object):
    """
    Base class for accessing CHIANTI data attached to a particular ion

    Examples
    --------

    Notes
    -----
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
    def _abundance(self):
        return DataIndexer(self.hdf5_dbase_root, '/'.join([self.atomic_symbol.lower(), 'abundance']))

    @property
    def _ip(self):
        return DataIndexer(self.hdf5_dbase_root, '/'.join([self.atomic_symbol.lower(), self._ion_name, 'ip']))

    @property
    def _ioneq(self):
        return DataIndexer(self.hdf5_dbase_root, '/'.join([self.atomic_symbol.lower(), self._ion_name, 'ioneq']))
      

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
    property_template.__name__ = '_{}'.format('_'.join(filetype.split('/')))
    setattr(cls, property_template.__name__, property(property_template))

# Collect the filetypes and add the methods
all_ext = [cls.filetype for cls in all_subclasses(GenericIonParser) if hasattr(cls, 'filetype')]
for filetype in all_ext:
    add_property(IonBase, filetype)
    add_property(IonBase, '/'.join(['dielectronic', filetype]))
