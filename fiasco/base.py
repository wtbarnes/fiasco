"""
Base classes for access to CHIANTI ion data
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

__all__ = ['DataIndexer', 'IonBase', 'ContinuumBase']


class DataIndexer(object):
    """
    Data access layer for each distinct CHIANTI dataset

    Acts as an interface layer between `Ion` and the CHIANTI data stored in the
    HDF5 database. All data that the user interacts with passes through this layer.

    .. warning:: This object is not meant to be instantiated directly by the user. Rather, instances
                 are created by higher-level objects in order to provide access to the CHIANTI data.

    Parameters
    ----------
    hdf5_path : `str`
    top_level_path : `str`
    """
    
    def __init__(self, hdf5_path, top_level_path):
        self.top_level_path = top_level_path
        self.hdf5_dbase_root = hdf5_path

    @classmethod
    def create_indexer(cls, hdf5_path, top_level_path):
        """
        Create an instance as long as the dataset exists
        """
        with h5py.File(hdf5_path, 'r') as hf:
            path_is_valid = True if top_level_path in hf else False
        return cls(hdf5_path, top_level_path) if path_is_valid else None

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
            if key not in self:
                return None
            ds = hf[self.top_level_path][key]
            if isinstance(ds, h5py.Group):
                data = DataIndexer.create_indexer(self.hdf5_dbase_root,
                                                  '/'.join([self.top_level_path, key]))
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
            return ('' if 'unit' not in x.attrs or x.attrs['unit'] == 'SKIP'
                    or x.attrs['unit'] == '' else '({})'.format(x.attrs['unit']))

        def dfilter(x):
            return '' if 'description' not in x.attrs else '{}'.format(x.attrs['description'])

        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            grp = hf[self.top_level_path]
            var_names = [(key, ufilter(grp[key]), dfilter(grp[key])) for key in grp]
            footer = '' if 'footer' not in grp.attrs else grp.attrs['footer']
            
        name_strs = '\n'.join(['{} {} -- {}'.format(*v) for v in var_names])
        version = '' if self.version is None else f'-- v{self.version}'
        return f"""{self.top_level_path} {version}

Fields
------
{name_strs}

Footer
------
{footer}"""


class Base(object):
    """
    Base class for setting up ion metadata and building database as needed
    """

    def __init__(self, ion_name, hdf5_path=None, **kwargs):
        element, ion = ion_name.split()
        if '+' in ion:
            ion = f"{int(ion.strip('+')) + 1}"
        self.atomic_number = plasmapy.atomic.atomic_number(element.capitalize())
        self.element_name = plasmapy.atomic.element_name(element.capitalize())
        self.atomic_symbol = plasmapy.atomic.atomic_symbol(element.capitalize())
        self.ionization_stage = int(ion)
        self.charge_state = self.ionization_stage - 1
        self.ion_name = f'{self.atomic_symbol} {self.ionization_stage}'
        # This is the preferred CHIANTI format and is only preserved for internal data access
        self._ion_name = f'{self.atomic_symbol.lower()}_{self.ionization_stage}'
        if hdf5_path is None:
            self.hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
        else:
            self.hdf5_dbase_root = hdf5_path
        ask_before = kwargs.get('ask_before', True)
        download_dbase(fiasco.defaults['ascii_dbase_root'], ask_before=ask_before)
        build_hdf5_dbase(fiasco.defaults['ascii_dbase_root'], self.hdf5_dbase_root,
                         ask_before=ask_before)


class ContinuumBase(Base):
    """
    Base class for retrieving continuum datasets.
    """

    @property
    def _gffgu(self):
        data_path = '/'.join(['continuum', 'gffgu'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _gffint(self):
        data_path = '/'.join(['continuum', 'gffint'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _klgfb(self):
        data_path = '/'.join(['continuum', 'klgfb'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _verner(self):
        data_path = '/'.join([self.atomic_symbol.lower(), self._ion_name, 'continuum',
                              'verner_short'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _itoh(self):
        data_path = '/'.join([self.atomic_symbol.lower(), 'continuum', 'itoh'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _hseq(self):
        data_path = '/'.join([self.atomic_symbol.lower(), 'continuum', 'hseq_2photon'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _heseq(self):
        data_path = '/'.join([self.atomic_symbol.lower(), 'continuum', 'heseq_2photon'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)


class IonBase(Base):
    """
    Base class for accessing CHIANTI data attached to a particular ion

    Parameters
    ----------
    ion_name : `str`
        Name of ion, e.g. for Fe V, 'Fe 5', 'iron 5', 'Fe 4+'
    hdf5_path : `str`, optional
    
    Examples
    --------
    """
       
    @property
    def _abundance(self):
        data_path = '/'.join([self.atomic_symbol.lower(), 'abundance'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _ip(self):
        data_path = '/'.join([self.atomic_symbol.lower(), self._ion_name, 'ip'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _ioneq(self):
        data_path = '/'.join([self.atomic_symbol.lower(), self._ion_name, 'ioneq'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)


def add_property(cls, filetype):
    """
    Dynamically add filetype properties to base data access class
    """
    def property_template(self):
        data_path = '/'.join([self.atomic_symbol.lower(), self._ion_name, filetype])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    property_template.__doc__ = f'Data in {filetype} type file'
    property_template.__name__ = '_{}'.format('_'.join(filetype.split('/')))
    setattr(cls, property_template.__name__, property(property_template))

# Collect the filetypes and add the methods
all_ext = [cls.filetype for cls in all_subclasses(GenericIonParser) if hasattr(cls, 'filetype')]
for filetype in all_ext:
    add_property(IonBase, filetype)
    add_property(IonBase, '/'.join(['dielectronic', filetype]))