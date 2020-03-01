"""
Access layer for interfacing with CHIANTI stored in HDF5
"""
import os

import numpy as np
import h5py
try:
    import h5pyd
except ImportError:
    pass
import astropy.units as u
from astropy.table import QTable

import fiasco
from fiasco.util.exceptions import MissingDatabaseError

__all__ = ['DataIndexer']


class DataIndexer(object):
    """
    Data access layer for CHIANTI data
    """

    def __new__(cls, *args):
        if fiasco.defaults['use_remote_data']:
            return DataIndexerRemote(
                fiasco.defaults['remote_domain'], fiasco.defaults['remote_endpoint'], args[1])
        else:
            return DataIndexerLocal(*args)

    @classmethod
    def create_indexer(cls, *args):
        if fiasco.defaults['use_remote_data']:
            return DataIndexerRemote.create_indexer(
                fiasco.defaults['remote_domain'], fiasco.defaults['remote_endpoint'], args[1])
        else:
            return DataIndexerLocal.create_indexer(*args)


class DataIndexerRemote(object):

    def __init__(self, domain, endpoint, top_level_path):
        self.domain = domain
        self.endpoint = endpoint
        self.top_level_path = top_level_path

    @classmethod
    def create_indexer(cls, domain, endpoint, top_level_path):
        """
        Create an instance as long as the dataset exists
        """
        with h5pyd.File(domain, 'r', endpoint=endpoint) as hf:
            path_is_valid = True if top_level_path in hf else False
        return cls(domain, endpoint, top_level_path) if path_is_valid else None

    @property
    def version(self):
        with h5pyd.File(self.domain, 'r', endpoint=self.endpoint) as hf:
            if 'chianti_version' in hf[self.top_level_path].attrs:
                version = hf[self.top_level_path].attrs['chianti_version']
            else:
                version = None
        return version

    @property
    def fields(self):
        with h5pyd.File(self.domain, 'r', endpoint=self.endpoint) as hf:
            fields = [k for k in hf[self.top_level_path]]
        return fields

    def as_table(self):
        qt = QTable()
        for field in self.fields:
            qt[field] = self[field]
        return qt

    def __contains__(self, key):
        with h5pyd.File(self.domain, 'r', endpoint=self.endpoint) as hf:
            key_in_grp = key in hf[self.top_level_path]
        return key_in_grp

    def __getitem__(self, key):
        """
        NOTE: There seems to be a weird in bug in h5pyd where if a dataset
        is returned directly to a numpy array, the slicing/indexing fails. Thus,
        all of the gymnastics around returning datasets and casting to types appropriately.
        """
        if type(key) is int:
            raise NotImplementedError('Iteration not supported.')
        with h5pyd.File(self.domain, 'r', endpoint=self.endpoint) as hf:
            if key not in self:
                return None
            ds = hf[self.top_level_path][key]
            if isinstance(ds, h5pyd.Group):
                data = DataIndexerRemote.create_indexer(self.domain, self.endpoint,
                                                        '/'.join([self.top_level_path, key]))
            else:
                # Scalars cannot be sliced
                if not ds.shape:
                    data = np.array(ds.value)
                else:
                    data = ds[:]
                # Some things are just arrays
                if ds.attrs['unit'] == 'SKIP' or ds.dtype == 'object':
                    data = data.astype(ds.dtype)
                else:
                    data = u.Quantity(data, ds.attrs['unit'], dtype=ds.dtype)
                if '|S' in data.dtype.str:
                    data = data.astype(str)
        return data

    def __repr__(self):
        with h5pyd.File(self.domain, 'r', endpoint=self.endpoint) as hf:
            grp = hf[self.top_level_path]
            var_names = [key for key in grp]
            footer = '' if 'footer' not in grp.attrs else grp.attrs['footer']

        name_strs = '\n'.join(var_names)
        version = '' if self.version is None else f'-- v{self.version}'
        return f"""{self.top_level_path} {version}

Fields
------
{name_strs}

Footer
------
{footer}"""


class DataIndexerLocal(object):
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
        Create an instance as long as the dataset exists. This class method
        exists so that None can be returned if the dataset specified by
        `top_level_path` does not exist.
        """
        if not os.path.isfile(hdf5_path):
            raise MissingDatabaseError(f'No HDF5 database found at {hdf5_path}')
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
    def footer(self):
        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            if 'footer' in hf[self.top_level_path].attrs:
                footer = hf[self.top_level_path].attrs['footer']
            else:
                footer = None
        return footer

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
        if key not in self:
            return None
        with h5py.File(self.hdf5_dbase_root, 'r') as hf:
            ds = hf[self.top_level_path][key]
            if isinstance(ds, h5py.Group):
                data = DataIndexer.create_indexer(
                    self.hdf5_dbase_root, '/'.join([self.top_level_path, key]))
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
