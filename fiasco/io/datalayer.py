"""
Access layer for interfacing with CHIANTI stored in HDF5.
"""
import astropy.units as u
import h5py
import numpy as np
import pathlib

from astropy.table import QTable

from fiasco.util.exceptions import MissingDatabaseError

__all__ = ['DataIndexer', 'DataIndexerHDF5']


class DataIndexer:
    """
    Data access layer for each distinct CHIANTI dataset

    Acts as an interface layer between `~fiasco.Ion` and the CHIANTI data. All data that the user interacts
    with passes through this layer.

    .. warning:: This object is not meant to be instantiated directly by the user. Rather, instances
                 are created by higher-level objects in order to provide access to the CHIANTI data.
    """

    def __new__(cls, *args):
        return DataIndexerHDF5(*args)

    @classmethod
    def create_indexer(cls, *args):
        return DataIndexerHDF5.create_indexer(*args)


class DataIndexerHDF5:
    """
    Interface layer for CHIANTI data stored in HDF5 format.

    Parameters
    ----------
    hdf5_path : `str`
    top_level_path : `str`
    """

    def __init__(self, hdf5_path, top_level_path):
        self.top_level_path = top_level_path
        self._hdf5_dbase_root = hdf5_path

    @property
    def hdf5_dbase_root(self):
        dbase_root = pathlib.Path(self._hdf5_dbase_root)
        if not dbase_root.is_file():
            raise MissingDatabaseError(f'No HDF5 database found at {dbase_root}')
        return dbase_root

    @classmethod
    def create_indexer(cls, hdf5_path, top_level_path):
        """
        Create an instance as long as the dataset exists. This class method
        exists so that None can be returned if the dataset specified by
        ``top_level_path`` does not exist.
        """
        hdf5_path = pathlib.Path(hdf5_path)
        if not hdf5_path.is_file():
            raise MissingDatabaseError(f'No HDF5 database found at {hdf5_path}')
        with h5py.File(hdf5_path, 'r') as hf:
            path_is_valid = top_level_path in hf
        if path_is_valid:
            return cls(hdf5_path, top_level_path)
        else:
            raise KeyError(f'{top_level_path} not found in {hdf5_path}')

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
        if isinstance(key, int):
            raise NotImplementedError('Iteration not supported.')
        if key not in self:
            raise KeyError(f'{key} not found in {self.top_level_path}')
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

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

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
