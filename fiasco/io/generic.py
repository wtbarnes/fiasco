"""
Base class for file parser.
"""
import astropy.units as u
import h5py
import numpy as np
import pathlib
import warnings

from astropy.table import QTable
from astropy.utils.exceptions import AstropyUserWarning

import fiasco

from fiasco.util.exceptions import MissingASCIIFileError
from fiasco.util.util import read_chianti_version

__all__ = ['GenericParser', 'GenericIonParser']

class GenericParser:
    """
    Base class for CHIANTI file parsers
    """
    dtypes = []
    units = []
    headings = []
    descriptions = []

    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.ascii_dbase_root = pathlib.Path(kwargs.get('ascii_dbase_root', fiasco.defaults['ascii_dbase_root']))

    @property
    def full_path(self):
        if self.standalone:
            return pathlib.Path(self.filename)
        elif hasattr(self, '_full_path'):
            return self._full_path
        else:
            return self.ascii_dbase_root / self.filename

    @full_path.setter
    def full_path(self, value):
        self._full_path = value

    @property
    def standalone(self):
        return pathlib.Path(self.filename).exists()

    @property
    def chianti_version(self):
        if not self.standalone:
            return read_chianti_version(self.ascii_dbase_root)

    def parse(self):
        """
        Generate Astropy QTable from a CHIANTI ion file
        """
        # NOTE: put this here and not in __init__ as __init__ may be overwritten in a subclass
        if not self.full_path.is_file():
            raise MissingASCIIFileError(f'Could not find file {self.full_path}')
        with self.full_path.open(encoding='UTF-8') as f:
            lines = f.readlines()
        table = []
        for i, line in enumerate(lines):
            # Footer denotes end of file and is fenced by -1's
            if line.strip() == '-1' or r'%file' in line:  # sometimes the first -1 is missing
                break
            else:
                self.preprocessor(table, line, i)

        df = QTable(data=list(map(list, zip(*table))), names=self.headings)
        # This catches a warning thrown when we convert a staggered array into
        # a unitful column. This happens in several of the scups files for the
        # bt_t, bt_type, bt_upsilon columns.
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=AstropyUserWarning)
            for name, unit, dtype in zip(self.headings, self.units, self.dtypes):
                df[name].unit = unit
                df[name] = df[name].astype(dtype)

        df.meta['footer'] = self.extract_footer(lines)
        df.meta['chianti_version'] = f'{self.chianti_version}' if self.chianti_version else ''

        df = self.postprocessor(df)

        return df

    def extract_footer(self, lines):
        """
        Extract metadata from raw text and format appropriately.
        """
        for i, line in enumerate(lines):
            if line.strip() == '-1' or r'%file' in line:  # sometimes the first -1 is missing
                i_start = i if r'%file' in line else i+1
                comment = lines[i_start:len(lines)]
                break

        footer = '\n'.join([l.strip('%').strip() for l in comment if l.strip('%').strip().strip('-1')])
        return footer

    def preprocessor(self, table, line, *args):
        """
        Default preprocessor method run on each line ingested.
        """
        if hasattr(self, 'fformat'):
            line = self.fformat.read(line)
        else:
            line = line.strip().split()
        line = [item.strip() if isinstance(item, str) else item for item in line]
        table.append(line)

    def postprocessor(self, df):
        """
        Default postprocessor method run on the whole dataframe
        """
        df.meta['filename'] = self.filename
        df.meta['descriptions'] = {h: d for h, d in zip(self.headings, self.descriptions)}
        return df

    def to_hdf5(self, *args, **kwargs):
        raise NotImplementedError('No method for converting QTable to HDF5')


class GenericIonParser(GenericParser):
    """
    Base class for CHIANTI files attached to a particular ion
    """
    def __init__(self, ion_filename, **kwargs):
        super().__init__(ion_filename, **kwargs)
        self.dielectronic = False
        self.ion_name = pathlib.Path(self.filename).stem
        if self.ion_name and self.ion_name[-1] == 'd':
            self.dielectronic = True
            self.ion_name = self.ion_name[:-1]
        self.element = self.ion_name.split('_')[0]

    @property
    def full_path(self):
        if self.standalone:
            return pathlib.Path(self.filename)
        else:
            return self.ascii_dbase_root / self.element / pathlib.Path(self.filename).stem / self.filename

    def postprocessor(self, df):
        df = super().postprocessor(df)
        df.meta['element'] = self.element
        df.meta['ion'] = self.ion_name
        df.meta['dielectronic'] = self.dielectronic
        return df

    def to_hdf5(self, hf, df, **kwargs):
        """
        Add datasets to a group for an HDF5 file handler
        """
        if self.dielectronic:
            grp_name = '/'.join([self.element, self.ion_name, 'dielectronic', self.filetype])
        else:
            grp_name = '/'.join([self.element, self.ion_name, self.filetype])
        if grp_name not in hf:
            grp = hf.create_group(grp_name)
            grp.attrs['chianti_version'] = df.meta['chianti_version']
            grp.attrs['footer'] = df.meta['footer']
        else:
            grp = hf[grp_name]
        hf['/'.join([self.element, self.ion_name])].attrs['element'] = self.element
        hf['/'.join([self.element, self.ion_name])].attrs['ion'] = self.ion_name
        for name in df.colnames:
            col = df[name]
            data = col.value if type(col) is u.Quantity else col.data
            if '<U' in data.dtype.str:
                numchar = data.dtype.str[2:]
                data = data.astype(f'|S{numchar}')
            if name in grp:
                ds = grp[name]
            else:
                if data.dtype == np.dtype('O'):
                    ragged_dtype = h5py.special_dtype(vlen=np.dtype('float64'))
                    ds = grp.create_dataset(name, data=data, dtype=ragged_dtype)
                else:
                    ds = grp.create_dataset(name, data=data, dtype=data.dtype)
            if col.unit is None:
                ds.attrs['unit'] = 'SKIP'
            else:
                ds.attrs['unit'] = col.unit.to_string()
            ds.attrs['description'] = df.meta['descriptions'][name]
