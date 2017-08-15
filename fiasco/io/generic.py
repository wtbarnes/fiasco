"""
Base class for file parser
"""
import os
import h5py
import numpy as np
from astropy.table import QTable
import astropy.units as u
import fiasco


class GenericParser(object):
    """
    Base class for CHIANTI file parsers
    """

    dtypes = []
    units = []
    headings = []
    
    def __init__(self,ion_filename,**kwargs):
        self.ion_filename = ion_filename
        self.filetype = self.ion_filename.split('.')[-1]
        self.ion_name = self.ion_filename.split('.')[0]
        self.element = self.ion_name.split('_')[0]
        self.full_path = os.path.join(fiasco.defaults['chianti_dbase_root'],
                                      self.element,self.ion_name,self.ion_filename)
        
    def parse(self):
        """
        Generate Astropy QTable from a CHIANTI ion file 
        """
        if not os.path.isfile(self.full_path):
            return None
        with open(self.full_path, 'r') as f:
            lines = f.readlines()
        table = []
        for i,line in enumerate(lines):
            if line.strip() == '-1':
                comment = ''.join(lines[i+1:len(lines)])
                break
            else:
                self.preprocessor(table,line,i)

        df = QTable(data=list(map(list,zip(*table))), names=self.headings)
        for name,unit,dtype in zip(self.headings, self.units, self.dtypes):
            df[name].unit = unit
            df[name] = df[name].astype(dtype)
        
        df.meta['footer'] = comment
        with open(os.path.join(fiasco.defaults['chianti_dbase_root'],'VERSION'),'r') as f:
            lines = f.readlines()
            version = lines[0].strip()
        df.meta['chianti_version'] = version

        df = self.postprocessor(df)
        
        return df
    
    def preprocessor(self,table,line,index):
        """
        Default preprocessor method run on each line ingested 
        """
        if hasattr(self,'fformat'):
            line = self.fformat.read(line)
        else:
            line = line.strip().split()
        line = [item.strip() if type(item) is str else item for item in line]
        table.append(line)
        
    def postprocessor(self,df):
        """
        Default postprocessor method run on the whole dataframe
        """
        df.meta['element'] = self.element
        df.meta['ion'] = self.ion_name
        return df
    
    def to_hdf5(self,hf,df):
        """
        Add datasets to a group for an HDF5 file handler
        """
        grp_name = '/'.join([self.element,self.ion_name,self.filetype])
        if grp_name not in hf:
            grp = hf.create_group(grp_name)
            grp.attrs['chianti_version'] = df.meta['chianti_version']
            grp.attrs['footer'] = df.meta['footer']
        else:
            grp = hf[grp_name]
        hf['/'.join([self.element,self.ion_name])].attrs['element'] = self.element
        hf['/'.join([self.element,self.ion_name])].attrs['ion'] = self.ion_name
        for name in df.colnames:
            col = df[name]
            if type(col) == u.Quantity:
                data = col.value
            else:
                data = col.data
            if '<U' in data.dtype.str:
                numchar = data.dtype.str[2:]
                data = data.astype('|S{}'.format(numchar))
            ds = self._write_to_hdf5(grp,name,data)
            if col.unit is None:
                ds.attrs['unit'] = ''
            else:
                ds.attrs['unit'] = col.unit.to_string()

    def _write_to_hdf5(self,grp,name,data):
        return grp.create_dataset(name,data=data,dtype=data.dtype)
        
        
        