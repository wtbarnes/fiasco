"""
Base class for access to CHIANTI ion data
"""
import os

import numpy as np
import h5py
import astropy.units as u

import fiasco
from .io.factory import all_subclasses
from .io.generic import GenericParser


class DataIndexer(object):
    
    def __init__(self,hdf5_path,top_level_path):
        self.top_level_path = top_level_path
        self.hdf5_path = hdf5_path
    
    def __getitem__(self,key):
        if type(key) is int:
            raise NotImplementedError('Iteration not supported.')
        with h5py.File(self.hdf5_path,'r') as hf:
            grp = hf[self.top_level_path]
            if key not in grp:
                raise IndexError('{} not found in {} filetype'.format(key,self.top_level_path))
            ds = grp[key]
            if isinstance(ds,h5py.Group):
                data = DataIndexer(self.hdf5_path,'/'.join([self.top_level_path,key]))
            else:
                if ds.attrs['unit'] == 'SKIP':
                    data = np.array(ds,dtype=ds.dtype)
                else:
                    data = u.Quantity(np.array(ds),ds.attrs['unit'],dtype=ds.dtype)
                if '|S' in data.dtype.str:
                    data = data.astype(str)
        return data
    
    def __repr__(self):
        ufilter = lambda x: 'unit' not in x.attrs or x.attrs['unit'] == 'SKIP' or x.attrs['unit'] == ''
        with h5py.File(self.hdf5_path,'r') as hf:
            grp = hf[self.top_level_path]
            var_names = [(key,'') if ufilter(grp[key]) else (key,'({})'.format(grp[key].attrs['unit'])) 
                         for key in grp]
            footer = grp.attrs['footer']
            
        name_strs = '\n'.join(['{} {} -- [description]'.format(v[0],v[1]) for v in var_names])
        return '''{top_level_path}
        
Fields
------
{vars_and_units}

Footer
------
{footer}
        '''.format(top_level_path=self.top_level_path,vars_and_units=name_strs,footer=footer)


class IonBase(object):
    """
    Base class for accessing CHIANTI data attached to a particular ion
    """
    
    def __init__(self,ion_name,hdf5_path=None):
        self.ion_name = ion_name
        self.element = ion_name.split('_')[0]
        self.stage = ion_name.split('_')[-1]
        if hdf5_path is None:
            self.hdf5_path = fiasco.defaults['chianti_hdf5_dbase_root']
        else:
            self.hdf5_path = hdf5_path
       
    @property
    def abundance(self):
        return DataIndexer(self.hdf5_path,'/'.join([self.element,'abundance']))
        

def add_property(cls,filetype):
    """
    Dynamically add filetype properties to base data access class
    """
    def property_template(self):
        with h5py.File(self.hdf5_path,'r') as hf:
            if '/'.join([self.element,self.ion_name,filetype]) not in hf:
                return None
        return DataIndexer(self.hdf5_path,'/'.join([self.element,self.ion_name,filetype]))

    property_template.__doc__ = 'Data in {} type file'.format(filetype)
    property_template.__name__ = filetype
    setattr(cls,property_template.__name__,property(property_template))
    
# Collect the filetypes and add the methods
all_ext = [cls.filetype for cls in all_subclasses(GenericParser) 
           if hasattr(cls,'filetype') and cls.filetype not in ['abund']]
for filetype in all_ext:
    add_property(IonBase,filetype)