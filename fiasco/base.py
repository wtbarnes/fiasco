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
    
    def __init__(self,top_level_path,**kwargs):
        self.top_level_path = top_level_path
        self.hdf5_path = kwargs.get('hdf5_path',
                                    fiasco.defaults['chianti_hdf5_dbase_root'])
    
    def __getitem__(self,key):
        with h5py.File(self.hdf5_path,'r') as hf:
            grp = hf[self.top_level_path]
            if key not in grp:
                raise IndexError('{} not a valid dataset for {}'.format(key,self.top_level_path))
            ds = grp[key]
            if ds.attrs['unit'] == 'SKIP':
                data = np.array(ds,dtype=ds.dtype)
            else:
                data = u.Quantity(np.array(ds),ds.attrs['unit'],dtype=ds.dtype)
            if '|S' in data.dtype.str:
                data = data.astype(str)
        return data
    
    def __repr__(self):
        with h5py.File(self.hdf5_path,'r') as hf:
            grp = hf[self.top_level_path]
            var_names = [(key,'')
                         if grp[key].attrs['unit'] == 'SKIP' or grp[key].attrs['unit'] == ''
                         else (key,'({})'.format(grp[key].attrs['unit'])) 
                         for key in grp]
            footer = grp.attrs['footer']
            
        name_strs = '\n'.join(['{} {} -- [description]'.format(v[0],v[1]) 
                               for v in var_names])
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
    
    def __init__(self,ion_name):
        self.ion_name = ion_name
        self.element = ion_name.split('_')[0]
        self.stage = ion_name.split('_')[-1]
       
    @property
    def abundance(self):
        return DataIndexer('/'.join([self.element,'abundance']))
   
    @property
    def ionization_potential(self):
        return DataIndexer('/'.join([self.element,'ionization_potential']))
        

def add_property(cls,filetype):
    """
    Dynamically add filetype properties to base data access class
    """
    def property_template(self):
        return DataIndexer('/'.join([self.element,self.ion_name,filetype]))
    property_template.__doc__ = 'Data in {} type file'.format(filetype)
    property_template.__name__ = filetype
    setattr(cls,property_template.__name__,property(property_template))
    
# Collect the filetypes and add the methods
all_ext = [cls.filetype for cls in all_subclasses(GenericParser) 
           if hasattr(cls,'filetype') and cls.filetype not in ['abund','ip']]
for filetype in all_ext:
    add_property(IonBase,filetype)