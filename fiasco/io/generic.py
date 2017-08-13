"""
Base class for file parser
"""
import os
from astropy.table import QTable
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
            line = list(filter(None, line.strip().split('  ')))
            if line[0] == '-1':
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
        table.append(line)
        
    def postprocessor(self,df):
        """
        Default postprocessor method run on the whole dataframe
        """
        df.meta['element'] = self.element
        df.meta['ion'] = self.ion_name
        return df
    
    def to_hdf5(self,hf):
        """
        Add datasets to a group for an HDF5 file handler
        """
        pass