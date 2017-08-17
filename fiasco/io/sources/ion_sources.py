"""
Source classes for CHIANTI filetypes attached to ions
"""
import os
import numpy as np
import h5py
import astropy.units as u
import fortranformat
import fiasco

from ..generic import GenericParser

__all__ = ['ElvlcParser', 'FblvlParser', 'ScupsParser',
           'PsplupsParser', 'EasplomParser', 'EasplupsParser',
           'WgfaParser', 'CilvlParser', 'ReclvlParser',
           'RrparamsParser', 'TrparamsParser', 'DrparamsParser',
           'DiparamsParser']


class ElvlcParser(GenericParser):
    filetype = 'elvlc'
    dtypes = [int,str,str,int,str,float,float,float]
    units = [None,None,None,None,None,u.dimensionless_unscaled,1/u.cm,1/u.cm]
    headings = ['level index','configuration','level label','multiplicity',
                'orbital angular momentum','total angular momentum',
                'observed energy','theoretical energy']
    fformat = fortranformat.FortranRecordReader('(I7,A30,A5,I5,A5,F5.1,F15.3,F15.3)')

    
class FblvlParser(GenericParser):
    filetype = 'fblvl'
    dtypes = [int,str,int,int,str,int,float,float]
    units = [None,None,None,None,None,None,1/u.cm,1/u.cm]
    headings = ['level index','configuration','principal quantum number',
                'azimuthal quantum number','orbital angular momentum',
                'multiplicity','observed energy','theoretical energy']
    fformat = fortranformat.FortranRecordReader('(I5,A20,2I5,A3,I5,2F20.3)')

    
class ScupsParser(GenericParser): 
    filetype = 'scups'
    dtypes = [int,int,float,float,float,int,int,float,'object','object']
    units = [None,None,u.Ry,u.dimensionless_unscaled,1/u.Ry,None,None,
             u.dimensionless_unscaled,u.dimensionless_unscaled,
             u.dimensionless_unscaled]
    headings = ['lower level index','upper level index','delta energy','oscillator strength',
                'high-temperature limit','number of scaled temperatures','Burgess-Tully scaling type',
                'Burgess-Tully scaling parameter','Burgess-Tully scaled temperatures',
                'Burgess-Tully scaled effective collision strengths']
    
    def preprocessor(self,table,line,index):
        if index%3 == 0:
            super().preprocessor(table,line,index)
        else:
            # scaled temperature or collision strengths
            scaled = np.array(line.strip().split(),dtype=float)
            table[-1].append(scaled)

    def _write_to_hdf5(self,grp,name,data):
        if data.dtype is np.dtype('O'):
            ragged_dtype = h5py.special_dtype(vlen=np.dtype(data[0].dtype.str))
            return grp.create_dataset(name,data=data,dtype=ragged_dtype)
        else:
            return grp.create_dataset(name,data=data,dtype=data.dtype)


class PsplupsParser(GenericParser):
    filetype = 'psplups'
    dtypes = [int,int,int,float,float,float,'object']
    units = [None,None,None,u.dimensionless_unscaled,u.Ry,u.dimensionless_unscaled,
             u.dimensionless_unscaled]
    headings = ['lower level index','upper level index','Burgess-Tully scaling type',
                'oscillator strength','delta energy','Burgess-Tully scaling parameter',
                'Burgess-Tully scaled effective collision strengths']
    
    def preprocessor(self,table,line,index):
        line = line.strip().split()
        new_line = []
        for i,item in enumerate(line):
            dot_splits = [j for j,char in enumerate(item) if char=='.']
            if len(dot_splits) > 1:
                dot_splits = np.hstack([dot_splits,len(item)+3])
                split_items = [item[dot_splits[j]-1:dot_splits[j+1]-2] 
                               for j,_ in enumerate(dot_splits[:-1])]
                new_line += split_items
            else:
                new_line.append(item)
        
        array = np.array(new_line[6:],dtype=float)
        new_line = new_line[:6] + [array]
        table.append(new_line)
        
            
class EasplomParser(GenericParser):
    filetype = 'easplom'
    dtypes = [int,int,int,float,float,float,float]
    units = [None,None,None,u.dimensionless_unscaled,u.Ry,u.dimensionless_unscaled,
             u.dimensionless_unscaled]
    headings = ['lower level index','upper level index','Burgess-Tully scaling type',
                'oscillator strength','delta energy','Burgess-Tully scaling parameter',
                'Burgess-Tully scaled cross-section']
    
    def preprocessor(self,table,line,index):
        line = line.strip().split()
        scaled_cs = np.array(line[8:],dtype=float)
        row = line[2:8] + [scaled_cs]
        table.append(row)
        
        
class EasplupsParser(EasplomParser):
    filetype = 'easplups'
    dtypes = [int,int,int,float,float,float,float]
    units = [None,None,None,u.dimensionless_unscaled,u.Ry,u.dimensionless_unscaled,
             u.dimensionless_unscaled]
    headings = ['lower level index','upper level index','Burgess-Tully scaling type',
                'oscillator strength','delta energy','upsilon coefficient',
                'excitation-autoionization rate coefficients']
    
    
class WgfaParser(GenericParser):
    filetype = 'wgfa'
    dtypes = [np.int,np.int,np.float,np.float,np.float,str,str]
    units = [None,None,u.angstrom,u.dimensionless_unscaled,1/u.s,None,None]
    headings = ['lower level index','upper level index',
                'transition wavelength','oscillator strength','radiative decay rate',
                'lower level label','upper level label']
    fformat = fortranformat.FortranRecordReader('(2I5,F15.3,2E15.3,A30,A30)')
    
    def preprocessor(self,table,line,index):
        super().preprocessor(table,line,index)
        table[-1][-2] = table[-1][-2].split('-')[0].strip() 
        
        
class CilvlParser(GenericParser):
    filetype = 'cilvl'
    dtypes = [int,int,float,float]
    units = [None,None,u.K,(u.cm**3)/u.s]
    headings = ['lower level index','upper level index','temperature',
                'ionization rate coefficient']
    
    def preprocessor(self,table,line,index):
        line = line.strip().split()
        if index%2 == 0:
            row = line[2:4]
            temperature = 10.**np.array(line[4:],dtype=float)
            row += [temperature]
            table.append(row)
        else:
            rate_coefficient = np.array(line[4:],dtype=float)
            table[-1].append(rate_coefficient)

            
class ReclvlParser(CilvlParser):
    filetype = 'reclvl'
    headings = ['lower level index','upper level index','temperature',
                'recombination rate coefficient']
    
    
class RrparamsParser(GenericParser):
    filetype = 'rrparams'
    
    def preprocessor(self,table,line,index):
        line = line.strip().split()
        if index == 0:
            filetype = int(line[0])
            table.append([filetype])
            if filetype == 1:
                self.dtypes = [int,float,float,float,float]
                self.units = [None,(u.cm**3)/u.s,None,u.K,u.K]
                self.headings = ['fit type','A fit parameter','B fit parameter',
                                 'T0 fit parameter','T1 fit parameter']
            elif filetype == 2:
                self.dtypes = [int,float,float,float,float,float,float]
                self.units = [None,(u.cm**3)/u.s,None,u.K,u.K,None,u.K]
                self.headings = ['fit type','A fit parameter','B fit parameter',
                                 'T0 fit parameter','T1 fit parameter','C fit parameter',
                                 'T2 fit parameter']
            elif filetype == 3:
                self.dtypes = [int,float,float]
                self.units = [None,(u.cm**3)/u.s,None]
                self.headings = ['fit type','A rad fit parameter','eta fit parameter']
            else:
                raise ValueError('Unrecognized .rrparams filetype {}'.format(filetype))
        else:
            if table[0][0] == 1 or table[0][0] == 2:
                table[0] += line[3:]
            else:
                table[0] += line[2:]
                
                
class TrparamsParser(GenericParser):
    filetype = 'trparams'
    dtypes = [float,float]
    units = [u.K,(u.cm**3)/u.s]
    headings = ['temperature','total recombination rate']
    
    def preprocessor(self,table,line,index):
        if index > 0:
            super().preprocessor(table,line,index)
            
            
class DrparamsParser(GenericParser):
    filetype = 'drparams'
      
    def preprocessor(self,table,line,index):
        raise NotImplementedError('Parser for .drparams files not yet implemented.')
        
        
class DiparamsParser(GenericParser):
    filetype = 'diparams'
    
    def preprocessor(self,table,line,index):
        raise NotImplementedError('Parser for .diparams files not yet implemented.')
