"""
Source classes for CHIANTI filetypes attached to ions
"""
import os

import numpy as np
import h5py
import astropy.units as u
import fortranformat

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

    def postprocessor(self,df):
        for cn in df.colnames:
            all_equal = np.all(np.array([row.size for row in df[cn]]) == df[cn][0].size)
            if df[cn].dtype == np.dtype('O') and all_equal:
                df[cn] = df[cn].astype(np.dtype('float64'))

        df = super().postprocessor(df)
        return df


class PsplupsParser(ScupsParser):
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
        
        array = np.array(new_line[6:],dtype=np.dtype('float64'))
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
    dtypes = [int,int,float,float,float,str,str]
    units = [None,None,u.angstrom,u.dimensionless_unscaled,1/u.s,None,None]
    headings = ['lower level index','upper level index',
                'transition wavelength','oscillator strength','radiative decay rate',
                'lower level label','upper level label']
    fformat = fortranformat.FortranRecordReader('(2I5,F15.3,2E15.3,A30,A30)')
    
    def preprocessor(self,table,line,index):
        super().preprocessor(table,line,index)
        # remove the dash in the second-to-last entry
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
        line = line.strip().split()
        if index == 0:
            self._drparams_filetype = int(line[0])
            if self._drparams_filetype == 1:
                # Badnell type
                self.dtypes = [float,float]
                self.units = [u.K,(u.cm**3)/u.s*(u.K**(3/2))]
                self.headings = ['E fit parameter','c fit parameter']
            elif self._drparams_filetype == 2:
                # Shull type
                self.dtypes = [float,float,float,float]
                self.units = [(u.cm**3)/u.s*(u.K**(3/2)),u.dimensionless_unscaled,u.K,u.K]
                self.headings = ['A fit coefficient','B fit coefficient',
                                 'T0 fit coefficient','T1 fit coefficient']
            else:
                raise ValueError('Unrecognized drparams filetype {}'.format(self._drparams_filetype))
        else:
            if self._drparams_filetype == 1:
                tmp = np.array(line[2:],dtype=float)
                if index%2 == 0:
                    table[-1].append(tmp)
                else:
                    table.append([tmp])
            else:
                table.append(line[2:])

        
class DiparamsParser(GenericParser):
    filetype = 'diparams'
    dtypes = [float,float,float,float]
    units = [1/u.cm,u.dimensionless_unscaled,u.dimensionless_unscaled,
             u.dimensionless_unscaled]
    headings = ['ionization potential','Burgess-Tully scaling factor',
                'scaled energy','scaled cross-section']
    
    def preprocessor(self,table,line,index):
        tmp = line.strip().split()
        if index == 0:
            self._num_fits = int(tmp[2])
            self._num_lines = int(tmp[3])
            self._has_excitation_autoionization = bool(int(tmp[4]))
        elif index==self._num_lines*2 + 1 and self._has_excitation_autoionization:
            self._excitation_autoionization = float(tmp[0])
        elif index%2 != 0:
            bt_factor = tmp[0]
            u_spline = np.array(tmp[1:],dtype=float)
            table.append([bt_factor,u_spline])
        else:
            ionization_potential = tmp[0]
            cs_spline = np.array(tmp[1:],dtype=float)
            table[-1] = [ionization_potential]+table[-1]+[cs_spline]
            
    def postprocessor(self,df):
        if hasattr(self,'_excitation_autoionization'):
            df['excitation autoionization'] = len(df['ionization potential'])*[self._excitation_autoionization]
        df = super().postprocessor(df)
        return df
