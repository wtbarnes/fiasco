"""
Source classes for CHIANTI filetypes
"""
import os
import numpy as np
import astropy.units as u
import fiasco

from .generic import GenericParser

__all__ = ['ElvlcParser', 'FblvlParser', 'ScupsParser',
           'PsplupsParser', 'EasplomParser', 'EasplupsParser',
           'WgfaParser', 'CilvlParser', 'ReclvlParser',
           'RrparamsParser', 'TrparamsParser', 'DrparamsParser',
           'DiparamsParser','AbundParser','IoneqParser']


class ElvlcParser(GenericParser):
    dtypes = [int,str,int,str,float,float,float]
    units = [None,None,None,None,u.dimensionless_unscaled,1/u.cm,1/u.cm]
    headings = ['level index','configuration','multiplicity',
                'orbital angular momentum','total angular momentum',
                'observed energy','theoretical energy']

    
class FblvlParser(GenericParser):
    dtypes = [int,str,int,int,str,int,float,float]
    units = [None,None,None,None,None,None,1/u.cm,1/u.cm]
    headings = ['level index','configuration','principal quantum number',
                'azimuthal quantum number','orbital angular momentum',
                'multiplicity','observed energy','theoretical energy']

    
class ScupsParser(GenericParser): 
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
            # main data
            table.append(line)
        else:
            # scaled temperature or collision strengths
            scaled = np.array(line,dtype=float)
            table[-1].append(scaled)


class PsplupsParser(GenericParser):
    dtypes = [int,int,int,float,float,float,float]
    units = [None,None,None,u.dimensionless_unscaled,u.Ry,u.dimensionless_unscaled,
             u.dimensionless_unscaled]
    headings = ['lower level index','upper level index','Burgess-Tully scaling type',
                'oscillator strength','delta energy','Burgess-Tully scaling parameter',
                'Burgess-Tully scaled effective collision strengths']
    
    def preprocessor(self,table,line,index):
        line = list(filter(None,('      '.join(line)).split()))
        row = line[:5]
        tmp = line[5]
        i_split = [i for i,char in enumerate(tmp) if char=='-' and tmp[i-1]!='e'][0]
        row += [tmp[:i_split]]
        scups = [tmp[i_split+1:]] 
        tmp_scups = line[6:]
        if len(tmp_scups[-1].split('-')) > 2:
            tmp = tmp_scups[-1]
            i_split = [i for i,char in enumerate(tmp) if char=='-' and tmp[i-1]!='e'][0]
            tmp_scups = tmp_scups[:-1] + [tmp[:i_split],tmp[i_split+1:]]
        scups = np.array(scups+tmp_scups,dtype=float)
        row += [scups]
        table.append(row)
        
            
class EasplomParser(GenericParser):
    dtypes = [int,int,int,float,float,float,float]
    units = [None,None,None,u.dimensionless_unscaled,u.Ry,u.dimensionless_unscaled,
             u.dimensionless_unscaled]
    headings = ['lower level index','upper level index','Burgess-Tully scaling type',
                'oscillator strength','delta energy','Burgess-Tully scaling parameter',
                'Burgess-Tully scaled cross-section']
    
    def preprocessor(self,table,line,index):
        line = list(filter(None,('      '.join(line)).split()))
        scaled_cs = np.array(line[8:],dtype=float)
        row = line[2:8] + [scaled_cs]
        table.append(row)
        
        
class EasplupsParser(EasplomParser):
    dtypes = [int,int,int,float,float,float,float]
    units = [None,None,None,u.dimensionless_unscaled,u.Ry,u.dimensionless_unscaled,
             u.dimensionless_unscaled]
    headings = ['lower level index','upper level index','Burgess-Tully scaling type',
                'oscillator strength','delta energy','upsilon coefficient',
                'excitation-autoionization rate coefficients']
    
    
class WgfaParser(GenericParser):
    dtypes = [np.int,np.int,np.float,np.float,np.float,
              str,np.int,str,np.float,
              str,np.int,str,np.float]
    units = [None,None,u.angstrom,u.dimensionless_unscaled,1/u.s,None,None,None,None,
             None,None,None,None]
    headings = ['lower level index','upper level index',
                'transition wavelength','oscillator strength','radiative decay rate',
                'lower level configuration','lower level multiplicity',
                'lower level orbital angular momentum',
                'lower level total angular momentum',
                'upper level configuration','upper level multiplicity',
                'upper level orbital angular momentum',
                'upper level total angular momentum']
    
    def preprocessor(self,table,line,index):
        ### lower ###
        tmp = line[-2].strip().split()
        del tmp[-1] # delete rogue dash
        tmp_pretty = tmp[-1]
        config = ' '.join(tmp[:-1])
        mult = tmp_pretty[0]
        orb = tmp_pretty[1]
        frac = tmp_pretty[2:]
        if len(frac) == 1:
            frac = frac[0]
        else:
            frac = float(frac.split('/')[0])/float(frac.split('/')[-1])
        lower = [config,mult,orb,frac] 
        ### upper ###
        tmp = line[-1].strip().split()
        tmp_pretty = tmp[-1]
        config = ' '.join(tmp[:-1])
        mult = tmp_pretty[0]
        orb = tmp_pretty[1]
        frac = tmp_pretty[2:]
        if len(frac) == 1:
            frac = frac[0]
        else:
            frac = float(frac.split('/')[0])/float(frac.split('/')[-1])
        upper = [config,mult,orb,frac] 
        ### recombine and assemble ###
        table.append(line[:-2] + lower + upper)
        
        
class CilvlParser(GenericParser):
    dtypes = [int,int,float,float]
    units = [None,None,u.K,(u.cm**3)/u.s]
    headings = ['lower level index','upper level index','temperature',
                'ionization rate coefficient']
    
    def preprocessor(self,table,line,index):
        line = (' '.join(line)).split()
        if index%2 == 0:
            row = line[2:4]
            temperature = 10.**np.array(line[4:],dtype=float)
            row += [temperature]
            table.append(row)
        else:
            rate_coefficient = np.array(line[4:],dtype=float)
            table[-1].append(rate_coefficient)

            
class ReclvlParser(CilvlParser):
    headings = ['lower level index','upper level index','temperature',
                'recombination rate coefficient']
    
    
class RrparamsParser(GenericParser):
    
    def preprocessor(self,table,line,index):
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
                table[0] += (' '.join(line[3:])).split()
            else:
                table[0] += (' '.join(line[2:])).split()
                
                
class TrparamsParser(GenericParser):
    dtypes = [float,float]
    units = [u.K,(u.cm**3)/u.s]
    headings = ['temperature','total recombination rate']
    
    def preprocessor(self,table,line,index):
        if index > 0:
            super().preprocessor(table,line,index)
            
            
class DrparamsParser(GenericParser):
    
    def preprocessor(self,table,line,index):
        raise NotImplementedError('Parser for .drparams files not yet implemented.')
        
        
class DiparamsParser(GenericParser):
    
    def preprocessor(self,table,line,index):
        raise NotImplementedError('Parser for .diparams files not yet implemented.')


class AbundParser(GenericParser):
    dtypes = [int,float,str]
    units = [None,u.dimensionless_unscaled,None]
    headings = ['atomic number','abundance relative to H','element']

    def __init__(self,abundance_filename):
        self.abundance_filename = abundance_filename
        self.full_path = os.path.join(fiasco.defaults['chianti_dbase_root'],
                                      'abundance', self.abundance_filename)

    def preprocessor(self,table,line,index):
        line[-1] = line[-1].strip()
        super().preprocessor(table,line,index)

    def postprocessor(self,df):
        df['abundance relative to H'] = 10.**(df['abundance relative to H'] 
                                              - df['abundance relative to H'][df['element']=='H'])
        df.meta['abundance_filename'] = self.abundance_filename
        return df


class IoneqParser(GenericParser):
    dtypes = [int,int,float,float]
    units = [None,None,u.K,u.dimensionless_unscaled]
    headings = ['atomic number','ion','temperature','ionization fraction']

    def __init__(self,ioneq_filename):
        self.ioneq_filename = ioneq_filename
        self.full_path = os.path.join(fiasco.defaults['chianti_dbase_root'],
                                      'ioneq', self.ioneq_filename)

    def preprocessor(self,table,line,index):
        if index == 0:
            pass
        elif index == 1:
            self.temperature = 10.**np.array(line,dtype=float)
        else:
            ioneq = np.array(line[2:],dtype=float)
            table.append(line[:2] + [self.temperature,ioneq])

    def postprocessor(self,df):
        df.meta['ioneq_filename'] = self.ioneq_filename
