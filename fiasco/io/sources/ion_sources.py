"""
Source classes for CHIANTI filetypes attached to ions
"""
import astropy.units as u
import fortranformat
import numpy as np

from fiasco.io.generic import GenericIonParser

__all__ = ['ElvlcParser', 'FblvlParser', 'ScupsParser',
           'PsplupsParser', 'EasplomParser', 'EasplupsParser',
           'WgfaParser', 'CilvlParser', 'ReclvlParser',
           'RrparamsParser', 'TrparamsParser', 'DrparamsParser',
           'DiparamsParser']


class ElvlcParser(GenericIonParser):
    """
    Energy levels and configurations for each level in an ion.
    """
    filetype = 'elvlc'
    dtypes = [int, str, str, int, str, float, float, float]
    units = [None, None, None, u.dimensionless_unscaled, None, u.dimensionless_unscaled, 1/u.cm, 1/u.cm]
    headings = ['level', 'config', 'label', 'multiplicity', 'L_label', 'J', 'E_obs', 'E_th']
    descriptions = [
        'level index',
        'configuration',
        'level label',
        'multiplicity, 2s+1',
        'orbital angular momentum',
        'total angular momentum',
        'observed energy',
        'theoretical energy'
    ]
    fformat = fortranformat.FortranRecordReader('(I7,A30,A5,I5,A5,F5.1,F15.3,F15.3)')


class FblvlParser(GenericIonParser):
    """
    Energy levels and configuration related to the calculation of the free-bound
    continuum. Only available for those ions and levels for which a free-bound
    continuum can be calculated.
    """
    filetype = 'fblvl'
    dtypes = [int, str, int, int, str, int, float, float]
    units = [None, None, None, None, None, None, 1/u.cm, 1/u.cm]
    headings = ['level', 'config', 'n', 'L', 'L_label', 'multiplicity', 'E_obs', 'E_th']
    descriptions = [
        'level index',
        'configuration',
        'principal quantum number',
        'azimuthal quantum number',
        'orbital angular momentum',
        'multiplicity, 2s+1',
        'observed energy',
        'theoretical energy'
    ]
    fformat = fortranformat.FortranRecordReader('(I5,A20,2I5,A3,I5,2F20.3)')


class ScupsParser(GenericIonParser):
    """
    Scaled collisions strengths (denoted by upsilon) between energy levels as described
    in :cite:t:`burgess_analysis_1992`.
    """
    filetype = 'scups'
    dtypes = [int, int, float, float, float, int, int, float, 'object', 'object']
    units = [
        None,
        None,
        u.Ry,
        u.dimensionless_unscaled,
        1/u.Ry,
        None,
        None,
        u.dimensionless_unscaled,
        u.dimensionless_unscaled,
        u.dimensionless_unscaled
    ]
    headings = [
        'lower_level',
        'upper_level',
        'delta_energy',
        'gf',
        'high_t_limit',
        'n_t',
        'bt_type',
        'bt_c',
        'bt_t',
        'bt_upsilon'
    ]
    descriptions = [
        'lower level index',
        'upper level index',
        'delta energy',
        'oscillator strength',
        'high-temperature limit',
        'number of scaled temperatures',
        'Burgess-Tully scaling type',
        'Burgess-Tully scaling parameter',
        'Burgess-Tully scaled temperatures',
        'Burgess-Tully scaled effective collision strengths'
    ]

    def preprocessor(self, table, line, index):
        if index % 3 == 0:
            super().preprocessor(table, line, index)
        else:
            # scaled temperature or collision strengths
            scaled = np.array(line.strip().split(), dtype=float)
            table[-1].append(scaled)

    def postprocessor(self, df):
        for cn in df.colnames:
            all_equal = np.all(np.array([row.size for row in df[cn]]) == df[cn][0].size)
            if df[cn].dtype == np.dtype('O') and all_equal:
                df[cn] = df[cn].astype(np.dtype('float64'))

        df = super().postprocessor(df)
        return df


class PsplupsParser(ScupsParser):
    """
    Spline fits to scaled collision rates for protons. These files are discussed in
    section 2.2 of :cite:t:`young_chianti-atomic_2003` and the details of how these
    quantities are scaled are given in :cite:t:`burgess_analysis_1992`.

    Notes
    -----
    * Unlike the electron "scups" and "splups" files which contain the collision strengths
      (upsilons), these files contain the scaled *rates*.
    * The number of spline points for the rates depends on the fit type, 5 points for type 6
      fits and 9 points for type 2.
    """
    filetype = 'psplups'
    dtypes = [int, int, int, float, float, float, 'object']
    units = [
        None,
        None,
        None,
        u.dimensionless_unscaled,
        u.Ry,
        u.dimensionless_unscaled,
        u.dimensionless_unscaled
    ]
    headings = ['lower_level', 'upper_level', 'bt_type', 'gf', 'delta_energy', 'bt_c', 'bt_rate']
    descriptions = [
        'lower level index',
        'upper level index',
        'Burgess-Tully scaling type',
        'oscillator strength',
        'delta energy',
        'Burgess-Tully scaling parameter',
        'Burgess-Tully scaled collision rate'
    ]

    def preprocessor(self, table, line, index):
        tmp = line.strip().split()
        # 5-point fit for type 6, 9-point fit for type 2
        n_spline = 5 if int(tmp[2]) == 6 else 9
        fformat = fortranformat.FortranRecordReader(f'(3I3,{3+n_spline}E10.3)')
        line = fformat.read(line)
        row = line[:6] + [np.array(line[6:])]
        table.append(row)


class EasplomParser(GenericIonParser):
    """
    Spline fits to the excitation-autoionization scaled cross-sections.
    See :cite:t:`burgess_analysis_1992` and :cite:t:`dere_ionization_2007`
    for more details.
    """
    filetype = 'easplom'
    dtypes = [int, int, int, float, float, float, float]
    units = [None, None, None, u.dimensionless_unscaled, u.Ry, u.dimensionless_unscaled, u.dimensionless_unscaled]
    headings = ['lower_level', 'upper_level', 'bt_type', 'gf', 'delta_energy', 'bt_c', 'bt_cross_section']
    descriptions = [
        'lower level index',
        'upper level index',
        'Burgess-Tully scaling type',
        'oscillator strength',
        'delta energy',
        'Burgess-Tully scaling parameter',
        'Burgess-Tully scaled cross-section'
    ]

    def preprocessor(self, table, line, index):
        line = line.strip().split()
        scaled_cs = np.array(line[8:], dtype=float)
        row = line[2:8] + [scaled_cs]
        table.append(row)


class EasplupsParser(EasplomParser):
    """
    Scaled collision strengths for calculating ionization rates due to excitation autoionization.
    """
    filetype = 'easplups'
    dtypes = [int, int, int, float, float, float, float]
    units = [None, None, None, u.dimensionless_unscaled, u.Ry, u.dimensionless_unscaled, u.dimensionless_unscaled]
    headings = ['lower_level', 'upper_level', 'bt_type', 'gf', 'delta_energy', 'bt_c', 'bt_upsilon']
    descriptions = [
        'lower level index',
        'upper level index',
        'Burgess-Tully scaling type',
        'oscillator strength',
        'delta energy',
        'upsilon coefficient',
        'Burgess-Tully scaled effective collision strength'
    ]


class WgfaParser(GenericIonParser):
    """
    Information about each possible transition in an ion, including level indices, wavelengths,
    and decay rates.
    """
    filetype = 'wgfa'
    dtypes = [int, int, float, float, float, str, str]
    units = [None, None, u.angstrom, u.dimensionless_unscaled, 1/u.s, None, None]
    headings = ['lower_level', 'upper_level', 'wavelength', 'gf', 'A', 'lower_label', 'upper_label']
    descriptions = [
        'lower level index',
        'upper level index',
        'transition wavelength',
        'oscillator strength',
        'radiative decay rate',
        'lower level label',
        'upper level label'
    ]
    fformat = fortranformat.FortranRecordReader('(2I5,F15.3,2E15.3,A30,A30)')

    def preprocessor(self, table, line, index):
        super().preprocessor(table, line, index)
        # remove the dash in the second-to-last entry
        table[-1][-2] = table[-1][-2].split('-')[0].strip()


class CilvlParser(GenericIonParser):
    filetype = 'cilvl'
    dtypes = [int, int, float, float]
    units = [None, None, u.K, (u.cm**3)/u.s]
    headings = ['lower_level', 'upper_level', 'temperature', 'ionization_rate']
    descriptions = ['lower level index', 'upper level index', 'temperature', 'ionization rate coefficient']

    def preprocessor(self, table, line, index):
        line = line.strip().split()
        if index % 2 == 0:
            row = line[2:4]
            temperature = 10.**np.array(line[4:], dtype=float)
            row += [temperature]
            table.append(row)
        else:
            rate_coefficient = np.array(line[4:], dtype=float)
            table[-1].append(rate_coefficient)


class ReclvlParser(CilvlParser):
    filetype = 'reclvl'
    headings = ['lower_level', 'upper_level', 'temperature', 'recombination_rate']
    descriptions = ['lower level index', 'upper level index', 'temperature', 'recombination rate coefficient']


class RrparamsParser(GenericIonParser):
    """
    Fit parameters for calculating radiative recombination rates. The first two fit types are
    given in Eqs. 1 and 2 of :cite:t:`badnell_radiative_2006` and the third fit type is given
    by Eq. 4 of :cite:t:`shull_ionization_1982`.
    """
    filetype = 'rrparams'

    def preprocessor(self, table, line, index):
        line = line.strip().split()
        if index == 0:
            filetype = int(line[0])
            table.append([filetype])
            if filetype == 1:
                self.dtypes = [int, float, float, float, float]
                self.units = [None, (u.cm**3)/u.s, None, u.K, u.K]
                self.headings = ['fit_type', 'A_fit', 'B_fit', 'T0_fit', 'T1_fit']
                self.descriptions = [
                    'fit type',
                    'A fit parameter',
                    'B fit parameter',
                    'T0 fit parameter',
                    'T1 fit parameter'
                ]
            elif filetype == 2:
                self.dtypes = [int, float, float, float, float, float, float]
                self.units = [None, (u.cm**3)/u.s, None, u.K, u.K, None, u.K]
                self.headings = ['fit_type', 'A_fit', 'B_fit', 'T0_fit', 'T1_fit', 'C_fit', 'T2_fit']
                self.descriptions = [
                    'fit type',
                    'A fit parameter',
                    'B fit parameter',
                    'T0 fit parameter',
                    'T1 fit parameter',
                    'C fit parameter',
                    'T2 fit parameter'
                ]
            elif filetype == 3:
                self.dtypes = [int, float, float]
                self.units = [None, (u.cm**3)/u.s, None]
                self.headings = ['fit_type', 'A_fit', 'eta_fit']
                self.descriptions = ['fit type', 'A rad fit parameter', 'eta fit parameter']
            else:
                raise ValueError(f'Unrecognized .rrparams filetype {filetype}')
        else:
            if table[0][0] == 1 or table[0][0] == 2:
                table[0] += line[3:]
            else:
                table[0] += line[2:]


class TrparamsParser(GenericIonParser):
    filetype = 'trparams'
    dtypes = [float, float]
    units = [u.K, (u.cm**3)/u.s]
    headings = ['temperature', 'recombination_rate']
    descriptions = ['temperature', 'total recombination rate']

    def preprocessor(self, table, line, index):
        if index > 0:
            super().preprocessor(table, line, index)


class DrparamsParser(GenericIonParser):
    """
    Fit parameters for calculating dielectronic recombination. The first fit type is given by Eq. 3
    of :cite:t:`zatsarinny_dielectronic_2003` and the second fit type is given by Eq. 5 of
    :cite:t:`shull_ionization_1982`.
    """
    filetype = 'drparams'

    def preprocessor(self, table, line, index):
        line = line.strip().split()
        if index == 0:
            self._drparams_filetype = int(line[0])
            if self._drparams_filetype == 1:
                # Badnell type
                self.dtypes = [int, float, float]
                self.units = [None, u.K, (u.cm**3)/u.s*(u.K**(3/2))]
                self.headings = ['fit_type', 'E_fit', 'c_fit']
                self.descriptions = ['fit type', 'E fit parameter', 'c fit parameter']
            elif self._drparams_filetype == 2:
                # Shull type
                self.dtypes = [int, float, float, float, float]
                self.units = [None, (u.cm**3)/u.s*(u.K**(3/2)), u.dimensionless_unscaled, u.K, u.K]
                self.headings = ['fit_type', 'A_fit', 'B_fit', 'T0_fit', 'T1_fit']
                self.descriptions = [
                    'fit type',
                    'A fit coefficient',
                    'B fit coefficient',
                    'T0 fit coefficient',
                    'T1 fit coefficient'
                ]
            else:
                raise ValueError(f'Unrecognized drparams filetype {self._drparams_filetype}')
        else:
            if self._drparams_filetype == 1:
                tmp = np.array(line[2:], dtype=float)
                if index % 2 == 0:
                    tmp_col = table[-1]
                    for i in range(tmp.shape[0]):
                        table.append([self._drparams_filetype, tmp_col[i], tmp[i]])
                    del table[0]
                else:
                    table.append(tmp)
            else:
                table.append([self._drparams_filetype]+line[2:])


class DiparamsParser(GenericIonParser):
    """
    Scaled cross-sections for calculating the ionization rate due to direct ionization.
    See :cite:t:`burgess_analysis_1992` and :cite:t:`dere_ionization_2007` for more details.

    Notes
    -----
    - The scaled cross-sections date have been multiplied by :math:`10^{14}`
    """
    filetype = 'diparams'
    dtypes = [float, float, float, float, float]
    units = [u.eV, u.dimensionless_unscaled, u.dimensionless_unscaled, u.cm**2*u.eV**2, None]
    headings = ['ip', 'bt_c', 'bt_e', 'bt_cross_section', 'ea']
    descriptions = [
        'ionization potential',
        'Burgess-Tully scaling factor',
        'Burgess-Tully scaled energy',
        'Burgess-Tully scaled cross-section',
        'excitation autoionization'
    ]

    def preprocessor(self, table, line, index):
        tmp = line.strip().split()
        if index == 0:
            self._num_fits = int(tmp[2])
            self._num_lines = int(tmp[3])
            self._has_excitation_autoionization = bool(int(tmp[4]))
        elif index == self._num_lines*2 + 1 and self._has_excitation_autoionization:
            for t in table:
                t[-1] = float(tmp[0])
        elif index % 2 != 0:
            bt_factor = tmp[0]
            u_spline = np.array(tmp[1:], dtype=float)
            table.append([bt_factor, u_spline])
        else:
            ionization_potential = tmp[0]
            cs_spline = np.array(tmp[1:], dtype=float)*1e-14
            table[-1] = [ionization_potential] + table[-1] + [cs_spline] + [0.0]
