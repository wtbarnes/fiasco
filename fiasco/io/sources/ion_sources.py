"""
Source classes for CHIANTI filetypes attached to ions
"""
import astropy.units as u
import fortranformat
import numpy as np

from fiasco.io.generic import GenericIonParser

__all__ = [
    'ElvlcParser',
    'FblvlParser',
    'ScupsParser',
    'SplupsParser',
    'PsplupsParser',
    'EasplomParser',
    'EasplupsParser',
    'WgfaParser',
    'CilvlParser',
    'ReclvlParser',
    'RrparamsParser',
    'TrparamsParser',
    'DrparamsParser',
    'DiparamsParser',
    'AutoParser',
    'RrlvlParser',
]


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


class SplupsParser(GenericIonParser):
    """
    Spline fits to scaled collisions strengths (denoted by upsilon) between energy levels as described
    in :cite:t:`burgess_analysis_1992`.  These files were used in CHIANTI versions prior to 8.0, and
    were replaced by ``.scups`` files in versions after that.

    Notes
    -----
    * The number of spline points for the rates depends on the fit type, 5 points for type 6
      fits and 9 points for type 2.
    """
    filetype = 'splups'
    dtypes = [int, int, int, int, int, float, float, float, 'object']
    units = [
        None,
        None,
        None,
        None,
        None,
        u.dimensionless_unscaled,
        u.Ry,
        u.dimensionless_unscaled,
        u.dimensionless_unscaled,
    ]
    headings = [
        'Z',
        'ion',
        'lower_level',
        'upper_level',
        'bt_type',
        'gf',
        'delta_energy',
        'bt_c',
        'bt_upsilon',
    ]
    descriptions = [
        'atomic number',
        'ionization state',
        'lower level index',
        'upper level index',
        'Burgess-Tully scaling type',
        'oscillator strength',
        'delta energy',
        'Burgess-Tully scaling parameter',
        'Burgess-Tully scaled effective collision strengths',
    ]

    def preprocessor(self, table, line, index):
        n_spline =  9  # Max number of spline points
        fformat = fortranformat.FortranRecordReader(f'(5I3,{3+n_spline}E10.3)')
        line = fformat.read(line)
        # NOTE: The first eight entries are fixed. The last entry is the scaled
        # spline fit to the array and can vary in length.
        # NOTE: Some spline fits only have 5 points and the scaling type is not
        # a reliable way to determine this so we have to filter these manually.
        # When fortranformat has missing entries, it fills them in as None. We
        # remove them here to avoid the undefined behavior of None in a ragged
        # array within an astropy Table.
        spline_fit = line[8:]
        spline_fit = [sf for sf in spline_fit if sf is not None]
        row = line[:8] + [np.array(spline_fit)]
        table.append(row)


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
        n_spline =  9  # Max number of spline points
        fformat = fortranformat.FortranRecordReader(f'(3I3,{3+n_spline}E10.3)')
        line = fformat.read(line)
        # NOTE: The first six entries are fixed. The last entry is the scaled
        # spline fit to the rate array and can vary in length.
        # NOTE: Some spline fits only have 5 points and the scaling type is not
        # a reliable way to determine this so we have to filter these manually.
        # When fortranformat has missing entries, it fills them in as None. We
        # remove them here to avoid the undefined behavior of None in a ragged
        # array within an astropy Table.
        spline_fit = line[6:]
        spline_fit = [sf for sf in spline_fit if sf is not None]
        row = line[:6] + [np.array(spline_fit)]
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
    dtypes = [int, int, 'object', 'object']
    headings = ['lower_level', 'upper_level', 'temperature', 'recombination_rate']
    descriptions = ['lower level index', 'upper level index', 'temperature', 'recombination rate coefficient']

    def postprocessor(self, df):
        # NOTE: For some versions of the database, not all of the temperatures and rates in a given file are
        # the same length. As such, the dtype of these columns must be set to "object". However, in cases
        # where they are all equal, we want to make sure they are preserved as floats.
        for cn in df.colnames:
            all_equal = np.all(np.array([row.size for row in df[cn]]) == df[cn][0].size)
            if df[cn].dtype == np.dtype('O') and all_equal:
                df[cn] = df[cn].astype(np.dtype('float64'))

        df = super().postprocessor(df)
        return df


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


class AutoParser(GenericIonParser):
    """
    Autoionization rates for each level in an ion.

    The autoionization rate is the rate of decay of atomic level through autoionization to a
    bound level. It is also needed to calculate the dielectronic recombination rate
    from the more highly ionized ions, by means of the principle of detailed-balance.
    For a full description of these files, see :cite:t:`dere_chianti_2017`.
    """
    filetype = 'auto'
    dtypes = [int, int, float, str, str]
    units = [None, None, 1/u.s, None, None]
    headings = [
        'lower_level',
        'upper_level',
        'autoionization_rate',
        'lower_label',
        'upper_label'
    ]
    descriptions = [
        'lower level index',
        'upper level index',
        'autoionization rate',
        'lower level label',
        'upper level label'
    ]
    fformat = fortranformat.FortranRecordReader('(2I7,E12.2,A30,A30)')

    def preprocessor(self, table, line, index):
        super().preprocessor(table, line, index)
        # remove the dash in the second-to-last entry
        table[-1][-2] = table[-1][-2].split('-')[0].strip()


class RrlvlParser(GenericIonParser):
    """
    Level-resolved recombination rates as a function of temperature.

    These files contain the *Direct* radiative recombination rates from
    the recombining ion to the recombined ion, whereas the ``.reclvl`` files
    contain the *effective* radiative recombination rates.
    A given ion should have either a ``.rrlvl`` file or a ``.reclvl`` file,
    but not both. For a full description of these files, see :cite:t:`young_chianti_2019`.
    """
    filetype = 'rrlvl'
    dtypes = [int, int, int, int, float, float]
    units = [None, None, None, None, u.K, (u.cm**3)/u.s]
    headings = [
        'Z',
        'ion',
        'initial_level',
        'final_level',
        'temperature',
        'rate',
    ]
    descriptions = [
        'atomic number',
        'ionization state',
        'level index of the recombining ion',
        'index of the final level of the transition in the recombined ion',
        'temperatures at which rates are tabulated',
        'direct radiative recombination rate coefficients',
    ]

    def preprocessor(self, table, line, index):
        # NOTE: Every pair of lines has the same first four entries. On the even
        # lines, the remaining entries contain the temperatures and on the odd
        # lines, the remaining entries are the rate coefficients. Thus, every
        # other line needs to be added to the one above it.
        line = line.strip().split()
        if index % 2 == 0:
            row = line[:4]
            temperature = np.array(line[4:], dtype=float)
            row += [temperature]
            table.append(row)
        else:
            rate_coefficient = np.array(line[4:], dtype=float)
            table[-1].append(rate_coefficient)
