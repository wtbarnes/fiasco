"""
Source classes for CHIANTI filetypes attached to ions
"""
import astropy.units as u
import fortranformat
import numpy as np

from packaging.version import Version

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
    """
    Total recombination rates as a function of temperature.
    """
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
    See :cite:t:`dere_ionization_2007` and :cite:t:`young_chianti_2025` for more details
    regarding the format of these files.

    .. note::

        The scaled cross-sections date have been multiplied by :math:`10^{14}`
    """
    filetype = 'diparams'
    dtypes = [float, float, float, float, float]
    units = [u.eV, u.dimensionless_unscaled, u.dimensionless_unscaled, u.cm**2*u.eV**2, u.dimensionless_unscaled]
    headings = ['ip', 'bt_c', 'bt_e', 'bt_cross_section', 'ea']
    descriptions = [
        'ionization potential',
        'Burgess-Tully scaling factor',
        'Burgess-Tully scaled energy',
        'Burgess-Tully scaled cross-section',
        'excitation autoionization scaling factor'
    ]

    def preprocessor(self, table, line, index):
        # NOTE: The order of these conditionals is important as the table is being modified in place
        if index == 0:
            # The first line has information about the number of fit points and number of transitions
            # included as well as the number of EA scaling coefficients.
            fformat = fortranformat.FortranRecordReader('(5I5)')
            _, _, self._n_fits, self._n_lines, self._n_ea = fformat.read(line)
        elif index == self._n_lines*2 + 1:
            # The last line is a standalone parameter with a length dependent on the number
            # of transitions included in the EA files. This is not necessarily present in all
            # files. In that case, this conditional is never met.
            ea_scaling = fortranformat.FortranRecordReader(f'({self._n_ea}E12.3)').read(line)
            for row in table:
                row[-1] = np.array(ea_scaling)
        elif index % 2 == 1:
            # The odd-numbered lines contain the scaling factor and energy array
            tmp = fortranformat.FortranRecordReader(f'({1+self._n_fits}F10.5)').read(line)
            table.append([tmp[0], np.array(tmp[1:])])
        else:
            # The even-numbered lines contain the ionization potential and the
            # scaled cross-section.
            tmp = fortranformat.FortranRecordReader(f'({1+self._n_fits}F10.5)').read(line)
            table[-1] = [tmp[0]] + table[-1] + [np.array(tmp[1:])*1e-14] + [1.0]


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


class DilvlParser(GenericIonParser):
    """
    Level-resolved direct ionization rates as a function of temperature.

    These files contain the direct ionization rates from the ion being ionized
    to the ionized ion. As an example, the ``c_2.dilvl`` file contains direct ionization
    rates connecting the levels in C II, given by the ``level_1``, to the levels in C III,
    given by the ``level_2`` column. For more information about these files, see
    :cite:t:`dufresne_chianti_2024-1`.
    """
    filetype = 'dilvl'
    dtypes = [int, int, float, float, float]
    units = [None, None, 'eV', 'cm3 s-1', 'K']
    headings = [
        'level_1',
        'level_2',
        'ip',
        'rate',
        'temperature',
    ]
    descriptions = [
        'initial level of the ion before ionization',
        'final level of the ion after ionization',
        'energy input required for the transition to take place',
        'rate coefficients for the transition',
        'temperatures at which the rate coefficients are calculated'
    ]

    def preprocessor(self, table, line, index):
        if index == 0:
            self._temperature = np.array(line.strip().split(), dtype=float)
        else:
            line = line.strip().split()
            rate = np.array(line[3:])
            temperature = np.array(self._temperature, dtype=rate.dtype)
            table.append(line[:3]+[rate,temperature])


class EalvlParser(DilvlParser):
    """
    Level-resolved indirect ionization rates as a function of temperature.

    These files contain the indirect ionization (also known as excitation autoionization)
    rates from the ion being ionized to the ionized ion. The format is identical to those
    files parsed by `DilvlParser`. For more information about these files, see
    :cite:t:`dufresne_chianti_2024-1`.
    """
    filetype = 'ealvl'


class RrcoeffsParser(GenericIonParser):
    """
    Level-resolved radiative recombination rate fitting coefficients as a function of temperature.

    These files contain the fitting coefficients for the level-resolved radiative recombination rates
    as a function of temperature. The coefficients in these files are analogous to those data in the files
    parsed by `RrparamsParser`. For more information about these files, see
    :cite:t:`dufresne_chianti_2024-1`.
    """
    filetype = 'rrcoeffs'
    dtypes = 3*[int] + 6*[float]
    units = [None, None, None, 'cm3 s-1', None, 'K', 'K', None, 'K']
    headings = [
        'level',
        'weight',
        'fit_type',
        'A_fit',
        'B_fit',
        'T0_fit',
        'T1_fit',
        'C_fit',
        'T2_fit',
    ]
    descriptions = [
        'Initial energy level of the ion prior to recombination',
        'Statistical weight (2J+1) of the initial level',
        'Type of fitting formula used',
        'A fit coefficient',
        'B fit coefficient',
        'T0 fit coefficient',
        'T1 fit coefficient',
        'C fit coefficient',
        'T2 fit coefficient',
    ]

    def preprocessor(self, table, line, index):
        if index == 0:
            return
        # NOTE: Format is dependent on fit type and there can be multiple fit types in a single
        # file. As such, we have to pull out the fit type first and adjust the format appropriately.
        # A type 1 fit has fewer parameters so we just pad those extra columns with NaN.
        _line = line.strip().split()
        fit_type = int(_line[2])
        if fit_type == 1:
            fformat = fortranformat.FortranRecordReader('(3I5,E12.4,F10.5,2E12.4)')
            line = fformat.read(line)
            line += 2*[np.nan]
        elif fit_type == 2:
            fformat = fortranformat.FortranRecordReader('(3I5,E12.4,F10.5,2E11.4,F10.5,E12.4)')
            line = fformat.read(line)
        else:
            raise ValueError(f'Unrecognized fit type {fit_type} for rrcoeffs file.')
        table.append(line)


class DrcoeffsParser(GenericIonParser):
    """
    Level-resolved dielectronic recombination rate fitting coefficients as a function of temperature.

    These files contain the fitting coefficients for level-resolved dielectronic recombination rates
    as a function of temperature. The coefficients in these files are analogous to the type 1 fitting
    coefficients in the files parsed by `DrparamsParser`. For more information about these files, see
    :cite:t:`dufresne_chianti_2024-1`.
    """
    filetype = 'drcoeffs'
    dtypes = [int, int, int, float, float]
    units = [None, None, None, 'K', u.cm**3/u.s*u.K**(3/2)]
    headings = [
        'level',
        'weight',
        'fit_type',
        'E_fit',
        'C_fit',
    ]
    descriptions = [
        'Initial energy level of the ion prior to recombination',
        'Statistical weight (2J+1) of the initial level',
        'Type of fitting formula used',
        'E fit parameter',
        'C fit parameter',
    ]

    def preprocessor(self, table, line, index):
        if self.chianti_version<=Version('11.0.2') and self.ion_name=='he_2':
            # Prior to and including v11.0.2, the he_2.drcoeffs file includes an erroneous "1"
            # in the second line. As such, we just need to skip this line. If this continues to persist
            # in future versions of the database, the version number needs to be incremented.
            index -= 1
        if index <= 0:
            return
        fformat = fortranformat.FortranRecordReader('(3I5,9E12.4)')
        line = fformat.read(line)
        # NOTE: This conditional is because the lines of this file come in pairs, with
        # the first three entries being repeated on each pair of lines and the remaining
        # entries being the different arrays of fitting coefficients.
        if (index-1)%2==0:
            table.append(line[:3]+[line[3:]])
        else:
            table[-1] += [line[3:]]


class CtilvlParser(GenericIonParser):
    """
    Level-resolved charge transfer ionization rate coefficients as a function of temperature.

    These files contain the rate coefficients for charge transfer ionization as a function
    of temperature between the levels of the ion before ionization and the ion after ionization.
    For more information about these files, see :cite:t:`dufresne_chianti_2024-1`.
    """
    filetype = 'ctilvl'
    dtypes = 4*[int] + 3*[float]
    units = 4*[None] + ['eV', 'cm3 s-1', 'K']
    headings = [
        'level_1',
        'level_2',
        'Z_perturber',
        'N_e_perturber',
        'delta_energy',
        'rate',
        'temperature',
    ]
    descriptions = [
        'Initial level of the ion before ionization.',
        'Final level of the ion after ionization.',
        'Atomic number of the perturber involved in the transition.',
        'Number of electrons in perturber before transition.',
        'Input energy required for the transition to take place.',
        'Ionization rate coefficients for the transition.',
        'Temperatures at which the rate coefficients are evaluated.'
    ]

    def preprocessor(self, table, line, index):
        if index == 0:
            self._temperature = np.array(line.strip().split(), dtype=self.dtypes[-1])
        else:
            n_T = self._temperature.shape[0]
            fformat = fortranformat.FortranRecordReader(f'(4I5,E11.3,{n_T}E11.3)')
            line = fformat.read(line)
            table.append(line[:5]+[line[5:]]+[self._temperature])


class CtrlvlParser(CtilvlParser):
    """
    Level-resolved charge transfer recombination rate coefficients as a function of temperature.

    These files contain the rate coefficients for charge transfer recombination as a function
    of temperature between the levels of the recombining and recombined ion. The format is identical
    to those parsed by `CtilvlParser`. For more information about these files, see
    :cite:t:`dufresne_chianti_2024-1`.
    """
    filetype = 'ctrlvl'
    descriptions = [
        'Initial level of the ion before recombination.',
        'Final level of the ion after recombination.',
        'Atomic number of the perturber involved in the transition.',
        'Number of electrons in perturber before transition.',
        'Input energy required for the transition to take place.',
        'Recombination rate coefficients for the transition.',
        'Temperatures at which the rate coefficients are evaluated.'
    ]
