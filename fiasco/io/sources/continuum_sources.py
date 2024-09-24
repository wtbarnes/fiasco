"""
Source classes for continuum files
"""
import astropy.units as u
import h5py
import numpy as np
import pathlib
import plasmapy

from fiasco.io.generic import GenericParser

__all__ = [
    'GffguParser',
    'GffintParser',
    'ItohIntegratedGauntParser',
    'ItohIntegratedGauntNonrelParser',
    'KlgfbParser',
    'VernerParser',
    'ItohParser',
    'HSeqParser',
    'HeSeqParser'
]


class GffguParser(GenericParser):
    """
    Free-free Gaunt factor as a function of scaled frequency and energy
    """
    filetype = 'gffgu'
    dtypes = [float, float, float]
    units = [u.dimensionless_unscaled, u.dimensionless_unscaled, u.dimensionless_unscaled]
    headings = ['u', 'gamma_squared', 'gaunt_factor']
    descriptions = ['scaled frequency', 'scaled temperature', 'free-free Gaunt factor']

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.full_path = pathlib.Path(kwargs.get('full_path',
                                    self.ascii_dbase_root / 'continuum' / filename))
        self.body_index = 5

    def preprocessor(self, table, line, index):
        if index >= self.body_index and '--' not in line:
            super().preprocessor(table, line, index)

    def extract_footer(self, lines):
        comment = []
        for i, l in enumerate(lines):
            if i <= self.body_index - 3:
                comment.append(l)
            else:
                break

        footer = '\n'.join([l.strip() for l in comment])
        return footer

    def to_hdf5(self, hf, df, **kwargs):
        grp_name = '/'.join(['continuum', self.filetype])
        if grp_name not in hf:
            grp = hf.create_group(grp_name)
            grp.attrs['chianti_version'] = df.meta['chianti_version']
            grp.attrs['footer'] = df.meta['footer']
        else:
            grp = hf[grp_name]

        for name in df.colnames:
            col = df[name]
            if type(col) == u.Quantity:
                data = col.value
            else:
                data = col.data
            if '<U' in data.dtype.str:
                numchar = data.dtype.str[2:]
                data = data.astype(f'|S{numchar}')
            if name in grp:
                ds = grp[name]
            else:
                if data.dtype == np.dtype('O'):
                    ragged_dtype = h5py.special_dtype(vlen=np.dtype('float64'))
                    ds = grp.create_dataset(name, data=data, dtype=ragged_dtype)
                else:
                    ds = grp.create_dataset(name, data=data, dtype=data.dtype)
            if col.unit is None:
                ds.attrs['unit'] = 'SKIP'
            else:
                ds.attrs['unit'] = col.unit.to_string()
            ds.attrs['description'] = df.meta['descriptions'][name]


class GffintParser(GffguParser):
    """
    Total free-free Gaunt factor as a function of scaled temperature.
    """
    filetype = 'gffint'
    dtypes = [float, float, float, float, float]
    units = [u.dimensionless_unscaled, u.dimensionless_unscaled, u.dimensionless_unscaled,
             u.dimensionless_unscaled, u.dimensionless_unscaled]
    headings = ['log_gamma_squared', 'gaunt_factor', 's1', 's2', 's3']
    descriptions = ['log scaled temperature', 'total free-free Gaunt factor', 'spline coefficient',
                    'spline coefficient', 'spline coefficient']

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.body_index = 4

class ItohIntegratedGauntParser(GenericParser):
    """
    Total (frequency-integrated) relativistic free-free Gaunt factor as a
    function of a scaled temperature and scaled atomic number.
    """
    filetype='itoh_integrated_gaunt'
    dtypes=[float]
    units = [u.dimensionless_unscaled]
    headings = ['a_ik']
    descriptions = ['fitting coefficient']

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.full_path = pathlib.Path(kwargs.get('full_path',
                                    self.ascii_dbase_root / 'continuum' / filename))

    def preprocessor(self, table, line, index):
        line = line.strip().split()
        gf = np.array(line, dtype=float)
        table.append([gf])

    def extract_footer(self, *args):
        return """Fit parameters for relativistic free-free Gaunt factors integrated over frequency.
Itoh et al., 2002, A&A, 382, 722
comment: These are the coefficients a_ik tabulated in Table 1."""

    def to_hdf5(self, hf, df, **kwargs):
        grp_name = '/'.join(['continuum', self.filetype])
        if grp_name not in hf:
            grp = hf.create_group(grp_name)
            grp.attrs['chianti_version'] = df.meta['chianti_version']
            grp.attrs['footer'] = df.meta['footer']
        else:
            grp = hf[grp_name]

        for name in df.colnames:
            col = df[name]
            if type(col) == u.Quantity:
                data = col.value
            else:
                data = col.data
            if '<U' in data.dtype.str:
                numchar = data.dtype.str[2:]
                data = data.astype(f'|S{numchar}')
            if name in grp:
                ds = grp[name]
            else:
                if data.dtype == np.dtype('O'):
                    ragged_dtype = h5py.special_dtype(vlen=np.dtype('float64'))
                    ds = grp.create_dataset(name, data=data, dtype=ragged_dtype)
                else:
                    ds = grp.create_dataset(name, data=data, dtype=data.dtype)
            ds.attrs['unit'] = col.unit.to_string()
            ds.attrs['description'] = df.meta['descriptions'][name]


class ItohIntegratedGauntNonrelParser(GenericParser):
    """
    Total (frequency-integrated) non-relativistic free-free Gaunt factor as a
    function of a scaled temperature.
    """
    filetype='itoh_integrated_gaunt_nonrel'
    dtypes=[float]
    units = [u.dimensionless_unscaled]
    headings = ['b_i']
    descriptions = ['fitting coefficient']

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.full_path = pathlib.Path(kwargs.get('full_path',
                                    self.ascii_dbase_root / 'continuum' / filename))

    def extract_footer(self, *args):
        return """Fit coefficients for non-relativistic, frequency integrated free-free Gaunt factor
Itoh et al. (2002, A&A, 382, 722)
Comment: Data taken from Table 2 of this work."""

    def to_hdf5(self, hf, df, **kwargs):
        grp_name = '/'.join(['continuum', self.filetype])
        if grp_name not in hf:
            grp = hf.create_group(grp_name)
            grp.attrs['chianti_version'] = df.meta['chianti_version']
            grp.attrs['footer'] = df.meta['footer']
        else:
            grp = hf[grp_name]

        for name in df.colnames:
            col = df[name]
            if type(col) == u.Quantity:
                data = col.value
            else:
                data = col.data
            if '<U' in data.dtype.str:
                numchar = data.dtype.str[2:]
                data = data.astype(f'|S{numchar}')
            if name in grp:
                ds = grp[name]
            else:
                if data.dtype == np.dtype('O'):
                    ragged_dtype = h5py.special_dtype(vlen=np.dtype('float64'))
                    ds = grp.create_dataset(name, data=data, dtype=ragged_dtype)
                else:
                    ds = grp.create_dataset(name, data=data, dtype=data.dtype)
            ds.attrs['unit'] = col.unit.to_string()
            ds.attrs['description'] = df.meta['descriptions'][name]


class KlgfbParser(GenericParser):
    """
    Free-bound gaunt factor as a function of photon energy for several different energy levels.
    """
    filetype = 'klgfb'
    dtypes = [int, int, float, float]
    units = [None, None, u.dimensionless_unscaled, u.dimensionless_unscaled]
    headings = ['n', 'l', 'log_pe', 'log_gaunt_factor']
    descriptions = ['principal quantum number', 'orbital angular momentum number',
                    'log photon energy divided by ionization potential',
                    'log free-bound Gaunt factor']

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.full_path = pathlib.Path(kwargs.get('full_path',
                                    self.ascii_dbase_root / 'continuum' / filename))

    def preprocessor(self, table, line, index):
        if index == 0:
            pass
        elif index == 1:
            self._photon_energy = np.array(line.strip().split(), dtype=float)
        else:
            line = line.strip().split()
            gf = np.array(line[2:], dtype=float)
            table.append(line[:2] + [self._photon_energy, gf])

    def extract_footer(self, *args):
        return """Log of free-bound Gaunt factors as a function of log of photon energy divided by ionization potential
From Karzas, W. J. and Latter, R., 1961, ApJS, 6, 167"""

    def to_hdf5(self, hf, df, **kwargs):
        grp_name = '/'.join(['continuum', self.filetype])
        if grp_name not in hf:
            grp = hf.create_group(grp_name)
            grp.attrs['chianti_version'] = df.meta['chianti_version']
            grp.attrs['footer'] = df.meta['footer']
        else:
            grp = hf[grp_name]

        for name in df.colnames:
            col = df[name]
            if type(col) == u.Quantity:
                data = col.value
            else:
                data = col.data
            if '<U' in data.dtype.str:
                numchar = data.dtype.str[2:]
                data = data.astype(f'|S{numchar}')
            if name in grp:
                ds = grp[name]
            else:
                if data.dtype == np.dtype('O'):
                    ragged_dtype = h5py.special_dtype(vlen=np.dtype('float64'))
                    ds = grp.create_dataset(name, data=data, dtype=ragged_dtype)
                else:
                    ds = grp.create_dataset(name, data=data, dtype=data.dtype)
            if col.unit is None:
                ds.attrs['unit'] = 'SKIP'
            else:
                ds.attrs['unit'] = col.unit.to_string()
            ds.attrs['description'] = df.meta['descriptions'][name]


class VernerParser(GenericParser):
    """
    Fit parameters for calculating partial photoionization cross-sections using
    the method of :cite:t:`verner_analytic_1995`.
    """
    filetype = 'verner_short'
    dtypes = [int, int, int, int, float, float, float, float, float, float]
    units = [None, None, None, None, u.eV, u.eV, u.megabarn, u.dimensionless_unscaled,
             u.dimensionless_unscaled, u.dimensionless_unscaled]
    headings = ['Z', 'n_electrons', 'n', 'l', 'E_thresh', 'E_0_fit', 'sigma_0', 'y_a_fit', 'P_fit',
                'y_w_fit']
    descriptions = ['atomic number', 'number of electrons', 'principal quantum number',
                    'orbital angular momentum number',
                    'threshold energy below which cross-section is 0',
                    'E_0 fit parameter', 'nominal value of cross-section', 'y_a fit parameter',
                    'P fit parameter', 'y_w fit parameter']

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.full_path = pathlib.Path(kwargs.get('full_path',
                                    self.ascii_dbase_root / 'continuum' / filename))

    def extract_footer(self, *args):
        return """Fit parameters for calculating partial photoionization cross-sections for individual ions
From Verner, D. A . and Yakovlev, D. G., 1995, A&AS, 109, 125"""

    def to_hdf5(self, hf, df, **kwargs):
        for row in df:
            el = plasmapy.particles.atomic_symbol(int(row['Z'])).lower()
            stage = row['Z'] - row['n_electrons'] + 1
            grp_name = f'{el}/{el}_{stage}/continuum/{self.filetype}'
            if grp_name not in hf:
                grp = hf.create_group(grp_name)
                grp.attrs['chianti_version'] = df.meta['chianti_version']
                grp.attrs['footer'] = df.meta['footer']
            else:
                grp = hf[grp_name]
            for col in row.colnames:
                if col == 'Z' or col == 'n_electrons':
                    continue
                ds = grp.create_dataset(col, data=row[col])
                ds.attrs['description'] = df.meta['descriptions'][col]
                if not hasattr(row[col], 'unit'):
                    ds.attrs['unit'] = 'SKIP'
                else:
                    ds.attrs['unit'] = row[col].unit.to_string()


class ItohParser(GenericParser):
    """
    Fit parameters for calculating relativistic free-free Gaunt factor using the method of
    :cite:t:`itoh_relativistic_2000`.
    """
    filetype = 'itoh'
    dtypes = [int, float]
    units = [None, u.dimensionless_unscaled]
    headings = ['Z', 'a']
    descriptions = ['atomic number', 'fit coefficient']

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.full_path = pathlib.Path(kwargs.get('full_path',
                                    self.ascii_dbase_root / 'continuum' / filename))

    def preprocessor(self, table, line, index):
        a_matrix = np.array(line.strip().split()).reshape((11, 11))
        table.append([index+1] + [a_matrix])

    def extract_footer(self, *args):
        return """Analytic fit coefficients as a function of scaled temperature and energy for calculating the relativistic free-free Gaunt factor
From Itoh, N., et al., ApJS, 2000, 128, 125"""

    def to_hdf5(self, hf, df, **kwargs):
        grp_name = '/'.join(['continuum', self.filetype])
        if grp_name not in hf:
            grp = hf.create_group(grp_name)
            grp.attrs['chianti_version'] = df.meta['chianti_version']
            grp.attrs['footer'] = df.meta['footer']
        else:
            grp = hf[grp_name]

        for name in df.colnames:
            col = df[name]
            if type(col) == u.Quantity:
                data = col.value
            else:
                data = col.data
            if '<U' in data.dtype.str:
                numchar = data.dtype.str[2:]
                data = data.astype(f'|S{numchar}')
            if name in grp:
                ds = grp[name]
            else:
                if data.dtype == np.dtype('O'):
                    ragged_dtype = h5py.special_dtype(vlen=np.dtype('float64'))
                    ds = grp.create_dataset(name, data=data, dtype=ragged_dtype)
                else:
                    ds = grp.create_dataset(name, data=data, dtype=data.dtype)
            if col.unit is None:
                ds.attrs['unit'] = 'SKIP'
            else:
                ds.attrs['unit'] = col.unit.to_string()
            ds.attrs['description'] = df.meta['descriptions'][name]


class HSeqParser(GenericParser):
    r"""
    Parameters for calculating two-photon continuum for hydrogen-like ions

    Notes
    -----
    * The parameter :math: `\psi_{\text{norm}}` (called :math: `A_{\text{sum}}` in CHIANTI) is a
      normalization factor of the integral of the spectral distribution function :math:`\psi(y)` from 0 to 1,
      such that :math: `\frac{1}{\psi_{\text{norm}}} \int_{0}^{1} \psi(y) dy = 2`.
      This normalization is only used for hydrogenic ions.
    """
    filetype = 'hseq_2photon'
    dtypes = [int, float, int, float, float, float]
    units = [None, u.dimensionless_unscaled, None, 1/u.s, u.dimensionless_unscaled, u.dimensionless_unscaled]
    headings = ['Z', 'y', 'Z_0', 'A', 'psi_norm', 'psi']
    descriptions = ['atomic number', 'fraction of energy carried by one of the two photons',
                    'nominal atomic number', 'radiative decay rate', 'normalization of the integral of psi from 0 to 1',
                    'spectral distribution function']

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.full_path = pathlib.Path(kwargs.get('full_path',
                                    self.ascii_dbase_root / 'continuum' / filename))

    def preprocessor(self, table, line, index):
        if index == 0:
            self._y0 = np.array(line.strip().split(), dtype=float)
        elif index == 1:
            self._z0 = np.array(line.strip().split(), dtype=float)
        else:
            line = line.strip().split()
            table.append([line[0]] + [self._y0] + [self._z0] + line[1:3] + [np.array(line[3:])])

    def extract_footer(self, *args):
        return """Information needed for calculating two-photon continuum emission for hydrogen isoelectronic sequence
Radiative decay rates from Parpia, F. A., and Johnson, W. R., 1982, Phys. Rev. A, 26, 1142
Spectral distribution function from Goldman, S.P. and Drake, G.W.F., 1981, Phys Rev A, 24, 183"""

    def to_hdf5(self, hf, df, **kwargs):
        for row in df:
            el = plasmapy.particles.atomic_symbol(int(row['Z'])).lower()
            grp_name = f'{el}/continuum/{self.filetype}'
            if grp_name not in hf:
                grp = hf.create_group(grp_name)
                grp.attrs['chianti_version'] = df.meta['chianti_version']
                grp.attrs['footer'] = df.meta['footer']
            else:
                grp = hf[grp_name]
            for col in row.colnames:
                if col == 'Z':
                    continue
                ds = grp.create_dataset(col, data=row[col])
                ds.attrs['description'] = df.meta['descriptions'][col]
                if not hasattr(row[col], 'unit'):
                    ds.attrs['unit'] = 'SKIP'
                else:
                    ds.attrs['unit'] = row[col].unit.to_string()


class HeSeqParser(GenericParser):
    """
    Parameters for calculating two-photon continuum for helium-like ions.
    """
    filetype = 'heseq_2photon'
    dtypes = [int, float, float, float]
    units = [None, u.dimensionless_unscaled, 1/u.s, u.dimensionless_unscaled]
    headings = ['Z', 'y', 'A', 'psi']
    descriptions = ['atomic number', 'fraction of energy carried by one of the two photons',
                    'radiative decay rate', 'spectral distribution function']

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.full_path = pathlib.Path(kwargs.get('full_path',
                                    self.ascii_dbase_root / 'continuum' / filename))

    def preprocessor(self, table, line, index):
        if index == 0:
            self._y0 = np.array(line.strip().split(), dtype=float)
        else:
            line = line.strip().split()
            table.append([line[0]] + [self._y0] + line[1:2] + [np.array(line[2:])])

    def extract_footer(self, *args):
        return """Information needed for calculating two-photon continuum emission for helium isoelectronic sequence
Radiative decay rates from Drake, G.W.F., 1986, Phys. Rev. A, 34, 2871
Spectral distribution function from Drake, G.W.F., Victor, G.A., Dalgarno, A., 1969, Phys. Rev. A, 180, 25."""

    def to_hdf5(self, hf, df, **kwargs):
        for row in df:
            el = plasmapy.particles.atomic_symbol(int(row['Z'])).lower()
            grp_name = f'{el}/continuum/{self.filetype}'
            if grp_name not in hf:
                grp = hf.create_group(grp_name)
                grp.attrs['chianti_version'] = df.meta['chianti_version']
                grp.attrs['footer'] = df.meta['footer']
            else:
                grp = hf[grp_name]
            for col in row.colnames:
                if col == 'Z':
                    continue
                ds = grp.create_dataset(col, data=row[col])
                ds.attrs['description'] = df.meta['descriptions'][col]
                if not hasattr(row[col], 'unit'):
                    ds.attrs['unit'] = 'SKIP'
                else:
                    ds.attrs['unit'] = row[col].unit.to_string()
