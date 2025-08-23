"""
Source classes for CHIANTI filetypes not attached to ions
"""
import astropy.units as u
import fortranformat
import numpy as np
import pathlib
import plasmapy

from astropy.table import Column

from fiasco.io.generic import GenericParser

__all__ = [
    'AbundParser',
    'IoneqParser',
    'IpParser',
    'DemParser',
    'AdvancedModelListParser',
    'ModelAtmosphereParser',
]


class AbundParser(GenericParser):
    filetype = 'abund'
    dtypes = [int, float, str]
    units = [None, u.dimensionless_unscaled, None]
    headings = ['Z', 'abundance', 'element']
    descriptions = ['atomic number', 'abundance relative to H', 'element']
    fformat = fortranformat.FortranRecordReader('(I2,F9.3,A2)')

    def __init__(self, abundance_filename, **kwargs):
        super().__init__(abundance_filename, **kwargs)
        self.full_path = pathlib.Path(self.ascii_dbase_root / 'abundance' / self.filename)

    def postprocessor(self, df):
        df['abundance'] = 10.**(df['abundance'] - df['abundance'][df['Z'] == 1])
        # repair missing data
        if df['element'][0] == '':
            col = []
            for atomic_number in df['Z']:
                col.append(plasmapy.particles.atomic_symbol(int(atomic_number)))
            df['element'] = Column(col)
        df = super().postprocessor(df)
        return df

    def to_hdf5(self, hf, df):
        dataset_name = pathlib.Path(self.filename).stem
        footer = f"""{dataset_name}
------------------
{df.meta['footer']}"""
        for row in df:
            grp_name = '/'.join([row['element'].lower(), 'abundance'])
            if grp_name not in hf:
                grp = hf.create_group(grp_name)
                grp.attrs['footer'] = ''
                grp.attrs['chianti_version'] = df.meta['chianti_version']
            else:
                grp = hf[grp_name]
            grp.attrs['footer'] += footer
            if dataset_name not in grp:
                ds = grp.create_dataset(dataset_name, data=row['abundance'])
                ds.attrs['unit'] = df['abundance'].unit.to_string()
                ds.attrs['description'] = df.meta['descriptions']['abundance']


class IoneqParser(GenericParser):
    filetype = 'ioneq'
    dtypes = [int, int, float, float]
    units = [None, None, u.K, u.dimensionless_unscaled]
    headings = ['Z', 'ion', 'temperature', 'ionization_fraction']
    descriptions = ['atomic number', 'ion', 'temperature', 'ionization fraction']

    def __init__(self, ioneq_filename, **kwargs):
        super().__init__(ioneq_filename, **kwargs)
        self.full_path = pathlib.Path(self.ascii_dbase_root / 'ioneq' / self.filename)

    def preprocessor(self, table, line, index):
        if index == 0:
            num_entries = int(line.strip().split()[0])
            self.fformat_temperature = fortranformat.FortranRecordReader(f'{num_entries}F6.2')
            self.fformat_ioneq = fortranformat.FortranRecordReader(f'2I3,{num_entries}E10.2')
        elif index == 1:
            self.temperature = 10.**np.array(self.fformat_temperature.read(line), dtype=float)
        else:
            line = self.fformat_ioneq.read(line)
            line = line[:2] + [self.temperature, np.array(line[2:], dtype=float)]
            table.append(line)

    def to_hdf5(self, hf, df):
        dataset_name = pathlib.Path(self.filename).stem
        for row in df:
            el = plasmapy.particles.atomic_symbol(int(row['Z'])).lower()
            ion = int(row['ion'])
            grp_name = '/'.join([el, f'{el}_{ion}', 'ioneq'])
            if grp_name not in hf:
                grp = hf.create_group(grp_name)
            else:
                grp = hf[grp_name]
            if dataset_name not in grp:
                sub_grp = grp.create_group(dataset_name)
                sub_grp.attrs['footer'] = df.meta['footer']
                sub_grp.attrs['chianti_version'] = df.meta['chianti_version']
                ds = sub_grp.create_dataset('temperature', data=row['temperature'])
                ds.attrs['unit'] = df['temperature'].unit.to_string()
                ds.attrs['description'] = df.meta['descriptions']['temperature']
                ds = sub_grp.create_dataset('ionization_fraction', data=row['ionization_fraction'])
                ds.attrs['unit'] = df['ionization_fraction'].unit.to_string()
                ds.attrs['description'] = df.meta['descriptions']['ionization_fraction']


class IpParser(GenericParser):
    filetype = 'ip'
    dtypes = [int, int, float]
    units = [None, None, 1/u.cm]
    headings = ['Z', 'ion', 'ip']
    descriptions = ['atomic number', 'ion', 'ionization potential']

    def __init__(self, ip_filename, **kwargs):
        super().__init__(ip_filename, **kwargs)
        self.full_path = pathlib.Path(self.ascii_dbase_root / 'ip' / self.filename)

    def to_hdf5(self, hf, df):
        dataset_name = pathlib.Path(self.filename).stem
        footer = f"""{dataset_name}
------------------
{df.meta['footer']}"""
        for row in df:
            el = plasmapy.particles.atomic_symbol(int(row['Z'])).lower()
            ion = int(row['ion'])
            grp_name = '/'.join([el, f'{el}_{ion}', 'ip'])
            if grp_name not in hf:
                grp = hf.create_group(grp_name)
                grp.attrs['footer'] = ''
                grp.attrs['chianti_version'] = df.meta['chianti_version']
            else:
                grp = hf[grp_name]
            grp.attrs['footer'] += footer
            if dataset_name not in grp:
                ds = grp.create_dataset(dataset_name, data=row['ip'])
                ds.attrs['unit'] = df['ip'].unit.to_string()
                ds.attrs['description'] = df.meta['descriptions']['ip']


class DemParser(GenericParser):
    filetype = 'dem'
    dtypes = [float, float]
    units = [u.K, u.cm**(-5)/u.K]
    headings = ['temperature_bin_center', 'dem']
    descriptions = ['center of the temperature bin', 'differential emission measure']

    def __init__(self, dem_filename, **kwargs):
        super().__init__(dem_filename, **kwargs)
        self.full_path = pathlib.Path(self.ascii_dbase_root / 'dem' / self.filename)

    def postprocessor(self, df):
        logT_center = df['temperature_bin_center'].value
        T_unit = df['temperature_bin_center'].unit
        df['temperature_bin_center'] = 10**logT_center * T_unit
        df['dem'] = 10**df['dem'].value * df['dem'].unit
        df = super().postprocessor(df)
        # NOTE: there are several DEM models which clearly do not have equal bin widths
        # so this conversion to EM and calculation of the edges is not valid
        # These are the 'flare_ext' and 'AU_mic' models.
        # For these models, the below steps are not valid so they are skipped.
        bins_equal = np.allclose(np.diff(logT_center), np.diff(logT_center)[0], atol=1e-15, rtol=0)
        if bins_equal:
            delta_logT = np.diff(logT_center)[0]
            logT_left = logT_center - delta_logT/2
            logT_right = logT_center + delta_logT/2
            df['temperature_bin_edge_left'] = 10**logT_left*T_unit
            df['temperature_bin_edge_right'] = 10**logT_right*T_unit
            df['temperature_bin_width'] = df['temperature_bin_edge_right'] - df['temperature_bin_edge_left']
            df['em'] = df['dem'] * df['temperature_bin_width']
            df.meta['descriptions'].update({
                'temperature_bin_edge_left': 'Left edge of the temperature bin',
                'temperature_bin_edge_right': 'Right edge of the temperature bin',
                'temperature_bin_width': 'Width of the temperature bin. Should be uniform in log10 space',
                'em': 'Emission measure; DEM integrated over each temperature bin',
            })
        else:
            from fiasco import log
            log.debug(f'''Cannot compute additional quantities for {self.filename}.
                          Temperature bins are not of equal width in log10 space.''')
        return df

    def to_hdf5(self, hf, df):
        dataset_name = pathlib.Path(self.filename).stem
        # NOTE: the following is necessary as there are sometimes subdirectories within the DEM
        # directory and some files within these subdirectories may have names that conflict with
        # DEM files in the top level DEM directory. Thus, we append the name of the subdirectory
        # to the dataset name to differentiate these DEM datasets.
        # NOTE: For relative paths, the first entry in 'parents' is always ".", the current
        # directory so we skip this.
        suffixes = [str(p) for p in pathlib.Path(self.filename).parents][:-1]
        dataset_name = '_'.join([dataset_name] + suffixes)
        footer = f"""{dataset_name}
------------------
{df.meta['footer']}"""
        grp_name = '/'.join(['dem', dataset_name])
        if grp_name not in hf:
            grp = hf.create_group(grp_name)
            grp.attrs['footer'] = footer
            grp.attrs['chianti_version'] = df.meta['chianti_version']
        else:
            grp = hf[grp_name]
        for col in df.colnames:
            ds = grp.create_dataset(col, data=df[col].value)
            ds.attrs['description'] = df.meta['descriptions'][col]
            ds.attrs['unit'] = df[col].unit.to_string()


class AdvancedModelListParser(GenericParser):
    filetype = 'advmodel_list'
    dtypes = [str, int]
    units = [None, None]
    headings = ['ion', 'n_levels']
    descriptions = ['ion name', 'number of included levels']
    fformat = fortranformat.FortranRecordReader('(A6,I5)')

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.full_path = pathlib.Path(self.ascii_dbase_root / 'ancillary_data' / 'advanced_models' / self.filename)

    def to_hdf5(self, hf, df):
        # The number of levels should be added as an attribute to the ion group
        # and then a property added to IonBase to access these number of advanced models.
        # This should default to 0.
        for row in df:
            ion_name = row['ion']
            element_name = ion_name.split('_')[0]
            grp_name = '/'.join([element_name, ion_name, 'advanced_model'])
            if grp_name not in hf:
                grp = hf.create_group(grp_name)
                grp.attrs['footer'] = df.meta['footer']
                grp.attrs['chianti_version'] = df.meta['chianti_version']
            else:
                grp = hf[grp_name]
            if 'n_levels' not in grp:
                ds = grp.create_dataset('n_levels', data=row['n_levels'])
                ds.attrs['description'] = df.meta['descriptions']['n_levels']
                ds.attrs['unit'] = 'SKIP'


class ModelAtmosphereParser(GenericParser):
    filetype = 'model_atmospheres'
    dtypes = 8*[float]
    units = [
        u.K,
        u.cm**(-3),
        u.km,
        u.K*u.cm**(-3),
        u.cm**(-3),
        u.dimensionless_unscaled,
        u.dimensionless_unscaled,
        u.dimensionless_unscaled,
    ]
    headings = [
        'temperature',
        'density_e',
        'height',
        'pressure',
        'density_H',
        'fraction_H_1',
        'fraction_He_1',
        'fraction_He_2',
    ]
    descriptions = [
        'temperature',
        'electron number density',
        'height',
        'pressure',
        'total hydrogen number density',
        'ionization fraction of neutral hydrogen',
        'ionization fraction of neutral helium',
        'ionization fraction of singly-ionized helium',
    ]
    fformat = fortranformat.FortranRecordReader('(E9.3,7E12.3)')

    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.full_path = pathlib.Path(
            self.ascii_dbase_root / 'ancillary_data' / 'advanced_models' / 'model_atmospheres' / pathlib.Path(self.filename).name
        )

    def to_hdf5(self, hf, df):
        dataset_name = pathlib.Path(self.filename).stem
        footer = f"""{dataset_name}
------------------
{df.meta['footer']}"""
        grp_name = '/'.join(['model_atmospheres', dataset_name])
        if grp_name not in hf:
            grp = hf.create_group(grp_name)
            grp.attrs['footer'] = footer
            grp.attrs['chianti_version'] = df.meta['chianti_version']
        else:
            grp = hf[grp_name]
        for col in df.colnames:
            ds = grp.create_dataset(col, data=df[col].value)
            ds.attrs['description'] = df.meta['descriptions'][col]
            ds.attrs['unit'] = df[col].unit.to_string()
