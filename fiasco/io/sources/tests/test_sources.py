"""
Test parsers for all files attached to specific ions
"""
import numpy as np
import pathlib
import plasmapy.particles
import pytest

from astropy.table import QTable

import fiasco.io


@pytest.mark.parametrize('filename', [
    'h_1.elvlc',
    'h_1.fblvl',
    pytest.param('h_1.scups', marks=pytest.mark.requires_dbase_version('>= 8')),
    'c_2.psplups',
    'be_2.easplom',
    'al_3.easplups',
    'o_2.wgfa',
    'fe_17.cilvl',
    'fe_17.reclvl',
    'al_3.rrparams',
    'fe_2.trparams',
    'fe_12.drparams',
    'al_3.diparams',
    'al_6.diparams',  # NOTE: This file does not have the row of EA scaling factors
    pytest.param('fe_23.auto', marks=pytest.mark.requires_dbase_version('>= 9')),
    pytest.param('fe_23.rrlvl', marks=pytest.mark.requires_dbase_version('>= 9')),
    pytest.param('c_5.splups', marks=pytest.mark.requires_dbase_version('< 8')),
    pytest.param('c_6.splups', marks=pytest.mark.requires_dbase_version('< 8')),
    pytest.param('c_3.dilvl', marks=pytest.mark.requires_dbase_version('>= 11')),
    pytest.param('c_3.ealvl', marks=pytest.mark.requires_dbase_version('>= 11')),
    pytest.param('c_3.rrcoeffs', marks=pytest.mark.requires_dbase_version('>= 11')),
    pytest.param('c_3.drcoeffs', marks=pytest.mark.requires_dbase_version('>= 11')),
    pytest.param('c_2.ctilvl', marks=pytest.mark.requires_dbase_version('>= 11')),
    pytest.param('c_2.ctrlvl', marks=pytest.mark.requires_dbase_version('>= 11')),
])
def test_ion_sources(ascii_dbase_root, filename,):
    parser = fiasco.io.Parser(filename, ascii_dbase_root=ascii_dbase_root)
    table = parser.parse()
    assert isinstance(table, QTable)
    assert all([k in table.meta for k in
                ('footer', 'filename', 'element', 'ion', 'chianti_version', 'descriptions')])
    assert table.colnames == parser.headings
    for h, unit in zip(parser.headings, parser.units):
        assert table[h].unit == unit
    assert table.meta['element'] == filename.split('_')[0]
    assert table.meta['ion'] == filename.split('.')[0]
    assert table.meta['filename'] == filename


@pytest.mark.parametrize('filename', [
    'chianti.ioneq',
    pytest.param('sun_coronal_1992_feldman.abund', marks=pytest.mark.requires_dbase_version('< 10')),
    pytest.param('sun_coronal_2021_chianti.abund', marks=pytest.mark.requires_dbase_version('>= 10')),
    'chianti.ip',
    'gffgu.dat',
    'gffint.dat',
    'hseq_2photon.dat',
    'itoh.dat',
    pytest.param('klgfb.dat', marks=pytest.mark.requires_dbase_version('< 10')),
    pytest.param('klgfb_6.dat', marks=pytest.mark.requires_dbase_version('>= 10')),
    'verner_short.txt',
    'flare.dem',
    pytest.param('itoh_integrated_gaunt.txt', marks=pytest.mark.requires_dbase_version('>= 9.0.1')),
    pytest.param('itoh_integrated_gaunt_nonrel.txt', marks=pytest.mark.requires_dbase_version('>= 9.0.1')),
    pytest.param('advmodel_list.ions', marks=pytest.mark.requires_dbase_version('>= 11')),
    pytest.param(pathlib.Path('model_atmospheres') / 'fontenla_plage.dat', marks=pytest.mark.requires_dbase_version('>= 11'))
])
def test_non_ion_sources(ascii_dbase_root, filename):
    parser = fiasco.io.Parser(filename, ascii_dbase_root=ascii_dbase_root)
    table = parser.parse()
    assert isinstance(table, QTable)
    assert all([k in table.meta for k in ('footer', 'filename', 'chianti_version', 'descriptions')])
    assert all([h in table.colnames for h in parser.headings])
    for h, unit in zip(parser.headings, parser.units):
        assert table[h].unit == unit
    assert table.meta['filename'] == filename


@pytest.mark.parametrize('filename', [
    # Many abundance files were moved in v11 of the database
    pytest.param('sun_coronal_1992_feldman.abund', marks=pytest.mark.requires_dbase_version('< 10')),
    pytest.param('archive/sun_coronal_1992_feldman.abund', marks=pytest.mark.requires_dbase_version('>= 10')),
    pytest.param('sun_coronal_2021_chianti.abund', marks=pytest.mark.requires_dbase_version('>= 10')),
    'version_3/allen.abund',
    'version_3/grevesse_anders.abund'
])
def test_abundance_parsing(ascii_dbase_root, filename):
    parser = fiasco.io.Parser(filename, ascii_dbase_root=ascii_dbase_root)
    table = parser.parse()
    assert all([h in table.colnames for h in ['Z', 'abundance', 'element']])
    for row in table:
        Z = row['Z'].item() if isinstance(row['Z'], np.int64) else row['Z']
        assert row['element'] == plasmapy.particles.atomic_symbol(Z)
        assert row['element'] == plasmapy.particles.atomic_symbol(row['element'])
