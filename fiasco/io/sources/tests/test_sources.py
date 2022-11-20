"""
Test parsers for all files attached to specific ions
"""
import pytest

from astropy.table import QTable

import fiasco.io


@pytest.mark.parametrize('filename', [
    'h_1.elvlc',
    'h_1.fblvl',
    'h_1.scups',
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
    'sun_coronal_1992_feldman.abund',
    'chianti.ip',
    'gffgu.dat',
    'gffint.dat',
    'hseq_2photon.dat',
    'itoh.dat',
    'klgfb.dat',
    'verner_short.txt',
])
def test_non_ion_sources(ascii_dbase_root, filename):
    parser = fiasco.io.Parser(filename, ascii_dbase_root=ascii_dbase_root)
    table = parser.parse()
    assert isinstance(table, QTable)
    assert all([k in table.meta for k in ('footer', 'filename', 'chianti_version', 'descriptions')])
    assert table.colnames == parser.headings
    for h, unit in zip(parser.headings, parser.units):
        assert table[h].unit == unit
    assert table.meta['filename'] == filename
