"""
Test parsers for all files attached to specific ions
"""
import os
import pytest
import numpy as np
from astropy.table import QTable

from fiasco.io.sources.ion_sources import (ElvlcParser, FblvlParser, ScupsParser, PsplupsParser, 
                                           EasplomParser)


def test_elvlcparser(tmpdir):
    """
    Test parser for energy level files
    """
    # Write temporary file
    f = tmpdir.mkdir('ion_sources').join('h_1.elvlc')
    f.write("""      1     1s                                2  S    0.5          0.000          0.000
-1
%filename: h_1.elvlc
%observed energy levels: Fuhr et al, 1999, NIST Atomic Spectra Database Version 2.0
%produced as part of the Arcetri/Cambridge/NRL 'CHIANTI' atomic data base collaboration
%
%  Ken Dere  May 3 2001""")
    # Make model dataframe
    table = QTable(data=[[1], ['1s'], [''], [2], ['S'], [0.5], [0.000], [0.000]],
                   names=ElvlcParser.headings,
                   meta={'element': 'h',
                         'ion': 'h_1',
                         'filename': os.path.join(f.dirname, f.basename),
                         'footer': "filename: h_1.elvlc\nobserved energy levels: Fuhr et al, 1999, NIST Atomic Spectra Database Version 2.0\nproduced as part of the Arcetri/Cambridge/NRL 'CHIANTI' atomic data base collaboration\nKen Dere  May 3 2001"})
    for h, unit in zip(ElvlcParser.headings, ElvlcParser.units):
        table[h].unit = unit
    # Read file
    p = ElvlcParser(os.path.join(f.dirname, f.basename), standalone=True)
    table_parsed = p.parse()
    # Check columns
    for c in table.colnames:
        assert all(table[c] == table_parsed[c])
    # Check metadata
    for k in ('footer', 'filename', 'element', 'ion'):
        assert table_parsed.meta[k] == table.meta[k]
    assert 'chianti_version' in table_parsed.meta
    assert 'descriptions' in table_parsed.meta


def test_fblvlparser(tmpdir):
    """
    Test parser for free-bound level files
    """
    f = tmpdir.mkdir('ion_sources').join('h_1.fblvl')
    f.write("""    1                  1s    1    0  s    2               0.000               0.000
 -1
%filename: h_1.fblvl
%observed energy levels: Fuhr et al, 1999, NIST Atomic Spectra Database Version 2.0
%additional energy levels: NIST Atomic Spectra Database Version 2.0
%comment: free-bound data file developed from .elvlc file
%produced as part of the Arcetri/Cambridge/NRL 'CHIANTI' atomic data base collaboration
%
% Ken Dere - Jul 17 2002
 -1""")
    table = QTable(data=[[1], ['1s'], [1], [0], ['s'], [2], [0.000], [0.000]],
                   names=FblvlParser.headings,
                   meta={'element': 'h',
                         'ion': 'h_1',
                         'filename': os.path.join(f.dirname, f.basename),
                         'footer': "filename: h_1.fblvl\nobserved energy levels: Fuhr et al, 1999, NIST Atomic Spectra Database Version 2.0\nadditional energy levels: NIST Atomic Spectra Database Version 2.0\ncomment: free-bound data file developed from .elvlc file\nproduced as part of the Arcetri/Cambridge/NRL 'CHIANTI' atomic data base collaboration\nKen Dere - Jul 17 2002"})
    for h, unit in zip(FblvlParser.headings, FblvlParser.units):
        table[h].unit = unit
    p = FblvlParser(os.path.join(f.dirname, f.basename), standalone=True)
    table_parsed = p.parse()
    for c in table.colnames:
        assert all(table[c] == table_parsed[c])
    for k in ('footer', 'filename', 'element', 'ion'):
        assert table_parsed.meta[k] == table.meta[k]
    assert 'chianti_version' in table_parsed.meta
    assert 'descriptions' in table_parsed.meta


def test_scupsparser(tmpdir):
    """
    Test parser for scaled upsilon collision strength files
    """
    f = tmpdir.mkdir('ion_sources').join('h_1.scups')
    f.write("""      1      2   7.500e-01   0.000e+00          -1    5    2   4.000e-01
   0.000e+00   2.500e-01   5.000e-01   7.500e-01   1.000e+00
   2.117e-01   3.088e-01   3.327e-01   3.891e-01   5.751e-01
-1
% produced as part of the Arcetri/Cambridge/NRL atomic data base collaboration
% A values:  Parpia, F. A., and Johnson, W. R., 1972, Phys. Rev. A, 26, 1142.
%produced as part of the Arcetri/Cambridge/NRL atomic data base
 collaboration 'CHIANTI' 
  Ken Dere   Wed Jul 26 08:24:54 2000
%energy levels:    NIST Atomic Spectra Database Version 2.0
  http://physics.nist.gov/cgi-bin/AtData/main_asd
%oscillator strengths:  Wiese, W. L., Smith, M. W., and Glennon, B. M., 1966,
 Atomic Transition Probabilities, NSRDS-NBS-4.
%note:  hydrogenic f values
%collision strengths:  Anderson, H., Ballance, C.P., Badnell, N.R., Summers, H.P., 2000, J. Phys. B, 33,1255, revised 2002
%note:  fine structure collision assumes LS coupling or distribution according to statistical weights
  compiled by Ken Dere (NRL)  Thu Feb  7 15:11:55 2002
-1""")
    table = QTable(data=[[1], [2], [0.75], [0.0], [-1.0], [5], [2], [0.4],
                         [np.array([0, 0.25, 0.5, 0.75, 1])],
                         [np.array([0.2117, 0.3088, 0.3327, 0.3891, 0.5751])]],
                   names=ScupsParser.headings,
                   meta={'element': 'h',
                         'ion': 'h_1',
                         'filename': os.path.join(f.dirname, f.basename),
                         'footer': "produced as part of the Arcetri/Cambridge/NRL atomic data base collaboration\nA values:  Parpia, F. A., and Johnson, W. R., 1972, Phys. Rev. A, 26, 1142.\nproduced as part of the Arcetri/Cambridge/NRL atomic data base\ncollaboration 'CHIANTI'\nKen Dere   Wed Jul 26 08:24:54 2000\nenergy levels:    NIST Atomic Spectra Database Version 2.0\nhttp://physics.nist.gov/cgi-bin/AtData/main_asd\noscillator strengths:  Wiese, W. L., Smith, M. W., and Glennon, B. M., 1966,\nAtomic Transition Probabilities, NSRDS-NBS-4.\nnote:  hydrogenic f values\ncollision strengths:  Anderson, H., Ballance, C.P., Badnell, N.R., Summers, H.P., 2000, J. Phys. B, 33,1255, revised 2002\nnote:  fine structure collision assumes LS coupling or distribution according to statistical weights\ncompiled by Ken Dere (NRL)  Thu Feb  7 15:11:55 2002"})
    for h, unit in zip(ScupsParser.headings, ScupsParser.units):
        table[h].unit = unit
    p = ScupsParser(os.path.join(f.dirname, f.basename), standalone=True)
    table_parsed = p.parse()
    for c in table.colnames:
        for row, row_parsed in zip(table[c], table_parsed[c]):
            test = row == row_parsed
            if isinstance(test, np.ndarray):
                assert np.all(test)
            else:
                assert test
    for k in ('footer', 'filename', 'element', 'ion'):
        assert table_parsed.meta[k] == table.meta[k]
    assert 'chianti_version' in table_parsed.meta
    assert 'descriptions' in table_parsed.meta


def test_psplupsparser(tmpdir):
    f = tmpdir.mkdir('ion_sources').join('c_2.psplups')
    f.write("""  1  2  2 0.000e+00 5.830e-04 2.900e+02 6.029e-11 5.015e-10 3.796e-09 7.844e-09 1.140e-08 1.433e-08 1.685e-08 1.906e-08 1.883e-08
 -1
%filename: c_2.psplups
%rates: Foster VJ, Keenan FP, Reid RHG, ADNDT 67, 99, 1997
%energies: From .elvlc file, experimental energies
%comment: The rate coefficients are in units cm^3/s.
%comment: Fits are valid for temperatures 2e3 to 4e5 K.
%produced as part of the Arcetri/Cambridge/NRL 'CHIANTI' atomic data base collaboration
%
% Peter Young   15-May-2001
 -1""")
    table = QTable(data=[[1], [2], [2], [0.0], [0.000583], [290.0],
                         [np.array([6.029e-11, 5.015e-10, 3.796e-09, 7.844e-09, 1.140e-08, 1.433e-08, 1.685e-08, 1.906e-08, 1.883e-08])]],
                   names=PsplupsParser.headings,
                   meta={'element': 'c',
                         'ion': 'c_2',
                         'filename': os.path.join(f.dirname, f.basename),
                         'footer': "filename: c_2.psplups\nrates: Foster VJ, Keenan FP, Reid RHG, ADNDT 67, 99, 1997\nenergies: From .elvlc file, experimental energies\ncomment: The rate coefficients are in units cm^3/s.\ncomment: Fits are valid for temperatures 2e3 to 4e5 K.\nproduced as part of the Arcetri/Cambridge/NRL 'CHIANTI' atomic data base collaboration\nPeter Young   15-May-2001"})
    for h, unit in zip(PsplupsParser.headings, PsplupsParser.units):
        table[h].unit = unit
    p = PsplupsParser(os.path.join(f.dirname, f.basename), standalone=True)
    table_parsed = p.parse()
    for c in table.colnames:
        for row, row_parsed in zip(table[c], table_parsed[c]):
            test = row == row_parsed
            if isinstance(test, np.ndarray):
                assert np.all(test)
            else:
                assert test
    for k in ('footer', 'filename', 'element', 'ion'):
        assert table_parsed.meta[k] == table.meta[k]
    assert 'chianti_version' in table_parsed.meta
    assert 'descriptions' in table_parsed.meta


def test_easplomparser(tmpdir):
    f = tmpdir.mkdir('ion_sources').join('be_2.easplom')
    f.write("""  4  2  1  3  1 7.551e-01 8.649e+00 1.700e+00 9.785e-02 1.007e-01 1.147e-01 1.415e-01 1.749e-01
-1
%file:  be_2.easplom
%excitation autoionization cross section parameter file
 derived from fits to experimental and theoretical data
Dere, K. P., 2007, A&A, 466, 771
ADS ref:  http://adsabs.harvard.edu/abs/2007A%26A...466..771D
 created for CHIANTI database for astrophysical spectroscopy
  created by Ken Dere  (GMU)  Fri Jan 26 12:40:26 2007
-1""")
    table = QTable(data=[[1], [3], [1], [0.7551], [8.649], [1.7],
                         [np.array([0.09785, 0.1007, 0.1147, 0.1415, 0.1749])]],
                   names=EasplomParser.headings,
                   meta={'element': 'be',
                         'ion': 'be_2',
                         'filename': os.path.join(f.dirname, f.basename),
                         'footer': "file:  be_2.easplom\nexcitation autoionization cross section parameter file\nderived from fits to experimental and theoretical data\nDere, K. P., 2007, A&A, 466, 771\nADS ref:  http://adsabs.harvard.edu/abs/2007A%26A...466..771D\ncreated for CHIANTI database for astrophysical spectroscopy\ncreated by Ken Dere  (GMU)  Fri Jan 26 12:40:26 2007"})
    for h, unit in zip(EasplomParser.headings, EasplomParser.units):
        table[h].unit = unit
    p = EasplomParser(os.path.join(f.dirname, f.basename), standalone=True)
    table_parsed = p.parse()
    for c in table.colnames:
        for row, row_parsed in zip(table[c], table_parsed[c]):
            test = row == row_parsed
            if isinstance(test, np.ndarray):
                assert np.all(test)
            else:
                assert test
    for k in ('footer', 'filename', 'element', 'ion'):
        assert table_parsed.meta[k] == table.meta[k]
    assert 'chianti_version' in table_parsed.meta
    assert 'descriptions' in table_parsed.meta
