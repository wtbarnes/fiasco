"""
IDL comparison tests for continuum calculations
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco


@pytest.fixture
def all_ions(hdf5_dbase_root):
    temperature = np.logspace(5, 8, 61) * u.K
    abundance_name = 'sun_coronal_1992_feldman_ext'
    ioneq_name = 'chianti'
    ion_kwargs = {
        'abundance_filename': abundance_name,
        'ioneq_filename': ioneq_name,
        'hdf5_dbase_root': hdf5_dbase_root,
    }
    all_ions = [fiasco.Ion(i, temperature, **ion_kwargs) for i in fiasco.list_ions(hdf5_dbase_root)]
    return fiasco.IonCollection(*all_ions)


@pytest.fixture
def wavelength():
    return np.arange(25, 414, 1) * u.Angstrom


def test_idl_compare_free_free(idl_env, ascii_dbase_root, all_ions, wavelength):
    ff_python = all_ions.free_free(wavelength).to_value('erg cm3 s-1 Angstrom-1')
    # Run IDL free-free calculation
    script = """
    ; set common block
    common elements, abund, abund_ref, ioneq, ioneq_logt, ioneq_ref

    ; read abundance and ionization equilibrium
    abundfile='{{ abundance_filename }}'
    read_abund, abundfile, abund, abund_ref
    ioneqfile="{{ ioneq_filename }}"
    read_ioneq, ioneqfile, ioneq_logt, ioneq, ioneq_ref

    ; set temperature and wavelength
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    wavelength = {{ wavelength | to_unit('Angstrom') | force_double_precision }}

    ; calculate free-free
    freefree, temperature, wavelength, ff, /no_setup

    ; The output of freefree is scaled by 10^40
    ff = ff/1d40
    """
    abundance_name = all_ions[0]._dset_names['abundance_filename']
    ioneq_name = all_ions[0]._dset_names['ioneq_filename']
    args = {
        'wavelength': wavelength,
        'temperature': all_ions.temperature,
        'abundance_filename': ascii_dbase_root / 'abundance' / f'{abundance_name}.abund',
        'ioneq_filename': ascii_dbase_root / 'ioneq' / f'{ioneq_name}.ioneq',
    }
    res_idl = idl_env.run(script, args=args, verbose=True)
    ff_idl = res_idl['ff']
    # Compare IDL and Python calculation
    assert u.allclose(ff_idl, ff_python, atol=0.0, rtol=0.005)


def test_idl_compare_free_bound(idl_env, ascii_dbase_root, all_ions, wavelength):
    fb_python = all_ions.free_bound(wavelength).to_value('erg cm3 s-1 Angstrom-1')
    # Run IDL free-free calculation
    script = """
    ; set common block
    common elements, abund, abund_ref, ioneq, ioneq_logt, ioneq_ref

    ; read abundance and ionization equilibrium
    abundfile='{{ abundance_filename }}'
    read_abund, abundfile, abund, abund_ref
    ioneqfile="{{ ioneq_filename }}"
    read_ioneq, ioneqfile, ioneq_logt, ioneq, ioneq_ref

    ; set temperature and wavelength
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    wavelength = {{ wavelength | to_unit('Angstrom') | force_double_precision }}

    ; calculate free-bound
    freebound, temperature, wavelength, fb, /no_setup

    ; The output of freefree is scaled by 10^40
    fb = fb/1d40
    """
    abundance_name = all_ions[0]._dset_names['abundance_filename']
    ioneq_name = all_ions[0]._dset_names['ioneq_filename']
    args = {
        'wavelength': wavelength,
        'temperature': all_ions[0].temperature,
        'abundance_filename': ascii_dbase_root / 'abundance' / f'{abundance_name}.abund',
        'ioneq_filename': ascii_dbase_root / 'ioneq' / f'{ioneq_name}.ioneq',
    }
    res_idl = idl_env.run(script, args=args, verbose=True)
    fb_idl = res_idl['fb']
    # Compare IDL and Python calculation
    assert u.allclose(fb_idl, fb_python, atol=0.0, rtol=0.005)
