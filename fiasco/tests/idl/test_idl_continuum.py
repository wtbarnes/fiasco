"""
IDL comparison tests for continuum calculations
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

from fiasco.tests.idl.helpers import run_idl_script
from fiasco.util.exceptions import MissingDatasetException


@pytest.fixture
def all_ions(hdf5_dbase_root):
    temperature = np.logspace(5, 8, 61) * u.K
    abundance_name = 'sun_coronal_1992_feldman_ext'
    ioneq_name = 'chianti'
    ion_kwargs = {
        'abundance': abundance_name,
        'ioneq_filename': ioneq_name,
        'hdf5_dbase_root': hdf5_dbase_root,
    }
    all_ions = [fiasco.Ion(i, temperature, **ion_kwargs) for i in fiasco.list_ions(hdf5_dbase_root)]
    return fiasco.IonCollection(*all_ions)


@pytest.fixture
def wavelength():
    return np.arange(25, 414, 1) * u.Angstrom


@pytest.fixture
def idl_input_args(all_ions, wavelength):
    return {
        'wavelength': wavelength,
        'temperature': all_ions.temperature,
        'abundance': all_ions[0]._dset_names['abundance'],
        'ioneq': all_ions[0]._dset_names['ioneq_filename'],
    }


def test_idl_compare_free_free(idl_env, all_ions, idl_input_args, dbase_version):
    script = """
    ; set common block
    common elements, abund, abund_ref, ioneq, ioneq_logt, ioneq_ref

    ; read abundance and ionization equilibrium
    abund_file = FILEPATH('{{abundance}}.abund', ROOT_DIR=!xuvtop, SUBDIR='abundance')
    ioneq_file = FILEPATH('{{ioneq}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
    read_abund, abund_file, abund, abund_ref
    read_ioneq, ioneq_file, ioneq_logt, ioneq, ioneq_ref

    ; set temperature and wavelength
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    wavelength = {{ wavelength | to_unit('Angstrom') | force_double_precision }}

    ; calculate free-free
    freefree, temperature, wavelength, free_free, /no_setup
    """
    idl_result = run_idl_script(idl_env,
                                script,
                                idl_input_args,
                                ['free_free'],
                                'freefree_all_ions',
                                dbase_version,
                                format_func={'free_free': lambda x: x/1e40*u.Unit('erg cm3 s-1 Angstrom-1')})

    free_free_python = all_ions.free_free(idl_result['wavelength'])
    # Compare IDL and Python calculation
    assert u.allclose(idl_result['free_free'], free_free_python, atol=None, rtol=0.005)


def test_idl_compare_free_bound(idl_env, all_ions, idl_input_args, dbase_version):
    script = """
    ; set common block
    common elements, abund, abund_ref, ioneq, ioneq_logt, ioneq_ref

    ; read abundance and ionization equilibrium
    abund_file = FILEPATH('{{abundance}}.abund', ROOT_DIR=!xuvtop, SUBDIR='abundance')
    ioneq_file = FILEPATH('{{ioneq}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
    read_abund, abund_file, abund, abund_ref
    read_ioneq, ioneq_file, ioneq_logt, ioneq, ioneq_ref

    ; set temperature and wavelength
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    wavelength = {{ wavelength | to_unit('Angstrom') | force_double_precision }}

    ; calculate free-bound
    freebound, temperature, wavelength, free_bound, /no_setup
    """
    idl_result = run_idl_script(idl_env,
                                script,
                                idl_input_args,
                                ['free_bound'],
                                'freebound_all_ions',
                                dbase_version,
                                format_func={'free_bound': lambda x: x*4*np.pi/1e40*u.Unit('erg cm3 s-1 Angstrom-1')})
    free_bound_python = all_ions.free_bound(idl_result['wavelength'])
    # Compare IDL and Python calculation
    assert u.allclose(idl_result['free_bound'], free_bound_python, atol=None, rtol=0.005)


def test_idl_compare_free_bound_ion(idl_env, all_ions, idl_input_args, dbase_version):
    script = """
    ; set common block
    common elements, abund, abund_ref, ioneq, ioneq_logt, ioneq_ref

    ; read abundance and ionization equilibrium
    abund_file = FILEPATH('{{abundance}}.abund', ROOT_DIR=!xuvtop, SUBDIR='abundance')
    ioneq_file = FILEPATH('{{ioneq}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
    read_abund, abund_file, abund, abund_ref
    read_ioneq, ioneq_file, ioneq_logt, ioneq, ioneq_ref

    ; set temperature and wavelength
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    wavelength = {{ wavelength | to_unit('Angstrom') | force_double_precision }}

    ; calculate free-bound
    freebound_ion, temperature, wavelength, free_bound, {{ atomic_number }}, {{ ionization_stage }}
    """
    # NOTE: this should be something we parametrize over such that there is a test for each
    # ion rather than all being lumped into a single test that takes a really long time
    for ion in all_ions:
        try:
            free_bound_python = ion.free_bound(idl_input_args['wavelength'])
        except MissingDatasetException:
            continue
        args = {**idl_input_args, 'atomic_number':ion.atomic_number, 'ionization_stage':ion.ionization_stage}
        idl_result = run_idl_script(idl_env,
                                    script,
                                    args,
                                    ['free_bound'],
                                    f'freebound_{ion.atomic_number}_{ion.ionization_stage}',
                                    dbase_version,
                                    format_func={'free_bound': lambda x: x*4*np.pi/1e40*u.Unit('erg cm3 s-1 Angstrom-1')},
                                    write_file=False)
        # Compare IDL and Python calculation
        assert u.allclose(idl_result['free_bound'], free_bound_python, atol=None, rtol=0.006)
