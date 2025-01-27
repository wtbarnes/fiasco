"""
IDL comparison tests for continuum calculations
"""
import astropy.units as u
import copy
import numpy as np
import pytest

import fiasco

from fiasco.tests.idl.helpers import run_idl_script
from fiasco.util.exceptions import MissingDatasetException


def build_ion_collection(root, temperature, **kwargs):
    all_ions = [fiasco.Ion(i, temperature, hdf5_dbase_root=root, **kwargs) for i in fiasco.list_ions(root)]
    return fiasco.IonCollection(*all_ions)


@pytest.fixture
def ion_input_args():
    return {
        'abundance': 'sun_coronal_1992_feldman_ext',
        'ionization_fraction': 'chianti',
    }


@pytest.fixture
def temperature():
    return np.logspace(5, 8, 61) * u.K


@pytest.fixture
def all_ions(ion_input_args, temperature, hdf5_dbase_root):
    # NOTE: the reason for separating the fixture and the function that generates
    # the collection is so that it is easier to generate a new collection for
    # different inputs
    return build_ion_collection(hdf5_dbase_root, temperature, **ion_input_args)


@pytest.fixture
def idl_input_args(ion_input_args, temperature):
    return {
        'wavelength': np.arange(25, 414, 1) * u.Angstrom,
        'temperature': temperature,
        'density': 1e9*u.cm**(-3),
        **ion_input_args,
    }


def test_idl_compare_free_free(idl_env, all_ions, idl_input_args, dbase_version):
    script = """
    ; set common block
    common elements, abund, abund_ref, ioneq, ioneq_logt, ioneq_ref

    ; read abundance and ionization equilibrium
    abund_file = FILEPATH('{{abundance}}.abund', ROOT_DIR=!xuvtop, SUBDIR='abundance')
    ioneq_file = FILEPATH('{{ionization_fraction}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
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
    ioneq_file = FILEPATH('{{ionization_fraction}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
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
    ioneq_file = FILEPATH('{{ionization_fraction}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
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


def test_idl_compare_free_free_radiative_loss(idl_env, ion_input_args, hdf5_dbase_root, dbase_version):
    script = """
    abund_file = FILEPATH('{{abundance}}.abund', ROOT_DIR=!xuvtop, SUBDIR='abundance')
    ioneq_file = FILEPATH('{{ionization_fraction}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
    ff_rad_loss, temperature, free_free_radiative_loss, abund_file=abund_file, ioneq_file=ioneq_file
    """
    idl_result = run_idl_script(idl_env,
                                script,
                                ion_input_args,
                                ['temperature', 'free_free_radiative_loss'],
                                'freefree_radiative_loss_all_ions',
                                dbase_version,
                                format_func={'free_free_radiative_loss': lambda x: x*u.Unit('erg cm3 s-1'),
                                             'temperature': lambda x: x*u.K})
    all_ions = build_ion_collection(hdf5_dbase_root, idl_result['temperature'], **ion_input_args)
    free_free_radiative_loss_python = all_ions.free_free_radiative_loss()
    assert u.allclose(idl_result['free_free_radiative_loss'],
                      free_free_radiative_loss_python,
                      atol=None,
                      rtol=0.005)


def test_idl_compare_free_bound_radiative_loss(idl_env, ion_input_args, hdf5_dbase_root, dbase_version):
    script = """
    abund_file = FILEPATH('{{abundance}}.abund', ROOT_DIR=!xuvtop, SUBDIR='abundance')
    ioneq_file = FILEPATH('{{ionization_fraction}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
    fb_rad_loss, temperature, free_bound_radiative_loss, abund_file=abund_file, ioneq_file=ioneq_file
    """
    idl_result = run_idl_script(idl_env,
                                script,
                                ion_input_args,
                                ['temperature', 'free_bound_radiative_loss'],
                                'freebound_radiative_loss_all_ions',
                                dbase_version,
                                format_func={'free_bound_radiative_loss': lambda x: x*u.Unit('erg cm3 s-1'),
                                             'temperature': lambda x: x*u.K})
    all_ions = build_ion_collection(hdf5_dbase_root, idl_result['temperature'], **ion_input_args)
    free_bound_radiative_loss_python = all_ions.free_bound_radiative_loss()
    assert u.allclose(idl_result['free_bound_radiative_loss'],
                      free_bound_radiative_loss_python,
                      atol=None,
                      rtol=0.01)


@pytest.mark.xfail()
def test_idl_compare_two_photon(idl_env, all_ions, idl_input_args, dbase_version):
    script = """
    ; set common block
    common elements, abund, abund_ref, ioneq, ioneq_logt, ioneq_ref

    ; read abundance and ionization equilibrium
    abund_file = FILEPATH('{{abundance}}.abund', ROOT_DIR=!xuvtop, SUBDIR='abundance')
    ioneq_file = FILEPATH('{{ionization_fraction}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
    read_abund, abund_file, abund, abund_ref
    read_ioneq, ioneq_file, ioneq_logt, ioneq, ioneq_ref

    ; set temperature and wavelength
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    density = {{ density | to_unit('cm-3') | force_double_precision }}
    wavelength = {{ wavelength | to_unit('Angstrom') | force_double_precision }}

    ; calculate two-photon
    two_photon,temperature,wavelength,two_photon_continuum,edensity=density,/no_setup
    """
    # NOTE: Extend wavelength range for the two-photon test
    new_input_args = copy.deepcopy(idl_input_args)
    new_input_args['temperature'] = 10**np.arange(4, 7.05, 0.05) * u.K
    new_input_args['wavelength'] = np.arange(1,2000,1) * u.Angstrom
    idl_result = run_idl_script(idl_env,
                                script,
                                new_input_args,
                                ['two_photon_continuum'],
                                'twophoton_all_ions',
                                dbase_version,
                                format_func={'two_photon_continuum': lambda x: x*4*np.pi/1e40*u.Unit('erg cm3 s-1 Angstrom-1')})
    two_photon_python = all_ions.two_photon(idl_result['wavelength'], idl_result['density']).squeeze()
    assert u.allclose(idl_result['two_photon_continuum'], two_photon_python, atol=None, rtol=0.005)
