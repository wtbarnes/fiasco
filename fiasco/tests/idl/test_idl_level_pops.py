"""
IDL comparison tests related to the level populations calculation.
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

from fiasco.tests.idl.helpers import run_idl_script
from fiasco.util import parse_ion_name


@pytest.fixture(params=['fe_23', 'ca_18'])
def idl_two_ion_model(request, idl_env, dbase_version, chianti_idl_version):
    ion_name = request.param
    temperature = np.logspace(5.5, 7.5, 10) * u.K
    abundance = 'sun_coronal_1992_feldman_ext'
    ionization_fraction = 'chianti'
    script = """
    temp = {{ temperature | to_unit('K') | force_double_precision }}
    {% if database_version | version_check('>=', '10.1') %}
    abundance_subdirs = ['abundance', 'archive']
    {% else %}
    abundance_subdirs = 'abundance'
    {% endif %}
    abund_file = FILEPATH('{{abundance}}.abund', ROOT_DIR=!xuvtop, SUBDIR=abundance_subdirs)
    ioneq_file = FILEPATH('{{ionization_fraction}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')

    ion_1 = ch_load_ion_rates('{{ ion_name }}',$
                              temp,$
                              abund_file=abund_file,$
                              ioneq_file=ioneq_file,$
                              /verbose)
    ion_2 = ch_load_ion_rates('{{ next_ion_name }}',$
                              temp,$
                              abund_file=abund_file,$
                              ioneq_file=ioneq_file,$
                              /verbose)
    two_ion = ch_load_2ion_rates(ion_1,ion_2)

    ionization = two_ion.ioniz
    radiative_recombination = two_ion.rr
    autoionization = two_ion.ai
    dielectronic_capture = two_ion.dc
    dielectronic_recombination = two_ion.dr
    """
    Z, iz = parse_ion_name(ion_name)
    next_ion_name = f'{ion_name.split("_")[0]}_{iz+1}'
    input_args = {
        'ion_name': ion_name,
        'next_ion_name': next_ion_name,
        'temperature': temperature,
        'abundance': abundance,
        'ionization_fraction': ionization_fraction,
    }
    outputs = [
        'ionization',
        'radiative_recombination',
        'autoionization',
        'dielectronic_capture',
        'dielectronic_recombination',
    ]
    formatters = {
        'ionization': lambda x: u.Quantity(np.swapaxes(x.T, 1, 2), 'cm3 s-1'),
        'radiative_recombination': lambda x: u.Quantity(np.swapaxes(x.T, 1, 2), 'cm3 s-1'),
        'autoionization': lambda x: u.Quantity(x, 's-1'),
        'dielectronic_capture': lambda x: u.Quantity(np.swapaxes(x.T, 1, 2), 'cm3 s-1'),
        'dielectronic_recombination': lambda x: u.Quantity(np.swapaxes(x.T, 1, 2), 'cm3 s-1'),
    }
    idl_result = run_idl_script(idl_env,
                                script,
                                input_args,
                                outputs,
                                f'two_ion_rate_matrices_{Z}_{iz}',
                                dbase_version,
                                chianti_idl_version,
                                format_func=formatters,
                                write_file=False)
    return idl_result


@pytest.fixture
def ion(idl_two_ion_model, hdf5_dbase_root):
    return fiasco.Ion(idl_two_ion_model['ion_name'],
                      idl_two_ion_model['temperature'],
                      abundance=idl_two_ion_model['abundance'],
                      ionization_fraction=idl_two_ion_model['ionization_fraction'],
                      hdf5_dbase_root=hdf5_dbase_root)


@pytest.mark.requires_dbase_version('>= 9')
def test_idl_compare_rate_matrix_ionization(ion, idl_two_ion_model):
    assert u.allclose(ion._rate_matrix_ionization,
                      idl_two_ion_model['ionization'],
                      rtol=0.05)


@pytest.mark.requires_dbase_version('>= 9')
def test_idl_compare_rate_matrix_radiative_recombination(ion, idl_two_ion_model):
    assert u.allclose(ion._rate_matrix_radiative_recombination,
                      idl_two_ion_model['radiative_recombination'],
                      rtol=1e-5)


@pytest.mark.requires_dbase_version('>= 9')
def test_idl_compare_rate_matrix_autoionization(ion, idl_two_ion_model):
    assert u.allclose(ion._rate_matrix_autoionization,
                      idl_two_ion_model['autoionization'],
                      rtol=1e-7)


@pytest.mark.requires_dbase_version('>= 9')
def test_idl_compare_rate_matrix_dielectronic_capture(ion, idl_two_ion_model):
    assert u.allclose(ion._rate_matrix_dielectronic_capture,
                      idl_two_ion_model['dielectronic_capture'],
                      rtol=2e-4)


@pytest.mark.requires_dbase_version('>= 9')
def test_idl_compare_rate_matrix_dielectronic_recombination(ion, idl_two_ion_model):
    assert u.allclose(ion._rate_matrix_dielectronic_recombination,
                      idl_two_ion_model['dielectronic_recombination'],
                      rtol=2e-3)
