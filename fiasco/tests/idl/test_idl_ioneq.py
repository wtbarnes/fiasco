"""
IDL comparison tests for ionization equilibria
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

from fiasco.tests.idl.helpers import run_idl_script
from fiasco.util import parse_ion_name


@pytest.mark.parametrize('ion_name', [
    'H 1',
    'Fe 5',
    'Fe 16',
    'Fe 18',
    'Fe 27',
    'C 1',
    'C 2',
    'C 3',
    'Ca 2',
    'Fe 20',
])
def test_ionization_fraction_from_idl(ion_name, idl_env, dbase_version, chianti_idl_version, hdf5_dbase_root):
    Z, iz = parse_ion_name(ion_name)
    script = """
    ioneq_file = FILEPATH('{{ ionization_fraction }}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
    read_ioneq, ioneq_file, temperature, ioneq, ioneq_ref
    ioneq = ioneq[*,{{Z-1}},{{iz-1}}]
    """
    formatters = {'temperature': lambda x: 10**x*u.K,
                  'ioneq': lambda x: x*u.dimensionless_unscaled}
    idl_result = run_idl_script(idl_env,
                                script,
                                {'ionization_fraction': 'chianti', 'Z': Z, 'iz': iz},
                                ['temperature', 'ioneq'],
                                f'ioneq_{Z}_{iz}',
                                dbase_version,
                                chianti_idl_version,
                                format_func=formatters)
    ion = fiasco.Ion(ion_name,
                     idl_result['temperature'],
                     hdf5_dbase_root=hdf5_dbase_root,
                     ionization_fraction=idl_result['ionization_fraction'])
    assert u.allclose(idl_result['ioneq'], ion.ionization_fraction, rtol=0.0, atol=1e-5)


@pytest.fixture
def temperature():
    return 10**np.arange(5,8,0.05) * u.K


@pytest.mark.parametrize('ion_name', [
    'Ne VII',
    'Al III',
    'Fe I',
    'Fe XI',
    'Fe XII',
    'Fe XVIII',
    'Fe XXIV',
])
def test_ionization_rate_from_idl(ion_name, temperature, idl_env, dbase_version, chianti_idl_version, hdf5_dbase_root):
    script = """
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    rate = ioniz_rate('{{ ion_name }}', temperature)
    """
    ion = fiasco.Ion(ion_name, temperature, hdf5_dbase_root=hdf5_dbase_root)
    idl_result = run_idl_script(idl_env,
                                script,
                                {'temperature': temperature, 'ion_name': ion._ion_name},
                                ['rate'],
                                f'ionization_rate_{ion.atomic_number}_{ion.ionization_stage}',
                                dbase_version,
                                chianti_idl_version,
                                format_func={'rate': lambda x: x*u.Unit('cm3 s-1')})
    assert u.allclose(idl_result['rate'], ion.ionization_rate, rtol=0.05)


@pytest.mark.parametrize('ion_name', [
    'Ne VII',
    'Al III',
    'Fe I',
    'Fe XI',
    'Fe XII',
    'Fe XVIII',
    'Fe XXIV',
])
def test_recombination_rate_from_idl(ion_name, temperature, idl_env, dbase_version, chianti_idl_version, hdf5_dbase_root):
    script = """
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    rate = recomb_rate('{{ ion_name }}', temperature)
    """
    ion = fiasco.Ion(ion_name, temperature, hdf5_dbase_root=hdf5_dbase_root)
    idl_result = run_idl_script(idl_env,
                                script,
                                {'temperature': temperature, 'ion_name': ion._ion_name},
                                ['rate'],
                                f'recombination_rate_{ion.atomic_number}_{ion.ionization_stage}',
                                dbase_version,
                                chianti_idl_version,
                                format_func={'rate': lambda x: x*u.Unit('cm3 s-1')})
    assert u.allclose(idl_result['rate'], ion.recombination_rate, rtol=0.02)
