"""
IDL comparison tests for ionization equilibria
"""
import astropy.units as u
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
def test_ionization_fraction_from_idl(ion_name, idl_env, dbase_version, hdf5_dbase_root):
    Z, iz = parse_ion_name(ion_name)
    script = """
    ioneq_file = FILEPATH('{{ioneq_filename}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
    read_ioneq, ioneq_file, temperature, ioneq, ioneq_ref
    ioneq = ioneq[*,{{Z-1}},{{iz-1}}]
    """
    formatters = {'temperature': lambda x: 10**x*u.K,
                  'ioneq': lambda x: x*u.dimensionless_unscaled}
    idl_result = run_idl_script(idl_env,
                                script,
                                {'ioneq_filename': 'chianti', 'Z': Z, 'iz': iz},
                                ['temperature', 'ioneq'],
                                f'ioneq_{Z}_{iz}',
                                dbase_version,
                                format_func=formatters)
    ion = fiasco.Ion(ion_name,
                     idl_result['temperature'],
                     hdf5_dbase_root=hdf5_dbase_root,
                     ionization_fraction='chianti')
    assert u.allclose(idl_result['ioneq'], ion.ionization_fraction, rtol=0.0, atol=1e-5)
