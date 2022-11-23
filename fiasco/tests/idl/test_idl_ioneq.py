"""
IDL comparison tests for ionization equilibria
"""
import astropy.units as u
import numpy as np
import pathlib
import pytest

import fiasco


@pytest.fixture
def ioneq_from_idl(idl_env, ascii_dbase_root):
    script = """
    ioneqfile="{{ ioneq_filename }}"
    read_ioneq, ioneqfile, ioneq_logt, ioneq, ioneq_ref
    """
    args = {
        'ioneq_filename': ascii_dbase_root / 'ioneq' / f'chianti.ioneq',
    }
    res_idl = idl_env.run(script, args=args, verbose=True)
    return res_idl


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
])
def test_ioneq_from_idl(ion_name, ioneq_from_idl, hdf5_dbase_root):
    temperature = 10**ioneq_from_idl['ioneq_logt'] * u.K
    ion_kwargs = {
        'ioneq_filename': pathlib.Path(str(ioneq_from_idl['ioneqfile'])).stem,
        'hdf5_dbase_root': hdf5_dbase_root,
    }
    ion = fiasco.Ion(ion_name, temperature, **ion_kwargs)
    ioneq_idl = ioneq_from_idl['ioneq'][ion.charge_state, ion.atomic_number-1, :]
    ioneq_python = ion.ioneq
    assert u.allclose(ioneq_idl, ioneq_python, rtol=0.0, atol=1e-5)
