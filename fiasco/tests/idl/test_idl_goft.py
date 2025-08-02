"""
IDL comparison tests for a few representative contribution functions.
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

from fiasco.tests.idl.helpers import run_idl_script, version_check
from fiasco.util import parse_ion_name


@pytest.mark.parametrize(('ion_name', 'wavelength'), [
    ('Ca XV', 200.972*u.Angstrom),
    ('Fe IX', 171.073*u.Angstrom),
    ('Fe IX', 188.496*u.Angstrom),
    ('Fe XI', 188.497*u.Angstrom),
    ('Fe XIV', 197.862*u.Angstrom),
    ('Fe XVI', 262.984*u.Angstrom),
])
def test_idl_compare_goft(idl_env, hdf5_dbase_root, dbase_version, chianti_idl_version, ion_name, wavelength):
    goft_script = """
    {% if database_version | version_check('>=', '10.1') %}
    abundance_subdirs = ['abundance', 'archive']
    {% else %}
    abundance_subdirs = 'abundance'
    {% endif %}
    abund_file = FILEPATH('{{abundance}}.abund', ROOT_DIR=!xuvtop, SUBDIR=abundance_subdirs)
    ioneq_file = FILEPATH('{{ionization_fraction}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')

    {% if chianti_idl_version | version_check('>=', 9) %}
    ; Set ioneq_file this way to get around a bug that always causes the GUI picker to pop up
    ; even when the file is specified. Seems to only happen in v9.
    defsysv,'!ioneq_file',ioneq_file
    {% endif %}

    density = {{ density | to_unit('cm-3') | log10 | force_double_precision }}
    wave_min = {{ (wavelength - wave_window) | to_unit('angstrom') | force_double_precision }}
    wave_max = {{ (wavelength + wave_window) | to_unit('angstrom') | force_double_precision }}

    contribution_function = g_of_t({{ Z }},$
                                   {{ iz }},$
                                   dens=density,$
                                   abund_file=abund_file,$
                                   ioneq_file=ioneq_file,$
                                   wrange=[wave_min, wave_max])
    ; Call this function to get the temperature array
    read_ioneq,ioneq_file,temperature,ioneq,ref

    {% if chianti_idl_version | version_check('>=', 9) %}
    defsysv,'!ioneq_file',''
    {% endif %}
    """
    # Setup IDl arguments
    Z, iz = parse_ion_name(ion_name)
    input_args = {
        'Z': Z,
        'iz': iz,
        'wavelength': wavelength,
        'wave_window': 1 * u.angstrom,
        'density': 1e+10 * u.cm**(-3),
        'abundance': 'sun_coronal_1992_feldman_ext',
        'ionization_fraction': 'chianti',
    }
    formatters = {'temperature': lambda x: 10**x*u.K,
                  'contribution_function': lambda x: x*u.Unit('erg cm3 s-1')}
    idl_result = run_idl_script(idl_env,
                                goft_script,
                                input_args,
                                ['temperature', 'contribution_function'],
                                f'goft_{Z}_{iz}_{wavelength.to_value("AA"):.3f}',
                                dbase_version,
                                chianti_idl_version,
                                format_func=formatters)
    # Run equivalent fiasco code
    ion = fiasco.Ion(ion_name,
                     idl_result['temperature'],
                     hdf5_dbase_root=hdf5_dbase_root,
                     abundance=idl_result['abundance'],
                     ionization_fraction=idl_result['ionization_fraction'])
    contribution_func = ion.contribution_function(idl_result['density'])
    idx = np.argmin(np.abs(ion.transitions.wavelength[ion.transitions.is_bound_bound] - idl_result['wavelength']))
    # NOTE: fiasco does not include the n_H/n_e ratio
    if version_check(idl_result['chianti_idl_version'], '<', 9):
        # Prior to v9, the CHIANTI IDL software assumed this ratio was a constant 0.83
        n_H_n_e = 0.83
    else:
        # Later versions use the actual temperature-dependent proton-to-electron ratio
        n_H_n_e = ion.proton_electron_ratio
    goft_python = contribution_func[:, 0, idx] * n_H_n_e
    # Find relevant range for comparison. The solutions may diverge many orders of magnitude below
    # the peak of the contribution function, but that is not relevant in assessing meaningful
    # differences between the IDL and Python approaches. Thus, we only assess the differences
    # where the contribution function is above some threshold relative to the peak of the
    # contribution function.
    goft_idl = idl_result['contribution_function']
    i_compare = np.where(goft_idl>=goft_idl.max()*1e-3)
    # Compare results
    assert u.allclose(goft_idl[i_compare], goft_python[i_compare], atol=None, rtol=0.05)
