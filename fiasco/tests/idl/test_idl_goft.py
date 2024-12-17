"""
IDL comparison tests for a few representative contribution functions.
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

from fiasco.tests.idl.helpers import run_idl_script
from fiasco.util import parse_ion_name

# NOTE: This is necessary because I cannot figure out how to select the contribution function
# by wavelength so I have to run the g_of_t function once with the GUI selector to work out
# what integer index corresponds to which wavelength. For any additional test that is added to
# the list below, this index will need to be computed and added to the list below.
# NOTE: These indices will change with database version number so they'll need to be updated
# or we will have to find a way to select them automatically.
INDEX_WAVE_MAPPING = {
    'Ca XV 200.972 Angstrom': 297,
    'Fe IX 171.073 Angstrom': 13946,
    'Fe IX 188.496 Angstrom': 18403,
    'Fe XI 188.497 Angstrom': 22527,
    'Fe XIV 197.862 Angstrom': 27964,
    'Fe XVI 262.984 Angstrom': 1415,
}


@pytest.mark.parametrize(('ion_name', 'wavelength'), [
    ('Ca XV', 200.972*u.Angstrom),
    ('Fe IX', 171.073*u.Angstrom),
    ('Fe IX', 188.496*u.Angstrom),
    ('Fe XI', 188.497*u.Angstrom),
    ('Fe XIV', 197.862*u.Angstrom),
    ('Fe XVI', 262.984*u.Angstrom),
])
def test_idl_compare_goft(idl_env, hdf5_dbase_root, dbase_version, ion_name, wavelength):
    goft_script = """
    abund_file = FILEPATH('{{abundance}}.abund', ROOT_DIR=!xuvtop, SUBDIR='abundance')
    ioneq_file = FILEPATH('{{ionization_fraction}}.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
    density = {{ density | to_unit('cm-3') | log10 | force_double_precision }}
    wave_min = {{ (wavelength - wave_window) | to_unit('angstrom') | force_double_precision }}
    wave_max = {{ (wavelength + wave_window) | to_unit('angstrom') | force_double_precision }}
    contribution_function = g_of_t({{ Z }},$
                                    {{ iz }},$
                                    dens=density,$
                                    abund_file=abund_file,$
                                    ioneq_file=ioneq_file,$
                                    {% if index %}index={{ index }},/quiet,${% endif %}
                                    wrange=[wave_min, wave_max])
    ; Call this function to get the temperature array
    read_ioneq,ioneq_file,temperature,ioneq,ref
    """
    # Setup IDl arguments
    Z, iz = parse_ion_name(ion_name)
    input_args = {
        'Z': Z,
        'iz': iz,
        'wavelength': wavelength,
        'wave_window': 1 * u.angstrom,
        'index': INDEX_WAVE_MAPPING[f'{ion_name} {wavelength}'],
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
                                format_func=formatters)
    # Run equivalent fiasco code
    ion = fiasco.Ion(ion_name,
                     idl_result['temperature'],
                     hdf5_dbase_root=hdf5_dbase_root,
                     abundance=idl_result['abundance'],
                     ionization_fraction=idl_result['ioneq'])
    contribution_func = ion.contribution_function(idl_result['density'])
    idx = np.argmin(np.abs(ion.transitions.wavelength[~ion.transitions.is_twophoton] - idl_result['wavelength']))
    # NOTE: Multiply by 0.83 because the fiasco calculation does not include the n_H/n_e ratio
    goft_python = contribution_func[:, 0, idx] * 0.83
    # Find relevant range for comparison. The solutions may diverge many orders of magnitude below
    # the peak of the contribution function, but that is not relevant in assessing meaningful
    # differences between the IDL and Python approaches. Thus, we only assess the differences
    # where the contribution function is above some threshold relative to the peak of the
    # contribution function.
    goft_idl = idl_result['contribution_function']
    i_compare = np.where(goft_idl>=goft_idl.max()*1e-10)
    # Compare results
    assert u.allclose(goft_idl[i_compare], goft_python[i_compare], atol=None, rtol=0.05)
