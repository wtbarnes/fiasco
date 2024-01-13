"""
IDL comparison tests for a few representative contribution functions.
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

from fiasco.tests.idl import save_idl_test_output
from fiasco.util import parse_ion_name

# NOTE: This is necessary because I cannot figure out how to select the contribution function
# by wavelength so I have to run the g_of_t function once with the GUI selector to work out
# what integer index corresponds to which wavelength. For any additional test that is added to
# the list below, this index will need to be computed and added to the list below.
INDEX_WAVE_MAPPING = {
    200.972*u.Angstrom: 297,
    171.073*u.Angstrom: 13946,
    188.497*u.Angstrom: 22527,
    197.862*u.Angstrom: 27964,
    262.984*u.Angstrom: 1415,
}


@pytest.mark.parametrize(('ion_name', 'wavelength'), [
    ('Ca XV', 200.972*u.Angstrom),
    ('Fe IX', 171.073*u.Angstrom),
    ('Fe XI', 188.497*u.Angstrom),
    ('Fe XIV', 197.862*u.Angstrom),
    ('Fe XVI', 262.984*u.Angstrom),
])
def test_idl_compare_goft(idl_env, ascii_dbase_root, hdf5_dbase_root, ion_name, wavelength):
    goft_script = """
    abund_file = '{{ [xuvtop, 'abundance', abundance_file] | join('/') }}'
    ioneq_file = '{{ [xuvtop, 'ioneq', ioneq_file] | join('/') }}'
    density = {{ density | to_unit('cm-3') | log10 | force_double_precision }}
    wave_min = {{ wave_min | to_unit('angstrom') | force_double_precision }}
    wave_max = {{ wave_max | to_unit('angstrom') | force_double_precision }}
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
    Z, iz = parse_ion_name(ion_name)
    # Setup IDl arguments
    abundance_set = 'sun_coronal_1992_feldman_ext'
    ioneq_set = 'chianti'
    density = 1e+10 * u.cm**(-3)
    wave_window = 1 * u.angstrom
    idl_args = {
        'Z': Z,
        'iz': iz,
        'wave_min': wavelength - wave_window,
        'wave_max': wavelength + wave_window,
        'index': INDEX_WAVE_MAPPING[wavelength],
        'density': density,
        'abundance_file': f'{abundance_set}.abund',
        'ioneq_file': f'{ioneq_set}.ioneq',
        'xuvtop': str(ascii_dbase_root),
    }
    # Run IDL code
    res = idl_env.run(goft_script, args=idl_args, save_vars=['temperature','contribution_function'])
    temperature = 10**res['temperature']*u.K
    goft_idl = res['contribution_function'] * u.Unit('erg cm3 s-1')
    # Save the the contribution function to the data directory
    save_idl_test_output({'temperature': temperature, 'goft': goft_idl, 'idl_script': goft_script, **idl_args},
                         f'goft_{Z}_{iz}_{wavelength.to_value("AA"):.3f}',
                         ascii_dbase_root)
    # Run equivalent fiasco code
    ion = fiasco.Ion((Z,iz),
                     temperature,
                     hdf5_dbase_root=hdf5_dbase_root,
                     abundance_filename=abundance_set,
                     ioneq_filename=ioneq_set)
    contribution_func = ion.contribution_function(density)
    idx = np.argmin(np.abs(ion.transitions.wavelength[~ion.transitions.is_twophoton] - wavelength))
    goft_python = contribution_func[:, 0, idx]
    # Find relevant range for comparison. The solutions may diverge many orders of magnitude below
    # the peak of the contribution function, but that is not relevant in assessing meaningful
    # differences between the IDL and Python approaches. Thus, we only assess the differences
    # where the contribution function is above some threshold relative to the peak of the
    # contribution function.
    i_compare = np.where(goft_idl>=goft_idl.max()*1e-10)
    # Compare results
    assert u.allclose(goft_idl[i_compare], goft_python[i_compare], atol=None, rtol=0.05)
