"""
IDL comparison tests for line ratio calculations.
"""
import astropy.units as u
import pytest

import fiasco

from fiasco.tests.idl.helpers import run_idl_script


@pytest.mark.requires_dbase_version('>= 10.1')
def test_eis_fe_xiii_density_ratio_from_idl(idl_env, dbase_version, chianti_idl_version, hdf5_dbase_root):
    script = """
    {% if database_version | version_check('>=', '10.1') %}
    abundance_subdirs = ['abundance', 'archive']
    {% else %}
    abundance_subdirs = 'abundance'
    {% endif %}
    abund_file = FILEPATH('sun_coronal_1992_feldman_ext.abund', ROOT_DIR=!xuvtop, SUBDIR=abundance_subdirs)
    density = findgen(21) * 0.2 + 8.0
    ioneq_file = FILEPATH('chianti.ioneq', ROOT_DIR=!xuvtop, SUBDIR='ioneq')
    defsysv,'!abund_file',abund_file
    defsysv,'!ioneq_file',ioneq_file
    temperature = ch_tmax('fe_13', /log, ioneqname=ioneq_file)
    em = emiss_calc(26, 13, dens=density, temp=temperature, ioneq_file=ioneq_file, abund_file=abund_file, /quiet)

    k = where(em.level1 EQ 1 AND em.level2 EQ 20, nk)
    denominator_index = k[0]
    k = where((em.level1 EQ 3 AND em.level2 EQ 24) OR (em.level1 EQ 3 AND em.level2 EQ 25), nk)
    numerator_index = k

    ratio = total(reform(em[numerator_index].em), 2) / reform(em[denominator_index].em)
    numerator = em[numerator_index].lambda
    denominator = em[denominator_index].lambda
    """
    formatters = {
        'density': lambda x: (10.0**x) * u.cm**(-3),
        'temperature': lambda x: (10.0**x) * u.K,
        'numerator': lambda x: x * u.angstrom,
        'denominator': lambda x: x * u.angstrom,
    }
    idl_result = run_idl_script(
        idl_env,
        script,
        {},
        ['density', 'ratio', 'numerator', 'denominator', 'temperature'],
        'eis_fe_xiii_density_ratio',
        dbase_version,
        chianti_idl_version,
        format_func=formatters,
    )
    ion = fiasco.Ion('Fe XIII', idl_result['temperature'], hdf5_dbase_root=hdf5_dbase_root)
    ratio = fiasco.line_ratio(
        ion,
        idl_result['numerator'],
        idl_result['denominator'],
        idl_result['density'],
        use_two_ion_model=False,
    )
    assert u.allclose(ratio.squeeze(), idl_result['ratio'], rtol=2e-3)
