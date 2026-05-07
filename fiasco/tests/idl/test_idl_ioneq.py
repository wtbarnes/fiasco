"""
IDL comparison tests for ionization equilibria
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

from fiasco.tests.idl.helpers import run_idl_script, version_check
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
@pytest.mark.parametrize(('idl_function', 'fiasco_function', 'rtol'), [
    ('ioniz_rate', 'ionization_rate', 0.05),
    ('recomb_rate', 'recombination_rate', 0.02),
])
def test_rate_from_idl(
    idl_function, fiasco_function, rtol, ion_name, temperature, idl_env, dbase_version, chianti_idl_version, hdf5_dbase_root
):
    script = """
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    rate = {{ rate_function }}('{{ ion_name }}', temperature)
    """
    ion = fiasco.Ion(ion_name, temperature, hdf5_dbase_root=hdf5_dbase_root)
    idl_result = run_idl_script(idl_env,
                                script,
                                {'temperature': temperature, 'ion_name': ion._ion_name, 'rate_function': idl_function},
                                ['rate'],
                                f'{fiasco_function}_{ion.atomic_number}_{ion.ionization_stage}',
                                dbase_version,
                                chianti_idl_version,
                                format_func={'rate': lambda x: x*u.Unit('cm3 s-1')})
    assert u.allclose(idl_result['rate'], getattr(ion, fiasco_function)(), rtol=rtol)


@pytest.mark.requires_dbase_version('>= 11')
@pytest.mark.parametrize('ion_name', ['C II',])
@pytest.mark.parametrize(('idl_function', 'fiasco_function'), [
    ('ch_diel_recomb', 'dielectronic_recombination_rate'),
    ('ch_rad_recomb', 'radiative_recombination_rate'),
])
def test_level_resolved_recombination_rate_from_idl(
    idl_function, fiasco_function, ion_name, temperature, idl_env, dbase_version, chianti_idl_version, hdf5_dbase_root
):
    script = """
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    rate = {{ rate_function }}('{{ ion_name }}', temperature, /level_resolved, /quiet)
    """
    ion = fiasco.Ion(ion_name, temperature, hdf5_dbase_root=hdf5_dbase_root)
    idl_result = run_idl_script(
        idl_env,
        script,
        {'temperature': temperature, 'ion_name': ion._ion_name, 'rate_function': idl_function},
        ['rate'],
        f'level_resolved_{fiasco_function}_{ion.atomic_number}_{ion.ionization_stage}',
        dbase_version,
        chianti_idl_version,
        format_func={'rate': lambda x: x*u.Unit('cm3 s-1')}
    )
    assert u.allclose(idl_result['rate'],
                      getattr(ion, fiasco_function)(level_resolved=True).T,
                      rtol=1e-6)


@pytest.mark.requires_dbase_version('>= 11')
@pytest.mark.parametrize('ion_name', ['C II',])
@pytest.mark.parametrize(('filetype', 'fiasco_function'), [
    ('dilvl', 'direct_ionization_rate'),
    ('ealvl', 'excitation_autoionization_rate'),
])
def test_level_resolved_ionization_rate_from_idl(
    filetype, fiasco_function, ion_name, temperature, idl_env, dbase_version, chianti_idl_version, hdf5_dbase_root
):
    script = """
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    rate = ch_ioniz_rate_lr(!xuvtop+'/{{ file_name }}', temp_in=temperature)
    rate = rate.rates
    """
    ion = fiasco.Ion(ion_name, temperature, hdf5_dbase_root=hdf5_dbase_root)
    idl_result = run_idl_script(
        idl_env,
        script,
        {'temperature': temperature, 'file_name': f'{ion.atomic_symbol.lower()}/{ion._ion_name}/{ion._ion_name}.{filetype}'},
        ['rate'],
        f'level_resolved_{fiasco_function}_{ion.atomic_number}_{ion.ionization_stage}',
        dbase_version,
        chianti_idl_version,
        format_func={'rate': lambda x: x*u.Unit('cm3 s-1')}
    )
    assert u.allclose(idl_result['rate'],
                      getattr(ion, fiasco_function)(level_resolved=True).T,
                      rtol=1e-6)


# NOTE: The list of ions here is motivated by the need to test the different cases for different
# isoelectronic sequences when calculating the suppression factor.
@pytest.mark.requires_dbase_version('>= 10')
@pytest.mark.parametrize('ion_name', [
    'Si IV',
    'S V',
    'Fe XII',
    'Fe XIII',
    'Fe XIV',
    'O II',
    'O III',
    'O IV',
    'O V',
    'O VI',
    'O VII',
    'O VIII',
])
def test_dielectronic_recombination_suppression_factor_from_idl(ion_name, idl_env, dbase_version, chianti_idl_version, hdf5_dbase_root):
    script = """
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    density = {{ density | to_unit('cm-3') | force_double_precision }}
    suppression = fltarr(n_elements(density))
    for i=0,n_elements(density)-1 do begin
        rate = ch_dr_suppress('{{ion_name}}', temperature, density=density[i], q0=q0, xa=xa, s=s)
        suppression[i] = s
    end
    """
    # Derived from the example suppression factor calculation shown in Fig. 2 of Nikolic et al. (2018)
    temperature = 10**4.5 * u.K
    ion = fiasco.Ion(ion_name, temperature, hdf5_dbase_root=hdf5_dbase_root)
    # There is a bug in the calculation of the suppression factors for H- and He-like ions that was
    # fixed in version 11.0.2 of the CHIANTI software
    density = np.logspace(0,15,30)*u.cm**(-3)
    suppression = ion._dielectronic_recombination_suppression(density, couple_density_to_temperature=False)
    idl_result = run_idl_script(idl_env,
                                script,
                                {'temperature': temperature, 'density': density, 'ion_name': ion._ion_name},
                                ['suppression'],
                                f'dielectronic_recombination_suppression_factor_{ion.atomic_number}_{ion.ionization_stage}',
                                dbase_version,
                                chianti_idl_version)
    if (ion.isoelectronic_sequence in ('H', 'He')) and version_check(idl_result['chianti_idl_version'], '<', '11.0.2'): # pragma: no cover
        pytest.skip('Skipping dielectronic recombination suppression test for H- and He-like ions due '
                    'to a bug in the IDL software in versions prior to v11.0.2.')
    u.allclose(idl_result['suppression'], suppression.squeeze(), rtol=0.01)


@pytest.mark.requires_dbase_version('>= 11')
@pytest.mark.parametrize('ion_name', ['C II',])
@pytest.mark.parametrize('atmosphere_model', [None, 'avrett_atmosphere'])
def test_advanced_model_rates_from_idl(
    ion_name, atmosphere_model, temperature, idl_env, dbase_version, chianti_idl_version, hdf5_dbase_root
):
    # Compare density-dependent rates for advanced model ions
    # This tests the output of ch_adv_model_rates.pro which outputs
    # both the ionization and recombination rates
    # should also test with and without a model atmosphere
    script = """
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    density = {{ density | to_unit('cm-3') | force_double_precision }}
    {% if atmosphere_model %}
    atmosphere_file = FILEPATH('{{ atmosphere_model }}.dat', $
                               ROOT_DIR=!xuvtop, $
                               SUBDIR=['ancillary_data', 'advanced_models', 'model_atmospheres'])
    params = ch_adv_model_setup(temperature, /ct, atmosphere_file=atmosphere_file)
    rate = ch_adv_model_rates('{{ ion_name }}', temperature, density, model_atm=params.ct_model)
    {% else %}
    rate = ch_adv_model_rates('{{ ion_name }}', temperature, density)
    {% endif %}
    ionization_rate = rate.final_ioniz
    recombination_rate = rate.final_recomb
    """
    ion = fiasco.Ion(ion_name, temperature, hdf5_dbase_root=hdf5_dbase_root)
    density = 1e15*u.K*u.cm**(-3)/temperature
    input_args = {
        'temperature': temperature,
        'density': density,
        'ion_name': ion._ion_name,
        'atmosphere_model': atmosphere_model,
    }
    atm_name = 'no_atm' if atmosphere_model is None else atmosphere_model
    file_name = f'advanced_model_rates_{atm_name}_{ion.atomic_number}_{ion.ionization_stage}'
    idl_result = run_idl_script(
        idl_env,
        script,
        input_args,
        ['ionization_rate', 'recombination_rate'],
        file_name,
        dbase_version,
        chianti_idl_version,
        format_func={'ionization_rate': lambda x: x*u.Unit('cm3 s-1'),
                     'recombination_rate': lambda x: x*u.Unit('cm3 s-1')}
    )
    assert u.allclose(
        idl_result['ionization_rate'], ion.ionization_rate(density=density).T, rtol=1e-6
    )
    assert u.allclose(
        idl_result['recombination_rate'], ion.recombination_rate(density=density).T, rtol=1e-6
    )


def test_density_dependent_ionization_fraction_from_idl():
    # Compare density-dependent ionization fractions for a all ions of a few different elements
    # included in the advanced model list.
    # This tests the output of ch_calc_ioneq.pro
    ...
