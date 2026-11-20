"""
IDL comparison tests related to the level populations calculation.
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

from fiasco.tests.idl.helpers import run_idl_script

# NOTE: This runs the same calculation regardless of which outputs are saved so that the
# level populations and the rate matrices used to compute them are always consistent.
POP_SOLVER_SCRIPT = """
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
density = {{ density | to_unit('cm-3') | force_double_precision }}
temperature = {{ temperature | to_unit('K') | force_double_precision }}
ion_obj = ch_setup_ion('{{ion_name}}', ioneq_file=ioneq_file, abund_file=abund_file, /noprot)
pop_solver, ion_obj, temperature, density, pops, /no_auto, /noionrec, out_rates=rates
elec_collision_rates = rates.qq
rad_decay_rates = rates.aa
"""
FORMATTERS = {
    'pops': lambda x: u.Quantity(x.squeeze().T),
    'elec_collision_rates': lambda x: u.Quantity(np.swapaxes(x.T, 1, 2), 'cm3 s-1'),
    'rad_decay_rates': lambda x: u.Quantity(x, 's-1'),
}


@pytest.fixture
def temperature():
    return 10**np.arange(4, 9, 0.1)*u.K


@pytest.fixture
def density():
    return 1e11 * u.cm**(-3)


# These ions are chosen specifically as their level populations have previously
# shown large differences compared to the IDL results.
@pytest.fixture(params=['C IV', 'Cr XXII', 'Na II', 'Si IV', 'O III', 'Ne IX', 'Mg VII', 'Fe XXVI'])
def ion(request, temperature, hdf5_dbase_root):
    ion = fiasco.Ion(request.param,
                     temperature,
                     abundance='sun_coronal_1992_feldman_ext',
                     ionization_fraction='chianti',
                     hdf5_dbase_root=hdf5_dbase_root)
    if ion.n_levels == 0:
        pytest.skip(f'The test database does not include an atomic model for {ion.ion_name_roman}.')
    return ion


def run_pop_solver(ion, density, name, outputs, idl_env, dbase_version, chianti_idl_version, **kwargs):
    input_args = {
        'density': density,
        'temperature': ion.temperature,
        'abundance': ion._dset_names['abundance'],
        'ionization_fraction': ion._dset_names['ionization_fraction'],
        'ion_name': ion._ion_name,
    }
    return run_idl_script(idl_env,
                          POP_SOLVER_SCRIPT,
                          input_args,
                          outputs,
                          f'{name}_{ion.atomic_number}_{ion.ionization_stage}',
                          dbase_version,
                          chianti_idl_version,
                          format_func={k: FORMATTERS[k] for k in outputs},
                          **kwargs)


@pytest.fixture
def idl_level_populations(ion, density, idl_env, dbase_version, chianti_idl_version):
    # NOTE: Only the level populations are small enough to store. This means that this
    # comparison, unlike those of the rate matrices, does not require an IDL installation.
    return run_pop_solver(ion,
                          density,
                          'level_populations',
                          ['pops'],
                          idl_env,
                          dbase_version,
                          chianti_idl_version)


@pytest.fixture
def idl_rate_matrices(ion, density, idl_env, dbase_version, chianti_idl_version):
    return run_pop_solver(ion,
                          density,
                          'level_populations_rate_matrices',
                          ['elec_collision_rates', 'rad_decay_rates'],
                          idl_env,
                          dbase_version,
                          chianti_idl_version,
                          write_file=False)


def test_idl_compare_level_populations(ion, density, idl_level_populations):
    level_pops = ion.level_populations(density,
                                       include_level_resolved_rate_correction=False,
                                       use_two_ion_model=False,
                                       include_protons=False).squeeze()
    assert u.allclose(level_pops, idl_level_populations['pops'], atol=1e-15, rtol=0.01)


def test_idl_compare_radiative_decay_matrices(ion, idl_rate_matrices):
    assert u.allclose(ion._rate_matrix_radiative_decay,
                      idl_rate_matrices['rad_decay_rates'],
                      rtol=1e-7)


def test_idl_compare_collisional_electron_rate_matrices(ion, idl_rate_matrices):
    # NOTE: IDL interpolates the scaled collision strengths in single precision. At high
    # temperatures, where the interpolated value can be many orders of magnitude smaller
    # than the spline nodes, this gives errors of a few percent in rates that are too small
    # to matter (the largest difference seen across all ions is 5e-18 cm3 s-1).
    assert u.allclose(ion._rate_matrix_collisional_electron,
                      idl_rate_matrices['elec_collision_rates'],
                      rtol=0.01,
                      atol=1e-16*u.cm**3/u.s)


@pytest.mark.parametrize('ion_name', ['Fe X', 'Fe XII', 'Ca XV'])
def test_idl_compare_proton_collision_excitation_rate(ion_name,
                                                      temperature,
                                                      idl_env,
                                                      dbase_version,
                                                      chianti_idl_version,
                                                      hdf5_dbase_root):
    script = """
    temperature = {{ temperature | to_unit('K') | force_double_precision }}
    zion2filename, {{ atomic_number }}, {{ ionization_stage }}, filename
    read_splups, filename+'.psplups', protons, /prot
    rate = dblarr(n_elements(protons), n_elements(temperature))
    FOR i=0, n_elements(protons)-1 DO BEGIN
        descale_all, temperature, protons, i, rate_i
        rate[i, *] = rate_i
    ENDFOR
    """
    ion = fiasco.Ion(ion_name, temperature, hdf5_dbase_root=hdf5_dbase_root)
    input_args = {
        'temperature': temperature,
        'atomic_number': ion.atomic_number,
        'ionization_stage': ion.ionization_stage,
    }
    idl_result = run_idl_script(idl_env,
                                script,
                                input_args,
                                ['rate'],
                                f'proton_excitation_rate_{ion.atomic_number}_{ion.ionization_stage}',
                                dbase_version,
                                chianti_idl_version,
                                format_func={'rate': lambda x: x*u.Unit('cm3 s-1')})
    assert u.allclose(ion.proton_collision_excitation_rate, idl_result['rate'], rtol=1e-3)
