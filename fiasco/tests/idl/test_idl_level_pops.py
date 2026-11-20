"""
IDL comparison tests related to the level populations calculation.
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

from fiasco.tests.idl.helpers import run_idl_script


@pytest.fixture
def density():
    return 1e11 * u.cm**(-3)


# These ions are chosen specifically as their level populations have previously
# shown large differences compared to the IDL results.
@pytest.fixture(params=['C IV', 'Cr XXII', 'Na II', 'Si IV', 'O III', 'Ne IX', 'Mg VII'])
def ion(request, hdf5_dbase_root):
    temperature = 10**np.arange(4, 9, 0.1)*u.K
    return fiasco.Ion(request.param,
                      temperature,
                      abundance='sun_coronal_1992_feldman_ext',
                      ionization_fraction='chianti',
                      hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def idl_result(ion, density, idl_env, dbase_version, chianti_idl_version):
    script = """
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
    input_args={
        'density': density,
        'temperature': ion.temperature,
        'abundance': ion._dset_names['abundance'],
        'ionization_fraction': ion._dset_names['ionization_fraction'],
        'ion_name': ion._ion_name,
    }
    outputs = ['pops', 'elec_collision_rates', 'rad_decay_rates']
    formatters = {
        'pops': lambda x: u.Quantity(x.squeeze().T),
        'elec_collision_rates': lambda x: u.Quantity(np.swapaxes(x.T, 1, 2), 'cm3 s-1'),
        'rad_decay_rates': lambda x: u.Quantity(x, 's-1'),
    }
    idl_result = run_idl_script(idl_env,
                                script,
                                input_args,
                                outputs,
                                f'level_populations_{ion.atomic_number}_{ion.ionization_stage}',
                                dbase_version,
                                chianti_idl_version,
                                format_func=formatters,
                                write_file=False)
    return idl_result


def test_idl_compare_level_populations(ion, density, idl_result):
    level_pops = ion.level_populations(density,
                                       include_level_resolved_rate_correction=False,
                                       use_two_ion_model=False,
                                       include_protons=False).squeeze()
    # Apply thresholds such that comparisons are focused around the relevant temperature range.
    ionfrac_threshold = 1e-10
    level_threshold = 1e-15
    above_pops_thresh = idl_result['pops'].max(axis=0) > level_threshold
    above_ionfrac_thresh = ion.ionization_fraction > ionfrac_threshold
    idl_pops = idl_result['pops'][above_ionfrac_thresh][:, above_pops_thresh]
    fiasco_pops = level_pops[above_ionfrac_thresh][:, above_pops_thresh]
    assert u.allclose(idl_pops, fiasco_pops, rtol=0.01)


def test_idl_compare_radiative_decay_matrices(ion, idl_result):
    assert u.allclose(ion._rate_matrix_radiative_decay,
                      idl_result['rad_decay_rates'],
                      rtol=1e-7)


def test_idl_compare_collisional_electron_rate_matrices(ion, idl_result):
    assert u.allclose(ion._rate_matrix_collisional_electron,
                      idl_result['elec_collision_rates'],
                      rtol=0.01)
