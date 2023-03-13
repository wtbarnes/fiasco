"""
Test ion functionality
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

from fiasco.util.exceptions import MissingDatasetException

temperature = np.logspace(5, 8, 100)*u.K


@pytest.fixture
def ion(hdf5_dbase_root):
    return fiasco.Ion('Fe 5', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def another_ion(hdf5_dbase_root):
    return fiasco.Ion('Fe 6', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def fe10(hdf5_dbase_root):
    return fiasco.Ion('Fe 10', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def c6(hdf5_dbase_root):
    return fiasco.Ion('C VI', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def fe20(hdf5_dbase_root):
    # NOTE: This ion was added because it has reclvl and cilvl files which
    # we need to test the level-resolved rate correction factor
    return fiasco.Ion('Fe XX', temperature, hdf5_dbase_root=hdf5_dbase_root)


def test_new_instance(ion):
    abundance_filename = ion._instance_kwargs['abundance_filename']
    new_ion = ion._new_instance()
    for k in new_ion._instance_kwargs:
        assert new_ion._instance_kwargs[k] == ion._instance_kwargs[k]
    assert u.allclose(new_ion.temperature, ion.temperature, rtol=0)
    new_ion = ion._new_instance(temperature=ion.temperature[:1])
    assert u.allclose(new_ion.temperature, ion.temperature[:1])
    new_ion = ion._new_instance(abundance_filename='sun_coronal_1992_feldman')
    assert new_ion._instance_kwargs['abundance_filename'] == 'sun_coronal_1992_feldman'
    assert ion._instance_kwargs['abundance_filename'] == abundance_filename


def test_level_indexing(ion):
    # Integer
    assert isinstance(ion[0], fiasco.Level)
    assert ion[0].__repr__() == fiasco.Level(0, ion._elvlc).__repr__()
    # Slice
    levels = ion[:5]
    assert len(levels) == 5
    assert isinstance(levels, list)
    assert isinstance(levels[0], fiasco.Level)
    assert levels[2].__repr__() == fiasco.Level(2, ion._elvlc).__repr__()
    # Fancy indexing
    levels = ion[[1, 5, 10]]
    assert len(levels) == 3
    assert isinstance(levels, list)
    assert isinstance(levels[0], fiasco.Level)
    assert levels[2].__repr__() == fiasco.Level(10, ion._elvlc).__repr__()


def test_level(ion):
    level = ion[0]
    assert isinstance(level, fiasco.Level)
    assert level.multiplicity == 5
    assert level.total_angular_momentum == 0
    assert level.orbital_angular_momentum_label == 'D'


def test_repr(ion):
    assert 'Fe 5' in ion.__repr__()


def test_repr_scalar_temp(hdf5_dbase_root):
    assert 'Fe 5' in fiasco.Ion('Fe 5', 1e6 * u.K, hdf5_dbase_root=hdf5_dbase_root).__repr__()


def test_ion_properties(ion):
    assert ion.atomic_number == 26
    assert ion.element_name == 'iron'
    assert ion.atomic_symbol == 'Fe'
    assert ion.ion_name == 'Fe 5'


def test_level_properties(ion):
    assert hasattr(ion[0], 'level')
    assert hasattr(ion[0], 'energy')
    assert hasattr(ion[0], 'configuration')


def test_scalar_temperature(hdf5_dbase_root):
    ion = fiasco.Ion('H 1', 1 * u.MK, hdf5_dbase_root=hdf5_dbase_root)
    ioneq = ion.ioneq
    assert ioneq.shape == (1,)
    t_data = ion._ioneq[ion._dset_names['ioneq_filename']]['temperature']
    ioneq_data = ion._ioneq[ion._dset_names['ioneq_filename']]['ionization_fraction']
    i_t = np.where(t_data == ion.temperature)
    assert u.allclose(ioneq, ioneq_data[i_t])


def test_no_elvlc_raises_index_error(hdf5_dbase_root):
    with pytest.raises(IndexError):
        fiasco.Ion('H 2', temperature, hdf5_dbase_root=hdf5_dbase_root)[0]


def test_ioneq(ion):
    t_data = ion._ioneq[ion._dset_names['ioneq_filename']]['temperature']
    ioneq_data = ion._ioneq[ion._dset_names['ioneq_filename']]['ionization_fraction']
    ion_at_nodes = ion._new_instance(temperature=t_data)
    assert u.allclose(ion_at_nodes.ioneq, ioneq_data, rtol=1e-6)


def test_ioneq_positive(ion):
    assert np.all(ion.ioneq >= 0)


def test_ioneq_out_bounds_is_nan(ion):
    t_data = ion._ioneq[ion._dset_names['ioneq_filename']]['temperature']
    t_out_of_bounds = t_data[[0,-1]] + [-100, 1e6] * u.K
    ion_out_of_bounds = ion._new_instance(temperature=t_out_of_bounds)
    assert np.isnan(ion_out_of_bounds.ioneq).all()


def test_formation_temeprature(ion):
    assert ion.formation_temperature == ion.temperature[np.argmax(ion.ioneq)]


def test_abundance(ion):
    assert ion.abundance.dtype == np.dtype('float64')
    # This value has not been tested for correctness
    assert u.allclose(ion.abundance, 0.0001258925411794166)


def test_proton_collision(fe10):
    rate = fe10.proton_collision_excitation_rate
    assert u.allclose(rate[0, 0], 4.69587161e-13 * u.cm**3 / u.s)

    rate = fe10.proton_collision_deexcitation_rate
    assert u.allclose(rate[0, 0], 1.17688025e-12 * u.cm**3 / u.s)


def test_missing_abundance(hdf5_dbase_root):
    ion = fiasco.Ion('Li 1',
                     temperature,
                     abundance_filename='sun_coronal_1992_feldman',
                     hdf5_dbase_root=hdf5_dbase_root)
    with pytest.raises(KeyError):
        _ = ion.abundance


def test_ip(ion):
    assert ion.ip.dtype == np.dtype('float64')
    # This value has not been tested for correctness
    assert u.allclose(ion.ip, 1.2017997435751017e-10 * u.erg)


def test_missing_ip(hdf5_dbase_root):
    ion = fiasco.Ion('Fe 27', temperature, hdf5_dbase_root=hdf5_dbase_root)
    with pytest.raises(MissingDatasetException):
        _ = ion.ip


def test_level_populations(ion):
    pop = ion.level_populations(1e8 * u.cm**-3)
    assert pop.shape == ion.temperature.shape + (1,) + ion._elvlc['level'].shape
    # This value has not been checked for correctness
    assert u.allclose(pop[0, 0, 0], 0.011643747849652244)
    # Check that the total populations are normalized to 1 for all temperatures
    assert u.allclose(pop.squeeze().sum(axis=1), 1, atol=None, rtol=1e-15)


def test_level_populations_proton_data_toggle(ion):
    # Fe V has no psplups data so the toggle should have no effect
    lp_protons = ion.level_populations(1e9*u.cm**(-3), include_protons=True)
    lp_no_protons = ion.level_populations(1e9*u.cm**(-3), include_protons=False)
    assert u.allclose(lp_protons, lp_no_protons, atol=0, rtol=0)


def test_contribution_function(ion):
    cont_func = ion.contribution_function(1e7 * u.cm**-3)
    assert cont_func.shape == ion.temperature.shape + (1, ) + ion._wgfa['wavelength'].shape
    # This value has not been tested for correctness
    assert u.allclose(cont_func[0, 0, 0], 2.08668713e-30 * u.cm**3 * u.erg / u.s)


def test_emissivity_shape(c6):
    # NOTE: Explicitly testing C VI here because it has a psplups file
    # and thus will compute the proton rates as well which have a specific
    # codepath for coupled density/temperature.
    # NOTE: Test that coupled temperature/density appropriately propagate through
    # and that resulting quantity has the right shape.
    # Using the emissivity quantity here because it is the highest level
    # product that needs to manipulate the density. This will implicitly test the
    # contribution function as well.
    #
    # Scalar, no coupling
    density = 1e9 * u.cm**(-3)
    emiss = c6.emissivity(density)
    wavelength = c6.transitions.wavelength[~c6.transitions.is_twophoton]
    assert emiss.shape == c6.temperature.shape + (1,) + wavelength.shape
    # Array, no coupling
    density = [1e8, 1e9, 1e10] * u.cm**(-3)
    emiss = c6.emissivity(density)
    wavelength = c6.transitions.wavelength[~c6.transitions.is_twophoton]
    assert emiss.shape == c6.temperature.shape + density.shape + wavelength.shape
    # Array, with coupling
    pressure = 1e15 * u.K * u.cm**(-3)
    density = pressure / c6.temperature
    emiss = c6.emissivity(density, couple_density_to_temperature=True)
    wavelength = c6.transitions.wavelength[~c6.transitions.is_twophoton]
    assert emiss.shape == c6.temperature.shape + (1,) + wavelength.shape


def test_coupling_unequal_dimensions_exception(ion):
    with pytest.raises(ValueError, match='Temperature and density must be of equal length'):
        _ = ion.level_populations([1e7, 1e8]*u.cm**(-3), couple_density_to_temperature=True)


@pytest.fixture
def pops_with_correction(fe20):
    return fe20.level_populations(1e9*u.cm**(-3)).squeeze()


@pytest.fixture
def pops_no_correction(fe20):
    return fe20.level_populations(1e9*u.cm**(-3),
                                  include_level_resolved_rate_correction=False).squeeze()


def test_level_populations_normalized(pops_no_correction, pops_with_correction):
    assert u.allclose(pops_with_correction.sum(axis=1), 1, atol=None, rtol=1e-15)
    assert u.allclose(pops_no_correction.sum(axis=1), 1, atol=None, rtol=1e-15)


def test_level_populations_correction(fe20, pops_no_correction, pops_with_correction):
    # Test level-resolved correction applied to correct levels
    i_corrected = np.unique(np.concatenate([fe20._cilvl['upper_level'], fe20._reclvl['upper_level']]))
    i_corrected -= 1
    # This tests that, for at least some portion of the temperature axis, the populations are
    # significantly different for each corrected level
    pops_equal = u.isclose(pops_with_correction[:, i_corrected], pops_no_correction[:, i_corrected],
                           atol=0.0, rtol=1e-5)
    assert ~np.all(np.all(pops_equal, axis=0))
    # All other levels should be unchanged (with some tolerance for renormalization)
    is_uncorrected = np.ones(pops_no_correction.shape[-1], dtype=bool)
    is_uncorrected[i_corrected] = False
    i_uncorrected = np.where(is_uncorrected)
    assert u.allclose(pops_with_correction[:, i_uncorrected], pops_no_correction[:, i_uncorrected],
                      atol=0.0, rtol=1e-5)


def test_emissivity(ion):
    emm = ion.emissivity(1e7 * u.cm**-3)
    assert emm.shape == ion.temperature.shape + (1, ) + ion._wgfa['wavelength'].shape
    # This value has not been tested for correctness
    assert u.allclose(emm[0, 0, 0], 2.08668713e-16 * u.erg / u.cm**3 / u.s)


@pytest.mark.parametrize('em', [
    1e29 * u.cm**-5,
    [1e29] * u.cm**-5,
    1e29 * np.ones(temperature.shape) * u.cm**-5,
])
def test_intensity(ion, em):
    wave_shape = ion._wgfa['wavelength'].shape
    intens = ion.intensity(1e7 * u.cm**-3, em)
    assert intens.shape == ion.temperature.shape + (1, ) + wave_shape
    # Test density varying along independent axis
    density = [1e7, 1e9, 1e10] * u.cm**(-3)
    intens = ion.intensity(density, em)
    assert intens.shape == ion.temperature.shape + density.shape + wave_shape
    # Test density varying along same axis as temperature
    density = 1e15 * u.K * u.cm**(-3) / ion.temperature
    intens = ion.intensity(density, em, couple_density_to_temperature=True)
    assert intens.shape == ion.temperature.shape + (1, ) + wave_shape


def test_excitation_autoionization_rate(ion):
    rate = ion.excitation_autoionization_rate
    assert rate.shape == ion.temperature.shape
    # This value has not been tested for correctness
    assert u.allclose(rate[0], 1.14821255e-12 * u.cm**3 / u.s)


def test_dielectronic_recombination_rate(ion):
    rate = ion.dielectronic_recombination_rate
    assert rate.shape == ion.temperature.shape
    # This value has not been tested for correctness
    assert u.allclose(rate[0], 1.60593802e-11 * u.cm**3 / u.s)


def test_free_free(ion):
    emission = ion.free_free(200 * u.Angstrom)
    assert emission.shape == ion.temperature.shape + (1, )
    # This value has not been tested for correctness
    assert u.allclose(emission[0], 1.72804216e-29 * u.cm**3 * u.erg / u.Angstrom / u.s)


def test_free_bound(ion):
    emission = ion.free_bound(200 * u.Angstrom)
    assert emission.shape == ion.temperature.shape + (1, )
    # This value has not been tested for correctness
    assert u.allclose(emission[0, 0], 9.7902609e-26 * u.cm**3 * u.erg / u.Angstrom / u.s)


def test_add_ions(ion, another_ion):
    collection = ion + another_ion
    assert isinstance(collection, fiasco.IonCollection)
    assert collection[0] == ion
    assert collection[1] == another_ion


def test_radd_ions(ion, another_ion):
    collection = another_ion + ion
    assert isinstance(collection, fiasco.IonCollection)
    assert collection[1] == ion
    assert collection[0] == another_ion


def test_transitions(ion):
    trans = ion.transitions
    assert isinstance(trans, fiasco.Transitions)
    assert len(trans) == 361
    # These values have not been tested for correctness
    assert not trans.is_twophoton[0]
    assert trans.is_observed[0]
    assert u.allclose(trans.A[0], 0.000155 / u.s)
    assert u.allclose(trans.wavelength[0], 703729.75 * u.Angstrom)
    assert u.allclose(trans.upper_level[0], 2)
    assert u.allclose(trans.lower_level[0], 1)
    assert u.allclose(trans.delta_energy[0], 2.82273956e-14 * u.erg)


def test_create_ion_without_units_raises_units_error(hdf5_dbase_root):
    with pytest.raises(TypeError):
        fiasco.Ion('Fe 5', temperature.value, hdf5_dbase_root=hdf5_dbase_root)


def test_create_ion_with_wrong_units_raises_unit_conversion_error(hdf5_dbase_root):
    with pytest.raises(u.UnitsError):
        fiasco.Ion('Fe 5', temperature.value*u.s, hdf5_dbase_root=hdf5_dbase_root)


def test_indexing_no_levels(hdf5_dbase_root):
    ion = fiasco.Ion('Fe 1', temperature, hdf5_dbase_root=hdf5_dbase_root)
    print(ion)
    assert [l for l in ion] == []
    with pytest.raises(IndexError, match='No energy levels available for Fe 1'):
        ion[0]


def test_repr_no_levels(hdf5_dbase_root):
    """
    Ensures the repr can be printed without errors even when
    no energy level or transition information is available.
    """
    assert fiasco.Ion('Fe 1', temperature, hdf5_dbase_root=hdf5_dbase_root).__repr__


def test_next_ion(ion):
    next_ion = ion.next_ion()
    assert next_ion.ionization_stage == ion.ionization_stage + 1
    assert next_ion.atomic_number == ion.atomic_number


def test_previous_ion(ion):
    prev_ion = ion.previous_ion()
    assert prev_ion.ionization_stage == ion.ionization_stage - 1
    assert prev_ion.atomic_number == ion.atomic_number
