"""
Test ion functionality
"""
import numpy as np
import astropy.units as u
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


def test_level_indexing(ion):
    assert isinstance(ion[0], fiasco.ion.Level)


def test_repr(ion):
    assert 'Fe 5' in ion.__repr__()


def test_repr_scalar_temp(ion, hdf5_dbase_root):
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
    np.testing.assert_allclose(ioneq, ioneq_data[i_t])


def test_scalar_density(hdf5_dbase_root):
    ion = fiasco.Ion('H 1', temperature, hdf5_dbase_root=hdf5_dbase_root)
    pop = ion.level_populations(1e8 * u.cm**-3)
    assert pop.shape == ion.temperature.shape + (1,) + ion._elvlc['level'].shape
    # This value has not been checked for correctness
    np.testing.assert_allclose(pop[0, 0, 0], 0.9965048292729177)


def test_no_elvlc_raises_index_error(hdf5_dbase_root):
    with pytest.raises(IndexError):
        fiasco.Ion('H 2', temperature, hdf5_dbase_root=hdf5_dbase_root)[0]


def test_ioneq(ion):
    assert ion.ioneq.shape == temperature.shape
    t_data = ion._ioneq[ion._dset_names['ioneq_filename']]['temperature']
    ioneq_data = ion._ioneq[ion._dset_names['ioneq_filename']]['ionization_fraction']
    i_t = np.where(t_data == ion.temperature[0])
    # Essentially test that we've done the interpolation to the data correctly
    # for a single value
    np.testing.assert_allclose(ion.ioneq[0], ioneq_data[i_t])


def test_abundance(ion):
    assert ion.abundance.dtype == np.dtype('float64')
    # This value has not been tested for correctness
    np.testing.assert_allclose(ion.abundance, 3.1622776601683795e-05)


def test_proton_collision(fe10):
    rate = fe10.proton_collision_excitation_rate()
    assert u.allclose(rate[0, 0], 4.69587161e-13 * u.cm**3 / u.s)

    rate = fe10.proton_collision_deexcitation_rate()
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


def test_contribution_function(ion):
    cont_func = ion.contribution_function(1e7 * u.cm**-3)
    assert cont_func.shape == ion.temperature.shape + (1, ) + ion._wgfa['wavelength'].shape
    # This value has not been tested for correctness
    assert u.allclose(cont_func[0, 0, 0], 5.24152109e-31 * u.cm**3 * u.erg / u.s)


def test_emissivity(ion):
    emm = ion.emissivity(1e7 * u.cm**-3)
    assert emm.shape == ion.temperature.shape + (1, ) + ion._wgfa['wavelength'].shape
    # This value has not been tested for correctness
    assert u.allclose(emm[0, 0, 0], 5.24152109e-17 * u.erg / u.cm**3 / u.s)


def test_excitation_autoionization_rate(ion):
    rate = ion.excitation_autoionization_rate()
    assert rate.shape == ion.temperature.shape
    # This value has not been tested for correctness
    assert u.allclose(rate[0], 1.14821255e-12 * u.cm**3 / u.s)


def test_dielectronic_recombination_rate(ion):
    rate = ion.dielectronic_recombination_rate()
    assert rate.shape == ion.temperature.shape
    # This value has not been tested for correctness
    assert u.allclose(rate[0], 1.60593802e-11 * u.cm**3 / u.s)


def test_free_free(ion):
    emission = ion.free_free(200 * u.Angstrom)
    assert emission.shape == ion.temperature.shape + (1, )
    # This value has not been tested for correctness
    assert u.allclose(emission[0], 6.81123745e-28 * u.cm**3 * u.erg / u.Angstrom / u.s)


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


def test_create_ion_without_units_raises_units_error(hdf5_dbase_root):
    with pytest.raises(TypeError):
        fiasco.Ion('Fe 5', temperature.value, hdf5_dbase_root=hdf5_dbase_root)


def test_create_ion_with_wrong_units_raises_unit_conversion_error(hdf5_dbase_root):
    with pytest.raises(u.UnitsError):
        fiasco.Ion('Fe 5', temperature.value*u.s, hdf5_dbase_root=hdf5_dbase_root)
