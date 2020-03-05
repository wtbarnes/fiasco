"""
Test ion functionality
"""
import numpy as np
import astropy.units as u
import pytest

import fiasco

temperature = np.logspace(5, 8, 100)*u.K


@pytest.fixture
def ion(hdf5_dbase_root):
    return fiasco.Ion('Fe 5', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def another_ion(hdf5_dbase_root):
    return fiasco.Ion('Fe 6', temperature, hdf5_dbase_root=hdf5_dbase_root)


def test_level_indexing(ion):
    assert isinstance(ion[0], fiasco.ion.Level)


def test_repr(ion):
    assert 'Fe 5' in ion.__repr__()


def test_repr_scalar_temp(ion, hdf5_dbase_root):
    assert 'Fe 5' in fiasco.Ion('Fe 5', 1e6 * u.K, hdf5_dbase_root=hdf5_dbase_root).__repr__()


def test_level_properties(ion):
    assert hasattr(ion[0], 'level')
    assert hasattr(ion[0], 'energy')
    assert hasattr(ion[0], 'configuration')


def test_scalar_temperature(hdf5_dbase_root):
    ion = fiasco.Ion('H 1', 1 * u.MK, hdf5_dbase_root=hdf5_dbase_root)
    ioneq = ion.ioneq
    assert ioneq.shape == (1,)


def test_scalar_density(hdf5_dbase_root):
    ion = fiasco.Ion('H 1', temperature, hdf5_dbase_root=hdf5_dbase_root)
    pop = ion.level_populations(1e8 * u.cm**-3)
    assert pop.shape == ion.temperature.shape + (1,) + ion._elvlc['level'].shape


def test_no_elvlc_raises_index_error(hdf5_dbase_root):
    with pytest.raises(IndexError):
        fiasco.Ion('H 2', temperature, hdf5_dbase_root=hdf5_dbase_root)[0]


def test_ioneq(ion):
    assert ion.ioneq.shape == temperature.shape


def test_abundance(ion):
    assert ion.abundance.dtype == np.dtype('float64')


def test_missing_abundance(hdf5_dbase_root):
    ion = fiasco.Ion('Li 1',
                     temperature,
                     abundance_filename='sun_coronal_1992_feldman',
                     hdf5_dbase_root=hdf5_dbase_root)
    assert ion.abundance is None


def test_ip(ion):
    assert ion.ip.dtype == np.dtype('float64')


def test_missing_ip(hdf5_dbase_root):
    ion = fiasco.Ion('Fe 27', temperature, hdf5_dbase_root=hdf5_dbase_root)
    assert ion.ip is None


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
