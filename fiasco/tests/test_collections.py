"""
Test collection functionality
"""
import warnings

import numpy as np
import astropy.units as u
import pytest

import fiasco


temperature = np.logspace(5, 8, 100)*u.K


@pytest.fixture
def ion(hdf5_dbase_root):
    return fiasco.Ion('Fe 2', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def another_ion(hdf5_dbase_root):
    return fiasco.Ion('Ca 2', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def element(hdf5_dbase_root):
    return fiasco.Element('hydrogen', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def another_element(hdf5_dbase_root):
    return fiasco.Element('helium', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def collection(hdf5_dbase_root):
    return fiasco.IonCollection(fiasco.Ion('H 1', temperature, hdf5_dbase_root=hdf5_dbase_root),
                                fiasco.Ion('He 2', temperature, hdf5_dbase_root=hdf5_dbase_root))


def test_create_collection_from_ions(ion, another_ion):
    assert isinstance(ion + another_ion, fiasco.IonCollection)
    assert isinstance(fiasco.IonCollection(ion, another_ion), fiasco.IonCollection)


def test_create_collection_from_elements(element, another_element):
    assert isinstance(element + another_element, fiasco.IonCollection)
    assert isinstance(fiasco.IonCollection(element, another_element), fiasco.IonCollection)


def test_create_collection_from_mixture(ion, element, collection):
    assert isinstance(ion + element + collection, fiasco.IonCollection)
    assert isinstance(fiasco.IonCollection(ion, element, collection), fiasco.IonCollection)


def test_create_collection_from_collection(collection):
    assert isinstance(fiasco.IonCollection(collection), fiasco.IonCollection)


def test_getitem(ion, another_ion, element, another_element):
    collection = fiasco.IonCollection(ion, another_ion, element, another_element)
    assert collection[0] == ion
    assert collection[1] == another_ion
    assert isinstance(collection[1:2], fiasco.IonCollection)
    assert collection[1:2][0] == another_ion
    assert isinstance(collection[:2], fiasco.IonCollection)
    for i in collection:
        assert isinstance(i, fiasco.Ion)


def test_contains(collection, hdf5_dbase_root):
    assert 'H 1' in collection
    assert 'hydrogen 1' in collection
    assert 'hydrogen +0' in collection
    ion = fiasco.Ion('H 1', temperature, hdf5_dbase_root=hdf5_dbase_root)
    assert ion in collection


def test_spectrum(hdf5_dbase_root):
    i1 = fiasco.Ion('H 1', 1 * u.MK, hdf5_dbase_root=hdf5_dbase_root)
    i2 = fiasco.Ion('Fe 5', 1 * u.MK, hdf5_dbase_root=hdf5_dbase_root)
    c = i1 + i2
    density = 1e9 * u.cm**-3
    em = 1e29 * u.cm**-5
    with pytest.warns(UserWarning, match='Using 0.1 Angstroms'):
        w, spec = c.spectrum(density, em)
    assert spec.shape == (1, 1, ) + w.shape
    # Add an ion with no spectral information
    i3 = fiasco.Ion('H 2', 1 * u.MK, hdf5_dbase_root=hdf5_dbase_root)
    c += i3
    with pytest.warns(UserWarning, match=f'No transition data available for {i3.ion_name}'):
        w2, spec2 = c.spectrum(density, em)
    assert spec2.shape == (1, 1, ) + w2.shape
    assert np.all(spec == spec2)


def test_spectrum_no_valid_ions(hdf5_dbase_root):
    # Consider the case of an collection with ions with no spectral information
    c2 = fiasco.IonCollection(fiasco.Ion('H 2', 1 * u.MK, hdf5_dbase_root=hdf5_dbase_root))
    with pytest.warns(UserWarning, match='No transition data available for H 2'):
        with pytest.raises(ValueError, match='No collision or transition data available for any ion in collection.'):
            c2.spectrum(1e9 * u.cm**-3, 1e29 * u.cm**-5)


def test_unequal_temperatures_raise_assertion_error(hdf5_dbase_root):
    first_ion = fiasco.Ion('Fe 12', [1e6, 1e7]*u.K, hdf5_dbase_root=hdf5_dbase_root)
    second_ion = fiasco.Ion('Fe 9', [1e4, 1e5]*u.K, hdf5_dbase_root=hdf5_dbase_root)
    with pytest.raises(AssertionError):
        fiasco.IonCollection(first_ion, second_ion)


def test_create_with_wrong_type_raise_type_error(ion, collection):
    with pytest.raises(TypeError):
        fiasco.IonCollection(ion, 0)
    with pytest.raises(TypeError):
        collection + 0
