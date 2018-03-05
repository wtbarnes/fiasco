"""
Test collection functionality
"""
import numpy as np
import astropy.units as u
import pytest

import fiasco


temperature = np.logspace(5, 8, 100)*u.K


@pytest.fixture
def ion():
    return fiasco.Ion('Fe 11', temperature)


@pytest.fixture
def another_ion():
    return fiasco.Ion('Ca 1', temperature)


@pytest.fixture
def element():
    return fiasco.Element('calcium', temperature)


@pytest.fixture
def collection():
    return fiasco.IonCollection(fiasco.Ion('H 1', temperature), fiasco.Ion('He 2', temperature))


def test_create_collection_from_ions(ion):
    another_ion = fiasco.Ion('Fe 12', temperature)
    assert isinstance(ion + another_ion, fiasco.IonCollection)
    assert isinstance(fiasco.IonCollection(ion, another_ion), fiasco.IonCollection)


def test_create_collection_from_elements(element):
    another_element = fiasco.Element('iron', temperature)
    assert isinstance(element + another_element, fiasco.IonCollection)
    assert isinstance(fiasco.IonCollection(element, another_element), fiasco.IonCollection)


def test_create_collection_from_mixture(ion, element, collection):
    assert isinstance(ion + element + collection, fiasco.IonCollection)
    assert isinstance(fiasco.IonCollection(ion, element, collection), fiasco.IonCollection)


def test_create_collection_from_collection(collection):
    assert isinstance(fiasco.IonCollection(collection), fiasco.IonCollection)


def test_getitem(ion, another_ion):
    collection = fiasco.IonCollection(ion, another_ion)
    assert collection[0] == ion
    assert collection[1] == another_ion


def test_contains(collection):
    assert 'H 1' in collection
    assert 'hydrogen 1' in collection
    assert 'hydrogen +0' in collection
    ion = fiasco.Ion('H 1', temperature)
    assert ion in collection


def test_unequal_temperatures_raise_assertion_error():
    first_ion = fiasco.Ion('Fe 12', [1e6, 1e7]*u.K)
    second_ion = fiasco.Ion('Fe 9', [1e4, 1e5]*u.K)
    with pytest.raises(AssertionError):
        fiasco.IonCollection(first_ion, second_ion)


def test_create_with_wrong_type_raise_type_error(ion, collection):
    with pytest.raises(TypeError):
        fiasco.IonCollection(ion, 0)
    with pytest.raises(TypeError):
        collection + 0
