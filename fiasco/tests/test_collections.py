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
    return fiasco.Ion('fe_11', temperature)


@pytest.fixture
def element():
    return fiasco.Element('calcium', temperature)


@pytest.fixture
def collection():
    return fiasco.IonCollection(fiasco.Ion('h_1',temperature), fiasco.Ion('he_2',temperature))


def test_create_collection_from_ions(ion):
    another_ion = fiasco.Ion('fe_12', temperature)
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


def test_unequal_temperatures_raise_assertion_error():
    first_ion = fiasco.Ion('fe_12', [1e6, 1e7]*u.K)
    second_ion = fiasco.Ion('fe_9', [1e4, 1e5]*u.K)
    with pytest.raises(AssertionError):
        fiasco.IonCollection(first_ion, second_ion)
