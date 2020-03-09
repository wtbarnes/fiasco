"""
Test element functionality
"""
import numpy as np
import astropy.units as u
import pytest

import fiasco


temperature = np.logspace(5, 8, 100)*u.K


@pytest.fixture
def element(hdf5_dbase_root):
    return fiasco.Element('H', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def another_element(hdf5_dbase_root):
    return fiasco.Element('He', temperature, hdf5_dbase_root=hdf5_dbase_root)


def test_atomic_symbol(element):
    assert element.atomic_symbol == 'H'


def test_atomic_number(element):
    assert element.atomic_number == 1


def test_element_name(element):
    assert element.element_name == 'hydrogen'


def test_add_elements(element, another_element):
    collection = element + another_element
    assert isinstance(collection, fiasco.IonCollection)
    assert collection[0].ion_name == element[0].ion_name


def test_radd_elements(element, another_element):
    collection = another_element + element
    assert isinstance(collection, fiasco.IonCollection)
    assert collection[0].ion_name == another_element[0].ion_name


def test_create_element_number(element, hdf5_dbase_root):
    other_element = fiasco.Element(element.atomic_number,
                                   temperature,
                                   hdf5_dbase_root=hdf5_dbase_root)
    for ion in other_element:
        assert ion in element
    assert other_element.atomic_symbol == element.atomic_symbol
    assert other_element.atomic_number == element.atomic_number
    assert other_element.element_name == element.element_name


def test_create_element_name(element, hdf5_dbase_root):
    other_element = fiasco.Element(element.element_name,
                                   temperature,
                                   hdf5_dbase_root=hdf5_dbase_root)
    for ion in other_element:
        assert ion in element
    assert other_element.atomic_symbol == element.atomic_symbol
    assert other_element.atomic_number == element.atomic_number
    assert other_element.element_name == element.element_name


def test_create_element_lowercase(element, hdf5_dbase_root):
    other_element = fiasco.Element(element.atomic_symbol.lower(),
                                   temperature,
                                   hdf5_dbase_root=hdf5_dbase_root)
    for ion in other_element:
        assert ion in element
    assert other_element.atomic_symbol == element.atomic_symbol
    assert other_element.atomic_number == element.atomic_number
    assert other_element.element_name == element.element_name


def test_getitem_ion_name(element):
    assert element['H 1'] == element[0]
    assert element['H 2'] == element[-1]
    

def test_create_element_without_units_raises_units_error(hdf5_dbase_root):
    with pytest.raises(TypeError):
        fiasco.Element('H', temperature.value, hdf5_dbase_root=hdf5_dbase_root)


def test_create_element_with_wrong_units_raises_unit_conversion_error(hdf5_dbase_root):
    with pytest.raises(u.UnitsError):
        fiasco.Element('H', temperature.value*u.s, hdf5_dbase_root=hdf5_dbase_root)


def test_equilibrium_ionization(element):
    ioneq = element.equilibrium_ionization()
    assert ioneq.shape == element.temperature.shape + (element.atomic_number + 1,)
