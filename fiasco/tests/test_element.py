"""
Test element functionality
"""
import astropy.units as u
import numpy as np
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


@pytest.mark.parametrize('symbol', [1, 'hydrogen', 'Hydrogen', 'H'])
def test_create_element_number(symbol, hdf5_dbase_root):
    other_element = fiasco.Element(symbol, temperature, hdf5_dbase_root=hdf5_dbase_root)
    assert other_element.atomic_symbol == 'H'
    assert other_element.atomic_number == 1
    assert other_element.element_name == 'hydrogen'


def test_getitem_ion_name(element):
    assert element['H 1'] == element[0]
    assert element['H 1'].ionization_stage == 1
    assert element['H 1'].atomic_number == 1
    assert element['H 2'] == element[1]
    assert element['H 2'].ionization_stage == 2
    assert element['H 2'].atomic_number == 1
    assert element['H +1'] == element[1]


def test_create_element_without_units_raises_units_error(hdf5_dbase_root):
    with pytest.raises(TypeError):
        fiasco.Element('H', temperature.value, hdf5_dbase_root=hdf5_dbase_root)


def test_create_element_with_wrong_units_raises_unit_conversion_error(hdf5_dbase_root):
    with pytest.raises(u.UnitsError):
        fiasco.Element('H', temperature.value*u.s, hdf5_dbase_root=hdf5_dbase_root)


def test_equilibrium_ionization(hdf5_dbase_root):
    # NOTE: Using an element with Z>1 so that we are testing the full rate matrix
    # computation. Using C here because we are already using this for the gallery
    carbon = fiasco.Element('C', temperature, hdf5_dbase_root=hdf5_dbase_root)
    ioneq = carbon.equilibrium_ionization
    assert ioneq.shape == carbon.temperature.shape + (carbon.atomic_number + 1,)
    assert ioneq[33, 5] == u.Quantity(0.5787345345914312)


def test_element_repr(element):
    assert element.element_name in element.__repr__()
