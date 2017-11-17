"""
Test element functionality
"""
import numpy as np
import astropy.units as u
import pytest

import fiasco


temperature = np.logspace(5, 8, 100)*u.K


@pytest.fixture
def element():
    return fiasco.Element('Fe', temperature)


@pytest.fixture
def another_element():
    return fiasco.Element('Ca', temperature)


def test_atomic_symbol(element):
    assert element.atomic_symbol == 'Fe'


def test_atomic_number(element):
    assert element.atomic_number == 26


def test_element_name(element):
    assert element.element_name == 'iron'


def test_add_elements(element, another_element):
    collection = element + another_element
    assert isinstance(collection, fiasco.IonCollection)
    assert collection[0].ion_name == element[0].ion_name


def test_radd_elements(element, another_element):
    collection = another_element + element
    assert isinstance(collection, fiasco.IonCollection)
    assert collection[0].ion_name == another_element[0].ion_name


def test_create_element_number(element):
    other_element = fiasco.Element(element.atomic_number, temperature)
    assert other_element.ions == element.ions
    assert other_element.atomic_symbol == element.atomic_symbol
    assert other_element.atomic_number == element.atomic_number
    assert other_element.element_name == element.element_name


def test_create_element_name(element):
    other_element = fiasco.Element(element.element_name, temperature)
    assert other_element.ions == element.ions
    assert other_element.atomic_symbol == element.atomic_symbol
    assert other_element.atomic_number == element.atomic_number
    assert other_element.element_name == element.element_name


def test_create_element_lowercase(element):
    other_element = fiasco.Element(element.atomic_symbol.lower(), temperature)
    assert other_element.ions == element.ions
    assert other_element.atomic_symbol == element.atomic_symbol
    assert other_element.atomic_number == element.atomic_number
    assert other_element.element_name == element.element_name


def test_create_element_without_units_raises_units_error():
    with pytest.raises(TypeError):
        fiasco.Element('Fe', temperature.value)


def test_create_element_with_wrong_units_raises_unit_conversion_error():
    with pytest.raises(u.UnitsError):
        fiasco.Element('Fe', temperature.value*u.s)

