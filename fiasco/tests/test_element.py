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


def test_atomic_symbol(element):
    assert element.atomic_symbol == 'Fe'


def test_atomic_number(element):
    assert element.atomic_number == 26


def test_element_name(element):
    assert element.element_name == 'iron'


def test_create_element_without_units_raises_units_error():
    with pytest.raises(TypeError):
        fiasco.Element('Fe', temperature.value)


def test_create_element_with_wrong_units_raises_unit_conversion_error():
    with pytest.raises(u.UnitsError):
        fiasco.Element('Fe', temperature.value*u.s)

