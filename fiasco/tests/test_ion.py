"""
Test ion functionality
"""
import numpy as np
import astropy.units as u
import pytest

import fiasco


temperature = np.logspace(5, 8, 100)*u.K


@pytest.fixture
def ion():
    return fiasco.Ion('fe_5', temperature)


@pytest.fixture
def another_ion():
    return fiasco.Ion('fe_6', temperature)


def test_ioneq(ion):
    assert ion.ioneq.shape == temperature.shape


def test_abundance(ion):
    assert ion.abundance.dtype == np.dtype('float64')


def test_missing_abundance():
    ion = fiasco.Ion('li_1', temperature, abundance_filename='sun_coronal_1992_feldman')
    assert ion.abundance is None


def test_ip(ion):
    assert ion.ip.dtype == np.dtype('float64')

def test_missing_ip():
    ion = fiasco.Ion('fe_27', temperature)
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


def test_create_ion_without_units_raises_units_error():
    with pytest.raises(TypeError):
        fiasco.Ion('fe_5', temperature.value)


def test_create_ion_with_wrong_units_raises_unit_conversion_error():
    with pytest.raises(u.UnitsError):
        fiasco.Ion('fe_5', temperature.value*u.s)
