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


def test_ioneq(ion):
    assert ion.ioneq.shape == temperature.shape


def test_abundance(ion):
    assert ion.abundance.dtype == np.dtype('float64')


def test_ip(ion):
    assert ion.ip.dtype == np.dtype('float64')


def test_create_ion_without_units_raises_units_error():
    with pytest.raises(TypeError):
        fiasco.Ion('fe_5', temperature.value)


def test_create_ion_with_wrong_units_raises_unit_conversion_error():
    with pytest.raises(u.UnitsError):
        fiasco.Ion('fe_5', temperature.value*u.s)
