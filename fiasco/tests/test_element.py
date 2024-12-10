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
    ionization_fraction = carbon.equilibrium_ionization
    assert ionization_fraction.shape == carbon.temperature.shape + (carbon.atomic_number + 1,)
    assert u.allclose(ionization_fraction[33, 5], u.Quantity(0.5787345345914312), atol=0.0, rtol=1e-10)


def test_element_repr(element):
    assert element.element_name in element.__repr__()


@pytest.mark.parametrize(('value', 'dset'),[
    (0.07943282347242822, 'sun_coronal_1992_feldman_ext'),
    (0.08511380382023759, 'sun_photospheric_2007_grevesse'),
    (1e-3, None),
])
def test_change_element_abundance(another_element, value, dset):
    another_element.abundance = value if dset is None else dset
    assert u.allclose(another_element.abundance, value)
    assert u.allclose(another_element[1].abundance, value)
