"""
Test base ion functionality
"""
import pytest

import fiasco


@pytest.fixture
def ionbase():
    return fiasco.IonBase('fe 5')


def test_ion_name(ionbase):
    assert ionbase.ion_name == 'Fe 5'
    assert ionbase._ion_name == 'fe_5'


def test_atomic_symbol(ionbase):
    assert ionbase.atomic_symbol == 'Fe'


def test_atomic_number(ionbase):
    assert ionbase.atomic_number == 26


def test_element_name(ionbase):
    assert ionbase.element_name == 'iron'


def test_ionization_stage(ionbase):
    assert ionbase.ionization_stage == 5


def test_charge_state(ionbase):
    assert ionbase.charge_state == 4


def test_create_ion_symbol_lower():
    ion = fiasco.IonBase('h 1')
    assert ion.element_name == 'hydrogen'
    assert ion.atomic_symbol == 'H'
    assert ion.atomic_number == 1
    assert ion.ionization_stage == 1
    assert ion.charge_state == 0
    assert ion._ion_name == 'h_1'


def test_create_ion_symbol_upper():
    ion = fiasco.IonBase('H 1')
    assert ion.element_name == 'hydrogen'
    assert ion.atomic_symbol == 'H'
    assert ion.atomic_number == 1
    assert ion.ionization_stage == 1
    assert ion.charge_state == 0
    assert ion._ion_name == 'h_1'


def test_create_ion_charge_state():
    ion = fiasco.IonBase('h +0')
    assert ion.element_name == 'hydrogen'
    assert ion.atomic_symbol == 'H'
    assert ion.atomic_number == 1
    assert ion.ionization_stage == 1
    assert ion.charge_state == 0
    assert ion._ion_name == 'h_1'


def test_create_ion_element_name():
    ion = fiasco.IonBase('hydrogen 1')
    assert ion.element_name == 'hydrogen'
    assert ion.atomic_symbol == 'H'
    assert ion.atomic_number == 1
    assert ion.ionization_stage == 1
    assert ion.charge_state == 0
    assert ion._ion_name == 'h_1'
