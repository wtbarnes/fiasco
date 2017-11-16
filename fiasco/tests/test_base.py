"""
Test base ion functionality
"""
import pytest

import fiasco


@pytest.fixture
def ionbase():
    return fiasco.IonBase('fe_5')


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
