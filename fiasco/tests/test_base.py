"""
Test base ion functionality
"""
import pytest

import fiasco
from fiasco.util.exceptions import MissingIonError


@pytest.fixture
def ionbase(hdf5_dbase_root):
    return fiasco.base.IonBase('fe 5', hdf5_dbase_root=hdf5_dbase_root)


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


def test_create_ion_symbol_lower(hdf5_dbase_root):
    ion = fiasco.base.IonBase('h 1', hdf5_dbase_root=hdf5_dbase_root)
    assert ion.element_name == 'hydrogen'
    assert ion.atomic_symbol == 'H'
    assert ion.atomic_number == 1
    assert ion.ionization_stage == 1
    assert ion.charge_state == 0
    assert ion._ion_name == 'h_1'


def test_create_ion_symbol_upper(hdf5_dbase_root):
    ion = fiasco.base.IonBase('H 1', hdf5_dbase_root=hdf5_dbase_root)
    assert ion.element_name == 'hydrogen'
    assert ion.atomic_symbol == 'H'
    assert ion.atomic_number == 1
    assert ion.ionization_stage == 1
    assert ion.charge_state == 0
    assert ion._ion_name == 'h_1'


def test_create_ion_charge_state(hdf5_dbase_root):
    ion = fiasco.base.IonBase('h +0', hdf5_dbase_root=hdf5_dbase_root)
    assert ion.element_name == 'hydrogen'
    assert ion.atomic_symbol == 'H'
    assert ion.atomic_number == 1
    assert ion.ionization_stage == 1
    assert ion.charge_state == 0
    assert ion._ion_name == 'h_1'


def test_create_ion_element_name(hdf5_dbase_root):
    ion = fiasco.base.IonBase('hydrogen 1', hdf5_dbase_root=hdf5_dbase_root)
    assert ion.element_name == 'hydrogen'
    assert ion.atomic_symbol == 'H'
    assert ion.atomic_number == 1
    assert ion.ionization_stage == 1
    assert ion.charge_state == 0
    assert ion._ion_name == 'h_1'


def test_create_invalid_ion_raises_missing_ion_error(hdf5_dbase_root):
    with pytest.raises(MissingIonError):
        fiasco.base.IonBase('hydrogen 3', hdf5_dbase_root=hdf5_dbase_root)
