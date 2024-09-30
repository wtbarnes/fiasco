"""
Test base ion functionality
"""
import pytest

import fiasco

from fiasco.util.exceptions import MissingIonError


@pytest.mark.parametrize('ion_name', [
    "fe 21",
    "Fe 21",
    "Iron 21",
    "iron XXI",
    "Fe xxi",
    "Fe 20+",
    (26, 21),
    ("Fe", "XXI"),
])
def test_create_ion_input_formats(hdf5_dbase_root, ion_name):
    ion = fiasco.base.IonBase(ion_name, hdf5_dbase_root=hdf5_dbase_root)
    assert ion.element_name == 'iron'
    assert ion.atomic_symbol == 'Fe'
    assert ion.atomic_number == 26
    assert ion.ionization_stage == 21
    assert ion.charge_state == 20
    assert ion.isoelectronic_sequence == 'C'
    assert ion._ion_name == 'fe_21'
    assert ion.ion_name_roman == 'Fe XXI'
    assert ion.ionization_stage_roman == 'XXI'


def test_create_invalid_ion_raises_missing_ion_error(hdf5_dbase_root):
    with pytest.raises(MissingIonError):
        fiasco.base.IonBase('hydrogen 3', hdf5_dbase_root=hdf5_dbase_root)
