"""
Tests for util functions
"""
import pytest

from fiasco.util import parse_ion_name

@pytest.mark.parametrize('ion_name', [
    "fe 21",
    "Fe 21",
    "Iron 21",
    "iron XXI",
    "Fe xxi",
    "Fe 20+",
    "fe_21",
    "26 21",
    (26, 21),
    ("26", "21"),
    ("iron", "XXI"),
])
def test_parse_ion_name(ion_name):
    element, ion = parse_ion_name(ion_name)
    assert element == 26
    assert ion == 21
