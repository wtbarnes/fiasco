"""
Tests for util functions
"""
import pytest

from fiasco.util import get_chianti_catalog, parse_ion_name


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


def test_get_chianti_catalog(ascii_dbase_root):
    catalog = get_chianti_catalog(ascii_dbase_root)
    keys = ['abundance_files',
            'ioneq_files',
            'ip_files',
            'continuum_files',
            'ion_files']
    for k in catalog:
        assert k in catalog
        assert isinstance(catalog[k], list)
