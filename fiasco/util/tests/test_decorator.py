"""
Tests for any decorators
"""
import astropy.units as u
import pytest

from fiasco import Ion
from fiasco.util import needs_dataset
from fiasco.util.exceptions import MissingDatasetException


@pytest.fixture()
def ion_like_object(hdf5_dbase_root):
    class NewIon(Ion):
        @needs_dataset('elvlc')
        def elvlc_method(self):
            pass
        @needs_dataset('fake')
        def fake_method(self):
            pass
        @property
        def _fake(self):
            raise KeyError
    return NewIon('Fe XI', 1*u.MK, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture()
def non_ion_like_object():
    class NonIon:
        @needs_dataset('elvlc')
        def elvlc_method(self):
            pass
        @property
        def _elvlc(self):
            pass
        @needs_dataset('fake')
        def fake_method(self):
            pass
        @property
        def _fake(self):
            raise KeyError
        def __str__(self):
            return 'non-ion'
    return NonIon()


def test_decorator_ion_like_has_data(ion_like_object):
    # Just test that this does not throw an exception
    assert ion_like_object.elvlc_method() is None


def test_decorator_ion_like_has_no_data(ion_like_object):
    with pytest.raises(MissingDatasetException, match='fake dataset missing for Fe 11'):
        ion_like_object.fake_method()


def test_decorator_non_ion_like_has_data(non_ion_like_object):
    # Just test that this does not throw an exception
    assert non_ion_like_object.elvlc_method() is None


def test_decorator_non_ion_like_has_no_data(non_ion_like_object):
    with pytest.raises(MissingDatasetException, match='fake dataset missing for non-ion'):
        non_ion_like_object.fake_method()
