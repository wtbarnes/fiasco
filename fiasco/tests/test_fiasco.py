"""
Test package level functions
"""
import numpy as np
import astropy.units as u
import pytest

import fiasco


def test_list_elements(hdf5_dbase_root):
    # FIXME: actually test the expected elements
    elements = fiasco.list_elements(hdf5_dbase_root)
    assert type(elements) is list


def test_list_ions(hdf5_dbase_root):
    # FIXME: actually test the expected ions
    ions = fiasco.list_ions(hdf5_dbase_root)
    assert type(ions) is list


def test_proton_electron_ratio(hdf5_dbase_root):
    t = np.logspace(4, 9, 100) * u.K
    # NOTE: this number will not be accurate as we are using only a subset of
    # the database
    pe_ratio = fiasco.proton_electron_ratio(t, hdf5_dbase_root=hdf5_dbase_root)
    assert type(pe_ratio) is u.Quantity
    assert pe_ratio.shape == t.shape


def test_proton_electron_ratio_missing_abundance_warning(hdf5_dbase_root):
    t = np.logspace(4, 9, 100) * u.K
    # NOTE: use the Feldman abundance set as we know it does not contain abundances
    # for all of the elements in CHIANTI. For the missing elements, a UserWarning
    # should be raised.
    with pytest.warns(UserWarning, match='Abundance not available'):
        _ = fiasco.proton_electron_ratio(t, hdf5_dbase_root=hdf5_dbase_root,
                                         abundance_filename='sun_coronal_1992_feldman')
