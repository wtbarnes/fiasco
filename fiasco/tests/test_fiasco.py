"""
Test package level functions
"""
import numpy as np
import pytest
import astropy.units as u

import fiasco


def test_list_elements():
    # FIXME: actually test the expected elements
    elements = fiasco.list_elements(fiasco.defaults['hdf5_dbase_root'])
    assert type(elements) is list


def test_list_ions():
    # FIXME: actually test the expected ions
    ions = fiasco.list_ions(fiasco.defaults['hdf5_dbase_root'])
    assert type(ions) is list


def test_proton_electron_ratio():
    t = np.logspace(4, 9, 100) * u.K
    pe_ratio = fiasco.proton_electron_ratio(t)
    assert type(pe_ratio) is u.Quantity
    assert pe_ratio.shape == t.shape
