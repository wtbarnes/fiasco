"""
Test package level functions
"""
import astropy.units as u
import numpy as np

import fiasco


def test_list_elements(hdf5_dbase_root):
    elements = fiasco.list_elements(hdf5_dbase_root)
    assert elements[:5] == ['H', 'He', 'Li', 'Be', 'B']


def test_list_ions(hdf5_dbase_root):
    ions = fiasco.list_ions(hdf5_dbase_root)
    assert ions[:5] == ['H 1', 'H 2', 'He 1', 'He 2', 'He 3']


def test_get_isoelectronic_sequence(hdf5_dbase_root):
    iso_seq = fiasco.get_isoelectronic_sequence('iron',
                                                hdf5_dbase_root=hdf5_dbase_root)
    assert iso_seq == ['Fe 1',
                       'Co 2',
                       'Ni 3',
                       'Cu 4',
                       'Zn 5',]


def test_proton_electron_ratio(hdf5_dbase_root):
    t = np.logspace(4, 9, 100) * u.K
    # NOTE: this number will not be accurate as we are using only a subset of
    # the database
    pe_ratio = fiasco.proton_electron_ratio(t, hdf5_dbase_root=hdf5_dbase_root)
    assert type(pe_ratio) is u.Quantity
    assert pe_ratio.shape == t.shape
    assert np.isfinite(pe_ratio).all()
