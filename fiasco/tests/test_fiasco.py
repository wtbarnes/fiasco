"""
Test package level functions
"""
import pytest

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


@pytest.fixture
def fe_13(hdf5_dbase_root):
    return fiasco.Ion('Fe XIII', 10**np.arange(5.5, 6.5, 0.05)*u.K, hdf5_dbase_root=hdf5_dbase_root)


@pytest.mark.parametrize(('numerator', 'denominator', 'result'),[
    (203.795*u.AA, 202.044*u.AA, [0.03940454, 0.30392223, 0.95424066]),
    ('3s2 3p 3d 3D2 -- 3s2 3p2 3P2', '3s2 3p 3d 3P1 -- 3s2 3p2 3P0', [0.03940454, 0.30392223, 0.95424066]),
    ([203.795, 203.826]*u.AA, 202.044*u.AA, [0.12535491, 1.05989413, 3.60720143]),
    (['3s2 3p 3d 3D2 -- 3s2 3p2 3P2', '3s2 3p 3d 3D3 -- 3s2 3p2 3P2'], '3s2 3p 3d 3P1 -- 3s2 3p2 3P0', [0.12535491, 1.05989413, 3.60720143]),
    ([203.795*u.AA, '3s2 3p 3d 3D3 -- 3s2 3p2 3P2'], '3s2 3p 3d 3P1 -- 3s2 3p2 3P0', [0.12535491, 1.05989413, 3.60720143]),
])
def test_line_ratio_wavelengths(fe_13, numerator, denominator, result):
    density = [1e8, 1e9, 1e10] * u.cm**(-3)
    ratio = fiasco.line_ratio(fe_13, numerator, denominator, density, use_two_ion_model=False)
    assert ratio.shape == fe_13.temperature.shape + density.shape
    assert u.allclose(ratio[np.argmax(fe_13.ionization_fraction), :], result)


def test_line_ratio_temperature(fe_13):
    density = 1e15 * u.K * u.cm**(-3) / fe_13.temperature
    ratio = fiasco.line_ratio(fe_13,
                              [203.795, 203.826]*u.AA,
                              202.044*u.AA,
                              density,
                              couple_density_to_temperature=True,
                              use_two_ion_model=False)
    assert ratio.shape == fe_13.temperature.shape + (1,)
    assert u.allclose(ratio.squeeze()[[0, 10, 19]], [3.04999452, 1.22415167, 0.3683284])
