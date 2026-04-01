"""
Test package level functions
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco


def test_list_elements(hdf5_dbase_root):
    # FIXME: actually test the expected elements
    elements = fiasco.list_elements(hdf5_dbase_root)
    assert isinstance(elements, list)


def test_list_ions(hdf5_dbase_root):
    # FIXME: actually test the expected ions
    ions = fiasco.list_ions(hdf5_dbase_root)
    assert isinstance(ions, list)


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


def test_map_ratio_to_quantity():
    density = [1e8, 1e9, 1e10] * u.cm**(-3)
    theoretical_ratio = [0.2, 0.5, 0.8] * u.dimensionless_unscaled
    observed_ratio = [0.2, 0.65, 1.0]

    mapped_density = fiasco.map_ratio_to_quantity(observed_ratio, density, theoretical_ratio)

    assert type(mapped_density) is u.Quantity
    assert mapped_density.shape == (3,)
    assert mapped_density.unit == density.unit
    assert u.allclose(mapped_density[:2], [1e8, 5.5e9] * u.cm**(-3))
    assert np.isnan(mapped_density[-1].value)

    mapped_density = fiasco.map_ratio_to_quantity(0.35, density, theoretical_ratio[::-1])
    assert u.allclose(mapped_density, 5.5e9 * u.cm**(-3))


def test_map_ratio_to_quantity_non_monotonic_curve():
    density = [1e8, 1e9, 1e10] * u.cm**(-3)
    theoretical_ratio = [0.2, 0.8, 0.5] * u.dimensionless_unscaled

    with pytest.raises(ValueError, match='must be monotonic'):
        fiasco.map_ratio_to_quantity(0.35, density, theoretical_ratio)


def test_map_ratio_to_quantity_requires_quantity_axis():
    theoretical_ratio = [0.2, 0.5, 0.8] * u.dimensionless_unscaled

    with pytest.raises(TypeError, match='quantity must be an astropy Quantity'):
        fiasco.map_ratio_to_quantity(0.35, [1e8, 1e9, 1e10], theoretical_ratio)
