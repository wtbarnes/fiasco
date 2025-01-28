"""
Test Gaunt factor functionality
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

from fiasco.gaunt import GauntFactor
from fiasco.util.exceptions import MissingDatasetException

temperature = np.logspace(5, 8, 100)*u.K

@pytest.fixture
def ion(hdf5_dbase_root):
    return fiasco.Ion('Fe 5', temperature, hdf5_dbase_root=hdf5_dbase_root)

@pytest.fixture
def gaunt_factor(hdf5_dbase_root):
    return GauntFactor(hdf5_dbase_root=hdf5_dbase_root)

def test_repr(gaunt_factor):
    assert 'Gaunt factor' in gaunt_factor.__repr__()

@pytest.mark.parametrize(('ionization_stage', 'zeta'), [
    (2, 32.0),
    (16, 18.0),
    (24, 8.0),
    (26, 2.0),
])
def test_zeta0(hdf5_dbase_root, ionization_stage, zeta):
    iron = fiasco.Element('Fe', temperature=temperature, hdf5_dbase_root=hdf5_dbase_root)
    z0 = iron[ionization_stage].gaunt_factor._zeta_0(iron.atomic_number, iron[ionization_stage].charge_state)
    assert u.isclose(z0, zeta)

@pytest.mark.parametrize(('charge_state', 'expected'),
                        [(0, 0.0), (1, 1.41409406), (2, 1.3265169), (3, 1.27165875)])
def test_gaunt_factor_free_free_integrated(gaunt_factor, charge_state, expected):
    gf = gaunt_factor.free_free_integrated(temperature, charge_state)
    assert gf.shape == temperature.shape
    # This value has not been tested for correctness
    assert u.allclose(gf[0], expected * u.dimensionless_unscaled)

@pytest.mark.requires_dbase_version('>= 9.0.1')
@pytest.mark.parametrize(('charge_state','index','expected'),
                        [(1, 0, 1.40965757), (1, 75, 1.20459895),
                         (2, 0, 1.32203268), (2, 75, 1.28680539),
                         (3, 0, 1.2673757), (3, 75, 1.34662286),
                        ])
def test_gaunt_factor_free_free_integrated_itoh(gaunt_factor, charge_state, index, expected):
    gf = gaunt_factor.free_free_integrated(temperature, charge_state, use_itoh=True)
    assert gf.shape == temperature.shape
    # This value has not been tested for correctness
    assert u.allclose(gf[index], expected * u.dimensionless_unscaled)

@pytest.mark.requires_dbase_version('< 9.0.1')
def test_free_free_integrated_itoh_missing_data(gaunt_factor):
    gf = gaunt_factor._free_free_itoh_integrated(temperature, 1)
    assert gf.shape == temperature.shape
    assert np.isnan(gf[0])
    with pytest.raises(MissingDatasetException):
        gf = gaunt_factor._free_free_itoh_integrated_relativistic(temperature, 1)
    with pytest.raises(MissingDatasetException):
        gf = gaunt_factor._free_free_itoh_integrated_nonrelativistic(temperature, 1)

@pytest.mark.requires_dbase_version('>= 9.0.1')
@pytest.mark.parametrize(('T'), [(1e3,1e9)])
def test_gaunt_factor_free_free_integrated_itoh_invalid_temperature(gaunt_factor, T):
    gf = gaunt_factor._free_free_itoh_integrated_relativistic(T*u.K, 1)
    assert np.isnan(gf[0])
    gf = gaunt_factor._free_free_itoh_integrated_nonrelativistic(T*u.K, 1)
    assert np.isnan(gf[0])

def test_gaunt_factor_free_bound_nl_missing(gaunt_factor):
    #test cases where n or l is not in the klgfb data
    assert u.isclose(gaunt_factor.free_bound(0.5, 10, 1), 1.0 * u.dimensionless_unscaled)

@pytest.mark.parametrize(('ground_state', 'expected'), [(True, 55.18573076316151), (False, 11.849092513590998)])
def test_gaunt_factor_free_bound_integrated(ion, ground_state, expected):
    gf = ion.gaunt_factor.free_bound_integrated(ion.temperature,
                                                ion.atomic_number,
                                                ion.charge_state,
                                                ion.previous_ion()._fblvl['n'][0],
                                                ion.previous_ion().ionization_potential,
                                                ground_state=ground_state)
    assert gf.shape == ion.temperature.shape
    # These values have not been tested for correctness
    assert u.isclose(gf[20], expected * u.dimensionless_unscaled)

def test_free_bound_integrated_zero_charge(gaunt_factor):
    assert u.allclose(gaunt_factor.free_bound_integrated(temperature, 1, 0, 1, 13.6*u.eV), 0.0 * u.dimensionless_unscaled)

@pytest.mark.parametrize('gs', [True, False])
def test_free_bound_gaunt_factor_low_temperature(gs, hdf5_dbase_root):
    # At low temperatures (~1e4 K), exponential terms in the gaunt factor used to compute the
    # free-bound radiative loss can blow up. This just tests to make sure those are handled correctly
    ion = fiasco.Ion('N 8', np.logspace(4,6,100)*u.K, hdf5_dbase_root=hdf5_dbase_root)
    gf_fb_int = ion.gaunt_factor.free_bound_integrated(ion.temperature,
                                                       ion.atomic_number,
                                                       ion.charge_state,
                                                       ion.previous_ion()._fblvl['n'][0],
                                                       ion.previous_ion().ionization_potential,
                                                       ground_state=gs)
    assert not np.isinf(gf_fb_int).any()
