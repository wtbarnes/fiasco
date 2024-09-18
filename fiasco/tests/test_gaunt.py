"""
Test Gaunt factor functionality
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

temperature = np.logspace(5, 8, 100)*u.K

@pytest.fixture
def ion(hdf5_dbase_root):
    return fiasco.Ion('Fe 5', temperature, hdf5_dbase_root=hdf5_dbase_root)

@pytest.fixture
def gaunt_factor(hdf5_dbase_root):
    return fiasco.GauntFactor()

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
                        [(1, 1.41409406), (2, 1.3265169), (3, 1.27165875)])
def test_gaunt_factor_free_free_total(gaunt_factor, charge_state, expected):
    gf = gaunt_factor.free_free_total(temperature, charge_state)
    assert gf.shape == temperature.shape
    # This value has not been tested for correctness
    assert u.allclose(gf[0], expected * u.dimensionless_unscaled)

@pytest.mark.requires_dbase_version('>= 9.0.1')
@pytest.mark.parametrize(('charge_state','itoh','relativistic','index','expected'),
                        [(3, True, True, 0, 1.2673757), (3, True, False, 0, 1.2673757),
                        (3, False, False, 0, 1.27165875), (3, True, True, 75, 1.3466916),
                        (3, True, False, 75, 1.34662286), (3, False, False, 75, 1.35061768)])
def test_gaunt_factor_free_free_total_itoh(gaunt_factor, charge_state, itoh, relativistic, index, expected):
    gf = gaunt_factor.free_free_total(temperature, charge_state, itoh, relativistic)
    assert gf.shape == temperature.shape
    # This value has not been tested for correctness
    assert u.allclose(gf[index], expected * u.dimensionless_unscaled)

def test_gaunt_factor_free_bound_total(ion):
    ion_gf_0 = ion.gaunt_factor.free_bound_total(ion.temperature, ion.atomic_number, ion.charge_state, ion.previous_ion()._fblvl['n'][0], ion.previous_ion().ip, ground_state=True)
    ion_gf_1 = ion.gaunt_factor.free_bound_total(ion.temperature, ion.atomic_number, ion.charge_state, ion.previous_ion()._fblvl['n'][0], ion.previous_ion().ip, ground_state=False)
    assert ion_gf_0.shape == ion.temperature.shape
    # These values have not been tested for correctness
    assert u.isclose(ion_gf_0[20], 55.18573076316151 * u.dimensionless_unscaled)
    assert u.isclose(ion_gf_1[20], 11.849092513590998 * u.dimensionless_unscaled)

@pytest.mark.parametrize('gs', [True, False])
def test_free_bound_gaunt_factor_low_temperature(gs, hdf5_dbase_root):
    # At low temperatures (~1e4 K), exponential terms in the gaunt factor used to compute the
    # free-bound radiative loss can blow up. This just tests to make sure those are handled correctly
    ion = fiasco.Ion('N 8', np.logspace(4,6,100)*u.K, hdf5_dbase_root=hdf5_dbase_root)
    gf_fb_total = ion.gaunt_factor.free_bound_total(ion.temperature, ion.atomic_number, ion.charge_state, ion.previous_ion()._fblvl['n'][0], ion.previous_ion().ip, ground_state=gs)
    assert not np.isinf(gf_fb_total).any()
