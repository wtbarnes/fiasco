"""
Test collection functionality
"""
import astropy.units as u
import numpy as np
import pytest

import fiasco

temperature = np.logspace(4, 8, 100)*u.K
wavelength = np.arange(1, 500, 1) * u.Angstrom


@pytest.fixture
def ion(hdf5_dbase_root):
    return fiasco.Ion('Fe 2', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def another_ion(hdf5_dbase_root):
    return fiasco.Ion('Ca 2', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def element(hdf5_dbase_root):
    return fiasco.Element('hydrogen', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def another_element(hdf5_dbase_root):
    return fiasco.Element('helium', temperature, hdf5_dbase_root=hdf5_dbase_root)


@pytest.fixture
def collection(hdf5_dbase_root):
    return fiasco.IonCollection(fiasco.Ion('H 1', temperature, hdf5_dbase_root=hdf5_dbase_root),
                                fiasco.Ion('He 2', temperature, hdf5_dbase_root=hdf5_dbase_root))


@pytest.fixture
def another_collection(hdf5_dbase_root):
    # This collection is for the continuum tests because we already have the needed
    # files.
    fe_5 = fiasco.Ion('Fe 5', temperature, hdf5_dbase_root=hdf5_dbase_root)
    return fe_5 + fe_5.next_ion()


def test_create_collection_from_ions(ion, another_ion):
    assert isinstance(ion + another_ion, fiasco.IonCollection)
    assert isinstance(fiasco.IonCollection(ion, another_ion), fiasco.IonCollection)


def test_create_collection_from_elements(element, another_element):
    assert isinstance(element + another_element, fiasco.IonCollection)
    assert isinstance(fiasco.IonCollection(element, another_element), fiasco.IonCollection)


def test_create_collection_from_mixture(ion, element, collection):
    assert isinstance(ion + element + collection, fiasco.IonCollection)
    # This tests the reflected operator
    assert isinstance(collection + ion, fiasco.IonCollection)
    assert isinstance(fiasco.IonCollection(ion, element, collection), fiasco.IonCollection)


def test_create_collection_from_collection(collection):
    assert isinstance(fiasco.IonCollection(collection), fiasco.IonCollection)


def test_getitem(ion, another_ion, element, another_element):
    collection = fiasco.IonCollection(ion, another_ion, element, another_element)
    assert collection[0] == ion
    assert collection[1] == another_ion
    assert isinstance(collection[1:2], fiasco.IonCollection)
    assert collection[1:2][0] == another_ion
    assert isinstance(collection[:2], fiasco.IonCollection)
    for i in collection:
        assert isinstance(i, fiasco.Ion)
    assert isinstance(collection[[1,2]], fiasco.IonCollection)


def test_contains(collection, hdf5_dbase_root):
    assert 'H 1' in collection
    assert 'hydrogen 1' in collection
    assert 'hydrogen +0' in collection
    ion = fiasco.Ion('H 1', temperature, hdf5_dbase_root=hdf5_dbase_root)
    assert ion in collection


def test_length(collection):
    assert len(collection) == len(collection._ion_list)


@pytest.mark.parametrize('wavelength', [wavelength, wavelength[50]])
def test_free_free(another_collection, wavelength):
    ff = another_collection.free_free(wavelength)
    assert ff.shape == temperature.shape + wavelength.shape if wavelength.shape else (1,)
    index_w = 50 if wavelength.shape else 0
    index_t = 24  # This is approximately where the ionization fraction for Fe V peaks
    assert u.allclose(ff[index_t, index_w], 3.2914969734961024e-42 * u.Unit('erg cm3 s-1 Angstrom-1'))


@pytest.mark.parametrize('wavelength', [wavelength, wavelength[50]])
def test_free_bound(another_collection, wavelength):
    fb = another_collection.free_bound(wavelength)
    assert fb.shape == temperature.shape + wavelength.shape if wavelength.shape else (1,)
    index_w = 50 if wavelength.shape else 0
    index_t = 24  # This is approximately where the ionization fraction for Fe V peaks
    assert u.allclose(fb[index_t, index_w], 3.057781475607237e-36 * u.Unit('erg cm3 s-1 Angstrom-1'))


@pytest.mark.requires_dbase_version('>= 8')
@pytest.mark.parametrize(('wavelength','index_w'), [(wavelength,450), (wavelength[450],0)])
def test_two_photon(collection, wavelength, index_w, hdf5_dbase_root):
    # add Li III to the test to include an ion that throws a MissingDatasetException
    collection = collection + fiasco.Ion('Li III', collection.temperature, hdf5_dbase_root=hdf5_dbase_root)
    tp = collection.two_photon(wavelength, electron_density = 1e10 * u.cm**(-3))
    wavelength = np.atleast_1d(wavelength)
    assert tp.shape == temperature.shape + (1, ) + wavelength.shape
    # This value has not been checked for correctness
    assert u.allclose(tp[30, 0, index_w], 3.4862925e-27 * u.Unit('erg cm3 s-1 Angstrom-1'))


@pytest.mark.requires_dbase_version('>= 8')
def test_radiative_loss(collection, hdf5_dbase_root):
    # add Li III to the test to include an ion that throws a MissingDatasetException
    collection = collection + fiasco.Ion('Li III', collection.temperature, hdf5_dbase_root=hdf5_dbase_root)
    density = [1e9,1e10,1e11] * u.cm**-3
    rl = collection.radiative_loss(density)
    assert rl.shape == (len(temperature), len(density))
    # These values have not been checked for correctness
    u.allclose(rl[0], [3.90235371e-24, 4.06540902e-24, 4.08411295e-24] * u.erg * u.cm**3 / u.s)


@pytest.mark.requires_dbase_version('>= 8')
def test_radiative_loss_bound_bound(collection, hdf5_dbase_root):
    # add Li III to the test to include an ion that throws a MissingDatasetException
    collection = collection + fiasco.Ion('Li III', collection.temperature, hdf5_dbase_root=hdf5_dbase_root)
    density = [1e9,1e10,1e11] * u.cm**-3
    rl = collection.bound_bound_radiative_loss(density)
    assert rl.shape == (len(temperature), len(density))
    # These values have not been checked for correctness
    u.allclose(rl[0], [3.90235371e-24, 4.06540902e-24, 4.08411295e-24] * u.erg * u.cm**3 / u.s)

@pytest.mark.requires_dbase_version('>= 8')
@pytest.mark.parametrize(('index','expected'),
                        [(0, 2.72706455e-35),
                        (75, 5.59153955e-31),])
def test_radiative_loss_free_free(collection, index, expected):
    rl = collection.free_free_radiative_loss()
    assert rl.shape == collection.temperature.shape
    # This value has not been checked for correctness
    u.isclose(rl[index], expected * u.erg * u.cm**3 / u.s)

@pytest.mark.requires_dbase_version('>= 9.0.1')
@pytest.mark.parametrize(('index','expected'),
                        [(0, 2.71800458e-35),
                        (75, 5.57346003e-31),])
def test_radiative_loss_free_free_itoh(collection, index, expected):
    """
    For database versions >= 9.0.1, the Itoh Gaunt factors give a different
    result for the free-free radiative loss.
    """
    rl = collection.free_free_radiative_loss(use_itoh=True)
    assert rl.shape == collection.temperature.shape
    # This value has not been checked for correctness
    u.isclose(rl[index], expected * u.erg * u.cm**3 / u.s)


@pytest.mark.requires_dbase_version('>= 8')
def test_radiative_loss_free_bound(collection):
    rl = collection.free_bound_radiative_loss()
    assert rl.shape == collection.temperature.shape
    # This value has not been checked for correctness
    u.isclose(rl[0], 1.13808317e-33 * u.erg * u.cm**3 / u.s)


def test_radiative_loss_temperature_density_coupling(another_collection):
    p0 = 1e15 * u.K * u.cm**(-3)
    density = p0 / another_collection.temperature
    rl = another_collection.radiative_loss(density, couple_density_to_temperature=True)
    assert rl.shape == another_collection.temperature.shape + (1,)


@pytest.mark.requires_dbase_version('>= 8')
def test_spectrum(hdf5_dbase_root):
    i1 = fiasco.Ion('H 1', 1 * u.MK, hdf5_dbase_root=hdf5_dbase_root)
    i2 = fiasco.Ion('Fe 5', 1 * u.MK, hdf5_dbase_root=hdf5_dbase_root)
    c = i1 + i2
    density = 1e9 * u.cm**-3
    em = 1e29 * u.cm**-5
    w, spec = c.spectrum(density, em)
    assert spec.shape == (1, 1, ) + w.shape
    # Add an ion with no spectral information
    i3 = fiasco.Ion('H 2', 1 * u.MK, hdf5_dbase_root=hdf5_dbase_root)
    c += i3
    w2, spec2 = c.spectrum(density, em)
    assert spec2.shape == (1, 1, ) + w2.shape
    assert np.all(spec == spec2)


@pytest.mark.requires_dbase_version('>= 8')
def test_build_line_list(hdf5_dbase_root):
    h1 = fiasco.Ion('H 1', temperature, hdf5_dbase_root=hdf5_dbase_root)
    fe5 = fiasco.Ion('Fe 5', temperature, hdf5_dbase_root=hdf5_dbase_root)
    # H 2 has no transition data and should be skipped
    h2 = fiasco.Ion('H 2', temperature, hdf5_dbase_root=hdf5_dbase_root)
    c = h1 + fe5 + h2
    density = 1e9 * u.cm**-3
    line_list = c.build_line_list(density)
    n_h1 = h1.transitions.is_bound_bound.sum()
    n_fe5 = fe5.transitions.is_bound_bound.sum()
    assert len(line_list) == n_h1 + n_fe5
    assert set(line_list['ion_name']) == {'H I', 'Fe V'}
    assert line_list['contribution_function'].shape == (n_h1 + n_fe5, temperature.shape[0], 1)
    assert 'intensity' not in line_list.colnames
    assert u.allclose(line_list.meta['temperature'], temperature)
    # Level indices point back at the levels of the ion
    is_fe5 = line_list['ion_name'] == 'Fe V'
    assert np.array_equal(line_list['upper_level'][is_fe5],
                          fe5.transitions.upper_level[fe5.transitions.is_bound_bound])
    assert np.all(line_list['upper_level'] > line_list['lower_level'])
    # Consistent with the contribution function of each ion
    g_fe5 = fe5.contribution_function(density)
    assert u.allclose(line_list['contribution_function'][line_list['ion_name'] == 'Fe V'],
                      np.moveaxis(g_fe5, -1, 0))
    # Selecting by abundance removes iron
    line_list = c.build_line_list(density, minimum_abundance=1e-3)
    assert set(line_list['ion_name']) == {'H I'}
    # Selecting by wavelength range keeps only the lines in range
    wavelength_range = [900, 1300] * u.angstrom
    line_list = c.build_line_list(density, wavelength_range=wavelength_range)
    assert np.all(line_list['wavelength'] >= wavelength_range[0])
    assert np.all(line_list['wavelength'] <= wavelength_range[1])
    assert 'H I' in set(line_list['ion_name'])


@pytest.mark.requires_dbase_version('>= 8')
def test_build_line_list_intensity(hdf5_dbase_root):
    fe5 = fiasco.Ion('Fe 5', temperature, hdf5_dbase_root=hdf5_dbase_root)
    c = fiasco.IonCollection(fe5)
    density = [1e9, 1e10] * u.cm**-3
    dem = fiasco.get_dem_model('quiet_sun', hdf5_dbase_root=hdf5_dbase_root)
    line_list = c.build_line_list(density, emission_measure=dem)
    assert line_list['intensity'].shape == (len(line_list), 2)
    assert line_list['intensity'].unit.is_equivalent(u.erg / u.cm**2 / u.s / u.steradian)
    assert np.all(line_list['intensity'] >= 0)


def test_build_line_list_no_valid_ions(hdf5_dbase_root):
    c = fiasco.IonCollection(fiasco.Ion('H 2', temperature, hdf5_dbase_root=hdf5_dbase_root))
    with pytest.raises(ValueError, match='No lines found for any ion in collection.'):
        c.build_line_list(1e9 * u.cm**-3)


def test_spectrum_no_valid_ions(hdf5_dbase_root):
    # Consider the case of an collection with ions with no spectral information
    c2 = fiasco.IonCollection(fiasco.Ion('H 2', 1 * u.MK, hdf5_dbase_root=hdf5_dbase_root))
    with pytest.raises(ValueError, match='No collision or transition data available for any ion in collection.'):
        c2.spectrum(1e9 * u.cm**-3, 1e29 * u.cm**-5)


def test_unequal_temperatures_raise_value_error(hdf5_dbase_root):
    first_ion = fiasco.Ion('Fe 12', [1e6, 1e7]*u.K, hdf5_dbase_root=hdf5_dbase_root)
    second_ion = fiasco.Ion('Fe 9', [1e4, 1e5]*u.K, hdf5_dbase_root=hdf5_dbase_root)
    with pytest.raises(ValueError):
        fiasco.IonCollection(first_ion, second_ion)


def test_create_with_wrong_type_raise_type_error(ion, collection):
    with pytest.raises(TypeError):
        fiasco.IonCollection(ion, 0)
    with pytest.raises(TypeError):
        collection + 0


def test_collections_repr(collection):
    assert isinstance(collection.__repr__(), str)
    for ion in collection:
        assert ion.ion_name in collection.__repr__()
