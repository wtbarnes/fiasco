"""
Package-level functions.
"""
import astropy.units as u
import numpy as np
import plasmapy.particles

from plasmapy.particles.exceptions import InvalidParticleError
from scipy.interpolate import interp1d

import fiasco

from fiasco.io import DataIndexer
from fiasco.util import parse_ion_name

__all__ = [
    'list_elements',
    'list_ions',
    'proton_electron_ratio',
    'get_isoelectronic_sequence',
    'map_ratio_to_quantity',
]


def list_elements(hdf5_dbase_root=None, sort=True):
    """
    List all available elements in the CHIANTI database.

    Parameters
    ----------
    hdf5_dbase_root: path-like, optional
        If not specified, will default to that specified in ``fiasco.defaults``.
    sort: `bool`, optional
        If True, sort the list of elements by increasing atomic number.
    """
    if hdf5_dbase_root is None:
        hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
    elements = []
    root = DataIndexer.create_indexer(hdf5_dbase_root, '/')
    for f in root.fields:
        try:
            elements.append(plasmapy.particles.atomic_symbol(f.capitalize()))
        except InvalidParticleError:
            continue
    if sort:
        elements = sorted(elements, key=lambda x: plasmapy.particles.atomic_number(x))
    return elements


def list_ions(hdf5_dbase_root=None, sort=True):
    """
    List all available ions in the CHIANTI database

    Parameters
    ----------
    hdf5_dbase_root: path-like, optional
        If not specified, will default to that specified in ``fiasco.defaults``.
    sort: `bool`, optional
        If True, sort the list of elements by increasing atomic number.
    """
    if hdf5_dbase_root is None:
        hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
    root = DataIndexer(hdf5_dbase_root, '/')
    # NOTE: get the list from the index if possible. This is ~30x faster
    try:
        ions = root['ion_index']
    except KeyError:
        ions = []
        for f in root.fields:
            try:
                el = plasmapy.particles.atomic_symbol(f.capitalize())
                for i in root[f].fields:
                    if f == i.split('_')[0]:
                        ions.append(f"{el} {i.split('_')[1]}")
            except InvalidParticleError:
                continue
    # Optional because adds significant overhead
    if sort:
        ions = sorted(ions, key=lambda x: (plasmapy.particles.atomic_number(x.split()[0]),
                                           int(x.split()[1])))
    # NOTE: when grabbing straight from the index and not sorting, the result will be
    # a numpy array. Cast to a list to make sure the return type is consistent for
    # all possible inputs
    return ions.tolist() if isinstance(ions, np.ndarray) else ions


def get_isoelectronic_sequence(element, hdf5_dbase_root=None):
    """
    List of ions in the isoelectronic sequence of ``element``.

    Ions in the same isoelectronic sequence, that is, ions that have the
    same number of bound electrons, :math:`Z - z`, often share common properties
    despite having different charge states and being from different elements.
    These so-called isoelectronic sequences are typically denoted by the element
    with an atomic number equal to the number of bound electrons, e.g. C II is in
    the boron isoelectronic sequence or equivalently may be said to be boron-like.
    Given the name of a sequence, as denoted by an element label, this function
    returns a list of all ions in that sequence.

    Parameters
    ----------
    element: `str`, `int`
        Name of sequence. Can be either the full name (e.g. "hydrogren"),
        the atomic symbol (e.g. "H") or the atomic number (e.g. 1)
    hdf5_dbase_root: path-like, optional
        If not specified, will default to that specified in ``fiasco.defaults``.
    """
    Z_iso = plasmapy.particles.atomic_number(element)
    all_ions = list_ions(hdf5_dbase_root=hdf5_dbase_root)

    def _is_in_sequence(ion):
        Z, z = parse_ion_name(ion)
        return Z_iso == (Z - z + 1)

    return [ion for ion in all_ions if _is_in_sequence(ion)]


@u.quantity_input
def proton_electron_ratio(temperature: u.K, **kwargs):
    """
    Calculate ratio between proton and electron densities as a function of temperature
    according to Eq. 7 of :cite:t:`young_chianti-atomic_2003`.

    Parameters
    ----------
    temperature : `~astropy.units.Quantity`

    See Also
    --------
    fiasco.Ion : Accepts same keyword arguments for setting database and dataset names
    """
    # Import here to avoid circular imports
    from fiasco import log
    h_2 = fiasco.Ion('H +1', temperature, **kwargs)
    numerator = h_2.abundance * h_2._ion_fraction[h_2._instance_kwargs['ionization_fraction']]['ionization_fraction']
    denominator = u.Quantity(np.zeros(numerator.shape))
    for el_name in list_elements(h_2.hdf5_dbase_root):
        el = fiasco.Element(el_name, temperature, **h_2._instance_kwargs)
        try:
            abundance = el.abundance
        except KeyError:
            abund_file = el[0]._instance_kwargs['abundance']
            log.warning(
                f'Not including {el.atomic_symbol}. Abundance not available from {abund_file}.')
            continue
        for ion in el:
            ionization_file = ion._instance_kwargs['ionization_fraction']
            # NOTE: We use ._ion_fraction here rather than .ionization_fraction to avoid
            # doing an interpolation to the temperature array every single time and instead only
            # interpolate once at the end.
            # It is assumed that the ionization_fraction temperature array for each ion is the same.
            try:
                ionization_fraction = ion._ion_fraction[ionization_file]['ionization_fraction']
                t_ionization_fraction = ion._ion_fraction[ionization_file]['temperature']
            except KeyError:
                log.warning(
                    f'Not including {ion.ion_name}. Ionization fraction not available from {ionization_file}.')
                continue
            denominator += ionization_fraction * abundance * ion.charge_state

    ratio = numerator / denominator
    f_interp = interp1d(t_ionization_fraction.to(temperature.unit).value,
                        ratio.value,
                        kind='linear',
                        bounds_error=False,
                        fill_value=(ratio[0], ratio[-1]))

    return u.Quantity(f_interp(temperature.value))


def map_ratio_to_quantity(observed_ratio,
                          quantity,
                          theoretical_ratio,
                          *,
                          bounds_error=False,
                          fill_value=np.nan):
    """
    Map an observed line ratio to the associated physical quantity using a theoretical curve.

    Parameters
    ----------
    observed_ratio : array-like or `~astropy.units.Quantity`
        Observed line ratio values. Must be dimensionless.
    quantity : `~astropy.units.Quantity`
        Quantity grid associated with ``theoretical_ratio``, e.g. an electron density grid.
    theoretical_ratio : array-like or `~astropy.units.Quantity`
        Theoretical line ratio sampled on ``quantity``. Must be dimensionless and monotonic.
    bounds_error : `bool`, optional
        If True, raise an exception when ``observed_ratio`` falls outside of the
        theoretical ratio range. By default, points outside the range are set by
        ``fill_value``.
    fill_value : scalar, tuple, or `~astropy.units.Quantity`, optional
        Value used outside of the interpolation interval. If a quantity is given,
        it must be convertible to the units of ``quantity``.

    Returns
    -------
    `~astropy.units.Quantity`
        Interpolated quantity with the same shape as ``observed_ratio``.
    """
    observed_ratio = u.Quantity(observed_ratio, u.dimensionless_unscaled)
    if not isinstance(quantity, u.Quantity):
        raise TypeError('quantity must be an astropy Quantity.')
    quantity = np.ravel(quantity)
    theoretical_ratio = np.ravel(
        u.Quantity(theoretical_ratio, u.dimensionless_unscaled).value
    )

    if quantity.shape != theoretical_ratio.shape:
        raise ValueError('quantity and theoretical_ratio must have the same shape.')

    is_finite = np.isfinite(quantity.value) & np.isfinite(theoretical_ratio)
    quantity = quantity[is_finite]
    theoretical_ratio = theoretical_ratio[is_finite]

    if quantity.size < 2:
        raise ValueError('At least two finite samples are required to interpolate a ratio curve.')

    diff = np.diff(theoretical_ratio)
    nonzero = diff != 0.0
    if not np.any(nonzero):
        raise ValueError('theoretical_ratio must vary over the provided quantity grid.')
    if np.all(diff[nonzero] < 0.0):
        quantity = quantity[::-1]
        theoretical_ratio = theoretical_ratio[::-1]
    elif not np.all(diff[nonzero] > 0.0):
        raise ValueError(
            'theoretical_ratio must be monotonic to map ratios back to a quantity. '
            'Restrict the grid to a monotonic interval first.'
        )

    theoretical_ratio, unique_index = np.unique(theoretical_ratio, return_index=True)
    quantity = quantity[unique_index]
    if quantity.size < 2:
        raise ValueError('Theoretical ratio must contain at least two unique samples.')

    if isinstance(fill_value, tuple):
        fill_value = tuple(
            value.to_value(quantity.unit) if isinstance(value, u.Quantity) else value
            for value in fill_value
        )
    elif isinstance(fill_value, u.Quantity):
        fill_value = fill_value.to_value(quantity.unit)

    interp = interp1d(
        theoretical_ratio,
        quantity.value,
        bounds_error=bounds_error,
        fill_value=fill_value,
    )
    return u.Quantity(interp(observed_ratio.to_value(u.dimensionless_unscaled)), quantity.unit)
