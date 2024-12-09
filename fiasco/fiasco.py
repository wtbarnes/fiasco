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

__all__ = ['list_elements', 'list_ions', 'proton_electron_ratio', 'get_isoelectronic_sequence']


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
