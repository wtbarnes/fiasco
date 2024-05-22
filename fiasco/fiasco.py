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

__all__ = ['list_elements', 'list_ions', 'proton_electron_ratio']


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
    return ions.tolist() if type(ions) == np.ndarray else ions


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
    numerator = h_2.abundance * h_2._ioneq[h_2._instance_kwargs['ioneq_filename']]['ionization_fraction']
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
            ioneq_file = ion._instance_kwargs['ioneq_filename']
            # NOTE: We use ._ioneq here rather than .ioneq to avoid doing an interpolation to the
            # temperature array every single time and instead only interpolate once at the end.
            # It is assumed that the ioneq temperature array for each ion is the same.
            try:
                ioneq = ion._ioneq[ioneq_file]['ionization_fraction']
                t_ioneq = ion._ioneq[ioneq_file]['temperature']
            except KeyError:
                log.warning(
                    f'Not including {ion.ion_name}. Ionization fraction not available from {ioneq_file}.')
                continue
            denominator += ioneq * abundance * ion.charge_state

    ratio = numerator / denominator
    f_interp = interp1d(t_ioneq.to(temperature.unit).value,
                        ratio.value,
                        kind='linear',
                        bounds_error=False,
                        fill_value=(ratio[0], ratio[-1]))

    return u.Quantity(f_interp(temperature.value))
