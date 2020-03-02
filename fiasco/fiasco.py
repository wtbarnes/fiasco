"""
Package-level functions
"""
import functools
import warnings

import numpy as np
from scipy.interpolate import interp1d
import astropy.units as u
import plasmapy.atomic
from plasmapy.atomic.exceptions import InvalidParticleError

import fiasco
from fiasco.io import DataIndexer

__all__ = ['list_elements', 'list_ions', 'proton_electron_ratio']


def list_elements(hdf5_dbase_root, sort=True):
    """
    List all available elements in the CHIANTI database.
    """
    elements = []
    root = DataIndexer.create_indexer(hdf5_dbase_root, '/')
    for f in root.fields:
        try:
            elements.append(plasmapy.atomic.atomic_symbol(f.capitalize()))
        except InvalidParticleError:
            continue
    if sort:
        elements = sorted(elements, key=lambda x: plasmapy.atomic.atomic_number(x))
    return elements


@functools.lru_cache()
def list_ions(hdf5_dbase_root, sort=True):
    """
    List all available ions in the CHIANTI database
    """
    ions = []
    root = DataIndexer.create_indexer(hdf5_dbase_root, '/')
    for f in root.fields:
        try:
            el = plasmapy.atomic.atomic_symbol(f.capitalize())
            for i in root[f].fields:
                if f == i.split('_')[0]:
                    ions.append(f"{el} {i.split('_')[1]}")
        except InvalidParticleError:
            continue
    # Optional because adds significant overhead
    if sort:
        ions = sorted(ions, key=lambda x: (plasmapy.atomic.atomic_number(x.split()[0]),
                                           int(x.split()[1])))
    return ions


@u.quantity_input
def proton_electron_ratio(temperature: u.K, **kwargs):
    """
    Calculate ratio between proton and electron densities as a function of temperature
    according to Eq. 7 of [1]_.

    Parameters
    ----------
    temperature : `~astropy.units.Quantity`

    See Also
    --------
    fiasco.Ion : Accepts same keyword arguments for setting database and dataset names

    References
    ----------
    .. [1] Young, P. et al., 2003, ApJS, `144 135 <http://adsabs.harvard.edu/abs/2003ApJS..144..135Y>`_
    """
    h_2 = fiasco.Ion('H +1', temperature, **kwargs)
    numerator = h_2.abundance * h_2._ioneq[h_2._dset_names['ioneq_filename']]['ionization_fraction']
    denominator = u.Quantity(np.zeros(numerator.shape))
    for el_name in list_elements(h_2.hdf5_dbase_root):
        el = fiasco.Element(el_name, temperature, **kwargs)
        abundance = el.abundance
        if abundance is None:
            warnings.warn(f'Not including {el.atomic_symbol}. Abundance not available.')
            continue
        for ion in el:
            ioneq = ion._ioneq[ion._dset_names['ioneq_filename']]['ionization_fraction']
            if ioneq is None:
                warnings.warn(f'Not including {ion.ion_name}. Ionization fraction not available.')
                continue
            denominator += ioneq * abundance * ion.charge_state

    ratio = numerator / denominator
    interp = interp1d(ion._ioneq[ion._dset_names['ioneq_filename']]['temperature'].value,
                      ratio.value,
                      kind='linear',
                      bounds_error=False,
                      fill_value=(ratio[0], ratio[-1]))

    return u.Quantity(interp(temperature))
