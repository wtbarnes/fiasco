"""
Package-level functions
"""
import warnings

import numpy as np
from scipy.interpolate import interp1d
import astropy.units as u
import plasmapy.atomic
from plasmapy.utils import InvalidParticleError

import fiasco

__all__ = ['list_elements', 'proton_electron_ratio']


def list_elements():
    """
    List all available elements in the CHIANTI database.
    """
    dl = fiasco.DataIndexer(fiasco.defaults['hdf5_dbase_root'], '/')
    elements = []
    for f in dl.fields:
        try:
            elements.append(plasmapy.atomic.atomic_symbol(f.capitalize()))
        except InvalidParticleError:
            continue
    elements = sorted(elements, key=lambda x: plasmapy.atomic.atomic_number(x))
    return elements


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
    fiasco.Ion : Accepts same keyword arguments for setting dataset names

    References
    ----------
    .. [1] Young, P. et al., 2003, ApJS, `144 135 <http://adsabs.harvard.edu/abs/2003ApJS..144..135Y>`_
    """
    h_2 = fiasco.Ion('H +1', temperature, **kwargs)
    numerator = h_2.abundance * h_2._ioneq[h_2._dset_names['ioneq_filename']]['ionization_fraction']
    denominator = u.Quantity(np.zeros(numerator.shape))
    for el_name in list_elements():
        el = fiasco.Element(el_name, temperature, ion_kwargs=kwargs)
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
                      ratio.value, kind='linear', bounds_error=False,
                      fill_value=(ratio[0], ratio[-1]))

    return u.Quantity(interp(temperature))
