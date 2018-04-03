"""
Package-level functions
"""
import numpy as np
import astropy.units as u
import plasmapy.atomic

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
        except plasmapy.atomic.names.InvalidParticleError:
            continue
    return elements


@u.quantity_input
def proton_electron_ratio(temperature: u.K, **kwargs):
    """
    Calculate ratio between proton and electron densities as a function of temperature
    according to Eq. 7 of [1]_.

    Parameters
    ----------
    temperature : `~astropy.units.Quantity`

    Other Parameters
    ----------------
    See `~fiasco.Ion`
    
    References
    ----------
    .. [1] Young, P. et al., 2003, ApJS, `144 135 <http://adsabs.harvard.edu/abs/2003ApJS..144..135Y>`_
    """
    h_2 = fiasco.Ion('H +1', temperature, **kwargs)
    numerator = h_2.abundance * h_2.ioneq
    denominator = u.Quantity(np.zeros(numerator.shape))
    for el_name in list_elements():
        el = fiasco.Element(el_name, temperature, ion_kwargs=kwargs)
        for ion in el:
            denominator += ion.ioneq * ion.abundance * ion.charge_state
            
    ratio = numerator / denominator
    # Set out of range values to closest valid values
    indices_valid, = np.where(~np.isnan(ratio))
    ratio[:indices_valid[0]] = ratio[indices_valid[0]]
    ratio[indices_valid[-1]+1:] = ratio[indices_valid[-1]]

    return ratio
