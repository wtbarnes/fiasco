"""
Numerical tools
"""
from functools import partial
import numpy as np
from scipy.interpolate import splrep, splev
import astropy.units as u

__all__ = ['vectorize_where', 'vectorize_where_sum', 'burgess_tully_descale']


def vectorize_where(x_1, x_2):
    """
    Find indices of one array in another

    Parameters
    ----------
    x_1 : array-like
        Array to search through
    x_2 : array-like
        Values to search for
    """
    return np.vectorize(lambda a, b: np.where(a == b)[0], excluded=[0])(x_1, x_2)


def vectorize_where_sum(x_1, x_2, y, axis=None):
    """
    Find all occurences of one array in another and sum over a third

    Parameters
    ----------
    x_1 : array-like
        Array to search through
    x_2 : array-like
        Values to search for
    y : array-like
    axis : `int`, optional
        Axis to sum over
    """
    unit = None
    if isinstance(y, u.Quantity):
        unit = y.unit
        y = y.value
    if len(y.shape) == 2:
        signature = '()->(n)'
    elif len(y.shape) == 1:
        signature = '()->()'
    else:
        raise ValueError('y cannot have dimension greater than 2')
    collect = np.vectorize(lambda a, b, c: c[np.where(a == b)].sum(axis=axis),
                           excluded=[0, 2], signature=signature)
    return u.Quantity(collect(x_1, x_2, y), unit)


def _xnew(energy_ratio, c, scaling_type):
    energy_ratio = energy_ratio.T
    if scaling_type in [1, 4]:
        return 1.0 - np.log(c) / np.log(energy_ratio + c)
    elif scaling_type in [2, 3, 5, 6]:
        return energy_ratio / (energy_ratio + c)


def burgess_tully_descale(x, y, energy_ratio, c, scaling_type):
    """
    Convert scaled Burgess-Tully [1]_ parameters to physical quantities.

    For a scaled temperature, :math:`x` and scaled effective collision strength
    :math:`y`, the effective collision strength can be calculated as a function
    of the scaled energy :math:`U=k_BT_e/\Delta E_{ij}` the ratio between the thermal
    energy and the energy of the transition :math:`ij`.

    There are 6 different scaling types, depending on the type of transition. This scaling
    is explained in detail in section 5 of [1]_. For types 1 and 4, the scaled temperatures
    and collision strengths are related to :math:`U` and :math:`\\Upsilon` by,

    * type 1
      
      .. math::

            x = 1 - \\frac{\ln C}{\ln{(U + C)}},\quad
            y = \\frac{\\Upsilon}{\log(U + e)}

    * type 2

      .. math::

            x = \\frac{U}{U + C},\quad
            y = \\Upsilon

    * type 3

      .. math::

            x = \\frac{U}{U + C},\quad
            y = (U + 1)\\Upsilon

    * type 4

      .. math::

            x = 1 - \\frac{\ln C}{\ln{(U + C)}},\quad
            y = \\frac{\\Upsilon}{\log(U + C)}

    * type 5

      .. math::

            x = \\frac{U}{U + C},\quad
            y = \\Upsilon U

    * type 6

      .. math::

            x = \\frac{U}{U + C},\quad
            y = \log_{10}\\Upsilon

    where :math:`C` is a scaling constant that is different for each transition. Note that [1]_
    only considered scaling types 1 through 4. Types 5 and 6 correspond to dielectron and proton
    transitions, respectively.

    To "descale" the scaled effective collision strengths that are stored in the database,
    a spline fit is computed to the new :math:`x` as computed from :math:`U` and then
    the relationship between :math:`\\Upsilon` and :math:`y` is inverted to get
    :math:`\\Upsilon` as a function of :math:`U`.
    
    Parameters
    ----------
    x : `array-like`
        Scaled temperature. First dimension should have length ``n``, the number of
        transitions. The second dimension will be the number of spline points, but may
        be different for each row. If each row has ``l`` spline points, `x` should 
        have shape ``(n,l)``. If they are not all equal, `x` will have shape ``(n,)``.
    y : `array-like`
        Scaled collision strength. Must have the same dimensions as `x`.
    energy_ratio : `array-like`
        Ratio between the thermal energy and that of each transition with shape ``(n,m)``,
        where ``m`` is the dimension of the temperature array.
    c : `array-like`
        Scaling constant for each transiton with shape ``(n,)``
    scaling_type : `array-like`
        The type of descaling to apply for each transition with shape ``(n,)``. Must be between
        1 and 6

    Returns
    -------
    upsilon : `array-like`
        Descaled collision strength or cross-section with the same shape as `energy_ratio`.

    References
    ----------
    .. [1] Burgess, A. and Tully, J. A., 1992, A&A, `254, 436 <http://adsabs.harvard.edu/abs/1992A%26A...254..436B>`_
    """
    # NOTE: Arrays with staggered number of columns, which have an 'object'
    # dtype (denoted by 'O') appear to be 1D, but should not be cast to 2D
    # as this will actually add an extra dimension and throw off the function
    # mapping
    x = np.asarray(x)
    y = np.asarray(y)
    if x.dtype != np.dtype('O'):
        x = np.atleast_2d(x)
        y = np.atleast_2d(y)
    energy_ratio = np.atleast_2d(u.Quantity(energy_ratio).to_value(u.dimensionless_unscaled))
    c = u.Quantity(c).to_value(u.dimensionless_unscaled)

    out = np.zeros(energy_ratio.shape)
    xnew = np.zeros(energy_ratio.shape)

    for type in np.unique(scaling_type):
        idxs = scaling_type == type
        xnew[idxs, :] = _xnew(energy_ratio[idxs, :], c[idxs], type).T

    # Use list(map()) here to allow varying shaped inputs for x, y
    splrep_szero = partial(splrep, s=0)
    nots = np.array(list(map(splrep_szero, x, y)))
    splev_derzero = partial(splev, der=0)
    out = np.array(list(map(splev_derzero, xnew, nots)))

    for type in np.unique(scaling_type):
        idxs = scaling_type == type
        if type == 1:
            out[idxs, ...] *= np.log(energy_ratio[idxs, ...] + np.e)
        elif type == 3:
            out[idxs, ...] /= (energy_ratio[idxs, ...] + 1.0)
        elif type == 4:
            out[idxs, ...] *= np.log(energy_ratio[idxs, ...].T + c[idxs]).T
        elif type == 5:
            out[idxs, ...] /= energy_ratio[idxs, ...]
        elif type == 6:
            out[idxs, ...] = 10**out[idxs]

    return out
