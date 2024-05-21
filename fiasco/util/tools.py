"""
Numerical tools
"""
import astropy.units as u
import numpy as np

from functools import partial
from scipy.interpolate import splev, splrep

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
    x_1 = np.atleast_1d(x_1)
    x_2 = np.atleast_1d(x_2)
    return np.array([np.where(x_1==x)[0] for x in x_2]).squeeze()


def vectorize_where_sum(x_1, x_2, y, axis=None):
    """
    Find all occurrences of one array in another and sum over a third

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
    r"""
    Convert scaled Burgess-Tully :cite:p:`burgess_analysis_1992` parameters to physical quantities.

    For a scaled temperature, :math:`x` and scaled effective collision strength
    :math:`y`, the effective collision strength can be calculated as a function
    of the scaled energy :math:`U=k_BT_e/\Delta E_{ij}` the ratio between the thermal
    energy and the energy of the transition :math:`ij`.

    There are 6 different scaling types, depending on the type of transition. This scaling
    is explained in detail in section 5 of :cite:t:`burgess_analysis_1992`.
    The scaled temperatures and collision strengths are related to :math:`U` and :math:`\Upsilon` by,

    * type 1

      .. math::

            x = 1 - \frac{\ln C}{\ln{(U + C)}},\quad
            y = \frac{\Upsilon}{\log(U + e)}

    * type 2

      .. math::

            x = \frac{U}{U + C},\quad
            y = \Upsilon

    * type 3

      .. math::

            x = \frac{U}{U + C},\quad
            y = (U + 1)\Upsilon

    * type 4

      .. math::

            x = 1 - \frac{\ln C}{\ln{(U + C)}},\quad
            y = \frac{\Upsilon}{\log(U + C)}

    * type 5

      .. math::

            x = \frac{U}{U + C},\quad
            y = \Upsilon U

    * type 6

      .. math::

            x = \frac{U}{U + C},\quad
            y = \log_{10}\Upsilon

    where :math:`C` is a scaling constant that is different for each transition. Note that :cite:t:`burgess_analysis_1992`
    only considered scaling types 1 through 4. Types 5 and 6 correspond to dielectron and proton
    transitions, respectively.

    To "descale" the scaled effective collision strengths that are stored in the database,
    a spline fit is computed to the new :math:`x` as computed from :math:`U` and then
    the relationship between :math:`\Upsilon` and :math:`y` is inverted to get
    :math:`\Upsilon` as a function of :math:`U`.

    Parameters
    ----------
    x : `array-like`
        Scaled temperature. First dimension should have length ``n``, the number of
        transitions. The second dimension will be the number of spline points, but may
        be different for each row. If each row has ``l`` spline points, ``x`` should
        have shape ``(n,l)``. If they are not all equal, ``x`` will have shape ``(n,)``.
    y : `array-like`
        Scaled collision strength. Must have the same dimensions as ``x``.
    energy_ratio : `array-like`
        Ratio between the thermal energy and that of each transition with shape ``(n,m)``,
        where ``m`` is the dimension of the temperature array.
    c : `array-like`
        Scaling constant for each transition with shape ``(n,)``
    scaling_type : `array-like`
        The type of descaling to apply for each transition with shape ``(n,)``. Must be between
        1 and 6

    Returns
    -------
    upsilon : `array-like`
        Descaled collision strength or cross-section with the same shape as ``energy_ratio``.
    """
    # NOTE: Arrays with staggered number of columns, which have an 'object'
    # dtype (denoted by 'O') appear to be 1D, but should not be cast to 2D
    # as this will actually add an extra dimension and throw off the function
    # mapping
    # NOTE: We need to first work out whether the array is ragged or not in order
    # to set the dtype of the resulting array. Not doing so is now deprecated in
    # numpy. See https://github.com/wtbarnes/fiasco/issues/120
    is_ragged = False
    if isinstance(x, list):
        if not all([_x.shape[0] == x[0].shape[0] for _x in x]):
            is_ragged = True
    elif isinstance(x, np.ndarray):
        if x.dtype == np.dtype('O'):
            is_ragged = True
    else:
        raise TypeError(f'x has unsupported type {type(x)}')
    if is_ragged:
        x = np.asarray(x, dtype='O')
        y = np.asarray(y, dtype='O')
    else:
        x = np.atleast_2d(x)
        y = np.atleast_2d(y)
    energy_ratio = np.atleast_2d(u.Quantity(energy_ratio).to_value(u.dimensionless_unscaled))
    c = u.Quantity(c).to_value(u.dimensionless_unscaled)

    out = np.zeros(energy_ratio.shape)
    xnew = np.zeros(energy_ratio.shape)

    for s_type in np.unique(scaling_type):
        idxs = scaling_type == s_type
        xnew[idxs, :] = _xnew(energy_ratio[idxs, :], c[idxs], s_type).T

    # Use list(map()) here to allow varying shaped inputs for x, y
    splrep_szero = partial(splrep, s=0)
    nots = np.array(list(map(splrep_szero, x, y)), dtype=object)
    splev_derzero = partial(splev, der=0)
    out = np.array(list(map(splev_derzero, xnew, nots)), dtype=object)

    for s_type in np.unique(scaling_type):
        idxs = scaling_type == s_type
        if s_type == 1:
            out[idxs, ...] *= np.log(energy_ratio[idxs, ...] + np.e)
        elif s_type == 3:
            out[idxs, ...] /= (energy_ratio[idxs, ...] + 1.0)
        elif s_type == 4:
            out[idxs, ...] *= np.log(energy_ratio[idxs, ...].T + c[idxs]).T
        elif s_type == 5:
            out[idxs, ...] /= energy_ratio[idxs, ...]
        elif s_type == 6:
            out[idxs, ...] = 10**out[idxs]

    return out
