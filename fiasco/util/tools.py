"""
Numerical tools
"""
import numpy as np
import astropy.units as u

__all__ = ['vectorize_where', 'vectorize_where_sum']


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
    axis : int, optional
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
