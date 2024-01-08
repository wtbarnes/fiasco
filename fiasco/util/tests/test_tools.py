import numpy as np
import pytest

from fiasco.util import burgess_tully_descale, vectorize_where


@pytest.mark.parametrize(('a', 'b', 'c'), [
    (np.arange(10), [3, 5, 9], np.array([3, 5, 9])),
    (np.arange(10), [6], np.array([6])),
    (np.arange(10), 6, np.array([6])),
])
def test_vectorize_where(a, b, c):
    indices = vectorize_where(a, b)
    assert (indices == c).all()


def test_vectorize_where_result_in_array():
    a = np.arange(20)
    b = np.array([5, 7, 17, 19])
    c = np.random.rand(*(15, 10,)+a.shape)
    indices = vectorize_where(a, b)
    assert c[:, :, indices].shape == (15, 10) + b.shape


@pytest.mark.parametrize('x', [
    np.linspace(0, 1, 10).reshape(1,10),  # one array
    np.tile(np.linspace(0, 1, 10), (2, 1)),  # 2D array
    [np.linspace(0,1,10), np.linspace(0,1,10)], # list with equal lengths
    np.array([np.linspace(0, 1, 5), np.linspace(0, 1, 9)], dtype=object),  # staggered array
    [np.linspace(0,1,5), np.linspace(0,1,9)], # list with unequal lengths
])
def test_burgess_tully(x):
    if isinstance(x, list):
        y = [_x**1.6 + 10 for _x in x]
    else:
        y = x**1.6 + 10
    energy_ratio = np.ones((len(x), 10))
    c = 1*np.ones((len(x),))
    for scaling_type in range(1, 7):
        st = scaling_type*np.ones(c.shape)
        ups = burgess_tully_descale(x, y, energy_ratio, c, st)
        assert ups.shape == energy_ratio.shape


def test_burgess_tully_different_scaling_types():
    x = np.tile(np.linspace(0, 1, 10), (6, 1))
    y = x**1.6 + 10
    energy_ratio = np.ones(x.shape[:1] + (10,))
    c = 1*np.ones(x.shape[:1])
    scaling_type = np.arange(1, 7)
    ups = burgess_tully_descale(x, y, energy_ratio, c, scaling_type)
    assert ups.shape == energy_ratio.shape
