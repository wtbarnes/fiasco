import numpy as np

from fiasco.util import burgess_tully_descale


def test_burgess_tully_one_array():
    x = np.linspace(0, 1, 10)
    y = x**1.6 + 10
    energy_ratio = np.ones((1,) + x.shape)
    c = 1
    for scaling_type in range(1, 7):
        ups = burgess_tully_descale(x, y, energy_ratio, c, scaling_type)
        assert ups.shape == energy_ratio.shape


def test_burgess_tully_2d_array():
    x = np.tile(np.linspace(0, 1, 10), (2, 1))
    y = x**1.6 + 10
    energy_ratio = np.ones(x.shape[:1] + (10,))
    c = 1*np.ones(x.shape[:1])
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


def test_burgess_tully_staggered_array():
    x = np.array([np.linspace(0, 1, 5), np.linspace(0, 1, 9)])
    y = x**1.6 + 10
    energy_ratio = np.ones(x.shape[:1] + (10,))
    c = 1*np.ones(x.shape)
    for scaling_type in range(1, 7):
        st = scaling_type*np.ones(c.shape)
        ups = burgess_tully_descale(x, y, energy_ratio, c, st)
        assert ups.shape == energy_ratio.shape
