import numpy as np

from fiasco.util import burgess_tully_descale_vectorize, burgess_tully_descale


def test_burgess_tully():
    # Check that normal and vectorized versions give the same result
    x = np.linspace(0, 1, 10)
    y = x**1.6 + 10
    energy_ratio = np.ones(x.shape)
    c = 1
    for scaling_type in range(1, 7):
        out = burgess_tully_descale(x, y, energy_ratio, c, scaling_type)
        out_vect = burgess_tully_descale_vectorize(x, y, energy_ratio, c, scaling_type)
        np.testing.assert_equal(out, out_vect[0])
