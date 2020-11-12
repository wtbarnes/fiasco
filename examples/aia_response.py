"""
Understanding AIA temperature response with `aiapy`
===================================================

In this example, we will use the `~aiapy` package to convolve
the Fe XVIII contribution function with the instrument response
function to understand the temperature sensitivity of the 94
angstrom channel.
"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
from aiapy.response import Channel
from astropy.visualization import quantity_support
quantity_support()

import fiasco

#############################################################
# First, create the `~aiapy.response.Channel` object for the
# 94 angstrom channel and compute the wavelength response.
ch = Channel(94*u.angstrom)
response = ch.wavelength_response() * ch.plate_scale

#############################################################
# Plot the wavelength response
plt.plot(ch.wavelength, response)
plt.xlim(ch.channel + [-4, 4]*u.angstrom)
plt.axvline(x=ch.channel, ls='--', color='k')

############################################################
# Next, we construct for the `~fiasco.Ion` object for Fe 18. Note
# that we choose to use the coronal abundances of
# `Feldman et al. (1992) <https://ui.adsabs.harvard.edu/abs/2012SoPh..275...41B/>`_.
temperature = 10.**(np.arange(4.5, 8, 0.05)) * u.K
fe18 = fiasco.Ion('Fe 18', temperature,
                  abundance_filename='sun_coronal_1992_feldman')

############################################################
# Compute contribution function,
#
# .. math:: G(n,T,\lambda) = 0.83\mathrm{Ab}\frac{hc}{\lambda}N_{\lambda}A_{\lambda}f\frac{1}{n}\quad \mathrm{[erg~cm^3~s^{-1}]}
#
# for each transition of Fe 18 at a single density.
# `~fiasco.Ion.contribution_function` also accepts an array of
# densities, but this requires significantly more computation.
g = fe18.contribution_function(1e9 * u.cm**(-3), include_protons=False)

############################################################
# Get the corresponding transition wavelengths
transitions = fe18.transitions.wavelength[~fe18.transitions.is_twophoton]
energy = const.h * const.c / transitions / u.photon

############################################################
# Interpolate the response function to the transition wavelengths
f = interp1d(ch.wavelength, response, fill_value='extrapolate')
response_transitions = f(transitions) * response.unit

############################################################
# In order to compute the temperature response function, we integrate the
# wavelength response function, weighted by the contribution function, over
# the given wavelength range. We divide by :math:`hc/\\lambda` in order to
# convert from units of energy to photons. The factor of :math:`0.83` is a
# relative scaling factor for the abundance of H and is not included in the
# temperature responses computed by `aia_get_response.pro`. For more
# more information on the AIA wavelength response calculation,
# see `Boerner et al. (2012) <https://ui.adsabs.harvard.edu/abs/2012SoPh..275...41B/>`_.
K = (g / energy * response_transitions).sum(axis=2) / (4*np.pi*u.steradian) / 0.83
K = K.squeeze().to('cm5 ct pix-1 s-1')

############################################################
# Plot the effective temperature response function for the 94
# channel.
plt.plot(temperature, K)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-29, 1e-26)
plt.xlim(1e5, 3e7)
plt.show()

############################################################
# Note that this represents only the contribution of Fe 18 to the 94
# angstrom bandpass for a single density. In order to construct the full
# response functions, one would need to repeat this procedure for every
# ion listed in the CHIANTI atomic database and for all AIA channels. The
# assumption of constant pressure is also often used such that the density
# varies inversely with the pressure.
