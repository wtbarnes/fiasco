"""
Calculating a contribution function
===================================

This example shows how to calculate the contribution function of a particular
transition of a particular ion, at given temperatures and densities.
"""
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import quantity_support

from fiasco import Ion

quantity_support()

###############################################################################
# Specify the plasma properties; note that an `~fiasco.Ion` has to be created with a
# range of temperatures, but the density is only used later in the contribution
# function calculation.
Te = np.geomspace(0.1, 100, 51) * u.MK
ne = 1e8 * u.cm**-3

###############################################################################
# Create the ion object
ion = Ion('O 5+', Te)
print(ion)

###############################################################################
# Calculate the contribution function, and get the index of the transition
# closest to the desired wavelength
contribution_func = ion.contribution_function(ne)
wlen = 1031.93 * u.Angstrom

###############################################################################
# Because the contribution function is calculated for all transitions at once,
# we need to get the index of the transition closest to the
# wavelength specified earlier.
transitions = ion.transitions.wavelength[~ion.transitions.is_twophoton]
idx = np.argmin(np.abs(transitions - wlen))

###############################################################################
# Plot the result
plt.plot(Te, contribution_func[:, 0, idx],
         label=f'{ion.atomic_symbol} {ion.charge_state}+ {wlen}')
plt.title('Contribution function')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
