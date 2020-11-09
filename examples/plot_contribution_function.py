"""
Calculating a contribution function
===================================

This example shows how to calculate the contribution function of a particular
transition of a particular ion, at given temperatures and densities.
"""
###############################################################################
# Import required modules, and set up support for plotting quantities
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import quantity_support
quantity_support()

from fiasco import Ion, IonCollection

###############################################################################
# Specify the plasma properties; note that an `Ion` has to be created with a
# range of temperatures, but the density is only used later in the contribution
# function calculation.
Te = np.geomspace(0.1, 100, 51) * u.MK
ne = 1e8 * u.cm**-3

###############################################################################
# Create the ion object
ion_name = 'O 5+'
ion = Ion(ion_name, Te)
print(ion)

###############################################################################
# Because the contribution function is calculated for all transitions at once,
# we need a function to get the index of the transition closest to the
# wavelength specified earlier.
def get_idx(ion, wlen):
    """
    Get the wavelength and index of that wavelength closest to *wlen*.
    """
    wlens = ion.transitions.wavelength[~ion.transitions.is_twophoton]
    idx = np.argmin(np.abs(wlens - wlen))
    return idx, wlens[idx]

###############################################################################
# Calculate the contribution function, and get the index of the transition
# closest to the desired wavelength
contribution_func = ion.contribution_function(ne)
wlen = 1031.93 * u.Angstrom
transition_idx, wlen = get_idx(ion, wlen)

###############################################################################
# Plot the result
fig, ax = plt.subplots(tight_layout=True)
ax.plot(Te, contribution_func[:, 0, transition_idx], label=f'{ion_name} {wlen}')

ax.set_title('Contribution function')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
plt.show()
