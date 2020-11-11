"""
Ionization and recombination rates
===================================

This example shows how to compute the recombination rate as a function of
temperature for Fe 16, including the radiative and dielectronic
components.
"""

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.visualization import quantity_support
quantity_support()

from fiasco import Ion

#################################################
# First, instantiate the Fe XVI ion object. This can be done in a number of
# ways, but we will use the spectroscopic roman numeral notation.
ion = Ion('Fe XVI', np.logspace(4, 8, 100) * u.K)

#################################################
# Calculate the recombination rate, which includes contributions from both
# radiative and dielectronic recombination. We can also compute these
# separately in order to understand over what temperature range each
# process dominates.
total_recombination = ion.recombination_rate()
dielectronic_recombination = ion.dielectronic_recombination_rate()
radiative_recombination = ion.radiative_recombination_rate()

#################################################
# Similarly, we can calculate the ionization rate which includes
# contributions from direction ionization and excitation
# autoionization
total_ionization = ion.ionization_rate()
direct_ionization = ion.direct_ionization_rate()
excitation_autoionization = ion.excitation_autoionization_rate()

#################################################
# Plot the recombination and ionization rates, including all
# components, as a function of temperature.
fig, ax = plt.subplots(tight_layout=True)
ax.plot(ion.temperature, total_recombination,
        label='Recombination', color='C0',)
ax.plot(ion.temperature, dielectronic_recombination,
        label='Dielectronic', color='C0', ls='--')
ax.plot(ion.temperature, radiative_recombination,
        label='Radiative', color='C0', ls=':')
ax.plot(ion.temperature, total_ionization,
        label='Ionization', color='C1')
ax.plot(ion.temperature, direct_ionization,
        label='Direct', color='C1', ls='--')
ax.plot(ion.temperature, excitation_autoionization,
        label='Excitation Autoionization', color='C1', ls=':')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1e-12, 1e-9)
ax.legend()
plt.show()
