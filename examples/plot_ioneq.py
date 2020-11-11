"""
Ionization fractions in equilibrium
===============================================

This example shows how to compute the ionization fraction as a function of
temperature, assuming equilibrium, for both a single ion as well as a whole
element.
"""
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.visualization import quantity_support
quantity_support()

from fiasco import Ion, Element

################################################
# Create an Fe XV ion for a temperature array between :math:`10^4` and
# :math:`10^8` K.
temperature = 10**np.arange(4, 8, 0.01) * u.K
fe15 = Ion('Fe 15', temperature)

################################################
# The CHIANTI database includes tabulated ionization equilibria for
# all ions in the database. The `ioneq` attribute on each
# `~fiasco.Ion` object returns the tabulated population
# fractions interpolated onto the `temperature` array.
plt.plot(fe15.temperature, fe15.ioneq)
plt.xscale('log')
plt.title(f'{fe15.roman_name} Equilibrium Ionization')
plt.show()

################################################
# Note that these population fractions returned by `~fiasco.Ion.ioneq` are
# stored in the CHIANTI database and therefore are set to NaN
# for temperatures outside of the tabulated temperature data given in CHIANTI.
# If you need to calculate the population fractions over a wider
# temperature range, it is better to do so by calculating the
# ionization and recombination rates as a function of temperature for every
# ion of the element and then solving the associated system of equations.
# This can be done by creating a `~fiasco.Element` object and then calling
# the `~fiasco.Element.equilibrium_ionization` method.
fe = Element('iron', temperature)
ioneq = fe.equilibrium_ionization()

################################################
# Plot the population fraction of each ion as a function of temperature.
for ion in fe:
    plt.plot(fe.temperature, ioneq[:, ion.charge_state])
plt.xscale('log')
plt.title(f'{fe.atomic_symbol} Equilibrium Ionization')
plt.show()

################################################
# We can then compare tabulated and calculated results. Note that the two
# may not be equal due to differences in the rates when the tabulated results
# were calculated or due to artifacts from the interpolation.
plt.plot(fe.temperature, ioneq[:, fe[14].charge_state], label='calculated')
plt.plot(fe15.temperature, fe15.ioneq, label='interpolated')
plt.xlim(1e6, 1e7)
plt.xscale('log')
plt.legend()
plt.show()
