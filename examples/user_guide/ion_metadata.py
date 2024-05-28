"""
Creating an ``Ion``
===================

This example shows how to create a `~fiasco.Ion` object and
access various pieces of metadata.
"""
import astropy.units as u
import numpy as np

from fiasco import Ion

######################################################
# The CHIANTI database is organized around individual ions, with
# multiple types of datafiles attached to each ion. In keeping with
# the organization of the database, the primary unit of the fiasco
# library is the `~fiasco.Ion` object which can be created in the
# following way,
temperature = np.logspace(5, 7, 100) * u.K
ion = Ion('Fe 15', temperature)
print(ion)

#######################################################
# This creates a `~fiasco.Ion` object for the Fe XV ion. Note also
# the same object can also be created in the following ways,
ion = Ion('iron 15', temperature)
ion = Ion('iron 14+', temperature)
ion = Ion('Fe XV', temperature)

#######################################################
# The `~fiasco.Ion` object holds several basic pieces of metadata
# related to the particular ion,
print(ion.element_name)
print(ion.atomic_symbol)
print(ion.atomic_number)
print(ion.ion_name)
print(ion.charge_state)
print(ion.ionization_stage)
print(ion.abundance)

#######################################################
# The `~fiasco.Ion` object can also be indexed like an array in
# order to get information about the energy levels.
print(ion[0])

#######################################################
# Each level also holds various bits of metadata
for i in range(5):
    lev = ion[i]
    print(f'Level {lev.level} {lev.configuration}, {lev.energy}')
