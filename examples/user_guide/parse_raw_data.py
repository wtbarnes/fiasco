"""
Parsing Raw CHIANTI Data
=========================

This example shows how to directly parse the raw CHIANTI database
files.
"""
from fiasco.io import Parser

##################################################################
# While the main advantages of fiasco lie in its high-level
# interface to the CHIANTI data, some users may wish to just parse
# the raw data directly. fiasco provides a convenient interface to
# parsing any of the raw CHIANTI datafiles and provides detailed
# metadata about each datafile. Specifically, fiasco returns a
# `~astropy.table.QTable` object with appropriate units and
# descriptive names attached to each column.
#
# For example, say we want to parse the energy level file for Fe V
# (i.e. iron with four electrons missing)
p = Parser('fe_5.elvlc')
table = p.parse()
print(table)

###################################################
# The individual columns can easily be accessed as well.
print(table.colnames)
print(table['E_obs'])

#####################################################
# Each above column is an `~astropy.units.Quantity` object with
# units attached to it if appropriate. Metadata, including the
# original footer from the raw CHIANTI data and detailed
# descriptions of each of the columns, is included with each table,
print(table.meta.keys())
print(table.meta['footer'])
