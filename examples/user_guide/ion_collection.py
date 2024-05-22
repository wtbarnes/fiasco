"""
Creating an ``IonCollection``
=============================

This example shows how to create an `~fiasco.IonCollection` object.
"""
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import quantity_support

import fiasco

# sphinx_gallery_thumbnail_number = 2

quantity_support()

#############################################################
# The `~fiasco.IonCollection` object allows you to create
# arbitrary collections of `~fiasco.Ion` objects. This
# provides an easy way to compute quantities that combine
# emission from multiple ions, such as a spectra or a
# radiative loss curve.
#
# We will create `~fiasco.IonCollection` containing both O VI
# and Fe XVIII
temperature = np.geomspace(1e4,1e8,100) * u.K
fe18 = fiasco.Ion('Fe XVIII', temperature)
o6 = fiasco.Ion('O VI', temperature)
col = fiasco.IonCollection(fe18, o6)
print(col)

#############################################################
# The only requirement for ions in the same collection is
# that they use the same ``temperature`` array. We can also
# index our collection in the same manner as a
# `~fiasco.Element` object to access the individual ions.
print(col[0])

#############################################################
# Or iterate through them.
for ion in col:
    print(ion.ion_name_roman)

#############################################################
# You can also create a collection using the addition
# operator. This is equivalent to using the
# `~fiasco.IonCollection` constructor as we did above.
col = fe18 + o6

#############################################################
# Furthermore, you can iteratively build a collection in this
# way as well. For example, to create a new collection that
# includes Fe XV from the collection we created above,
new_col = col + fiasco.Ion('Fe XV', temperature)
print(new_col)

#############################################################
# You can even add `~fiasco.Element` objects as well, which
# results in every ion of that element being added to your
# collection.
new_col = col + fiasco.Element('carbon', temperature)
print(new_col)

#############################################################
# There are several methods on `~fiasco.IonCollection` for
# easily computing quantities with contributions from multiple
# ions. As an example, let's compute the radiative loss for
# our collection containing Fe XVIII and O VI.
#
# First, let's create an `~fiasco.IonCollection` containing
# just Fe XVIII.
col = fiasco.IonCollection(fe18)

#############################################################
# Now, let's compute the radiative loss as a function of
# temperature for just Fe XVIII
density = 1e9*u.cm**(-3)
rl = col.radiative_loss(density)
plt.plot(col.temperature, rl)
plt.yscale('log')
plt.xscale('log')
plt.ylim(1e-30, 1e-20)

#############################################################
# Next, let's add O VI and again compute the radiative loss.
new_col = fe18 + o6
rl_new = new_col.radiative_loss(density)
plt.plot(col.temperature, rl, label=','.join([i.ion_name_roman for i in col]))
plt.plot(new_col.temperature, rl_new, label=','.join([i.ion_name_roman for i in new_col]))
plt.yscale('log')
plt.xscale('log')
plt.ylim(1e-30, 1e-20)
plt.legend()

#############################################################
# By comparing the radiative losses, we can clearly see at
# what temperatures the respective ions are dominating.
