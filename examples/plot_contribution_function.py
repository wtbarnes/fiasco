import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import quantity_support
quantity_support()

from fiasco import Ion, IonCollection

ion_name = 'C 2+'
wlen = 977.03 * u.Angstrom

Te = np.logspace(4, 8, 51) * u.K
ne = 1e8 * u.cm**-3

ion = Ion(ion_name, Te)

def get_idx(ion, wlen):
    """
    Get the wavelength and index of that wavelength closest to *wlen*.
    """
    wlens = ion.transitions.wavelength[~ion.transitions.is_twophoton]
    idx = np.argmin(np.abs(wlens - wlen))
    return idx, wlens[idx]


contribution_func = ion.contribution_function(ne)
idx, _ = get_idx(ion, wlen)

fig, ax = plt.subplots(tight_layout=True)
ax.plot(Te, contribution_func[:, 0, idx], label=f'{ion_name} {wlen}')

ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
plt.show()
