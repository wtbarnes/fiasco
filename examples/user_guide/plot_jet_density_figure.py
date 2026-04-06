"""
========================================
Fe XII Density Diagnostics with EIS Data
========================================

This example shows how to compute a density diagnostic from Hinode/EIS data
using `fiasco`.

The goal is to reproduce the main elements of Fig. 7 from
`Temperature and density structure of a recurring active region jet <https://www.aanda.org/articles/aa/full_html/2017/02/aa28796-16/aa28796-16.html>`__
by Mulay et al. (2017).

The exact density values differ from the paper because this example uses the
current CHIANTI database through `fiasco`, together with precomputed EIS
products generated from reprocessed data.
"""
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from astropy.utils.data import download_file
from astropy.visualization import quantity_support

import fiasco

quantity_support()

###############################################################################
# First download the precomputed EIS products. These two FITS files are hosted
# in a separate data repository so that this example does not need to run the
# EIS fitting step with `eispac <https://eispac.readthedocs.io/en/latest/>`__.
data_base_url = 'https://media.githubusercontent.com/media/nabobalis/data/main/fiasco'
observed_path = download_file(
        f'{data_base_url}/jet_footpoint_fe12_observed.fits',
        cache=True,
        pkgname='fiasco',
        show_progress=False,
)

intensity_path = download_file(
        f'{data_base_url}/jet_footpoint_fe12_intensities.fits',
        cache=True,
        pkgname='fiasco',
        show_progress=False,
)


###############################################################################
# Load the observed Fe XII 195 map and the fitted Fe XII 186 and 195 intensity
# maps from the downloaded FITS files.
with fits.open(observed_path) as observed_hdul:
    observed_195 = observed_hdul['OBSERVED_195'].data

with fits.open(intensity_path) as intensity_hdul:
    intensity_186 = intensity_hdul['INTENSITY_186'].data
    intensity_195 = intensity_hdul['INTENSITY_195'].data
    intensity_header = intensity_hdul['INTENSITY_195'].header

###############################################################################
# Compute the theoretical Fe XII density-sensitive ratio curve.
temperature = 10**6.2 * u.K
temperature_grid = np.array([temperature.value]) * temperature.unit
density = np.logspace(8, 12, 80) * u.cm**-3

fe12 = fiasco.Ion('Fe XII', temperature_grid)
ratio_curve = fiasco.line_ratio(
    fe12,
    [186.854, 186.887] * u.angstrom,
    [195.119, 195.179] * u.angstrom,
    density,
    use_two_ion_model=False,
).squeeze()

###############################################################################
# Map the observed Fe XII intensity ratio onto the theoretical curve to derive
# the density map.
observed_ratio = np.divide(
    intensity_186,
    intensity_195,
    out=np.full_like(intensity_186, np.nan),
    where=intensity_195 > 0,
) * u.dimensionless_unscaled
is_finite = np.isfinite(ratio_curve.value)
sort_index = np.argsort(ratio_curve.value[is_finite])
density_map = u.Quantity(
    np.interp(
        observed_ratio.to_value(u.dimensionless_unscaled),
        ratio_curve.value[is_finite][sort_index],
        density.to_value('cm-3')[is_finite][sort_index],
        left=np.nan,
        right=np.nan,
    ),
    'cm-3',
)

###############################################################################
# Plot the theoretical ratio curve, the observed Fe XII map, and the derived
# density map for the jet footpoint.
fig, axes = plt.subplot_mosaic(
    [['curve', 'curve'],
     ['observed', 'density']],
    figsize=(10, 8),
    layout='constrained',
)

axes['curve'].plot(density, ratio_curve, color='black', lw=2)
axes['curve'].set_title('Fe XII $(186.854 + 186.887) / (195.119 + 195.179)$')
axes['curve'].set_xlabel(r'Electron density [$\mathrm{cm^{-3}}$]')
axes['curve'].set_ylabel('Ratio')
axes['curve'].set_xscale('log')

observed_image = axes['observed'].imshow(observed_195, origin='lower', cmap='magma', aspect="auto")
axes['observed'].set_title('Observed Fe XII 195')
axes['observed'].set_xlabel('X [Pixels]')
axes['observed'].set_ylabel('Y [Pixels]')

density_image = axes['density'].imshow(np.log10(density_map.to_value('cm-3')), origin='lower', cmap='viridis', aspect="auto")
axes['density'].set_title('Derived Electron Density')
axes['density'].set_xlabel('X [Pixels]')
axes['density'].set_ylabel('Y [Pixels]')

fig.colorbar(
    observed_image,
    cax=axes['observed'].inset_axes([1.02, 0.0, 0.04, 1.0]),
    label='Wavelength-integrated Fe XII 195',
)

fig.colorbar(
    density_image,
    cax=axes['density'].inset_axes([1.02, 0.0, 0.04, 1.0]),
    label=r'$\log_{10}(n_e / \mathrm{cm^{-3}})$',
)

plt.show()
