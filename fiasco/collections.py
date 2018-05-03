"""
Multi-ion container
"""
import numpy as np
import astropy.units as u
from astropy.convolution import convolve, Model1DKernel
from astropy.modeling.models import Gaussian1D
import plasmapy

import fiasco

__all__ = ['IonCollection']


class IonCollection(object):
    """
    Container for holding multiple ions. Instantiate with ions, elements, or another
    ion collection.

    Examples
    --------
    """
    
    def __init__(self, *args, **kwargs):
        self._ion_list = []
        for item in args:
            if isinstance(item, fiasco.Ion):
                self._ion_list.append(item)
            elif isinstance(item, fiasco.IonCollection):
                self._ion_list += item._ion_list
            else:
                raise TypeError(f'{item} has unrecognized type {type(item)}',
                                'and cannot be added to collection.')
        # TODO: check for duplicates
        assert all([all(self[0].temperature == ion.temperature) for ion in self]), (
            'Temperatures for all ions in collection must be the same.')
        
    def __getitem__(self, value):
        return self._ion_list[value]
    
    def __contains__(self, value):
        if type(value) is str:
            el, ion = value.split()
            if '+' in ion:
                ion = int(ion.strip('+')) + 1
            value = f'{plasmapy.atomic.atomic_symbol(el)} {ion}'
        elif isinstance(value, fiasco.Ion):
            value = value.ion_name
        return value in [i.ion_name for i in self._ion_list]
    
    def __add__(self, value):
        return IonCollection(self, value)
    
    def __radd__(self, value):
        return IonCollection(value, self)

    @u.quantity_input
    def spectrum(self, density: u.cm**(-3), emission_measure: u.cm**(-5), wavelength_range=None,
                 bin_width=None, kernel=None, **kwargs):
        """
        Calculate spectrum for multiple ions

        Parameters
        ----------
        density : `~astropy.units.Quantity`
            Electron number density
        emission_measure : `~astropy.units.Quantity`
            Column emission measure
        wavelength_range : `~astropy.units.Quantity`, optional
            Tuple of bounds on which transitions to include. Default includes all
        bin_width : `~astropy.units.Quantity`, optional
            Wavelength resolution to bin intensity values. Default to 1/10 of range
        kernel : `~astropy.convolution.kernels.Model1DKernel`, optional
            Convolution kernel for computing spectrum. Default is gaussian kernel with thermal width

        Returns
        -------
        wavelength : `~astropy.units.Quantity`
            Continuous wavelength
        spectrum : `~astropy.units.Quantity`
            Continuous intensity distribution as a function of wavelength

        See Also
        --------
        fiasco.Ion.spectrum : Compute spectrum for a single ion
        """
        if wavelength_range is None:
            wavelength_range = u.Quantity([0, np.inf], 'angstrom')
        # Compute all intensities
        intensity, wavelength = None, None
        for ion in self:
            wave = ion.transitions.wavelength[~ion.transitions.is_twophoton]
            i_wavelength, = np.where(np.logical_and(wave >= wavelength_range[0],
                                                    wave <= wavelength_range[1]))
            # Skip if no transitions in this range
            if i_wavelength.shape[0] == 0:
                continue
            intens = ion.intensity(density, emission_measure, **kwargs)
            if wavelength is None:
                wavelength = wave[i_wavelength].value
                intensity = intens[:, :, i_wavelength].value
            else:
                wavelength = np.concatenate((wavelength, wave[i_wavelength].value))
                intensity = np.concatenate((intensity, intens[:, :, i_wavelength].value), axis=2)
        if np.any(np.isinf(wavelength_range)):
            wavelength_range = u.Quantity([wavelength.min(), wavelength.max()], wave.unit)
        # Setup bins
        if bin_width is None:
            bin_width = np.diff(wavelength_range)[0]/100.
        num_bins = int((np.diff(wavelength_range)[0]/bin_width).value)
        wavelength_edges = np.linspace(*wavelength_range.value, num_bins+1)
        # Setup convolution kernel
        if kernel is None:
            std = 0.1*u.angstrom  # FIXME: should default to thermal width
            std_eff = (std/bin_width).value  # Scale sigma by bin width
            # Kernel size must be odd
            x_size = int(8*std_eff)+1 if (int(8*std_eff) % 2) == 0 else int(8*std_eff)
            m = Gaussian1D(amplitude=1./np.sqrt(2.*np.pi)/std.value,
                           mean=0., stddev=std_eff)
            kernel = Model1DKernel(m,  x_size=x_size)
        # FIXME: This is very inefficient! Vectorize somehow...
        spectrum = np.zeros(intensity.shape[:2]+(num_bins,))
        for i in range(spectrum.shape[0]):
            for j in range(spectrum.shape[1]):
                tmp, _ = np.histogram(wavelength, bins=wavelength_edges, weights=intensity[i, j, :])
                spectrum[i, j, :] = convolve(tmp, kernel, normalize_kernel=False)
        wavelength = (wavelength_edges[1:] + wavelength_edges[:-1])/2.
        spectrum = spectrum * intens.unit / bin_width.unit

        return wavelength, spectrum

    @u.quantity_input
    def radiative_loss(self, density: u.cm**(-3), **kwargs):
        """
        Calculate radiative loss curve which includes multiple ions
        """
        rad_loss = u.Quantity(
            np.zeros(self[0].temperature.shape + density.shape), 'erg cm^3 s^-1')
        for ion in self:
            g = ion.contribution_function(density, **kwargs)
            if g is None:
                continue
            rad_loss += g.sum(axis=2)

        return rad_loss
