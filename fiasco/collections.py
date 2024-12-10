"""
Multi-ion container
"""
import astropy.units as u
import numpy as np

from astropy.convolution import convolve, Model1DKernel
from astropy.modeling.models import Gaussian1D

import fiasco

from fiasco.util import parse_ion_name
from fiasco.util.exceptions import MissingDatasetException

__all__ = ['IonCollection']


class IonCollection:
    """
    Container for holding multiple `~fiasco.Ion` instances.

    This container is most useful when needing to group many ions together in
    order to perform some aggregate calculation like a radiative loss curve or
    a composite spectrum made up of ions from many different elements.

    Parameters
    ----------
    *ions: `fiasco.Ion` or `fiasco.IonCollection`
        Entries can be either ions or collections of ion.
    """

    def __init__(self, *args):
        # Import here to avoid circular imports
        from fiasco import log
        self.log = log
        self._ion_list = []
        for item in args:
            if isinstance(item, fiasco.Ion):
                self._ion_list.append(item)
            elif isinstance(item, fiasco.IonCollection):
                self._ion_list += item._ion_list
            else:
                raise TypeError(f'{item} has unrecognized type {type(item)}',
                                'and cannot be added to collection.')
        if not all([all(self[0].temperature == ion.temperature) for ion in self]):
            raise ValueError('Temperatures for all ions in collection must be the same.')

    def __getitem__(self, value):
        ions = np.array(self._ion_list)[value]
        if isinstance(ions, fiasco.Ion):
            return ions
        else:
            return IonCollection(*ions)

    def __contains__(self, value):
        if isinstance(value, (str, tuple)):  # NOQA: UP038
            pair = parse_ion_name(value)
        elif isinstance(value, fiasco.Ion):
            pair = value._base_rep
        return pair in [i._base_rep for i in self._ion_list]

    def __add__(self, value):
        return IonCollection(self, value)

    def __radd__(self, value):
        return IonCollection(value, self)

    def __len__(self):
        return len(self._ion_list)

    def __repr__(self):
        ion_name_list = '\n'.join([i.ion_name for i in self._ion_list])
        return f"""Ion Collection
--------------
Number of ions: {len(self)}
Temperature range: [{self.temperature[0].to(u.MK):.3f}, {self.temperature[-1].to(u.MK):.3f}]

Available Ions
--------------
{ion_name_list}"""

    @property
    def temperature(self,):
        # Temperatures for all ions must be the same
        return self[0].temperature

    @u.quantity_input
    def free_free(self, wavelength: u.angstrom):
        r"""
        Compute combined free-free continuum emission (bremsstrahlung).

        .. note:: Both abundance and ionization fraction are included here

        The combined free-free continuum is given by,

        .. math::

            P_{ff}(\lambda,T) = \sum_{X,k}\mathrm{Ab}(X)f(X_{k})C_{ff, X_k}(\lambda,T)

        where :math:`\mathrm{Ab}(X)` is the abundance of element :math:`X`,
        :math:`f(X_{k})` is the ionization equilibrium of the ion,
        and :math:`C_{ff, X_k}(\lambda,T)` is the free-free emission of the ion
        as computed by `fiasco.Ion.free_free`.
        The sum is taken over all ions in the collection.

        Parameters
        ----------
        wavelength : `~astropy.units.Quantity`

        See Also
        --------
        fiasco.Ion.free_free
        """
        wavelength = np.atleast_1d(wavelength)
        free_free = u.Quantity(np.zeros(self.temperature.shape + wavelength.shape),
                               'erg cm^3 s^-1 Angstrom^-1')
        for ion in self:
            ff = ion.free_free(wavelength)
            try:
                abundance = ion.abundance
                ionization_fraction = ion.ionization_fraction
            except MissingDatasetException as e:
                self.log.warning(f'{ion.ion_name} not included in free-free emission. {e}')
                continue
            else:
                free_free += ff * abundance * ionization_fraction[:, np.newaxis]
        return free_free

    @u.quantity_input
    def free_bound(self, wavelength: u.angstrom, **kwargs):
        r"""
        Compute combined free-bound continuum emission.

        .. note:: Both abundance and ionization fraction are included here.

        The combined free-bound continuum is given by,

        .. math::

            P_{fb}(\lambda,T) = \sum_{X,k}\mathrm{Ab}(X)f(X_{k+1})C_{fb, X_k}(\lambda,T)

        where :math:`\mathrm{Ab}(X)` is the abundance of element :math:`X`,
        :math:`f(X_{k+1})` is the ionization equilibrium of the recombining ion
        :math:`X_{k+1}`,
        and :math:`C_{fb, X_k}(\lambda,T)` is the free-bound emission of the recombined
        ion :math:`X_k` as computed by `fiasco.Ion.free_bound`.
        The sum is taken over all ions in the collection.

        Parameters
        ----------
        wavelength : `~astropy.units.Quantity`

        See Also
        --------
        fiasco.Ion.free_bound
        """
        wavelength = np.atleast_1d(wavelength)
        free_bound = u.Quantity(np.zeros(self.temperature.shape + wavelength.shape),
                                'erg cm^3 s^-1 Angstrom^-1')
        for ion in self:
            try:
                fb = ion.free_bound(wavelength, **kwargs)
                abundance = ion.abundance
                # NOTE: the free-bound emissivity gets multiplied by the population
                # fraction of the recombining ion, that is, the ion with one higher
                # charge state.
                ionization_fraction = ion.next_ion().ionization_fraction
            except MissingDatasetException as e:
                self.log.warning(f'{ion.ion_name} not included in free-bound emission. {e}')
                continue
            else:
                free_bound += fb * abundance * ionization_fraction[:, np.newaxis]
        return free_bound

    @u.quantity_input
    def two_photon(self, wavelength: u.angstrom, electron_density: u.cm**-3, **kwargs):
        r"""
        Compute the two-photon continuum emission.

        .. note:: Both abundance and ionization equilibrium are included here.

        The combined two-photon continuum is given by

        .. math::

            P_{2p}(\lambda,T,n_{e}) = \sum_{X,k}\mathrm{Ab}(X)f(X_{k})C_{2p, X_k}(\lambda,T,n_{e})

        where :math:`\mathrm{Ab}(X)` is the abundance of element :math:`X`,
        :math:`f(X_{k})` is the ionization equilibrium of the emitting ion :math:`X_{k}`,
        and :math:`C_{fb, X_k}(\lambda,T)` is the two-photon emission of the
        ion :math:`X_k` as computed by `fiasco.Ion.two_photon`.
        The sum is taken over all ions in the collection.

        Parameters
        ----------
        wavelength : `~astropy.units.Quantity`
        electron_density: `~astropy.units.Quantity`

        See Also
        --------
        fiasco.Ion.two_photon
        """
        wavelength = np.atleast_1d(wavelength)
        electron_density = np.atleast_1d(electron_density)
        two_photon = u.Quantity(np.zeros(self.temperature.shape + electron_density.shape + wavelength.shape),
                                'erg cm^3 s^-1 Angstrom^-1')
        for ion in self:
            try:
                tp = ion.two_photon(wavelength, electron_density, **kwargs)
            except MissingDatasetException as e:
                self.log.warning(f'{ion.ion_name} not included in two-photon emission. {e}')
                continue
            else:
                two_photon += tp * ion.abundance * ion.ionization_fraction[:, np.newaxis, np.newaxis]
        return two_photon

    @u.quantity_input
    def spectrum(self, density: u.cm**(-3), emission_measure: u.cm**(-5), wavelength_range=None,
                 bin_width=None, kernel=None, **kwargs):
        """
        Calculate spectrum for multiple ions

        .. warning:: This function is still experimental and may be removed or significantly
                     refactored in future releases.

        Parameters
        ----------
        density : `~astropy.units.Quantity`
            Electron number density
        emission_measure : `~astropy.units.Quantity`
            Column emission measure
        wavelength_range : `~astropy.units.Quantity`, optional
            Tuple of bounds on which transitions to include. Default includes all
        bin_width : `~astropy.units.Quantity`, optional
            Wavelength resolution to bin intensity values. Default to 1/100 of range
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
            # If no elvlc and wgfa data, cannot calculate spectra
            try:
                wave = ion.transitions.wavelength[~ion.transitions.is_twophoton]
            except MissingDatasetException:
                self.log.warning(f'No transition data available for {ion.ion_name}')
                continue
            else:
                i_wavelength, = np.where(np.logical_and(wave >= wavelength_range[0],
                                                        wave <= wavelength_range[1]))
            # Skip if no transitions in this range
            if i_wavelength.shape[0] == 0:
                continue
            # If no scups data, cannot calculate spectra
            try:
                intens = ion.intensity(density, emission_measure, **kwargs)
            except MissingDatasetException:
                self.log.warning(f'No collision data available for {ion.ion_name}')
                continue
            if wavelength is None:
                wavelength = wave[i_wavelength].value
                intensity = intens[:, :, i_wavelength].value
            else:
                wavelength = np.concatenate((wavelength, wave[i_wavelength].value))
                intensity = np.concatenate((intensity, intens[:, :, i_wavelength].value), axis=2)

        if wavelength is None:
            raise ValueError('No collision or transition data available for any ion in collection.')

        if np.any(np.isinf(wavelength_range)):
            wavelength_range = u.Quantity([wavelength.min(), wavelength.max()], wave.unit)

        # Setup bins
        if bin_width is None:
            bin_width = np.diff(wavelength_range)[0]/100.
        num_bins = int((np.diff(wavelength_range)[0]/bin_width).value)
        wavelength_edges = np.linspace(*wavelength_range.value, num_bins+1)
        # Setup convolution kernel
        if kernel is None:
            self.log.warning('Using 0.1 Angstroms (ie. not the actual thermal width) for thermal broadening')
            std = 0.1*u.angstrom  # FIXME: should default to thermal width
            std_eff = (std/bin_width).value  # Scale sigma by bin width
            # Kernel size must be odd
            x_size = int(8*std_eff)+1 if (int(8*std_eff) % 2) == 0 else int(8*std_eff)
            m = Gaussian1D(amplitude=1./np.sqrt(2.*np.pi)/std.value, mean=0., stddev=std_eff)
            kernel = Model1DKernel(m,  x_size=x_size)

        # FIXME: This is very inefficient! Vectorize somehow...
        spectrum = np.zeros(intensity.shape[:2]+(num_bins,))
        for i in range(spectrum.shape[0]):
            for j in range(spectrum.shape[1]):
                tmp, _ = np.histogram(wavelength, bins=wavelength_edges, weights=intensity[i, j, :])
                spectrum[i, j, :] = convolve(tmp, kernel, normalize_kernel=False)

        wavelength = (wavelength_edges[1:] + wavelength_edges[:-1])/2. * wave.unit
        spectrum = spectrum * intens.unit / bin_width.unit
        return wavelength, spectrum

    @u.quantity_input
    def radiative_loss(self, density: u.cm**(-3), use_itoh=False, **kwargs) -> u.Unit('erg cm3 s-1'):
        r"""
        Calculate the total wavelength-integrated radiative loss rate including the
        bound-bound, free-bound, and free-free emission contributions

        .. note::  The calculation does not include two-photon continuum emission, which is also
                   neglected in the CHIANTI IDL routines.

        Parameters
        ----------
        density : `~astropy.units.Quantity`
            Electron number density
        use_itoh : `bool`, optional
            Whether to use Gaunt factors taken from :cite:t:`itoh_radiative_2002` for the calculation
            of free-free emission.  Defaults to false.

        Returns
        -------
        rad_loss_total : `~astropy.units.Quantity`
            The total bolometric radiative loss rate
        """
        rad_loss_bound_bound = self.bound_bound_radiative_loss(density, **kwargs)
        rad_loss_free_free = self.free_free_radiative_loss(use_itoh=use_itoh)
        rad_loss_free_bound = self.free_bound_radiative_loss()

        rad_loss_total = (rad_loss_bound_bound
                          + rad_loss_free_free[:, np.newaxis]
                          + rad_loss_free_bound[:, np.newaxis])

        return rad_loss_total

    @u.quantity_input
    def bound_bound_radiative_loss(self, density, **kwargs) -> u.Unit('erg cm3 s-1'):
        r"""
        Calculate the radiative loss rate from bound-bound emission (line emission)
        integrated over wavelength.

        Parameters
        ----------
        density : `~astropy.units.Quantity`
            Electron number density

        Returns
        -------
        rad_loss : `~astropy.units.Quantity`
            The bolometric bound-bound radiative loss rate per unit emission measure
        """
        density = np.atleast_1d(density)
        if kwargs.get('couple_density_to_temperature', False):
            shape = self.temperature.shape + (1,)
        else:
            shape = self.temperature.shape + density.shape
        rad_loss = u.Quantity(np.zeros(shape), 'erg cm^3 s^-1')
        for ion in self:
            try:
                g = ion.contribution_function(density, **kwargs)
            except MissingDatasetException as e:
                self.log.warning(f'{ion.ion_name} not included in the bound-bound emission. {e}')
                continue
            rad_loss += g.sum(axis=2)
        return rad_loss

    @u.quantity_input
    def free_free_radiative_loss(self, use_itoh=False) -> u.Unit('erg cm3 s-1'):
        r"""
        Calculate the radiative loss rate from free-free emission (bremsstrahlung)
        integrated over wavelength.

        Parameters
        ----------
        use_itoh : `bool`, optional
            Whether to use Gaunt factors taken from :cite:t:`itoh_radiative_2002`.
            Defaults to false.

        Returns
        -------
        rad_loss : `~astropy.units.Quantity`
            The bolometric free-free radiative loss rate per unit emission measure

        See Also
        -------
        fiasco.GauntFactor.free_free_integrated
        """
        free_free = u.Quantity(np.zeros(self.temperature.shape), 'erg cm^3 s^-1')
        for ion in self:
            try:
                ff = ion.free_free_radiative_loss(use_itoh=use_itoh)
                abundance = ion.abundance
                ionization_fraction = ion.ionization_fraction
            except MissingDatasetException as e:
                self.log.warning(f'{ion.ion_name} not included in free-free radiative loss. {e}')
                continue
            else:
                free_free += ff * abundance * ionization_fraction
        return free_free

    @u.quantity_input
    def free_bound_radiative_loss(self) ->  u.Unit('erg cm3 s-1'):
        r"""
        Calculate the radiative loss rate from free-bound emission (collisional recombination)
        integrated over wavelength.

        Returns
        -------
        rad_loss : `~astropy.units.Quantity`
            The bolometric free-bound radiative loss rate per unit emission measure
        """
        free_bound = u.Quantity(np.zeros(self.temperature.shape), 'erg cm^3 s^-1')
        for ion in self:
            try:
                fb = ion.free_bound_radiative_loss()
                abundance = ion.abundance
                ionization_fraction = ion.ionization_fraction
            except MissingDatasetException as e:
                self.log.warning(f'{ion.ion_name} not included in free-bound radiative loss. {e}')
                continue
            else:
                free_bound += fb * abundance * ionization_fraction
        return free_bound
