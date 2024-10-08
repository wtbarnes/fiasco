"""
Gaunt factor object.
"""
import astropy.constants as const
import astropy.units as u
import numpy as np

from scipy.interpolate import splev, splrep
from scipy.ndimage import map_coordinates
from scipy.special import polygamma

import fiasco

from fiasco.io.datalayer import DataIndexer
from fiasco.util import check_database, needs_dataset
from fiasco.util.exceptions import MissingDatasetException

__all__ = ['GauntFactor']

class GauntFactor:
    """
    Class for calculating the Gaunt factor for various continuum processes.

    The Gaunt factor is defined as the ratio of the true cross-section to the
    semi-classical Kramers cross-section, and thus is essentially a multiplicative
    correction for quantum mechanical effects.  It is a unitless quantity.

    Parameters
    ------------
    hdf5_dbase_root: path-like, optional
        Path to built database
    kwargs:
        All keyword arguments to `fiasco.util.check_database` are also supported here
    """
    def __init__(self, hdf5_dbase_root=None, **kwargs):
        if hdf5_dbase_root is None:
            self.hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
        else:
            self.hdf5_dbase_root = hdf5_dbase_root
        check_database(self.hdf5_dbase_root, **kwargs)

    def __repr__(self):
        return f"""Gaunt factor object
---------------------
HDF5 Database: {self.hdf5_dbase_root}"""

    @property
    def _gffgu(self):
        data_path = '/'.join(['continuum', 'gffgu'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _gffint(self):
        data_path = '/'.join(['continuum', 'gffint'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _itoh(self):
        data_path = '/'.join(['continuum', 'itoh'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _itoh_integrated_gaunt(self):
        data_path = '/'.join(['continuum', 'itoh_integrated_gaunt'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _itoh_integrated_gaunt_nonrel(self):
        data_path = '/'.join(['continuum', 'itoh_integrated_gaunt_nonrel'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _klgfb(self):
        data_path = '/'.join(['continuum', 'klgfb'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @u.quantity_input
    def free_free(self, temperature: u.K, wavelength: u.angstrom, atomic_number, charge_state) -> u.dimensionless_unscaled:
        r"""
        The Gaunt factor for free-free emission as a function of temperature and wavelength.

        The free-free Gaunt factor is calculated from a lookup table of temperature averaged
        free-free Gaunt factors from Table 2 of :cite:t:`sutherland_accurate_1998` as a function
        of :math:`\log{\gamma^2},\log{u}`, where :math:`\gamma^2=Z^2\mathrm{Ry}/k_BT`
        and :math:`u=hc/\lambda k_BT`.

        For the regime, :math:`6<\log_{10}(T)< 8.5` and :math:`-4<\log_{10}(u)<1`, the above
        prescription is replaced with the fitting formula of :cite:t:`itoh_relativistic_2000`
        for the relativistic free-free Gaunt factor. This is given by Eq. 4 of
        :cite:t:`itoh_relativistic_2000`,

        .. math::

            g_{ff} = \sum_{i,j=0}^{10} a_{ij}t^iU^j,

        where :math:`t=(\log{T} - 7.25)/1.25` and :math:`U=(\log{u} + 1.5)/2.5`.

        Parameters
        ----------
        temperature : `~astropy.units.Quantity`
            The temperature(s) for which to calculate the Gaunt factor
        wavelength : `~astropy.units.Quantity`
            The wavelength(s) at which to calculate the Gaunt factor
        atomic_number : `int`
            The atomic number of the emitting element
        charge_state : `int`
            The charge state of the emitting ion

        See Also
        --------
        fiasco.Ion.free_free
        """
        gf_itoh = self._free_free_itoh(temperature, wavelength, atomic_number)
        gf_sutherland = self._free_free_sutherland(temperature, wavelength, charge_state)
        gf = np.where(np.isnan(gf_itoh), gf_sutherland, gf_itoh)
        return gf

    @needs_dataset('itoh')
    @u.quantity_input
    def _free_free_itoh(self, temperature: u.K, wavelength: u.angstrom, atomic_number) -> u.dimensionless_unscaled:
        log10_temperature = np.log10(temperature.to(u.K).value)
        # calculate scaled energy and temperature
        tmp = np.outer(temperature, wavelength)
        lower_u = const.h * const.c / const.k_B / tmp
        upper_u = 1. / 2.5 * (np.log10(lower_u) + 1.5)
        t = 1. / 1.25 * (log10_temperature - 7.25)
        itoh_coefficients = self._itoh['a'][self._itoh['Z']==atomic_number].squeeze()
        # calculate Gaunt factor
        gf = u.Quantity(np.zeros(upper_u.shape))
        for i in range(itoh_coefficients.shape[0]):
            for j in range(itoh_coefficients.shape[1]):
                gf += (itoh_coefficients[i, j] * (t**i))[:, np.newaxis] * (upper_u**j)
        # apply NaNs where Itoh approximation is not valid
        gf = np.where(np.logical_and(np.log10(lower_u) >= -4., np.log10(lower_u) <= 1.0), gf, np.nan)
        gf[np.where(np.logical_or(log10_temperature <= 6.0, log10_temperature >= 8.5)), :] = np.nan
        return gf

    @needs_dataset('gffgu')
    @u.quantity_input
    def _free_free_sutherland(self, temperature: u.K, wavelength: u.angstrom, charge_state) -> u.dimensionless_unscaled:
        Ry = const.h * const.c * const.Ryd
        tmp = np.outer(temperature, wavelength)
        lower_u = const.h * const.c / const.k_B / tmp
        gamma_squared = ((charge_state**2) * Ry / const.k_B / temperature[:, np.newaxis]
                         * np.ones(lower_u.shape))
        # NOTE: This escape hatch avoids a divide-by-zero warning as we cannot take log10
        # of 0. This does not matter as the free-free continuum will be 0 for zero charge
        # state anyway.
        if charge_state == 0:
            return u.Quantity(np.zeros(lower_u.shape))
        # convert to index coordinates
        i_lower_u = (np.log10(lower_u) + 4.) * 10.
        i_gamma_squared = (np.log10(gamma_squared) + 4.) * 5.
        indices = [i_gamma_squared.flatten(), i_lower_u.flatten()]
        # interpolate data to scaled quantities
        # FIXME: interpolate without reshaping?
        gf_data = self._gffgu['gaunt_factor'].reshape(
            np.unique(self._gffgu['u']).shape[0],
            np.unique(self._gffgu['gamma_squared']).shape[0],
        )
        gf = map_coordinates(gf_data, indices, order=1, mode='nearest').reshape(lower_u.shape)
        return u.Quantity(np.where(gf < 0., 0., gf))

    @u.quantity_input
    def free_free_integrated(self, temperature: u.K, charge_state, use_itoh=False) -> u.dimensionless_unscaled:
        r"""
        The wavelength-integrated Gaunt factor for free-free emission as a function of temperature.

        The wavelength-integrated Gaunt factor is primarily used for calculating the total radiative losses from
        free-free emission.
        By default, this calculation is done with the form specified in :cite:t:`sutherland_accurate_1998`,
        which is valid over a wide range of temperatures.
        The ``use_itoh`` option substitutes the form specified by :cite:t:`itoh_radiative_2002`, which is more
        accurate but has a more limited range of validity.
        The difference between the two forms is small, as shown in :cite:t:`young_chianti_2019-1`.
        The CHIANTI atomic database only uses the :cite:t:`sutherland_accurate_1998` form as a result, but
        includes the data sets for both forms.

        .. note:: The Gaunt factor calculation of :cite:t:`itoh_radiative_2002` includes both a relativistic
                  (Eq. 5) and non-relativistic (Eq. 13) form.
                  The relativistic form is valid over the temperature range :math:`6.0\leq\log_{10}T\leq8.5`
                  and for charge states :math:`1\le z\le 28`.
                  The nonrelativistic form is valid over :math:`-3\leq\log_{10}\gamma^{2}\leq 2` where
                  :math:`\gamma^2=z^2\mathrm{Ry}/k_BT`.
                  Outside of these ranges, the form of :cite:t:`sutherland_accurate_1998` is used.

        Parameters
        ----------
        temperature : `~astropy.units.Quantity`
            The temperature(s) for which to calculate the Gaunt factor
        charge_state : `int`
            The charge state of the ion
        use_itoh : `bool`, optional
            If true, use the :cite:t:`itoh_radiative_2002` Gaunt factors over valid ranges.
            If false (default), use the :cite:t:`sutherland_accurate_1998` Gaunt factors instead.

        See Also
        --------
        fiasco.Ion.free_free_radiative_loss
        """
        if charge_state == 0:
            return u.Quantity(np.zeros(temperature.shape))
        gf = self._free_free_sutherland_integrated(temperature, charge_state)
        if use_itoh:
            gf_itoh = self._free_free_itoh_integrated(temperature, charge_state)
            gf = np.where(np.isnan(gf_itoh), gf, gf_itoh)
        return gf

    @needs_dataset('gffint')
    @u.quantity_input
    def _free_free_sutherland_integrated(self, temperature: u.K, charge_state) -> u.dimensionless_unscaled:
        """
        The wavelength-integrated free-free Gaunt factor, as specified by :cite:t:`sutherland_accurate_1998`,
        in Section 2.4 of that work.

        Parameters
        ----------
        temperature : `~astropy.units.Quantity`
            The temperature(s) for which to calculate the Gaunt factor
        charge_state : `int`
            The charge state of the ion
        """
        temperature = np.atleast_1d(temperature)
        Ry = const.h * const.c * const.Ryd
        log_gamma_squared = np.log10((charge_state**2 * Ry) / (const.k_B * temperature))
        index = [np.abs(self._gffint['log_gamma_squared'] - x).argmin() for x in log_gamma_squared]
        delta = log_gamma_squared - self._gffint['log_gamma_squared'][index]
        # The spline fit was pre-calculated by Sutherland 1998:
        return self._gffint['gaunt_factor'][index] + delta * (self._gffint['s1'][index] + delta * (self._gffint['s2'][index] + delta * self._gffint['s3'][index]))

    @u.quantity_input
    def _free_free_itoh_integrated(self, temperature: u.K, charge_state) -> u.dimensionless_unscaled:
        r"""
        The wavelength-integrated free-free Gaunt factor, as specified by :cite:t:`itoh_radiative_2002`.

        Parameters
        ----------
        temperature : `~astropy.units.Quantity`
            The temperature(s) for which to calculate the Gaunt factor
        charge_state : `int`
            The charge state of the ion
        """
        temperature = np.atleast_1d(temperature)
        try:
            gf_relativistic = self._free_free_itoh_integrated_relativistic(temperature, charge_state)
        except MissingDatasetException:
            gf_relativistic = np.full(len(temperature), np.nan)
        try:
            gf_nonrelativistic = self._free_free_itoh_integrated_nonrelativistic(temperature, charge_state)
        except MissingDatasetException:
            gf_nonrelativistic = np.full(len(temperature), np.nan)
        return np.where(np.isnan(gf_nonrelativistic), gf_relativistic, gf_nonrelativistic)

    @needs_dataset('itoh_integrated_gaunt_nonrel')
    @u.quantity_input
    def _free_free_itoh_integrated_nonrelativistic(self, temperature: u.K, charge_state) -> u.dimensionless_unscaled:
        r"""
        The wavelength-integrated non-relativistic free-free Gaunt factor, as specified by
        :cite:t:`itoh_radiative_2002`.

        This form is only valid between :math:`-3.0 \leq \log_{10} \gamma^{2} \leq 2.0`.  We use the form
        specified by :cite:t:`sutherland_accurate_1998` outside of this range.

        Parameters
        ----------
        temperature : `~astropy.units.Quantity`
            The temperature(s) for which to calculate the Gaunt factor
        charge_state : `int`
            The charge state of the ion
        """
        temperature = np.atleast_1d(temperature)
        Ry = const.h * const.c * const.Ryd
        gamma_squared = (charge_state**2) * Ry / (const.k_B * temperature)
        gf = u.Quantity(np.zeros(temperature.shape))
        Gamma = (np.log10(gamma_squared) + 0.5) / 2.5
        b_array = self._itoh_integrated_gaunt_nonrel['b_i']
        for i in range(b_array.shape[0]):
            gf += b_array[i] * Gamma**i
        not_valid = np.logical_or(np.log10(gamma_squared) < -3.0, np.log10(gamma_squared) > 2.0)
        return np.where(not_valid, np.nan, gf)

    @needs_dataset('itoh_integrated_gaunt')
    @u.quantity_input
    def _free_free_itoh_integrated_relativistic(self, temperature: u.K, charge_state) -> u.dimensionless_unscaled:
        r"""
        The wavelength-integrated relativistic free-free Gaunt factor, as specified by
        :cite:t:`itoh_radiative_2002`.

        The relativistic approximation is only valid between :math:`6.0 \leq \log_{10} T_{e} \leq 8.5`, and charges between 1 and 28.
        At cooler temperatures, the calculation uses the non-relativistic form, while at higher temperatures it defaults back to the
        expressions from :cite:t:`sutherland_accurate_1998`.

        Parameters
        ----------
        temperature : `~astropy.units.Quantity`
            The temperature(s) for which to calculate the Gaunt factor
        charge_state : `int`
            The charge state of the ion
        """
        temperature = np.atleast_1d(temperature)
        log_temperature = np.log10(temperature.to_value('K'))
        z = (charge_state - 14.5) / 13.5
        t = (log_temperature-7.25)/1.25
        gf = u.Quantity(np.zeros(temperature.shape))
        a_matrix = self._itoh_integrated_gaunt['a_ik']
        for i in range(a_matrix.shape[0]):
            for k in range(a_matrix.shape[1]):
                gf += a_matrix[i,k] * z**i * t**k
        not_valid = np.logical_or(log_temperature < 6.0, log_temperature > 8.5)
        return np.where(not_valid, np.nan, gf)

    @needs_dataset('klgfb')
    @u.quantity_input
    def free_bound(self, E_scaled, n, l) -> u.dimensionless_unscaled:
        r"""
        The Gaunt factor for free-bound emission as a function of scaled energy.

        The empirical fits are taken from Table 1 of :cite:t:`karzas_electron_1961`.
        In CHIANTI, this is used to compute the cross-sections in the free-bound continuum.

        Parameters
        ----------
        E_scaled : `float`
            Rratio of photon energy to the ionization energy.
        n : `int`
            The principal quantum number
        l : `int`
            The azimuthal quantum number

        See Also
        --------
        fiasco.Ion.free_bound
        """
        E_scaled = np.atleast_1d(E_scaled)
        index_nl = np.where(np.logical_and(self._klgfb['n'] == n, self._klgfb['l'] == l))[0]
        # If there is no Gaunt factor for n, l, set it to 1
        if index_nl.shape == (0,):
            gf = np.ones(E_scaled.shape)
        else:
            gf_interp = splrep(self._klgfb['log_pe'][index_nl, :].squeeze(),
                               self._klgfb['log_gaunt_factor'][index_nl, :].squeeze())
            gf = np.exp(splev(E_scaled, gf_interp))

        return gf


    @u.quantity_input
    def free_bound_integrated(self, temperature: u.K, atomic_number, charge_state, n_0,
                            ionization_potential: u.eV, ground_state=True) -> u.dimensionless_unscaled:
        r"""
        The wavelength-integrated Gaunt factor for free-bound emission as a function of temperature.

        The wavelength-integrated free-bound Gaunt factor is calculated using the approach of
        :cite:t:`mewe_calculated_1986`.
        The Gaunt factor is not calculated for individual levels, except that the ground state has
        been specified to be :math:`g_{fb}(n_{0}) = 0.9` following :cite:t:`mewe_calculated_1986`.
        For more details on this calculation, see :ref:`fiasco-topic-guide-freebound-gaunt-factor`.

        Parameters
        ----------
        temperature : `~astropy.units.Quantity`
            The temperature(s) for which to calculate the Gaunt factor
        atomic_number : `int`
            The atomic number of the element
        charge_state : `int`
            The charge state of the ion
        n_0 : `int`
            The principal quantum number n of the ground state of the recombined ion
        ionization_potential : `~astropy.units.Quantity`
            The ionization potential of the recombined ion
        ground_state : `bool`, optional
            If True (default), calculate the Gaunt factor for recombination onto the ground state :math:`n = 0`.
            Otherwise, calculate for recombination onto higher levels with :math:`n > 1`.  See Equation 16 of
            :cite:t:`mewe_calculated_1986`.

        See Also
        --------
        fiasco.Ion.free_bound_radiative_loss
        """
        z = charge_state
        if z == 0:
            return u.Quantity(np.zeros(temperature.shape))
        else:
            Ry = const.h * const.c * const.Ryd
            prefactor = (Ry / (const.k_B * temperature)).decompose()
            if ground_state:
                g_n_0 = 0.9 # The ground state Gaunt factor, g_fb(n_0), prescribed by Mewe et al 1986
                z_0 = n_0 * np.sqrt(ionization_potential / Ry).decompose()
                # NOTE: Low temperatures can lead to very large terms in the exponential and in some
                # cases, the exponential term can exceed the maximum number expressible in double precision.
                # These terms are eventually set to zero anyway since the ionization fraction is so small
                # at these temperatures.
                with np.errstate(over='ignore'):
                    f_2 = g_n_0 * self._zeta_0(atomic_number, charge_state) * (z_0**4 / n_0**5) * np.exp((prefactor * z_0**2)/(n_0**2))
            else:
                f_1 = -0.5 * polygamma(2, n_0 + 1)
                with np.errstate(over='ignore'):
                    f_2 = 2.0 * f_1 * z**4 * np.exp((prefactor * z**2)/((n_0+1)**2))
            # Large terms in the exponential can lead to infs in the f_2 term. These will be zero'd
            # out when multiplied by the ionization equilibrium anyway
            f_2 = np.where(np.isinf(f_2), 0.0, f_2)

            return (prefactor * f_2)

    def _zeta_0(self, atomic_number, charge_state) -> u.dimensionless_unscaled:
        r"""
        :math:`\zeta_{0}`, the number of vacancies in an ion, which is used to calculate
        the wavelength-integrated free-bound Gaunt factor of that ion.

        See Section 2.2 and Table 1 of :cite:t:`mewe_calculated_1986`.
        """
        difference = atomic_number - (charge_state + 1)
        if difference <= 0:
            max_vacancies = 1
        elif difference <= 8 and difference > 0:
            max_vacancies = 9
        elif difference <= 22 and difference > 8:
            max_vacancies = 27
        else:
            max_vacancies = 55
        return max_vacancies - difference
