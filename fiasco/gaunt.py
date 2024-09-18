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

__all__ = ['GauntFactor']

class GauntFactor:
    """
    Class for calculating the Gaunt factor for various continuum processes.

    The Gaunt factor is defined as the ratio of the true cross-section to the
    semi-classical Kramers cross-section, and thus is essentially a multiplicative
    correction for quantum mechanical effects.  It is a unitless quantity.
    """
    def __init__(self, hdf5_dbase_root=None, *args, **kwargs):
        if hdf5_dbase_root is None:
            self.hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
        else:
            self.hdf5_dbase_root = hdf5_dbase_root
        check_database(self.hdf5_dbase_root, **kwargs)

    def __str__(self):
        return "Gaunt factor object"

    def __repr__(self):
        return "Gaunt factor object"

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
    def _itohintrel(self):
        data_path = '/'.join(['continuum', 'itoh_integrated_gaunt'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _itohintnonrel(self):
        data_path = '/'.join(['continuum', 'itoh_integrated_gaunt_nonrel'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _klgfb(self):
        data_path = '/'.join(['continuum', 'klgfb'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @u.quantity_input
    def free_free(self, temperature: u.K, atomic_number, charge_state, wavelength: u.angstrom) -> u.dimensionless_unscaled:
        r"""
        Free-free Gaunt factor as a function of wavelength.

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
        atomic_number : `int`
            The atomic number of the emitting element
        charge_state : `int`,
            The charge state of the emitting ion
        wavelength : `~astropy.units.Quantity`
            The wavelength(s) at which to calculate the Gaunt factor
        """
        gf_itoh = self._free_free_itoh(temperature, atomic_number, wavelength)
        gf_sutherland = self._free_free_sutherland(temperature, charge_state, wavelength)
        gf = np.where(np.isnan(gf_itoh), gf_sutherland, gf_itoh)

        return gf

    @u.quantity_input
    def _free_free_itoh(self, temperature: u.K, atomic_number, wavelength: u.angstrom) -> u.dimensionless_unscaled:
        log10_temperature = np.log10(temperature.to(u.K).value)
        # calculate scaled energy and temperature
        tmp = np.outer(temperature, wavelength)
        lower_u = const.h * const.c / const.k_B / tmp
        upper_u = 1. / 2.5 * (np.log10(lower_u) + 1.5)
        t = 1. / 1.25 * (log10_temperature - 7.25)
        itoh_coefficients = self._itoh['a'][atomic_number-1]
        # calculate Gaunt factor
        gf = u.Quantity(np.zeros(upper_u.shape))
        for j in range(11):
            for i in range(11):
                gf += (itoh_coefficients[i, j] * (t**i))[:, np.newaxis] * (upper_u**j)
        # apply NaNs where Itoh approximation is not valid
        gf = np.where(np.logical_and(np.log10(lower_u) >= -4., np.log10(lower_u) <= 1.0),
                      gf, np.nan)
        gf[np.where(np.logical_or(log10_temperature <= 6.0, log10_temperature >= 8.5)), :] = np.nan

        return gf

    @needs_dataset('gffgu')
    @u.quantity_input
    def _free_free_sutherland(self, temperature: u.K, charge_state, wavelength: u.angstrom) -> u.dimensionless_unscaled:
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
    def free_free_total(self, temperature: u.K, charge_state, itoh=False, relativistic=True) -> u.dimensionless_unscaled:
        """
        The total (wavelength-averaged) free-free Gaunt factor, used for calculating
        the total radiative losses from free-free emission.

        Parameters
        ----------
        temperature : `~astropy.units.Quantity`
            The temperature(s) for which to calculate the Gaunt factor
        charge_state : `int`,
            The charge state of the ion
        itoh : `bool`, optional
            Specify whether to use the approximations specified by :cite:t:`itoh_radiative_2002`.  
            If true, use the forms by :cite:t:`itoh_radiative_2002`.  If false, use the forms by
            :cite:t:`sutherland_accurate_1998`.
        relativistic : `bool`, optional
            If using the :cite:t:`itoh_radiative_2002` approximations, use the relativistic form
            instead of the non-relativistic form.
        """
        if charge_state == 0:
            return u.Quantity(np.zeros(temperature.shape))
        else:
            if not itoh:
                return self._free_free_sutherland_integrated(temperature, charge_state)
            else:
                if relativistic:
                    return self._free_free_itoh_integrated_relativistic(temperature, charge_state)
                else:
                    return self._free_free_itoh_integrated_nonrelativistic(temperature, charge_state)

    @needs_dataset('gffint')
    @u.quantity_input
    def _free_free_sutherland_integrated(self, temperature: u.K, charge_state) -> u.dimensionless_unscaled:
        """
        The total (wavelength-averaged) free-free Gaunt factor, as specified by :cite:t:`sutherland_accurate_1998`,
        in Section 2.4 of that work.
        
        This is the default option used by CHIANTI for integrated free-free Gaunt factor, and we also use it 
        where the other versions are not valid.
        
        Parameters
        ----------
        temperature : `~astropy.units.Quantity`
            The temperature(s) for which to calculate the Gaunt factor
        charge_state : `int`,
            The charge state of the ion
        """
        temperature = np.atleast_1d(temperature)
        Ry = const.h * const.c * const.Ryd
        log_gamma_squared = np.log10((charge_state**2 * Ry) / (const.k_B * temperature))
        index = [np.abs(self._gffint['log_gamma_squared'] - x).argmin() for x in log_gamma_squared]
        delta = log_gamma_squared - self._gffint['log_gamma_squared'][index]
        # The spline fit was pre-calculated by Sutherland 1998:
        return self._gffint['gaunt_factor'][index] + delta * (self._gffint['s1'][index] + delta * (self._gffint['s2'][index] + delta * self._gffint['s3'][index]))


    @needs_dataset('itohintnonrel')
    @u.quantity_input
    def _free_free_itoh_integrated_nonrelativistic(self, temperature: u.K, charge_state) -> u.dimensionless_unscaled:
        """
        The total (wavelength-averaged) non-relativistic free-free Gaunt factor, as specified by
        :cite:t:`itoh_radiative_2002`.
        
        This form is only valid between :math:`-3.0 \leq \log_{10} \gamma^{2} \leq 2.0`.  We use the form
        specified by :cite:t:`sutherland_accurate_1998` outside of this range.  
        
        Parameters
        ----------
        temperature : `~astropy.units.Quantity`
            The temperature(s) for which to calculate the Gaunt factor
        charge_state : `int`,
            The charge state of the ion
        """
        temperature = np.atleast_1d(temperature)
        Ry = const.h * const.c * const.Ryd
        gamma_squared = (charge_state**2) * Ry / (const.k_B * temperature)
        summation = u.Quantity(np.zeros(temperature.shape))
        Gamma = (np.log10(gamma_squared) + 0.5) / 2.5
        for j in range(len(summation)):
            if np.log10(gamma_squared[j]) < -3.0 or np.log10(gamma_squared[j]) > 2.0:
                summation[j] = self._free_free_sutherland_integrated(temperature[j], charge_state)
            else:
                for i in range(len(self._itohintnonrel['b_i'])):
                    summation[j] += self._itohintnonrel['b_i'][i] * Gamma[j]**i
        return summation

    @needs_dataset('itohintrel')
    @u.quantity_input
    def _free_free_itoh_integrated_relativistic(self, temperature: u.K, charge_state) -> u.dimensionless_unscaled:
        """
        The total (wavelength-averaged) relativistic free-free Gaunt factor, as specified by
        :cite:t:`itoh_radiative_2002`.

        The relativistic approximation is only valid between :math:`6.0 \leq \log_{10} T_{e} \leq 8.5`, and charges between 1 and 28.
        At cooler temperatures, the calculation uses the non-relativistic form, while at higher temperatures it defaults back to the 
        expressions from :cite:t:`sutherland_accurate_1998`.

        Parameters
        ----------
        temperature : `~astropy.units.Quantity`
            The temperature(s) for which to calculate the Gaunt factor
        charge_state : `int`,
            The charge state of the ion
        """
        temperature = np.atleast_1d(temperature)
        z = (charge_state - 14.5) / 13.5
        t = (np.log10(temperature.data)-7.25)/1.25
        summation = u.Quantity(np.zeros(temperature.shape))
        for j in range(len(summation)):
            if np.log10(temperature[j].data) < 6.0:
                summation[j] = self._free_free_itoh_integrated_nonrelativistic(temperature[j], charge_state)
            elif np.log10(temperature[j].data) > 8.5:
                summation[j] = self._free_free_sutherland_integrated(temperature[j], charge_state)
            else:
                for i in range(len(self._itohintrel['a_ik'][:][0])):
                    for k in range(len(self._itohintrel['a_ik'][0][:])):
                        summation[j] += self._itohintrel['a_ik'][i][k] * z**i * t[j]**k
        return summation

    @needs_dataset('klgfb')
    @u.quantity_input
    def free_bound(self, E_scaled: u.dimensionless_unscaled, n, l) -> u.dimensionless_unscaled:
        r"""
        Free-bound Gaunt factor as a function of scaled energy.

        The empirical fits are taken from Table 1 of :cite:t:`karzas_electron_1961`.

        Parameters
        ----------
        E_scaled : `~astropy.units.Quantity`
            A scaled energy, the ratio of photon energy divided by ionization energy.
        n : `int`
            The principal quantum number
        l : `int`
            The azimuthal quantum number
        """
        index_nl = np.where(np.logical_and(self._klgfb['n'] == n, self._klgfb['l'] == l))[0]
        # If there is no Gaunt factor for n, l, set it to 1
        if index_nl.shape == (0,):
            gf = 1
        else:
            gf_interp = splrep(self._klgfb['log_pe'][index_nl, :].squeeze(),
                               self._klgfb['log_gaunt_factor'][index_nl, :].squeeze())
            gf = np.exp(splev(E_scaled, gf_interp))

        return gf


    @u.quantity_input
    def free_bound_total(self, temperature: u.K, atomic_number, charge_state, n_0,
                            ionization_potential: u.eV, ground_state=True) -> u.dimensionless_unscaled:
        r"""
        The total Gaunt factor for free-bound emission, using the expressions from :cite:t:`mewe_calculated_1986`.

        The Gaunt factor is not calculated for individual levels, except that the ground state has
        been specified to be :math: `g_{fb}(n_{0}) = 0.9` following :cite:t:`mewe_calculated_1986`.

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
        ionization_potential :
            The ionization potential of the recombined ion
        ground_state : `bool`, optional
            If True (default), calculate the Gaunt factor for recombination onto the ground state :math: `n = 0`.
            Otherwise, calculate for recombination onto higher levels with :math: `n > 1`.  See Equation 16 of
            :cite:t:`mewe_calculated_1986`.

        Notes
        -----
        Equation 14 of :cite:t:`mewe_calculated_1986` has a simple
        analytic solution.  They approximate
        .. math::
            f_{1}(Z, n, n_{0} ) = \sum_{1}^{\infty} n^{-3} - \sum_{1}^{n_{0}} n^{-3} = \zeta(3) - \sum_{1}^{n_{0}} n^{-3} \approx 0.21 n_{0}^{-1.5}

        where :math: `\zeta(x)` is the Riemann zeta function.

        However, the second sum is analytic, :math: `\sum_{1}^{n_{0}} n^{-3} = \zeta(3) + \frac{1}{2}\psi^{(2)}(n_{0}+1)`
        where :math: `\psi^{n}(x)` is the n-th derivative of the digamma function (a.k.a. the polygamma function).
        So, we can write the full solution as:
        .. math::
            f_{1}(Z, n, n_{0}) = \zeta(3) - \sum_{1}^{n_{0}} n^{-3} = - \frac{1}{2}\psi^{(2)}(n_{0}+1)

        The final expression is therefore simplified and more accurate than :cite:t:`mewe_calculated_1986`.
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
