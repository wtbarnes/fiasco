"""
Ion object. Holds all methods and properties of a CHIANTI ion.
"""
import numpy as np
from scipy.interpolate import splrep, splev, interp1d
from scipy.ndimage import map_coordinates
import astropy.units as u
import astropy.constants as const

from .base import IonBase, ContinuumBase
from .collections import IonCollection
from .level import Level, Transitions
from fiasco import proton_electron_ratio
from fiasco.util import (needs_dataset, vectorize_where, vectorize_where_sum,
                         burgess_tully_descale)

__all__ = ['Ion']


class Ion(IonBase, ContinuumBase):
    """
    Ion class

    The ion object is the fundamental unit of the fiasco library. An Ion object contains
    all of the properties and methods needed to access important information about each ion
    from the CHIANTI database.

    Parameters
    ----------
    ion_name : `str`
    temperature : `~astropy.units.Quantity`

    Other Parameters
    ----------------
    ioneq_filename : `str`, optional
        Ionization equilibrium dataset
    abundance_filename : `str`, optional
        Abundance dataset
    ip_filename : `str`, optional
        Ionization potential dataset

    Examples
    --------
    """

    @u.quantity_input
    def __init__(self, ion_name, temperature: u.K, *args, **kwargs):
        super().__init__(ion_name, *args, **kwargs)
        self.temperature = np.atleast_1d(temperature)
        # Get selected datasets
        # TODO: do not hardcode defaults, pull from rc file
        self._dset_names = {}
        self._dset_names['ioneq_filename'] = kwargs.get('ioneq_filename', 'chianti')
        self._dset_names['abundance_filename'] = kwargs.get('abundance_filename',
                                                            'sun_photospheric_1998_grevesse')
        self._dset_names['ip_filename'] = kwargs.get('ip_filename', 'chianti')

    def __repr__(self):
        n_levels = self._elvlc['level'].shape[0] if self._elvlc else 0
        n_transitions = self._wgfa['lower_level'].shape[0] if self._wgfa else 0
        return f"""CHIANTI Database Ion
---------------------
Name: {self.ion_name}
Element: {self.element_name} ({self.atomic_number})
Charge: +{self.charge_state}
Number of Levels: {n_levels}
Number of Transitions: {n_transitions}

Temperature range: [{self.temperature[0].to(u.MK)}, {self.temperature[-1].to(u.MK)}]

HDF5 Database: {self.hdf5_dbase_root}
Using Datasets:
  ioneq: {self._dset_names['ioneq_filename']}
  abundance: {self._dset_names['abundance_filename']}
  ip: {self._dset_names['ip_filename']}"""

    def __getitem__(self, key):
        if self._elvlc is None:
            raise IndexError(f'No energy levels available for {self.ion_name}')
        # Throw an index error to stop iteration
        _ = self._elvlc['level'][key]
        return Level(key, self._elvlc)

    def __add__(self, value):
        return IonCollection(self, value)

    def __radd__(self, value):
        return IonCollection(value, self)

    @property
    @needs_dataset('elvlc', 'wgfa')
    def transitions(self):
        return Transitions(self._elvlc, self._wgfa)

    @property
    def ioneq(self):
        """
        Ionization equilibrium data interpolated to the given temperature

        Interpolated the pre-computed ionization fractions stored in CHIANTI to the temperature
        of the ion. Returns NaN where interpolation is out of range of the data. For computing
        ionization equilibrium outside of this temperature range, it is better to use the ionization
        and recombination rates.

        See Also
        --------
        fiasco.Element.equilibrium_ionization
        """
        f = interp1d(self._ioneq[self._dset_names['ioneq_filename']]['temperature'],
                     self._ioneq[self._dset_names['ioneq_filename']]['ionization_fraction'],
                     kind='linear',
                     bounds_error=False,
                     fill_value=np.nan)
        ioneq = f(self.temperature)
        isfinite = np.isfinite(ioneq)
        ioneq[isfinite] = np.where(ioneq[isfinite] < 0., 0., ioneq[isfinite])
        return u.Quantity(ioneq)

    @property
    def abundance(self):
        """
        Elemental abundance relative to H
        """
        return self._abundance[self._dset_names['abundance_filename']]

    @property
    @needs_dataset('ip')
    @u.quantity_input
    def ip(self) -> u.erg:
        """
        Ionization potential with reasonable units
        """
        return self._ip[self._dset_names['ip_filename']] * const.h * const.c

    @property
    def hydrogenic(self):
        """
        Is the ion hydrogen-like or not
        """
        return (self.atomic_number - self.charge_state == 1) and (self.atomic_number >= 6)

    @property
    def helium_like(self):
        """
        Is the ion helium like or not
        """
        return (self.atomic_number - self.charge_state == 2) and (self.atomic_number >= 10)

    @needs_dataset('scups')
    @u.quantity_input
    def effective_collision_strength(self) -> u.dimensionless_unscaled:
        """
        Maxwellian-averaged collision strength, typically denoted by :math:`\\Upsilon`

        See Also
        --------
        fiasco.util.burgess_tully_descale : Descale and interpolate :math:`\\Upsilon`
        """
        kBTE = np.outer(const.k_B*self.temperature, 1.0/self._scups['delta_energy'])
        upsilon = burgess_tully_descale(self._scups['bt_t'],
                                        self._scups['bt_upsilon'],
                                        kBTE.T,
                                        self._scups['bt_c'],
                                        self._scups['bt_type'])
        upsilon = u.Quantity(np.where(upsilon > 0., upsilon, 0.))
        return upsilon.T

    @needs_dataset('elvlc', 'scups')
    @u.quantity_input
    def electron_collision_deexcitation_rate(self) -> u.cm**3 / u.s:
        """
        Collisional de-excitation rate coefficient for electrons.

        According to Eq. (4.12) of [1]_, the rate coefficient for collisional de-excitation
        is given by,

        .. math::

           C^d_{ji} = I_Ha_0^2\sqrt{\\frac{8\pi}{mk_B}}\\frac{\\Upsilon}{\omega_jT^{1/2}},

        where :math:`j,i` are the upper and lower level indices, respectively, :math:`I_H` is the
        ionization potential for H, :math:`a_0` is the Bohr radius, :math:`\\Upsilon` is the
        effective collision strength, and :math:`\omega_j` is the statistical weight of the
        level :math:`j`.

        References
        ----------
        .. [1] Phillips, K., et al., 2008, `Ultraviolet and X-ray Spectroscopy of the Solar Atmosphere <http://adsabs.harvard.edu/abs/2008uxss.book.....P>`_

        See Also
        --------
        electron_collision_excitation_rate : Excitation rate due to collisions
        effective_collision_strength : Maxwellian-averaged collision strength, :math:`\\Upsilon`
        """
        c = (const.h**2)/((2. * np.pi * const.m_e)**(1.5) * np.sqrt(const.k_B))
        upsilon = self.effective_collision_strength()
        omega_upper = 2. * self._elvlc['J'][self._scups['upper_level'] - 1] + 1.
        return c * upsilon / np.sqrt(self.temperature[:, np.newaxis]) / omega_upper

    @needs_dataset('elvlc', 'scups')
    @u.quantity_input
    def electron_collision_excitation_rate(self, deexcitation_rate=None) -> u.cm**3 / u.s:
        """
        Collisional excitation rate coefficient for electrons.

        The rate coefficient for collisional excitation is given by,

        .. math::

            C^e_{ij} = \\frac{\omega_j}{\omega_i}C^d_{ji}e^{-k_BT_e/\Delta E_{ij}}

        where :math:`j,i` are the upper and lower level indices, respectively, :math:`\omega_j,\omega_i` 
        are the statistical weights of the upper and lower levels, respectively, and :math:`\Delta E_{ij}`
        is the energy of the transition.

        Parameters
        ----------
        deexcitation_rate : `~astropy.units.Quantity`, optional
            Optionally specify deexcitation rate to speedup calculation

        See Also
        --------
        electron_collision_deexcitation_rate : De-excitation rate due to collisions
        """
        if deexcitation_rate is None:
            deexcitation_rate = self.electron_collision_deexcitation_rate()
        omega_upper = 2. * self._elvlc['J'][self._scups['upper_level'] - 1] + 1.
        omega_lower = 2. * self._elvlc['J'][self._scups['lower_level'] - 1] + 1.
        kBTE = np.outer(1./const.k_B/self.temperature, self._scups['delta_energy'])
        return omega_upper / omega_lower * deexcitation_rate * np.exp(-kBTE)

    @needs_dataset('psplups')
    @u.quantity_input
    def proton_collision_excitation_rate(self) -> u.cm**3 / u.s:
        """
        Collisional excitation rate coefficient for protons
        """
        # Create scaled temperature--these are not stored in the file
        bt_t = [np.linspace(0, 1, ups.shape[0]) for ups in self._psplups['bt_rate']]
        # Get excitation rates directly from scaled data
        kBTE = np.outer(const.k_B*self.temperature, 1.0/self._psplups['delta_energy'])
        ex_rate = burgess_tully_descale(bt_t,
                                        self._psplups['bt_rate'],
                                        kBTE.T,
                                        self._psplups['bt_c'],
                                        self._psplups['bt_type'])
        return u.Quantity(np.where(ex_rate > 0., ex_rate, 0.), u.cm**3/u.s).T

    @needs_dataset('elvlc', 'psplups')
    @u.quantity_input
    def proton_collision_deexcitation_rate(self, excitation_rate=None) -> u.cm**3 / u.s:
        """
        Collisional de-excitation rate coefficient for protons
        """
        kBTE = np.outer(const.k_B*self.temperature, 1.0/self._psplups['delta_energy'])
        if excitation_rate is None:
            excitation_rate = self.proton_collision_excitation_rate()
        omega_upper = 2.*self._elvlc['J'][self._psplups['upper_level'] - 1] + 1.
        omega_lower = 2.*self._elvlc['J'][self._psplups['lower_level'] - 1] + 1.
        dex_rate = (omega_lower / omega_upper) * excitation_rate * np.exp(1. / kBTE)

        return dex_rate

    @needs_dataset('elvlc', 'wgfa', 'scups')
    @u.quantity_input
    def level_populations(self,
                          density: u.cm**(-3),
                          include_protons=True) -> u.dimensionless_unscaled:
        """
        Compute energy level populations as a function of temperature and density

        Parameters
        ----------
        density : `~astropy.units.Quantity`
        include_protons : `bool`, optional
            If True (default), include proton excitation and de-excitation rates

        Returns
        -------
        `~astropy.units.Quantity`
            A ``(l, m, n)`` shaped quantity, where ``l`` is the number of
            temperatures, ``m`` is the number of densities, and ``n``
            is the number of energy levels.
        """
        level = self._elvlc['level']
        lower_level = self._scups['lower_level']
        upper_level = self._scups['upper_level']
        coeff_matrix = np.zeros(self.temperature.shape + (level.max(), level.max(),))/u.s

        # Radiative decay out of current level
        coeff_matrix[:, level-1, level-1] -= vectorize_where_sum(
            self.transitions.upper_level, level, self.transitions.A.value) * self.transitions.A.unit
        # Radiative decay into current level from upper levels
        coeff_matrix[:, self.transitions.lower_level-1, self.transitions.upper_level-1] += (
            self.transitions.A)

        # Collisional--electrons
        dex_rate_e = self.electron_collision_deexcitation_rate()
        ex_rate_e = self.electron_collision_excitation_rate(deexcitation_rate=dex_rate_e)
        ex_diagonal_e = vectorize_where_sum(
            lower_level, level, ex_rate_e.value.T, 0).T * ex_rate_e.unit
        dex_diagonal_e = vectorize_where_sum(
            upper_level, level, dex_rate_e.value.T, 0).T * dex_rate_e.unit
        # Collisional--protons
        if include_protons and self._psplups is not None:
            lower_level_p = self._psplups['lower_level']
            upper_level_p = self._psplups['upper_level']
            pe_ratio = proton_electron_ratio(self.temperature,
                                             **self._dset_names,
                                             hdf5_dbase_root=self.hdf5_dbase_root)
            proton_density = np.outer(pe_ratio, density)[:, :, np.newaxis]
            ex_rate_p = self.proton_collision_excitation_rate()
            dex_rate_p = self.proton_collision_deexcitation_rate(excitation_rate=ex_rate_p)
            ex_diagonal_p = vectorize_where_sum(
                lower_level_p, level, ex_rate_p.value.T, 0).T * ex_rate_p.unit
            dex_diagonal_p = vectorize_where_sum(
                upper_level_p, level, dex_rate_p.value.T, 0).T * dex_rate_p.unit

        # Populate density dependent terms and solve matrix equation for each density value
        density = np.atleast_1d(density)
        populations = np.zeros(self.temperature.shape + density.shape + (level.max(),))
        for i, d in enumerate(density):
            c_matrix = coeff_matrix.copy()
            # Collisional excitation and de-excitation out of current state
            c_matrix[:, level-1, level-1] -= d*(ex_diagonal_e + dex_diagonal_e)
            # De-excitation from upper states
            c_matrix[:, lower_level-1, upper_level-1] += d*dex_rate_e
            # Excitation from lower states
            c_matrix[:, upper_level-1, lower_level-1] += d*ex_rate_e
            # Same processes as above, but for protons
            if include_protons and self._psplups is not None:
                d_p = proton_density[:, i, :]
                c_matrix[:, level-1, level-1] -= d_p*(ex_diagonal_p + dex_diagonal_p)
                c_matrix[:, lower_level_p-1, upper_level_p-1] += d_p * dex_rate_p
                c_matrix[:, upper_level_p-1, lower_level_p-1] += d_p * ex_rate_p
            # Invert matrix
            val, vec = np.linalg.eig(c_matrix.value)
            # Eigenvectors with eigenvalues closest to zero are the solutions to the homogeneous
            # system of linear equations
            # NOTE: Sometimes eigenvalues may have complex component due to numerical stability.
            # We will take only the real component as our rate matrix is purely real
            i_min = np.argmin(np.fabs(np.real(val)), axis=1)
            pop = np.take(np.real(vec), i_min, axis=2)[range(vec.shape[0]), :, range(vec.shape[0])]
            # NOTE: The eigenvectors can only be determined up to a sign so we must enforce
            # positivity
            np.fabs(pop, out=pop)
            np.divide(pop, pop.sum(axis=1)[:, np.newaxis], out=pop)
            populations[:, i, :] = pop

        return u.Quantity(populations)

    @needs_dataset('scups', 'elvlc', 'wgfa')
    @u.quantity_input
    def contribution_function(self, density: u.cm**(-3), **kwargs) -> u.cm**3 * u.erg / u.s:
        """
        Contribution function :math:`G(n,T)` for all transitions

        The contribution function for ion :math:`k` of element :math:`X` for a
        particular transition :math:`ij` is given by,

        .. math::

           G_{ij} = \\frac{n_H}{n_e}\mathrm{Ab}(X)f_{X,k}N_jA_{ij}\Delta E_{ij}\\frac{1}{n_e},

        Note that the contribution function is often defined in differing ways by different authors.
        The contribution function is defined as above in [1]_.

        The corresponding wavelengths can be retrieved with,

        .. code-block:: python

           ion.transitions.wavelength[~ion.transitions.is_twophoton]

        Parameters
        ----------
        density : `~astropy.units.Quantity`
            Electron number density

        References
        ----------
        .. [1] Young, P. et al., 2016, J. Phys. B: At. Mol. Opt. Phys., `49, 7 <http://iopscience.iop.org/article/10.1088/0953-4075/49/7/074009/meta>`_
        """
        populations = self.level_populations(density, **kwargs)
        term = np.outer(self.ioneq, 1./density) * self.abundance * 0.83
        # Exclude two-photon transitions
        upper_level = self.transitions.upper_level[~self.transitions.is_twophoton]
        wavelength = self.transitions.wavelength[~self.transitions.is_twophoton]
        A = self.transitions.A[~self.transitions.is_twophoton]
        energy = const.h * const.c / wavelength
        i_upper = vectorize_where(self._elvlc['level'], upper_level)
        g = term[:, :, np.newaxis] * populations[:, :, i_upper] * (A * energy)
        return g

    @needs_dataset('scups', 'elvlc', 'wgfa')
    @u.quantity_input
    def emissivity(self, density: u.cm**(-3), **kwargs) -> u.erg * u.cm**(-3) / u.s:
        """
        Emissivity as a function of temperature and density for all transitions

        The emissivity is given by the expression,

        .. math::

           \epsilon(n,T) = G(n,T)n^2

        Note that, like the contribution function, emissivity is often defined in
        in differing ways by different authors.

        Parameters
        ----------
        density : `~astropy.units.Quantity`
            Electron number density

        See Also
        --------
        contribution_function : Calculate contribution function, :math:`G(n,T)`
        """
        g = self.contribution_function(density, **kwargs)
        return g * (density**2)[np.newaxis, :, np.newaxis]

    @needs_dataset('scups', 'elvlc', 'wgfa')
    @u.quantity_input
    def intensity(self,
                  density: u.cm**(-3),
                  emission_measure: u.cm**(-5), **kwargs) -> u.erg / u.cm**2 / u.s / u.steradian:
        """
        Line-of-sight intensity computed assuming a particular column emission measure

        The intensity along the line-of-sight can be written as,

        .. math::

           I = \\frac{1}{4\pi}\int\mathrm{d}T,G(n,T)n^2\\frac{dh}{dT}

        which, in the isothermal approximation, can be simplified to,

        .. math::

           I(T_0) \\approx \\frac{1}{4\pi}G(n,T_0)\mathrm{EM}(T_0)

        where,

        .. math::

           \mathrm{EM}(T) = \int\mathrm{d}h\,n^2

        is the column emission measure.

        Parameters
        ----------
        density : `~astropy.units.Quantity`
            Electron number density
        emission_measure : `~astropy.units.Quantity`
            Column emission measure
        """
        g = self.contribution_function(density, **kwargs)
        return 1/(4.*np.pi*u.steradian) * g * emission_measure[:, np.newaxis, np.newaxis]

    def spectrum(self, *args, **kwargs):
        """
        Construct the spectrum using a given filter over a specified wavelength range.

        See Also
        --------
        fiasco.IonCollection.spectrum : Compute spectrum for multiple ions
        intensity : Compute LOS intensity for all transitions
        """
        return IonCollection(self).spectrum(*args, **kwargs)

    @needs_dataset('ip')
    @u.quantity_input
    def direct_ionization_rate(self) -> u.cm**3 / u.s:
        """
        Calculate direct ionization rate

        Needs an equation reference or explanation
        """
        xgl, wgl = np.polynomial.laguerre.laggauss(12)
        kBT = const.k_B*self.temperature
        energy = np.outer(xgl, kBT) + self.ip
        cross_section = self.direct_ionization_cross_section(energy)
        if cross_section is None:
            return None
        term1 = np.sqrt(8./np.pi/const.m_e)*np.sqrt(kBT)*np.exp(-self.ip/kBT)
        term2 = ((wgl*xgl)[:, np.newaxis]*cross_section).sum(axis=0)
        term3 = (wgl[:, np.newaxis]*cross_section).sum(axis=0)*self.ip/kBT
        return term1*(term2 + term3)

    @u.quantity_input
    def direct_ionization_cross_section(self, energy: u.erg) -> u.cm**2:
        """
        Calculate direct ionization cross-section.

        The cross-sections are calculated one of two ways:

        - Using the method of [1]_ for hydrogenic and He-like ions
        - Using the scaled cross-sections of [2]_ for all other ions

        References
        ----------
        .. [1] Fontes, C. J., et al., 1999, Phys. Rev. A., `59 1329 <https://journals.aps.org/pra/abstract/10.1103/PhysRevA.59.1329>`_
        .. [2] Dere, K. P., 2007, A&A, `466, 771 <http://adsabs.harvard.edu/abs/2007A%26A...466..771D>`_
        """
        if self.hydrogenic or self.helium_like:
            return self._fontes_cross_section(energy)
        else:
            return self._dere_cross_section(energy)

    @needs_dataset('diparams')
    @u.quantity_input
    def _dere_cross_section(self, energy: u.erg) -> u.cm**2:
        """
        Calculate direct ionization cross-sections according to [1]_.

        References
        ----------
        .. [1] Dere, K. P., 2007, A&A, `466, 771 <http://adsabs.harvard.edu/abs/2007A%26A...466..771D>`_
        """
        # Cross-sections from diparams file
        cross_section_total = np.zeros(energy.shape)
        for ip, bt_c, bt_e, bt_cross_section in zip(self._diparams['ip'],
                                                    self._diparams['bt_c'],
                                                    self._diparams['bt_e'],
                                                    self._diparams['bt_cross_section']):
            U = energy/ip
            scaled_energy = 1. - np.log(bt_c)/np.log(U - 1. + bt_c)
            f_interp = interp1d(bt_e, bt_cross_section, kind='cubic', fill_value='extrapolate')
            scaled_cross_section = f_interp(scaled_energy) * bt_cross_section.unit
            # Only nonzero at energies above the ionization potential
            scaled_cross_section *= (U > 1.)
            cross_section = scaled_cross_section * (np.log(U) + 1.) / U / (ip**2)
            if not hasattr(cross_section_total, 'unit'):
                cross_section_total = cross_section_total*cross_section.unit
            cross_section_total += cross_section

        return cross_section_total

    @needs_dataset('ip')
    @u.quantity_input
    def _fontes_cross_section(self, energy: u.erg) -> u.cm**2:
        """
        Calculate direct ionization cross-section according to [1]_.

        References
        ----------
        .. [1] Fontes, C. J., et al., 1999, Phys. Rev. A., `59 1329 <https://journals.aps.org/pra/abstract/10.1103/PhysRevA.59.1329>`_
        """
        U = energy/self.ip
        A = 1.13
        B = 1 if self.hydrogenic else 2
        F = 1 if self.atomic_number < 20 else (140 + (self.atomic_number/20)**3.2)/141
        if self.atomic_number >= 16:
            c, d, C, D = -0.28394, 1.95270, 0.20594, 3.70590
            if self.atomic_number > 20:
                C += ((self.atomic_number - 20)/50.5)**1.11
        else:
            c, d, C, D = -0.80414, 2.32431, 0.14424, 3.82652

        Qrp = 1./U * (A * np.log(U) + D * (1. - 1./U)**2 + C * U * (1. - 1./U)**4
                      + (c / U + d / U**2) * (1. - 1. / U))

        # NOTE: conversion to Rydbergs equivalent to scaling to the ionization energy
        # of hydrogen such that it is effectively unitless
        return B * (np.pi * const.a0**2) * F * Qrp / (self.ip.to(u.Ry).value**2)

    @needs_dataset('easplups')
    @u.quantity_input
    def excitation_autoionization_rate(self) -> u.cm**3 / u.s:
        """
        Calculate ionization rate due to excitation autoionization
        """
        c = (const.h**2)/((2. * np.pi * const.m_e)**(1.5) * np.sqrt(const.k_B))
        kBTE = np.outer(const.k_B*self.temperature, 1.0/self._easplups['delta_energy'])
        # NOTE: Transpose here to make final dimensions compatible with multiplication with
        # temperature when computing rate
        kBTE = kBTE.T
        xs = [np.linspace(0, 1, ups.shape[0]) for ups in self._easplups['bt_upsilon']]
        upsilon = burgess_tully_descale(xs,
                                        self._easplups['bt_upsilon'].value,
                                        kBTE,
                                        self._easplups['bt_c'].value,
                                        self._easplups['bt_type'])
        # NOTE: The 1/omega multiplicity factor is already included in the scaled upsilon
        # values provided by CHIANTI
        rate = c * upsilon * np.exp(-1 / kBTE) / np.sqrt(self.temperature)

        return rate.sum(axis=0)

    @u.quantity_input
    def ionization_rate(self) -> u.cm**3 / u.s:
        """
        Total ionization rate.

        Includes contributions from both direct ionization and excitation-autoionization

        See Also
        --------
        direct_ionization_rate
        excitation_autoionization_rate
        """
        di_rate = self.direct_ionization_rate()
        di_rate = np.zeros(self.temperature.shape)*u.cm**3/u.s if di_rate is None else di_rate
        ea_rate = self.excitation_autoionization_rate()
        ea_rate = np.zeros(self.temperature.shape)*u.cm**3/u.s if ea_rate is None else ea_rate
        return di_rate + ea_rate

    @needs_dataset('rrparams')
    @u.quantity_input
    def radiative_recombination_rate(self) -> u.cm**3 / u.s:
        """
        Radiative recombination rate

        The recombination rate due to interaction with the ambient radiation field
        is calculated using a set of fit parameters using one of two methods:

        - Method of [1]_, (show expression)
        - Method of [2]_, (show expression)

        References
        ----------
        .. [1] Badnell, N. R., 2006, APJS, `167 334 <http://adsabs.harvard.edu/abs/2006ApJS..167..334B>`_
        .. [2] Shull, J. M., Van Steenberg, M., 1982, `48 95 <http://adsabs.harvard.edu/abs/1982ApJS...48...95S>`_
        """
        if self._rrparams['fit_type'][0] == 1 or self._rrparams['fit_type'][0] == 2:
            A = self._rrparams['A_fit']
            B = self._rrparams['B_fit']
            if self._rrparams['fit_type'] == 2:
                B = B + self._rrparams['C_fit']*np.exp(-self._rrparams['T2_fit']/self.temperature)
            T0 = self._rrparams['T0_fit']
            T1 = self._rrparams['T1_fit']

            return A/(np.sqrt(self.temperature/T0) * (1 + np.sqrt(self.temperature/T0))**(1. - B)
                      * (1. + np.sqrt(self.temperature/T1))**(1. + B))
        elif self._rrparams['fit_type'][0] == 3:
            return self._rrparams['A_fit'] * (
                    (self.temperature/(1e4*u.K))**(-self._rrparams['eta_fit']))
        else:
            raise ValueError(f"Unrecognized fit type {self._rrparams['fit_type']}")

    @needs_dataset('drparams')
    @u.quantity_input
    def dielectronic_recombination_rate(self) -> u.cm**3 / u.s:
        """
        Dielectronic recombination rate

        Calculated according to one of two methods,

        - Method of [1]_, (show expression)
        - Method of [2]_, (show expression)

        References
        ----------
        """
        if self._drparams['fit_type'][0] == 1:
            E_over_T = np.outer(self._drparams['E_fit'], 1./self.temperature)
            return self.temperature**(-1.5)*(
                    self._drparams['c_fit'][:, np.newaxis]*np.exp(-E_over_T)).sum(axis=0)
        elif self._drparams['fit_type'][0] == 2:
            A = self._drparams['A_fit']
            B = self._drparams['B_fit']
            T0 = self._drparams['T0_fit']
            T1 = self._drparams['T1_fit']
            return A * self.temperature**(-1.5) * np.exp(-T0/self.temperature) * (
                    1. + B * np.exp(-T1/self.temperature))
        else:
            raise ValueError(f"Unrecognized fit type {self._drparams['fit_type']}")

    @u.quantity_input
    def recombination_rate(self) -> u.cm**3 / u.s:
        """
        Total recombination rate.

        Includes contributions from dielectronic recombination and radiative recombination.

        See Also
        --------
        radiative_recombination_rate
        dielectronic_recombination_rate
        """
        rr_rate = self.radiative_recombination_rate()
        rr_rate = np.zeros(self.temperature.shape)*u.cm**3/u.s if rr_rate is None else rr_rate
        dr_rate = self.dielectronic_recombination_rate()
        dr_rate = np.zeros(self.temperature.shape)*u.cm**3/u.s if dr_rate is None else dr_rate
        return rr_rate + dr_rate

    @u.quantity_input
    def free_free(self, wavelength: u.angstrom) -> u.erg * u.cm**3 / u.s / u.angstrom:
        """
        Free-free continuum emission or bremsstrahlung

        .. note:: Does not include ionization equilibrium or abundance

        Parameters
        ----------
        wavelength : `~astropy.units.Quantity`
        """
        prefactor = (const.c / 3. / const.m_e
                     * (const.alpha * const.h / np.pi)**3
                     * np.sqrt(2. * np.pi / 3. / const.m_e / const.k_B))
        tmp = np.outer(self.temperature, wavelength)
        exp_factor = np.exp(-const.h * const.c / const.k_B / tmp) / (wavelength**2)
        gf = self._gaunt_factor_free_free(wavelength)

        return (prefactor * self.atomic_number**2 * exp_factor * gf
                / np.sqrt(self.temperature)[:, np.newaxis])

    @needs_dataset('fblvl', 'ip')
    @u.quantity_input
    def free_bound(self,
                   wavelength: u.angstrom,
                   use_verner=True) -> u.erg * u.cm**3 / u.s / u.angstrom:
        """
        Free-bound continuum emission of the recombined ion.

        .. note:: Does not include ionization equilibrium or abundance

        Parameters
        ----------
        wavelength : `~astropy.units.Quantity`
        use_verner : `bool`, optional
            If True, evaluate ground-state cross-sections using method of Verner and Yakovlev

        See Also
        --------
        _verner_cross_section
        """
        prefactor = (2/np.sqrt(2*np.pi)/(4*np.pi)/(
            const.h*(const.c**3) * (const.m_e * const.k_B)**(3/2)))
        recombining = Ion(f'{self.element_name} {self.ionization_stage + 1}', self.temperature,
                          **self._dset_names)
        omega_0 = 1. if recombining._fblvl is None else recombining._fblvl['multiplicity'][0]
        E_photon = const.h * const.c / wavelength
        energy_temperature_factor = np.outer(self.temperature**(-3/2), E_photon**5)
        # Sum over levels of recombined ion
        sum_factor = u.Quantity(np.zeros(self.temperature.shape+wavelength.shape), 'cm^2')
        for omega, E, n, L, level in zip(self._fblvl['multiplicity'],
                                         self._fblvl['E_obs']*const.h*const.c,
                                         self._fblvl['n'],
                                         self._fblvl['L'],
                                         self._fblvl['level']):
            # Energy required to ionize ion from level i
            E_ionize = self.ip - E
            # Check if ionization potential and photon energy sufficient
            if (E_ionize < 0*u.erg) or (E_photon.max() < E):
                continue
            # Only use Verner cross-section for ground state
            if level == 1 and use_verner:
                cross_section = self._verner_cross_section(E_photon)
            else:
                cross_section = self._karzas_cross_section(E_photon, E_ionize, n, L)
            E_scaled = np.outer(1/(const.k_B*self.temperature), E_photon - E_ionize)
            # Scaled energy can blow up at low temperatures; not an issue when cross-section is 0
            E_scaled[:, np.where(cross_section == 0*cross_section.unit)] = 0.0
            sum_factor += omega / omega_0 * np.exp(-E_scaled) * cross_section

        return (prefactor * energy_temperature_factor * sum_factor)

    def free_free_loss(self):
        """
        Wavelength-integrated radiative losses due to free-free emission
        """
        raise NotImplementedError

    def free_bound_loss(self):
        """
        Wavelength-integrated radiative losses due to free-bound emission
        """
        raise NotImplementedError

    @u.quantity_input
    def _gaunt_factor_free_free(self, wavelength: u.angstrom) -> u.dimensionless_unscaled:
        """
        Compute free-free gaunt factor using the approaches of [1]_
        and [2]_ where appropriate.

        References
        ----------
        .. [1] Itoh, N. et al., 2000, ApJS, `128, 125
            <http://adsabs.harvard.edu/abs/2000ApJS..128..125I>`_
        .. [2] Sutherland, R. S., 1998, MNRAS, `300, 321
            <http://adsabs.harvard.edu/abs/1998MNRAS.300..321S>`_
        """
        gf_itoh = self._gaunt_factor_free_free_itoh(wavelength)
        gf_sutherland = self._gaunt_factor_free_free_sutherland(wavelength)
        gf = np.where(np.isnan(gf_itoh), gf_sutherland, gf_itoh)

        return gf

    @u.quantity_input
    def _gaunt_factor_free_free_itoh(self, wavelength: u.angstrom) -> u.dimensionless_unscaled:
        """
        Calculates the free-free gaunt factor of [1]_.

        Need some equations here...

        Notes
        -----
        The relativistic values are valid for :math:`6<\log_{10}(T)< 8.5` and
        :math:`-4<\log_{10}(u)<1`

        References
        ----------
        .. [1] Itoh, N. et al., 2000, ApJS, `128, 125
            <http://adsabs.harvard.edu/abs/2000ApJS..128..125I>`_
        """
        log10_temperature = np.log10(self.temperature.to(u.K).value)
        # calculate scaled energy and temperature
        tmp = np.outer(self.temperature, wavelength)
        lower_u = const.h * const.c / const.k_B / tmp
        upper_u = 1. / 2.5 * (np.log10(lower_u) + 1.5)
        t = 1. / 1.25 * (log10_temperature - 7.25)
        itoh_coefficients = self._itoh['a']
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

    @u.quantity_input
    def _gaunt_factor_free_free_sutherland(self,
                                           wavelength: u.angstrom) -> u.dimensionless_unscaled:
        """
        Calculates the free-free gaunt factor calculations of [1]_.

        Need some equations here.

        References
        ----------
        .. [1] Sutherland, R. S., 1998, MNRAS, `300, 321
            <http://adsabs.harvard.edu/abs/1998MNRAS.300..321S>`_
        """
        Ry = const.h * const.c * const.Ryd
        tmp = np.outer(self.temperature, wavelength)
        lower_u = const.h * const.c / const.k_B / tmp
        gamma_squared = ((self.atomic_number**2) * Ry / const.k_B / self.temperature[:, np.newaxis]
                         * np.ones(lower_u.shape))
        # convert to index coordinates
        i_lower_u = (np.log10(lower_u) + 4.) * 10.
        i_gamma_squared = (np.log10(gamma_squared) + 4.) * 5.
        # interpolate data to scaled quantities
        # FIXME: interpolate without reshaping?
        gf_data = self._gffgu['gaunt_factor'].reshape(
            np.unique(self._gffgu['u']).shape[0],
            np.unique(self._gffgu['gamma_squared']).shape[0],
        )
        gf = map_coordinates(
            gf_data, [i_gamma_squared.flatten(), i_lower_u.flatten()]).reshape(lower_u.shape)

        return u.Quantity(np.where(gf < 0., 0., gf))

    @needs_dataset('verner')
    @u.quantity_input
    def _verner_cross_section(self, energy: u.erg) -> u.cm**2:
        """
        Ground state photoionization cross-section using the method of [1]_.

        Parameters
        ----------
        energy : `~astropy.units.Quantity`
            Photon energy

        References
        ----------
        .. [1] Verner & Yakovlev, 1995, A&AS, `109, 125
            <http://adsabs.harvard.edu/abs/1995A%26AS..109..125V>`_
        """
        # decompose simplifies units and makes sure y is unitless
        y = (energy / self._verner['E_0_fit']).decompose()
        Q = 5.5 + self._verner['l'] - 0.5*self._verner['P_fit']
        F = ((y - 1)**2 + self._verner['y_w_fit']**2) * (y**(-Q))*(
            1. + np.sqrt(y / self._verner['y_a_fit']))**(-self._verner['P_fit'])
        return np.where(energy < self._verner['E_thresh'], 0.,
                        F.decompose().value) * self._verner['sigma_0']

    @u.quantity_input
    def _karzas_cross_section(self,
                              photon_energy: u.erg,
                              ionization_energy: u.erg,
                              n,
                              l) -> u.cm**2:
        """
        Photoionization cross-section using the method of [1]_.

        Parameters
        ----------
        photon_energy : `~astropy.units.Quantity`
            Energy of emitted photon
        ionization_energy : `~astropy.units.Quantity`
            Ionization potential of recombined ion for level `n`
        n : `int`
            Principal quantum number
        l : `int`
            Orbital angular momentum number

        References
        ----------
        .. [1] Karzas and Latter, 1961, ApJSS, `6, 167
            <http://adsabs.harvard.edu/abs/1961ApJS....6..167K>`_
        """
        prefactor = (2**4)*const.h*(const.e.gauss**2)/(3*np.sqrt(3)*const.m_e*const.c)
        index_nl = np.where(np.logical_and(self._klgfb['n'] == n, self._klgfb['l'] == l))[0]
        # If there is no Gaunt factor for n, l, set it to 1
        if index_nl.shape == (0,):
            gaunt_factor = 1
        else:
            E_scaled = np.log(photon_energy/ionization_energy)
            gf_interp = splrep(self._klgfb['log_pe'][index_nl, :].squeeze(),
                               self._klgfb['log_gaunt_factor'][index_nl, :].squeeze())
            gaunt_factor = np.exp(splev(E_scaled, gf_interp))
        cross_section = prefactor * ionization_energy**2 * photon_energy**(-3) * gaunt_factor / n
        cross_section[np.where(photon_energy < ionization_energy)] = 0.*cross_section.unit
        return cross_section
