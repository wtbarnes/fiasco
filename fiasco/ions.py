"""
Ion object. Holds all methods and properties of a CHIANTI ion.
"""
import astropy.constants as const
import astropy.units as u
import numpy as np

from functools import cached_property
from scipy.interpolate import interp1d, splev, splrep
from scipy.ndimage import map_coordinates

from fiasco import proton_electron_ratio
from fiasco.base import ContinuumBase, IonBase
from fiasco.collections import IonCollection
from fiasco.levels import Level, Transitions
from fiasco.util import (
    burgess_tully_descale,
    needs_dataset,
    vectorize_where,
    vectorize_where_sum,
)
from fiasco.util.exceptions import MissingDatasetException

__all__ = ['Ion']


class Ion(IonBase, ContinuumBase):
    """
    Class for representing a CHIANTI ion.

    The ion object is the fundamental unit of `fiasco`. This object contains
    all of the properties and methods needed to access important information about each ion
    from the CHIANTI database as well as compute common derived quantities.

    Parameters
    ----------
    ion_name : `str` or `tuple`
        Name of the ion. This can be either a string denoting the name or a tuple containing the
        atomic number and ionization stage. See `~fiasco.util.parse_ion` for a list of all possible
        input formats.
    temperature : `~astropy.units.Quantity`
        Temperature array over which to evaluate temperature dependent quantities.
    ioneq_filename : `str`, optional
        Ionization equilibrium dataset
    abundance_filename : `str`, optional
        Abundance dataset
    ip_filename : `str`, optional
        Ionization potential dataset
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
                                                            'sun_coronal_1992_feldman_ext')
        self._dset_names['ip_filename'] = kwargs.get('ip_filename', 'chianti')

    def _new_instance(self, temperature=None, **kwargs):
        """
        Convenience method for creating an ion of the same type with
        possibly different arguments. If different arguments are not
        specified, this will just create a copy of itself.
        """
        if temperature is None:
            temperature = self.temperature.copy()
        new_kwargs = self._instance_kwargs
        new_kwargs.update(kwargs)
        return type(self)(self.ion_name, temperature, **new_kwargs)

    def __repr__(self):
        try:
            n_levels = self._elvlc['level'].shape[0]
        except KeyError:
            n_levels = 0
        try:
            n_transitions = self._wgfa['lower_level'].shape[0]
        except KeyError:
            n_transitions = 0
        return f"""CHIANTI Database Ion
---------------------
Name: {self.ion_name}
Element: {self.element_name} ({self.atomic_number})
Charge: +{self.charge_state}
Number of Levels: {n_levels}
Number of Transitions: {n_transitions}

Temperature range: [{self.temperature[0].to(u.MK):.3f}, {self.temperature[-1].to(u.MK)}:.3f]

HDF5 Database: {self.hdf5_dbase_root}
Using Datasets:
  ioneq: {self._dset_names['ioneq_filename']}
  abundance: {self._dset_names['abundance_filename']}
  ip: {self._dset_names['ip_filename']}"""

    @cached_property
    def _all_levels(self):
        try:
            _ = self._elvlc
        except KeyError:
            return None
        else:
            n_levels = self._elvlc['level'].shape[0]
            return [Level(i, self._elvlc) for i in range(n_levels)]

    def __getitem__(self, key):
        if self._all_levels is None:
            raise IndexError(f'No energy levels available for {self.ion_name}')
        else:
            # NOTE: casting to array first so that "fancy indexing" can
            # be used to index the energy levels.
            indexed_levels = np.array(self._all_levels)[key]
            if isinstance(indexed_levels, Level):
                return indexed_levels
            else:
                return indexed_levels.tolist()

    def __add__(self, value):
        return IonCollection(self, value)

    def __radd__(self, value):
        return IonCollection(value, self)

    @property
    def _instance_kwargs(self):
        # Keyword arguments used to istantiate this Ion. These are useful when
        # constructing a new Ion instance that pulls from exactly the same
        # data sources.
        kwargs = {
            'hdf5_dbase_root': self.hdf5_dbase_root,
            **self._dset_names,
        }
        return kwargs

    def next_ion(self):
        """
        Return an `~fiasco.Ion` instance with the next highest ionization stage.

        For example, if the current instance is Fe XII (+11), this method returns
        an instance of Fe XIII (+12). All other input arguments remain the same.
        """
        return type(self)((self.atomic_number, self.ionization_stage+1),
                           self.temperature,
                           **self._instance_kwargs)

    def previous_ion(self):
        """
        Return an `~fiasco.Ion` instance with the next lowest ionization stage.

        For example, if the current instance is Fe XII (+11), this method returns
        an instance of Fe XI (+10). All other input arguments remain the same.
        """
        return type(self)((self.atomic_number, self.ionization_stage-1),
                           self.temperature,
                           **self._instance_kwargs)

    @property
    @needs_dataset('elvlc', 'wgfa')
    def transitions(self):
        return Transitions(self._elvlc, self._wgfa)

    @cached_property
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
        f = interp1d(self._ioneq[self._dset_names['ioneq_filename']]['temperature'].to('MK').value,
                     self._ioneq[self._dset_names['ioneq_filename']]['ionization_fraction'],
                     kind='linear',
                     bounds_error=False,
                     fill_value=np.nan)
        ioneq = f(self.temperature.to('MK').value)
        isfinite = np.isfinite(ioneq)
        ioneq[isfinite] = np.where(ioneq[isfinite] < 0., 0., ioneq[isfinite])
        return u.Quantity(ioneq)

    @property
    def abundance(self):
        """
        Elemental abundance relative to H.
        """
        return self._abundance[self._dset_names['abundance_filename']]

    @property
    @needs_dataset('ip')
    @u.quantity_input
    def ip(self) -> u.erg:
        """
        Ionization potential.
        """
        return self._ip[self._dset_names['ip_filename']] * const.h * const.c

    @property
    def hydrogenic(self):
        r"""
        Is the ion hydrogen-like or not.

        Notes
        -----
        This is `True` if :math:`Z - z = 1` and :math:`Z\ge6`.
        """
        return (self.atomic_number - self.charge_state == 1) and (self.atomic_number >= 6)

    @property
    def helium_like(self):
        r"""
        Is the ion helium like or not.

        Notes
        -----
        This is `True` if :math:`Z - z = 2` and :math:`Z\ge10`.
        """
        return (self.atomic_number - self.charge_state == 2) and (self.atomic_number >= 10)

    @property
    @u.quantity_input
    def formation_temperature(self) -> u.K:
        """
        Temperature at which `~fiasco.Ion.ioneq` is maximum. This is a useful proxy for
        the temperature at which lines for this ion are formed.
        """
        return self.temperature[np.argmax(self.ioneq)]

    @cached_property
    @needs_dataset('scups')
    @u.quantity_input
    def effective_collision_strength(self) -> u.dimensionless_unscaled:
        r"""
        Maxwellian-averaged collision strength, typically denoted by :math:`\Upsilon`, as a function of temperature.

        According to Eq. 4.11 of :cite:t:`phillips_ultraviolet_2008`,
        :math:`\Upsilon` is given by,

        .. math::

            \Upsilon = \int_0^\infty\mathrm{d}\left(\frac{E}{k_BT}\right)\,\Omega_{ji}\exp{\left(-\frac{E}{k_BT}\right)}

        where :math:`\Omega_{ji}` is the collision strength.
        These Maxwellian-averaged collision strengths are stored in
        dimensionless form in CHIANTI and are rescaled to the appropriate
        temperature.

        See Also
        --------
        fiasco.util.burgess_tully_descale : Descale and interpolate :math:`\Upsilon`.
        """
        kBTE = np.outer(const.k_B * self.temperature, 1.0 / self._scups['delta_energy'])
        upsilon = burgess_tully_descale(self._scups['bt_t'],
                                        self._scups['bt_upsilon'],
                                        kBTE.T,
                                        self._scups['bt_c'],
                                        self._scups['bt_type'])
        upsilon = u.Quantity(np.where(upsilon > 0., upsilon, 0.))
        return upsilon.T

    @cached_property
    @needs_dataset('elvlc', 'scups')
    @u.quantity_input
    def electron_collision_deexcitation_rate(self) -> u.cm**3 / u.s:
        r"""
        Collisional de-excitation rate coefficient for electrons.

        According to Eq. 4.12 of :cite:t:`phillips_ultraviolet_2008`, the rate coefficient for collisional de-excitation
        is given by,

        .. math::

           C^d_{ji} = I_Ha_0^2\sqrt{\frac{8\pi}{mk_B}}\frac{\Upsilon}{\omega_jT^{1/2}},

        where :math:`j,i` are the upper and lower level indices, respectively, :math:`I_H` is the
        ionization potential for H, :math:`a_0` is the Bohr radius, :math:`\Upsilon` is the
        effective collision strength, and :math:`\omega_j` is the statistical weight of the
        level :math:`j`.

        See Also
        --------
        electron_collision_excitation_rate : Excitation rate due to collisions
        effective_collision_strength : Maxwellian-averaged collision strength, :math:`\Upsilon`
        """
        c = (const.h**2) / ((2. * np.pi * const.m_e)**(1.5) * np.sqrt(const.k_B))
        upsilon = self.effective_collision_strength
        omega_upper = 2. * self._elvlc['J'][self._scups['upper_level'] - 1] + 1.
        return c * upsilon / np.sqrt(self.temperature[:, np.newaxis]) / omega_upper

    @cached_property
    @needs_dataset('elvlc', 'scups')
    @u.quantity_input
    def electron_collision_excitation_rate(self) -> u.cm**3 / u.s:
        r"""
        Collisional excitation rate coefficient for electrons.

        The rate coefficient for collisional excitation is given by,

        .. math::

            C^e_{ij} = \frac{\omega_j}{\omega_i}C^d_{ji}\exp{\left(-\frac{k_BT_e}{\Delta E_{ij}}\right)}

        where :math:`j,i` are the upper and lower level indices, respectively, :math:`\omega_j,\omega_i`
        are the statistical weights of the upper and lower levels, respectively, and :math:`\Delta E_{ij}`
        is the energy of the transition :cite:p:`phillips_ultraviolet_2008`.

        Parameters
        ----------
        deexcitation_rate : `~astropy.units.Quantity`, optional
            Optionally specify deexcitation rate to speedup calculation

        See Also
        --------
        electron_collision_deexcitation_rate : De-excitation rate due to collisions
        """
        omega_upper = 2. * self._elvlc['J'][self._scups['upper_level'] - 1] + 1.
        omega_lower = 2. * self._elvlc['J'][self._scups['lower_level'] - 1] + 1.
        kBTE = np.outer(1./const.k_B/self.temperature, self._scups['delta_energy'])
        return omega_upper / omega_lower * self.electron_collision_deexcitation_rate * np.exp(-kBTE)

    @cached_property
    @needs_dataset('psplups')
    @u.quantity_input
    def proton_collision_excitation_rate(self) -> u.cm**3 / u.s:
        """
        Collisional excitation rate coefficient for protons.

        These excitation rates are stored in CHIANTI and then rescaled
        to the appropriate temperatures using the method of
        :cite:t:`burgess_analysis_1992`.

        See Also
        --------
        electron_collision_excitation_rate
        """
        # Create scaled temperature--these are not stored in the file
        bt_t = [np.linspace(0, 1, ups.shape[0]) for ups in self._psplups['bt_rate']]
        # Get excitation rates directly from scaled data
        kBTE = np.outer(const.k_B * self.temperature, 1.0 / self._psplups['delta_energy'])
        ex_rate = burgess_tully_descale(bt_t,
                                        self._psplups['bt_rate'],
                                        kBTE.T,
                                        self._psplups['bt_c'],
                                        self._psplups['bt_type'])
        with np.errstate(invalid='ignore'):
            return u.Quantity(np.where(ex_rate > 0., ex_rate, 0.), u.cm**3/u.s).T

    @cached_property
    @needs_dataset('elvlc', 'psplups')
    @u.quantity_input
    def proton_collision_deexcitation_rate(self) -> u.cm**3 / u.s:
        r"""
        Collisional de-excitation rate coefficient for protons.

        As in the electron case, the proton collision de-excitation
        rate is given by,

        .. math::

            C^{d,p}_{ji} = \frac{\omega_i}{\omega_j}\exp{\left(\frac{E}{k_BT}\right)}C^{e,p}_{ij}

        where :math:`C^{e,p}_{ji}` is the excitation rate due to collisions
        with protons.

        Note that :math:`T` is technically the proton temperature. In the case of a thermal plasma, the electron and proton temperatures are equal, :math:`T_e=T_p`.
        See Section 4.9.4 of :cite:t:`phillips_ultraviolet_2008` for additional information on proton collision rates.

        See Also
        --------
        proton_collision_excitation_rate : Excitation rate due to collisions with protons
        """
        kBTE = np.outer(const.k_B * self.temperature, 1.0 / self._psplups['delta_energy'])
        omega_upper = 2. * self._elvlc['J'][self._psplups['upper_level'] - 1] + 1.
        omega_lower = 2. * self._elvlc['J'][self._psplups['lower_level'] - 1] + 1.
        dex_rate = (omega_lower / omega_upper) * self.proton_collision_excitation_rate * np.exp(1. / kBTE)

        return dex_rate

    @needs_dataset('elvlc', 'wgfa', 'scups')
    @u.quantity_input
    def level_populations(self,
                          density: u.cm**(-3),
                          include_protons=True) -> u.dimensionless_unscaled:
        """
        Energy level populations as a function of temperature and density.

        Parameters
        ----------
        density : `~astropy.units.Quantity`
        include_protons : `bool`, optional
            If True (default), include proton excitation and de-excitation rates.

        Returns
        -------
        `~astropy.units.Quantity`
            A ``(l, m, n)`` shaped quantity, where ``l`` is the number of
            temperatures, ``m`` is the number of densities, and ``n``
            is the number of energy levels.
        """
        # NOTE: Cannot include protons if psplups data not available
        try:
            _ = self._psplups
        except KeyError:
            self.log.warning(
                f'No proton data available for {self.ion_name}. '
                'Not including proton excitation and de-excitation in level populations calculation.')
            include_protons = False

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
        dex_rate_e = self.electron_collision_deexcitation_rate
        ex_rate_e = self.electron_collision_excitation_rate
        ex_diagonal_e = vectorize_where_sum(
            lower_level, level, ex_rate_e.value.T, 0).T * ex_rate_e.unit
        dex_diagonal_e = vectorize_where_sum(
            upper_level, level, dex_rate_e.value.T, 0).T * dex_rate_e.unit
        # Collisional--protons
        if include_protons:
            lower_level_p = self._psplups['lower_level']
            upper_level_p = self._psplups['upper_level']
            pe_ratio = proton_electron_ratio(self.temperature, **self._instance_kwargs)
            proton_density = np.outer(pe_ratio, density)[:, :, np.newaxis]
            ex_rate_p = self.proton_collision_excitation_rate
            dex_rate_p = self.proton_collision_deexcitation_rate
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

    @needs_dataset('abundance', 'elvlc')
    @u.quantity_input
    def contribution_function(self, density: u.cm**(-3), **kwargs) -> u.cm**3 * u.erg / u.s:
        r"""
        Contribution function :math:`G(n_e,T)` for all transitions.

        The contribution function for ion :math:`k` of element :math:`X` for a
        particular transition :math:`ij` is given by,

        .. math::

           G_{ij} = \frac{n_H}{n_e}\mathrm{Ab}(X)f_{X,k}N_jA_{ij}\Delta E_{ij}\frac{1}{n_e},

        Note that the contribution function is often defined in differing ways by different authors.
        The contribution function is defined as above in :cite:t:`young_chianti_2016`.

        The corresponding wavelengths can be retrieved with,

        .. code-block:: python

           ion.transitions.wavelength[~ion.transitions.is_twophoton]

        Parameters
        ----------
        density : `~astropy.units.Quantity`
            Electron number density
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

    @u.quantity_input
    def emissivity(self, density: u.cm**(-3), **kwargs) -> u.erg * u.cm**(-3) / u.s:
        r"""
        Emissivity as a function of temperature and density for all transitions.

        The emissivity is given by the expression,

        .. math::

           \epsilon(n_e,T) = G(n_e,T)n_e^2

        where :math:`G` is the contribution function, :math:`n_e` is the electron
        density, and :math:`T` is the temperature.
        Note that, like the contribution function, emissivity is often defined in
        in differing ways by different authors.
        Here, we use the definition of the emissivity as given by Eq. 3 of
        :cite:t:`young_chianti_2016`.

        Parameters
        ----------
        density : `~astropy.units.Quantity`
            Electron number density

        See Also
        --------
        contribution_function : Calculate contribution function, :math:`G(n,T)`
        """
        density = np.atleast_1d(density)
        g = self.contribution_function(density, **kwargs)
        return g * (density**2)[np.newaxis, :, np.newaxis]

    @u.quantity_input
    def intensity(self,
                  density: u.cm**(-3),
                  emission_measure: u.cm**(-5), **kwargs) -> u.erg / u.cm**2 / u.s / u.steradian:
        r"""
        Line-of-sight intensity computed assuming a particular column emission measure.

        The intensity along the line-of-sight can be written as,

        .. math::

           I = \frac{1}{4\pi}\int\mathrm{d}T,G(n,T)n^2\frac{dh}{dT}

        which, in the isothermal approximation, can be simplified to,

        .. math::

           I(T_0) \approx \frac{1}{4\pi}G(n,T_0)\mathrm{EM}(T_0)

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
        emission_measure = np.atleast_1d(emission_measure)
        g = self.contribution_function(density, **kwargs)
        return 1/(4.*np.pi*u.steradian) * g * emission_measure[:, np.newaxis, np.newaxis]

    def spectrum(self, *args, **kwargs):
        """
        Construct the spectrum using a given filter over a specified wavelength range.

        All arguments are passed directly to `fiasco.IonCollection.spectrum`.

        See Also
        --------
        fiasco.IonCollection.spectrum : Compute spectrum for multiple ions
        intensity : Compute LOS intensity for all transitions
        """
        return IonCollection(self).spectrum(*args, **kwargs)

    @cached_property
    @needs_dataset('ip')
    @u.quantity_input
    def direct_ionization_rate(self) -> u.cm**3 / u.s:
        r"""
        Total ionization rate due to collisions as a function of temperature.

        The ionization rate due to the collisions with free electrons can
        be written as the integral of the velocity-weighted collisional
        cross-section over the Maxwell-Boltzmann distribution.
        Following Section 3.5.1 of :cite:t:`del_zanna_solar_2018`, this can be
        written as,

        .. math::

            C^I = \sqrt{\frac{8}{\pi m_e}}(k_BT)^{-3/2}\int_I^{\infty}\mathrm{d}E\,E\sigma_I(E)\exp{\left(-\frac{E}{k_BT}\right)}

        where :math:`E` is the energy of the incident electron,
        :math:`I` is the ionization energy of the initially bound electron,
        and :math:`\sigma_I` is the ionization cross-section.

        Making the substitution :math:`x=(E-I)/k_BT`, the above integral can be
        rewritten as,

        .. math::

            \begin{aligned}
                C^I = \sqrt{\frac{8k_BT}{\pi m_e}}\exp{\left(-\frac{I}{k_BT}\right)}&\left(\int_0^{\infty}\mathrm{d}x\,x\sigma_{I}(k_BTx+I)e^{-x} \right. \\
                                                                                    &\left. + \frac{I}{k_BT}\int_0^{\infty}\mathrm{d}x\,\sigma_{I}(k_BTx+I)e^{-x}\right).
            \end{aligned}

        Each of these integrals is of the form such that they can be evaluated using Gauss-Laguerre quadrature.
        Note that there is a typo in the expression for the ionization rate integral in Eq. 32 of :cite:t:`del_zanna_solar_2018`.
        The details of the ionization cross-section calculation can be found in `direct_ionization_cross_section`.

        See Also
        --------
        direct_ionization_cross_section : Calculation of :math:`\sigma_I` as a function of :math:`E`.
        """
        xgl, wgl = np.polynomial.laguerre.laggauss(12)
        kBT = const.k_B * self.temperature
        energy = np.outer(xgl, kBT) + self.ip
        cross_section = self.direct_ionization_cross_section(energy)
        term1 = np.sqrt(8./np.pi/const.m_e)*np.sqrt(kBT)*np.exp(-self.ip/kBT)
        term2 = ((wgl*xgl)[:, np.newaxis]*cross_section).sum(axis=0)
        term3 = (wgl[:, np.newaxis]*cross_section).sum(axis=0)*self.ip/kBT
        return term1*(term2 + term3)

    @u.quantity_input
    def direct_ionization_cross_section(self, energy: u.erg) -> u.cm**2:
        r"""
        Direct ionization cross-section as a function of energy.

        The direction ionization cross-section is calculated one of two ways.
        For H and He like ions, the cross-section is computed according to
        the method of :cite:t:`fontes_fully_1999`,

        .. math::

            \sigma_I = B\frac{\pi a_0^2}{I^2}Q_R

        where :math:`B=1` for H-like ions and :math:`B=2` for He-like ions,
        :math:`I` is the ionization energy (expressed in Ry),
        :math:`a_0` is the Bohr radius,
        and :math:`Q_R` is a reduced cross-section which can be approximated by
        the fitting formula given in Eqs. 2.10, 2.11, and 2.12 of
        :cite:t:`fontes_fully_1999`.

        For all other ions, the cross-section is computed according to the method
        of :cite:t:`dere_ionization_2007` which employs a scaling similar to that
        used by :cite:t:`burgess_analysis_1992`.
        Rearranging Eq. 3 of :cite:t:`dere_ionization_2007`,

        .. math::

            \sigma_I = \frac{\Sigma (\log{u} + 1)}{uI^2}

        where :math:`u=E/I` is the energy of the incident electron scaled by ionization
        potential and :math:`\Sigma` is the scaled cross-section which is defined over,

        .. math::

            U = 1 - \frac{\log{f}}{\log{u - 1 + f}}

        where :math:`f` is a fitting parameter.
        :math:`U,f,\Sigma` are all stored in the CHIANTI database such that :math:`\sigma_I`
        can be computed for a given :math:`E`.
        """
        if self.hydrogenic or self.helium_like:
            return self._fontes_cross_section(energy)
        else:
            return self._dere_cross_section(energy)

    @needs_dataset('diparams')
    @u.quantity_input
    def _dere_cross_section(self, energy: u.erg) -> u.cm**2:
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

    @cached_property
    @needs_dataset('easplups')
    @u.quantity_input
    def excitation_autoionization_rate(self) -> u.cm**3 / u.s:
        r"""
        Ionization rate due to excitation autoionization.

        Following Eq. 4.74 of :cite:t:`phillips_ultraviolet_2008`, the excitation
        autoionization rate is given by,

        .. math::

            \alpha_{EA} = \frac{h^2}{(2\pi m_e)^{3/2}}(k_BT)^{-1/2}\sum_{lj}\Upsilon^{EA}_{lj}\exp{\left(-\frac{\Delta E_{lj}}{k_BT}\right)}

        where :math:`\Upsilon^{EA}` is the thermally-averaged excitation autoionization
        cross-section as stored in CHIANTI and includes the additional :math:`\omega_j`
        multiplicity factor compared to the expression in :cite:t:`phillips_ultraviolet_2008`.
        The sum is taken over inelastic collisions to level :math:`j` from a level :math:`l`
        below the ionization threshold.
        Additionally, note that the constant has been rewritten in terms of :math:`h`
        rather than :math:`I_H` and :math:`a_0`.
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

    @cached_property
    @u.quantity_input
    def ionization_rate(self) -> u.cm**3 / u.s:
        r"""
        Total ionization rate as a function of temperature.

        The total ionization rate, as a function of temperature, for a given ion
        is the sum of the direct ionization and excitation autoionization rates such that,

        .. math::

            \alpha_{I} = \alpha_{DI} + \alpha_{EA}

        See Also
        --------
        direct_ionization_rate
        excitation_autoionization_rate
        """
        try:
            di_rate = self.direct_ionization_rate
        except MissingDatasetException:
            di_rate = u.Quantity(np.zeros(self.temperature.shape), 'cm3 s-1')
        try:
            ea_rate = self.excitation_autoionization_rate
        except MissingDatasetException:
            ea_rate = u.Quantity(np.zeros(self.temperature.shape), 'cm3 s-1')
        return di_rate + ea_rate

    @cached_property
    @needs_dataset('rrparams')
    @u.quantity_input
    def radiative_recombination_rate(self) -> u.cm**3 / u.s:
        r"""
        Radiative recombination rate as a function of temperature.

        The recombination rate due to interaction with the ambient radiation field
        is calculated using a set of fit parameters using one of two methods.
        The methodology used depends on the type of radiative recombination
        rate fitting coefficients available for the particular ion in the CHIANTI atomic
        database.

        The first method is given in Eq. 4 of :cite:t:`verner_atomic_1996` and
        Eq. 1 of :cite:t:`badnell_radiative_2006`,

        .. math::

            \alpha_{RR} = A(\sqrt{T/T_0}(1 + \sqrt{T/T_0})^{1-B}(1 + \sqrt{T/T_1})^{1+B})^{-1}

        where :math:`A,B,T_0,T_1` are fitting coefficients provided for each ion in the CHIANTI
        atomic database.
        In some cases, the fitting coefficient :math:`B` is also modified as,

        .. math::

            B \to B + Ce^{-T_2/T}

        where :math:`C` and :math:`T_2` are additional fitting coefficients
        (see Eq. 2 of :cite:t:`badnell_radiative_2006`).

        The second method is given by Eq. 4 of :cite:t:`shull_ionization_1982`
        and Eq. 1 of :cite:t:`verner_atomic_1996`,

        .. math::

            \alpha_{RR} = A(T/T_0)^{-\eta}

        where :math:`A` and :math:`\eta` are fitting parameters provided in the
        CHIANTI atomic database and :math:`T_0=10^4` K.
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

    @cached_property
    @needs_dataset('drparams')
    @u.quantity_input
    def dielectronic_recombination_rate(self) -> u.cm**3 / u.s:
        r"""
        Dielectronic recombination rate as a function of temperature.

        The dielectronic recombination rate, as a function of :math:`T`, is computed
        using one of two methods.
        The methodology used depends on the type of dielectronic recombination
        rate fitting coefficients available for the particular ion in the CHIANTI atomic
        database.

        The first method is given in Eq. 3 of :cite:t:`zatsarinny_dielectronic_2003`,

        .. math::

            \alpha_{DR} = T^{-3/2}\sum_ic_ie^{-E_i/T}

        where :math:`c_i` and :math:`E_i` are fitting coefficients stored in the CHIANTI
        database.

        The second method is given by Eq. 5 of :cite:t:`shull_ionization_1982`,

        .. math::

            \alpha_{DR} = A T^{-3/2}e^{-T_0/T}(1 + B e^{-T_1/T})

        where :math:`A,B,T_0,T_1` are fitting coefficients stored in the CHIANTI database.
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

    @cached_property
    @u.quantity_input
    def recombination_rate(self) -> u.cm**3 / u.s:
        r"""
        Total recombination rate as a function of temperature.

        The total recombination rate, as a function of temperature, for a given ion
        is the sum of the radiative and dielectronic recombination rates such that,

        .. math::

            \alpha_{R} = \alpha_{RR} + \alpha_{DR}

        See Also
        --------
        radiative_recombination_rate
        dielectronic_recombination_rate
        """
        try:
            rr_rate = self.radiative_recombination_rate
        except MissingDatasetException:
            rr_rate = u.Quantity(np.zeros(self.temperature.shape), 'cm3 s-1')
        try:
            dr_rate = self.dielectronic_recombination_rate
        except MissingDatasetException:
            dr_rate = u.Quantity(np.zeros(self.temperature.shape), 'cm3 s-1')
        return rr_rate + dr_rate

    @u.quantity_input
    def free_free(self, wavelength: u.angstrom) -> u.erg * u.cm**3 / u.s / u.angstrom:
        r"""
        Free-free continuum emission as a function of temperature and wavelength.

        Free-free emission, also known as *bremsstrahlung* (or “braking radiation”),
        is produced when an ion interacts with a free electron, reduces the momentum
        of the free electron, and, by conservation of energy and momentum, produces
        a photon. According to Eq. 4.114 of :cite:t:`phillips_ultraviolet_2008` the
        free-free emission produced by a thermal distribution of electrons as a function
        of temperature and wavelength is given by,

        .. math::

            C_{ff}(\lambda,T_e) = \frac{c}{3m_e}\left(\frac{\alpha h}{\pi}\right)^3\sqrt{\frac{2\pi}{3m_ek_B}}\frac{z^2}{\lambda^2T_e^{1/2}}\exp{\left(-\frac{hc}{\lambda k_BT_e}\right)}\langle g_{ff}\rangle,

        where :math:`\alpha` is the fine-structure constant, :math:`z` is the charge of the ion, and
        :math:`\langle g_{ff}\rangle` is the velocity-averaged free-free Gaunt factor.

        Parameters
        ----------
        wavelength : `~astropy.units.Quantity`

        Notes
        -----
        The result does not include ionization equilibrium or abundance factors.

        See Also
        --------
        gaunt_factor_free_free
        fiasco.IonCollection.free_free: Includes abundance and ionization equilibrium
        """
        prefactor = (const.c / 3. / const.m_e * (const.alpha * const.h / np.pi)**3
                     * np.sqrt(2. * np.pi / 3. / const.m_e / const.k_B))
        tmp = np.outer(self.temperature, wavelength)
        exp_factor = np.exp(-const.h * const.c / const.k_B / tmp) / (wavelength**2)
        gf = self.gaunt_factor_free_free(wavelength)

        return (prefactor * self.charge_state**2 * exp_factor * gf
                / np.sqrt(self.temperature)[:, np.newaxis])

    @u.quantity_input
    def gaunt_factor_free_free(self, wavelength: u.angstrom) -> u.dimensionless_unscaled:
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
        """
        gf_itoh = self._gaunt_factor_free_free_itoh(wavelength)
        gf_sutherland = self._gaunt_factor_free_free_sutherland(wavelength)
        gf = np.where(np.isnan(gf_itoh), gf_sutherland, gf_itoh)

        return gf

    @needs_dataset('itoh')
    @u.quantity_input
    def _gaunt_factor_free_free_itoh(self, wavelength: u.angstrom) -> u.dimensionless_unscaled:
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

    @needs_dataset('gffgu')
    @u.quantity_input
    def _gaunt_factor_free_free_sutherland(self,
                                           wavelength: u.angstrom) -> u.dimensionless_unscaled:
        Ry = const.h * const.c * const.Ryd
        tmp = np.outer(self.temperature, wavelength)
        lower_u = const.h * const.c / const.k_B / tmp
        gamma_squared = ((self.charge_state**2) * Ry / const.k_B / self.temperature[:, np.newaxis]
                         * np.ones(lower_u.shape))
        # NOTE: This escape hatch avoids a divide-by-zero warning as we cannot take log10
        # of 0. This does not matter as the free-free continuum will be 0 for zero charge
        # state anyway.
        if self.charge_state == 0:
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

    @needs_dataset('fblvl', 'ip')
    @u.quantity_input
    def free_bound(self,
                   wavelength: u.angstrom,
                   use_verner=True) -> u.Unit('erg cm3 s-1 Angstrom-1'):
        r"""
        Free-bound continuum emission of the recombined ion.

        .. note:: Does not include ionization equilibrium or abundance.
                  Unlike the equivalent IDL routine, the output here is not
                  expressed per steradian and as such the factor of
                  :math:`1/4\pi` is not included.

        When an electron is captured by an ion of charge :math:`z+1`
        (the recombining ion), it creates a an ion of charge :math:`z`
        (the recombined ion) and produces a continuum of emission
        called the free-bound continuum. The emission of the
        recombined ion is given by,

        .. math::

            C_{fb}(\lambda, T) = \frac{2}{hc^3(k_B m_e)^{3/2}\sqrt{2\pi}}\frac{E^5}{T^{3/2}}\sum_i\frac{\omega_i}{\omega_0}\sigma_i^{\mathrm{bf}}\exp{\left(-\frac{E-I_i}{k_BT}\right)}

        where :math:`E` is the energy of the outgoing photon,
        :math:`\omega_i,\omega_0` are the statastical weights of the
        :math:`i`-th level of the recombined ion and the ground level of the recombining ion, respectively,
        :math:`\sigma_i^{\mathrm{bf}}` is the free-bound cross-section,
        and :math:`I_i` is the energy required to ionize the recombined ion from level :math:`i`.
        A detailed derivation of this formula can be found in
        `CHIANTI Technical Report No. 12 <http://www.chiantidatabase.org/tech_reports/12_freebound/chianti_report_12.pdf>`_.

        For ground state transitions, the photoionization cross-section :math:`\sigma_i^{\mathrm{bf}}` is evaluated
        using Eq. 1 of :cite:t:`verner_analytic_1995` if ``use_verner`` is set to True.
        For all other transitions, and in all cases if ``use_verner`` is set to False, :math:`\sigma_i^{\mathrm{bf}}`
        is evaluated using the method of :cite:t:`karzas_electron_1961`.

        Parameters
        ----------
        wavelength : `~astropy.units.Quantity`
        use_verner : `bool`, optional
            If True, evaluate ground-state cross-sections using method of
            :cite:t:`verner_analytic_1995`.
        """
        wavelength = np.atleast_1d(wavelength)
        prefactor = (2/np.sqrt(2*np.pi)/(const.h*(const.c**3) * (const.m_e * const.k_B)**(3/2)))
        recombining = self.next_ion()
        try:
            # NOTE: This checks whether the fblvl data is available for the
            # recombining ion
            needs_dataset('fblvl')(lambda _: None)(recombining)
            omega_0 = recombining._fblvl['multiplicity'][0]
        except MissingDatasetException:
            omega_0 = 1.0
        E_photon = const.h * const.c / wavelength
        energy_temperature_factor = np.outer(self.temperature**(-3/2), E_photon**5)
        # Fill in observed energies with theoretical energies
        E_obs = self._fblvl['E_obs']*const.h*const.c
        E_th = self._fblvl['E_th']*const.h*const.c
        level_fb = self._fblvl['level']
        use_theoretical = np.logical_and(E_obs==0*u.erg, level_fb!=1)
        E_fb = np.where(use_theoretical, E_th, E_obs)
        # Sum over levels of recombined ion
        sum_factor = u.Quantity(np.zeros(self.temperature.shape + wavelength.shape), 'cm^2')
        for omega, E, n, L, level in zip(self._fblvl['multiplicity'],
                                         E_fb,
                                         self._fblvl['n'],
                                         self._fblvl['L'],
                                         level_fb):
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

    @needs_dataset('verner')
    @u.quantity_input
    def _verner_cross_section(self, energy: u.erg) -> u.cm**2:
        """
        Ground state photoionization cross-section using the method of :cite:t:`verner_analytic_1995`.

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

    @needs_dataset('klgfb')
    @u.quantity_input
    def _karzas_cross_section(self, photon_energy: u.erg, ionization_energy: u.erg, n, l) -> u.cm**2:
        """
        Photoionization cross-section using the method of :cite:t:`karzas_electron_1961`.

        Parameters
        ----------
        photon_energy : `~astropy.units.Quantity`
            Energy of emitted photon
        ionization_energy : `~astropy.units.Quantity`
            Ionization potential of recombined ion for level ``n``
        n : `int`
            Principal quantum number
        l : `int`
            Orbital angular momentum number
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
