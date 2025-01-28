"""
Ion object. Holds all methods and properties of a CHIANTI ion.
"""
import astropy.constants as const
import astropy.units as u
import numpy as np

from functools import cached_property
from scipy.interpolate import CubicSpline, interp1d, PchipInterpolator

from fiasco import proton_electron_ratio
from fiasco.base import IonBase
from fiasco.collections import IonCollection
from fiasco.gaunt import GauntFactor
from fiasco.levels import Level, Transitions
from fiasco.util import (
    burgess_tully_descale,
    needs_dataset,
    vectorize_where,
    vectorize_where_sum,
)
from fiasco.util.exceptions import MissingDatasetException

__all__ = ['Ion']


class Ion(IonBase):
    """
    Class for representing a CHIANTI ion.

    The ion object is the fundamental unit of `fiasco`. This object contains
    all of the properties and methods needed to access important information about each ion
    from the CHIANTI database as well as compute common derived quantities.

    Parameters
    ----------
    ion_name : `str` or `tuple`
        Name of the ion. This can be either a string denoting the name or a tuple containing the
        atomic number and ionization stage. See `~fiasco.util.parse_ion_name` for a list of all possible
        input formats.
    temperature : `~astropy.units.Quantity`
        Temperature array over which to evaluate temperature dependent quantities.
    abundance : `str` or `float`, optional
        If a string is provided, use the appropriate abundance dataset.
        If a float is provided, use that value as the abundance.
    ionization_fraction : `str` or `float` or array-like, optional
        If a string is provided, use the appropriate "ioneq" dataset. If an array is provided, it must be the same shape
        as ``temperature``. If a scalar value is passed in, the ionization fraction is assumed constant at all temperatures.
    ionization_potential : `str` or `~astropy.units.Quantity`, optional
        If a string is provided, use the appropriate "ip" dataset.
        If a scalar value is provided, use that value for the ionization potential. This value should be convertible to eV.
    """

    @u.quantity_input
    def __init__(self,
                 ion_name,
                 temperature: u.K,
                 abundance='sun_coronal_1992_feldman_ext',
                 ionization_fraction='chianti',
                 ionization_potential='chianti',
                 *args,
                 **kwargs):
        super().__init__(ion_name, *args, **kwargs)
        self.temperature = np.atleast_1d(temperature)
        self._dset_names = {}
        self.abundance = abundance
        self.ionization_fraction = ionization_fraction
        self.ionization_potential = ionization_potential
        self.gaunt_factor = GauntFactor(hdf5_dbase_root=self.hdf5_dbase_root)

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
Isoelectronic Sequence: {self.isoelectronic_sequence}
Number of Levels: {n_levels}
Number of Transitions: {n_transitions}

Temperature range: [{self.temperature[0].to(u.MK):.3f}, {self.temperature[-1].to(u.MK):.3f}]

HDF5 Database: {self.hdf5_dbase_root}
Using Datasets:
  ionization_fraction: {self._dset_names['ionization_fraction']}
  abundance: {self._dset_names['abundance']}
  ip: {self._dset_names['ionization_potential']}"""

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
        # Keyword arguments used to instantiate this Ion. These are useful when
        # constructing a new Ion instance that pulls from exactly the same
        # data sources.
        kwargs = {
            'hdf5_dbase_root': self.hdf5_dbase_root,
            **self._dset_names,
        }
        # If any of the datasets are set using a string specifying the name of the dataset,
        # the dataset name is in _dset_names. We want to pass this to the new instance
        # so that the new instance knows that the dataset was specified using a
        # dataset name. Otherwise, we can just pass the actual value.
        if kwargs['abundance'] is None:
            kwargs['abundance'] = self.abundance
        if kwargs['ionization_fraction'] is None:
            kwargs['ionization_fraction'] = self.ionization_fraction
        if kwargs['ionization_potential'] is None:
            kwargs['ionization_potential'] = self.ionization_potential
        return kwargs

    def _has_dataset(self, dset_name):
        # There are some cases where we need to check for the existence of a dataset
        # within a function as opposed to checking for the existence of that dataset
        # before entering the function using the decorator approach.
        try:
            needs_dataset(dset_name)(lambda _: None)(self)
        except MissingDatasetException:
            return False
        else:
            return True

    @property
    @u.quantity_input
    def thermal_energy(self) -> u.erg:
        """
        Thermal energy, :math:`k_BT`, as a function of temperature.
        """
        return self.temperature.to('erg', equivalencies=u.equivalencies.temperature_energy())

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
        "A `~fiasco.Transitions` object holding the information about transitions for this ion."
        return Transitions(self._elvlc, self._wgfa)

    @property
    @u.quantity_input
    def ionization_fraction(self) -> u.dimensionless_unscaled:
        """
        Ionization fraction of an ion
        """
        if self._ionization_fraction is None:
            raise MissingDatasetException(
                f"{self._dset_names['ionization_fraction']} ionization fraction data missing for {self.ion_name}"
            )
        return self._ionization_fraction

    @ionization_fraction.setter
    def ionization_fraction(self, ionization_fraction):
        if isinstance(ionization_fraction, str):
            self._dset_names['ionization_fraction'] = ionization_fraction
            ionization_fraction = None
            if self._has_dataset('ion_fraction') and (ionization_fraction := self._ion_fraction.get(self._dset_names['ionization_fraction'])):
                ionization_fraction = self._interpolate_ionization_fraction(
                    self.temperature,
                    ionization_fraction['temperature'],
                    ionization_fraction['ionization_fraction']
                )
            self._ionization_fraction = ionization_fraction
        else:
           # Multiplying by np.ones allows for passing in scalar values
            ionization_fraction = np.atleast_1d(ionization_fraction) * np.ones(self.temperature.shape)
            self._dset_names['ionization_fraction'] = None
            self._ionization_fraction = ionization_fraction

    @staticmethod
    def _interpolate_ionization_fraction(temperature, temperature_data, ionization_data):
        """
        Ionization equilibrium data interpolated to the given temperature

        Interpolated the pre-computed ionization fractions stored in CHIANTI to a new temperature
        array. Returns NaN where interpolation is out of range of the data. For computing
        ionization equilibrium outside of this temperature range, it is better to use the ionization
        and recombination rates.

        .. note::

            The cubic interpolation is performed in log-log space using a Piecewise Cubic Hermite
            Interpolating Polynomial with `~scipy.interpolate.PchipInterpolator`. This helps to
            ensure smoothness while reducing oscillations in the interpolated ionization fractions.

        Parameters
        ----------
        temperature: `~astropy.units.Quantity`
            Temperature array to interpolation onto.
        temperature_data: `~astropy.units.Quantity`
            Temperature array on which the ionization fraction is defined
        ionization_data: `~astropy.units.Quantity`
            Ionization fraction as a function of temperature.

        See Also
        --------
        fiasco.Element.equilibrium_ionization
        """
        temperature = temperature.to_value('K')
        temperature_data = temperature_data.to_value('K')
        ionization_data = ionization_data.to_value()
        # Perform PCHIP interpolation in log-space on only the non-zero ionization fractions.
        # See https://github.com/wtbarnes/fiasco/pull/223 for additional discussion.
        is_nonzero = ionization_data > 0.0
        f_interp = PchipInterpolator(np.log10(temperature_data[is_nonzero]),
                                     np.log10(ionization_data[is_nonzero]),
                                     extrapolate=False)
        ionization_fraction = f_interp(np.log10(temperature))
        ionization_fraction = 10**ionization_fraction
        # This sets all entries that would have interpolated to zero ionization fraction to zero
        ionization_fraction = np.where(np.isnan(ionization_fraction), 0.0, ionization_fraction)
        # Set entries that are truly out of bounds of the original temperature data back to NaN
        out_of_bounds = np.logical_or(temperature<temperature_data.min(), temperature>temperature_data.max())
        ionization_fraction = np.where(out_of_bounds, np.nan, ionization_fraction)
        is_finite = np.isfinite(ionization_fraction)
        ionization_fraction[is_finite] = np.where(ionization_fraction[is_finite] < 0., 0., ionization_fraction[is_finite])
        return u.Quantity(ionization_fraction)

    @property
    @u.quantity_input
    def abundance(self) -> u.dimensionless_unscaled:
        """
        Elemental abundance relative to H.
        """
        if self._abundance is None:
            raise MissingDatasetException(
                f"{self._dset_names['abundance']} abundance data missing for {self.ion_name}"
            )
        return self._abundance

    @abundance.setter
    def abundance(self, abundance):
        """
        Sets the abundance of an ion (relative to H).
        If the abundance is given as a string, use the matching abundance set.
        If the abundance is given as a float, use that value directly.
        """
        self._dset_names['abundance'] = None
        if isinstance(abundance, str):
            self._dset_names['abundance'] = abundance
            abundance = None
            if self._has_dataset('abund'):
                abundance = self._abund.get(self._dset_names['abundance'])
        self._abundance = abundance

    @property
    @u.quantity_input
    def ionization_potential(self) -> u.eV:
        """
        Ionization potential.
        """
        if self._ionization_potential is None:
            raise MissingDatasetException(
                f"{self._dset_names['ionization_potential']} ionization potential data missing for {self.ion_name}"
            )
        # NOTE: Ionization potentials in CHIANTI are stored in units of cm^-1
        # Using this here also means that ionization potentials can be passed
        # in wavenumber units as well.
        return self._ionization_potential.to('eV', equivalencies=u.spectral())

    @ionization_potential.setter
    def ionization_potential(self, ionization_potential):
        """
        Sets the ionization potential of an ion.
        If the ionization potential is given as a string, use the matching ionization potential set.
        if the ionization potential is given as a float, use that value directly.
        """
        self._dset_names['ionization_potential'] = None
        if isinstance(ionization_potential, str):
            self._dset_names['ionization_potential'] = ionization_potential
            ionization_potential = None
            if self._has_dataset('ip'):
                ionization_potential = self._ip.get(self._dset_names['ionization_potential'])
        self._ionization_potential = ionization_potential

    @property
    def hydrogenic(self):
        r"""
        Is the ion in the hydrogen isoelectronic sequence.
        """
        return self.isoelectronic_sequence == 'H'

    @property
    def helium_like(self):
        r"""
        Is the ion in the helium isoelectronic sequence.
        """
        return self.isoelectronic_sequence == 'He'

    @property
    @u.quantity_input
    def formation_temperature(self) -> u.K:
        """
        Temperature at which `~fiasco.Ion.ionization_fraction` is maximum. This is a useful proxy for
        the temperature at which lines for this ion are formed.
        """
        return self.temperature[np.argmax(self.ionization_fraction)]

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
        kBTE = np.outer(self.thermal_energy, 1.0 / self._scups['delta_energy'])
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
        c = const.h**2 / (2. * np.pi * const.m_e)**(1.5)
        upsilon = self.effective_collision_strength
        omega_upper = 2. * self._elvlc['J'][self._scups['upper_level'] - 1] + 1.
        return c * upsilon / np.sqrt(self.thermal_energy[:, np.newaxis]) / omega_upper

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
        kBTE = np.outer(1./self.thermal_energy, self._scups['delta_energy'])
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
        proton_collision_deexcitation_rate
        electron_collision_excitation_rate
        """
        # Create scaled temperature--these are not stored in the file
        bt_t = [np.linspace(0, 1, ups.shape[0]) for ups in self._psplups['bt_rate']]
        # Get excitation rates directly from scaled data
        kBTE = np.outer(self.thermal_energy, 1.0 / self._psplups['delta_energy'])
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
        kBTE = np.outer(self.thermal_energy, 1.0 / self._psplups['delta_energy'])
        omega_upper = 2. * self._elvlc['J'][self._psplups['upper_level'] - 1] + 1.
        omega_lower = 2. * self._elvlc['J'][self._psplups['lower_level'] - 1] + 1.
        dex_rate = (omega_lower / omega_upper) * self.proton_collision_excitation_rate * np.exp(1. / kBTE)

        return dex_rate

    @needs_dataset('elvlc', 'wgfa', 'scups')
    @u.quantity_input
    def level_populations(self,
                          density: u.cm**(-3),
                          include_protons=True,
                          include_level_resolved_rate_correction=True,
                          couple_density_to_temperature=False) -> u.dimensionless_unscaled:
        """
        Energy level populations as a function of temperature and density.

        Compute the level populations of the given ion as a function of temperature and
        density. This is done by solving the homogeneous linear system of equations
        describing the processes that populate and depopulate each energy level of each
        ion. Section 3 of :cite:t:`young_chianti_2016` provides a brief description of
        this set of equations.

        Parameters
        ----------
        density : `~astropy.units.Quantity`
        include_protons : `bool`, optional
            If True (default), include proton excitation and de-excitation rates.
        include_level_resolved_rate_correction: `bool`, optional
            If True (default), include the level-resolved ionization and recombination rate
            correction in the resulting level populations as described in Section 2.3 of
            :cite:t:`landi_chianti-atomic_2006`.
        couple_density_to_temperature: `bool`, optional
            If True, the density will vary along the same axis as temperature
            in the computed level populations and the number of densities must be the same as
            the number of temperatures. This is useful, for example, when computing the level
            populations at constant pressure and is also much faster than computing the level
            populations along an independent density axis. By default, this is set to False.

        Returns
        -------
        `~astropy.units.Quantity`
            A ``(l, m, n)`` shaped quantity, where ``l`` is the number of
            temperatures, ``m`` is the number of densities, and ``n`` is the number of energy
            levels. If ``couple_density_to_temperature=True``, then ``m=1`` and ``l``
            represents the number of temperatures and densities.
        """
        density = np.atleast_1d(density)
        if couple_density_to_temperature:
            if density.shape != self.temperature.shape:
                raise ValueError('Temperature and density must be of equal length if density is '
                                 'coupled to the temperature axis.')

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
            # NOTE: Cannot include protons if psplups data not available for this ion
            try:
                ex_rate_p = self.proton_collision_excitation_rate
                dex_rate_p = self.proton_collision_deexcitation_rate
            except MissingDatasetException:
                self.log.warning(
                    f'No proton data available for {self.ion_name}. '
                    'Not including proton excitation and de-excitation in level populations calculation.')
                include_protons = False
        # NOTE: Having two of the same if blocks here is ugly, but necessary. We cannot continue
        # with the proton rate calculation if the data is not available.
        if include_protons:
            lower_level_p = self._psplups['lower_level']
            upper_level_p = self._psplups['upper_level']
            pe_ratio = proton_electron_ratio(self.temperature, **self._instance_kwargs)
            if couple_density_to_temperature:
                proton_density = (pe_ratio * density)[:, np.newaxis, np.newaxis]
            else:
                proton_density = np.outer(pe_ratio, density)[:, :, np.newaxis]
            ex_diagonal_p = vectorize_where_sum(
                lower_level_p, level, ex_rate_p.value.T, 0).T * ex_rate_p.unit
            dex_diagonal_p = vectorize_where_sum(
                upper_level_p, level, dex_rate_p.value.T, 0).T * dex_rate_p.unit

        # NOTE: The reason for this conditional is so that the loop below only
        # performs one iteration and the density value at that one iteration is
        # the entire density array such that density and temperature vary together
        if couple_density_to_temperature:
            density = [density]
            density_shape = (1,)
        else:
            density_shape = density.shape
        populations = np.zeros(self.temperature.shape + density_shape + (level.max(),))
        # Populate density dependent terms and solve matrix equation for each density value
        for i, _d in enumerate(density):
            c_matrix = coeff_matrix.copy()
            # NOTE: the following manipulation is needed such that both
            # scalar densities (in the case of no n-T coupling) and arrays
            # (in the case of n-T coupling) can be multiplied by the multi-
            # dimensional excitation rates, whose leading dimension
            # corresponds to the temperature axis.
            d = np.atleast_1d(_d)[:, np.newaxis]
            # Collisional excitation and de-excitation out of current state
            c_matrix[:, level-1, level-1] -= d*(ex_diagonal_e + dex_diagonal_e)
            # De-excitation from upper states
            c_matrix[:, lower_level-1, upper_level-1] += d*dex_rate_e
            # Excitation from lower states
            c_matrix[:, upper_level-1, lower_level-1] += d*ex_rate_e
            # Same processes as above, but for protons
            if include_protons:
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
            # Apply ionization/recombination correction
            if include_level_resolved_rate_correction:
                correction = self._population_correction(pop, d, c_matrix)
                pop *= correction
                np.divide(pop, pop.sum(axis=1)[:, np.newaxis], out=pop)
            populations[:, i, :] = pop

        return u.Quantity(populations)

    def _level_resolved_rates_interpolation(self, temperature_table, rate_table,
                                            extrapolate_above=False,
                                            extrapolate_below=False):
        # NOTE: According to CHIANTI Technical Report No. 20, Section 5,
        # the interpolation of the level resolved recombination,
        # the rates should be zero below the temperature range and above
        # the temperature range, the last two points should be used to perform
        # a linear extrapolation. For the ionization rates, the rates should be
        # zero above the temperature range and below the temperature range, the
        # last two points should be used. Thus, we need to perform two interpolations
        # for each level.
        # NOTE: In the CHIANTI IDL code, the interpolation is done using a cubic spline.
        # Here, the rates are interpolated using a Piecewise Cubic Hermite Interpolating
        # Polynomial (PCHIP) which balances smoothness and also reduces the oscillations
        # that occur with higher order spline fits. This is needed mostly due to the wide
        # range over which this data is fit.
        temperature = self.temperature.to_value('K')
        rates = []
        for t, r in zip(temperature_table.to_value('K'), rate_table.to_value('cm3 s-1')):
            rate_interp = PchipInterpolator(t, r, extrapolate=False)(temperature)
            # NOTE: Anything outside of the temperature range will be set to NaN by the
            # interpolation but we want these to be 0.
            rate_interp = np.where(np.isnan(rate_interp), 0, rate_interp)
            if extrapolate_above:
                f_extrapolate = interp1d(t[-2:], r[-2:], kind='linear', fill_value='extrapolate')
                i_extrapolate = np.where(temperature > t[-1])
                rate_interp[i_extrapolate] = f_extrapolate(temperature[i_extrapolate])
            if extrapolate_below:
                f_extrapolate = interp1d(t[:2], r[:2], kind='linear', fill_value='extrapolate')
                i_extrapolate = np.where(temperature < t[0])
                rate_interp[i_extrapolate] = f_extrapolate(temperature[i_extrapolate])
            rates.append(rate_interp)
        # NOTE: Take transpose to maintain consistent ordering of temperature in the leading
        # dimension and levels in the last dimension
        rates =  u.Quantity(rates, 'cm3 s-1').T
        # NOTE: The linear extrapolation at either end may return rates < 0 so we set these
        # to zero.
        rates = np.where(rates<0, 0, rates)
        return rates

    @cached_property
    @needs_dataset('cilvl')
    @u.quantity_input
    def _level_resolved_ionization_rate(self):
        ionization_rates = self._level_resolved_rates_interpolation(
            self._cilvl['temperature'],
            self._cilvl['ionization_rate'],
            extrapolate_below=True,
            extrapolate_above=False,
        )
        return self._cilvl['upper_level'], ionization_rates

    @cached_property
    @needs_dataset('reclvl')
    @u.quantity_input
    def _level_resolved_recombination_rate(self):
        recombination_rates = self._level_resolved_rates_interpolation(
            self._reclvl['temperature'],
            self._reclvl['recombination_rate'],
            extrapolate_below=False,
            extrapolate_above=True,
        )
        return self._reclvl['upper_level'], recombination_rates

    @u.quantity_input
    def _population_correction(self, population, density, rate_matrix):
        """
        Correct level population to account for ionization and
        recombination processes.

        Parameters
        ----------
        population: `np.ndarray`
        density: `~astropy.units.Quantity`
        rate_matrix: `~astropy.units.Quantity`

        Returns
        -------
        correction: `np.ndarray`
            Correction factor to multiply populations by
        """
        # NOTE: These are done in separate try/except blocks because some ions have just a cilvl file,
        # some have just a reclvl file, and some have both.
        # NOTE: Ionization fraction values for surrounding ions are retrieved afterwards because first and last ions do
        # not have previous or next ions but also do not have reclvl or cilvl files.
        # NOTE: stripping the units off and adding them at the end because of some strange astropy
        # Quantity behavior that does not allow for adding these two compatible shapes together.
        numerator = np.zeros(population.shape)
        try:
            upper_level_ionization, ionization_rate = self._level_resolved_ionization_rate
            ionization_fraction_previous = self.previous_ion().ionization_fraction.value[:, np.newaxis]
            numerator[:, upper_level_ionization-1] += (ionization_rate * ionization_fraction_previous).to_value('cm3 s-1')
        except MissingDatasetException:
            pass
        try:
            upper_level_recombination, recombination_rate = self._level_resolved_recombination_rate
            ionization_fraction_next = self.next_ion().ionization_fraction.value[:, np.newaxis]
            numerator[:, upper_level_recombination-1] += (recombination_rate * ionization_fraction_next).to_value('cm3 s-1')
        except MissingDatasetException:
            pass
        numerator *= density.to_value('cm-3')

        c = rate_matrix.to_value('s-1').copy()
        # This excludes processes that depopulate the level
        i_diag, j_diag = np.diag_indices(c.shape[1])
        c[:, i_diag, j_diag] = 0.0
        # Sum of the population-weighted excitations from lower levels
        # and cascades from higher levels
        denominator = np.einsum('ijk,ik->ij', c, population)
        denominator *= self.ionization_fraction.value[:, np.newaxis]
        # Set any zero entries to NaN to avoid divide by zero warnings
        denominator = np.where(denominator==0.0, np.nan, denominator)

        ratio = numerator / denominator
        # Set ratio to zero where denominator is zero. This also covers the
        # case of out-of-bounds ionization fractions (which will be NaN)
        ratio = np.where(np.isfinite(ratio), ratio, 0.0)
        # NOTE: Correction should not affect the ground state populations
        ratio[:, 0] = 0.0
        return 1.0 + ratio

    @needs_dataset('elvlc')
    @u.quantity_input
    def contribution_function(self, density: u.cm**(-3), **kwargs) -> u.cm**3 * u.erg / u.s:
        r"""
        Contribution function :math:`G(n_e,T)` for all transitions.

        The contribution function for ion :math:`k` of element :math:`X` for a
        particular transition :math:`ij` is given by,

        .. math::

           G_{ij} = \mathrm{Ab}(X)f_{X,k}N_jA_{ij}\Delta E_{ij}\frac{1}{n_e},

        Note that the contribution function is often defined in differing ways by different authors.
        The contribution function is defined as above in :cite:t:`young_chianti_2016`.

        The corresponding wavelengths can be retrieved with,

        .. code-block:: python

           ion.transitions.wavelength[~ion.transitions.is_twophoton]

        .. important::

            The ratio :math:`n_H/n_e`, which is often approximated as :math:`n_H/n_e\approx0.83`, is
            explicitly not included here. This means that when computing an intensity with the
            result of this function, the accompanying emission measure is
            :math:`\mathrm{EM}=\mathrm{d}hn_Hn_e` rather than :math:`n_e^2`.

        Parameters
        ----------
        density: `~astropy.units.Quantity`
            Electron number density
        couple_density_to_temperature: `bool`, optional
            If True, the density will vary along the same axis as temperature
            in the computed level populations. The number of densities must be the same as the number of temperatures. This is useful, for
            example, when computing the level populations at constant
            pressure and is also much faster than computing the level
            populations along an independent density axis. By default, this
            is set to False.

        Returns
        -------
        g: `~astropy.units.Quantity`
            A ``(l, m, k)`` shaped quantity, where ``l`` is the number of
            temperatures, ``m`` is the number of densities, and ``k``
            is the number of transitions corresponding to the transition
            wavelengths described above. If ``couple_density_to_temperature=True``,
            then ``m=1`` and ``l`` represents the number of temperatures and densities.

        See Also
        --------
        level_populations
        """
        couple_density_to_temperature = kwargs.get('couple_density_to_temperature', False)
        populations = self.level_populations(density, **kwargs)
        if couple_density_to_temperature:
            term = self.ionization_fraction / density
            term = term[:, np.newaxis, np.newaxis]
        else:
            term = np.outer(self.ionization_fraction, 1./density)
            term = term[:, :, np.newaxis]
        term *= self.abundance
        # Exclude two-photon transitions
        upper_level = self.transitions.upper_level[~self.transitions.is_twophoton]
        wavelength = self.transitions.wavelength[~self.transitions.is_twophoton]
        A = self.transitions.A[~self.transitions.is_twophoton]
        energy = const.h * const.c / wavelength
        i_upper = vectorize_where(self._elvlc['level'], upper_level)
        g = term * populations[:, :, i_upper] * (A * energy)
        return g

    @u.quantity_input
    def emissivity(self, density: u.cm**(-3), **kwargs) -> u.erg * u.cm**(-3) / u.s:
        r"""
        Emissivity as a function of temperature and density for all transitions.

        The emissivity is given by the expression,

        .. math::

           \epsilon(n_e,T) = G(n_e,T)n_Hn_e

        where :math:`G` is the contribution function, :math:`n_H` is the H (or proton) density,
        :math:`n_e` is the electron density, and :math:`T` is the temperature.
        Note that, like the contribution function, emissivity is often defined in
        in differing ways by different authors.
        Here, we use the definition of the emissivity as given by Eq. 3 of
        :cite:t:`young_chianti_2016`.

        .. note::

            The H number density, :math:`n_H`, is computed using ``density`` combined with the
            output of `~fiasco.proton_electron_ratio`.

        Parameters
        ----------
        density : `~astropy.units.Quantity`
            Electron number density.
        couple_density_to_temperature: `bool`, optional
            If True, the density will vary along the same axis as temperature
            in the computed level populations. The number of densities must be the same as the number of temperatures. This is useful, for
            example, when computing the level populations at constant
            pressure and is also much faster than computing the level
            populations along an independent density axis. By default, this
            is set to False.

        Returns
        -------
        `~astropy.units.Quantity`
            A ``(l, m, k)`` shaped quantity, where ``l`` is the number of
            temperatures, ``m`` is the number of densities, and ``k``
            is the number of transitions corresponding to the transition
            wavelengths described in `contribution_function`.
            If ``couple_density_to_temperature=True``, then ``m=1`` and
            ``l`` represents the number of temperatures and densities.

        See Also
        --------
        contribution_function : Calculate contribution function, :math:`G(n,T)`
        """
        density = np.atleast_1d(density)
        pe_ratio = proton_electron_ratio(self.temperature, **self._instance_kwargs)
        pe_ratio = pe_ratio[:, np.newaxis, np.newaxis]
        g = self.contribution_function(density, **kwargs)
        density_squared = density**2
        couple_density_to_temperature = kwargs.get('couple_density_to_temperature', False)
        if couple_density_to_temperature:
            density_squared = density_squared[:, np.newaxis, np.newaxis]
        else:
            density_squared = density_squared[np.newaxis, :, np.newaxis]
        return g * pe_ratio * density_squared

    @u.quantity_input
    def intensity(self,
                  density: u.cm**(-3),
                  emission_measure: u.cm**(-5), **kwargs) -> u.erg / u.cm**2 / u.s / u.steradian:
        r"""
        Line-of-sight intensity computed assuming a particular column emission measure.

        The intensity along the line-of-sight can be written as,

        .. math::

           I = \frac{1}{4\pi}\int\mathrm{d}T,G(n,T)n_Hn_e\frac{dh}{dT}

        which, in the isothermal approximation, can be simplified to,

        .. math::

           I(T_0) \approx \frac{1}{4\pi}G(n,T_0)\mathrm{EM}(T_0)

        where,

        .. math::

           \mathrm{EM}(T) = \int\mathrm{d}h\,n_Hn_e

        is the column emission measure.

        Parameters
        ----------
        density : `~astropy.units.Quantity`
            Electron number density
        emission_measure : `~astropy.units.Quantity`
            Column emission measure. Must be either a scalar, an array of length 1, or
            an array with the same length as ``temperature``. Note that it is assumed
            that the emission measure is the product of the H and electron density.
        couple_density_to_temperature: `bool`, optional
            If True, the density will vary along the same axis as temperature.
            The number of densities must be the same as the number of temperatures.
            This is useful, for example, when computing the intensities at constant
            pressure and is also much faster than computing the intensity
            along an independent density axis. By default, this is set to False.

        Returns
        -------
        `~astropy.units.Quantity`
            A ``(l, m, k)`` shaped quantity, where ``l`` is the number of
            temperatures, ``m`` is the number of densities, and ``k``
            is the number of transitions corresponding to the transition
            wavelengths described in `contribution_function`.
            If ``couple_density_to_temperature=True``, then ``m=1`` and
            ``l`` represents the number of temperatures and densities.
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
        kBT = self.thermal_energy
        energy = np.outer(xgl, kBT) + self.ionization_potential
        cross_section = self.direct_ionization_cross_section(energy)
        term1 = np.sqrt(8./np.pi/const.m_e)*np.sqrt(kBT)*np.exp(-self.ionization_potential/kBT)
        term2 = ((wgl*xgl)[:, np.newaxis]*cross_section).sum(axis=0)
        term3 = (wgl[:, np.newaxis]*cross_section).sum(axis=0)*self.ionization_potential/kBT
        return term1*(term2 + term3)

    @u.quantity_input
    def direct_ionization_cross_section(self, energy: u.erg) -> u.cm**2:
        r"""
        Direct ionization cross-section as a function of energy.

        The direction ionization cross-section is calculated one of two ways.
        See :cite:t:`dere_ionization_2007`, Sections 3.1 and 3.2 for details.
        For H-like ions with :math:`Z\ge6` and He-like ions with :math:`Z\ge10`,
        the cross-section is computed according to the method of
        :cite:t:`fontes_fully_1999`,

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
        if ((self.hydrogenic and self.atomic_number >= 6) or
            (self.helium_like and self.atomic_number >= 10)):
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

    @u.quantity_input
    def _fontes_cross_section(self, energy: u.erg) -> u.cm**2:
        U = energy/self.ionization_potential
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
        return B * (np.pi * const.a0**2) * F * Qrp / (self.ionization_potential.to_value(u.Ry)**2)

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
        c = const.h**2/(2. * np.pi * const.m_e)**(1.5)
        kBTE = np.outer(self.thermal_energy, 1.0/self._easplups['delta_energy'])
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
        rate = c * upsilon * np.exp(-1 / kBTE) / np.sqrt(self.thermal_energy)

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
    @needs_dataset('trparams')
    @u.quantity_input
    def _total_recombination_rate(self) -> u.cm**3 / u.s:
        temperature_data = self._trparams['temperature'].to_value('K')
        rate_data = self._trparams['recombination_rate'].to_value('cm3 s-1')
        f_interp = interp1d(temperature_data, rate_data, fill_value='extrapolate', kind='cubic')
        f_interp = PchipInterpolator(np.log10(temperature_data), np.log10(rate_data), extrapolate=True)
        rate_interp = 10**f_interp(np.log10(self.temperature.to_value('K')))
        return u.Quantity(rate_interp, 'cm3 s-1')

    @cached_property
    @u.quantity_input
    def recombination_rate(self) -> u.cm**3 / u.s:
        r"""
        Total recombination rate as a function of temperature.

        The total recombination rate, as a function of temperature, for a given ion
        is the sum of the radiative and dielectronic recombination rates such that,

        .. math::

            \alpha_{R} = \alpha_{RR} + \alpha_{DR}

        .. important::

            For most ions, this total recombination rate is computed by summing the
            outputs of the `radiative_recombination_rate` and `dielectronic_recombination_rate` methods.
            However, for some ions, total recombination rate data is available in the
            so-called ``.trparams`` files. For these ions, the output of this method
            will *not* be equal to the sum of the `dielectronic_recombination_rate` and
            `radiative_recombination_rate` method. As such, when computing the total
            recombination rate, this method should always be used.

        See Also
        --------
        radiative_recombination_rate
        dielectronic_recombination_rate
        """
        # NOTE: If the trparams data is available, then it is prioritized over the sum
        # of the dielectronic and radiative recombination rates. This is also how the
        # total recombination rates are computed in IDL. The reasoning here is that the
        # total recombination rate data, if available, is more reliable than the sum of
        # the radiative and dielectronic recombination rates. According to P. Young, there
        # is some slight controversy over this within some communities, but CHIANTI has chosen
        # to prioritize this data if it exists.
        try:
            tr_rate = self._total_recombination_rate
        except MissingDatasetException:
            self.log.debug(f'No total recombination data available for {self.ion_name}.')
        else:
            return tr_rate
        try:
            rr_rate = self.radiative_recombination_rate
        except MissingDatasetException:
            self.log.debug(f'No radiative recombination data available for {self.ion_name}.')
            rr_rate = u.Quantity(np.zeros(self.temperature.shape), 'cm3 s-1')
        try:
            dr_rate = self.dielectronic_recombination_rate
        except MissingDatasetException:
            self.log.debug(f'No dielectronic recombination data available for {self.ion_name}.')
            dr_rate = u.Quantity(np.zeros(self.temperature.shape), 'cm3 s-1')
        return rr_rate + dr_rate

    @u.quantity_input
    def free_free(self, wavelength: u.angstrom) -> u.erg * u.cm**3 / u.s / u.angstrom:
        r"""
        Free-free continuum emission as a function of temperature and wavelength.

        .. important::

            This does not include ionization fraction or abundance factors.

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

        See Also
        --------
        fiasco.GauntFactor.free_free: Calculation of :math:`\langle g_{ff}\rangle`.
        fiasco.IonCollection.free_free: Includes abundance and ionization equilibrium.
        """
        prefactor = (const.c / 3. / const.m_e * (const.alpha * const.h / np.pi)**3
                     * np.sqrt(2. * np.pi / 3. / const.m_e))
        tmp = np.outer(self.thermal_energy, wavelength)
        exp_factor = np.exp(-const.h * const.c / tmp) / (wavelength**2)
        gf = self.gaunt_factor.free_free(self.temperature, wavelength, self.atomic_number, self.charge_state, )

        return (prefactor * self.charge_state**2 * exp_factor * gf
                / np.sqrt(self.thermal_energy)[:, np.newaxis])

    @u.quantity_input
    def free_free_radiative_loss(self, use_itoh=False) -> u.erg * u.cm**3 / u.s:
        r"""
        Free-free continuum radiative losses as a function of temperature.

        .. important::

            This does not include the ionization fraction or abundance factors.

        The total free-free radiative loss is given by integrating the emissivity over
        all wavelengths.  The total losses per unit emission measure are then given by
        Equation 18 of :cite:`sutherland_accurate_1998`,

        .. math::

            R_{ff}(T_e) = F_{k} \sqrt{(T_{e})} z^{2} \langle g_{t,ff}\rangle

        where :math:`T_{e}` is the electron temperature, :math:`F_{k}` is a constant,
        :math:`z` is the charge state, and :math:`\langle g_{t,ff}\rangle` is the
        wavelength-integrated free-free Gaunt factor.
        The prefactor :math:`F_{k}` is defined in Equation 19 of :cite:t:`sutherland_accurate_1998`,

        .. math::

            F_k =& \frac{16e^6}{3^{3/2}c^3}\sqrt{\frac{2\pi k_B}{\hbar^2m_e^3}}\\
                \approx& 1.42555669\times10^{-27}\,\mathrm{cm}^{5}\,\mathrm{g}\,\mathrm{K}^{-1/2}\,\mathrm{s}^{3}.

        Parameters
        ----------
        use_itoh : `bool`, optional
            Whether to use Gaunt factors taken from :cite:t:`itoh_radiative_2002`.
            Defaults to false.

        See Also
        --------
        fiasco.GauntFactor.free_free_integrated: Calculation of :math:`\langle g_{t,ff}\rangle`.
        """
        prefactor = (16./3**1.5) * np.sqrt(2. * np.pi / (const.hbar**2 * const.m_e**3)) * (const.e.esu**6 / const.c**3)
        gf = self.gaunt_factor.free_free_integrated(self.temperature, self.charge_state, use_itoh=use_itoh)
        return (prefactor * self.charge_state**2 * gf * np.sqrt(self.thermal_energy))

    @needs_dataset('fblvl', 'ip')
    @u.quantity_input
    def free_bound(self,
                   wavelength: u.angstrom,
                   use_verner=True) -> u.Unit('erg cm3 s-1 Angstrom-1'):
        r"""
        Free-bound continuum emission of the recombined ion.

        .. important::

            This does not include the ionization fraction or abundance factors.

        .. important::

            Unlike the equivalent IDL routine, the output here is not
            expressed per steradian and as such the factor of :math:`1/4\pi` is not included.

        When an electron is captured by an ion of charge :math:`z+1`
        (the recombining ion), it creates a an ion of charge :math:`z`
        (the recombined ion) and produces a continuum of emission
        called the free-bound continuum. The emission of the
        recombined ion is given by,

        .. math::

            C_{fb}(\lambda, T) = \frac{2}{hc^3(k_B m_e)^{3/2}\sqrt{2\pi}}\frac{E^5}{T^{3/2}}\sum_i\frac{\omega_i}{\omega_0}\sigma_i^{\mathrm{bf}}\exp{\left(-\frac{E-I_i}{k_BT}\right)}

        where :math:`E` is the energy of the outgoing photon,
        :math:`\omega_i,\omega_0` are the statistical weights of the
        :math:`i`-th level of the recombined ion and the ground level of the recombining ion, respectively,
        :math:`\sigma_i^{\mathrm{bf}}` is the free-bound cross-section,
        and :math:`I_i` is the energy required to ionize the recombined ion from level :math:`i`.
        A detailed derivation of this formula can be found in :cite:t:`young_chianti_2021`.

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
        prefactor = 2/np.sqrt(2*np.pi)/(const.h*(const.c**3) * const.m_e**(3/2))
        recombining = self.next_ion()
        omega_0 = recombining._fblvl['multiplicity'][0] if recombining._has_dataset('fblvl') else 1.0
        E_photon = const.h * const.c / wavelength
        # Precompute this here to avoid repeated outer product calculations
        exp_energy_ratio = np.exp(-np.outer(1/self.thermal_energy, E_photon))
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
            E_ionize = self.ionization_potential - E
            # Check if ionization potential and photon energy sufficient
            if (E_ionize < 0*u.erg) or (E_photon.max() < E):
                continue
            # Only use Verner cross-section for ground state
            if level == 1 and use_verner:
                cross_section = self._verner_cross_section(E_photon)
            else:
                cross_section = self._karzas_cross_section(E_photon, E_ionize, n, L)
            # NOTE: Scaled energy can blow up at low temperatures such that taking an
            # exponential yields numbers too high to be expressed with double precision.
            # At these temperatures, the cross-section is 0 anyway so we can just zero
            # these terms. Just multiplying by 0 is not sufficient because 0*inf=inf
            with np.errstate(over='ignore', invalid='ignore'):
                exp_ip_ratio = np.exp(E_ionize/self.thermal_energy)
                xs_exp_ip_ratio = np.outer(exp_ip_ratio, cross_section)
            xs_exp_ip_ratio[:,cross_section==0.0*u.cm**2] = 0.0 * u.cm**2
            sum_factor += omega * xs_exp_ip_ratio

        return (prefactor
                * np.outer(self.thermal_energy**(-3/2), E_photon**5)
                * exp_energy_ratio
                * sum_factor / omega_0)

    @u.quantity_input
    def free_bound_radiative_loss(self) -> u.erg * u.cm**3 / u.s:
        r"""
        The radiative loss rate for free-bound emission as a function of temperature,
        integrated over all wavelengths.

        .. important::

            This does not include the ionization fraction or abundance factors.

        .. note::

            This ion, for which the free-bound radiative loss is being calculated,
            is taken to be the recombining ion. The ion one ionization stage lower
            is taken to be the recombined ion.

        The calculation integrates Equation 1a of :cite:t:`mewe_calculated_1986`, where the
        Gaunt factor is summed only for free-bound emission :cite:p:`young_chianti_2019-1`.
        Since the form of the Gaunt factor used by :cite:t:`mewe_calculated_1986` does not
        depend on wavelength, the integral is straightforward.

        The continuum intensity per unit emission measure is given by:

        .. math::

            C_{fb}(\lambda, T) = \frac{F g_{fb}}{\lambda^{2}\ T^{1/2}} \exp{\Big(\frac{-h c}{\lambda k_{B} T}\Big)}

        where

        .. math::

            F = \frac{64 \pi}{3} \sqrt{\frac{\pi}{6}} \frac{q_{e}^{6}}{c^{2} m_{e}^{2} k_{B}^{1/2}}

        is a constant :cite:p:`gronenschild_calculated_1978`.
        Integrating in wavelength space gives the free-bound loss rate,

        .. math::

            R_{fb} = \frac{F k_{B} g_{fb} T^{1/2}}{h c} \exp{\Big(\frac{-h c}{\lambda k_{B} T}\Big)}

        We have dropped the factor of :math:`n_{e}^{2}` here to make the loss rate per unit emission measure.

        .. note:: The form of :math:`C_{fb}` used by :cite:t:`mewe_calculated_1986` and given above is slightly
                  different than the form used in `~fiasco.Ion.free_bound` and as such the two approaches are
                  not entirely self-consistent. This particular form is used, rather than calling
                  `~fiasco.Ion.free_bound` and integrating the result, for the sake of efficiency.

        See Also
        --------
        fiasco.GauntFactor.free_bound_integrated: Calculation of :math:`g_{fb}`
        """
        if self.charge_state == 0:
            return u.Quantity(np.zeros(self.temperature.shape) * u.erg * u.cm**3 / u.s)

        recombined = self.previous_ion()
        if not recombined._has_dataset('fblvl'):
            return u.Quantity(np.zeros(self.temperature.shape) * u.erg * u.cm**3 / u.s)
        C_ff = 64 * np.pi / 3.0 * np.sqrt(np.pi/6.) * (const.e.esu**6)/(const.c**2 * const.m_e**1.5)
        prefactor = C_ff * np.sqrt(self.thermal_energy) / (const.h*const.c)

        E_obs = recombined._fblvl['E_obs']*const.h*const.c
        E_th = recombined._fblvl['E_th']*const.h*const.c
        n0 = recombined._fblvl['n'][0]
        E_fb = np.where(E_obs==0*u.erg, E_th, E_obs)
        wvl_n0 = 1 / (recombined.ionization_potential - E_fb[0]).to('cm-1', equivalencies=u.spectral())
        wvl_n1 = (n0 + 1)**2 /(const.Ryd * self.charge_state**2)
        g_fb0 = self.gaunt_factor.free_bound_integrated(self.temperature,
                                                        self.atomic_number,
                                                        self.charge_state,
                                                        n0,
                                                        recombined.ionization_potential,
                                                        ground_state=True)
        g_fb1 = self.gaunt_factor.free_bound_integrated(self.temperature,
                                                        self.atomic_number,
                                                        self.charge_state,
                                                        n0,
                                                        recombined.ionization_potential,
                                                        ground_state=False)
        term1 = g_fb0 * np.exp(-const.h*const.c/(self.thermal_energy * wvl_n0))
        term2 = g_fb1 * np.exp(-const.h*const.c/(self.thermal_energy * wvl_n1))

        return prefactor * (term1 + term2)

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
        E_scaled = np.log(photon_energy/ionization_energy)
        gaunt_factor = self.gaunt_factor.free_bound(E_scaled, n, l)
        cross_section = prefactor * ionization_energy**2 * photon_energy**(-3) * gaunt_factor / n
        cross_section[np.where(photon_energy < ionization_energy)] = 0.*cross_section.unit
        return cross_section

    @needs_dataset('elvlc')
    @u.quantity_input
    def two_photon(self,
                   wavelength: u.angstrom,
                   electron_density: u.cm**(-3),
                   include_protons=False,
                   couple_density_to_temperature=False) -> u.Unit('erg cm3 s-1 Angstrom-1'):
        r"""
        Two-photon continuum emission of a hydrogenic or helium-like ion.

        .. important::

            This does not include the ionization fraction or abundance factors.

        .. important::

            Unlike the equivalent IDL routine, the output here is not
            expressed per steradian and as such the factor of :math:`1/4\pi` is not included.

        For more details regarding this calculation, see :ref:`fiasco-topic-guide-two-photon-continuum`.

        Parameters
        ----------
        wavelength : `~astropy.units.Quantity`
        electron_density : `~astropy.units.Quantity`
        include_protons : `bool`, optional
            If True, use proton excitation and de-excitation rates in the level population calculation.
        couple_density_to_temperature: `bool`, optional
            If True, the density will vary along the same axis as temperature
        """
        wavelength = np.atleast_1d(wavelength)
        electron_density = np.atleast_1d(electron_density)

        final_shape = self.temperature.shape + electron_density.shape + wavelength.shape
        if couple_density_to_temperature:
            final_shape = self.temperature.shape + (1,) + wavelength.shape
        if self.hydrogenic:
            A_ji = self._hseq['A']
            psi_norm = self._hseq['psi_norm']
            x_interp, y_interp = self._hseq['y'], self._hseq['psi']
            config = '2s'  # Get the index of the 2S1/2 state for H-like
            J = 0.5
        elif self.helium_like:
            A_ji = self._heseq['A']
            psi_norm = 1.0 * u.dimensionless_unscaled
            x_interp, y_interp = self._heseq['y'], self._heseq['psi']
            config = '1s.2s'  # Get the index of the 1s2s 1S0 state for He-like:
            J = 0
        else:
            return u.Quantity(np.zeros(final_shape),  'erg cm^3 s^-1 Angstrom^-1')
        level_index = np.where((self._elvlc['config'] == config) & (np.isclose(self._elvlc['J'], J)) )[0][0]

        E_obs = self._elvlc['E_obs'][level_index]
        E_th = self._elvlc['E_th'][level_index]
        E_2p = E_obs if E_obs > 0.0 else E_th
        rest_wavelength = 1 / E_2p

        # NOTE: Explicitly setting the boundary condition type here to match the behavior of the
        # IDL spline interpolation functions. See https://github.com/wtbarnes/fiasco/pull/297 for
        # additional details.
        cubic_spline = CubicSpline(x_interp, y_interp, bc_type='natural')
        x_new = (rest_wavelength / wavelength).decompose().to_value(u.dimensionless_unscaled)
        psi_interp = cubic_spline(x_new)
        psi_interp = np.where(x_new>1.0, 0.0, psi_interp)

        energy_dist = (A_ji * rest_wavelength * psi_interp) / (psi_norm * wavelength**3)

        level_population = self.level_populations(electron_density, include_protons=include_protons)
        level_population = level_population[..., level_index]

        if couple_density_to_temperature:
            electron_density = electron_density[:, np.newaxis]
        level_density = level_population / electron_density
        matrix = np.outer(level_density, energy_dist).reshape(final_shape)

        return const.h * const.c * matrix
