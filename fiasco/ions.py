"""
Ion object. Holds all methods and properties of a CHIANTI ion.
"""
import astropy.constants as const
import astropy.table
import astropy.units as u
import numpy as np
import pathlib
import plasmapy.particles
import scipy.special
import warnings

from astropy.utils.data import get_pkg_data_path
from functools import cached_property
from scipy.interpolate import CubicSpline, interp1d, PchipInterpolator

from fiasco import proton_electron_ratio
from fiasco.base import IonBase
from fiasco.collections import IonCollection
from fiasco.gaunt import GauntFactor
from fiasco.levels import Levels, Transitions
from fiasco.util import (
    burgess_tully_descale,
    needs_dataset,
    periodic_table_period,
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
        return f"""CHIANTI Database Ion
---------------------
Name: {self.ion_name}
Element: {self.element_name} ({self.atomic_number})
Charge: +{self.charge_state}
Isoelectronic Sequence: {self.isoelectronic_sequence}
Number of Levels: {self.n_levels}
Number of Transitions: {self.n_transitions}

Temperature range: [{self.temperature[0].to(u.MK):.3f}, {self.temperature[-1].to(u.MK):.3f}]

HDF5 Database: {self.hdf5_dbase_root}
Using Datasets:
  ionization_fraction: {self._dset_names['ionization_fraction']}
  abundance: {self._dset_names['abundance']}
  ip: {self._dset_names['ionization_potential']}"""

    @cached_property
    @needs_dataset('elvlc')
    def levels(self):
        return Levels(self._elvlc)

    def __getitem__(self, key):
        try:
            indexed_levels = self.levels[key]
        except MissingDatasetException:
            raise IndexError(f'No energy levels available for {self.ion_name}')
        else:
            return indexed_levels

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
    def n_levels(self):
        """
        Number of energy levels in the CHIANTI model
        """
        try:
            return len(self.levels)
        except MissingDatasetException:
            return 0

    @property
    def n_transitions(self):
        """
        Number of transitions in the CHIANTI model
        """
        try:
            return len(self.transitions)
        except MissingDatasetException:
            return 0

    @property
    @u.quantity_input
    def thermal_energy(self) -> u.erg:
        """
        Thermal energy, :math:`k_BT`, as a function of temperature.
        """
        return self.temperature.to('erg', equivalencies=u.equivalencies.temperature_energy())

    @cached_property
    @u.quantity_input
    def proton_electron_ratio(self) -> u.dimensionless_unscaled:
        return proton_electron_ratio(self.temperature, **self._instance_kwargs)

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
        return Transitions(self.levels, self._wgfa)

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

    @u.quantity_input
    def level_populations(self,
                          density: u.cm**(-3),
                          include_protons=True,
                          include_level_resolved_rate_correction=True,
                          couple_density_to_temperature=False,
                          use_two_ion_model=True) -> u.dimensionless_unscaled:
        """
        Energy level populations as a function of temperature and density.

        Compute the level populations of the given ion as a function of temperature and
        density. This is done by solving the homogeneous linear system of equations
        describing the processes that populate and depopulate each energy level of each
        ion. Section 3 of :cite:t:`young_chianti_2016` provides a brief description of
        this set of equations.

        Parameters
        ----------
        density: `~astropy.units.Quantity`
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
        use_two_ion_model: `bool`, optional
            If True, include processes that connect the ion to the adjacent ionization stage
            :math:`z+1`. This only makes a difference for CHIANTI database v9 and later. Note
            that this will likely increase the compute time for ions that have a two-ion model.

        Returns
        -------
        `~astropy.units.Quantity`
            A ``(l, m, n)`` shaped quantity, where ``l`` is the number of
            temperatures, ``m`` is the number of densities, and ``n`` is the number of energy
            levels. If ``couple_density_to_temperature=True``, then ``m=1`` and ``l``
            represents the number of temperatures and densities.
        """
        if use_two_ion_model:
            # This avoids running the two-ion model when the data to connect the
            # two ionization stages is not available.
            if not self._has_dataset('auto') and not self._has_dataset('rrlvl'):
                self.log.warning(
                    'No autoionization or level-resolved radiative recombination data '
                    f'available for {self.ion_name}. Using single-ion model for level '
                    'populations calculation. '
                    'To silence this warning, set use_two_ion_model=False.'
                )
                use_two_ion_model = False
        density = np.atleast_1d(density)
        if couple_density_to_temperature:
            if density.shape != self.temperature.shape:
                raise ValueError('Temperature and density must be of equal length if density is '
                                 'coupled to the temperature axis.')
            # NOTE: The reason for this conditional is so that the loop below only
            # performs one iteration and the density value at that one iteration is
            # the entire density array such that density and temperature vary together
            density = [density]
            density_shape = (1,)
        else:
            density_shape = density.shape
        populations = np.zeros(self.temperature.shape + density_shape + (self.n_levels,))
        # Populate density dependent terms and solve matrix equation for each density value
        for i, _d in enumerate(density):
            # NOTE: the following manipulation is needed such that both
            # scalar densities (in the case of no n-T coupling) and arrays
            # (in the case of n-T coupling) can be multiplied by the multi-
            # dimensional excitation rates, whose leading dimension
            # corresponds to the temperature axis.
            d = np.atleast_1d(_d)
            # Compute rate matrix
            if use_two_ion_model:
                c_matrix = self._build_two_ion_coefficient_matrix(d, include_protons=include_protons)
            else:
                c_matrix = self._build_coefficient_matrix(d, include_protons=include_protons)
            # Invert matrix
            c_matrix[:, -1, :] = 1.*c_matrix.unit
            b = np.zeros(c_matrix.shape[2:])
            b[-1] = 1.0
            pop = np.linalg.solve(c_matrix.value, b)
            pop[pop<0] = 0.0
            np.divide(pop, pop.sum(axis=1)[:, np.newaxis], out=pop)
            # Apply ionization/recombination correction
            if include_level_resolved_rate_correction:
                correction = self._population_correction(pop, d, c_matrix)
                pop *= correction
                np.divide(pop, pop.sum(axis=1)[:, np.newaxis], out=pop)
            # NOTE: In cases where the two-ion model is used, this selects only the populations of the
            # recombined ion and renormalizes
            pop = pop[:, :self.n_levels]
            np.divide(pop, pop.sum(axis=1)[:, np.newaxis], out=pop)
            populations[:, i, :] = pop

        return u.Quantity(populations)

    def _build_coefficient_matrix(self, electron_density, include_protons=False):
        d_e = electron_density[:, np.newaxis, np.newaxis]
        rate_matrix = self._rate_matrix_radiative_decay + d_e*self._rate_matrix_collisional_electron
        if include_protons:
            try:
                d_p = self.proton_electron_ratio[:,np.newaxis,np.newaxis]*d_e
                rate_matrix += d_p*self._rate_matrix_collisional_proton
            except MissingDatasetException:
                self.log.warning(
                    f'No proton data available for {self.ion_name}. '
                    'Not including proton excitation and de-excitation in level populations calculation.'
                )
        # Add depopulating terms
        # NOTE: By summing over the rows, we are computing the processes that depopulate
        # that level by summing up all of the processes that populate *from* that level.
        idx = np.diag_indices_from(rate_matrix[0,...])
        rate_matrix[:, idx[0], idx[1]] = -rate_matrix.sum(axis=1)
        return rate_matrix

    def _build_two_ion_coefficient_matrix(self, electron_density, include_protons=False):
        # Get coefficient matrix of recombined ion
        c_matrix_recombined = self._build_coefficient_matrix(electron_density, include_protons=include_protons)
        # Get coefficient matrix of recombining ion
        try:
            c_matrix_recombining = self.next_ion()._build_coefficient_matrix(electron_density,
                                                                             include_protons=include_protons)
        except MissingDatasetException:
            self.log.warning(
                f'No rate data available for recombining ion {self.next_ion().ion_name}. '
                f'Using single-ion model for {self.ion_name}.'
            )
            return c_matrix_recombined
        rate_matrix_total = self._empty_rate_matrix(unit='s-1')
        # Add terms that include both ions
        d_e = electron_density[:, np.newaxis, np.newaxis]
        rate_matrix_total += self._rate_matrix_autoionization
        rate_matrix_total += d_e * (self._rate_matrix_ionization
                                    + self._rate_matrix_radiative_recombination
                                    + self._rate_matrix_dielectronic_capture
                                    + self._rate_matrix_dielectronic_recombination)
        # Add depopulating terms
        # NOTE: By summing over the rows, we are computing the processes that depopulate
        # that level by summing up all of the processes that populate *from* that level.
        idx = np.diag_indices_from(rate_matrix_total[0,...])
        rate_matrix_total[:, idx[0], idx[1]] = -rate_matrix_total.sum(axis=1)
        # Add single-ion terms
        # NOTE: These are added after the two-ion terms because they already include the depopulating
        # terms along the diagonal
        rate_matrix_total[:, :self.n_levels, :self.n_levels] += c_matrix_recombined
        rate_matrix_total[:, self.n_levels:, self.n_levels:] += c_matrix_recombining
        return rate_matrix_total

    @cached_property
    @u.quantity_input
    def _rate_matrix_radiative_decay(self) -> u.Unit('s-1'):
        rate_matrix = u.Quantity(np.zeros((self.n_levels, self.n_levels)), 's-1')
        # Radiative decay into current level from upper levels
        lower_index = self.transitions.lower_level-1
        upper_index = self.transitions.upper_level-1
        rate_matrix[lower_index, upper_index] += self.transitions.A
        return rate_matrix

    @cached_property
    @needs_dataset('scups')
    @u.quantity_input
    def _rate_matrix_collisional_electron(self) -> u.Unit('cm3 s-1'):
        rate_matrix = u.Quantity(np.zeros(self.temperature.shape + (self.n_levels, self.n_levels,)), 'cm3 s-1')
        lower_index = self._scups['lower_level'] - 1
        upper_index = self._scups['upper_level'] - 1
        # De-excitation from upper states
        rate_matrix[:, lower_index, upper_index] += self.electron_collision_deexcitation_rate
        # Excitation from lower states
        rate_matrix[:, upper_index, lower_index] += self.electron_collision_excitation_rate
        return rate_matrix

    @cached_property
    @needs_dataset('psplups')
    @u.quantity_input
    def _rate_matrix_collisional_proton(self) -> u.Unit('cm3 s-1'):
        rate_matrix = u.Quantity(np.zeros(self.temperature.shape + (self.n_levels, self.n_levels,)), 'cm3 s-1')
        lower_index = self._psplups['lower_level'] - 1
        upper_index = self._psplups['upper_level'] - 1
        rate_matrix[:, lower_index, upper_index] += self.proton_collision_deexcitation_rate
        rate_matrix[:, upper_index, lower_index] += self.proton_collision_excitation_rate
        return rate_matrix

    def _empty_rate_matrix(self, temperature_dependent=True, unit='cm3 s-1'):
        n_levels = self.n_levels + self.next_ion().n_levels
        shape = (n_levels, n_levels)
        if temperature_dependent:
            shape = self.temperature.shape + shape
        return u.Quantity(np.zeros(shape), unit)

    @cached_property
    @u.quantity_input
    def _rate_matrix_ionization(self) -> u.Unit('cm3 s-1'):
        rate_matrix = self._empty_rate_matrix()
        # Ionization from ground state of recombined to ground state of recombining
        rate_matrix[:, self.n_levels, 0] = self.ionization_rate
        return rate_matrix

    @cached_property
    @u.quantity_input
    def _rate_matrix_radiative_recombination(self) -> u.Unit('cm3 s-1'):
        rate_matrix = self._empty_rate_matrix()
        try:
            # NOTE: Using copy to avoid in-place modification of cached property
            rr_rate_ground = self.next_ion().radiative_recombination_rate.copy()
        except MissingDatasetException:
            rr_rate_ground = u.Quantity(np.zeros(self.temperature.shape), 'cm3 s-1')
        if self._has_dataset('rrlvl'):
            rr_rate_interp = self._level_resolved_rates_interpolation(self._rrlvl['temperature'],
                                                                      self._rrlvl['rate'],
                                                                      interpolator=interp1d,
                                                                      interpolator_kwargs={'fill_value': np.nan},
                                                                      log_space=True)
            level_final = self._rrlvl['final_level']
            level_initial = self._rrlvl['initial_level']
            # TODO: understand whether we need to sum over repeated level combinations
            rate_matrix[:, level_final-1, level_initial+self.n_levels-1] = rr_rate_interp
            # Subtract total level-resolved rate for the ground level from ground state rate
            # but excluding the ground transition from the sum
            idx_ground = np.logical_and(level_final>1, level_initial==1)
            rr_rate_ground -= rr_rate_interp[:, idx_ground].sum(axis=1)
        # NOTE: The total of the level-resolved rates may sometimes be larger than the ground state rate
        # NOTE: Explicitly setting (and overriding) this value rather than adding to it as this is what
        # is done in the IDL code. See https://github.com/chianti-atomic/chianti-idl/issues/11.
        rate_matrix[:, 0, self.n_levels] = np.where(rr_rate_ground<0.0, 0.0, rr_rate_ground)
        return rate_matrix

    @cached_property
    @u.quantity_input
    def _rate_matrix_autoionization(self) -> u.Unit('s-1'):
        rate_matrix = self._empty_rate_matrix(temperature_dependent=False, unit='s-1')
        # NOTE: Explicitly not using a decorator here in order to return an empty matrix
        # and avoid repeated exception handling later on.
        if not self._has_dataset('auto'):
            self.log.debug(
                f'No .auto data available for {self.ion_name}. '
                'Not including autoionization rates in two-ion rate matrix.'
            )
            return rate_matrix
        # NOTE: Only include those transitions with an upper level below or equal to that of the highest
        # energy level of the recombined ion
        idx = np.where(self._auto['upper_level']<=self.n_levels)
        lower_index = self._auto['lower_level'][idx] + self.n_levels - 1
        upper_index = self._auto['upper_level'][idx] - 1
        rate_matrix[lower_index, upper_index] += self._auto['autoionization_rate'][idx]
        return rate_matrix

    def _dielectronic_capture_rate(self, level_lower, level_upper, A_auto):
        # See Eq. A4 of Dere et al. (2019)
        next_ion = self.next_ion()
        levels_recombined = self[vectorize_where(self.levels.level, level_upper)]
        levels_recombining = next_ion[vectorize_where(next_ion.levels.level, level_lower)]
        g_ratio = levels_recombined.weight / levels_recombining.weight
        delta_energy = levels_recombined.energy - self.ionization_potential - levels_recombining.energy
        kBTE = np.outer(1/self.thermal_energy, delta_energy)
        prefactor = const.h**3 / 2 / (2*np.pi*const.m_e*self.thermal_energy)**(3/2)
        return prefactor[:,np.newaxis] * g_ratio * A_auto * np.exp(-kBTE)

    @cached_property
    @u.quantity_input
    def _rate_matrix_dielectronic_capture(self) -> u.Unit('cm3 s-1'):
        rate_matrix = self._empty_rate_matrix()
        # NOTE: Explicitly not using a decorator here in order to return an empty matrix
        # and avoid repeated exception handling later on.
        if not self._has_dataset('auto'):
            self.log.debug(
                f'No .auto data available for {self.ion_name}. '
                'Not including level-resolved dielectronic capture rates in two-ion rate matrix.'
            )
            return rate_matrix
        # NOTE: Only include those transitions with an upper level below or equal to that of the highest
        # energy level of the recombined ion
        idx = np.where(self._auto['upper_level']<=self.n_levels)
        level_lower = self._auto['lower_level'][idx]
        level_upper = self._auto['upper_level'][idx]
        A_auto = self._auto['autoionization_rate'][idx]
        dc_rate = self._dielectronic_capture_rate(level_lower, level_upper, A_auto)
        upper_index = level_upper - 1
        lower_index = level_lower + self.n_levels - 1
        rate_matrix[:, upper_index, lower_index] += dc_rate
        return rate_matrix

    @cached_property
    @u.quantity_input
    def _rate_matrix_dielectronic_recombination(self) -> u.Unit('cm3 s-1'):
        rate_matrix = self._empty_rate_matrix()
        # Compute ground-ground dielectronic recombination rate
        try:
            dr_rate_ground = self.next_ion().dielectronic_recombination_rate
        except MissingDatasetException:
            dr_rate_ground = u.Quantity(np.zeros(self.temperature.shape), 'cm3 s-1')
        # NOTE: Explicitly not using a decorator here in order to return an empty matrix
        # and avoid repeated exception handling later on.
        if self._has_dataset('auto'):
            # See Eq. A5 of Dere et al. (2019)
            # Compute total of level-resolved dielectronic recombination rates
            # Select only those transitions which autoionize to the ground state of the
            # recombining ion
            idx_ground = np.where(self._auto['lower_level']==1)
            # Sum autoionization rates between upper level and all states
            A_auto_sum = vectorize_where_sum(self._auto['upper_level'],
                                             self._auto['upper_level'][idx_ground],
                                             self._auto['autoionization_rate'])
            # Sum radiative decay rates between upper levels and lower bound levels
            idx_bound = self._wgfa['lower_level'] < self._auto['upper_level'].min()
            A_rad_sum = vectorize_where_sum(self._wgfa['upper_level'][idx_bound],
                                            self._auto['upper_level'][idx_ground],
                                            self._wgfa['A'][idx_bound],)
            branching_ratio = A_rad_sum / (A_rad_sum + A_auto_sum)
            # Get needed levels for recombined and recombining ions
            dc_rate = self._dielectronic_capture_rate(self._auto['lower_level'][idx_ground],
                                                      self._auto['upper_level'][idx_ground],
                                                      self._auto['autoionization_rate'][idx_ground])
            dr_rate_ground -= (dc_rate * branching_ratio).sum(axis=1)
        # NOTE: In some cases, the summed dielectronic capture rates may be larger than the
        # ground-ground dielectronic recombination rates
        rate_matrix[:, 0, self.n_levels] += np.where(dr_rate_ground<0, 0, dr_rate_ground)
        return rate_matrix

    @u.quantity_input
    def _level_resolved_rates_interpolation(self,
                                            temperature_table: u.K,
                                            rate_table: u.cm**3/u.s,
                                            log_space=False,
                                            fill_above=None,
                                            fill_below=None,
                                            interpolator=None,
                                            interpolator_kwargs=None) -> u.cm**3/u.s:
        """
        Extrapolate tables of level-resolved rates as a function of temperature.

        Extrapolate table of level-resolved ionization or recombination rates to
        the temperature array of the ion. Values within the bounds of the temperature
        data are interpolated using ``interpolator``. Values outside of the interpolation
        range are either filled with a scalare value or interpolated using the two points
        on the boundary.

        .. note:: The reason this is a separate function is that in the CHIANTI IDL
                  code the values above and below the temperature data are handled
                  in very particular ways and the level-resolved rates are interpolated
                  multiple times throughout the codebase.

        Parameters
        ----------
        temperature table: `~astropy.units.Quantity`
            Temperature array corresponding to each level-resolved rate
        rate_table: `~astropy.units.Quantity`
            Temperature-dependent, level-resolved rate. The first axis must correspond
            to level and the second axis to temperature.
        log_space: `bool`, optional
            If True, take the base-10 logarithm of the temperature and the rates before
            interpolating. This can be useful when interpolating very small numbers.
        fill_above: `str` or `float`, optional
            If "extrapolate", use the last two points to extrapolate above the temperature
            range. If a `float`, fill in all values above the temperature range using that
            value. If `None` (default), use the rate at the upper temperature boundary as
            the fill value.
        fill_below: `str` or `float`, optional
            If "extrapolate", use the first two points to extrapolate below the temperature
            range. If a `float`, fill in all values below the temperature range using that
            value. If `None` (default), use the rate at the lower temperature boundary as
            the fill value.
        interpolator: callable, optional
            Interpolator to use. By default, this is `~scipy.interpolation.PchipInterpolator`.
        interpolator_kwargs: `dict`, optional
            Keyword arguments to be passed to ``interpolator``.

        Returns
        -------
        rates: `~astropy.units.Quantity`
            Array of rates where the first axis corresponds to temperature and the
            second axis corresponds to level.
        """
        if interpolator is None:
            interpolator = PchipInterpolator
        if interpolator_kwargs is None:
            interpolator_kwargs = {'extrapolate': False}
        # NOTE: According to CHIANTI Technical Report No. 20, Section 5,
        # the interpolation of the level resolved recombination,
        # the rates should be zero below the temperature range and above
        # the temperature range, the last two points should be used to perform
        # a linear extrapolation. For the ionization rates, the rates should be
        # zero above the temperature range and below the temperature range, the
        # last two points should be used. Thus, we need to perform two interpolations
        # for each level.
        # NOTE: In the case of the level-resolved radiative recombination rates in the
        # rrlvl files, we choose to fill the values outside of the interpolation range
        # using the minimum and maximum data values as appropriate.
        # NOTE: In the CHIANTI IDL code, the interpolation is done using a cubic spline.
        # Here, by default, the rates are interpolated using a Piecewise Cubic Hermite Interpolating
        # Polynomial (PCHIP) which balances smoothness and also reduces the oscillations
        # that occur with higher order spline fits. This is needed mostly due to the wide
        # range over which this data is fit.
        temperature = self.temperature.to_value('K')
        temperature_table = temperature_table.to_value('K')
        rate_table = rate_table.to_value('cm3 s-1')
        if log_space:
            temperature = np.log10(temperature)
            temperature_table = np.log10(temperature_table)
            rate_table = np.log10(rate_table)
        rates = []
        for t, r in zip(temperature_table, rate_table):
            # NOTE: Values outside of the temperature data range are set to NaN
            rate_interp = interpolator(t, r, **interpolator_kwargs)(temperature)
            # Extrapolate above temperature range
            f_extrapolate = interp1d(t[-2:],
                                     r[-2:],
                                     kind='linear',
                                     bounds_error=False,
                                     fill_value=r[-1] if fill_above is None else fill_above)
            i_extrapolate = np.where(temperature > t[-1])
            rate_interp[i_extrapolate] = f_extrapolate(temperature[i_extrapolate])
            # Extrapolate below temperature range
            f_extrapolate = interp1d(t[:2],
                                     r[:2],
                                     kind='linear',
                                     bounds_error=False,
                                     fill_value=r[0] if fill_below is None else fill_below)
            i_extrapolate = np.where(temperature < t[0])
            rate_interp[i_extrapolate] = f_extrapolate(temperature[i_extrapolate])
            rates.append(rate_interp)
        if log_space:
            rates = 10**np.array(rates)
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
            fill_below='extrapolate',
            fill_above=0.0,
        )
        return self._cilvl['upper_level'], ionization_rates

    @cached_property
    @needs_dataset('reclvl')
    @u.quantity_input
    def _level_resolved_recombination_rate(self):
        recombination_rates = self._level_resolved_rates_interpolation(
            self._reclvl['temperature'],
            self._reclvl['recombination_rate'],
            fill_below=0.0,
            fill_above='extrapolate',
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
            upper_index_ionization = upper_level_ionization-1
            numerator[:, upper_index_ionization] += (ionization_rate * ionization_fraction_previous).to_value('cm3 s-1')
        except MissingDatasetException:
            pass
        try:
            upper_level_recombination, recombination_rate = self._level_resolved_recombination_rate
            ionization_fraction_next = self.next_ion().ionization_fraction.value[:, np.newaxis]
            upper_index_recombination = upper_level_recombination-1
            numerator[:, upper_index_recombination] += (recombination_rate * ionization_fraction_next).to_value('cm3 s-1')
        except MissingDatasetException:
            pass
        numerator *= density.to_value('cm-3')[:,np.newaxis]

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

           ion.transitions.wavelength[ion.transitions.is_bound_bound]

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
        upper_level = self.transitions.upper_level[self.transitions.is_bound_bound]
        wavelength = self.transitions.wavelength[self.transitions.is_bound_bound]
        A = self.transitions.A[self.transitions.is_bound_bound]
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
    @needs_dataset('diparams')
    @u.quantity_input
    def direct_ionization_rate(self) -> u.cm**3 / u.s:
        r"""
        Ionization rate due to collisions as a function of temperature.

        The ionization rate due to collisions with free electrons assuming a Maxwell-Boltzmann
        distribution. At a minimum, this represents the contribution from the outer-shell electron
        though contributions from inner-shell electrons are also considered for some ions.
        For more details, see the topic guide on :ref:`fiasco-topic-guide-direct-ionization-rate`
        as well as :cite:t:`young_chianti_2025`.
        """
        xgl, wgl = np.polynomial.laguerre.laggauss(12)
        kBT = self.thermal_energy
        cross_section = self._direct_ionization_cross_section(np.outer(xgl, kBT))
        rate_total = u.Quantity(np.zeros(self.temperature.shape), 'cm3 s-1')
        for ip, xs in zip(self._diparams['ip'], cross_section):
            term1 = np.sqrt(8./np.pi/const.m_e)*np.sqrt(kBT)*np.exp(-ip/kBT)
            term2 = ((wgl*xgl)[:, np.newaxis]*xs).sum(axis=0)
            term3 = (wgl[:, np.newaxis]*xs).sum(axis=0)*ip/kBT
            rate_total += term1*(term2 + term3)
        return rate_total

    @needs_dataset('diparams')
    @u.quantity_input
    def _direct_ionization_cross_section(self, energy: u.erg) -> u.cm**2:
        cross_section_all = []
        for ip, bt_c, bt_e, bt_cross_section in zip(self._diparams['ip'],
                                                    self._diparams['bt_c'],
                                                    self._diparams['bt_e'],
                                                    self._diparams['bt_cross_section']):
            U = (energy + ip)/ip
            scaled_energy = 1. - np.log(bt_c)/np.log(U - 1. + bt_c)
            f_interp = PchipInterpolator(bt_e, bt_cross_section, extrapolate=True)
            scaled_cross_section = f_interp(scaled_energy) * bt_cross_section.unit
            cross_section = scaled_cross_section * (np.log(U) + 1.) / U / (ip**2)
            cross_section_all.append(cross_section)

        return u.Quantity(cross_section_all)

    @cached_property
    @needs_dataset('easplups', 'diparams')
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
        # NOTE: Use just the first row as the EA scaling from the diparams files is the same
        # for all rows as it is not related to the number of lines included in the DI calculation.
        # They are contained in this datastructure as a result of the quirk of them being stored in
        # the diparams file in the database.
        scaling = self._diparams['ea'][0][:, np.newaxis]
        # NOTE: The 1/omega multiplicity factor is already included in the scaled upsilon
        # values provided by CHIANTI
        rate = c * scaling * upsilon * np.exp(-1 / kBTE) / np.sqrt(self.thermal_energy)

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

    @u.quantity_input
    def _dielectronic_recombination_suppression(self, density:u.Unit('cm-3'), couple_density_to_temperature=True):
        """
        Density-dependent suppression factor for dielectronic recombination.

        Calculates the density-dependent suppression factor for dielectronic
        recombination following the formulation in Eq. 2 of
        :cite:t:`nikolic_suppression_2018`, including the lower-temperature
        correction given in Eq. 14.

        Parameters
        ----------
        density: `~astropy.units.Quantity`
        """
        if self.isoelectronic_sequence is None:
            return 1
        density = np.atleast_1d(density)
        if couple_density_to_temperature:
            if density.shape != (1,) and density.shape != self.temperature.shape:
                raise ValueError('Temperature and density must be of equal length if density is '
                                 'coupled to the temperature axis.')
        # "A" factor
        A_N = self._nikolic_a_factor()
        q_0 = (1 - np.sqrt(2/3/self.charge_state))*A_N/np.sqrt(self.charge_state)
        T_0 = 5e4*u.K * q_0**2
        # Activation log density (Eq. 3 of Nikolic et al. 2018)
        x_a0 = 10.1821
        x_a = x_a0 + np.log10((self.charge_state/q_0)**7*np.sqrt(self.temperature/T_0))
        # Suppression factor (Eq. 2 of Nikolic et al. 2018)
        width = 5.64586
        x = np.log10(density.to_value('cm-3'))
        if not couple_density_to_temperature:
            x = np.tile(x[:,np.newaxis], self.temperature.shape)
        suppression = np.exp(-((x-x_a)/width*np.sqrt(np.log(2)))**2)
        suppression = np.where(x<=x_a, 1, suppression)
        # Low-temperature correction (Eq. 14 of Nikolic et al. 2018)
        filename = pathlib.Path(get_pkg_data_path('data', package='fiasco')) / 'nikolic_table_5.dat'
        coefficient_table = astropy.table.QTable.read(filename, format='ascii.mrt')
        if self.isoelectronic_sequence in coefficient_table['Sequence']:
            row = coefficient_table[coefficient_table['Sequence']==self.isoelectronic_sequence]
        else:
            # NOTE: if all coefficients are 0, exp_factor evaluates to 1
            row = {f'p_{i}':0*u.eV for i in range(6)}
        # NOTE: Per the footnote to Table 5 of Nikolic et al. (2018), there are two special cases for
        # the p_0 coefficient for H-,He-,and Ne-like ions and for Si-like S III
        if self.ion_name == 'S III':
            row['p_0'] = 17.6874 * u.eV
        if self.isoelectronic_sequence in ['H', 'He', 'Ne']:
            row['p_0'] = 20*scipy.special.erfc(2*(x-x_a0)) * u.eV
            # NOTE: This loop allows for broadcasting later on
            for i in range(1,6):
                row[f'p_{i}'] = row[f'p_{i}']*np.ones(row['p_0'].shape)
        eps_energies = u.Quantity([row[f'p_{i}']*(self.charge_state/10)**i for i in range(6)]).sum(axis=0)
        exp_factor = np.exp(-eps_energies/10/self.thermal_energy)
        return 1 - (1 - suppression)*exp_factor

    def _nikolic_a_factor(self):
        """
        Compute :math:`A(N)` according to Equations 6 and 9 of :cite:t:`nikolic_suppression_2018`.
        """
        Z_iso = plasmapy.particles.atomic_number(self.isoelectronic_sequence)
        # Compute nominal A value according to Eq. 6 and 7 or Table 1
        if Z_iso <= 5:
            # NOTE: According to the paragraph below Eq. 7 of Nikolic et al. (2018), "...the given
            # parameterization was not flexible enough to provide an adequate fit to the
            # Summers (1974 & 1979) data for the lower isoelectronic sequences N<=5.
            # Instead, we explicitly list the optimal values for A(N), for lower ionization
            # stages, in Table 1."
            # NOTE: These values comes from the leftmost columns of Table 1 in Nikolic et al. (2018).
            A_N = {1: 16, 2: 18, 3: 66, 4: 66, 5: 52}[Z_iso]
        else:
            # NOTE: This lookup table comes from Eq. 7 of Nikolic et al. (2018). This is dependent
            # on the "period" (or row on the periodic table) of the isolectronic sequence to which
            # the given ion belongs.
            period_iso = periodic_table_period(self.isoelectronic_sequence)
            N_1, N_2 = {
                2: (3,10), 3: (11,18), 4: (19,36), 5: (37,54), 6: (55,86), 7: (87,118)
            }[period_iso]
            A_N = 12 + 10*N_1 + (10*N_1 - 2*N_2)/(N_1 - N_2)*(Z_iso - N_1)
        # Compute additional modifications according to Eqs. 9, 10, and 11
        filename = pathlib.Path(get_pkg_data_path('data', package='fiasco')) / 'nikolic_table_2.dat'
        coefficient_table = astropy.table.QTable.read(filename, format='ascii.mrt')
        if Z_iso not in coefficient_table['N']:
            return A_N*np.ones(self.temperature.shape)
        # Calculate pis/gammas. Relabel as c_i as the formula is the same
        c_i = []
        for i in range(1,7):
            row = coefficient_table[np.logical_and(coefficient_table['N'] == Z_iso, coefficient_table['i'] == i)]
            c_i.append(
                row['c_1'] + row['c_2']*self.charge_state**row['c_3']*np.exp(-self.charge_state/row['c_4'])
            )
        c_i = np.array(c_i)
        # Calculate psi term According to Eqs. 10 and 11
        logT = np.log10(self.temperature.to_value('K'))
        psi = 1 + c_i[2]*np.exp(-((logT-c_i[0])/np.sqrt(2)/c_i[1])**2) + c_i[5]*np.exp(-((logT-c_i[3])/np.sqrt(2)/c_i[4])**2)
        if Z_iso < 5:
            psi = 2*psi/(1 + np.exp(-2.5e4*u.K*self.charge_state**2/self.temperature))
        return A_N*psi

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

        emission = (prefactor
                  * np.outer(self.thermal_energy**(-3/2), E_photon**5)
                  * exp_energy_ratio
                  / omega_0)
        # NOTE: Necessary because ratio of ionization energy to thermal energy can blow
        # up to infinity at low temperatures for some ionization potentials. Simple multiplication
        # will not sufficiently deal with these as 0*infinity=infinity.
        with warnings.catch_warnings(action='ignore', category=RuntimeWarning):
            emission = np.where(
                np.logical_and(emission==0, np.isinf(sum_factor)), 0, emission*sum_factor
            )
        return emission

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
        gaunt_factor = self.gaunt_factor.free_bound(photon_energy/ionization_energy, n, l)
        cross_section = prefactor * ionization_energy**2 * photon_energy**(-3) * gaunt_factor / n
        cross_section[np.where(photon_energy < ionization_energy)] = 0.*cross_section.unit
        return cross_section

    @needs_dataset('elvlc')
    @u.quantity_input
    def two_photon(self,
                   wavelength: u.angstrom,
                   electron_density: u.cm**(-3),
                   **kwargs) -> u.Unit('erg cm3 s-1 Angstrom-1'):
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
        kwargs : `dict`, optional
            All valid keyword arguments to `level_populations` can also be passed here. Note that in this
            method, proton rates are *not* included by default.
        """
        wavelength = np.atleast_1d(wavelength)
        electron_density = np.atleast_1d(electron_density)

        final_shape = self.temperature.shape + electron_density.shape + wavelength.shape
        couple_density_to_temperature = kwargs.setdefault('couple_density_to_temperature', False)
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
        # NOTE: There are known issues when including the proton rates here for some ions so these
        # are excluded by default. See https://github.com/wtbarnes/fiasco/pull/260#issuecomment-1955770878
        # for more details.
        kwargs.setdefault('include_protons', False)
        level_population = self.level_populations(electron_density, **kwargs)
        level_population = level_population[..., level_index]

        if couple_density_to_temperature:
            electron_density = electron_density[:, np.newaxis]
        level_density = level_population / electron_density
        matrix = np.outer(level_density, energy_dist).reshape(final_shape)

        return const.h * const.c * matrix
