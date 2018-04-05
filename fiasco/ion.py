"""
Ion object. Holds all methods and properties of a CHIANTI ion.
"""
import warnings

import numpy as np
from scipy.interpolate import splrep, splev, interp1d
import astropy.units as u
import astropy.constants as const

from .base import IonBase, ContinuumBase
from .collections import IonCollection
from fiasco.util import needs_dataset


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
        self.temperature = temperature
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
        _ = self._elvlc['level'][key]
        return Level(key, self._elvlc)
    
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
                     kind='linear', bounds_error=False, fill_value=np.nan)
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
    def ip(self):
        """
        Ionization potential with reasonable units
        """
        if self._ip is not None:
            return (self._ip[self._dset_names['ip_filename']]
                    * const.h.cgs * const.c.cgs).decompose().cgs
        else:
            return None

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

    def __add__(self, value):
        return IonCollection(self, value)

    def __radd__(self, value):
        return IonCollection(value, self)
        
    @staticmethod
    def burgess_tully_descale(x, y, energy_ratio, c, scaling_type):
        """
        Convert scaled Burgess-Tully parameters to physical quantities. For more details see
        [1]_.

        Parameters
        ----------
        x : `~astropy.units.Quantity`
        y : `~astropy.units.Quantity`
        energy_ratio : `~astropy.units.Quantity`
            Ratio of temperature to photon energy
        c : `~astropy.units.Quantity`
            Scaling constant
        scaling_type : `int`

        Returns
        -------
        upsilon : `~numpy.NDArray`
            Descaled collision strength or cross-section

        References
        ----------
        .. [1] Burgess, A. and Tully, J. A., 1992, A&A, `254, 436 <http://adsabs.harvard.edu/abs/1992A%26A...254..436B>`_ 
        """
        nots = splrep(x, y, s=0)
        if scaling_type == 1:
            x_new = 1.0 - np.log(c) / np.log(energy_ratio + c)
            upsilon = splev(x_new, nots, der=0) * np.log(energy_ratio + np.e)
        elif scaling_type == 2:
            x_new = energy_ratio / (energy_ratio + c)
            upsilon = splev(x_new, nots, der=0)
        elif scaling_type == 3:
            x_new = energy_ratio / (energy_ratio + c)
            upsilon = splev(x_new, nots, der=0) / (energy_ratio + 1.0)
        elif scaling_type == 4:
            x_new = 1.0 - np.log(c) / np.log(energy_ratio + c)
            upsilon = splev(x_new, nots, der=0) * np.log(energy_ratio + c)
        elif scaling_type == 5:
            # dielectronic
            x_new = energy_ratio / (energy_ratio + c)
            upsilon = splev(x_new, nots, der=0) / energy_ratio
        elif scaling_type == 6:
            # protons
            x_new = energy_ratio / (energy_ratio + c)
            upsilon = 10**splev(x_new, nots, der=0)
        else:
            raise ValueError('Unrecognized BT92 scaling option.')

        return upsilon

    @needs_dataset('scups')
    def effective_collision_strength(self):
        """
        Maxwellian-averaged collision strength, typically denoted by :math:`\\Upsilon`
        
        Note
        ----
        Need a more efficient way of calculating upsilon for all transitions. Current
        method is slow for ions with many transitions, e.g. Fe IX and Fe XI

        See Also
        --------
        burgess_tully_descale : Descale and interpolate :math:`\\Upsilon`
        """
        energy_ratio = np.outer(const.k_B.cgs*self.temperature,
                                1.0/self._scups['delta_energy'].to(u.erg))
        upsilon = np.array(list(map(self.burgess_tully_descale, self._scups['bt_t'],
                                    self._scups['bt_upsilon'], energy_ratio.T, self._scups['bt_c'],
                                    self._scups['bt_type'])))
        upsilon = u.Quantity(np.where(upsilon > 0., upsilon, 0.))
        return upsilon.T

    @needs_dataset('elvlc', 'scups')
    def electron_collision_deexcitation_rate(self):
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

        Reference
        ---------
        .. [1] Phillips, K., et al., 2008, `Ultraviolet and X-ray Spectroscopy of the Solar Atmosphere <http://adsabs.harvard.edu/abs/2008uxss.book.....P>`_

        See Also
        --------
        electron_collision_excitation_rate : Excitation rate due to collisions
        effective_collision_strength : Maxwellian-averaged collision strength, :math:`\\Upsilon`
        """
        c = (const.h.cgs**2)/((2. * np.pi * const.m_e.cgs)**(1.5) * np.sqrt(const.k_B.cgs))
        upsilon = self.effective_collision_strength()
        omega_upper = 2. * self._elvlc['J'][self._scups['upper_level'] - 1] + 1.
        return c * upsilon / np.sqrt(self.temperature[:, np.newaxis]) / omega_upper

    @needs_dataset('elvlc', 'scups')
    def electron_collision_excitation_rate(self):
        """
        Collisional excitation rate coefficient for electrons.

        See Also
        --------
        electron_collision_deexcitation_rate : De-excitation rate due to collisions
        """
        dex_rate = self.electron_collision_deexcitation_rate()
        omega_upper = 2. * self._elvlc['J'][self._scups['upper_level'] - 1] + 1.
        omega_lower = 2. * self._elvlc['J'][self._scups['lower_level'] - 1] + 1.
        energy_ratio = np.outer(1./const.k_B.cgs/self.temperature,
                                self._scups['delta_energy'].to(u.erg))
        return omega_upper / omega_lower * dex_rate * np.exp(-energy_ratio)

    @needs_dataset('ip')
    def direct_ionization_rate(self):
        """
        Calculate direct ionization rate in cm3/s
        
        Needs an equation reference or explanation
        """
        xgl, wgl = np.polynomial.laguerre.laggauss(12)
        kBT = const.k_B.cgs*self.temperature
        energy = np.outer(xgl, kBT) * kBT.unit + self.ip
        cross_section = self.direct_ionization_cross_section(energy)
        if cross_section is None:
            return None
        term1 = np.sqrt(8./np.pi/const.m_e.cgs)*np.sqrt(kBT)*np.exp(-self.ip/kBT)
        term2 = ((wgl*xgl)[:, np.newaxis]*cross_section).sum(axis=0)
        term3 = (wgl[:, np.newaxis]*cross_section).sum(axis=0)*self.ip/kBT
        
        return term1*(term2 + term3)
    
    @u.quantity_input
    def direct_ionization_cross_section(self, energy: u.erg):
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
    def _dere_cross_section(self, energy: u.erg):
        """
        Calculate direct ionization cross-sections according to [1]_.

        References
        ----------
        .. [1] Dere, K. P., 2007, A&A, `466, 771 <http://adsabs.harvard.edu/abs/2007A%26A...466..771D>`_
        """
        # Cross-sections from diparams file
        cross_section_total = np.zeros(energy.shape)
        for ip, bt_c, bt_e, bt_cross_section in zip(self._diparams['ip'], self._diparams['bt_c'],
                                                    self._diparams['bt_e'],
                                                    self._diparams['bt_cross_section']):
            U = energy/(ip.to(u.erg))
            scaled_energy = 1. - np.log(bt_c)/np.log(U - 1. + bt_c)
            f_interp = interp1d(bt_e.value, bt_cross_section.value, kind='cubic',
                                fill_value='extrapolate')
            scaled_cross_section = f_interp(scaled_energy.value)*bt_cross_section.unit
            # Only nonzero at energies above the ionization potential
            scaled_cross_section *= (U.value > 1.)
            cross_section = scaled_cross_section * (np.log(U) + 1.) / U / (ip**2)
            if not hasattr(cross_section_total, 'unit'):
                cross_section_total = cross_section_total*cross_section.unit
            cross_section_total += cross_section
            
        return cross_section_total

    @needs_dataset('ip')
    @u.quantity_input
    def _fontes_cross_section(self, energy: u.erg):
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
        
        return B * (np.pi * const.a0.cgs**2) * F * Qrp / (self.ip.to(u.Ry).value**2)

    @needs_dataset('easplups')
    def excitation_autoionization_rate(self):
        """
        Calculate ionization rate due to excitation autoionization
        """
        # Collision constant
        c = (const.h.cgs**2)/((2. * np.pi * const.m_e.cgs)**(1.5) * np.sqrt(const.k_B.cgs))
        kBTE = u.Quantity(np.outer(
            const.k_B.cgs*self.temperature, 1.0/self._easplups['delta_energy'].to(u.erg)))
        # Descale upsilon
        shape = self._easplups['bt_upsilon'].shape
        xs = np.tile(np.linspace(0, 1, shape[1]), shape[0]).reshape(shape)
        args = [xs, self._easplups['bt_upsilon'].value, kBTE.value, self._easplups['bt_c'].value, 
                self._easplups['bt_type']]
        upsilon = u.Quantity(list(map(self.burgess_tully_descale, *args)))
        # Rate coefficient
        rate = c * upsilon * np.exp(-1 / kBTE) / np.sqrt(self.temperature[np.newaxis, :])
        
        return rate.sum(axis=0)
    
    def ionization_rate(self):
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
    def radiative_recombination_rate(self):
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
    def dielectronic_recombination_rate(self):
        """
        Dielectronic recombination rate

        Calculated according to one of two methods,

        - Method of [1]_, (show expression)
        - Method of [2]_, (show expression)

        References
        ----------
        """
        if self._drparams['fit_type'][0] == 1:
            E_over_T = (np.outer(self._drparams['E_fit'], 1./self.temperature)
                        * (self._drparams['E_fit'].unit/self.temperature.unit))
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
    
    def recombination_rate(self):
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


class Level(object):

    def __init__(self, index, elvlc):
        self._index = index
        self._elvlc = elvlc

    def __repr__(self):
        return f"""Level: {self.level}
Configuration: {self.configuration}
Orbital Angular Momentum: {self.orbital_angular_momentum_label}
Energy: {self.energy}"""

    @property
    def level(self):
        return self._elvlc['level'][self._index]

    @property
    def configuration(self):
        return self._elvlc['config'][self._index]

    @property
    def multiplicity(self):
        return self._elvlc['mult'][self._index]

    @property
    def total_angular_momentum(self):
        return self._elvlc['J'][self._index]

    @property
    def orbital_angular_momentum_label(self):
        return self._elvlc['L_label'][self._index]

    @property
    def energy(self):
        if self._elvlc['E_obs'][self._index] < 0:
            return (self._elvlc['E_th'][self._index]*const.h.cgs*const.c.cgs).decompose().cgs
        else:
            return (self._elvlc['E_obs'][self._index]*const.h.cgs*const.c.cgs).decompose().cgs
