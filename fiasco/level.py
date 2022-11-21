"""
Energy level and transitions classes
"""
import astropy.constants as const
import astropy.units as u
import numpy as np

from fiasco.util import vectorize_where

__all__ = ['Level', 'Transitions']


class Level:

    def __init__(self, index, elvlc):
        self._index = index
        self._elvlc = elvlc

    def __repr__(self):
        return f"""Level: {self.level}
Configuration: {self.configuration}
Orbital Angular Momentum: {self.orbital_angular_momentum_label}
Energy: {self.energy.to(u.eV)}"""

    @property
    def level(self):
        return self._elvlc['level'][self._index]

    @property
    def configuration(self):
        return self._elvlc['config'][self._index]

    @property
    def multiplicity(self):
        return self._elvlc['multiplicity'][self._index]

    @property
    def total_angular_momentum(self):
        return self._elvlc['J'][self._index]

    @property
    def orbital_angular_momentum_label(self):
        return self._elvlc['L_label'][self._index]

    @property
    @u.quantity_input
    def energy(self) -> u.erg:
        key = 'E_th' if self._elvlc['E_obs'][self._index] < 0 else 'E_obs'
        return self._elvlc[key][self._index]*const.h*const.c


class Transitions:

    def __init__(self, elvlc, wgfa):
        self._elvlc = elvlc
        self._wgfa = wgfa

    def __len__(self):
        return len(self._wgfa['wavelength'])

    @property
    def is_twophoton(self):
        """
        True if the transition is a two-photon decay
        """
        return self._wgfa['wavelength'] == 0.*u.angstrom

    @property
    def is_observed(self):
        """
        True for transitions that connect two observed energy levels
        """
        return self._wgfa['wavelength'] > 0.*u.angstrom

    @property
    @u.quantity_input
    def A(self) -> u.s**(-1):
        """
        Spontaneous transition probability due to radiative decay
        """
        return self._wgfa['A']

    @property
    @u.quantity_input
    def wavelength(self) -> u.angstrom:
        return np.fabs(self._wgfa['wavelength'])

    @property
    def upper_level(self):
        return self._wgfa['upper_level']

    @property
    def lower_level(self):
        return self._wgfa['lower_level']

    @property
    @u.quantity_input
    def delta_energy(self) -> u.erg:
        energy = u.Quantity(np.where(
            self._elvlc['E_obs'].value == -1, self._elvlc['E_th'].value,
            self._elvlc['E_obs'].value), self._elvlc['E_obs'].unit)
        indices = np.vstack([vectorize_where(self._elvlc['level'], self.lower_level),
                             vectorize_where(self._elvlc['level'], self.upper_level)])
        return np.diff(energy[indices], axis=0).flatten() * const.h * const.c
