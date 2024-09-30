"""
Energy level and transitions classes
"""
import astropy.constants as const
import astropy.units as u
import numpy as np

from fiasco.util import vectorize_where

__all__ = ['Level', 'Transitions']


class Level:
    """
    An object for holding atomic energy level data from CHIANTI.

    .. warning:: This is not meant to be instantiated directly,
                 but rather accessed by indexing a `fiasco.Ion`
                 object.

    Parameters
    ----------
    index: `int`
        The index of the energy level.
    elvlc: `~fiasco.io.datalayer.DataIndexer`
        Pointer to the energy level information for a given ion in
        the CHIANTI database.
    """

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
        "Index of the level."
        return self._elvlc['level'][self._index]

    @property
    def configuration(self):
        "Label denoting the electronic configuration of the level"
        return self._elvlc['config'][self._index]

    @property
    def multiplicity(self):
        "The multiplicity of the level"
        return self._elvlc['multiplicity'][self._index]

    @property
    def total_angular_momentum(self):
        "The total angular momentum number :math:`J`."
        return self._elvlc['J'][self._index]

    @property
    def orbital_angular_momentum_label(self):
        "The orbital angular momentum number :math:`L`."
        return self._elvlc['L_label'][self._index]

    @property
    def is_observed(self) -> u.erg:
        "True if the energy of the level is from laboratory measurements."
        return self._elvlc['E_obs'][self._index].to_value('cm-1') != -1

    @property
    @u.quantity_input
    def energy(self) -> u.erg:
        """
        Energy of level. Defaults to observed energy and falls back to
        theoretical energy if no measured energy is available.
        """
        if self.is_observed:
            return self._elvlc['E_obs'][self._index].to('erg', equivalencies=u.equivalencies.spectral())
        return self._elvlc['E_th'][self._index].to('erg', equivalencies=u.equivalencies.spectral())


class Transitions:
    """
    An object for holding atomic transition data from CHIANTI.

    .. warning:: This is not meant to be instantiated directly,
                 but rather accessed through `fiasco.Ion.transitions`.

    Parameters
    ----------
    elvlc: `~fiasco.io.datalayer.DataIndexer`
        Pointer to the energy level information for a given ion in
        the CHIANTI database.
    wgfa: `~fiasco.io.datalayer.DataIndexer`
        Pointer to the transition information for a given ion in
        the CHIANTI database.
    """

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
        "Wavelength of each transition."
        return np.fabs(self._wgfa['wavelength'])

    @property
    def upper_level(self):
        "Index of the upper level of the transition."
        return self._wgfa['upper_level']

    @property
    def lower_level(self):
        "Index of the lower level of the transition."
        return self._wgfa['lower_level']

    @property
    @u.quantity_input
    def delta_energy(self) -> u.erg:
        "Energy spacing between the upper and lower level of the transition."
        energy = u.Quantity(np.where(
            self._elvlc['E_obs'].value == -1, self._elvlc['E_th'].value,
            self._elvlc['E_obs'].value), self._elvlc['E_obs'].unit)
        indices = np.vstack([vectorize_where(self._elvlc['level'], self.lower_level),
                             vectorize_where(self._elvlc['level'], self.upper_level)])
        return np.diff(energy[indices], axis=0).flatten() * const.h * const.c
