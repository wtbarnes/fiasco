"""
Energy level and transitions classes
"""
import astropy.units as u
import numpy as np

from fiasco.util import vectorize_where

__all__ = ['Levels', 'Transitions']


_ORBITAL_ANGULAR_MOMENTUM_LABELS = [
    'S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U','V','W','X','Y'
]
_ORBITAL_ANGULAR_MOMENTUM_LOOKUP = {
    l:i for i,l in enumerate(_ORBITAL_ANGULAR_MOMENTUM_LABELS)
}


class Levels:
    """
    Data associated with all energy levels of an ion.

    .. warning:: This is not meant to be instantiated directly,
                 but rather accessed by indexing a `fiasco.Ion`
                 object.

    Parameters
    ----------
    elvlc: `~fiasco.io.datalayer.DataIndexer`
        Pointer to the energy level information for a given ion in
        the CHIANTI database.
    index: `int`, `slice`, or array-like, optional
        Index or slice for which to extract energy level data.
    """

    def __init__(self, elvlc, index=None):
        self._elvlc = elvlc
        self._index = np.s_[:] if index is None else index

    def __len__(self):
        try:
            return self.level.shape[0]
        except IndexError:
            return 0

    def __getitem__(self, index):
        # NOTE: This throws an IndexError to stop iteration
        _ = self.level[index]
        return type(self)(self._elvlc, index=index)

    def __repr__(self):
        return f"""Level: {self.level}
Configuration: {self.configuration}
Spin: {self.spin}
Total Angular Momentum: {self.total_angular_momentum}
Orbital Angular Momentum: {self.orbital_angular_momentum_label}
Energy: {self.energy}"""

    @property
    def level(self):
        "Index of each level."
        return self._elvlc['level'][self._index]

    @property
    def configuration(self):
        "Label denoting the electronic configuration."
        return self._elvlc['config'][self._index]

    @property
    def multiplicity(self):
        "Multiplicity, :math:`2S+1`"
        return self._elvlc['multiplicity'][self._index]

    @property
    @u.quantity_input
    def spin(self) -> u.dimensionless_unscaled:
        "Spin, :math:`S`, of the electronic configuration."
        return (self.multiplicity - 1)/2

    @property
    def total_angular_momentum(self):
        "Total angular momentum number :math:`J`."
        return self._elvlc['J'][self._index]

    @property
    def weight(self):
        "Statistical weight, :math:`2J + 1`."
        return 2*self.total_angular_momentum + 1

    @property
    def orbital_angular_momentum_label(self):
        "Orbital angular momentum label."
        return self._elvlc['L_label'][self._index]

    @property
    @u.quantity_input
    def orbital_angular_momentum(self) -> u.dimensionless_unscaled:
        "Orbital angular momentum number :math:`L`."
        return np.array(
            [_ORBITAL_ANGULAR_MOMENTUM_LOOKUP[l]
            for l in self.orbital_angular_momentum_label]
        ).squeeze()

    @property
    def is_observed(self):
        "True if the energy of the level is from laboratory measurements."
        return self._elvlc['E_obs'][self._index].to_value('cm-1') != -1

    @property
    @u.quantity_input
    def energy(self) -> u.eV:
        """
        Energy of each level. Defaults to observed energy and falls back to
        theoretical energy if no measured energy is available.
        """
        energy = np.where(self.is_observed,
                          self._elvlc['E_obs'][self._index],
                          self._elvlc['E_th'][self._index])
        return energy.to('eV', equivalencies=u.equivalencies.spectral())


class Transitions:
    """
    An object for holding atomic transition data from CHIANTI.

    .. warning:: This is not meant to be instantiated directly,
                 but rather accessed through `fiasco.Ion.transitions`.

    Parameters
    ----------
    levels: `~fiasco.Levels`
        Data structure holding information about all of the energy levels
        for a given ion in the CHIANTI database.
    wgfa: `~fiasco.io.datalayer.DataIndexer`
        Pointer to the transition information for a given ion in
        the CHIANTI database.
    """

    def __init__(self, levels, wgfa):
        self._levels = levels
        self._wgfa = wgfa

    def __len__(self):
        return self.wavelength.shape[0]

    @property
    def is_twophoton(self):
        """
        True if the transition is a two-photon decay
        """
        return np.logical_and(self.wavelength == 0,
                              self.upper_level<=10)

    @property
    def is_autoionization(self):
        """
        True if the transition corresponds to an autoionization.
        """
        return np.logical_and(self.wavelength==0, ~self.is_twophoton)

    @property
    def is_bound_bound(self):
        """
        True for bound-bound transitions.
        """
        return self._wgfa['wavelength'] != 0

    @property
    def is_observed(self):
        """
        True for transitions that connect two observed energy levels
        """
        return self._wgfa['wavelength'] > 0

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
    def delta_energy(self) -> u.eV:
        "Energy spacing between the upper and lower level of the transition."
        indices = np.vstack([vectorize_where(self._levels.level, self.lower_level),
                             vectorize_where(self._levels.level, self.upper_level)])
        return np.diff(self._levels.energy[indices], axis=0).squeeze()
