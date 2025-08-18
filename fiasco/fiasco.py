"""
Package-level functions.
"""
import astropy.table
import astropy.units as u
import mendeleev
import numpy as np
import pathlib
import plasmapy.particles

from astropy.utils.data import get_pkg_data_path
from plasmapy.particles.exceptions import InvalidParticleError
from scipy.interpolate import interp1d

import fiasco

from fiasco.io import DataIndexer
from fiasco.util import parse_ion_name

__all__ = [
    'list_elements',
    'list_ions',
    'proton_electron_ratio',
    'get_isoelectronic_sequence',
    'dielectronic_recombination_suppression',
]


def list_elements(hdf5_dbase_root=None, sort=True):
    """
    List all available elements in the CHIANTI database.

    Parameters
    ----------
    hdf5_dbase_root: path-like, optional
        If not specified, will default to that specified in ``fiasco.defaults``.
    sort: `bool`, optional
        If True, sort the list of elements by increasing atomic number.
    """
    if hdf5_dbase_root is None:
        hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
    elements = []
    root = DataIndexer.create_indexer(hdf5_dbase_root, '/')
    for f in root.fields:
        try:
            elements.append(plasmapy.particles.atomic_symbol(f.capitalize()))
        except InvalidParticleError:
            continue
    if sort:
        elements = sorted(elements, key=lambda x: plasmapy.particles.atomic_number(x))
    return elements


def list_ions(hdf5_dbase_root=None, sort=True):
    """
    List all available ions in the CHIANTI database

    Parameters
    ----------
    hdf5_dbase_root: path-like, optional
        If not specified, will default to that specified in ``fiasco.defaults``.
    sort: `bool`, optional
        If True, sort the list of elements by increasing atomic number.
    """
    if hdf5_dbase_root is None:
        hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
    root = DataIndexer(hdf5_dbase_root, '/')
    # NOTE: get the list from the index if possible. This is ~30x faster
    try:
        ions = root['ion_index']
    except KeyError:
        ions = []
        for f in root.fields:
            try:
                el = plasmapy.particles.atomic_symbol(f.capitalize())
                for i in root[f].fields:
                    if f == i.split('_')[0]:
                        ions.append(f"{el} {i.split('_')[1]}")
            except InvalidParticleError:
                continue
    # Optional because adds significant overhead
    if sort:
        ions = sorted(ions, key=lambda x: (plasmapy.particles.atomic_number(x.split()[0]),
                                           int(x.split()[1])))
    # NOTE: when grabbing straight from the index and not sorting, the result will be
    # a numpy array. Cast to a list to make sure the return type is consistent for
    # all possible inputs
    return ions.tolist() if isinstance(ions, np.ndarray) else ions


def get_isoelectronic_sequence(element, hdf5_dbase_root=None):
    """
    List of ions in the isoelectronic sequence of ``element``.

    Ions in the same isoelectronic sequence, that is, ions that have the
    same number of bound electrons, :math:`Z - z`, often share common properties
    despite having different charge states and being from different elements.
    These so-called isoelectronic sequences are typically denoted by the element
    with an atomic number equal to the number of bound electrons, e.g. C II is in
    the boron isoelectronic sequence or equivalently may be said to be boron-like.
    Given the name of a sequence, as denoted by an element label, this function
    returns a list of all ions in that sequence.

    Parameters
    ----------
    element: `str`, `int`
        Name of sequence. Can be either the full name (e.g. "hydrogren"),
        the atomic symbol (e.g. "H") or the atomic number (e.g. 1)
    hdf5_dbase_root: path-like, optional
        If not specified, will default to that specified in ``fiasco.defaults``.
    """
    Z_iso = plasmapy.particles.atomic_number(element)
    all_ions = list_ions(hdf5_dbase_root=hdf5_dbase_root)

    def _is_in_sequence(ion):
        Z, z = parse_ion_name(ion)
        return Z_iso == (Z - z + 1)

    return [ion for ion in all_ions if _is_in_sequence(ion)]


@u.quantity_input
def proton_electron_ratio(temperature: u.K, **kwargs):
    """
    Calculate ratio between proton and electron densities as a function of temperature
    according to Eq. 7 of :cite:t:`young_chianti-atomic_2003`.

    Parameters
    ----------
    temperature : `~astropy.units.Quantity`

    See Also
    --------
    fiasco.Ion : Accepts same keyword arguments for setting database and dataset names
    """
    # Import here to avoid circular imports
    from fiasco import log
    h_2 = fiasco.Ion('H +1', temperature, **kwargs)
    numerator = h_2.abundance * h_2._ion_fraction[h_2._instance_kwargs['ionization_fraction']]['ionization_fraction']
    denominator = u.Quantity(np.zeros(numerator.shape))
    for el_name in list_elements(h_2.hdf5_dbase_root):
        el = fiasco.Element(el_name, temperature, **h_2._instance_kwargs)
        try:
            abundance = el.abundance
        except KeyError:
            abund_file = el[0]._instance_kwargs['abundance']
            log.warning(
                f'Not including {el.atomic_symbol}. Abundance not available from {abund_file}.')
            continue
        for ion in el:
            ionization_file = ion._instance_kwargs['ionization_fraction']
            # NOTE: We use ._ion_fraction here rather than .ionization_fraction to avoid
            # doing an interpolation to the temperature array every single time and instead only
            # interpolate once at the end.
            # It is assumed that the ionization_fraction temperature array for each ion is the same.
            try:
                ionization_fraction = ion._ion_fraction[ionization_file]['ionization_fraction']
                t_ionization_fraction = ion._ion_fraction[ionization_file]['temperature']
            except KeyError:
                log.warning(
                    f'Not including {ion.ion_name}. Ionization fraction not available from {ionization_file}.')
                continue
            denominator += ionization_fraction * abundance * ion.charge_state

    ratio = numerator / denominator
    f_interp = interp1d(t_ionization_fraction.to(temperature.unit).value,
                        ratio.value,
                        kind='linear',
                        bounds_error=False,
                        fill_value=(ratio[0], ratio[-1]))

    return u.Quantity(f_interp(temperature.value))


@u.quantity_input
def dielectronic_recombination_suppression(ion, density:u.Unit('cm-3')):
    """
    Density-dependent suppression factor for dielectronic recombination.

    Calculates the density-dependent suppression factor for dielectronic
    recombination following the formulation of :cite:t:`nikolic_suppression_2018`.

    Parameters
    ----------
    ion
    density: `~astropy.units.Quantity`
    """
    if ion.isoelectronic_sequence is None:
        return 1
    # "A" factor
    A_N = _nikolic_a_factor(ion)
    # Activation log density (Eq. 3 of Nikolic et al. 2018)
    x_a0 = 10.1821
    q_0 = (1 - np.sqrt(2/3/ion.charge_state))*A_N/np.sqrt(ion.charge_state)
    T_0 = 5e4*u.K * q_0**2
    x_a = x_a0 + np.log((ion.charge_state/q_0)**7*np.sqrt(ion.temperature/T_0))
    # Suppression factor (Eq. 2 of Nikolic et al. 2018)
    width = 5.64586
    x = np.log10(density.to_value('cm-3'))
    suppression = np.exp(-((x-x_a)/width*np.sqrt(np.log(2)))**2)
    suppression = np.where(x<=x_a, 1, suppression)
    # Low-temperature correction (Eq. 14 of Nikolic et al. 2018)
    filename = pathlib.Path(get_pkg_data_path('data', package='fiasco')) / 'nikolic_table_5.dat'
    coefficient_table = astropy.table.QTable.read(filename, format='ascii.mrt')
    if ion.isoelectronic_sequence not in coefficient_table['Sequence']:
        return suppression
    row = coefficient_table[coefficient_table['Sequence']==ion.isoelectronic_sequence]
    eps_energies = u.Quantity([row[f'p_{i}']*(ion.charge_state/10)**i for i in range(6)]).sum()
    exp_factor = np.exp(-eps_energies/10/ion.thermal_energy)
    return 1 - (1 - suppression)*exp_factor


def _nikolic_a_factor(ion):
    """
    Compute :math:`A(N)` according to Equations 6 and 9 of :cite:t:`nikolic_suppression_2018`.
    """
    Z_iso = plasmapy.particles.atomic_number(ion.isoelectronic_sequence)
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
        period_iso = mendeleev.element(ion.isoelectronic_sequence).period
        N_1, N_2 = {
            2: (3,10), 3: (11,18), 4: (19,36), 5: (37,54), 6: (55,86), 7: (87,118)
        }[period_iso]
        A_N = 12 + 10*N_1 + (10*N_1 - 2*N_2)/(N_1 - N_2)*(Z_iso - N_1)
    # Compute additional modifications according to Eqs. 9, 10, and 11
    filename = pathlib.Path(get_pkg_data_path('data', package='fiasco')) / 'nikolic_table_2.dat'
    coefficient_table = astropy.table.QTable.read(filename, format='ascii.mrt')
    if Z_iso not in coefficient_table['N']:
        return A_N
    # Calculate pis/gammas. Relabel as c_i as the formula is the same
    c_i = []
    for i in range(1,7):
        row = coefficient_table[np.logical_and(coefficient_table['N'] == Z_iso, coefficient_table['i'] == i)]
        c_i.append(
            row['c_1'] + row['c_2']*ion.charge_state**row['c_3']*np.exp(-ion.charge_state/row['c_4'])
        )
    c_i = np.array(c_i)
    # Calculate psi term According to Eqs. 10 and 11
    logT = np.log10(ion.temperature.to_value('K'))
    psi = 1 + c_i[2]*np.exp(-((logT-c_i[0])/np.sqrt(2)/c_i[1])**2) + c_i[5]*np.exp(-((logT-c_i[3])/np.sqrt(2)/c_i[4])**2)
    if Z_iso < 5:
        psi = 2*psi/(1 + np.exp(-2.5e4*u.K*ion.charge_state**2/ion.temperature))
    return A_N*psi
