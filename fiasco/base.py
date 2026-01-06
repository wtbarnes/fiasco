"""
Base classes for access to CHIANTI ion data.

These classes are not meant to be instantiated directly by the user.
"""
import astropy.units as u
import plasmapy.particles

from packaging.version import Version
from plasmapy.utils import roman

import fiasco

from fiasco.io.datalayer import DataIndexer
from fiasco.io.factory import all_subclasses
from fiasco.io.generic import GenericIonParser
from fiasco.util import check_database, parse_ion_name, periodic_table_period
from fiasco.util.exceptions import MissingIonError

__all__ = ['IonBase']


class IonBase:
    """
    Base class for accessing data attached to a particular ion.

    .. note:: This is not meant to be instantiated directly by the user
              and primarily serves as a base class for `~fiasco.Ion`.

    Parameters
    ----------
    ion_name : str
        Name of ion, e.g. for Fe V, 'Fe 5', 'iron 5', 'Fe 4+'
    hdf5_dbase_root : str, optional
        Path to the root of the CHIANTI HDF5 database.
        If not provided, the default path is used.
    kwargs :
        Further keyword arguments are passed to `fiasco.util.check_database`.
    """

    def __init__(self, ion_name, hdf5_dbase_root=None, **kwargs):
        # base rep is a tuple of integers (atomic_number, ionization_stage)
        self._base_rep = parse_ion_name(ion_name)
        if hdf5_dbase_root is None:
            self.hdf5_dbase_root = fiasco.defaults['hdf5_dbase_root']
        else:
            self.hdf5_dbase_root = hdf5_dbase_root
        check_database(self.hdf5_dbase_root, **kwargs)
        if self.ion_name not in fiasco.list_ions(self.hdf5_dbase_root, sort=False):
            raise MissingIonError(f'{self.ion_name} not found in {self.hdf5_dbase_root}')
        # Put import here to avoid circular imports
        from fiasco import log
        self.log = log
        # Warn users if the database they are using is potentially stale.
        self._check_dbase_fiasco_version()

    def _check_dbase_fiasco_version(self):
        "Warn if database was generated with an earlier version of fiasco."
        dbase_version = DataIndexer(self.hdf5_dbase_root, '/').fiasco_version
        if dbase_version is not None:
            current_version = Version(fiasco.__version__)
            if dbase_version < current_version:
                warning_msg = (f'{self.hdf5_dbase_root} was produced with an earlier version of '
                               f'fiasco ({dbase_version}) than the current version ({current_version}). '
                               'You may need to rebuild the HDF5 database. '
                               'See https://fiasco.readthedocs.io/en/stable/how_to_guides.html#re-building-the-hdf5-database '
                               'for more information.')
                self.log.warning(warning_msg)

    @property
    def atomic_number(self):
        """The atomic number of the element, :math:`Z`."""
        return plasmapy.particles.atomic_number(self._base_rep[0])

    @property
    def element_name(self):
        """The full name of the element, e.g. "hydrogen"."""
        return plasmapy.particles.element_name(self.atomic_number)

    @property
    def atomic_symbol(self):
        """The standard atomic symbol for the element, e.g. "H" for hydrogen."""
        return plasmapy.particles.atomic_symbol(self.atomic_number)

    @property
    def ion_name(self):
        """The name of the element and ionization stage, e.g. "Fe 11"."""
        return f'{self.atomic_symbol} {self.ionization_stage}'

    @property
    def ionization_stage(self):
        """Number denoting the degree of ionization, with 1 denoting the neutral stage."""
        return self._base_rep[1]

    @property
    def charge_state(self):
        "Total number of electrons removed, :math:`z`."
        return self.ionization_stage - 1

    @property
    def isoelectronic_sequence(self):
        "Atomic symbol denoting to which isoelectronic sequence this ion belongs."
        if (Z_iso := self.atomic_number - self.charge_state) > 0:
            return plasmapy.particles.atomic_symbol(Z_iso)

    @property
    @u.quantity_input
    def mass(self) -> u.g:
        "Atomic mass of the ion."
        return plasmapy.particles.Particle(self.ion_name_roman).mass

    @property
    def period(self):
        "Row of the periodic table to which the element belongs"
        return periodic_table_period(self.atomic_number)

    @property
    def _ion_name(self):
        # Old CHIANTI format, only preserved for internal data access
        return f'{self.atomic_symbol.lower()}_{self.ionization_stage}'

    @property
    def ionization_stage_roman(self):
        "Ionization stage in roman numeral format."
        return roman.to_roman(int(self.ionization_stage))

    @property
    def ion_name_roman(self):
        "Name of the element and ionization stage in roman numeral format."
        return f'{self.atomic_symbol} {self.ionization_stage_roman}'

    @property
    def _advanced_model_levels(self):
        dl = DataIndexer.create_indexer(self.hdf5_dbase_root,
                                        f'{self.atomic_symbol.lower()}/{self._ion_name}')
        if (advanced_model:=dl.get('advanced_model', None)):
            return int(advanced_model.get('n_levels', 0))
        return 0

    @property
    def _verner(self):
        data_path = '/'.join([self.atomic_symbol.lower(), self._ion_name, 'continuum',
                              'verner_short'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _hseq(self):
        data_path = '/'.join([self.atomic_symbol.lower(), 'continuum', 'hseq_2photon'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _heseq(self):
        data_path = '/'.join([self.atomic_symbol.lower(), 'continuum', 'heseq_2photon'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _abund(self):
        data_path = '/'.join([self.atomic_symbol.lower(), 'abundance'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _ion_fraction(self):
        data_path = '/'.join([self.atomic_symbol.lower(), self._ion_name, 'ioneq'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _ip(self):
        data_path = '/'.join([self.atomic_symbol.lower(), self._ion_name, 'ip'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)


def add_property(cls, filetype):
    """
    Dynamically add filetype properties to base data access class
    """
    def property_template(self):
        data_path = '/'.join([self.atomic_symbol.lower(), self._ion_name, filetype])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    property_template.__doc__ = f'Data in {filetype} type file'
    property_template.__name__ = f'_{"_".join(filetype.split("/"))}'
    setattr(cls, property_template.__name__, property(property_template))


# Collect the filetypes and add the methods
all_ext = [cls.filetype for cls in all_subclasses(GenericIonParser) if hasattr(cls, 'filetype')]
for filetype in all_ext:
    add_property(IonBase, filetype)
    add_property(IonBase, '/'.join(['dielectronic', filetype]))
