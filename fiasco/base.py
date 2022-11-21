"""
Base classes for access to CHIANTI ion data
"""
import plasmapy.particles

from plasmapy.utils import roman

import fiasco

from fiasco.io.datalayer import DataIndexer
from fiasco.io.factory import all_subclasses
from fiasco.io.generic import GenericIonParser
from fiasco.util import check_database, parse_ion_name
from fiasco.util.exceptions import MissingIonError

__all__ = ['IonBase', 'ContinuumBase']


class Base:
    """
    Base class for setting up ion metadata and building database if necessary.

    Parameters
    ----------
    ion_name : str
    hdf5_dbase_root : str, optional
    kwargs :
        kwargs are passed to `check_database`.
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

    @property
    def atomic_number(self):
        return plasmapy.particles.atomic_number(self._base_rep[0])

    @property
    def element_name(self):
        return plasmapy.particles.element_name(self.atomic_number)

    @property
    def atomic_symbol(self):
        return plasmapy.particles.atomic_symbol(self.atomic_number)

    @property
    def ion_name(self):
        return f'{self.atomic_symbol} {self.ionization_stage}'

    @property
    def ionization_stage(self):
        return self._base_rep[1]

    @property
    def charge_state(self):
        return self.ionization_stage - 1

    @property
    def _ion_name(self):
        # Old CHIANTI format, only preserved for internal data access
        return f'{self.atomic_symbol.lower()}_{self.ionization_stage}'

    @property
    def ionization_stage_roman(self):
        return roman.to_roman(int(self.ionization_stage))

    @property
    def ion_name_roman(self):
        return f'{self.atomic_symbol} {self.ionization_stage_roman}'


class ContinuumBase(Base):
    """
    Base class for retrieving continuum datasets.

    .. note:: This is not meant to be instantiated directly by the user
              and primarily serves as a base class for `~fiasco.Ion`.
    """

    @property
    def _gffgu(self):
        data_path = '/'.join(['continuum', 'gffgu'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _gffint(self):
        data_path = '/'.join(['continuum', 'gffint'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _klgfb(self):
        data_path = '/'.join(['continuum', 'klgfb'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _verner(self):
        data_path = '/'.join([self.atomic_symbol.lower(), self._ion_name, 'continuum',
                              'verner_short'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _itoh(self):
        data_path = '/'.join([self.atomic_symbol.lower(), 'continuum', 'itoh'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _hseq(self):
        data_path = '/'.join([self.atomic_symbol.lower(), 'continuum', 'hseq_2photon'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)

    @property
    def _heseq(self):
        data_path = '/'.join([self.atomic_symbol.lower(), 'continuum', 'heseq_2photon'])
        return DataIndexer.create_indexer(self.hdf5_dbase_root, data_path)


class IonBase(Base):
    """
    Base class for accessing CHIANTI data attached to a particular ion

    .. note:: This is not meant to be instantiated directly by the user
              and primarily serves as a base class for `~fiasco.Ion`.

    Parameters
    ----------
    ion_name : `str`
        Name of ion, e.g. for Fe V, 'Fe 5', 'iron 5', 'Fe 4+'
    hdf5_path : `str`, optional

    Examples
    --------
    """

    @property
    def _abundance(self):
        data_path = '/'.join([self.atomic_symbol.lower(), 'abundance'])
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
add_property(IonBase, 'ip')
add_property(IonBase, 'ioneq')
