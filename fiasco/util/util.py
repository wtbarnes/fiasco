"""
Basic utilities
"""
import configparser
import os
import pathlib
import plasmapy.particles
import sys

from packaging.version import Version
from plasmapy.utils import roman

FIASCO_HOME = pathlib.Path.home() / '.fiasco'
FIASCO_RC = FIASCO_HOME / 'fiascorc'

__all__ = ['setup_paths', 'get_chianti_catalog', 'read_chianti_version', 'parse_ion_name']


def parse_ion_name(ion_name):
    """
    Parse the atomic number and ionization stage from representation of ion.

    This function can take a number of formats for the ion name. As an example, all
    of the following representations of Fe 18, that is, an iron ion
    with 17 electrons removed and a total charge of +17, will return (26,18):

    1. ``'Fe 18'`, ``'fe 18'`` (atomic symbol and ionization stage)
    2. ``'Fe 17+'`` (atomic symbol and charge state)
    3. ``'Iron 18'``, ``'iron 18'`` (element name and ionization stage)
    4. ``'Fe XVIII'``, ``'fe xviii'`` (atomic symbol and ionization stage in spectroscopic notation)
    5. ``'26 18'`` (atomic number and ionization stage)
    6. ``(26, 18)`` (tuple of any combination of the above)
    """
    if isinstance(ion_name, tuple):
        element, ion = ion_name
    elif isinstance(ion_name, str):
        element, ion = ion_name.split('_' if '_' in ion_name else None)
    else:
        raise TypeError(f'Unrecognized type {type(ion_name)}. ion_name must be either a string or tuple')
    # Parse element string
    if isinstance(element, str):
        element = element.capitalize()
    element = plasmapy.particles.atomic_number(element)
    # Parse ion string
    if isinstance(ion, str):
        if '+' in ion:
            ion = f"{int(ion.strip('+')) + 1}"
        if roman.is_roman_numeral(ion.upper()):
            ion = roman.from_roman(ion.upper())
    ion = int(ion)
    return (element, ion)


def setup_paths():
    """
    Parse .rc file and set ASCII and HDF5 database paths.
    """
    paths = {}
    if FIASCO_RC.is_file():
        config = configparser.ConfigParser()
        config.read(FIASCO_RC)
        if 'database' in config:
            paths = dict(config['database'])

    if 'ascii_dbase_root' not in paths:
        paths['ascii_dbase_root'] = FIASCO_HOME / 'chianti_dbase'
    if 'hdf5_dbase_root' not in paths:
        paths['hdf5_dbase_root'] = FIASCO_HOME / 'chianti_dbase.h5'
    paths['ascii_dbase_root'] = pathlib.Path(paths['ascii_dbase_root'])
    paths['hdf5_dbase_root'] = pathlib.Path(paths['hdf5_dbase_root'])

    return paths


def get_chianti_catalog(ascii_dbase_root):
    """
    Return a dictionary of all CHIANTI data files, separated by category.

    Parse CHIANTI filetree and return list of all files, separated by category. This will be only
    be useful when dealing with the raw ASCII data.

    Parameters
    ----------
    ascii_dbase_root: path-like
        Path to the top of the CHIANTI filetree

    Returns
    -------
    : `dict`
        All CHIANTI files, separated by category. The resulting dictionary should have the
        following keys: 'abundance_files', 'ioneq_files', 'ip_files', 'continuum_files',
        'ion_files'.
    """
    ascii_dbase_root = pathlib.Path(ascii_dbase_root)
    # TODO: Replace usage with pathlib, noting that pathlib does not
    # have a direct equivalent to os.walk

    skip_dirs = ['version_3', 'deprecated', 'masterlist', 'ioneq', 'dem', 'ancillary_data', 'ip',
                 'abundance', 'continuum', 'instrument_responses']
    # List of all files associated with ions
    ion_files = []
    for root, sub, files in os.walk(ascii_dbase_root):
        if all(sd not in root for sd in skip_dirs) and all(sd not in sub for sd in skip_dirs):
            ion_files += [f for f in files if f[0] != '.']

    # List all of the non-ion files, excluding any "dot"/hidden files
    def walk_sub_dir(subdir):
        subdir_files = []
        subdir_root = ascii_dbase_root / subdir
        for root, _, files in os.walk(subdir_root):
            subdir_files += [os.path.relpath(os.path.join(root, f), subdir_root) for f in files
                             if f[0] != '.']

        return subdir_files

    non_ion_subdirs = ['abundance', 'ip', 'ioneq', 'continuum', 'dem']
    all_files = {f'{sd}_files': walk_sub_dir(sd) for sd in non_ion_subdirs}
    all_files['ion_files'] = ion_files

    return all_files


def read_chianti_version(ascii_dbase_root):
    """
    Read the CHIANTI version number from the ASCII database.

    Parameters
    ----------
    asciii_dbase_root: `str` or `pathlib.Path`

    Returns
    -------
    : `packaging.version.Version`
    """
    version_file = pathlib.Path(ascii_dbase_root) / 'VERSION'
    with version_file.open() as f:
        lines = f.readlines()
    version = lines[0].strip()
    return Version(version)


def query_yes_no(question, default="yes"):
    """
    Ask a yes/no question via raw_input() and return their answer.
    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).
    The "answer" return value is one of "yes" or "no".

    See `this gist <https://gist.github.com/hrouault/1358474>`_
    """
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        msg = f"invalid default answer: {default}"
        raise ValueError(msg)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")
