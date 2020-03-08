"""
Basic utilities
"""
import os
import sys
import configparser
from builtins import input

FIASCO_HOME = os.path.join(os.environ['HOME'], '.fiasco')

__all__ = ['setup_paths', 'get_masterlist']


def setup_paths():
    """
    Parse .rc file and set ASCII and HDF5 database paths.
    """
    paths = {}
    if os.path.isfile(os.path.join(FIASCO_HOME, 'fiascorc')):
        config = configparser.ConfigParser()
        config.read(os.path.join(FIASCO_HOME, 'fiascorc'))
        if 'database' in config:
            paths = dict(config['database'])

    if 'ascii_dbase_root' not in paths:
        paths['ascii_dbase_root'] = os.path.join(FIASCO_HOME, 'chianti_dbase')
    if 'hdf5_dbase_root' not in paths:
        paths['hdf5_dbase_root'] = os.path.join(FIASCO_HOME, 'chianti_dbase.h5')

    return paths


def get_masterlist(ascii_dbase_root):
    """
    Parse CHIANTI filetree and return list of all files, separated by category. This will be only
    be useful when dealing with the raw ASCII data.
    """
    skip_dirs = ['version_3', 'deprecated', 'masterlist', 'ioneq', 'dem', 'ancillary_data', 'ip',
                 'abundance', 'continuum', 'instrument_responses']
    # List of all files associated with ions
    ion_files = []
    for root, sub, files in os.walk(ascii_dbase_root):
        if not any([sd in root for sd in skip_dirs]) and not any([sd in sub for sd in skip_dirs]):
            ion_files += [f for f in files if f[0] != '.']

    # List all of the non-ion files, excluding any "dot"/hidden files
    def walk_sub_dir(subdir):
        subdir_files = []
        subdir_root = os.path.join(ascii_dbase_root, subdir)
        for root, _, files in os.walk(subdir_root):
            subdir_files += [os.path.relpath(os.path.join(root, f), subdir_root) for f in files
                             if f[0] != '.']

        return subdir_files

    non_ion_subdirs = ['abundance', 'ioneq', 'ip', 'continuum']
    all_files = {f'{sd}_files': walk_sub_dir(sd) for sd in non_ion_subdirs}
    all_files['ion_files'] = ion_files

    return all_files


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
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")
