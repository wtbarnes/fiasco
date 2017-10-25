"""
Various functions for downloading and setting up the database
"""

import os
import warnings
import tarfile
import configparser
    
import h5py
from astropy.config import set_temp_cache
from astropy.utils.data import download_file
from astropy.utils.console import ProgressBar

from .yes_no import query_yes_no
import fiasco.io

__all__ = ['setup_paths', 'download_dbase', 'get_masterlist', 'build_hdf5_dbase']

FIASCO_HOME = os.path.join(os.environ['HOME'], '.fiasco')
CHIANTI_URL = 'http://www.chiantidatabase.org/download/CHIANTI_{version}_data.tar.gz'
LATEST_VERSION = '8.0.6'


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


def download_dbase(ascii_dbase_root, version=None, ask_before=True):
    """
    Download the database if it does not exist locally.
    """
    # VERSION file as a proxy for whether whole dbase exists
    # TODO: version checking, download newest version if installed version is out of date
    if os.path.isfile(os.path.join(ascii_dbase_root, 'VERSION')):
        return None
    if version is None:
        dbase_url = CHIANTI_URL.format(version=LATEST_VERSION)
    else:
        dbase_url = CHIANTI_URL.format(version=version)
    if ask_before:
        question = "No CHIANTI database found at {}. Download it from {}?"
        answer = query_yes_no(question.format(ascii_dbase_root, dbase_url), default='no')
        if not answer:
            return None
    
    # Download and extract
    tar_tmp_dir = os.path.join(FIASCO_HOME, 'tmp')
    if not os.path.exists(tar_tmp_dir):
        os.makedirs(tar_tmp_dir)
    with set_temp_cache(path=tar_tmp_dir, delete=True):
        tmp_tar = download_file(dbase_url, cache=True, show_progress=True)
        with tarfile.open(tmp_tar) as tar:
            tar.extractall(path=ascii_dbase_root)


def build_hdf5_dbase(ascii_dbase_root, hdf5_dbase_root, ask_before=True):
    """
    Assemble HDF5 file from raw ASCII CHIANTI database
    """
    # Check and ask
    if os.path.isfile(hdf5_dbase_root):
        return None
    if not os.path.isfile(os.path.join(ascii_dbase_root, 'VERSION')):
        raise FileNotFoundError('No CHIANTI database found at {}'.format(ascii_dbase_root))
    if ask_before:
        question = """No HDF5 database found at {}. Build it now?"""
        answer = query_yes_no(question.format(hdf5_dbase_root), default='yes')
        if not answer:
            return None
    
    # Build database
    all_files = []
    tmp = get_masterlist(ascii_dbase_root)
    for k in tmp:
        all_files += tmp[k]
    with ProgressBar(len(all_files)) as progress:
        with h5py.File(hdf5_dbase_root, 'a') as hf:
            for af in all_files:
                parser = fiasco.io.Parser(af)
                df = parser.parse()
                if df is None:
                    warnings.warn('Not including {} in {}'.format(af, hdf5_dbase_root), stacklevel=2)
                else:
                    parser.to_hdf5(hf, df)
                progress.update()


def get_masterlist(ascii_dbase_root):
    """
    Parse CHIANTI filetree and return several useful lists for indexing the database.

    Note
    -----
    This will be only be useful when dealing with the raw ASCII data.
    """
    skip_dirs = ['version_3', 'deprecated', 'masterlist', 'ioneq', 'dem', 'ancillary_data', 'ip', 'abundance',
                 'continuum', 'instrument_responses']
    # List of all files associated with ions
    ion_files = []
    for root, sub, files in os.walk(ascii_dbase_root):
        if not any([sd in root for sd in skip_dirs]) and not any([sd in sub for sd in skip_dirs]):
            ion_files += files

    # List all of the non-ion files
    def walk_sub_dir(subdir):
        subdir_files = []
        subdir_root = os.path.join(ascii_dbase_root, subdir)
        for root, _, files in os.walk(subdir_root):
            subdir_files += [os.path.relpath(os.path.join(root, f), subdir_root) for f in files]
        
        return subdir_files

    non_ion_subdirs = ['abundance', 'ioneq', 'ip']
    all_files = {'{}_files'.format(sd): walk_sub_dir(sd) for sd in non_ion_subdirs}
    all_files['ion_files'] = ion_files

    return all_files
    