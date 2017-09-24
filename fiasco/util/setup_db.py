"""
Various functions for downloading and setting up the database
"""

import os
import tarfile
try:
    # Python 3
    import configparser
    from urllib.request import urlretrieve
except ImportError:
    # Python 2
    import ConfigParser as configparser
    from urllib import urlretrieve
from astropy.config import set_temp_cache
from astropy.utils.data import download_file

from .yes_no import query_yes_no

__all__ = ['setup_paths', 'download_dbase']

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


def download_dbase(ascii_dbase_path, version=None):
    """
    Download the database if it does not exist locally.
    """
    # VERSION file as a proxy for whether whole dbase exists
    # TODO: version checking, download newest version if installed version is out of date
    if os.path.isfile(os.path.join(ascii_dbase_path, 'VERSION')):
        return None
    if version is None:
        dbase_url = CHIANTI_URL.format(version=LATEST_VERSION)
    else:
        dbase_url = CHIANTI_URL.format(version=version)
    # TODO: need a way to override this in "headless" situations, e.g. Travis CI
    question = "No CHIANTI database found at {}. Download it from the internet?"
    answer = query_yes_no(question.format(ascii_dbase_path), default='no')
    if not answer:
        return None
    tar_tmp_dir = os.path.join(FIASCO_HOME, 'tmp')
    if not os.path.exists(tar_tmp_dir):
        os.makedirs(tar_tmp_dir)
    with set_temp_cache(path=tar_tmp_dir, delete=True):
        tmp_tar = download_file(dbase_url, cache=True, show_progress=True)
        with tarfile.open(tmp_tar) as tar:
            tar.extractall(path=ascii_dbase_path)


def build_hdf5_dbase():
    """
    Assemble HDF5 file from raw ASCII CHIANTI database
    """
    pass