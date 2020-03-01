"""
Various functions for downloading and setting up the database
"""

import os
import warnings
import tarfile

import h5py
from astropy.config import set_temp_cache
from astropy.utils.data import download_file
from astropy.utils.console import ProgressBar

import fiasco.io
from .util import get_masterlist, query_yes_no
from .exceptions import MissingASCIIFileError

FIASCO_HOME = os.path.join(os.environ['HOME'], '.fiasco')
CHIANTI_URL = 'http://www.chiantidatabase.org/download/CHIANTI_{version}_data.tar.gz'
LATEST_VERSION = '8.0.7'

__all__ = ['check_database', 'download_dbase', 'build_hdf5_dbase']


def check_database(hdf5_dbase_root, **kwargs):
    """
    Check if the HDF5 database exists, download the ASCII files and build the
    HDF5 database.

    Check if the HDF5 database exists. If it does not, check if the ASCII database
    exists. If it does not, download the raw ASCII database from the CHIANTI webpage.
    Finally, build the HDF5 database from the ASCII database.

    Parameters
    ----------
    hdf5_dbase_root : `str`

    See Also
    --------
    download_dbase
    build_hdf5_dbase
    """
    ascii_dbase_root = kwargs.get('ascii_dbase_root', fiasco.defaults['ascii_dbase_root'])
    # Useful for building, downloading non-interactively
    ask_before = kwargs.get('ask_before', True)
    if os.path.isfile(hdf5_dbase_root):
        return None
    if ask_before:
        question = f"No HDF5 database found at {hdf5_dbase_root}. Build it now?"
        answer = query_yes_no(question, default='yes')
        if not answer:
            # FIXME: how to gracefully handle this exit?
            return None
    # VERSION file as a proxy for whether whole dbase exists
    # TODO: version checking, download newest version if installed version is out of date
    if not os.path.isfile(os.path.join(ascii_dbase_root, 'VERSION')):
        ascii_dbase_url = kwargs.get('ascii_dbase_url', CHIANTI_URL.format(version=LATEST_VERSION))
        if ask_before:
            question = f"No CHIANTI database found at {ascii_dbase_root}. Download it from {ascii_dbase_url}?"
            answer = query_yes_no(question, default='no')
            if not answer:
                # FIXME: how to gracefully handle this exit?
                return None
        download_dbase(ascii_dbase_url, ascii_dbase_root)
    # If we made it this far, build the HDF5 database
    build_hdf5_dbase(ascii_dbase_root, hdf5_dbase_root)


def download_dbase(ascii_dbase_url, ascii_dbase_root):
    """
    Download the CHIANTI database in ASCII format
    """
    tar_tmp_dir = os.path.join(FIASCO_HOME, 'tmp')
    if not os.path.exists(tar_tmp_dir):
        os.makedirs(tar_tmp_dir)
    with set_temp_cache(path=tar_tmp_dir, delete=True):
        tmp_tar = download_file(ascii_dbase_url, cache=True, show_progress=True)
        with tarfile.open(tmp_tar) as tar:
            tar.extractall(path=ascii_dbase_root)


def build_hdf5_dbase(ascii_dbase_root, hdf5_dbase_root, files=None):
    """
    Assemble HDF5 file from raw ASCII CHIANTI database.

    Parameters
    ----------
    ascii_dbase_root : `str`
        Path to top of CHIANTI database tree
    hdf5_dbase_root : `str`
        Path to HDF5 file
    files : `list`, optional
        A list of files to update in the HDF5 database. By default,
        this is all of the files in `ascii_dbase_root`
    """
    if files is None:
        files = []
        tmp = get_masterlist(ascii_dbase_root)
        for k in tmp:
            files += tmp[k]
    with ProgressBar(len(files)) as progress:
        with h5py.File(hdf5_dbase_root, 'a') as hf:
            for f in files:
                parser = fiasco.io.Parser(f, ascii_dbase_root=ascii_dbase_root)
                try:
                    df = parser.parse()
                except MissingASCIIFileError as e:
                    # FIXME: use the logger here
                    warnings.warn(f'{e}. Not including {f} in {hdf5_dbase_root}')
                else:
                    parser.to_hdf5(hf, df)
                progress.update()
