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

import fiasco.io
from .util import get_masterlist, query_yes_no

FIASCO_HOME = os.path.join(os.environ['HOME'], '.fiasco')
CHIANTI_URL = 'http://www.chiantidatabase.org/download/CHIANTI_{version}_data.tar.gz'
LATEST_VERSION = '8.0.6'

__all__ = ['download_dbase', 'build_hdf5_dbase']


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
        question = f"No CHIANTI database found at {ascii_dbase_root}. Download it from {dbase_url}?"
        answer = query_yes_no(question, default='no')
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
        raise FileNotFoundError(f'No CHIANTI database found at {ascii_dbase_root}')
    if ask_before:
        question = f"No HDF5 database found at {hdf5_dbase_root}. Build it now?"
        answer = query_yes_no(question, default='yes')
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
                parser = fiasco.io.Parser(af, ascii_dbase_root=ascii_dbase_root)
                df = parser.parse()
                if df is None:
                    warnings.warn(f'Not including {af} in {hdf5_dbase_root}')
                else:
                    parser.to_hdf5(hf, df)
                progress.update()
