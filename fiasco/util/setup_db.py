"""
Various functions for downloading and setting up the database
"""
import h5py
import hashlib
import numpy as np
import pathlib
import tarfile

from astropy.config import set_temp_cache
from astropy.utils.console import ProgressBar
from astropy.utils.data import download_file

import fiasco.io

from fiasco.io import DataIndexer
from fiasco.util.exceptions import MissingASCIIFileError, UnsupportedVersionError
from fiasco.util.util import get_chianti_catalog, query_yes_no

FIASCO_HOME = pathlib.Path.home() / '.fiasco'
CHIANTI_URL = 'http://download.chiantidatabase.org/CHIANTI_v{version}_database.tar.gz'
# List in order (oldest to newest) the supported versions of the database
SUPPORTED_VERSIONS = [
    '8.0',
    '8.0.2',
    '8.0.6',
    '8.0.7',
]
LATEST_VERSION = SUPPORTED_VERSIONS[-1]

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
    hdf5_dbase_root : `str` or `~pathlib.Path`
    ascii_dbase_root : `str` or `~pathlib.Path`, optional
    ask_before : bool, optional
    ascii_dbase_url : str, optional

    See Also
    --------
    download_dbase
    build_hdf5_dbase
    """
    ascii_dbase_root = kwargs.get('ascii_dbase_root',
                                  fiasco.defaults['ascii_dbase_root'])

    hdf5_dbase_root = pathlib.Path(hdf5_dbase_root)
    ascii_dbase_root = pathlib.Path(ascii_dbase_root)

    # Useful for building, downloading non-interactively
    ask_before = kwargs.get('ask_before', True)
    if hdf5_dbase_root.is_file():
        return None
    if ask_before:
        question = f"No HDF5 database found at {hdf5_dbase_root}. Build it now?"
        answer = query_yes_no(question, default='yes')
        if not answer:
            # FIXME: how to gracefully handle this exit?
            return None
    # VERSION file as a proxy for whether whole dbase exists
    # TODO: version checking, download newest version if installed version is out of date
    if not (ascii_dbase_root / "VERSION").is_file():
        ascii_dbase_url = kwargs.get('ascii_dbase_url')
        if ascii_dbase_url is None:
            ascii_dbase_url = CHIANTI_URL.format(version=LATEST_VERSION)
        if ask_before:
            question = f"No CHIANTI database found at {ascii_dbase_root}. Download it from {ascii_dbase_url}?"
            answer = query_yes_no(question, default='no')
            if not answer:
                # FIXME: how to gracefully handle this exit?
                return None
        download_dbase(ascii_dbase_url, ascii_dbase_root)
    # If we made it this far, build the HDF5 database
    files = kwargs.get('files')
    build_hdf5_dbase(ascii_dbase_root, hdf5_dbase_root, files=files)
    # Ensure that an unsupported version is not being passed to fiasco
    # NOTE: this check is only meant to be bypassed when testing new
    # versions. Hence, this kwarg is not documented
    if kwargs.get('check_chianti_version', True):
        check_database_version(hdf5_dbase_root)


def check_database_version(hdf5_dbase_root):
    dl = DataIndexer.create_indexer(hdf5_dbase_root, '/')
    if dl.version not in SUPPORTED_VERSIONS:
        raise UnsupportedVersionError(
            f'CHIANTI {dl.version} is not in the list of supported versions {SUPPORTED_VERSIONS}.')


def download_dbase(ascii_dbase_url, ascii_dbase_root):
    """
    Download the CHIANTI database in ASCII format
    """
    from fiasco import log
    log.debug(f'Downloading database from {ascii_dbase_url}')
    log.debug(f'Downloading database to {ascii_dbase_root}')
    tar_tmp_dir = FIASCO_HOME / 'tmp'
    tar_tmp_dir.mkdir(exist_ok=True, parents=True)
    with set_temp_cache(path=tar_tmp_dir, delete=True):
        tmp_tar = download_file(ascii_dbase_url, cache=True, show_progress=True)
        with tarfile.open(tmp_tar) as tar:
            tar.extractall(path=ascii_dbase_root)


def md5hash(path):
    path = pathlib.Path(path)
    with path.open('rb') as f:
        return hashlib.md5(f.read()).hexdigest()


def build_hdf5_dbase(ascii_dbase_root, hdf5_dbase_root, files=None):
    """
    Assemble HDF5 file from raw ASCII CHIANTI database.

    Parameters
    ----------
    ascii_dbase_root : `str` or `~pathlib.Path`
        Path to top of CHIANTI database tree
    hdf5_dbase_root : `str` or `~pathlib.Path`
        Path to HDF5 file
    files : `list` or `dict`, optional
        A list of files to update in the HDF5 database. By default,
        this is all of the files in `ascii_dbase_root`. If a `dict`, the
        dictionary keys must contain filenames and the items corresponding
        expected md5 hash of the file. Building the database will fail if any
        of the md5 hashes is not as expected.
    """
    # Import the logger here to avoid circular imports
    from fiasco import log

    if files is None:
        log.debug('Adding all files to CHIANTI database')
        files = []
        tmp = get_chianti_catalog(ascii_dbase_root)
        for k in tmp:
            files += tmp[k]
    log.debug(f'Building HDF5 database in {hdf5_dbase_root}')
    with ProgressBar(len(files)) as progress:
        with h5py.File(hdf5_dbase_root, 'a') as hf:
            for f in files:
                parser = fiasco.io.Parser(f, ascii_dbase_root=ascii_dbase_root)
                if isinstance(files, dict):
                    expected = files[f]
                    actual = md5hash(parser.full_path)
                    if expected != actual:
                        raise RuntimeError(f'Hash of {parser.full_path} ({actual}) did not match expected hash ({expected})')
                try:
                    df = parser.parse()
                except MissingASCIIFileError as e:
                    log.debug(f'{e}. Not including {f} in {hdf5_dbase_root}')
                else:
                    parser.to_hdf5(hf, df)
                progress.update()
            # Build an index for quick lookup of all ions in database
            from fiasco import list_ions  # import here to avoid circular imports

            # Delete it if it already exists to ensure the index is rebuilt
            if 'ion_index' in hf:
                del hf['ion_index']
            ion_list = list_ions(hdf5_dbase_root)
            ds = hf.create_dataset('ion_index', data=np.array(ion_list).astype(np.string_))
            ds.attrs['unit'] = 'SKIP'
