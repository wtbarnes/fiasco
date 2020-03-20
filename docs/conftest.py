# NOTE: this file is duplicated here because doctest cannot pick up custom
# fixtures when the conftest file is anywhere else for some reason
import os

import pytest

import fiasco
from fiasco.util import build_hdf5_dbase, download_dbase
from fiasco.util.setup_db import CHIANTI_URL, LATEST_VERSION

# Minimal set of CHIANTI files needed to run the doctests
TEST_FILES = {
    'sun_photospheric_1998_grevesse.abund': '9b175ee91f80fbe01967321a0fb051a8',
    'chianti.ip':                           'a5a5071535f14590210957f8783f2843',
    'chianti.ioneq':                        '81d8d24bb09cb2da63befd30e6c8767c',
    'fe_5.elvlc':                           'f7dd75c6bcea77584c5e00ddf2684d9c',
    'fe_15.elvlc':                          'cb8e6291bac17325130ea8564e118d48',
    'fe_18.rrparams':                       '6a36e5db1e377d0c30652960eb77b93d',
    'fe_18.drparams':                       '1e8b0194aeb099b2653508110e388200',
}


@pytest.fixture(scope='session', autouse=True)
def ascii_dbase_root(tmpdir_factory):
    # If we already have a local copy, just return the path to that
    path = fiasco.defaults.get('test_ascii_dbase_root')
    if path is not None and os.path.exists(path):
        return path

    # Otherwise download the database
    path = tmpdir_factory.mktemp('chianti_dbase')
    download_dbase(CHIANTI_URL.format(version=LATEST_VERSION), path)
    return path


@pytest.fixture(scope='session', autouse=True)
def hdf5_dbase_root(ascii_dbase_root, tmpdir_factory):
    path = tmpdir_factory.mktemp('fiasco').join('chianti_dbase.h5')
    build_hdf5_dbase(ascii_dbase_root, path, files=TEST_FILES)
    return path
