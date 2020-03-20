# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.
import os
import pytest
from astropy.tests.plugins import *

import fiasco
from fiasco.util import build_hdf5_dbase, download_dbase
from fiasco.util.setup_db import CHIANTI_URL, LATEST_VERSION, FIASCO_HOME

# Minimal set of CHIANTI files needed to run the tests
TEST_FILES = {
    'sun_photospheric_1998_grevesse.abund': '9b175ee91f80fbe01967321a0fb051a8',
    'chianti.ip':                           'a5a5071535f14590210957f8783f2843',
    'chianti.ioneq':                        '81d8d24bb09cb2da63befd30e6c8767c',
    'h_1.elvlc':                            'd31620aaf26a14486e635d14bcf7a6c1',
    'h_1.wgfa':                             'ed2a561ecdba5ee4b05ea56c297724ba',
    'h_1.scups':                            'de180c9e1b4f50a503efad8c83e714ab',
    'h_1.diparams':                         'f9be79794cc092c75733ece7aba9f29f',
    'h_2.rrparams':                         '05eb5044dc1ad070338d1ba3745dc4de',
    'he_1.elvlc':                           '577245da46cfc4d27a05f40147f17610',
    'he_2.elvlc':                           '58eee740f6842850f1f35a4f95b3e12c',
    'he_3.rrparams':                        '67e2e7a5e86c8eaa5db4e95ac3d35989',
    'li_1.diparams':                        '2a4b60bc799a1ea3ee3e911b5968733e',
    'ca_2.elvlc':                           'c8a9dbdc622d0a15eea6cf0eb1631d3b',
    'fe_2.elvlc':                           '9306d594a66b4649363a778a650853f9',
    'fe_5.elvlc':                           'f7dd75c6bcea77584c5e00ddf2684d9c',
    'fe_6.elvlc':                           '081519f986b8a8ed99a34ecf813f1358',
    'fe_9.elvlc':                           'b7d04f65a87a8de1c2bfc77331e438f3',
    'fe_12.elvlc':                          'ee9beb8b2ff03ba8fc046bd722992b21',
    'fe_15.elvlc':                          'cb8e6291bac17325130ea8564e118d48',
    'fe_18.rrparams':                       '6a36e5db1e377d0c30652960eb77b93d',
    'fe_18.drparams':                       '1e8b0194aeb099b2653508110e388200',
    'fe_27.rrparams':                       '75383b0f1b167f862cfd26bbadd2a029',
}


@pytest.fixture(scope='session')
def ascii_dbase_root(tmpdir_factory):
    # If we already have a local copy, just return the path to that
    path = fiasco.defaults.get('test_ascii_dbase_root')
    if path is not None and os.path.exists(path):
        return path

    # Otherwise download the database
    path = tmpdir_factory.mktemp('chianti_dbase')
    download_dbase(CHIANTI_URL.format(version=LATEST_VERSION), path)
    return path


@pytest.fixture(scope='session')
def hdf5_dbase_root(ascii_dbase_root, tmpdir_factory):
    path = tmpdir_factory.mktemp('fiasco').join('chianti_dbase.h5')
    build_hdf5_dbase(ascii_dbase_root, path, files=TEST_FILES)
    return path
