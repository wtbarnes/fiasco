# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.
import os
import pytest
from astropy.tests.plugins import *

from fiasco.util import build_hdf5_dbase, download_dbase
from fiasco.util.setup_db import CHIANTI_URL, LATEST_VERSION

# Minimal set of CHIANTI files needed to run the tests
TEST_FILES = [
    'sun_photospheric_1998_grevesse.abund',
    'chianti.ip',
    'chianti.ioneq',
    'h_1.elvlc',
    'h_1.wgfa',
    'h_1.scups',
    'h_2.rrparams',
    'he_1.elvlc',
    'he_2.elvlc',
    'he_3.rrparams',
    'li_1.diparams',
    'ca_2.elvlc',
    'fe_2.elvlc',
    'fe_5.elvlc',
    'fe_6.elvlc',
    'fe_9.elvlc',
    'fe_12.elvlc',
    'fe_15.elvlc',
    'fe_18.rrparams',
    'fe_18.drparams',
    'fe_27.rrparams',
]


@pytest.fixture(scope='session')
def ascii_dbase_root(tmpdir_factory):
    path = tmpdir_factory.mktemp('chianti_dbase')
    download_dbase(CHIANTI_URL.format(version=LATEST_VERSION), path)
    return path


@pytest.fixture(scope='session')
def hdf5_dbase_root(ascii_dbase_root, tmpdir_factory):
    path = tmpdir_factory.mktemp('fiasco').join('chianti_dbase.h5')
    build_hdf5_dbase(ascii_dbase_root, path, files=TEST_FILES)
    return path
