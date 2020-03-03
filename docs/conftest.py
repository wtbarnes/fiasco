# NOTE: this file is duplicated here because doctest cannot pick up custom
# fixtures when the conftest file is anywhere else for some reason
import pytest

from fiasco.util import build_hdf5_dbase, download_dbase
from fiasco.util.setup_db import CHIANTI_URL, LATEST_VERSION

# Minimal set of CHIANTI files needed to run the doctests
TEST_FILES = [
    'sun_photospheric_1998_grevesse.abund',
    'chianti.ip',
    'chianti.ioneq',
    'fe_5.elvlc',
    'fe_15.elvlc',
    'fe_18.rrparams',
    'fe_18.drparams',
]


@pytest.fixture(scope='session', autouse=True)
def ascii_dbase_root(tmpdir_factory):
    path = tmpdir_factory.mktemp('chianti_dbase')
    download_dbase(CHIANTI_URL.format(version=LATEST_VERSION), path)
    return path


@pytest.fixture(scope='session', autouse=True)
def hdf5_dbase_root(ascii_dbase_root, tmpdir_factory):
    path = tmpdir_factory.mktemp('fiasco').join('chianti_dbase.h5')
    build_hdf5_dbase(ascii_dbase_root, path, files=TEST_FILES)
    return path
