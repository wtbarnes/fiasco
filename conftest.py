# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.
import os
from urllib.request import urlopen

import pytest

from fiasco.util import build_hdf5_dbase, download_dbase
from fiasco.util.setup_db import CHIANTI_URL, LATEST_VERSION

# Minimal set of CHIANTI files needed to run the tests
# NOTE: need some way for this to be flexible depending on the supplied database version
TEST_FILES = {
    'sun_photospheric_1998_grevesse.abund': '9b175ee91f80fbe01967321a0fb051a8',
    'chianti.ip':                           'a5a5071535f14590210957f8783f2843',
    'chianti.ioneq':                        '81d8d24bb09cb2da63befd30e6c8767c',
    'gffgu.dat':                            '5895e1f211e7baa8293478df014b9ff1',
    'gffint.dat':                           '96381825257793805ca42eba03ae83e2',
    'klgfb.dat':                            '5496daa096dd4f7b935ec0d0e6417084',
    'itoh.dat':                             'e40aacc9d06889f9c4ba63a92673a398',
    'hseq_2photon.dat':                     '48656984fbcbe38f883ff9a4460c790a',
    'heseq_2photon.dat':                    '9e42ac8c37d67ba3109aaa56aab9e736',
    'verner_short.txt':                     'f05296ffc9c5306846ac9caf136fad26',

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
    'fe_5.fblvl':                           'e13a9becff489ea044fbd3cb1e02af13',
    'fe_5.drparams':                        '0f47e2fbf4be9eb5564a8adc19531999',
    'fe_5.easplups':                        'ec5443d429722ccdd9903743ca4718ab',
    'fe_5.scups':                           '0a863212d59c5eda5936e0664b513fb9',
    'fe_5.wgfa':                            '6f8e4d41760d5a0540008e15aca1038c',

    'fe_6.elvlc':                           '081519f986b8a8ed99a34ecf813f1358',
    'fe_6.fblvl':                           '9478aeab2479e9eefebfd1e552e27ddf',
    'fe_9.elvlc':                           'b7d04f65a87a8de1c2bfc77331e438f3',
    'fe_12.elvlc':                          'ee9beb8b2ff03ba8fc046bd722992b21',
    'fe_15.elvlc':                          'cb8e6291bac17325130ea8564e118d48',
    'fe_18.rrparams':                       '6a36e5db1e377d0c30652960eb77b93d',
    'fe_18.drparams':                       '1e8b0194aeb099b2653508110e388200',
    'fe_27.rrparams':                       '75383b0f1b167f862cfd26bbadd2a029',
    'fe_10.psplups':                        'dd34363f6daa81dbf106fbeb211b457d',
    'fe_10.elvlc':                          'f221d4c7167336556d57378ac368afc1',
}


def pytest_addoption(parser):
    parser.addoption('--ascii-dbase-root', action='store', default=None)


@pytest.fixture(scope='session')
def ascii_dbase_root(tmpdir_factory, request):
    path = request.config.getoption('--ascii-dbase-root')
    if path is None:
        path = CHIANTI_URL.format(version=LATEST_VERSION)
    try:
        _ = urlopen(path)
    except ValueError:
        if not os.path.exists(path):
            raise ValueError(f'{path} is not a valid URL or file path')
    else:
        _path = tmpdir_factory.mktemp('chianti_dbase')
        download_dbase(path, _path)
        path = _path
    return path


@pytest.fixture(scope='session')
def hdf5_dbase_root(ascii_dbase_root, tmpdir_factory):
    path = tmpdir_factory.mktemp('fiasco').join('chianti_dbase.h5')
    # FIXME: Disable hash checking for now
    build_hdf5_dbase(ascii_dbase_root, path, files=[k for k in TEST_FILES])
    return path
