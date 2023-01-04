import pathlib
import pytest

from fiasco.util import check_database

# Force MPL to use non-gui backends for testing.
try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')


# Minimal set of CHIANTI files needed to run the tests
# NOTE: need some way for this to be flexible depending on the supplied database version
TEST_FILES = {
    'sun_coronal_1992_feldman_ext.abund':   '8b18d62e03528e3b3800806a1fc9391e',
    'sun_coronal_1992_feldman.abund':       '75fde4f73bae8fdce136f65b736bc0c8',
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
    'c_1.diparams':                         '18e24ef4119abccd7226871045ebf769',
    'c_2.diparams':                         '499271a2634a66e6d0717f97cf1f6e7a',
    'c_2.drparams':                         'b9adf6474864b7d480c2132a95099358',
    'c_2.rrparams':                         '04c6cab3710f847f40e751e23fe29a54',
    'c_3.diparams':                         'e573d3ad1b34f86ef2896813a8452d4e',
    'c_3.drparams':                         '3d1d3b00866e77d51afd1206a6c63989',
    'c_3.easplups':                         '8281d99bfd4d9165302ca9b27c219d41',
    'c_3.rrparams':                         '8122c689559acd7715bd12a86672fb8f',
    'c_4.diparams':                         'a6f236513ffd130bfbab9a222c52e451',
    'c_4.drparams':                         'd2435409c5a86c4ff4e9687255d2d9f4',
    'c_4.easplups':                         'e9574be2d6434adcebbead5efdd6b638',
    'c_4.rrparams':                         '1bcd165f366b6e48c550ff4b841eb4d7',
    'c_5.diparams':                         '8087926a840a5da8c20a9eb36645316a',
    'c_5.drparams':                         'ff50ded77ec19c7a22f057663fff5c22',
    'c_5.rrparams':                         'e35b12f4fedc6326c91f31a22e037e88',
    'c_6.diparams':                         '6bdf83507b516f50f9561cc253d735d4',
    'c_6.drparams':                         '56503cbd1ed960101eefeaff77b41392',
    'c_6.rrparams':                         '811014c5b1ee782e12bab6ea53f52f04',
    'c_7.rrparams':                         '686156392d4cfabb3f0ff0280970f562',
    'ca_4.rrparams':                        '35f3713b90e6be7de6ca36c558546578',
    'ca_4.drparams':                        '7cbb34c0041aba833c389256863b66cc',
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
    'fe_16.diparams':                       '8c7733d443eceab3e2fd4073bacfa9a5',
    'fe_16.drparams':                       '5067a1e2ed714c673681ac8e6d8325f7',
    'fe_16.easplups':                       '9de5fcdb0814dbb8e51a4aad18273229',
    'fe_16.rrparams':                       'f5adaae054fd083a3a41d6bdddb06fa8',
    'fe_18.elvlc':                          'e5384e1f58d689bd29f83eede58bc290',
    'fe_18.scups':                          'ec57a94f781a3b6b553de90679e0a07c',
    'fe_18.wgfa':                           'f58accd58b7de6342ee174bb844642b8',
    'fe_27.rrparams':                       '75383b0f1b167f862cfd26bbadd2a029',
    'fe_10.psplups':                        'dd34363f6daa81dbf106fbeb211b457d',
    'fe_10.elvlc':                          'f221d4c7167336556d57378ac368afc1',
}


@pytest.fixture(scope='session')
def ascii_dbase_tree(tmpdir_factory, request):
    path = request.config.getoption('--ascii-dbase-root')
    if path is None:
        path = tmpdir_factory.mktemp('chianti_dbase')
    return pathlib.Path(path)


@pytest.fixture(scope='session')
def hdf5_dbase_root(ascii_dbase_tree, tmpdir_factory, request):
    # If specifying the HDF5 database explicitly, do not do any setup
    # as it is assumed that this database has already been built
    # Note that this test setup explicitly does not build the test HDF5
    # database in a custom location. It is always in a temporary directory
    # unless you explicitly build it ahead of time.
    path = request.config.getoption('--hdf5-dbase-root')
    if path is not None:
        return path
    # Otherwise, set it up here
    path = tmpdir_factory.mktemp('fiasco').join('chianti_dbase.h5')
    # Setup the test files. By default, only a limited number of files
    # are included in the test database to make running the tests more
    # efficient and a hash is checked to ensure the correct tests are
    # being checked
    if request.config.getoption('--include-all-files'):
        test_files = None
    else:
        if request.config.getoption('--disable-file-hash'):
            test_files = [k for k in TEST_FILES]
        else:
            test_files = TEST_FILES
    # Optionally use a different URL for the database (e.g. for testing different versions)
    ascii_dbase_url = request.config.getoption('--ascii-dbase-url')
    # Finally, run the database setup
    check_database(path,
                   ascii_dbase_root=ascii_dbase_tree,
                   ask_before=False,
                   ascii_dbase_url=ascii_dbase_url,
                   check_chianti_version=False,
                   files=test_files)
    return path


@pytest.fixture(scope='session')
def ascii_dbase_root(ascii_dbase_tree, hdf5_dbase_root):
    # The reason this exists is to ensure that the ASCII database is downloaded (if needed)
    # which happens in the above fixture. The downloading of the database was originally
    # in this fixture but that resulted in just recreating the check_database function
    # but in a worse way. What this means is that the HDF5 database is always built, even
    # when testing the ASCII parsing code, but that is not a huge price to pay.
    return ascii_dbase_tree
