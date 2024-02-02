import numpy as np
import pathlib
import pytest

from packaging.version import Version

from fiasco.util import check_database, read_chianti_version

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
    'sun_coronal_1992_feldman_ext.abund',
    'sun_coronal_1992_feldman.abund',
    'sun_photospheric_2007_grevesse.abund',
    'chianti.ip',
    'chianti.ioneq',
    'gffgu.dat',
    'gffint.dat',
    'klgfb.dat',
    'itoh.dat',
    'hseq_2photon.dat',
    'heseq_2photon.dat',
    'verner_short.txt',
    'flare.dem',
    'h_1.elvlc',
    'h_1.wgfa',
    'h_1.scups',
    'h_1.diparams',
    'h_2.rrparams',
    'he_1.elvlc',
    'he_2.elvlc',
    'he_3.rrparams',
    'c_1.diparams',
    'c_2.diparams',
    'c_2.drparams',
    'c_2.rrparams',
    'c_3.diparams',
    'c_3.drparams',
    'c_3.easplups',
    'c_3.rrparams',
    'c_4.diparams',
    'c_4.drparams',
    'c_4.easplups',
    'c_4.rrparams',
    'c_5.diparams',
    'c_5.drparams',
    'c_5.rrparams',
    'c_6.diparams',
    'c_6.drparams',
    'c_6.rrparams',
    'c_6.elvlc',
    'c_6.wgfa',
    'c_6.scups',
    'c_6.psplups',
    'c_7.rrparams',
    'ca_4.rrparams',
    'ca_4.drparams',
    'ca_15.elvlc',
    'ca_15.wgfa',
    'ca_15.psplups',
    'ca_15.scups',
    'li_1.diparams',
    'ca_2.elvlc',
    'fe_2.elvlc',
    'fe_5.elvlc',
    'fe_5.diparams',
    'fe_5.drparams',
    'fe_5.rrparams',
    'fe_5.trparams',
    'fe_5.easplups',
    'fe_5.scups',
    'fe_5.wgfa',
    'fe_6.elvlc',
    'fe_9.scups',
    'fe_9.wgfa',
    'fe_9.elvlc',
    'fe_9.fblvl',
    'fe_11.elvlc',
    'fe_11.wgfa',
    'fe_11.scups',
    'fe_11.psplups',
    'fe_12.elvlc',
    'fe_14.elvlc',
    'fe_14.wgfa',
    'fe_14.scups',
    'fe_14.psplups',
    'fe_15.elvlc',
    'fe_16.diparams',
    'fe_16.drparams',
    'fe_16.easplups',
    'fe_16.elvlc',
    'fe_16.scups',
    'fe_16.wgfa',
    'fe_16.rrparams',
    'fe_18.elvlc',
    'fe_18.scups',
    'fe_18.wgfa',
    'fe_27.rrparams',
    'fe_10.psplups',
    'fe_10.elvlc',
    'fe_20.elvlc',
    'fe_20.wgfa',
    'fe_20.scups',
    'fe_20.cilvl',
    'fe_20.reclvl',
    'he_2.fblvl',
    'he_1.fblvl',
    'ni_28.fblvl',
    'ni_17.fblvl',
    'ni_21.fblvl',
    'ni_19.fblvl',
    'ni_26.fblvl',
    'ni_18.fblvl',
    'ni_27.fblvl',
    'ni_20.fblvl',
    'ni_16.fblvl',
    'ni_11.fblvl',
    'ni_13.fblvl',
    'ni_25.fblvl',
    'ni_22.fblvl',
    'ni_23.fblvl',
    'ni_24.fblvl',
    'ni_12.fblvl',
    'ni_15.fblvl',
    'al_13.fblvl',
    'al_12.fblvl',
    'al_9.fblvl',
    'al_7.fblvl',
    'al_1.fblvl',
    'al_6.fblvl',
    'al_8.fblvl',
    'al_10.fblvl',
    'al_11.fblvl',
    'al_3.fblvl',
    'al_4.fblvl',
    'al_5.fblvl',
    'mg_7.fblvl',
    'mg_9.fblvl',
    'mg_8.fblvl',
    'mg_1.fblvl',
    'mg_6.fblvl',
    'mg_10.fblvl',
    'mg_11.fblvl',
    'mg_3.fblvl',
    'mg_4.fblvl',
    'mg_5.fblvl',
    'mg_2.fblvl',
    'mg_12.fblvl',
    'ca_17.fblvl',
    'ca_10.fblvl',
    'ca_19.fblvl',
    'ca_20.fblvl',
    'ca_18.fblvl',
    'ca_11.fblvl',
    'ca_16.fblvl',
    'ca_13.fblvl',
    'ca_14.fblvl',
    'ca_15.fblvl',
    'ca_12.fblvl',
    'ca_9.fblvl',
    'fe_10.fblvl',
    'fe_17.fblvl',
    'fe_21.fblvl',
    'fe_26.fblvl',
    'fe_19.fblvl',
    'fe_18.fblvl',
    'fe_20.fblvl',
    'fe_16.fblvl',
    'fe_11.fblvl',
    'fe_5.fblvl',
    'fe_4.fblvl',
    'fe_14.fblvl',
    'fe_13.fblvl',
    'fe_25.fblvl',
    'fe_22.fblvl',
    'fe_23.fblvl',
    'fe_24.fblvl',
    'fe_12.fblvl',
    'fe_15.fblvl',
    'fe_8.fblvl',
    'fe_6.fblvl',
    'fe_7.fblvl',
    'fe_9.fblvl',
    'n_3.fblvl',
    'n_4.fblvl',
    'n_5.fblvl',
    'n_2.fblvl',
    'n_7.fblvl',
    'n_1.fblvl',
    'n_6.fblvl',
    's_6.fblvl',
    's_1.fblvl',
    's_8.fblvl',
    's_9.fblvl',
    's_7.fblvl',
    's_12.fblvl',
    's_15.fblvl',
    's_14.fblvl',
    's_13.fblvl',
    's_2.fblvl',
    's_5.fblvl',
    's_4.fblvl',
    's_3.fblvl',
    's_16.fblvl',
    's_11.fblvl',
    's_10.fblvl',
    'o_7.fblvl',
    'o_6.fblvl',
    'o_1.fblvl',
    'o_8.fblvl',
    'o_4.fblvl',
    'o_3.fblvl',
    'o_2.fblvl',
    'o_5.fblvl',
    'h_1.fblvl',
    'si_12.fblvl',
    'si_14.fblvl',
    'si_13.fblvl',
    'si_2.fblvl',
    'si_5.fblvl',
    'si_4.fblvl',
    'si_3.fblvl',
    'si_11.fblvl',
    'si_10.fblvl',
    'si_6.fblvl',
    'si_1.fblvl',
    'si_8.fblvl',
    'si_9.fblvl',
    'si_7.fblvl',
    'ar_10.fblvl',
    'ar_17.fblvl',
    'ar_16.fblvl',
    'ar_11.fblvl',
    'ar_18.fblvl',
    'ar_4.fblvl',
    'ar_8.fblvl',
    'ar_1.fblvl',
    'ar_14.fblvl',
    'ar_13.fblvl',
    'ar_7.fblvl',
    'ar_9.fblvl',
    'ar_12.fblvl',
    'ar_15.fblvl',
    'ne_9.fblvl',
    'ne_7.fblvl',
    'ne_6.fblvl',
    'ne_1.fblvl',
    'ne_8.fblvl',
    'ne_4.fblvl',
    'ne_10.fblvl',
    'ne_3.fblvl',
    'ne_2.fblvl',
    'ne_5.fblvl',
    'c_5.fblvl',
    'c_2.fblvl',
    'c_3.fblvl',
    'c_4.fblvl',
    'c_1.fblvl',
    'c_6.fblvl',
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
        test_files = TEST_FILES
    check_hash = not request.config.getoption('--disable-file-hash')
    # Optionally use a different URL for the database (e.g. for testing different versions)
    ascii_dbase_url = request.config.getoption('--ascii-dbase-url')
    # Optionally disable version checking. This is useful when testing newer/older versions
    # of CHIANTI
    check_chianti_version = not request.config.getoption('--skip-version-check')
    # Finally, run the database setup
    check_database(path,
                   ascii_dbase_root=ascii_dbase_tree,
                   ask_before=False,
                   ascii_dbase_url=ascii_dbase_url,
                   check_chianti_version=check_chianti_version,
                   files=test_files,
                   check_hash=check_hash)
    return path


@pytest.fixture(scope='session')
def ascii_dbase_root(ascii_dbase_tree, hdf5_dbase_root):
    # The reason this exists is to ensure that the ASCII database is downloaded (if needed)
    # which happens in the above fixture. The downloading of the database was originally
    # in this fixture but that resulted in just recreating the check_database function
    # but in a worse way. What this means is that the HDF5 database is always built, even
    # when testing the ASCII parsing code, but that is not a huge price to pay.
    return ascii_dbase_tree


@pytest.fixture(scope="session")
def dbase_version(ascii_dbase_root):
    return read_chianti_version(ascii_dbase_root)


@pytest.fixture(autouse=True)
def requires_dbase_version(request, dbase_version):
    # NOTE: Fixtures that depend on other fixtures are awkward to implement.
    # See this SO answer: https://stackoverflow.com/a/28198398
    if marker := request.node.get_closest_marker('requires_dbase_version'):
        # NOTE: This has to have a space between the operator and the target
        if  len(marker.args) != 2:
            raise ValueError("Arguments must contain a condition and a version number, e.g. '<', '8.0.7'")
        operator, target_version = marker.args
        op_dict = {'<': np.less,
                   '<=': np.less_equal,
                   '>': np.greater,
                   '>=': np.greater_equal,
                   '=': np.equal,
                   '==': np.equal,
                   '!=': np.not_equal}
        if operator not in op_dict:
            raise ValueError(f'''{operator} is not a supported comparison operation.
                                 Must be one of {list(op_dict.keys())}.''')
        target_version = Version(target_version)
        allowed_dbase_version = op_dict[operator](dbase_version, target_version)
        if not allowed_dbase_version:
            pytest.skip(f'Skip because database version {dbase_version} is not {operator} {target_version}.')


def pytest_configure(config):
  config.addinivalue_line(
        "markers", "requires_dbase_version(dbase_version): Skip tests based on CHIANTI database version requirements.",
  )
