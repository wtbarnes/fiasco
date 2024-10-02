import numpy as np
import pathlib
import pytest

from packaging.version import Version

from fiasco.tests import get_test_file_list
from fiasco.util import check_database, read_chianti_version

# Force MPL to use non-gui backends for testing.
try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')

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
        test_files = get_test_file_list()
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
    """
    Skip a test if not all version requirements are met. Multiple requirements are joined by AND.
    """
    # NOTE: Fixtures that depend on other fixtures are awkward to implement.
    # See this SO answer: https://stackoverflow.com/a/28198398
    if marker := request.node.get_closest_marker('requires_dbase_version'):
        op_dict = {'<': np.less,
                   '<=': np.less_equal,
                   '>': np.greater,
                   '>=': np.greater_equal,
                   '=': np.equal,
                   '==': np.equal,
                   '!=': np.not_equal}

        def _evaluate_condtion(condition_string):
            condition_array = condition_string.split()
            if  len(condition_array) != 2:
                raise ValueError("Arguments must contain a condition and a version number with a space, e.g. '< 8.0.7'")
            operator, target_version = condition_array
            if operator not in op_dict:
                raise ValueError(f'''{operator} is not a supported comparison operation.
                                    Must be one of {list(op_dict.keys())}.''')
            target_version = Version(target_version)
            allowed_dbase_version = op_dict[operator](dbase_version, target_version)
            return allowed_dbase_version, operator, target_version

        conditions = np.atleast_1d(marker.args)
        for is_met, operator, target_version in list(map(_evaluate_condtion, conditions)):
            if not is_met:
                pytest.skip(f'Skipping because database version {dbase_version} is not {operator} {target_version}.')


def pytest_configure(config):
  config.addinivalue_line(
        "markers", "requires_dbase_version(dbase_version): Skip tests based on CHIANTI database version requirements.",
  )
