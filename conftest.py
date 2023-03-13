"""
Define command line options for running tests
"""


def pytest_addoption(parser):
    parser.addoption('--ascii-dbase-root', action='store', default=None)
    parser.addoption('--ascii-dbase-url', action='store', default=None)
    parser.addoption('--hdf5-dbase-root', action='store', default=None)
    parser.addoption('--disable-file-hash', action='store_true', default=False,
                     help='Disable MD5 hash checks on test files')
    parser.addoption('--idl-executable', action='store', default=None)
    parser.addoption('--idl-codebase-root', action='store', default=None)
    parser.addoption('--include-all-files', action='store_true', default=False)
    parser.addoption('--skip-version-check', action='store_true', default=False,
                     help='Do not check CHIANTI version')
