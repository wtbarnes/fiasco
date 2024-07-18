import pytest

from .helpers import setup_idl_environment


@pytest.fixture(scope="session")
def idl_env(ascii_dbase_root, request):
    # NOTE: The reason that we return None here rather than calling pytest.skip
    # is that the IDL tests can still be run if there is a file containing the IDL
    # results available. Thus, we have to wait until later to decide if we want to
    # skip the test.
    idl_executable = request.config.getoption('--idl-executable')
    idl_codebase_root = request.config.getoption('--idl-codebase-root')
    if idl_executable is None or idl_codebase_root is None:
        return None
    return setup_idl_environment(ascii_dbase_root, idl_codebase_root, idl_executable)
