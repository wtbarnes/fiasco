import pytest

from .helpers import get_chianti_idl_version, setup_idl_environment


@pytest.fixture(scope="session")
def idl_env(ascii_dbase_root, request):
    idl_executable = request.config.getoption('--idl-executable', default=None)
    idl_codebase_root = request.config.getoption('--idl-codebase-root', default=None)
    return setup_idl_environment(ascii_dbase_root, idl_codebase_root, idl_executable)


@pytest.fixture(scope="session")
def chianti_idl_version(request):
    idl_codebase_root = request.config.getoption('--idl-codebase-root', default=None)
    return get_chianti_idl_version(idl_codebase_root)
