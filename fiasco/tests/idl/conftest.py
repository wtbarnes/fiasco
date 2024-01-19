import os
import pytest

from astropy.utils.data import get_pkg_data_path


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
    try:
        import hissw
    except ImportError:
        return None
    extra_paths = [d for d, _, _ in os.walk(idl_codebase_root)]
    header = f'''
    defsysv,'!xuvtop','{ascii_dbase_root}'
    defsysv,'!abund_file',''
    defsysv,'!ioneq_file',''
    '''
    extra_paths += [get_pkg_data_path('ssw_gen_functions', package='fiasco.tests.idl')]
    env = hissw.Environment(
        idl_only=True,
        idl_home=idl_executable,
        header=header,
        extra_paths=extra_paths
    )
    try:
        _ = env.run('print,!xuvtop')
    except (hissw.util.SSWIDLError, hissw.util.IDLLicenseError):
        return None
    else:
        return env
