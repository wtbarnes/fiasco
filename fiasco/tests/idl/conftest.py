import os
import pytest

from astropy.utils.data import get_pkg_data_path


@pytest.fixture(scope="session")
def idl_env(ascii_dbase_root, request):
    idl_executable = request.config.getoption('--idl-executable')
    idl_codebase_root = request.config.getoption('--idl-codebase-root')
    if idl_codebase_root is None or idl_codebase_root is None:
        pytest.skip('''Must specify path to IDL executable and CHIANTI IDL code in order to
                       run the IDL comparison tests.''')
    try:
        import hissw
    except ImportError:
        pytest.skip('''The hissw package is required to install the IDL comparison tests.
                       See https://github.com/wtbarnes/hissw for installation information.''')
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
        pytest.skip('''Could not configure CHIANTI IDL environment.
                       You will not be able to run portions of the test suite.''')
    else:
        return env
