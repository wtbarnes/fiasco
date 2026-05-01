"""
Helpers for tests related to comparing with IDL output
"""
import asdf
import numpy as np
import os
import pathlib

from astropy.utils.data import get_pkg_data_path
from packaging.version import Version

import fiasco

from fiasco.util import read_chianti_version

__all__ = [
    'get_idl_test_output_filepath',
    'read_idl_test_output',
    'version_check',
    'setup_idl_environment',
    'run_idl_script',
    'get_chianti_idl_version',
]


def get_idl_test_output_filepath(name, version):
    if not isinstance(version, Version):
        version = Version(version)
    data_dir = pathlib.Path(get_pkg_data_path('data', package='fiasco.tests.idl'))
    return data_dir / f"{name}_v{version.major}.{version.minor}.{version.micro}.asdf"


def read_idl_test_output(name, version, keys=None):
    file_path = get_idl_test_output_filepath(name, version)
    with asdf.open(file_path, memmap=False, lazy_load=False) as af:
        if keys is None:
            keys = af.tree.keys()
        output = {k: af.tree[k] for k in keys}
    return output


def version_check(current_version, conditional, version):
    """
    Filter for easily evaluating conditional on version
    """
    op_dict = {'<': np.less,
               '<=': np.less_equal,
               '>': np.greater,
               '>=': np.greater_equal,
               '=': np.equal,
               '==': np.equal,
               '!=': np.not_equal}
    return op_dict[conditional](
        Version(str(current_version)),
        Version(str(version)),
    )


def setup_idl_environment(ascii_dbase_root,
                          idl_codebase_root,
                          idl_executable,
                          ssw_gen_root=None):
    """
    Setup an IDL environment for CHIANTI

    .. note:: If `hissw` is not installed or IDL or SSW are not installed, this
              function will return None.

    Parameters
    ----------
    ascii_dbase_root: path-like
        Path to top of IDL database tree
    idl_codebase_root: path-like
        Path to top of CHIANTI IDL software tree
    idl_executable: path-like
        Path to IDL executable
    """
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
    if ssw_gen_root is None:
        extra_paths += [get_pkg_data_path('ssw_gen_functions', package='fiasco.tests.idl')]
    else:
        extra_paths += [d for d, _, _ in os.walk(ssw_gen_root)]
    env = hissw.Environment(
        idl_only=True,
        idl_home=idl_executable,
        header=header,
        extra_paths=extra_paths,
        filters={'version_check': version_check},
    )
    try:
        _ = env.run('print,!xuvtop', verbose=False)
    except (hissw.util.SSWIDLError, hissw.util.IDLLicenseError):
        return None
    else:
        return env


def run_idl_script(idl_env,
                   script,
                   input_args,
                   save_vars,
                   file_name,
                   dbase_version,
                   chianti_idl_version,
                   format_func=None,
                   write_file=True):
    """
    Helper function for running CHIANTI IDL via hissw in tests

    This function runs the IDL script with hissw and saves the output to an asdf file
    such that tests can be run to compare against the IDL results without actually rerunning
    the IDL code. If the output file already exists, then the results are read from that file.

    Parameters
    ----------
    idl_env: `hissw.Environment`
        IDL environment to use to run the test script
    script: `str`
        IDL script to run to generate test comparison data
    input_args: `dict`
        All inputs to the IDL script
    save_vars: `list`
        List of variables to return from the IDL calculation
    file_name: `str`
        Name of the IDL results file
    dbase_version: `packaging.version.Version`, `str`
        Version of CHIANTI database used to generate these test results
    chianti_idl_version: `packaging.version.Version`, `str`
        Version of CHIANTI IDL software used to generate these test results
    format_func: `dict`, optional
        Functions to use to format output from the IDL function.
        This is most useful for adding the necessary units to any of the outputs.
    write_file: `bool`, optional
        If True, writes the result to a file so that the test can be run later without
        an IDL installation. If False, running this test will always require an IDL
        installation. Setting this to False may be needed to reduce the amount of test
        data.
    """
    file_path = get_idl_test_output_filepath(file_name, dbase_version)
    if not file_path.is_file():
        if idl_env is None:
            # Import here so that this can be used without a hard pytest dependency
            import pytest
            pytest.skip("""To run the IDL comparison tests, you must:
                            1. Specify a path to a working IDL executable,
                            2. Specify a path to the CHIANTI IDL code,
                            3. Install the hissw package (pip install hissw)
                           Without the following, you will not be able to generate new IDL comparison test results.""")
        # Add versions so they are accessible from within the script if needed
        input_args = {
            **input_args,
            'database_version': str(dbase_version),
            'chianti_idl_version': str(chianti_idl_version)
        }
        result = idl_env.run(script, args=input_args, save_vars=save_vars)
        if format_func is not None:
            for k in format_func:
                result[k] = format_func[k](result[k])
        if any([k in result for k in input_args]):
            fiasco.log.warning(
                'Conflict between input variables and result variables. Outputs to be unintentionally overwritten'
            )
        variables = {
            **result,
            **input_args,
            'idl_script': idl_env.render_script(script, input_args),
        }
        if write_file:
            with asdf.AsdfFile(variables) as af:
                af.write_to(file_path)
        else:
            return variables
    return read_idl_test_output(file_name, dbase_version)


def get_chianti_idl_version(idl_codebase_root):
    """
    Get the version of the CHIANTI IDL code being used.

    Get the current version of the CHIANTI IDL code used
    for generating test data so that it can be recorded in
    the results files.Try reading version from the current
    git tag as this is most reliable. If that fails, try to
    read it from the ``VERSION`` file. If that fails, raise a
    warning and return None.
    """
    # First try to get the version from the current git tag
    try:
        import git
    except ImportError:
        fiasco.log.warning('Cannot determine CHIANTI IDL version from git tag.')
    else:
        # NOTE: Nesting these exceptions like this because catching both
        # exceptions at once is not impossible as the second is imported
        # from a package that may be missing.
        try:
            tag = git.Repo(idl_codebase_root).git.describe('--tags')
        except git.exc.InvalidGitRepositoryError:
            fiasco.log.warning('Cannot determine CHIANTI IDL version from git tag.')
        else:
            return Version(tag)
    # Next try to get the version from the version file. This is secondary
    # because this file is not always present or accurate, particularly for
    # older versions. This file follows the same naming conventions and format
    # as the database version file.
    try:
        version = read_chianti_version(idl_codebase_root)
    except FileNotFoundError:
        fiasco.log.warning('Cannot determine CHIANTI IDL version from VERSION file')
    else:
        return version
    fiasco.log.warning('Cannot determine CHIANTI IDL version. Returning None.')
    return None
