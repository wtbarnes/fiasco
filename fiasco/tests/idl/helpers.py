"""
Helpers for tests related to comparing with IDL output
"""
import pathlib
import pytest

from astropy.utils.data import get_pkg_data_path

__all__ = ['get_idl_test_output_filepath', 'read_idl_test_output', 'run_idl_script']


def get_idl_test_output_filepath(name, version):
    data_dir = pathlib.Path(get_pkg_data_path('data', package='fiasco.tests.idl'))
    return data_dir / f"{name}_v{version}.asdf"


def read_idl_test_output(name, version, keys=None):
    # NOTE: Importing here to avoid it as a hard dependency for running tests
    import asdf
    file_path = get_idl_test_output_filepath(name, version)
    with asdf.open(file_path, copy_arrays=True) as af:
        if keys is None:
            keys = af.tree.keys()
        output = {k: af.tree[k] for k in keys}
    return output


def run_idl_script(idl_env, script, input_args, save_vars, file_name, version, format_func=None, write_file=True):
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
    version: `packaging.version.Version`, `str`
        Version of IDL used to generate these test results
    format_func: `dict`, optional
        Functions to use to format output from the IDL function.
        This is most useful for adding the necessary units to any of the outputs.
    write_file: `bool`, optional
        If True, writes the result to a file so that the test can be run later without
        an IDL installation. If False, running this test will always require an IDL
        installation. Setting this to False may be needed to reduce the amount of test
        data.
    """
    # NOTE: Importing here to avoid it as a hard dependency for running tests
    import asdf
    file_path = get_idl_test_output_filepath(file_name, version)
    if not file_path.is_file():
        if idl_env is None:
            pytest.skip("""To run the IDL comparison tests, you must:
                            1. Specify a path to a working IDL executable,
                            2. Specify a path to the CHIANTI IDL code,
                            3. Install the hissw package (pip install hissw)
                           Without the following, you will not be able to generate new IDL comparison test results.""")
        result = idl_env.run(script, args=input_args, save_vars=save_vars)
        if format_func is not None:
            for k in format_func:
                result[k] = format_func[k](result[k])
        # Save the the contribution function to the data directory
        variables = {**result, **input_args, 'idl_script': script}
        if write_file:
            with asdf.AsdfFile(variables) as af:
                af.write_to(file_path)
        else:
            return variables
    return read_idl_test_output(file_name, version)
