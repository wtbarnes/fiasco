"""
Helpers for tests related to comparing with IDL output
"""
import pathlib

from astropy.utils.data import get_pkg_data_path

from fiasco.util import read_chianti_version

__all__ = ['save_idl_test_output', 'get_idl_test_output']


def save_idl_test_output(variables, name, ascii_dbase_root):
    # NOTE: Importing here to avoid it as a hard dependency for running tests
    import asdf
    version = read_chianti_version(ascii_dbase_root)
    data_dir = pathlib.Path(get_pkg_data_path('data', package='fiasco.tests.idl'))
    file_path =  data_dir / f"{name}_v{version}.asdf"
    with asdf.AsdfFile(variables) as af:
        af.write_to(file_path)


def get_idl_test_output(name, version, keys=None):
    # NOTE: Importing here to avoid it as a hard dependency for running tests
    import asdf
    data_dir = pathlib.Path(get_pkg_data_path('data', package='fiasco.tests.idl'))
    file_path = data_dir / f"{name}_v{version}.asdf"
    with asdf.open(file_path, copy_arrays=True) as af:
        if keys is None:
            keys = af.tree.keys()
        output = {k: af.tree[k] for k in keys}
    return output
