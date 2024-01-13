"""
Helpers for tests related to comparing with IDL output
"""
import pathlib

from astropy.utils.data import get_pkg_data_path

from fiasco.util import read_chianti_version

__all__ = ['save_idl_test_output']


def save_idl_test_output(variables, name, ascii_dbase_root):
    # NOTE: Importing here to avoid it as a hard dependency for running tests
    import asdf
    version = read_chianti_version(ascii_dbase_root)
    fname = f"{name}_v{version['major']}.{version['minor']}.{version['patch']}.asdf"
    data_dir = pathlib.Path(get_pkg_data_path('data', package='fiasco.tests.idl'))
    file_path =  data_dir / fname
    with asdf.AsdfFile(variables) as af:
        af.write_to(file_path)
