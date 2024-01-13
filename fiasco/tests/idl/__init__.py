"""
Tests that directly compare the output of fiasco and CHIANTI IDL
"""
import asdf
import pathlib

from astropy.utils.data import get_pkg_data_path

from fiasco.util import read_chianti_version


def save_idl_test_output(variables, name, ascii_dbase_root):
    version = read_chianti_version(ascii_dbase_root)
    fname = f"{name}_v{version['major']}.{version['minor']}.{version['patch']}.asdf"
    data_dir = pathlib.Path(get_pkg_data_path('data', package='fiasco.tests.idl'))
    file_path =  data_dir / fname
    with asdf.AsdfFile(variables) as af:
        af.write_to(file_path)
