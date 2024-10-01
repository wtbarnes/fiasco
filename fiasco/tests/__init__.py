# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module contains package tests.
"""
import json
import pathlib

from astropy.utils.data import get_pkg_data_path

__all__ = ['get_test_file_list']


def get_test_file_list():
    data_dir = pathlib.Path(get_pkg_data_path('data', package='fiasco.tests'))
    file_path = data_dir / 'test_file_list.json'
    with open(file_path) as f:
        hash_table = json.load(f)
    return hash_table['test_files']
