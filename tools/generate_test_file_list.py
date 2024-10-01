"""
This script is for generating the file that contains the list of files that go into version of the database used for testing
"""
import json
import pathlib

from astropy.utils.data import get_pkg_data_path

__all__ = ['write_test_file_list']


def write_test_file_list(files):
    """
    Write a unique, sorted list of files to use in the test database.
    These are written to fiasco/util/data/test_files.json
    """
    # Remove duplicates
    files_unique = list(set(files))

    # Sort by element name and ionization stage
    def sort_func(x):
        components = x.split('.')
        # Case 0: There are multiple '.' characters in the filename
        if len(components) > 2:
            return ('.'.join(components[:-1]), 0, components[-1])
        name, ext = components
        components = name.split('_')
        # Case 1: Filenames are not ion files
        if len(components) != 2:
            return (name, 0, ext)
        element, ion = components
        # Dielectronic files have a "_d" appended to them
        if ion[-1] == 'd':
            ion = ion[:-1]
            element += '_d'
        try:
            ion = int(ion)
        except ValueError:
            # Case 2: Filenames have only 1 underscore and are not ion files
            return (name, 0, ext)
        else:
            # Case 3: Filenames are ion files
            return (element, ion, ext)

    files_sorted = sorted(files_unique, key=sort_func)

    # Write to JSON
    data_dir = pathlib.Path(get_pkg_data_path('data', package='fiasco.util'))
    with open(data_dir / 'test_file_list.json', mode='w') as f:
        json.dump({'test_files': files_sorted}, f, indent=2)


if __name__ == '__main__':
    # An example of how you might use this function to update the test file list
    from fiasco.tests import get_test_file_list

    test_files = get_test_file_list()  # Read current files
    test_files += ...  # Add new files here
    write_test_file_list(test_files)
