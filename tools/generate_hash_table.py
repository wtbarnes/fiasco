"""
Script for generating hash table for all files in the CHIANTI database
"""
import click
import json
import pathlib

from astropy.utils.data import get_pkg_data_path
from itertools import chain

from fiasco.util.setup_db import _md5hash
from fiasco.util.util import get_chianti_catalog, read_chianti_version


def build_ion_path(name):
    id = name.split('.')[0].split('_')
    return pathlib.Path(id[0]) / f'{id[0]}_{id[1]}' / name


def build_hash_table(dbase_root):
    dbase_root = pathlib.Path(dbase_root)
    catalogue = get_chianti_catalog(dbase_root)
    filepaths = chain(
        map(build_ion_path, catalogue['ion_files']),
        map(lambda x: pathlib.Path('abundance') / x, catalogue['abundance_files']),
        map(lambda x: pathlib.Path('ioneq') / x, catalogue['ioneq_files']),
        map(lambda x: pathlib.Path('ip') / x, catalogue['ip_files']),
        map(lambda x: pathlib.Path('continuum') / x, catalogue['continuum_files']),
        map(lambda x: pathlib.Path('dem') / x, catalogue['dem_files']),
    )
    filepaths = map(lambda x: dbase_root / x, filepaths)
    return {'_'.join(f.relative_to(dbase_root).parts): _md5hash(f) for f in filepaths}


@click.command()
@click.option('--database', required=True, type=str)
def write_hash_table(database):
    data_dir = pathlib.Path(get_pkg_data_path('data', package='fiasco.tests'))
    dbase_version = read_chianti_version(database)
    file_path = data_dir / f'file_hashes_v{dbase_version}.json'
    hash_table = build_hash_table(database)
    with open(file_path, 'w') as f:
        json.dump(hash_table, f)


if __name__ == '__main__':
    write_hash_table()
