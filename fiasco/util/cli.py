"""CLI entrypoints for building the CHIANTI database."""

import click
import pathlib

from fiasco.util import check_database
from fiasco.util.setup_db import _get_chianti_dbase_url


@click.command()
@click.argument(
    'output',
    type=click.Path(path_type=pathlib.Path, dir_okay=False, file_okay=True),
    help='Path to HDF5 file that will store the CHIANTI database.'
)
@click.option(
    '--ascii-dbase',
    'ascii_dbase',
    type=click.Path(path_type=pathlib.Path, file_okay=False, dir_okay=True),
    default=None,
    help=('Path to the ASCII CHIANTI database tree. If not specified, a new copy is downloaded to '
          'a new directory dbase in the same directory as the HDF5 database.'),
)
@click.option(
    '--url',
    type=str,
    default=None,
    help='URL to download the ASCII CHIANTI database from.',
)
@click.option(
    '--version',
    type=str,
    default=None,
    help='CHIANTI database version to use when downloading the ASCII database.',
)
def build_and_download_database(output, ascii_dbase, url, version):
    """Build the HDF5 CHIANTI database from an ASCII database tree."""
    if url is not None and version is not None:
        raise click.UsageError('Specify either --url or --version, not both.')
    hdf5_dbase_root = pathlib.Path(output)
    if hdf5_dbase_root.exists():
        click.echo(f'HDF5 database already exists at {hdf5_dbase_root}')
        return
    if url is None:
        url = _get_chianti_dbase_url(version=version)
    ascii_dbase_root = ascii_dbase
    if ascii_dbase_root is None:
        ascii_dbase_root = hdf5_dbase_root.parent / 'dbase'
    if (url is not None or version is not None) and ascii_dbase_root.exists():
        raise click.UsageError(f'ASCII database already exists at {ascii_dbase_root}. Specify URL, version, or existing path.')
    kwargs = {
        'ask_before': False,
        'check_chianti_version': False,
        'show_progress': True,
        'ascii_dbase_url': url,
        'ascii_dbase_root': ascii_dbase_root,
    }
    check_database(hdf5_dbase_root, **kwargs)
