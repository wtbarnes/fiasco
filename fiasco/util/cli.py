"""
CLI entrypoints for various tasks.
"""

import astropy.units as u
import click
import matplotlib.pyplot as plt
import numpy as np
import pathlib

import fiasco

from fiasco.util import check_database
from fiasco.util.setup_db import _get_chianti_dbase_url


@click.command()
@click.argument(
    'output',
    type=click.Path(path_type=pathlib.Path, dir_okay=False, file_okay=True),
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
        raise click.UsageError(f'ASCII database already exists at {ascii_dbase_root}. '
                               'Specify URL, version, or existing path.')
    kwargs = {
        'ask_before': False,
        'check_chianti_version': False,
        'show_progress': True,
        'ascii_dbase_url': url,
        'ascii_dbase_root': ascii_dbase_root,
    }
    check_database(hdf5_dbase_root, **kwargs)


@click.command()
@click.argument(
    'element',
    type=str,
)
@click.option(
    '--rates',
    'calculate_rates',
    is_flag=True,
    help=('Calculate ionization fraction using the ionization and recombination rates. '
          'Default is to use the precomputed ionization fractions in the database.')
)
@click.option(
    '--temperature',
    type=(float, float, int),
    default=(4,9,100),
    help='Temperature minimum, maximum, and number of points in log space.'
)
@click.option(
    '--ion',
    'ion_names',
    multiple=True,
    type=str,
    help='Plot only certain ions. By default, all ions of the given element are plotted.',
)
@click.option(
    '--filename',
    type=click.Path(path_type=pathlib.Path, file_okay=True, dir_okay=False),
    default=None,
    help='Path to save plot to. By default, plot appears in an interactive window.'
)
@click.option(
    '--database',
    type=click.Path(path_type=pathlib.Path, file_okay=True, dir_okay=False),
    default=None,
    help='Path to HDF5 CHIANTI database. If not specified, use default database.'
)
def plot_ionization_equilibrium(element, calculate_rates, temperature, ion_names, filename, database):
    """Make a quick-look plot of the ionization fraction in equilibrium as a function of temperature."""
    temperature = np.logspace(*temperature) * u.K
    el = fiasco.Element(element, temperature, hdf5_dbase_root=database)
    if not ion_names:
        ion_names = [ion.ion_name for ion in el]
    fig = plt.figure(figsize=(10,4), layout='constrained')
    ax = fig.add_subplot()
    for name in ion_names:
        ion = el[name]
        if calculate_rates:
            ionization_fraction = el.equilibrium_ionization[:, ion.charge_state]
        else:
            ionization_fraction = ion.ionization_fraction
        imax = np.argmax(ionization_fraction)
        ax.plot(el.temperature, ionization_fraction)
        ax.text(el.temperature[imax].to_value('K'),
                ionization_fraction[imax].to_value(),
                ion.ionization_stage_roman,
                verticalalignment='bottom',
                horizontalalignment='center')
    ax.set_xscale('log')
    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Ionization Fraction')
    ax.set_title(f'{el.element_name.capitalize()} Ionization Fractions in Equilibrium')
    ax.set_xlim(el.temperature[[0,-1]].to_value('K'))
    if filename is None:
        plt.show()
    else:
        fig.savefig(filename)
