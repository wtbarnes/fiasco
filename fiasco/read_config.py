"""
Parse .rc file
"""

import os
try:
    # Python 3
    import configparser
except ImportError:
    # Python 2
    import ConfigParser as configparser

__all__ = ['defaults']

# parse config file
defaults = {}
if os.path.isfile(os.path.join(os.environ['HOME'], '.fiasco', 'fiascorc')):
    config = configparser.ConfigParser()
    config.read(os.path.join(os.environ['HOME'], '.fiasco', 'fiascorc'))
    if 'database' in config:
        if 'dbase_root' in config['database']:
            defaults['chianti_dbase_root'] = config['database']['dbase_root']
        if 'hdf5_dbase_root' in config['database']:
            defaults['chianti_hdf5_dbase_root'] = config['database']['hdf5_dbase_root']
    
# set defaults
if 'chianti_dbase_root' not in defaults:
    try:
        defaults['chianti_dbase_root'] = os.environ['XUVTOP']
    except KeyError:
        defaults['chianti_dbase_root'] = ''
if 'chianti_hdf5_dbase_root' not in defaults:
    defaults['chianti_hdf5_dbase_root'] = os.path.join(os.environ['HOME'],'.fiasco','chianti_dbase.h5')
