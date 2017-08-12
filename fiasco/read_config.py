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
    if 'database' in config and 'dbase_root' in config['database']:
        defaults['chianti_dbase_root'] = config['database']['dbase_root']
    else:
        defaults['chianti_dbase_root'] = None

# set defaults
if defaults['chianti_dbase_root'] is None:
    try:
        defaults['chianti_dbase_root'] = os.environ['XUVTOP']
    except KeyError:
        defaults['chianti_dbase_root'] = ''
