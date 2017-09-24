"""
fiasco
======

a prototype package for interfacing with the CHIANTI atomic database
"""
from .base import IonBase
from fiasco.util import setup_paths, download_dbase
defaults = setup_paths()
download_dbase(defaults['ascii_dbase_root'])