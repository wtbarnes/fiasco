"""
fiasco
======

a prototype package for interfacing with the CHIANTI atomic database
"""
from .base import IonBase, ElementBase

from .util import setup_paths
defaults = setup_paths()
