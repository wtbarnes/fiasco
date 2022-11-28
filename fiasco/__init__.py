"""
A Python interface to the CHIANTI atomic database
"""
from fiasco.collections import *
from fiasco.elements import *
from fiasco.fiasco import *
from fiasco.ions import *
from fiasco.levels import *
from fiasco.util.util import setup_paths

try:
    from fiasco.version import __version__
except ImportError:
    __version__ = "unknown"

defaults = setup_paths()


from fiasco.util.logger import _init_log

log = _init_log()
