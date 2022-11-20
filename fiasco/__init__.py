"""
A Python interface to the CHIANTI atomic database
"""
from fiasco.collections import *
from fiasco.element import *
from fiasco.fiasco import *
from fiasco.ion import *
from fiasco.level import *
from fiasco.util.util import setup_paths

try:
    from fiasco.version import __version__
except ImportError:
    __version__ = "unknown"

defaults = setup_paths()


from fiasco.util.logger import _init_log

log = _init_log()
