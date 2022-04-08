"""
A Python interface to the CHIANTI atomic database
"""
from .util.util import setup_paths
from .fiasco import *
from .level import *
from .ion import *
from .collections import *
from .element import *

try:
    from .version import __version__
except ImportError:
    __version__ = "unknown"

defaults = setup_paths()


from fiasco.util.logger import _init_log
log = _init_log()
