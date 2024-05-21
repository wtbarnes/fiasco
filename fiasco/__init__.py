"""
fiasco: A Python interface to the CHIANTI atomic database
"""
from fiasco.collections import IonCollection
from fiasco.elements import Element
from fiasco.fiasco import list_elements, list_ions, proton_electron_ratio
from fiasco.ions import Ion
from fiasco.levels import Level, Transitions
from fiasco.util.util import setup_paths

try:
    from fiasco.version import __version__
except ImportError:
    __version__ = "unknown"

defaults = setup_paths()


from fiasco.util.logger import _init_log

log = _init_log()

__all__ = ["IonCollection", "Element", "list_elements", "list_ions", "proton_electron_ratio", "Ion", "Level", "Transitions", "defaults", "log", "__version__"]
