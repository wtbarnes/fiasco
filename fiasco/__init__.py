"""
A Python interface to the CHIANTI atomic database
"""

# Licensed under a 3-clause BSD style license - see LICENSE.rst
# ----------------------------------------------------------------------------
from ._sunpy_init import *
# ----------------------------------------------------------------------------

if not _ASTROPY_SETUP_:
    from .util.util import setup_paths
    defaults = setup_paths()

    from .fiasco import *
    from .level import *
    from .ion import *
    from .collections import *
    from .element import *
