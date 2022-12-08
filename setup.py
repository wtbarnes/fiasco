#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pathlib

from itertools import chain
from setuptools import setup
from setuptools.config.setupcfg import read_configuration

################################################################################
# Programmatically generate some extras combos.
################################################################################
extras = read_configuration("setup.cfg")['options']['extras_require']

# Dev is everything
extras['dev'] = list(chain(*extras.values()))

# All is everything but tests and docs
exclude_keys = ("tests", "docs", "dev")
ex_extras = dict(filter(lambda i: i[0] not in exclude_keys, extras.items()))
# Concatenate all the values together for 'all'
extras['all'] = list(chain.from_iterable(ex_extras.values()))

################################################################################
# Version configuration and setup call
################################################################################

VERSION_TEMPLATE = """
# Note that we need to fall back to the hard-coded version if either
# setuptools_scm can't be imported or setuptools_scm can't determine the
# version, so we catch the generic 'Exception'.
try:
    from setuptools_scm import get_version
    __version__ = get_version(root='..', relative_to=__file__)
except Exception:
    __version__ = '{version}'
""".lstrip()

setup(
    extras_require=extras,
    use_scm_version={'write_to': pathlib.Path('fiasco') / 'version.py',
                     'write_to_template': VERSION_TEMPLATE},

)
