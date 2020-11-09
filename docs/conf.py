# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
import os
import configparser

# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config


# -- Project information -----------------------------------------------------

project = 'fiasco'
copyright = '2020, Will Barnes'
author = 'Will Barnes'

# The full version, including alpha/beta/rc tags
from fiasco import __version__
release = __version__
is_development = '.dev' in __version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'matplotlib.sphinxext.plot_directive',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.inheritance_diagram',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'sphinx.ext.mathjax',
    'sphinx_automodapi.automodapi',
    'sphinx_automodapi.smart_resolver',
    'sunpy.util.sphinx.changelog',
]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The reST default role (used for this markup: `text`) to use for all
# documents. Set to the "smart" one.
default_role = 'obj'

# Disable having a separate return type row
napoleon_use_rtype = False

# Disable google style docstrings
napoleon_google_docstring = False

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/',
               (None, 'http://data.astropy.org/intersphinx/python3.inv')),
    'numpy': ('https://docs.scipy.org/doc/numpy/',
              (None, 'http://data.astropy.org/intersphinx/numpy.inv')),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/',
              (None, 'http://data.astropy.org/intersphinx/scipy.inv')),
    'matplotlib': ('https://matplotlib.org/',
                   (None, 'http://data.astropy.org/intersphinx/matplotlib.inv')),
    'astropy': ('http://docs.astropy.org/en/stable/', None),
    'sunpy': ('https://docs.sunpy.org/en/stable/', None)}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

try:
    html_theme = "sphinx_rtd_theme"
    import sphinx_rtd_theme
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
except ImportError:
    html_theme = 'default'


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# Render inheritance diagrams in SVG
graphviz_output_format = "svg"

graphviz_dot_args = [
    '-Nfontsize=10',
    '-Nfontname=Helvetica Neue, Helvetica, Arial, sans-serif',
    '-Efontsize=10',
    '-Efontname=Helvetica Neue, Helvetica, Arial, sans-serif',
    '-Gfontsize=10',
    '-Gfontname=Helvetica Neue, Helvetica, Arial, sans-serif'
]


ON_RTD = os.environ.get('READTHEDOCS').lower() == 'true'
ON_GHA = os.environ.get('CI').lower() == 'true'

# On Read the Docs and CI, download the database and build a minimal HDF5 version
if ON_RTD or ON_GHA:
    from fiasco.util import download_dbase, build_hdf5_dbase
    from fiasco.util.setup_db import CHIANTI_URL, LATEST_VERSION
    if ON_RTD:
        os.environ['HOME'] = '/home/docs'  # RTD does not set HOME?
    FIASCO_HOME = os.path.join(os.environ['HOME'], '.fiasco')
    if not os.path.exists(FIASCO_HOME):
        os.makedirs(FIASCO_HOME)
    ascii_dbase_root = os.path.join(FIASCO_HOME, 'ascii_dbase')
    hdf5_dbase_root = os.path.join(FIASCO_HOME, 'chianti_dbase.h5')
    download_dbase(CHIANTI_URL.format(version=LATEST_VERSION), ascii_dbase_root)
    build_hdf5_dbase(
        ascii_dbase_root,
        hdf5_dbase_root,
        files=['chianti.ip',
               'chianti.ioneq',
               'sun_photospheric_1998_grevesse.abund',
               'fe_5.elvlc',
               'fe_15.elvlc',
               'fe_18.rrparams',
               'fe_18.drparams',
               'o_6.scups',
               'o_6.elvlc',
               'o_6.wgfa',]
    )
    with open(os.path.join(FIASCO_HOME, 'fiascorc'), 'w') as f:
        c = configparser.ConfigParser()
        c.add_section('database')
        c.set('database', 'ascii_dbase_root', ascii_dbase_root)
        c.set('database', 'hdf5_dbase_root', hdf5_dbase_root)
        c.write(f)

# -- Sphinx gallery -----------------------------------------------------------
extensions += ['sphinx_gallery.gen_gallery']
sphinx_gallery_conf = {
     'examples_dirs': '../examples',   # path to your example scripts
     'gallery_dirs': 'auto_examples',  # path to where to save gallery generated output
}
