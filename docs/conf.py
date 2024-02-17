#
# Configuration file for the Sphinx documentation builder.
#
import configparser
import os
import pathlib

# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config


# -- Project information -----------------------------------------------------

project = 'fiasco'
copyright = '2022, Will Barnes'
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
    'sphinxcontrib.bibtex',
    'sphinx_design',
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
    'numpy': ('https://numpy.org/doc/stable/',
              (None, 'http://data.astropy.org/intersphinx/numpy.inv')),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/',
              (None, 'http://data.astropy.org/intersphinx/scipy.inv')),
    'matplotlib': ('https://matplotlib.org/',
                   (None, 'http://data.astropy.org/intersphinx/matplotlib.inv')),
    'astropy': ('https://docs.astropy.org/en/stable', None),
    'sunpy': ('https://docs.sunpy.org/en/stable/', None),
    'aiapy': ('https://aiapy.readthedocs.io/en/stable/', None),
    'plasmapy': ('https://docs.plasmapy.org/en/stable', None),
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'pydata_sphinx_theme'

# -- Sphinx Book Theme Options -----------------------------------------------------
html_theme_options = {
    "use_edit_page_button": True,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/wtbarnes/fiasco",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/fiasco/",
            "icon": "fa-brands fa-python",
        },
        {
            "name": "CHIANTI",
            "url": "http://chiantidatabase.org/",
            "icon": "fa-solid fa-wine-glass",
        }
    ],
    #"secondary_sidebar_items": {
    #    "index": []  # Remove secondary sidebar on landing page
    #},
    "announcement": "fiasco currently only supports version 8 of the CHIANTI database.",
}
html_context = {
    "github_user": "wtbarnes",
    "github_repo": "fiasco",
    "github_version": "main",
    "doc_path": "docs",
}
html_logo = '_static/fiasco-logo.png'


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

# Set path for BibTeX file for all of our references
bibtex_bibfiles = ['references.bib']

ON_RTD = os.environ.get('READTHEDOCS') == 'True'
ON_GHA = os.environ.get('CI') == 'true'

# On Read the Docs and CI, download the database and build a minimal HDF5 version
if ON_RTD or ON_GHA:
    from fiasco.util import build_hdf5_dbase, download_dbase, get_test_file_list
    from fiasco.util.setup_db import CHIANTI_URL, LATEST_VERSION
    FIASCO_HOME = pathlib.Path.home() / '.fiasco'
    FIASCO_HOME.mkdir(exist_ok=True, parents=True)
    ascii_dbase_root = FIASCO_HOME / 'chianti_dbase'
    hdf5_dbase_root = FIASCO_HOME / 'chianti_dbase.h5'
    download_dbase(CHIANTI_URL.format(version=LATEST_VERSION), ascii_dbase_root)
    build_hdf5_dbase(
        ascii_dbase_root,
        hdf5_dbase_root,
        files=get_test_file_list(),
    )
    with (FIASCO_HOME / 'fiascorc').open('w') as f:
        c = configparser.ConfigParser()
        c.add_section('database')
        c.set('database', 'ascii_dbase_root', str(ascii_dbase_root))
        c.set('database', 'hdf5_dbase_root', str(hdf5_dbase_root))
        c.write(f)

# -- Sphinx gallery -----------------------------------------------------------
from sphinx_gallery.sorting import ExampleTitleSortKey, ExplicitOrder  # NOQA: E402

extensions += ['sphinx_gallery.gen_gallery']
sphinx_gallery_conf = {
    'subsection_order': ExplicitOrder([
        '../examples/user_guide/',
        '../examples/idl_comparisons/',
    ]),
    'within_subsection_order': ExampleTitleSortKey,
    'examples_dirs': '../examples',   # path to your example scripts
    'gallery_dirs': 'generated/gallery',  # path to where to save gallery generated output
    'filename_pattern': '^((?!skip_).)*$',
    'default_thumb_file': '_static/fiasco-logo.png'
}
