#
# Configuration file for the Sphinx documentation builder.
#
import configparser
import datetime
import os

from packaging.version import Version
from sphinx_gallery.sorting import ExplicitOrder

# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config


# -- Project information -----------------------------------------------------
project = 'fiasco'
author = 'Will Barnes'
copyright = f'{datetime.datetime.utcnow().year}, {author}'

from fiasco import __version__

_version_ = Version(__version__)
# NOTE: Avoid "post" appearing in version string in rendered docs
if _version_.is_postrelease:
    version = release = f'{_version_.major}.{_version_.minor}.{_version_.micro}'
else:
    version = release = str(_version_)
is_development = _version_.is_devrelease

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinxcontrib.bibtex',
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
    'sphinx_design',
    'sphinx_gallery.gen_gallery',
]
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
source_suffix = '.rst'
master_doc = 'index'
default_role = 'obj'
napoleon_use_rtype = False
napoleon_google_docstring = False

# Enable nitpicky mode, which forces links to be non-broken
nitpicky = True
# This is not used. See docs/nitpick-exceptions file for the actual listing.
nitpick_ignore = []
for line in open('nitpick-exceptions'):
    if line.strip() == "" or line.startswith("#"):
        continue
    dtype, target = line.split(None, 1)
    target = target.strip()
    nitpick_ignore.append((dtype, target))

# -- Options for intersphinx extension ---------------------------------------
intersphinx_mapping = {
    "python": (
        "https://docs.python.org/3/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/python3.inv"),
    ),
    "numpy": (
        "https://numpy.org/doc/stable/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/numpy.inv"),
    ),
    "scipy": (
        "https://docs.scipy.org/doc/scipy/reference/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/scipy.inv"),
    ),
    'plasmapy': ('https://docs.plasmapy.org/en/stable', None),
    'sunpy': ('https://docs.sunpy.org/en/stable/', None),
    "aiapy": ("https://aiapy.readthedocs.io/en/stable/", None),
    "asdf": ("https://asdf.readthedocs.io/en/stable/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "packaging": ("https://packaging.pypa.io/en/stable/", None),
}

# -- Options for HTML output -------------------------------------------------
html_theme = 'pydata_sphinx_theme'

# -- Sphinx Book Theme Options -----------------------------------------------
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
    "announcement": "fiasco currently only supports version 8 of the CHIANTI database.",
}
html_context = {
    "github_user": "wtbarnes",
    "github_repo": "fiasco",
    "github_version": "main",
    "doc_path": "docs",
}
html_logo = '_static/fiasco-logo.png'
# Set path for BibTeX file for all of our references
bibtex_bibfiles = ['references.bib']
# Sidebar removal
html_sidebars = {
    "quick_start*": [],
    "how_to_guides*": [],
    "citation*": [],
}
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

# -- Database on RTD -----------------------------------------------
ON_RTD = os.environ.get('READTHEDOCS') == 'True'
ON_GHA = os.environ.get('CI') == 'true'

# On Read the Docs and CI, download the database and build a minimal HDF5 version
if (ON_RTD or ON_GHA):
    from fiasco.tests import get_test_file_list
    from fiasco.util import check_database
    from fiasco.util.setup_db import CHIANTI_URL, LATEST_VERSION
    from fiasco.util.util import FIASCO_HOME, FIASCO_RC
    FIASCO_HOME.mkdir(exist_ok=True, parents=True)
    ascii_dbase_root = FIASCO_HOME / 'chianti_dbase'
    hdf5_dbase_root = FIASCO_HOME / 'chianti_dbase.h5'
    check_database(
        hdf5_dbase_root=hdf5_dbase_root,
        ascii_dbase_root=ascii_dbase_root,
        ascii_dbase_url = CHIANTI_URL.format(version=LATEST_VERSION),
        ask_before=False,
        files=get_test_file_list(),
    )
    with FIASCO_RC.open(mode='w') as f:
        c = configparser.ConfigParser()
        c.add_section('database')
        c.set('database', 'ascii_dbase_root', str(ascii_dbase_root))
        c.set('database', 'hdf5_dbase_root', str(hdf5_dbase_root))
        c.write(f)

# -- Sphinx gallery -----------------------------------------------------------
sphinx_gallery_conf = {
    'backreferences_dir': os.path.join('generated', 'modules'),
    'filename_pattern': '^((?!skip_).)*$',
    'examples_dirs': os.path.join('..', 'examples'),
    'subsection_order': ExplicitOrder([
        '../examples/user_guide/',
        '../examples/idl_comparisons/',
    ]),
    'within_subsection_order': 'ExampleTitleSortKey',
    'gallery_dirs': os.path.join('generated', 'gallery'),
    'matplotlib_animations': True,
    "default_thumb_file": '_static/fiasco-logo.png',
    'abort_on_example_error': False,
    'plot_gallery': 'True',
    'remove_config_comments': True,
    'doc_module': ('fiasco',),
    'only_warn_on_example_error': True,
}
