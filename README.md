# fiasco

[![fiasco CI status](https://github.com/wtbarnes/fiasco/actions/workflows/ci.yml/badge.svg)](https://github.com/wtbarnes/fiasco/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/fiasco/badge/?version=stable)](https://fiasco.readthedocs.io/en/stable/?badge=stable)
[![PyPI](https://img.shields.io/pypi/v/fiasco.svg)](https://pypi.python.org/pypi)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7504785.svg)](https://doi.org/10.5281/zenodo.7504785)
[![codecov](https://codecov.io/gh/wtbarnes/fiasco/branch/main/graph/badge.svg?token=damCmTyRUN)](https://codecov.io/gh/wtbarnes/fiasco)
[![matrix](https://img.shields.io/matrix/atomic-data:openastronomy.org.svg?colorB=%23FE7900&label=Chat&logo=matrix&server_fqdn=openastronomy.modular.im)](https://openastronomy.element.io/#/room/#atomic-data:openastronomy.org)

`fiasco` provides a Python interface to the [CHIANTI](http://www.chiantidatabase.org/), an atomic database used primarily for astrophysical spectroscopy.
In addition to several high-level abstractions of the atomic data, fiasco also provides many common atomic physics calculations.

## Install

```shell
pip install fiasco
```

The [CHIANTI atomic database](http://www.chiantidatabase.org/chianti_download.html) is also required.
See [the docs](https://fiasco.readthedocs.io/en/latest/quick_start.html#fiasco-quick-start) for more details.

## Usage

The primary interface in `fiasco` is the `Ion` object:

```python
>>> import fiasco
>>> import astropy.units as u
>>> fe_18 = fiasco.Ion('Fe XVIII', 1*u.MK)
>>> fe_18
CHIANTI Database Ion
---------------------
Name: Fe 18
Element: iron (26)
Charge: +17
Number of Levels: 337
Number of Transitions: 7712

Temperature range: [1.000 MK, 1.000 MK]

HDF5 Database: ...chianti_dbase.h5
Using Datasets:
    ioneq: chianti
    abundance: sun_coronal_1992_feldman_ext
    ip: chianti
```

For a quick start guide to using `fiasco`, see [this page of the documentation](https://fiasco.readthedocs.io/en/stable/quick_start.html#fiasco-quick-start).
For more advanced examples, see [the example gallery](https://fiasco.readthedocs.io/en/stable/generated/gallery/index.html).

## Acknowledging or Citing fiasco

If you use `fiasco` in any published work, please cite the appropriate version of the software as well as the CHIANTI atomic database.
See [this page](https://fiasco.readthedocs.io/en/stable/citation.html) for additional details.

## Why *fiasco*?

A *fiasco*, or flask, is [the typical style of bottle](https://en.wikipedia.org/wiki/Fiasco_(bottle)) used to serve the *Chianti Classico* wine. It is typically larger and rounder at the bottom and is covered by a straw basket. In the same way, the `fiasco` package serves up the CHIANTI atomic database.
