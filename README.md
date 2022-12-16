# fiasco

[![Powered by SunPy Badge]( http://img.shields.io/badge/powered%20by-SunPy-orange.svg?style=flat)](http://www.sunpy.org)
[![fiasco CI status](https://github.com/wtbarnes/fiasco/workflows/CI/badge.svg?branch=main)](https://github.com/wtbarnes/fiasco/actions)
[![Documentation Status](https://readthedocs.org/projects/fiasco/badge/?version=latest)](http://fiasco.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/wtbarnes/fiasco/branch/main/graph/badge.svg?token=damCmTyRUN)](https://codecov.io/gh/wtbarnes/fiasco)
[![matrix](https://img.shields.io/matrix/atomic-data:openastronomy.org.svg?colorB=%23FE7900&label=Chat&logo=matrix&server_fqdn=openastronomy.modular.im)](https://openastronomy.element.io/#/room/#atomic-data:openastronomy.org)

A Python interface to the [CHIANTI atomic database](http://www.chiantidatabase.org/). For a high level
overview of the package, have a look at my talk ([slides](https://zenodo.org/record/1249002), [video](https://youtu.be/7_Nr700kBME)) from
the [2018 Python in Astronomy](http://openastronomy.org/pyastro/2018/) conference.

## Install

```shell
$ pip install fiasco
```

The [CHIANTI atomic database](http://www.chiantidatabase.org/chianti_download.html) is also required.
See [the docs](https://fiasco.readthedocs.io/en/latest/chianti.html#acquiring-the-atomic-data) for more details.

## Example

```python
>>> import astropy.units as u
>>> import fiasco
>>> iron = fiasco.Element('iron', [1e4, 1e6, 1e8]*u.K)
# Print some information about the element
>>> iron.atomic_number
26
>>> iron.atomic_symbol
'Fe'
>>> iron.abundance
<Quantity 3.16227766e-05>
# Select the Fe 16 ion
>>> iron[15].ion_name
'Fe 16'
>>> iron[15].charge_state
15
# Ionization fraction
>>> iron[15].ioneq
<Quantity [0.000e+00, 2.377e-08, 4.163e-18]>
```

## Why *fiasco*?

A *fiasco*, or flask, is [the typical style of bottle](https://en.wikipedia.org/wiki/Fiasco_(bottle)) used to serve the *Chianti Classico* wine. It is typically larger and rounder at the bottom and is covered by a straw basket. In the same way, the `fiasco` package serves up the CHIANTI atomic database.
