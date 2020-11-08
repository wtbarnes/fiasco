# fiasco
[![Powered by SunPy Badge]( http://img.shields.io/badge/powered%20by-SunPy-orange.svg?style=flat)](http://www.sunpy.org)
[![Build Status](https://travis-ci.org/wtbarnes/fiasco.svg?branch=master)](https://travis-ci.org/wtbarnes/fiasco)
[![Documentation Status](https://readthedocs.org/projects/fiasco/badge/?version=latest)](http://fiasco.readthedocs.io/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/github/wtbarnes/fiasco/badge.svg?branch=master&service=github)](https://coveralls.io/github/wtbarnes/fiasco?branch=master&service=github)

A Python interface to the [CHIANTI atomic database](http://www.chiantidatabase.org/). For a high level
overview of the package, have a look at my talk ([slides](https://zenodo.org/record/1249002), [video](https://youtu.be/7_Nr700kBME)) from
the [2018 Python in Astronomy](http://openastronomy.org/pyastro/2018/) conference.

**DISCLAIMER: fiasco is still in the very early stages of development. As such, the API is changing very frequently and drastically.**

## Dependencies
All packages can be installed with conda (recommended)

```shell
$ conda install {package_name}
```

or pip

```
$ pip install {package name}
```

* NumPy
* scipy
* Astropy
* h5py
* PlasmaPy (pip only, conda soon)
* fortranformat (pip only)

The [CHIANTI atomic database](http://www.chiantidatabase.org/chianti_download.html) is also required.

## Install
```shell
$ git clone https://github.com/wtbarnes/fiasco.git
$ cd fiasco && python setup.py install
```

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

## Develop
See the [wiki](https://github.com/wtbarnes/fiasco/wiki) for notes on the development of this package.

## License
This project is Copyright (c) Will Barnes and licensed under the terms of the BSD 3-Clause license. See the licenses folder for more information.
