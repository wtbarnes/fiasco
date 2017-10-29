# fiasco
[![Powered by SunPy Badge]( http://img.shields.io/badge/powered%20by-SunPy-orange.svg?style=flat)](http://www.sunpy.org)

A Python interface to the [CHIANTI atomic database](http://www.chiantidatabase.org/)

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
See [this notebook](https://gist.github.com/wtbarnes/1a9acc43381a33d016843c666150e8e9)

## Why *fiasco*?
A *fiasco*, or flask, is [the typical style of bottle](https://en.wikipedia.org/wiki/Fiasco_(bottle)) used to serve the *Chianti Classico* wine. It is typically larger and rounder at the bottom and is covered by a straw basket. In the same way, the `fiasco` package serves up the CHIANTI atomic database.

## Develop
See the [wiki](https://github.com/wtbarnes/fiasco/wiki) for notes on the development of this package.

## License
This project is Copyright (c) Will Barnes and licensed under the terms of the BSD 3-Clause license. See the licenses folder for more information.
