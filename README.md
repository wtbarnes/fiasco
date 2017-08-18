# fiasco
Prototype for Python interface to the [CHIANTI atomic database](http://www.chiantidatabase.org/)

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
* fortranformat (pip only)
* periodictable (pip only)

The [CHIANTI atomic database](http://www.chiantidatabase.org/chianti_download.html) is also required.

## Install
```shell
$ git clone https://github.com/wtbarnes/fiasco.git
$ cd fiasco && python setup.py install
```

## Example

## Why *fiasco*?
A *fiasco*, or flask, is [the typical style of bottle](https://en.wikipedia.org/wiki/Fiasco_(bottle)) used to serve the *Chianti Classico* wine. It is typically larger and rounder at the bottom and is covered by a straw basket. In the same way, the `fiasco` package serves up the CHIANTI atomic database.

## Develop
See the [wiki](https://github.com/wtbarnes/fiasco/wiki) for notes on the development of this package.
