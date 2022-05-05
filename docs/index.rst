fiasco Documentation
=====================

.. figure:: _static/fiasco-logo.png
   :figwidth: 50 %
   :align: center

Welcome to the documentation for fiasco. fiasco provides a Python interface to the `CHIANTI atomic database`_. In addition to several high-level abstractions of the atomic data, fiasco also provides many common atomic physics calculations. fiasco takes much of its inspiration from the `ChiantiPy`_ package which provides similar functionality.

.. warning:: fiasco is still in the early stages of development and as such frequent changes will be
             made to the API. Use at your own risk.

.. warning:: fiasco currently only supports v8 of the CHIANTI database. In
             the future, older versions of the database will be supported.

.. toctree::
  :maxdepth: 1

  install
  generated/gallery/index
  chianti
  code_ref/index
  references
  whatsnew/index

.. _CHIANTI atomic database: http://www.chiantidatabase.org/
.. _ChiantiPy: https://github.com/chianti-atomic/ChiantiPy
