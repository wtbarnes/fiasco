Getting Started
================

Installation
------------
Before installing fiasco, you'll need to install Python and several other packages. The recommended Python distribution is the `Anaconda distribution`_. fiasco requires Python 3.6+. For additional instructions on installing scientific Python, see `this page`_.

Next, install the following packages,

- numpy
- scipy
- astropy
- h5py
- plasmapy (pip only)
- fortranformat (pip only)

These packages can be installed either using conda (recommended),

.. code-block:: shell

   $ conda install {package-name}

or pip,

.. code-block:: shell

   $ pip install {package-name}

Finally, install fiasco from GitHub,

.. code-block:: shell

   $ git clone https://github.com/wtbarnes/fiasco.git
   $ cd fiasco
   $ python setup.py install

.. note:: Currently, fiasco is only available via GitHub. We plan to make it available on both
          `conda forge`_ and pip in the near future.

Acquiring the Atomic Data
-------------------------


.. _Anaconda distribution: https://docs.anaconda.com/anaconda/install/
.. _this page: http://docs.sunpy.org/en/stable/guide/installation/index.html#installing-scientific-python-and-sunpy
.. _conda forge: https://conda-forge.org/