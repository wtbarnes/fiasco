Getting Started
================

Dependencies
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


Installation
------------

Finally, install fiasco from GitHub, either via pip

.. code-block:: shell

   $ pip install git+git://github.com/wtbarnes/fiasco.git

or from source using `setup.py`,

.. code-block:: shell

   $ git clone https://github.com/wtbarnes/fiasco.git
   $ cd fiasco
   $ python setup.py install

.. note:: Currently, fiasco is only available via GitHub. We plan to make it available as a package
          on both `conda forge`_ and pip in the near future.

Acquiring the Atomic Data
-------------------------
The CHIANTI data is distributed by the CHIANTI team as a collection of ASCII text files in a series of subdirectories. Rather than interact with these raw text files directly, fiasco first builds the entire CHIANTI database into a single `HDF5`_ file to allow for easier and faster access to the data.

There are two ways of downloading and setting up the database to be used by fiasco,

1. Allow fiasco to download the latest release of the database for you **(recommended)**
2. Use an existing install of the CHIANTI database. This may be the best option for users who have already installed the CHIANTI package from `SSW`_.

**If you choose option 1, no action is required. This is the recommended option, particularly for users new to CHIANTI.** The first time that you instantiate a subclass of `~fiasco.IonBase`, you will be prompted to download the raw CHIANTI data and then build the HDF5 database from it. Both the raw data and the HDF5 data will be placed in `$HOME/.fiasco/`.

Option 2 requires you to tell fiasco where to find an existing version of the CHIANTI atomic data. You can do this by setting the path to the CHIANTI data in `$HOME/.fiasco/fiascorc`. For users who have installed the CHIANTI package with SSW, this file might look like,

.. code-block:: ini

   [database]
   ascii_dbase_root=/path/to/ssw/packages/chianti/dbase

Those who have manually installed the CHIANTI database can just provide the path directly to the (untarred) directory.

Additionally, you may also choose to place the HDF5 database in a location besides `$HOME/.fiasco/`. For example, a user who wants to use a manually installed version of the CHIANTI database and place the HDF5 file in a custom folder would place the following in their `$HOME/.fiasco/fiascorc` file,

.. code-block:: ini

   [database]
   ascii_dbase_root=/path/to/custom/untarred/chianti
   hdf5_dbase_root=/path/to/chianti/hdf5/file/chianti.h5

To double check that these presets have been set correctly,

.. code-block:: python

   >>> import fiasco
   >>> fiasco.defaults # doctest: +SKIP

should show the paths set in the configuration file.

.. _Anaconda distribution: https://docs.anaconda.com/anaconda/install/
.. _this page: http://docs.sunpy.org/en/stable/guide/installation/index.html#installing-scientific-python-and-sunpy
.. _conda forge: https://conda-forge.org/
.. _SSW: http://www.lmsal.com/solarsoft/
.. _HDF5: https://en.wikipedia.org/wiki/Hierarchical_Data_Format
