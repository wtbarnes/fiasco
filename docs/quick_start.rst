.. _fiasco-quick-start:

Quick Start
===========

In this quick start guide, you will learn how to install fiasco, download the CHIANTI atomic database, build an HDF5 version of the database required by fiasco, and create an `~fiasco.Ion` object.
To install fiasco, run the following command in the terminal,

.. code-block:: shell

   $ pip install fiasco

This will install latest version of the package and all needed dependencies from the `Python Package Index <https://pypi.org/project/fiasco/>`_ or PyPI.

Next, let's create an `~fiasco.Ion` object representing Fe XVIII:

.. code-block:: python

   >>> import fiasco
   >>> import astropy.units as u
   >>> fe_18 = fiasco.Ion('Fe XVIII', 1*u.MK)  # doctest: +SKIP

After running this last line, you should see the following message:

.. code-block:: console

   No HDF5 database found at $HOME/.fiasco/chianti_dbase.h5. Build it now? [Y/n]

where ``$HOME`` will be replaced by the path to your home directory.
Press "y" to confirm that you want the built database to have this filename.

.. note::

   If you want to build the database in a different directory, see :ref:`fiasco-how-to-hdf5-location`.

Next, you should see the following message:

.. code-block:: console

   No CHIANTI database found at $HOME/.fiasco/chianti_dbase. Download it from http://download.chiantidatabase.org/CHIANTI_v8.0.7_database.tar.gz? [y/N]

Press "y" to download this version of the CHIANTI database to this directory.
Once you confirm, you should see a progress bar indicating the progress of the download of the tar file containing all the files in the CHIANTI database.
This step should take about a minute or so.

.. note::

   If you want to download the database to a different directory or use an existing copy of the database, see :ref:`fiasco-how-to-download-chianti`.

You should then see a second progress bar indicating the progress of building the HDF5 database.
This step may take a few minutes.
The progress bar output should look something like the following:

.. code-block:: console

   |================================================| 265M/265M (100.00%)      1m 9s
   |=======>----------------------------------------| 496 /3.0k (16.06%) ETA  1m34s

These two steps together may take several minutes, but you will only need to do this once.
fiasco will now look in ``$HOME/.fiasco/chianti_dbase.h5`` for the built HDF5 version of the database.
You can confirm this by looking at the default values read from the configuration file:

.. code-block:: python

   >>> fiasco.defaults  # doctest: +SKIP
   {'ascii_dbase_root': PosixPath('.../.fiasco/chianti_dbase'),
    'hdf5_dbase_root': PosixPath('.../.fiasco/chianti_dbase.h5')}

Now that you have your database, you can use your ion object that you created above to access information about Fe XVIII:

.. code-block:: python

   >>> fe_18  # doctest: +SKIP
   CHIANTI Database Ion
   ---------------------
   Name: Fe 18
   Element: iron (26)
   Charge: +17
   Number of Levels: 337
   Number of Transitions: 7712
   <BLANKLINE>
   Temperature range: [1.000 MK, 1.000 MK]
   <BLANKLINE>
   HDF5 Database: ...chianti_dbase.h5
   Using Datasets:
      ionization_fraction: chianti
      abundance: sun_coronal_1992_feldman_ext
      ip: chianti
      <BLANKLINE>

You are now ready to use fiasco to access and compute derived quantities from the CHIANTI atomic database!
