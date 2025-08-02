.. _fiasco-how-to-guide:

How-To Guides
=============

These how-to guides provide examples of how to perform specific tasks with fiasco.

.. _fiasco-how-to-download-chianti:

Download the CHIANTI Database
-----------------------------

The CHIANTI data is distributed by the CHIANTI team as a collection of ASCII text files in a series of subdirectories. Rather than interact with these raw text files directly, fiasco first builds the entire CHIANTI database into a single `HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ file to allow for easier and faster access to the data.

There are two ways of downloading and setting up the database to be used by fiasco,

1. Allow fiasco to download the latest release of the database for you **(recommended)**
2. Use an existing install of the CHIANTI database. This may be the best option for users who have already installed the CHIANTI package from `SolarSoft (SSW) <http://www.lmsal.com/solarsoft/>`_.

**If you choose option 1, no action is required.
This is the recommended option, particularly for users new to CHIANTI.**
For more information on this approach, see :ref:`fiasco-quick-start`.

Option 2 requires you to tell fiasco where to find an existing version of the CHIANTI atomic data.
You can do this by setting the path to the CHIANTI data in :file:`$HOME/.fiasco/fiascorc`.
For users who have installed the CHIANTI package with SSW, this file might look like:

.. code-block:: ini

   [database]
   ascii_dbase_root=/path/to/ssw/packages/chianti/dbase

Those who are using a standalone version of the CHIANTI database can just provide the path directly to the (untarred) directory.
For users who are familiar with the CHIANTI IDL tools, ``ascii_dbase_root`` should point to the same directory as ``XUVTOP`` system variable.

.. _fiasco-how-to-hdf5-location:

Build the HDF5 Database in a Different Location
-----------------------------------------------

Additionally, you may also choose to place the HDF5 database in a location besides :file:`$HOME/.fiasco/`.
For example, a user who wants to place the HDF5 file in a custom folder would place the following in their :file:`$HOME/.fiasco/fiascorc` file:

.. code-block:: ini

   [database]
   hdf5_dbase_root=/path/to/chianti/hdf5/file/chianti.h5

.. _fiasco-how-to-rebuild-hdf5:

Re-build the HDF5 Database
--------------------------

There are occasional changes to fiasco that require the HDF5 database to be rebuilt because of changes to the data structure, such as when there are
updates to the data parsing code, or when additional data files have been added to the database.  In these cases, the HDF5 database can be
rebuilt by running the following code in python:

.. code-block:: pycon

   >>> import fiasco
   >>> fiasco.util.build_hdf5_dbase(fiasco.defaults['ascii_dbase_root'],
   ...                              fiasco.defaults['hdf5_dbase_root'],
   ...                              overwrite=True)  # doctest: +SKIP

which will overwrite the old HDF5 file (if any) with a new version.  The arguments of the function can alternatively point to your preferred location(s),
rather than just using the default values.
