.. _fiasco-how-to-guide:

How-To Guides
=============

These how-to guides provide examples of how to perform specific tasks with fiasco.

.. _fiasco-how-to-dev-install:

Install the Development Version of fiasco
-----------------------------------------

To install the development release of the package from GitHub, run:

.. code-block:: shell

   $ git clone https://github.com/wtbarnes/fiasco.git
   $ cd fiasco
   $ pip install -e .

If you also want the dependencies needed to run the tests and build the documentation, run:

.. code-block:: shell

   $ pip install -e .[dev]

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

Re-building the HDF5 Database
-----------------------------------------------

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

.. _fiasco-how-to-run-tests:

Run the Test Suite
------------------

To execute the test suite, run:

.. code-block:: shell

   $ pytest fiasco

This will download a copy of CHIANTI to a temporary directory, build a temporary version of the HDF5 database, run the tests, and then delete all of these files once the tests have been run.
When running the test suite many times locally, it is often preferable to use a copy of the database you have already downloaded to avoid repeated downloads.
To do this, run:

.. code-block:: shell

   $ pytest fiasco --ascii-dbase-root /path/to/chianti/dbase

By default, the test suite will only include a minimal set of files in the built database for testing to reduce the time needed to build the database.
To instead include all files, you can pass the following flag:

.. code-block:: shell

   $ pytest fiasco --ascii-dbase-root /path/to/chianti/dbase --include-all-files


If you would also like to avoid rebuilding the HDF5 database each time, you can pass the following flag:

.. code-block:: shell

   $ pytest fiasco --ascii-dbase-root /path/to/chianti/dbase --hdf5-dbase-root /path/to/chianti/chianti.h5

.. _fiasco-how-to-run-tests-idl:

Testing Against IDL Routines
----------------------------

The `fiasco` test suite includes a set of tests that automatically compare against the equivalent routines in the
CHIANTI IDL software.
The purpose of these tests is to provide a systematic way to assess any deviations from the original IDL software.
By default, these tests are run by comparing a set of stored results in :file:`fiasco/tests/idl/data` since it is not possible to run the accompanying IDL code in a continuous integration environment due to the licensing restrictions imposed by IDL.
However, it may be necessary to instead run the IDL commands locally in order to directly compare results, e.g. if you want to compare against a new version of CHIANTI for which there are not yet cached results.

To run these tests and also run the accompanying IDL code, first install the needed dependencies:

.. code-block:: shell

   $ pip install -e .[test-idl]

and then run pytest with the following flags:

.. code-block:: shell

   $ pytest fiasco/tests/idl/ \
            --idl-executable=/path/to/idl \
            --idl-codebase-root=/path/to/chianti/idl \
            --ascii-dbase-root=/path/to/chianti/dbase \
            --include-all-files

where :file:`/path/to/idl/` is the path to the directory containing :file:`bin/idl` (where ``idl`` is the IDL executable),
:file:`/path/to/chianti/idl` is the path to the directory containing all of the CHIANTI IDL routines,
and :file:`/path/to/chianti/dbase` is the path to the top directory of the CHIANTI atomic database (what you would usually
set as ``!XUVTOP`` in CHIANTI IDL).
Note that the IDL and database files should be from the same version.

This command will run the equivalent IDL commands inside of an isolated IDL environment using only those CHIANTI files.
Installing SSW is not required.
Note that these tests will also generate a new set of cached results in :file:`fiasco/tests/idl/data`.
If you're running these tests for a version of the database for which there are already cached results, the IDL code will not be executed and the cached results will be used instead.
To force the IDL code to run, you can delete the cached result files and they will be regenerated the next time you run the test.
