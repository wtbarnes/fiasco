The CHIANTI Atomic Database
===========================

Analyzing astrophysical observations requires detailed knowledge of the underlying physics responsible for the formation of the
observed spectra. The CHIANTI atomic database provides a comprehensive set of atomic data for relevant elements, ions, and
transitions in astrophysical plasmas. CHIANTI includes energy levels, transition wavelengths, as well as collision strength
and cross-section data for computing ionization and recombination rates. The papers :ref:`listed below <additional-resources>`
give an extensive description of each version of the database.

Acquiring the Atomic Data
-------------------------

The CHIANTI data is distributed by the CHIANTI team as a collection of ASCII text files in a series of subdirectories. Rather than interact with these raw text files directly, fiasco first builds the entire CHIANTI database into a single `HDF5`_ file to allow for easier and faster access to the data.

There are two ways of downloading and setting up the database to be used by fiasco,

1. Allow fiasco to download the latest release of the database for you **(recommended)**
2. Use an existing install of the CHIANTI database. This may be the best option for users who have already installed the CHIANTI package from `SSW`_.

**If you choose option 1, no action is required.
This is the recommended option, particularly for users new to CHIANTI.**
The first time that you create an instance of `~fiasco.Ion`, you will be prompted to download the raw CHIANTI data and then build the HDF5 database from it.
Both the raw ASCII data and the HDF5 data will be placed in `$HOME/.fiasco/`.
Alternatively, you can also download and set up the database using `~fiasco.util.check_database`.

Option 2 requires you to tell fiasco where to find an existing version of the CHIANTI atomic data. You can do this by setting the path to the CHIANTI data in `$HOME/.fiasco/fiascorc`.
For users who have installed the CHIANTI package with SSW, this file might look like,

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

.. _conda forge: https://conda-forge.org/
.. _SSW: http://www.lmsal.com/solarsoft/
.. _HDF5: https://en.wikipedia.org/wiki/Hierarchical_Data_Format

Testing Against IDL Routines
----------------------------

The `fiasco` test suite includes a set of tests that automatically compare against the equivalent routines in the
CHIANTI IDL software.
The purpose of these tests is to provide a systematic way to assess any deviations from the original IDL software.
To run these tests,

.. code-block:: shell

   $ pytest --idl-executable=/path/to/idl \
            --idl-codebase-root=/path/to/chianti/idl \
            --ascii-dbase-root=/path/to/chianti/dbase \
            --include-all-files \
            fiasco/tests/idl/

where :file:`/path/to/idl/` is the path to the directory containing :file:`bin/idl` (where ``idl`` is the IDL executable),
:file:`/path/to/chianti/idl` is the path to the directory containing all of the CHIANTI IDL routines,
and :file:`/path/to/chianti/dbase` is the path to the top directory of the CHIANTI atomic database (what you would usually
set as ``!XUVTOP`` in CHIANTI IDL).
Note that the IDL and database files should be from the same version.
This command will run the equivalent IDL commands inside of an isolated IDL environment using only those CHIANTI files.
Installing SSW is not required.
Note also that the `hissw <https://wtbarnes.github.io/hissw/>`_ package is required to run the IDL comparison tests.

It is not currently possible to run these tests in a continuous integration environment because of the licensing
restrictions imposed by IDL. However, this test suite will be run periodically, particularly as new versions of
CHIANTI are supported, to ensure that any differences between fiasco and the IDL software are well understood.

.. _additional-resources:

Additional Resources
--------------------

Each version of the CHIANTI database is described in detail in a set of publications. These are listed below:

- Version 1: :cite:t:`dere_chianti_1997`
- Version 2: :cite:t:`landi_chianti_1999`
- Version 3: :cite:t:`dere_chianti-atomic_2001`
- Version 4: :cite:t:`young_chianti-atomic_2003`
- Version 5: :cite:t:`landi_chianti-atomic_2006`
- Version 6: :cite:t:`dere_chianti_2009`
- Version 7: :cite:t:`landi_chiantiatomic_2012`
- Version 7.1: :cite:t:`landi_chiantiatomic_2013`
- Version 8: :cite:t:`del_zanna_chianti_2015`
- Version 9: :cite:t:`dere_chiantiatomic_2019`
- Version 10: :cite:t:`del_zanna_chiantiatomic_2021`

.. _CHIANTI Atomic Database: http://www.chiantidatabase.org/
