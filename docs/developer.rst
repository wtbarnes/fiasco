.. _fiasco-developer-guide:

Developer Guide
===============

This page contains useful information for those developing `fiasco`.

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

Test Against IDL Routines
-------------------------

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

.. _fiasco-how-to-support-new-database-version:

Add Support for a New Version of the CHIANTI Database
-----------------------------------------------------

Adding support for a new version of the CHIANTI database is a multi-step process.
Often, new versions of the database (and software) introduce new features for modeling the solar spectra.
These kinds of changes have to be dealt with on a case-by-case basis as new versions of the database and software introduce features of varying complexity.

At a minimum, compatibility with a new version of the database means being able to successfully parse the files within the database in order to provide equivalent functionality to the previous version of the database.
This is how we will define "support" for a new version of the database.
Note that "support" for a particular version of the database does not necessarily mean complete parity with the accompanying IDL routines for that version of the database.
Furthermore, as support for each new version of the CHIANTI database is added, fiasco must guarantee backwards compatibility with all past versions of the database for which it has provided support.

Below is a list of tasks that should be performed when supporting a new version of the database.
Note that this list is non-exhaustive and is the minimum number of steps required to support a new version.
Support for new versions may require additional steps depending on the complexity of the changes from one database version to the next.
It is probably easiest to do most of this in a single pull request though it is not required.

- Add the new version number(s) to the ``SUPPORTED_VERSIONS`` list in :file:`fiasco/util/setup_db.py`.
  The number of new supported versions should include all versions between the last supported version and the most current version (including minor or bugfix versions as appropriate) for which support is being added.
  As an example, if support for version 10 is being added and the most recent release of v10 of the database is 10.1.1, all versions matching the pattern ``10.0.x`` and ``10.1.x`` should be added to the supported versions list.
  Note that the last entry in ``SUPPORTED_VERSIONS`` will be treated as the latest version and the default version downloaded.
  As such, this list should be ordered according to increasing version number.
- Use :file:`tools/generate_test_file_list.py` to add any files needed for testing to the list of files included in the test database.
  Note that you should avoid editing :file:`fiasco/tests/data/test_file_list.json` by hand.
- Generate the list of file hashes for the new version of the database using :file:`tools/generate_hash_table.py`.
  Note that a new list needs to be generated only for the latest version for which support is being added.
  For example, if ``10.0.0``, ``10.0.1``, and ``10.1.0`` are all being added to the list of supported versions, a list of file hashes only needs to be generated for ``10.1.0``.
  This is because the file hashes are only used for the test database and automated testing is only run on the latest version of each major release of the database.
- Run the IDL comparison tests locally for the new database version following the instructions in :ref:`fiasco-how-to-run-tests-idl`.
  It is necessary to run the tests locally for two reasons:

    1. Running the IDL comparison tests generates the cached results for the new version of the database and stores them in the :file:`fiasco/tests/idl/data` directory.
       These cached results are then used to run these comparison tests on the CI such that an IDL installation is not required.
    2. There are a few IDL comparison tests which always require an IDL installation to run as their output is too large to store in the `fiasco` repository.
       As such, these must be run locally to ensure they pass with the latest version of the database.

  A few additional points regarding running these tests with a new version of the database:

  - When running these tests, you should use the latest version of the CHIANTI IDL software that is compatible with the version of the database you are testing.
    The versions of the CHIANTI IDL software are most easily accessed `here <https://github.com/chianti-atomic/chianti-idl/tags>`__.
  - The CHIANTI IDL software relies on some components of SSW.
    To avoid a dependency on SSW, the minimum set of routines that the IDL software depends on is provided in :file:`fiasco/tests/idl/ssw_gen_functions`.
    When running these tests using a new version of the software, you may find that additional SSW routines are required.
    These should be added to :file:`fiasco/tests/idl/ssw_gen_functions` as needed.
  - The contribution function comparison tests in :file:`fiasco/tests/idl/test_idl_goft.py` require manually selecting the wavelength of the line of interest as there is no direct way to specify a wavelength in physical units when computing the contribution function in the IDL software.
    When running the test, if there is no cached file present, the IDL GUI wavelength selector should automatically be displayed.
    You should select the wavelength appropriate for that test as shown in the test parametrization in :file:`fiasco/tests/idl/test_idl_goft.py`.
    Once the result is cached to a file, this will no longer be necessary to do when rerunning the tests.

- Add a job to the CI configuration in :file:`.github/workflows/ci.yml` that explicitly tests the previous version of the database.
  This should follow the example of the jobs that are titled following the pattern ``test_database_v<version_number>``.
  For example, if you are adding support for version 10 of the database such that version 10 is now the default version used when running the test suite, you should add a job titled ``test_database_v9``.
  This job should include the argument ``--ascii-dbase-url`` which points to the URL of the version of the database that was previously supported by default.
  The purpose of these jobs is to ensure backwards compatibility with previous releases of the database.
  Note that these jobs are only run on release branches, tags, and PRs with the "Run v<version_number> tests" label.
- Update the banner at the top of the documentation page to indicate that support has been added for a new version.
  This is configured in :file:`docs/conf.py`.
