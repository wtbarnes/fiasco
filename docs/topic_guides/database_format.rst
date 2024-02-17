.. _fiasco-topic-guide-database-format:

The Database Format Used by fiasco
==================================

Analyzing astrophysical observations requires detailed knowledge of the underlying physics responsible for the formation of the observed spectra.
The CHIANTI atomic database provides a comprehensive set of atomic data for relevant elements, ions, and transitions in astrophysical plasmas.
CHIANTI includes energy levels, transition wavelengths, as well as collision strength and cross-section data for computing ionization and recombination rates.
The papers :ref:`fiasco-reference-chianti-papers` give an extensive description of each version of the database.

The CHIANTI data is distributed by the CHIANTI team as a collection of ASCII text files in a series of subdirectories.
Rather than interact with these raw text files directly, fiasco first builds the entire CHIANTI database into a single `HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ file to allow for easier and faster access to the data.
