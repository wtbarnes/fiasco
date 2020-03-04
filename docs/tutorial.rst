Tutorial
=========

Parsing Raw CHIANTI Data
------------------------
While the main advantages of fiasco lie in its high-level interface to the CHIANTI data, some users may wish to just parse the raw data directly. fiasco provides a convenient interface to parsing any of the raw CHIANTI datafiles and provides detailed metadata about each datafile. Specifically, fiasco returns a `~astropy.table.QTable` object with appropriate units and descriptive names attached to each column.

For example, say we want to parse the energy level file for Fe V (i.e. iron with four electrons missing)

    >>> from fiasco.io import Parser
    >>> p = Parser('fe_5.elvlc', ascii_dbase_root=getfixture('ascii_dbase_root'))
    >>> table = p.parse()
    >>> table # doctest: +SKIP
    <QTable length=34>
    level  config label multiplicity L_label    J     E_obs      E_th
                                                      1 / cm    1 / cm
    int64   str7   str1    int64       str1  float64 float64   float64
    ----- ------- ----- ------------ ------- ------- -------- ----------
        1     3d4                  5       D     0.0      0.0        0.0
        2     3d4                  5       D     1.0    142.1    153.632
        3     3d4                  5       D     2.0    417.3    449.923
        4     3d4                  5       D     3.0    803.1    873.509
        5     3d4                  5       D     4.0   1282.8    1407.93
        6 3d4 (2)                  3       P     0.0  24055.4  24398.996
        7     3d4                  3       H     4.0  24932.5  28643.637
        8 3d4 (2)                  3       P     1.0  24972.9  25314.205
        9     3d4                  3       H     5.0  25225.9  29013.449
       10     3d4                  3       H     6.0  25528.5  29396.432
       11 3d4 (2)                  3       P     2.0  26468.3  26887.838
      ...     ...   ...          ...     ...     ...      ...        ...
       24 3d4 (2)                  1       D     2.0  46291.2  49168.906
       25     3d4                  1       F     3.0  52732.7  59302.051
       26 3d4 (1)                  3       P     2.0  61854.4  67320.555
       27 3d4 (1)                  3       F     4.0  62238.1  67484.062
       28 3d4 (1)                  3       F     2.0  62321.1  67574.047
       29 3d4 (1)                  3       F     3.0  62364.4  67616.844
       30 3d4 (1)                  3       P     1.0  62914.2  68356.477
       31 3d4 (1)                  3       P     0.0  63420.0  68879.922
       32 3d4 (1)                  1       G     4.0  71280.3  77939.836
       33 3d4 (1)                  1       D     2.0  93832.3 103239.773
       34 3d4 (1)                  1       S     0.0 121130.2 128805.281

The individual columns can easily be accessed as,

    >>> table['E_obs']
    <Quantity [     0. ,    142.1,    417.3,    803.1,   1282.8,  24055.4,
                24932.5,  24972.9,  25225.9,  25528.5,  26468.3,  26760.7,
                26842.3,  26974. ,  29817.1,  30147. ,  30430.1,  36586.3,
                36630.1,  36758.5,  36925.4,  37511.7,  39633.4,  46291.2,
                52732.7,  61854.4,  62238.1,  62321.1,  62364.4,  62914.2,
                63420. ,  71280.3,  93832.3, 121130.2] 1 / cm>

Each above column is an `~astropy.units.Quantity` object with units attached to it if appropriate. Metadata, including the original footer from the raw CHIANTI data and detailed descriptions of each of the columns, is included with each table,

    >>> table.meta.keys()
    odict_keys(['footer', 'chianti_version', 'filename', 'descriptions', 'element', 'ion', 'dielectronic'])
    >>> print(table.meta['footer'])
    filename:  fe_5.elvlc
    Observed energies: Ralchenko, Yu., Kramida, A.E., Reader, J., and NIST ASD Team (2008).
    NIST Atomic Spectra Database (version 3.1.5), [Online]. Available: http://physics.nist.gov/asd3
    [2009, September 1]. National Institute of Standards and Technology, Gaithersburg, MD.
    Theoretical energies: Ballance, C.P., Griffin, D.C., & McLaughlin, B.M. 2007, J.Phys.B, 40, F327
    produced as part of the Arcetri/Cambridge/NRL 'CHIANTI' atomic data base collaboration
    Peter Young, 3-Sep-2009

More details for how to access data from these tables can be found in the `Astropy Table documentation`_.

Note that in the above example we only provided the filename because fiasco is able to infer the path to the file using the information in `~fiasco.defaults`. However, the parser can be used on standalone files as well, i.e. CHIANTI data files that are outside the directory specificied in `fiasco.defaults['ascii_dbase_root']`,

    >>> p = Parser('/path/to/standalone/chianti/data/fe_5.elvlc') # doctest: +SKIP
    >>> table = p.parse() # doctest: +SKIP

Ions
------
The CHIANTI database is organized around individual ions, with multiple types of datafiles attached to each ion. In keeping with the organization of the database, the primary unit of the fiasco library is the `~fiasco.Ion` object which can be created in the following way,

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from fiasco import Ion
    >>> temperature = np.logspace(5, 7, 100) * u.K
    >>> ion = Ion('Fe 15', temperature, hdf5_dbase_root=getfixture('hdf5_dbase_root'))

This creates a `~fiasco.Ion` object for the Fe XV ion. Note also the same object can also be created in the following ways,

    >>> ion = Ion('iron 15', temperature, hdf5_dbase_root=getfixture('hdf5_dbase_root'))
    >>> ion = Ion('iron 14+', temperature, hdf5_dbase_root=getfixture('hdf5_dbase_root'))


Basic Metadata
***************
The `~fiasco.Ion` object holds several basic pieces of metadata related to the particular ion,

    >>> ion.element_name
    'iron'
    >>> ion.atomic_symbol
    'Fe'
    >>> ion.atomic_number
    26
    >>> ion.ion_name
    'Fe 15'
    >>> ion.charge_state
    14
    >>> ion.ionization_stage
    15
    >>> ion.abundance
    <Quantity 3.16227766e-05>

In the cases of the abundance and ionization potential (`ip`), specific datasets available in CHIANTI can be specified using the `abundance_filename` and `ip_filename` keywords, respectively, when instantiating the ion object. For more information on accessing specific datasets, see :ref:`ionbase` section.

The equilibrium population fractions are interpolated to `temperature` and can be accessed using the `ioneq` keyword.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from fiasco import Ion
    ion = Ion('Fe 15', np.logspace(5, 7, 100) * u.K)
    plt.plot(ion.temperature, ion.ioneq)
    plt.xlabel(r'$T [K]$')
    plt.ylabel(r'Population Fraction')
    plt.xscale('log')
    plt.show()

Note that these population fractions returned by `Ion.ioneq` are stored in the CHIANTI database and and therefore are set to NaN for temperatures outside of the temperature data given in CHIANTI. If you need to calculate the population fractions over a wider temperature range, it is better to do so by calculating the ionization and recombination rates. See the :ref:`elements` section for more info.

Derived Quantities
******************
In addition to providing an API to the CHIANTI data, `Ion` also provides several methods for computing derived quantities from the data. These include the ionization and recombination rates. 

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from fiasco import Ion
    ion = Ion('Fe 18', np.logspace(4, 8, 100) * u.K)
    plt.plot(ion.temperature, ion.recombination_rate(), label='Total')
    plt.plot(ion.temperature, ion.dielectronic_recombination_rate(), label='Dielectronic')
    plt.plot(ion.temperature, ion.radiative_recombination_rate(), label='Radiative')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel(f'Recombination Rate')
    plt.xlabel(f'Temperature [K]')
    plt.legend()
    plt.show()

Accessing Raw Data
******************

Working with Multiple Ions
--------------------------

.. _elements:

Elements
---------

.. _ionbase:

The IonBase Object
------------------


.. _Astropy Table documentation: http://docs.astropy.org/en/stable/table/access_table.html