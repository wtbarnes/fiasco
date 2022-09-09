"""
IDL comparison tests for continuum calculations
"""
import os

import numpy as np
import astropy.units as u

import fiasco


def test_idl_compare_free_free(idl_env, ascii_dbase_root, hdf5_dbase_root):
    # Inputs
    temperature = np.logspace(5, 8, 61) * u.K
    wavelength = np.arange(25, 414, 1) * u.Angstrom
    abundance_name =  'sun_coronal_1992_feldman.abund'
    ioneq_name = 'chianti'
    # Run fiasco free-free calculation
    ion_kwargs = {
        'abundance_filename': abundance_name,
        'ioneq_filename': ioneq_name,
        'hdf5_dbase_root': hdf5_dbase_root,
    }
    all_ions = [fiasco.Ion(i, temperature, **ion_kwargs) for i in fiasco.list_ions(hdf5_dbase_root)]
    all_ions = fiasco.IonCollection(*all_ions)
    # Run IDL free-free calculation
    script = """
    ; set common block
    common elements, abund, abund_ref, ioneq, ioneq_logt, ioneq_ref

    ; read abundance and ionization equilibrium
    abundfile='{{ abundance_filename }}'
    read_abund, abundfile, abund, abund_ref
    ioneqfile="{{ ioneq_filename }}"
    read_ioneq, ioneqfile, ioneq_logt, ioneq, ioneq_ref

    ; set temperature and wavelength
    temperature = {{ temperature | to_unit('K') | list }}
    wavelength = {{ wavelength | to_unit('Angstrom') | list }}

    ; calculate free-free
    freefree, temperature, wavelength, ff, /no_setup

    ; The output of freefree is scaled by 10^40
    ff = ff/1d40
    """
    temperature = np.logspace(5, 8, 61) * u.K
    wavelength = np.arange(25, 414, 1) * u.Angstrom
    abundance_filename =  ascii_dbase_root / 'abundance' / f'{abundance_name}.abund'
    ioneq_filename = ascii_dbase_root / 'ioneq' / f'{ioneq_name}.ioneq'
    args = {
        'wavelength': wavelength,
        'temperature': temperature,
        'abundance_filename': abundance_filename,
        'ioneq_filename': ioneq_filename,
    }
    res = idl_env.run(script, args=args, verbose=True)
    assert res['ff']
