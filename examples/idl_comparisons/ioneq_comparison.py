"""
CHIANTI IDL Comparison: Ionization Fraction
===============================================

Compare the ionization fraction calculation to that in the CHIANTI IDL routines for a few select ions.
"""
import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import quantity_support
from jinja2 import Template

import fiasco

from fiasco.tests.idl.helpers import read_idl_test_output

quantity_support()

################################################
# Define a function for comparing and plotting
# the outputs from fiasco and IDL. Note that we have precomputed
# the IDL outputs.
def plot_idl_comparison(x, y_idl, y_python, fig, n_rows, i_row, quantity_name, thresh=1e-6):
    # Highlight region where value of goft is within some value of the peak
    idx = np.where(y_idl>=y_idl.max()*thresh)
    x_thresh_min, x_thresh_max = x[idx[0][[0,-1]]]
    # Direct comparison
    ax = fig.add_subplot(n_rows, 3, i_row+1)
    ax.plot(x, y_python, label='fiasco')
    ax.plot(x, y_idl, color='k', marker='o', ls='--', markevery=2, label='IDL')
    ax.axvspan(x[0], x_thresh_min, color='r', alpha=0.25)
    ax.axvspan(x_thresh_max, x[-1], color='r', alpha=0.25)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if i_row == 0:
        ax.set_title('Ionization Fraction')
    ax.set_xlim(x[[0,-1]])
    ax.set_ylim(y_python.max()*np.array([thresh,2]))
    ax.set_ylabel(quantity_name)
    plt.legend()
    # Ratio
    ax = fig.add_subplot(n_rows, 3, i_row+2, sharex=ax)
    plt.plot(x, (y_python/y_idl).decompose())
    ax.axvspan(x[0], x_thresh_min, color='r', alpha=0.25)
    ax.axvspan(x_thresh_max, x[-1], color='r', alpha=0.25)
    ax.axhline(y=1, color='k', ls=':')
    if i_row == 0:
        ax.set_title('fiasco / IDL')
    diff_limit = .001
    ax.set_ylim(1-diff_limit, 1+diff_limit)
    # Normalized difference
    ax = fig.add_subplot(n_rows, 3, i_row+3, sharex=ax)
    ax.plot(x,  ((y_python - y_idl)/y_idl).decompose())
    ax.axvspan(x[0], x_thresh_min, color='r', alpha=0.25)
    ax.axvspan(x_thresh_max, x[-1], color='r', alpha=0.25)
    ax.axhline(y=0, color='k', ls=':')
    if i_row == 0:
        ax.set_title('(fiasco - IDL)/IDL')
    ax.set_ylim(-diff_limit, diff_limit)


################################################
# Next, plot the comparison between the ionization fraction results
# for a selected number of ions. In the case of fiasco, we will also
# compute the ionization fraction in addition to reading the tabulated
# values.
ioneq_files = [
    'ioneq_1_1',
    'ioneq_6_1',
    'ioneq_6_2',
    'ioneq_6_3',
    'ioneq_20_2',
    'ioneq_26_5',
    'ioneq_26_16',
    'ioneq_26_18',
    'ioneq_26_20',
    'ioneq_26_27',
]
fig = plt.figure(figsize=(9,3*len(ioneq_files)), layout='constrained')
for i, name in enumerate(ioneq_files):
    idl_result = read_idl_test_output(name, '8.0.7')
    ion = fiasco.Ion((idl_result['Z'], idl_result['iz']),
                     idl_result['temperature'],
                     ioneq_filename=idl_result['ioneq_filename'])
    print(f'IDL code to produce ioneq result for {ion.ion_name_roman}:')
    print(Template(idl_result['idl_script']).render(**idl_result))
    plot_idl_comparison(ion.temperature, idl_result['ioneq'], ion.ioneq,
                        fig, len(ioneq_files), 3*i, f'{ion.ion_name_roman}')
plt.show()
