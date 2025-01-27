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
    ax1 = fig.add_subplot(n_rows, 3, i_row+1)
    ax1.plot(x, y_python, label='fiasco')
    ax1.plot(x, y_idl, color='k', marker='o', ls='--', markevery=2, label='IDL')
    ax1.axvspan(x[0], x_thresh_min, color='r', alpha=0.25)
    ax1.axvspan(x_thresh_max, x[-1], color='r', alpha=0.25)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    if i_row == 0:
        ax1.set_title('Ionization Fraction')
    ax1.set_xlim(x[[0,-1]])
    ax1.set_ylim(y_python.max()*np.array([thresh,2]))
    ax1.set_ylabel(quantity_name)
    # Ratio
    ax2 = fig.add_subplot(n_rows, 3, i_row+2, sharex=ax1)
    ax2.plot(x, (y_python/y_idl).decompose())
    ax2.axvspan(x[0], x_thresh_min, color='r', alpha=0.25)
    ax2.axvspan(x_thresh_max, x[-1], color='r', alpha=0.25)
    ax2.axhline(y=1, color='k', ls=':')
    if i_row == 0:
        ax2.set_title('fiasco / IDL')
    diff_limit = .001
    ax2.set_ylim(1-diff_limit, 1+diff_limit)
    # Normalized difference
    ax3 = fig.add_subplot(n_rows, 3, i_row+3, sharex=ax1)
    ax3.plot(x,  ((y_python - y_idl)/y_idl).decompose())
    ax3.axvspan(x[0], x_thresh_min, color='r', alpha=0.25)
    ax3.axvspan(x_thresh_max, x[-1], color='r', alpha=0.25)
    ax3.axhline(y=0, color='k', ls=':')
    if i_row == 0:
        ax3.set_title('(fiasco - IDL)/IDL')
    ax3.set_ylim(-diff_limit, diff_limit)
    return ax1, ax2, ax3


################################################
# Next, plot the comparison between the ionization fraction results
# for a selected number of ions. In the case of fiasco, we will also
# compute the ionization fraction in addition to reading the tabulated
# values. The regions highlighted in red denote places where the ionization
# fraction is less than :math:`10^{-6}` times the peak of the ionization fraction.
# In these regions, the comparison is between the two results is not as critical
# as differences due to the interpolation schemes may show large deviations between
# the two approaches. However, the ionization fraction in these regions does not
# meaningfully contribute to any other derived quantities.
ionization_files = [
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
fig = plt.figure(figsize=(9,3*len(ionization_files)), layout='constrained')
for i, name in enumerate(ionization_files):
    idl_result = read_idl_test_output(name, '8.0.7')
    ion = fiasco.Ion((idl_result['Z'], idl_result['iz']),
                     idl_result['temperature'],
                     ionization_filename=idl_result['ioneq_filename'])
    element = fiasco.Element(ion.atomic_symbol, ion.temperature)
    ionization_fraction = element.equilibrium_ionization
    print(f'IDL code to produce ionization_fraction result for {ion.ion_name_roman}:')
    print(Template(idl_result['idl_script']).render(**idl_result))
    axes = plot_idl_comparison(ion.temperature, idl_result['ioneq'], ion.ionization_fraction,
                               fig, len(ionization_files), 3*i, f'{ion.ion_name_roman}')
    axes[0].plot(element.temperature, ionization_fraction[:, ion.charge_state],
                 label='fiasco (rates)', color='C1', ls='-')
    axes[0].legend()
plt.show()
