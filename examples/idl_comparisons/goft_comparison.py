"""
CHIANTI IDL Comparison: Contribution Functions
===============================================

Compare the contribution function calculation to that in the CHIANTI IDL routines for a few select ions.
"""
import hissw
import matplotlib.pyplot as plt
import numpy as np

from astropy.visualization import quantity_support

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
    ax1.set_xlim(1e5, 1e8)
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
    diff_limit = .05
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
# Next, plot the comparison between the contribution function results
# for a selected number of ions. The regions highlighted in red denote
# places where the contribution function is less than :math:`10^{-6}` times
# the peak of the contribution function. In these regions, the comparison is
# between the two results is not as critical as the contribution function
# is comparatively small in these regions.
goft_files = [
    'goft_20_15_200.972',
    'goft_26_9_171.073',
    'goft_26_9_188.496',
    'goft_26_11_188.497',
    'goft_26_14_197.862',
    'goft_26_16_262.984',
]
fig = plt.figure(figsize=(9,3*len(goft_files)), layout='constrained')
template_env = hissw.Environment(ssw_home='', idl_home='').env
for i, name in enumerate(goft_files):
    idl_result = read_idl_test_output(name, '8.0.7')
    ion = fiasco.Ion((idl_result['Z'], idl_result['iz']),
                     idl_result['temperature'],
                     abundance_filename=idl_result['abundance'],
                     ionization_fraction=idl_result['ioneq'])
    contribution_func = ion.contribution_function(idl_result['density'])
    transitions = ion.transitions.wavelength[~ion.transitions.is_twophoton]
    idx = np.argmin(np.abs(transitions - idl_result['wavelength']))
    # NOTE: Multiply by 0.83 because the fiasco calculation does not include the n_H/n_e ratio
    goft = contribution_func[:, 0, idx] * 0.83
    line_label = f'{ion.ion_name_roman} {idl_result["wavelength"]:latex_inline}'
    axes = plot_idl_comparison(
        ion.temperature,
        idl_result['contribution_function'],
        goft,
        fig,
        len(goft_files),
        3*i,
        line_label,
    )
    axes[0].legend()
    print(f'IDL code to produce {line_label} contribution function result:')
    print(template_env.from_string(idl_result['idl_script']).render(**idl_result))
