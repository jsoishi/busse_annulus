import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from dedalus.extras import plot_tools
from dedalus.tools.general import natural_sort

import pathlib

plt.ioff()
plt.style.use('prl')

def calc_dns_power(x, y, data):
    power = (data*data.conj()).real
    y_out = np.round(y/np.pi + 0.5).astype(np.int)

    return x, y_out, np.log10(power)

def calc_ce2_power(x, y, data):
    zero_data = data[:,0]
    data[:,0] = (zero_data*zero_data.conj()).real
    y_out = (y/np.pi + 0.5).astype(np.int)
    
    return x, y_out, np.log10(np.abs(data))


output = '../../figs'
field = 'zeta'
climits=(-1,7)
xlimits=(0,63)
ylimits=(1,63)
dns_tidx=-1
ce2_tidx=-1
label=None
dns_filename="../scratch/busse_annulus_ra4.00e+07_beta3.16e+05_C3.16e-01_Pr1.00e+00_filter5.00e-01_nx512_ny128_CFL/snapshots/snapshots_s2.h5"
ce2_filename="run_R_explicit_powerspec/data_snapshots/data_snapshots_s121.h5"
ce2_b_filename="run_S_explicit_no_mean_subtr/data_snapshots/data_snapshots_s68.h5"

with h5py.File(dns_filename,'r') as dns_data, h5py.File(ce2_filename,'r') as ce2_data, h5py.File(ce2_b_filename,'r') as ce2_b_data:
    nrows, ncols = 1, 3
    image = plot_tools.Box(1,1)
    pad = plot_tools.Frame(0.2, 0.1, 0.1, 0.1)
    margin = plot_tools.Frame(0.1, 0.1, 0.1, 0.1)

    scale = 4
    mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
    fig = mfig.figure

    dns_axes = mfig.add_axes(0, 0, [0, 0, 1, 1])
    if field == 'zeta':
        power_axis_label = r"$\log_{10} |\hat{\zeta}|^2$"
    elif field == 'theta':
        power_axis_label = r"$\log_{10} |\hat{\theta}|^2$"
    else:
        raise NotImplementedError("Field can only be theta or zeta, not {}".format(field))

    # DNS
    dns_image_axes = (2,1)
    dns_image_slices = (dns_tidx,slice(None), slice(None))
    dns_kspace_data = dns_data['tasks/{}_kspace'.format(field)]
    dns_time = dns_data['scales/sim_time'][dns_tidx]
    pax, cax = plot_tools.plot_bot(dns_kspace_data,dns_image_axes, dns_image_slices,func=calc_dns_power,clim=climits, cmap='viridis',axes=dns_axes)
    pax.xaxis.set_major_locator(MaxNLocator(integer=True))
    #pax.yaxis.set_major_locator(MaxNLocator(integer=True))
    pax.set_xlim(xlimits)
    pax.set_ylim(ylimits)
    cax.set_xlabel(power_axis_label)
    pax.set_xlabel(r"$k_x$")
    pax.set_ylabel(r"$k_{y}$")
    pax.text(0.9,0.85, 'a',
             bbox={'facecolor': 'grey', 'alpha': 0.5, 'boxstyle':'round'},transform=pax.transAxes, fontsize=18, color='white')
    pax.set_rasterized(True)
    # CE2
    ce2_axes = mfig.add_axes(0, 1, [0, 0, 1, 1])
    ce2_image_axes = (1,3)
    ce2_image_slices = (ce2_tidx,slice(None), 0,slice(None))
    ce2_kspace_data = ce2_data['tasks/{}_power'.format(field)]
    ce2_time  = ce2_data['scales/sim_time'][ce2_tidx]
    ce2_write = ce2_data['scales/write_number'][ce2_tidx]
    pax, cax = plot_tools.plot_bot(ce2_kspace_data,ce2_image_axes, ce2_image_slices,func=calc_ce2_power, clim=climits,cmap='viridis', axes=ce2_axes)

    cax.set_xlabel(power_axis_label)
    cax.set_visible(False)
    pax.xaxis.set_major_locator(MaxNLocator(integer=True))

    pax.set_xlim(xlimits)
    pax.set_ylim(ylimits)
    pax.set_xlabel(r"$k_x$")
    pax.set_ylabel(r"$k_{y}$")
    pax.text(0.9,0.85, 'b',
             bbox={'facecolor': 'grey', 'alpha': 0.5, 'boxstyle':'round'},transform=pax.transAxes, fontsize=18, color='white')

    pax.set_rasterized(True)

    
    # CE2
    ce2_b_axes = mfig.add_axes(0, 2, [0, 0, 1, 1])
    ce2_image_axes = (1,3)
    ce2_image_slices = (ce2_tidx,slice(None), 0,slice(None))
    ce2_b_kspace_data = ce2_b_data['tasks/{}_power'.format(field)]
    ce2_b_time  = ce2_b_data['scales/sim_time'][ce2_tidx]
    ce2_b_write = ce2_b_data['scales/write_number'][ce2_tidx]
    pax, cax = plot_tools.plot_bot(ce2_b_kspace_data,ce2_image_axes, ce2_image_slices,func=calc_ce2_power, clim=climits,cmap='viridis', axes=ce2_b_axes)

    cax.set_xlabel(power_axis_label)
    cax.set_visible(False)
    pax.xaxis.set_major_locator(MaxNLocator(integer=True))

    pax.set_xlim(xlimits)
    pax.set_ylim(ylimits)
    pax.set_xlabel(r"$k_x$")
    pax.set_ylabel(r"$k_{y}$")
    pax.text(0.9,0.85, 'c',
             bbox={'facecolor': 'grey', 'alpha': 0.5, 'boxstyle':'round'},transform=pax.transAxes, fontsize=18, color='white')

    pax.set_rasterized(True)

    outputdir = pathlib.Path(output)
    for ext in ['png']:#,'pdf']:
        outfilen = "power_spectra_{}_dns_run_R.{}".format(field,ext)
        fig.savefig(str(outputdir/outfilen), dpi=400)
    
