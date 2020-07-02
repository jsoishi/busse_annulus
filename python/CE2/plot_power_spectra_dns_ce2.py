"""plot_power_spectra_dns_ce2.py
Plot power spectra from DNS and CE2

Must be run serially

Usage:
    plot_power_spectra_dns_ce2.py <dns_file_path> <ce2_file_path> [--output=<dir> --tavg=<tavg> --label=<label> --climits=<climits> --xlimits=<xlimits> --ylimits=<ylimits> --dns_tidx=<dns_tidx> --ce2_tidx=<ce2_tidx>]

Options:
    --output=<dir>          Output directory [default: ../../figs]
    --label=<label>         label for plot 
    --climits=<limits>      color limits  [default: 4,6]
    --xlimits=<limits>      x limits  [default: 0,15]
    --ylimits=<limits>      y limits  [default: 1,15]
    --dns_tidx=<dns_tidx>   dns time index [default: -1] 
    --ce2_tidx=<ce2_tidx>   ce2 time index [default: -1] 
    --tavg=<tavg>           time to average over [default: 1]
"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

plt.ioff()

from dedalus.extras import plot_tools
from dedalus.tools.general import natural_sort

import pathlib
from docopt import docopt

plt.style.use('prl')
args = docopt(__doc__)
tavg = float(args['--tavg'])
label = args['--label']
climits = args['--climits']
xlimits = args['--xlimits']
ylimits = args['--ylimits']
dns_tidx = int(args['--dns_tidx'])
ce2_tidx = int(args['--ce2_tidx'])

climits = tuple([int(i) for i in climits.split(',')])
xlimits = tuple([int(i) for i in xlimits.split(',')])
ylimits = tuple([int(i) for i in ylimits.split(',')])

def calc_dns_power(x, y, data):
    power = (data*data.conj()).real
    y_out = np.round(y/np.pi + 0.5).astype(np.int)

    return x, y_out, np.log10(power)

def calc_ce2_power(x, y, data):
    zero_data = data[:,0]
    data[:,0] = (zero_data*zero_data.conj()).real
    y_out = (y/np.pi + 0.5).astype(np.int)
    
    return x, y_out, np.log10(np.abs(data))

# Create output directory if needed
output_path = pathlib.Path(args['--output']).absolute()

if not output_path.exists():
    output_path.mkdir()

dns_data = h5py.File(args['<dns_file_path>'],'r')
ce2_data = h5py.File(args['<ce2_file_path>'],'r')

nrows, ncols = 1, 3
image = plot_tools.Box(1,1)
pad = plot_tools.Frame(0.2, 0.1, 0.1, 0.1)
margin = plot_tools.Frame(0.1, 0.1, 0.1, 0.1)

scale = 4
mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
fig = mfig.figure

dns_axes = mfig.add_axes(0, 0, [0, 0, 1, 1])

# DNS
dns_image_axes = (2,1)
dns_image_slices = (dns_tidx,slice(None), slice(None))
dns_kspace_data = dns_data['tasks/zeta_kspace']
dns_time = dns_data['scales/sim_time'][dns_tidx]
pax, cax = plot_tools.plot_bot(dns_kspace_data,dns_image_axes, dns_image_slices,func=calc_dns_power,clim=climits, cmap='viridis',axes=dns_axes)
pax.xaxis.set_major_locator(MaxNLocator(integer=True))
#pax.yaxis.set_major_locator(MaxNLocator(integer=True))
pax.set_xlim(xlimits)
pax.set_ylim(ylimits)
cax.set_xlabel(r"$\log_{10} |\hat{\zeta}|^2$")
pax.set_xlabel(r"$k_x$")
pax.set_ylabel(r"$k_{y}$")
pax.text(0.9,0.85, 'a',
         bbox={'facecolor': 'grey', 'alpha': 0.5, 'boxstyle':'round'},transform=pax.transAxes, fontsize=18, color='white')
pax.set_rasterized(True)
# CE2
ce2_axes = mfig.add_axes(0, 1, [0, 0, 1, 1])
ce2_image_axes = (1,3)
ce2_image_slices = (ce2_tidx,slice(None), 0,slice(None))
ce2_kspace_data = ce2_data['tasks/zeta_power']
ce2_time  = ce2_data['scales/sim_time'][ce2_tidx]
ce2_write = ce2_data['scales/write_number'][ce2_tidx]
pax, cax = plot_tools.plot_bot(ce2_kspace_data,ce2_image_axes, ce2_image_slices,func=calc_ce2_power, clim=climits,cmap='viridis', axes=ce2_axes)
cax.set_xlabel(r"$\log_{10} |\hat{\zeta}|^2$")
pax.xaxis.set_major_locator(MaxNLocator(integer=True))
#pax.yaxis.set_major_locator(MaxNLocator(integer=True))
pax.set_xlim(xlimits)
pax.set_ylim(ylimits)
pax.set_xlabel(r"$k_x$")
pax.set_ylabel(r"$k_{y}$")
pax.text(0.9,0.85, 'b',
         bbox={'facecolor': 'grey', 'alpha': 0.5, 'boxstyle':'round'},transform=pax.transAxes, fontsize=18, color='white')

pax.set_rasterized(True)

# plot only kx = 0 mode
kx0_axes = mfig.add_axes(0, 2, [0.03, 0, 0.94, 0.94])
kx, ky, data = plot_tools.get_plane(ce2_kspace_data, ce2_image_axes[0], ce2_image_axes[1], ce2_image_slices)
ce2_kx, ce2_ky, ce2_power = calc_ce2_power(kx, ky, data)

kx, ky, data = plot_tools.get_plane(dns_kspace_data, dns_image_axes[0], dns_image_axes[1], dns_image_slices)
dns_kx, dns_ky, dns_power = calc_dns_power(kx, ky, data)

kx0_axes.plot(dns_power[:,0], label='DNS t = {:5.2e}'.format(dns_time))
kx0_axes.plot(ce2_power[:,0], label='CE2 t = {:5.2e}'.format(ce2_time))

kx0_axes.set_xlabel(r"$k_{y}$")
kx0_axes.set_ylabel(r"$\log_{10} |\hat{\zeta}|^2$")
kx0_axes.set_xlim(1,15)
kx0_axes.text(0.9,0.85, 'c',
         bbox={'facecolor': 'grey', 'alpha': 0.5, 'boxstyle':'round'},transform=kx0_axes.transAxes, fontsize=18, color='black')
kx0_axes.xaxis.set_major_locator(MaxNLocator(integer=True))
kx0_axes.legend(loc='lower right')

outputdir = pathlib.Path(args['--output'])
for ext in ['png','pdf']:
    outfilen = "power_spectra_dns_ce2_{}_write_{}.{}".format(label,ce2_write,ext)
    fig.savefig(str(outputdir/outfilen), dpi=400)
