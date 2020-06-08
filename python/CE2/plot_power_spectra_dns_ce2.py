"""plot_power_spectra_dns_ce2.py
Plot power spectra from DNS and CE2

Must be run serially

Usage:
    plot_power_spectra_dns_ce2.py <dns_file_path> <ce2_file_path> [--output=<dir> --tavg=<tavg> --label=<label> --climits=<climits> --xlimits=<xlimits> --ylimits=<ylimits>]

Options:
    --output=<dir>      Output directory [default: ../../figs]
    --label=<label>     label for plot 
    --climits=<limits>  color limits  [default: 4,6]
    --xlimits=<limits>  x limits  [default: 0,15]
    --ylimits=<limits>  y limits  [default: 1,15]
    --tavg=<tavg>       time to average over [default: 1]
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

nrows, ncols = 1, 2
image = plot_tools.Box(1,1)
pad = plot_tools.Frame(0.2, 0.1, 0.1, 0.1)
margin = plot_tools.Frame(0.1, 0.1, 0.1, 0.1)

scale = 4
mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
fig = mfig.figure

dns_axes = mfig.add_axes(0, 0, [0, 0, 1, 1])

# DNS
dns_image_axes = (2,1)
dns_image_slices = (-1,slice(None), slice(None))
pax, cax = plot_tools.plot_bot(dns_data['tasks/zeta_kspace'],dns_image_axes, dns_image_slices,func=calc_dns_power,clim=climits, cmap='viridis',axes=dns_axes)
pax.xaxis.set_major_locator(MaxNLocator(integer=True))
#pax.yaxis.set_major_locator(MaxNLocator(integer=True))
pax.set_xlim(xlimits)
pax.set_ylim(ylimits)
cax.set_xlabel(r"$\log_{10} |\hat{\zeta}|^2$")
pax.set_xlabel(r"$k_x$")
pax.set_ylabel(r"$k_{y}$")

# CE2
ce2_axes = mfig.add_axes(0, 1, [0, 0, 1, 1])
image_axes = (1,3)
image_slices = (-1,slice(None), 0,slice(None))
pax, cax = plot_tools.plot_bot(ce2_data['tasks/zeta_power'],image_axes, image_slices,func=calc_ce2_power, clim=climits,cmap='viridis', axes=ce2_axes)
cax.set_xlabel(r"$\log_{10} |\hat{\zeta}|^2$")
pax.xaxis.set_major_locator(MaxNLocator(integer=True))
#pax.yaxis.set_major_locator(MaxNLocator(integer=True))
pax.set_xlim(xlimits)
pax.set_ylim(ylimits)
pax.set_xlabel(r"$k_x$")
pax.set_ylabel(r"$k_{y}$")

#pax.text(0.95,0.85, 'b',
#         bbox={'facecolor': 'grey', 'alpha': 0.5, 'boxstyle':'round'},transform=pax.transAxes, fontsize=18)
pax.set_rasterized(True)
outputdir = pathlib.Path(args['--output'])
for ext in ['png','pdf']:
    outfilen = "power_spectra_dns_ce2_{}.{}".format(label,ext)
    fig.savefig(str(outputdir/outfilen), dpi=400)
