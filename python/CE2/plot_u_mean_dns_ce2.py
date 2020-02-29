"""
Plot

Usage:
    plot_u_mean_dns_ce2.py <dns_file_path> <ce2_file_path> [--output=<dir> --tavg=<tavg> --label=<label> --tavg_scale=<tavg_scale> --xlim=<xlim>]

Options:
    --output=<dir>               Output directory [default: ../../figs]
    --tavg=<tavg>                Time to average over (counted backwards from end of simulation) [default: 1]
    --label=<label>              label for plot 
    --tavg_scale=<tavg_scale>    scale for CE2 time [default: 0.1]
    --xlim=<xlim>                limits for x axis
"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from mpi4py import MPI

from dedalus.extras import plot_tools
from dedalus.tools.general import natural_sort

import pathlib
from docopt import docopt
from dedalus.tools.parallel import Sync

plt.style.use('prl')

def time_avg(data, task, scales_name='y/1',tavg_scale=1):
    time = data['scales/sim_time'][:]
    tmax = time.max()
    y = data['scales/{}'.format(scales_name)][:]
    time_mask = (time.max() - time) < tavg*tavg_scale
    umean = data['tasks'][task][time_mask,:,0].mean(axis = 0)

    return y, umean

args = docopt(__doc__)
tavg = float(args['--tavg'])
tavg_scale = float(args['--tavg_scale'])
label = args['--label']
xlim = args['--xlim']
if xlim:
    xlim = xlim.split(',')
    xlim = [float(i) for i in xlim]


# Create output directory if needed
output_path = pathlib.Path(args['--output']).absolute()

with Sync() as sync:
    if sync.comm.rank == 0:
        if not output_path.exists():
            output_path.mkdir()

dns_data = h5py.File(args['<dns_file_path>'],'r')
ce2_data = h5py.File(args['<ce2_file_path>'],'r')

y_dns, dns_umean = time_avg(dns_data, '<u>_x')
y_ce2, ce2_umean = time_avg(ce2_data, 'cu', scales_name='y1/1.0',tavg_scale=tavg_scale)

plt.plot(dns_umean, y_dns, label='DNS')
plt.plot(ce2_umean[0,:], y_ce2, '--',label='CE2')
plt.legend(loc='upper left')
plt.ylabel(r"$y$",rotation='horizontal')
plt.xlabel(r"$\left< u \right>_{x}$")
plt.ylim(0,1)
plt.xlim(xlim[0],xlim[1])
plt.tight_layout()

outputdir = pathlib.Path(args['--output'])
for ext in ['png','pdf']:
    outfilen = "umean_dns_ce2_{}.{}".format(label,ext)
    plt.savefig(str(outputdir/outfilen), dpi=400)
