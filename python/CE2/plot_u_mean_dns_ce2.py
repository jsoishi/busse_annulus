"""
Plot

Usage:
    plot_u_mean_dns_ce2.py <dns_file_path> <ce2_file_path> [--output=<dir> --tavg=<tavg> --label=<label>]

Options:
    --output=<dir>  Output directory [default: ../../figs]
    --tavg=<tavg>   Time to average over (counted backwards from end of simulation) [default: 1]
    --label=<label> label for plot 

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

def time_avg(data, task, scales_name='y/1'):
    time = data['scales/sim_time'][:]
    tmax = time.max()
    y = data['scales/{}'.format(scales_name)][:]
    time_mask = (time.max() - time) < tavg
    umean = data['tasks'][task][time_mask,:,0].mean(axis = 0)

    return y, umean

args = docopt(__doc__)
tavg = float(args['--tavg'])
label = args['--label']

# Create output directory if needed
output_path = pathlib.Path(args['--output']).absolute()

with Sync() as sync:
    if sync.comm.rank == 0:
        if not output_path.exists():
            output_path.mkdir()

dns_data = h5py.File(args['<dns_file_path>'],'r')
ce2_data = h5py.File(args['<ce2_file_path>'],'r')

y_dns, dns_umean = time_avg(dns_data, '<u>_x')
y_ce2, ce2_umean = time_avg(ce2_data, 'cu', scales_name='y1/1.0')

plt.plot(y_dns, dns_umean, label='DNS')
plt.plot(y_ce2, ce2_umean[0,:], '--',label='CE2')
plt.legend(loc='upper left')
plt.xlabel(r"$y$", fontsize=18)
plt.ylabel(r"$\left< u \right>_{x}$", fontsize=18)

outputdir = pathlib.Path(args['--output'])
outfilen = "umean_dns_ce2_{}.png".format(label)
plt.savefig(str(outputdir/outfilen))
