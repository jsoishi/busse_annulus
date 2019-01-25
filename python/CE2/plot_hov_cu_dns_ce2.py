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

ce2_image_axes = [0, 3]
dns_image_axes = [0, 1]
dns_data_slices = (slice(None), slice(None), 0)
ce2_data_slices = (slice(None), 0, 0, slice(None))

image_scales = ['sim_time', 0]

nrows, ncols = 2, 1
image = plot_tools.Box(2, 1)
pad = plot_tools.Frame(0.2, 0.2, 0.1, 0.1)
margin = plot_tools.Frame(0.3, 0.2, 0.1, 0.1)

scale = 4
mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
fig = mfig.figure


axes = mfig.add_axes(0, 0, [0, 0, 1, 1])
# DNS
print("dns_data shape", dns_data['tasks/<u>_x'].shape)
print('ce2 data shape', ce2_data['tasks/cu'].shape)
plot_tools.plot_bot(dns_data['tasks/<u>_x'], dns_image_axes, dns_data_slices, image_scales, axes=axes, title='<u> DNS', even_scale=True)

axes = mfig.add_axes(1, 0, [0, 0, 1, 1])
# CE2
plot_tools.plot_bot(ce2_data['tasks/cu'], ce2_image_axes, ce2_data_slices, image_scales, axes=axes, title='<u> CE2', even_scale=True)
axes.set_xlim(0,2)

outputdir = pathlib.Path(args['--output'])
outfilen = "hov_cu_dns_ce2_{}.png".format(label)
fig.savefig(str(outputdir/outfilen), dpi=400)