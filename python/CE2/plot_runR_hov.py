import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['xtick.major.pad'] = 8
#matplotlib.rcParams['ytick.major.pad'] = 8
import matplotlib.pyplot as plt
import colorcet as cc
plt.ioff()
from mpi4py import MPI

from dedalus.extras import plot_tools
from dedalus.tools.general import natural_sort

import pathlib
from dedalus.tools.parallel import Sync

plt.style.use('prl')
limits = (-400,400)
xy = False

dns_file = "../scratch/busse_annulus_ra4.00e+07_beta3.16e+05_C3.16e-01_Pr1.00e+00_filter5.00e-01_nx512_ny128_CFL/integrals.h5"
ce2_file = "run_R_dt1.25e-6/data_profiles.h5"
ce2_bias_file = "run_S_explicit_no_mean_subtr/data_profiles.h5"

dns_data = h5py.File(dns_file,'r')
ce2_data = h5py.File(ce2_file,'r')
ce2_bias_data = h5py.File(ce2_bias_file,'r')

ce2_image_axes = [0, 3]
if xy:
    dns_image_axes = [0, 2]
    dns_data_slices = (slice(None), 0, slice(None))
else:
    dns_image_axes = [0, 1]
    dns_data_slices = (slice(None), slice(None), 0)

ce2_data_slices = (slice(None), 0, 0, slice(None))

image_scales = ['sim_time', 0]

nrows, ncols = 1,3
image = plot_tools.Box(2,1)
pad = plot_tools.Frame(0.2, 0.1, 0.125, 0.125)
margin = plot_tools.Frame(0.1, 0.12, 0.1, 0.05)

y_ticks = (0.2,0.4,0.6,0.8,1.0)

scale = 2
mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
fig = mfig.figure

axes = mfig.add_axes(0, 0, [0, 0, 1, 1])
# DNS
print("dns_data shape", dns_data['tasks/<u>_x'].shape)
print('ce2 data shape', ce2_data['tasks/cu'].shape)
pax, cax = plot_tools.plot_bot(dns_data['tasks/<u>_x'], dns_image_axes, dns_data_slices, image_scales, axes=axes, title=r'$c_u$', even_scale=True, clim=limits, cmap=cc.cm['bwy'])
#pax.set_xlim(0,0.125)
pax.xaxis.set_label_text(r"$t$")
pax.yaxis.set_label_text("$y$")
pax.set_yticks(y_ticks)
pax.yaxis.label.set_rotation('horizontal')
pax.text(0.95,0.85, 'a',
         bbox={'facecolor': 'grey', 'alpha': 0.5, 'boxstyle':'round'},transform=pax.transAxes, fontsize=18)
pax.set_rasterized(True)

axes = mfig.add_axes(0, 1, [0, 0, 1, 1])
# CE2
pax, cax = plot_tools.plot_bot(ce2_data['tasks/cu'], ce2_image_axes, ce2_data_slices, image_scales, axes=axes, title=r'$\left<u\right>$ CE2', even_scale=True, clim=limits, cmap=cc.cm['bwy'])
pax.xaxis.set_label_text(r"$t$")
cax.set_visible(False)
pax.yaxis.set_label_text("$y_1$")
pax.set_yticks(y_ticks)
pax.yaxis.label.set_rotation('horizontal')

pax.text(0.95,0.85, 'b',
         bbox={'facecolor': 'grey', 'alpha': 0.5, 'boxstyle':'round'},transform=pax.transAxes, fontsize=18)
pax.set_rasterized(True)

# CE2 biased
axes = mfig.add_axes(0, 2, [0, 0, 1, 1])
# CE2
pax, cax = plot_tools.plot_bot(ce2_bias_data['tasks/cu'], ce2_image_axes, ce2_data_slices, image_scales, axes=axes, title=r'$\left<u\right>$ CE2 biased', even_scale=True, clim=limits, cmap=cc.cm['bwy'])
pax.xaxis.set_label_text(r"$t$")
cax.set_visible(False)
pax.yaxis.set_label_text("$y_1$")
pax.set_yticks(y_ticks)
pax.yaxis.label.set_rotation('horizontal')

pax.text(0.95,0.85, 'c',
         bbox={'facecolor': 'grey', 'alpha': 0.5, 'boxstyle':'round'},transform=pax.transAxes, fontsize=18)
pax.set_rasterized(True)

outputdir = pathlib.Path('../../figs/')
for ext in ['png','pdf']:
    outfilen = "run_R_S_hov_cu_dns_ce2.{}".format(ext)
    fig.savefig(str(outputdir/outfilen), dpi=400)
