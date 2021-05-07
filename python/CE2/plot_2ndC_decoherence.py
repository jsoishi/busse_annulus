"""Plot large decoherence figure

"""
import pathlib
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8
import matplotlib.pyplot as plt
plt.ioff()

plt.style.use('prl')
from dedalus.extras import plot_tools

output = pathlib.Path('../../figs/')
scale = 4

limits = (-0.0008,0.0008)
# Layout
nrows, ncols = 2,3 
image = plot_tools.Box(2, 1)
pad = plot_tools.Frame(0.2, 0.2, 0.1, 0.1)
margin = plot_tools.Frame(0.3, 0.2, 0.1, 0.1)

# Create multifigure
mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
fig = mfig.figure
datadir = pathlib.Path('run_B_DNS_restart/data_snapshots')

filenames = [datadir.joinpath(x) for x in ['data_snapshots_s1.h5', 'data_snapshots_s2.h5', 'data_snapshots_s3.h5','data_snapshots_s4.h5','data_snapshots_s5.h5','data_snapshots_s117.h5']]

title = r'$c_{\theta\theta}(\xi, y_1, 0.5)$'
for n,fn in enumerate(filenames):
    i, j = divmod(n, ncols)
    axes = mfig.add_axes(i, j, [0, 0, 1, 1])
    with h5py.File(fn, 'r') as df:
        dset = df['tasks/interp(ctt, y1=0.500)']
        image_axes = [1, 2]
        data_slices = [0, slice(None), slice(None), 0]
        pax, cax = plot_tools.plot_bot(dset, image_axes, data_slices, axes=axes, title=title, even_scale=True,clim=limits)
        pax.set_rasterized(True)
        time = df['scales/sim_time'][0]
        pax.text(0.1,0.8,'t = {:3.2e}'.format(time),bbox={'facecolor': 'grey', 'alpha': 0.5, 'boxstyle':'round'},transform=axes.transAxes,fontsize=18)

dpi=300
for ext in ['png','pdf']:
    savename = "run_B_decoherence.{}".format(ext)
    savepath = output.joinpath(savename)
    fig.savefig(str(savepath), dpi=dpi)
