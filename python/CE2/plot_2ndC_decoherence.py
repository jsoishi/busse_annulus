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
scale = 2

limits = (-0.01,0.01)
# Layout
nrows, ncols = 1,3
image = plot_tools.Box(2, 1)
pad = plot_tools.Frame(0.05, 0.02, 0.2, 0.005)
margin = plot_tools.Frame(0.3, 0.25, 0.1, 0.05)

# Create multifigure
mfig = plot_tools.MultiFigure(nrows, ncols, image, pad, margin, scale)
fig = mfig.figure
#datadir = pathlib.Path('run_B_DNS_restart/data_snapshots')
datadir = pathlib.Path('run_A_explicit_noproj_nomean_max_knowledge_Lx_adjust/data_snapshots')

filenames = [datadir.joinpath(x) for x in ['data_snapshots_s1.h5', 'data_snapshots_s3.h5','data_snapshots_s101.h5']]

title = r'$c_{\theta\theta}$'#(\xi, y_1, y_2=0.5)$'
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
        pax.set_xlabel(r"$\xi$")
        pax.set_ylabel(r"$y_1$")
        if n != 0:
            cax.set_visible(False)
        #cax.ticklabel_format(style='sci', axis='x',scilimits=(0,0))
        pax.collections[-1].colorbar.set_ticks((-0.01,0,0.01))

dpi=300
for ext in ['png','pdf']:
    savename = "run_A_decoherence.{}".format(ext)
    savepath = output.joinpath(savename)
    fig.savefig(str(savepath), dpi=dpi)

