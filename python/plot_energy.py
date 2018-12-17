"""
Plot a spacetime diagram

Usage:
    plot_spacetime.py [--var=<var>] <datapath>

Options:
    --var=<var>             variable to plot [default: <x kin en density>_x]

"""
import os
import sys
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import h5py
import numpy as np
plt.style.use('ggplot')
mpl.rcParams['font.size'] = 18

import dedalus.extras.plot_tools as dpt

from docopt import docopt

# parse arguments
args = docopt(__doc__)
var = args['--var']

datadir = args['<datapath>']
ts = h5py.File(datadir+"/timeseries/timeseries_s1.h5")

fig = plt.figure(figsize=(14,8))
ts_ax = fig.add_axes([0.1,0.1,0.85,0.85])

ts_ax.semilogy(ts['scales/sim_time'][:],ts['tasks/Ekin'][:,0], '-x')
#ts_ax.set_xlim(0.,1.)
#ts_ax.set_ylim(100,1500)
ts_ax.set_xlabel('$t$')
ts_ax.set_ylabel('$E_{kin}$')

savefilename = datadir.rstrip('/')
savefilename = savefilename.split('/')[1]
fig.savefig('../figs/{}_kin_energy.png'.format(savefilename))
