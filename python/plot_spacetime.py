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

integrals = h5py.File(datadir+"/integrals/integrals_s1.h5","r")

fig = plt.figure(figsize=(14,8))
st_ax = fig.add_axes([0.1,0.1,0.8,0.8])

data = integrals['tasks/{}'.format(var)]
dpt.plot_bot_3d(data,2,0,axes=st_ax)

savefilename = datadir.rstrip('/')
savefilename = savefilename.split('/')[1]
fig.savefig('../figs/{}_var_{}_spacetime.png'.format(savefilename,var))
