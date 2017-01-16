import os
import sys
import matplotlib as mpl


import matplotlib.pyplot as plt
import h5py
import numpy as np
plt.style.use('ggplot')
mpl.rcParams['font.size'] = 18
datadir = sys.argv[-1]

snapshots = h5py.File(datadir+"/snapshots/snapshots_s1.h5")

fig = plt.figure(figsize=(14,14))

ax = fig.add_axes([0.09,0.1,0.8,0.4])
plt_ax = fig.add_axes([0.09,0.5,0.8,0.4])
cb_ax = fig.add_axes([0.89,0.1,0.02,0.4])
data = snapshots['tasks/zeta_kspace'][-1,:,:]
dshape = data.shape

im =ax.imshow(np.log10((data*data.conj()).real),interpolation='none',cmap='viridis',origin='lower',aspect='auto', extent=[0,dshape[1]-1,0,dshape[0]-1])
ax.grid(b=False)
fig.colorbar(im,cax=cb_ax,label=r'$\log_{10} |\hat{\zeta}|^2$')

plt_ax.semilogy(data[1,:]*data[1,:].conj())
plt_ax.set_ylabel(r'$|\hat{\zeta}|^2$')
plt_ax.axes.get_xaxis().set_visible(False)
plt_ax.set_xlim(0,dshape[1]-1)
ax.set_xlabel('m')
ax.set_ylabel('n')
savefilename = datadir.rstrip('/')
savefilename = savefilename.split('/')[1]
fig.savefig('../figs/{}_zeta_kspace.png'.format(savefilename))

