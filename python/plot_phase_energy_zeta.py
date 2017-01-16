"""
Plot phase space, energy, and a snapshot of the vorticity (zeta)

"""
import os
import sys
import matplotlib as mpl


import matplotlib.pyplot as plt
import h5py
import numpy as np
plt.style.use('ggplot')
mpl.rcParams['font.size'] = 18
datadir = sys.argv[-1]

ts = h5py.File(datadir+"/timeseries/timeseries_s1.h5")
snapshots = h5py.File(datadir+"/snapshots/snapshots_s1.h5")

fig = plt.figure(figsize=(14,8))

aspect = 1/(2*np.pi)
vort_ax = fig.add_axes([0.08,0.1,0.8,0.4])
cb_ax = fig.add_axes([0.88,0.1,0.01,0.4])


phase_ax = fig.add_axes([0.08,0.6,0.35,0.35])
ts_ax = fig.add_axes([0.55,0.6,0.35,0.35])

ts_ax.plot(ts['scales/sim_time'][:],ts['tasks/E_zonal'][:,0])
ts_ax.set_xlabel('time')
ts_ax.set_ylabel('Zonal Kinetic Energy')

phase_ax.plot(ts['scales/sim_time'][:],ts['tasks/dEdt'][:,0])
phase_ax.set_xlabel('time')
phase_ax.set_ylabel('dE/dt')

im = vort_ax.imshow(snapshots['tasks/zeta'][-1,:,:],extent=[0,2*np.pi,0,1],aspect='auto',cmap='viridis',interpolation='none')
vort_ax.set_xlabel('x')
vort_ax.set_ylabel('y')
vort_ax.grid(b=False)
fig.colorbar(im,cax=cb_ax,label=r'$\zeta$')
savefilename = datadir.rstrip('/')
savefilename = savefilename.split('/')[1]
fig.savefig('../figs/{}_phase_energy_zeta.png'.format(savefilename))
