"""
Plot phase space, energy, and a snapshot of the vorticity (zeta)

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
datadir = sys.argv[-1]

def central_diff(t, data):
    diff = (data[2:] - data[0:-2])/(t[2:] - t[0:-2])
    return diff[:,0]

ts = h5py.File(datadir+"/timeseries/timeseries_s1.h5")
snapshots = h5py.File(datadir+"/snapshots/snapshots_s1.h5")

fig = plt.figure(figsize=(14,8))

aspect = 1/(2*np.pi)
vort_ax = fig.add_axes([0.08,0.1,0.8,0.4])
cb_ax = fig.add_axes([0.88,0.1,0.01,0.4])


phase_ax = fig.add_axes([0.1,0.6,0.35,0.35])
ts_ax = fig.add_axes([0.55,0.6,0.35,0.35])

ts_ax.plot(ts['scales/sim_time'][:],2*ts['tasks/E_zonal2'][:,0])
ts_ax.set_xlim(0.,1.)
#ts_ax.set_ylim(100,1500)
ts_ax.set_xlabel('$t$')
ts_ax.set_ylabel('$E_Z$')

dedt = central_diff(ts['scales/sim_time'][:],2*ts['tasks/E_zonal2'][:,0])

t_trim = ts['scales/sim_time'][1:-1]
ez = ts['tasks/E_zonal2'][1:-1,0]
t_select = (t_trim >0.5) & (t_trim < 0.7)
phase_ax.scatter(2*ez[t_select],2*dedt[t_select]/1000)
#phase_ax.set_xlim(750,850)
phase_ax.set_xlabel('$E_Z$')
phase_ax.set_ylabel('$10^{-3}\ dE_Z/dt$')
phase_ax.get_xaxis().get_major_formatter().set_useOffset(True)

im = vort_ax.imshow(snapshots['tasks/zeta'][-1,:,:],extent=[0,2*np.pi,0,1],aspect='auto',cmap='viridis',interpolation='none')
vort_ax.set_xlabel('x')
vort_ax.set_ylabel('y')
vort_ax.grid(b=False)
fig.colorbar(im,cax=cb_ax,label=r'$\zeta$')
savefilename = datadir.rstrip('/')
savefilename = savefilename.split('/')[1]
fig.savefig('../figs/{}_phase_energy_zeta.png'.format(savefilename))
