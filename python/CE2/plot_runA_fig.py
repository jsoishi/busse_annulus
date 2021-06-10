"""Plot large run A composite figure.

"""

import h5py
import numpy as np
import matplotlib
import colorcet as cc
matplotlib.use('Agg')
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8
import matplotlib.pyplot as plt
plt.ioff()

from dedalus.extras import plot_tools
from plot_joint_hov_cu_dns_ce2 import concatenate_dns_ce2_hov
plt.style.use('prl')

dns = h5py.File("../scratch/busse_annulus_ra7.60e+04_beta2.80e+03_C0.00e+00_Pr1.00e+00_filter5.00e-01_nx256_ny64_CFL_xy/integrals.h5","r")
ce2 = h5py.File("run_A_explicit_noproj_nomean_max_knowledge_Lx_adjust/data_profiles.h5","r")

dns_ts_file = h5py.File("../scratch/busse_annulus_ra7.60e+04_beta2.80e+03_C0.00e+00_Pr1.00e+00_filter5.00e-01_nx256_ny64_CFL_xy/timeseries/timeseries_s1.h5", "r")
ce2_ts_file = h5py.File("run_A_explicit_noproj_nomean_max_knowledge_Lx_adjust/data_scalars.h5","r")

aspect_ratio = 26/13.
fig_W = 10
fig_H = fig_W/aspect_ratio
fig = plt.figure(figsize=[fig_W,fig_H])

dns_ce2_cu_hov_ax = fig.add_axes([0.1,0.65,0.5,0.25])
dns_ce2_cu_hov_cb_ax = fig.add_axes([0.1,0.9,0.1,0.02])
dns_ce2_cu_y_ax = fig.add_axes([0.6,0.65,0.07,0.25])

dns_ce2_ke_ax = fig.add_axes([0.75,0.65,0.22,0.25])

dns_ce2_ct_hov_ax = fig.add_axes([0.1,0.15,0.5,0.25])
dns_ce2_ct_hov_cb_ax = fig.add_axes([0.1,0.4,0.1,0.02])
dns_ce2_ct_y_ax = fig.add_axes([0.6,0.15,0.07,0.25])




time_transition = dns['scales/sim_time'][-1]
ce2_time = ce2['scales/sim_time'][:]
if ce2['scales/sim_time'][0] == 0:
    ce2_time += time_transition

t_total = np.concatenate((dns['scales/sim_time'][:],ce2_time))
xx,yy = np.meshgrid(t_total, dns['scales/y/1'][:])

cu = concatenate_dns_ce2_hov(dns['tasks/<u>_x'][:,0,:], ce2['tasks/cu'][:,0,0,:])
ct = concatenate_dns_ce2_hov(dns['tasks/<theta>_x'][:,0,:], ce2['tasks/ct'][:,0,0,:])
cu_img = dns_ce2_cu_hov_ax.pcolormesh(xx,yy,cu.T,shading='nearest', rasterized=True, cmap=cc.cm['bwy'])#[0,3,0,1])
cb = fig.colorbar(cu_img, orientation='horizontal',cax=dns_ce2_cu_hov_cb_ax)
cb.set_label(label=r'$c_u$',fontsize=14)
dns_ce2_cu_hov_cb_ax.xaxis.tick_top()
dns_ce2_cu_hov_cb_ax.xaxis.set_label_position('top')
dns_ce2_cu_hov_cb_ax.tick_params(axis='x', which='major',labelsize=10, pad=2)
dns_ce2_cu_hov_ax.axvline(2.0,color='k',alpha=0.4,linewidth=1.)

dns_t_avg_start = 1000
dns_cu_tavg = dns['tasks/<u>_x'][dns_t_avg_start:,0,:].mean(axis=0)
ce2_cu_tavg = ce2['tasks/cu'][:,0,0,:].mean(axis=0)
dns_ce2_cu_y_ax.plot(dns_cu_tavg, dns['scales/y/1'][:], label='DNS')
dns_ce2_cu_y_ax.plot(ce2_cu_tavg, ce2['scales/y0/1.0'][:], label='CE2')
dns_ce2_cu_y_ax.set_ylim(0,1)
dns_ce2_cu_y_ax.set_xlabel(r"$c_u$")


dns_ce2_cu_hov_ax.set_xticks([0,0.5,1.,1.5,2., 2.5])
dns_ce2_cu_hov_ax.set_xlabel("t")
dns_ce2_cu_hov_ax.set_ylabel("y")

# Kinetic energy

#ke_total_time = np.concatenate(dns_ts_file['scales/sim_time'][:], time_transition+ce2_ts_file['scales/sim_time'][:])
#ke_total = np.concatenate(dns_ts_file['scales/sim_time'][:], time_transition+ce2_ts_file['scales/sim_time'][:])
dns_ce2_ke_ax.plot(dns_ts_file['scales/sim_time'][:], dns_ts_file['tasks/Ekin'][:,0],linewidth=1)
dns_ce2_ke_ax.plot(ce2_ts_file['scales/sim_time'][:]+time_transition, ce2_ts_file['tasks/KE'][:,0,0]/(2*np.pi),linewidth=1)
dns_ce2_ke_ax.set_xlabel("t")
dns_ce2_ke_ax.axvline(2.0,color='k',alpha=0.4,linewidth=1.)
dns_ce2_ke_ax.ticklabel_format(scilimits=[-3, 3])
dns_ce2_ke_ax.set_ylabel("Kinetic energy",fontsize=14)


dns_ct_tavg = dns['tasks/<theta>_x'][dns_t_avg_start:,0,:].mean(axis=0)
ce2_ct_tavg = ce2['tasks/ct'][:,0,0,:].mean(axis=0)

dns_ce2_cu_y_ax.get_yaxis().set_visible(False)
ct_img = dns_ce2_ct_hov_ax.pcolormesh(xx,yy,ct.T,shading='nearest', rasterized=True, cmap='RdBu_r')#[0,3,0,1])
cb = fig.colorbar(ct_img, orientation='horizontal',cax=dns_ce2_ct_hov_cb_ax)
cb.set_label(label=r'$c_{\theta}$',fontsize=14)
dns_ce2_ct_hov_ax.set_xticks([0,0.5,1.,1.5,2., 2.5])
dns_ce2_ct_hov_cb_ax.xaxis.tick_top()
dns_ce2_ct_hov_cb_ax.xaxis.set_label_position('top')
dns_ce2_ct_hov_cb_ax.tick_params(axis='x', which='major',labelsize=10, pad=2)
dns_ce2_ct_hov_ax.set_xlabel("t")
dns_ce2_ct_hov_ax.set_ylabel("y")
dns_ce2_ct_hov_ax.axvline(2.0,color='k',alpha=0.4,linewidth=1.)

dns_ce2_ct_y_ax.plot(dns_ct_tavg, dns['scales/y/1'][:], label='DNS')
dns_ce2_ct_y_ax.plot(ce2_ct_tavg, ce2['scales/y0/1.0'][:], label='CE2')
dns_ce2_ct_y_ax.get_yaxis().set_visible(False)
dns_ce2_ct_y_ax.set_ylim(0,1)
dns_ce2_ct_y_ax.set_xlim(-.15,.15)
dns_ce2_ct_y_ax.set_xticks([-.1,.1])
dns_ce2_ct_y_ax.set_xlabel(r"$c_{\theta}$")
dns_ce2_ct_y_ax.legend(bbox_to_anchor=(2.0,1.0))


fig.savefig("../../figs/run_A_fig.png",dpi=300)
fig.savefig("../../figs/run_A_fig.pdf",dpi=300)
