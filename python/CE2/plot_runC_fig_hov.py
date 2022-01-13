"""Plot large run C composite figure.

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

dns_dir = "../scratch/busse_annulus_ra8.00e+07_beta5.00e+05_C0.00e+00_Pr1.00e+00_filter5.00e-01_nx512_ny256_CFL_xy/"
ce2_dir = "run_C_DNS_burst_Lx_adjust/"
ce2_q_dir = "run_C_DNS_quiet/"
ce2_mi_dir = "run_C/" # maximum ignorance run
dns = h5py.File(dns_dir+"integrals.h5","r")
ce2_b = h5py.File(ce2_dir+"data_profiles.h5","r")
ce2_q = h5py.File(ce2_q_dir+"data_profiles.h5","r")
ce2_mi = h5py.File(ce2_mi_dir+"data_profiles.h5","r")
dns_ts_file = h5py.File(dns_dir+"timeseries/timeseries_s1.h5", "r")
ce2_b_ts_file = h5py.File(ce2_dir+"data_scalars.h5","r")
ce2_q_ts_file = h5py.File(ce2_q_dir+"data_scalars.h5","r")
ce2_mi_ts_file = h5py.File(ce2_mi_dir+"data_scalars.h5","r")

aspect_ratio = 26/(13.*1.5)
fig_W = 10
fig_H = fig_W/aspect_ratio
fig = plt.figure(figsize=[fig_W,fig_H])

# in subfig fractions
subfig_h = 0.5   #height
subfig_hm = 0.25 #h margin
subfig_h_cb = 0.04
subfig_wm = 0.15
subfig_wm_r = 0.1
subfig_w_hov = 0.5
subfig_w_cb = 0.1
subfig_w_ts = 0.22

npanels = 3
subfig_h /= npanels
subfig_hm /= npanels
subfig_h_cb /= npanels
# top to bottom

dns_ce2_cu_q_hov_ax    = fig.add_axes([subfig_wm,(2*npanels-1)*subfig_hm+(npanels-1)*subfig_h,subfig_w_hov,subfig_h])
dns_ce2_cu_q_hov_cb_ax = fig.add_axes([subfig_wm,(2*npanels-1)*subfig_hm+npanels*subfig_h,subfig_w_cb,subfig_h_cb])
dns_ce2_ke_q_ax = fig.add_axes([subfig_wm+subfig_w_hov+subfig_wm_r,(2*npanels-1)*subfig_hm+(npanels-1)*subfig_h,subfig_w_ts,subfig_h])

dns_ce2_cu_b_hov_ax = fig.add_axes([subfig_wm,(2*(npanels-1) -1)*subfig_hm+(npanels-2)*subfig_h,subfig_w_hov,subfig_h])
dns_ce2_cu_b_hov_cb_ax = fig.add_axes([subfig_wm,(2*(npanels-1) -1)*subfig_hm+(npanels-1)*subfig_h,0.1,subfig_h_cb])
dns_ce2_ke_b_ax = fig.add_axes([subfig_wm+subfig_w_hov+subfig_wm_r,(2*(npanels-1) -1)*subfig_hm+(npanels-2)*subfig_h,subfig_w_ts,subfig_h])

dns_ce2_cu_mi_hov_ax = fig.add_axes([subfig_wm,subfig_hm,subfig_w_hov,subfig_h])
dns_ce2_cu_mi_hov_cb_ax = fig.add_axes([subfig_wm,subfig_hm+subfig_h,0.1,subfig_h_cb])
dns_ce2_ke_mi_ax = fig.add_axes([subfig_wm+subfig_w_hov+subfig_wm_r,subfig_hm,subfig_w_ts,subfig_h])

cu_vmin = -1200
cu_vmax = 1200
dns_start_index = 1500
trans_index = 1666
quiet_index = 1690
burst_time_transition = dns['scales/sim_time'][trans_index]
quiet_time_transition = dns['scales/sim_time'][quiet_index]
ce2_b_time = ce2_b['scales/sim_time'][:]
if ce2_b['scales/sim_time'][0] == 0:
    ce2_b_time += burst_time_transition
ce2_q_time = ce2_q['scales/sim_time'][:] + quiet_time_transition
t_total_b = ce2_b_time
t_total_q = dns['scales/sim_time'][dns_start_index:quiet_index]
xx_b,yy_b = np.meshgrid(ce2_b_time, ce2_b['scales/y0/1.0'][:])
xx_q,yy_q = np.meshgrid(t_total_q, dns['scales/y/1'][:])
xx_mi,yy_mi = np.meshgrid(ce2_mi['scales/sim_time'][:], ce2_mi['scales/y0/1.0'][:])

#cu_q = concatenate_dns_ce2_hov(dns['tasks/<u>_x'][dns_start_index:quiet_index,0,:], ce2_q['tasks/cu'][:,0,0,:])
cu_b = ce2_b['tasks/cu'][:,0,0,:]
cu_mi = ce2_mi['tasks/cu'][:,0,0,:]

# burst 
cu_b_img = dns_ce2_cu_b_hov_ax.pcolormesh(xx_b,yy_b,cu_b.T,shading='nearest', rasterized=True,cmap=cc.cm['bwy'],vmin=cu_vmin, vmax=cu_vmax)
cb = fig.colorbar(cu_b_img, orientation='horizontal',cax=dns_ce2_cu_b_hov_cb_ax)
cb.set_label(label=r'$c_u$',fontsize=14)
dns_ce2_cu_b_hov_cb_ax.xaxis.tick_top()
dns_ce2_cu_b_hov_cb_ax.xaxis.set_label_position('top')
dns_ce2_cu_b_hov_cb_ax.tick_params(axis='x', which='major',labelsize=10, pad=2)
#dns_ce2_cu_b_hov_ax.axvline(burst_time_transition,color='k',alpha=0.4,linewidth=1.)

dns_t_avg_start = 1000
dns_cu_tavg = dns['tasks/<u>_x'][dns_t_avg_start:,0,:].mean(axis=0)
ce2_cu_b_tavg = ce2_b['tasks/cu'][:,0,0,:].mean(axis=0)
ce2_cu_q_tavg = ce2_q['tasks/cu'][:,0,0,:].mean(axis=0)
# dns_ce2_cu_b_y_ax.plot(dns_cu_tavg, dns['scales/y/1'][:], label='DNS')
# dns_ce2_cu_b_y_ax.plot(ce2_cu_b_tavg, ce2_b['scales/y0/1.0'][:], label='CE2')
# dns_ce2_cu_b_y_ax.set_ylim(0,1)
# dns_ce2_cu_b_y_ax.set_xlabel(r"$c_u$")

dns_ce2_cu_b_hov_ax.set_xticks([burst_time_transition,1.75,1.85])
dns_ce2_cu_b_hov_ax.set_xlabel("t")
dns_ce2_cu_b_hov_ax.set_ylabel("y")

# Kinetic energy
Lx = 2*np.pi
ce2_b_ts_time   = ce2_b_ts_file['scales/sim_time'][:]+burst_time_transition
print("b KE duration {}".format(ce2_b_ts_time[-1] - dns_ts_file['scales/sim_time'][trans_index]))
ce2_b_ts_ekin   = ce2_b_ts_file['tasks/KE'][:,0,0,0]/Lx
ce2_b_ts_ezonal = ce2_b_ts_file['tasks/KE_mean'][:,0,0,0]/Lx
ce2_b_ts_nz     = ce2_b_ts_ekin - ce2_b_ts_ezonal

ce2_q_ts_time = ce2_q_ts_file['scales/sim_time'][:]+quiet_time_transition
print("q KE duration {}".format(ce2_q_ts_time[-1] - dns_ts_file['scales/sim_time'][quiet_index]))
ce2_q_ts_ekin   = ce2_q_ts_file['tasks/KE'][:,0,0,0]/Lx
ce2_q_ts_ezonal = ce2_q_ts_file['tasks/KE_mean'][:,0,0,0]/Lx
ce2_q_ts_nz     = ce2_q_ts_ekin - ce2_q_ts_ezonal

dns_ts_ekin = dns_ts_file['tasks/Ekin'][:,0]
dns_ts_ezonal = dns_ts_file['tasks/E_zonal'][:,0]
dns_ts_nz = dns_ts_ekin - dns_ts_ezonal

ce2_mi_ts_time = ce2_mi_ts_file['scales/sim_time'][:]
ce2_mi_ts_ekin   = ce2_mi_ts_file['tasks/KE'][:,0,0,0]/Lx
ce2_mi_ts_ezonal = ce2_mi_ts_file['tasks/KE_mean'][:,0,0,0]/Lx
ce2_mi_ts_nz     = ce2_mi_ts_ekin - ce2_mi_ts_ezonal

#ke_total_time = np.concatenate(dns_ts_file['scales/sim_time'][:], time_transition+ce2_ts_file['scales/sim_time'][:])
#ke_total = np.concatenate(dns_ts_file['scales/sim_time'][:], time_transition+ce2_ts_file['scales/sim_time'][:])
#dns_ce2_ke_b_ax.plot(dns_ts_file['scales/sim_time'][:trans_index],dns_ts_ekin[:trans_index],linewidth=1)
#dns_ce2_ke_b_ax.plot(dns_ts_file['scales/sim_time'][:trans_index],dns_ts_ezonal[:trans_index],linewidth=1,alpha=0.4)
#dns_ce2_ke_b_ax.plot(dns_ts_file['scales/sim_time'][:trans_index],dns_ts_nz[:trans_index],linewidth=1)
#dns_ce2_ke_b_ax.set_prop_cycle(None)
# dns_ce2_ke_b_ax.plot(ce2_b_ts_time, ce2_b_ts_ekin,linewidth=1,ls=':')
# dns_ce2_ke_b_ax.plot(ce2_b_ts_time, ce2_b_ts_ezonal,linewidth=1,alpha=0.4,ls=':')
# dns_ce2_ke_b_ax.plot(ce2_b_ts_time, ce2_b_ts_nz,linewidth=1,ls=':')
dns_ce2_ke_b_ax.plot(ce2_b_ts_time, ce2_b_ts_ekin,linewidth=1)
dns_ce2_ke_b_ax.plot(ce2_b_ts_time, ce2_b_ts_ezonal,linewidth=1,alpha=0.4)
dns_ce2_ke_b_ax.plot(ce2_b_ts_time, ce2_b_ts_nz,linewidth=1)

#dns_ce2_ke_b_ax.plot(ce2_q_ts_time, ce2_q_ts_file['tasks/KE'][:,0,0]/(2*np.pi),linewidth=1)
dns_ce2_ke_b_ax.set_yscale('log')
dns_ce2_ke_b_ax.set_xlabel("t")
dns_ce2_ke_b_ax.set_xlim(burst_time_transition,ce2_b_ts_time[-1])
dns_ce2_ke_b_ax.set_ylim(3e3,1e7)
dns_ce2_ke_b_ax.set_yticks([1e4,1e5,1e6])
#dns_ce2_ke_b_ax.axvline(burst_time_transition,color='k',alpha=0.4,linewidth=1.)
#dns_ce2_ke_b_ax.ticklabel_format(scilimits=[-3, 3])
dns_ce2_ke_b_ax.set_ylabel("Kinetic energy",fontsize=14)


# DNS
#dns_ce2_cu_b_y_ax.get_yaxis().set_visible(False)
cu_q_img = dns_ce2_cu_q_hov_ax.pcolormesh(xx_q,yy_q,dns['tasks/<u>_x'][dns_start_index:quiet_index,0,:].T,shading='nearest', rasterized=True, cmap=cc.cm['bwy'], vmin=cu_vmin, vmax=cu_vmax)#[0,3,0,1])
cb = fig.colorbar(cu_q_img, orientation='horizontal',cax=dns_ce2_cu_q_hov_cb_ax)
cb.set_label(label=r'$c_{u}$',fontsize=14)
#dns_ce2_cu_b_hov_ax.set_xticks([1.5,1.6,1.7,1.8])
dns_ce2_cu_q_hov_cb_ax.xaxis.tick_top()
dns_ce2_cu_q_hov_cb_ax.xaxis.set_label_position('top')
dns_ce2_cu_q_hov_cb_ax.tick_params(axis='x', which='major',labelsize=10, pad=2)
dns_ce2_cu_q_hov_ax.set_xlabel("t")
dns_ce2_cu_q_hov_ax.set_ylabel("y")

# dns_ce2_cu_q_y_ax.plot(dns_cu_tavg, dns['scales/y/1'][:], label='DNS')
# dns_ce2_cu_q_y_ax.plot(ce2_cu_q_tavg, ce2_q['scales/y0/1.0'][:], label='CE2')
# dns_ce2_cu_q_y_ax.get_yaxis().set_visible(False)
# dns_ce2_cu_q_y_ax.set_ylim(0,1)
# #dns_ce2_cu_q_y_ax.set_xlim(-.15,.15)
# #dns_ce2_cu_q_y_ax.set_xticks([-.1,.1])
# dns_ce2_cu_q_y_ax.set_xlabel(r"$c_{u}$")
# #dns_ce2_cu_b_y_ax.legend(bbox_to_anchor=(2.0,1.0))

dns_ce2_ke_q_ax.plot(dns_ts_file['scales/sim_time'][:quiet_index],dns_ts_ekin[:quiet_index],linewidth=1)
dns_ce2_ke_q_ax.plot(dns_ts_file['scales/sim_time'][:quiet_index],dns_ts_ezonal[:quiet_index],linewidth=1,alpha=0.4)
dns_ce2_ke_q_ax.plot(dns_ts_file['scales/sim_time'][:quiet_index],dns_ts_nz[:quiet_index],linewidth=1)
#dns_ce2_ke_q_ax.set_prop_cycle(None)

#dns_ce2_ke_q_ax.plot(ce2_q_ts_time, ce2_q_ts_ekin,linewidth=1,ls=':')
#dns_ce2_ke_q_ax.plot(ce2_q_ts_time, ce2_q_ts_ezonal,linewidth=1,alpha=0.4,ls=':')
#dns_ce2_ke_q_ax.plot(ce2_q_ts_time, ce2_q_ts_nz,linewidth=1,ls=':')
dns_ce2_ke_q_ax.set_yscale('log')
dns_ce2_ke_q_ax.set_xlabel("t")
dns_ce2_ke_q_ax.set_ylabel("Kinetic energy",fontsize=14)
dns_ce2_ke_q_ax.set_xlim(1.5,dns_ts_file['scales/sim_time'][quiet_index])
dns_ce2_ke_q_ax.set_ylim(3e3,1e7)
dns_ce2_ke_q_ax.set_yticks([1e4,1e5,1e6])
#dns_ce2_ke_q_ax.axvline(quiet_time_transition,color='k',alpha=0.4,linewidth=1.)


# maximum ignorance
cu_mi_img = dns_ce2_cu_mi_hov_ax.pcolormesh(xx_mi,yy_mi,cu_mi.T,shading='nearest', rasterized=True, cmap=cc.cm['bwy'], vmin=cu_vmin, vmax=cu_vmax)#[0,3,0,1])
cb = fig.colorbar(cu_mi_img, orientation='horizontal',cax=dns_ce2_cu_mi_hov_cb_ax)
cb.set_label(label=r'$c_{u}$',fontsize=14)
#dns_ce2_cu_b_hov_ax.set_xticks([1.5,1.6,1.7,1.8])
dns_ce2_cu_mi_hov_cb_ax.xaxis.tick_top()
dns_ce2_cu_mi_hov_cb_ax.xaxis.set_label_position('top')
dns_ce2_cu_mi_hov_cb_ax.tick_params(axis='x', which='major',labelsize=10, pad=2)
dns_ce2_cu_mi_hov_ax.set_xlabel("t")
dns_ce2_cu_mi_hov_ax.set_ylabel("y")

ke_t_select = (ce2_mi_ts_time < 0.35) & (ce2_mi_ts_time > 0.3)
dns_ce2_ke_mi_ax.plot(ce2_mi_ts_time[ke_t_select], ce2_mi_ts_ekin[ke_t_select],linewidth=1)
dns_ce2_ke_mi_ax.plot(ce2_mi_ts_time[ke_t_select], ce2_mi_ts_ezonal[ke_t_select],linewidth=1,alpha=0.4)
dns_ce2_ke_mi_ax.plot(ce2_mi_ts_time[ke_t_select], ce2_mi_ts_nz[ke_t_select],linewidth=1)
dns_ce2_ke_mi_ax.set_yscale('log')
dns_ce2_ke_mi_ax.set_xlabel("t")
dns_ce2_ke_mi_ax.set_ylabel("Kinetic energy",fontsize=14)
dns_ce2_ke_mi_ax.set_xlim(ce2_mi_ts_time[ke_t_select][0],ce2_mi_ts_time[ke_t_select][-1])
dns_ce2_ke_mi_ax.set_ylim(3e3,1e7)
dns_ce2_ke_mi_ax.set_yticks([1e4,1e5,1e6])
dns_ce2_ke_mi_ax.set_xticks([0.3,0.35])

fig.savefig("../../figs/run_C_fig.png",dpi=300)
fig.savefig("../../figs/run_C_fig.pdf",dpi=300)
