"""
Plot

Usage:
    plot_joint_hov_cu_dns_ce2.py <dns_file_path> <ce2_file_path> [--output=<dir> --label_delta=<label_delta> --label=<label> --limits=<limits>]

Options:
    --output=<dir>                 Output directory [default: ../../figs]
    --label_delta=<label_delta>    offset for DNS/CE2 labels [default: 0.05]   
    --label=<label>                label for plot 
    --limits=<limits>              plot limits
"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8
import matplotlib.pyplot as plt
plt.ioff()
from mpi4py import MPI

from dedalus.extras import plot_tools
from dedalus.tools.general import natural_sort
 
import pathlib
from docopt import docopt
from dedalus.tools.parallel import Sync

def concatenate_dns_ce2_hov(dns, ce2):
    y_ratio = dns.shape[-1]//ce2.shape[-1]
    if dns.shape[-1]%ce2.shape[-1] != 0:
        raise ValueError("DNS and CE2 y ratio is not an integer.")
    if y_ratio != 1:
        ce2_expand = np.repeat(ce2,y_ratio,axis=-1)

    return np.concatenate((dns, ce2_expand))

if __name__ == "__main__":
    args = docopt(__doc__)

    label = args['--label']
    label_delta = float(args['--label_delta'])
    limits = args['--limits']

    if limits:
        limits = limits.split(',')

    dns_filename = args['<dns_file_path>']
    ce2_filename = args['<ce2_file_path>']
    dns = h5py.File(dns_filename,"r")
    ce2 = h5py.File(ce2_filename,"r")

    time_transition = dns['scales/sim_time'][-1]
    ce2_time = ce2['scales/sim_time'][:]
    if ce2['scales/sim_time'][0] == 0:
        ce2_time += time_transition

    t_total = np.concatenate((dns['scales/sim_time'][:],ce2_time))
    xx,yy = np.meshgrid(t_total, dns['scales/y/1'][:])

    cu = concatenate_dns_ce2_hov(dns['tasks/<u>_x'][:,0,:], ce2['tasks/cu'][:,0,0,:])
    ct = concatenate_dns_ce2_hov(dns['tasks/<theta>_x'][:,0,:], ce2['tasks/ct'][:,0,0,:])

    plt.style.use('prl')
    fig = plt.figure(figsize=[12,8])
    ct_axes = fig.add_axes([0.08,0.1,0.8,0.4])
    ct_cb_ax= fig.add_axes([0.89,0.1,0.02,0.4])
    cu_axes = fig.add_axes([0.08,0.55,0.8,0.4])
    cu_cb_ax= fig.add_axes([0.89,0.55,0.02,0.4])

    img = ct_axes.pcolormesh(xx,yy,ct.T,shading='nearest', rasterized=True)#[0,3,0,1])
    fig.colorbar(img, label=r'$c_\theta$',cax=ct_cb_ax)
    ct_axes.text(time_transition-label_delta,1.02,'DNS',fontsize=14, horizontalalignment='right')
    ct_axes.text(time_transition+label_delta,1.02,'CE2',fontsize=14)
    ct_axes.annotate("",(2,1),(2,1.15),arrowprops=dict(arrowstyle='-',shrinkB=0))
    ct_axes.set_xlabel('time')
    ct_axes.set_ylabel('y')

    cu_img = cu_axes.pcolormesh(xx,yy,cu.T,shading='nearest', rasterized=True)#[0,3,0,1])
    fig.colorbar(cu_img, label=r'$c_u$',cax=cu_cb_ax)
    cu_axes.set_ylabel('y')
    cu_axes.set_xticks([])

    
    for ax in [cu_axes, ct_axes]:
        ax.set_yticks([0,0.5,1.0])

    for ext in ['png','pdf']:
        fig.savefig('../../figs/{}_hov_dns_ce2.{}'.format(label, ext))
