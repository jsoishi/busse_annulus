import numpy as np
import matplotlib.pyplot as plt
import re
import sys

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['xtick.major.pad'] = 8
matplotlib.rcParams['ytick.major.pad'] = 8
import matplotlib.pyplot as plt
plt.ioff()

plt.style.use('prl')

niter = 10
dt = float(sys.argv[-1])
filename = sys.argv[-2]
name = sys.argv[-3]
kx_to_plot = [1,2,3,6,8]
rank = []
with open(filename,'r') as df:
    for l in df.readlines():
        match = re.search("rank: \[([ \d]+)\]", l)
        if match:
            rank.append([int(i) for i in match.group(1).split()])
rank = np.array(rank)

t_max = niter*dt*rank.shape[0]
time = niter*dt*np.arange(rank.shape[0])
print(rank.shape)
print("end time = {}, {}".format(t_max,time[-1]))

fig = plt.figure(figsize=[16,8])
im_ax = fig.add_axes([0.2/3.,0.1,0.4,0.8])
cb_ax = fig.add_axes([0.2/3.,0.9,0.4,0.02])
plt_ax = fig.add_axes([0.4+0.4/3.,0.1,0.4,0.8])


nrank = 15 #rank.max()+1
cmap = plt.get_cmap('viridis', nrank)
im = im_ax.imshow(rank.T, cmap=cmap, extent=[0,t_max,0,15], origin='lower', aspect='auto', interpolation='none', vmin=0, vmax=nrank-1)
cbar = plt.colorbar(mappable=im, cax=cb_ax, orientation='horizontal')
tick_locs = (np.arange(nrank) + 0.5)*(nrank-1)/nrank
cbar.set_ticks(tick_locs)

# set tick labels (as before)
cbar.set_ticklabels(np.arange(nrank))

cb_ax.set_xlabel("rank")
cb_ax.xaxis.set_ticks_position('top')
cb_ax.xaxis.set_label_position('top')
im_ax.set_xlabel("time")
im_ax.set_ylabel(r"$k_\xi$")

for i in kx_to_plot:
    plt_ax.plot(time,rank[:,i], label=r'$k_\xi = {}$'.format(i))

plt_ax.axhline(3, alpha=0.5, color='k')
plt_ax.set_xlim(0,1)
plt_ax.set_ylim(0,15)
plt_ax.legend()
plt_ax.set_xlabel("time")
plt_ax.set_ylabel("rank")

for ext in ['png', 'pdf']:
    fig.savefig("rank_vs_t_{}.{}".format(name, ext), dpi=300)

