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
dt = float(sys.argv[-1])
filename = sys.argv[-2]
rank = []
with open(filename,'r') as df:
    for l in df.readlines():
        match = re.search("rank: \[([ \d]+)\]", l)
        if match:
            rank.append([int(i) for i in match.group(1).split()])
rank = np.array(rank)

t_max = 100*dt*rank.shape[0]
time = 100*dt*np.arange(rank.shape[0])
print(rank.shape)
print("end time = {}, {}".format(t_max,time[-1]))

fig = plt.figure(figsize=[14,7])
im_ax = fig.add_axes([0.2/3.,0.1,0.4,0.8])
plt_ax = fig.add_axes([0.4+0.4/3.,0.1,0.4,0.8])
im_ax.imshow(rank.T, extent=[0,t_max,0,15], origin='lower', aspect='auto', interpolation='none')
im_ax.set_xlabel("time")
im_ax.set_ylabel(r"$k_\xi$")

for i in [1, 2, 3, 5, 10]:
    plt_ax.plot(time,rank[:,i], label=r'$k_\xi = {}$'.format(i))
plt_ax.legend()
plt_ax.set_xlabel("time")
plt_ax.set_ylabel("rank")


fig.savefig("rank_vs_t.png",dpi=300)
