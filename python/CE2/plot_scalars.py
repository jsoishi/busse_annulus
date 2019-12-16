"""
Plot scalars from joint analysis files.

Usage:
    plot_scalars.py <file>

"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()


def main(filename):
    """Save plot of specified tasks for given range of analysis writes."""

    # Plot settings
    tasks = ['KE', 'EN']
    figsize = (12, 8)
    dpi = 100

    # Plot integrals
    fig = plt.figure(figsize=figsize)
    with h5py.File(filename, mode='r') as file:
        for task in tasks:
            dset = file['tasks'][task]
            time = dset.dims[0]['sim_time'][:]
            data = dset[:,:,:,0].ravel()
            data = data
            plt.semilogy(time, data, label=task)
            print('Final %s: %.8f' %(task, data[-1]))
        plt.xlabel('time')
        plt.legend(loc="upper right")
        # Save figure
        fig.savefig('integrals.png', dpi=dpi)

        ke = file['tasks/KE'][:,0,0,0]
        thing = (ke == ke.max())
        t_peak = time[ke == ke.max()]
        t_window = (time > time[1]) & (time < t_peak)
        gamma_en, log_w0 = np.polyfit(time[t_window], np.log(ke[t_window]),1)
        print("growth_rate = {}".format(gamma_en))
    fig.clear()
    plt.close(fig)


if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging

    args = docopt(__doc__)
    main(args['<file>'])

