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
            data = dset[:].ravel()
            data = data / data[0]
            plt.plot(time, data, label=task)
            print('Final %s: %.8f' %(task, data[-1]))
        plt.ylim(0.8, 1.2)
        plt.xlabel('time')
        plt.legend(loc="upper right")
        # Save figure
        fig.savefig('integrals.png', dpi=dpi)
    fig.clear()
    plt.close(fig)


if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging

    args = docopt(__doc__)
    main(args['<file>'])

