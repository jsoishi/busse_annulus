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
plt.style.use('prl')

def calc_growth(time, data, t_start=0, t_stop=1):
    """calculate exponential growth rate for data
    return A, gamma such that

    data_model = A*exp(gamma*t)

    """
    t_window = (time > t_start) & (time < t_stop)
    gamma, log_A = np.polyfit(time[t_window], np.log(data[t_window]),1)
    
    return gamma, np.exp(log_A)

def main(filename):
    """Save plot of specified tasks for given range of analysis writes."""

    # Plot settings
    #tasks = ['KE', 'EN']
    tasks = ['cs_integ', 'ct_integ']
    names = [r'$\left< c_\psi \right>$',r'$\left< c_\theta \right>$']
    figsize = (12, 8)
    dpi = 100

    # Plot integrals
    fig = plt.figure(figsize=figsize)
    with h5py.File(filename, mode='r') as file:
        for name,task in zip(names,tasks):
            dset = file['tasks'][task]
            time = dset.dims[0]['sim_time'][:]
            data = dset[:,:,:,0].ravel()
            data = data
            plt.semilogy(time, data**2, label=name)
            print('Final %s: %.8f' %(task, data[-1]))
        plt.xlabel('time')
        plt.ylabel('Integral')
        plt.legend(loc="lower right")
        # Save figure
        fig.savefig('integrals.png', dpi=dpi)
        fig.savefig('integrals.pdf')
        ke = file['tasks/KE'][:,0,0,0]
        #gamma, A0 = calc_growth(time, file['tasks'][tasks[0]][:,0,0,0]**2,t_start=1,t_stop=2)
        #print("growth_rate = {}".format(gamma))
    fig.clear()
    plt.close(fig)


if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging

    args = docopt(__doc__)
    main(args['<file>'])

