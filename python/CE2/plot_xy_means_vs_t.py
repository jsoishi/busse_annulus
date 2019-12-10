"""
Plot time slice from joint analysis files.

Usage:
    plot_xy_means_vs_t.py <file_path> [--output=<dir>]

Options:
    --output=<dir>  Output directory [default: ./img_profiles]

"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('prl')
plt.ioff()
from mpi4py import MPI

from dedalus.extras import plot_tools
from dedalus.tools.general import natural_sort
#import parameters as param


def main(file_path, output_path):
    """Save plot of specified tasks for given range of analysis writes."""

    # Plot settings
    # with h5py.File(file_path, mode='r') as file:
    #     tasks = sorted(file['tasks'].keys())
    #     tasks = tasks[MPI.COMM_WORLD.rank::MPI.COMM_WORLD.size]
    #tasks = ["P1(cz)"]
    dpi = 128
    savename = lambda task: str(output_path.joinpath('hov_{:}.png'.format(task)))

    # Selection settings
    keep_axis = A = 3
    index = 0

    # Plot Hovmoller
    image_axes = [0, A]
    image_scales = ['sim_time', 0]
    data_slices = (slice(None), 0, 0, slice(None))

    with h5py.File(file_path, 'r') as file:
        sim_time = file['scales']['sim_time'][:]
        fig = plt.figure(figsize=(16,8))
        data_slices = (slice(None), 0, 0, slice(None))
        axes = fig.add_axes([0.1, 0.1, 0.85, 0.85])

        for task_name in file['tasks']:
            dset = file['tasks'][task_name]
            # Plot average vs time
            x = sim_time[data_slices[0]]
            y = np.mean(dset[data_slices], axis=1)
            plt.semilogy(x, np.abs(y), label=task_name)
            plt.xlabel(image_scales[0])
            plt.ylabel('x-y mean')
            #plt.xlim(min(x), max(x))
        plt.legend()
        fig.savefig(output_path.joinpath("xy_means_vs_t.png"), dpi=dpi)
        plt.close(fig)


if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    # Create output directory if needed
    output_path = pathlib.Path(args['--output']).absolute()
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_path.exists():
                output_path.mkdir()

    main(args['<file_path>'], output_path)


