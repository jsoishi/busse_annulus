"""
Plot time slice from joint analysis files.

Usage:
    plot_profiles.py <file_path> [--output=<dir>]

Options:
    --output=<dir>  Output directory [default: ./img_profiles]

"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
        for task_name in file['tasks']:
            dset = file['tasks'][task_name]
            # Plot Hovmoller
            axes = fig.add_axes([0.05, 0.4, 0.75, 0.45])
            plot_tools.plot_bot(dset, image_axes, data_slices, image_scales, even_scale=True, axes=axes, title=task_name)
            # Plot average vs time
            axes = fig.add_axes([0.05+0.03*0.75, 0.1, 0.75*0.94, 0.2])
            x = sim_time[data_slices[0]]
            y = np.mean(dset[data_slices], axis=1)
            plt.plot(x, y, '-k')
            plt.xlabel(image_scales[0])
            plt.title('<'+task_name+'>')
            plt.xlim(min(x), max(x))
            # Plot average vs height
            y = file['scales']['y1']['1.0']
            #t1 = np.max(sim_time) * 0.25
            #x1 = np.mean(dset[data_slices][sim_time < t1], axis=0)
            t1 = sim_time[data_slices[0]][0]
            x1 = dset[data_slices][0]
            #t2 = np.max(sim_time) * 0.75
            #x2 = np.mean(dset[data_slices][sim_time > t2], axis=0)
            t2 = sim_time[data_slices[0]][-1]
            x2 = dset[data_slices][-1]
            axes = fig.add_axes([0.85, 0.4, 0.1, 0.45*0.94])
            plt.plot(x1, y, '-k')
            plt.plot(x2, y, '--k')
            #plt.title(' t < %.1f (solid) \n t > %.1f (dashed)' %(t1, t2))
            plt.title('t = %.1f (solid) \n t = %.1f (dashed)' %(t1, t2))
            plt.ylabel('y')
            plt.ylim(min(y), max(y))
            plt.setp(plt.gca().get_xticklabels(), rotation=90)
            # Print double average
            fig.text(0.8, 0.2, ' <%s> \n t < %.1f = %.3e +- %.3e \n t > %.1f = %.3e +- %.3e ' %(task_name, t1, np.mean(x1), np.std(x1), t2, np.mean(x2), np.std(x2)))
            # Save figure
            #fig.suptitle('Reynolds = %.2e, Prandtl = %.2e, Schmidt = %.2e' %(param.Reynolds, param.Prandtl, param.Schmidt))
            fig.savefig(savename(task_name), dpi=dpi)
            fig.clear()
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


