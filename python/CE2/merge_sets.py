"""
Merge analysis sets from a FileHandler.

Usage:
    merge_sets.py <joint_path> <set_paths>...

"""


import pathlib
import h5py
import numpy as np

import logging
logger = logging.getLogger(__name__.split('.')[-1])


def merge_sets(joint_path, *set_paths, trigger_name='iteration'):

    joint_path = pathlib.Path(joint_path)
    set_paths = [pathlib.Path(sp) for sp in set_paths]

    # Find trigger cadence
    with h5py.File(str(set_paths[0]), mode='r') as set_file:
        trigger = set_file['scales'][trigger_name]
        cadence = trigger[1] - trigger[0]

    # Find min and max set indices
    min_index = 0
    max_index = 0
    for set_path in set_paths:
        logger.info("Examining indices {}".format(set_path))
        with h5py.File(str(set_path), mode='r') as set_file:
            # Get set indices
            trigger = set_file['scales'][trigger_name]
            #cadence = trigger[1] - trigger[0]
            indices = (trigger // cadence).astype(int)
            min_index = min(min_index, min(indices))
            max_index = max(max_index, max(indices))

    logger.info("Creating joint file {}".format(joint_path))
    with h5py.File(str(joint_path), mode='w') as joint_file:

        # Setup file
        logger.info("Merging setup from {}".format(set_paths[0]))
        with h5py.File(str(set_paths[0]), mode='r') as set_file:
            # File metadata
            joint_file.attrs['handler_name'] = set_file.attrs['handler_name']
            joint_file.attrs['writes'] = writes = max_index + 1 - min_index
            # Copy scales
            set_file.copy('scales', joint_file)
            # Expand time scales
            for scale_name in ['sim_time', 'wall_time', 'iteration', 'write_number']:
                joint_dset = joint_file['scales'][scale_name]
                joint_dset.resize(writes, axis=0)
                joint_dset[:] = 0
            # # Copy tasks
            # set_file.copy('tasks', joint_file)
            # # Expand time axes
            # for task_name in joint_file['tasks']:
            #     joint_dset = joint_file['tasks'][task_name]
            #     joint_dset.resize(writes, axis=0)
            #     joint_dset[:] = 0
            # Tasks
            joint_tasks = joint_file.create_group('tasks')
            set_tasks = set_file['tasks']
            for task_name in set_tasks:
                # Setup dataset with automatic chunking
                set_dset = set_tasks[task_name]
                spatial_shape = set_dset.shape[1:]
                joint_shape = (writes,) + tuple(spatial_shape)
                joint_dset = joint_tasks.create_dataset(name=set_dset.name,
                                                        shape=joint_shape,
                                                        dtype=set_dset.dtype,
                                                        chunks=True)
                # Dataset metadata
                joint_dset.attrs['task_number'] = set_dset.attrs['task_number']
                joint_dset.attrs['constant'] = set_dset.attrs['constant']
                joint_dset.attrs['grid_space'] = set_dset.attrs['grid_space']
                joint_dset.attrs['scales'] = set_dset.attrs['scales']
                # Dimension scales
                for i, set_dim in enumerate(set_dset.dims):
                    joint_dset.dims[i].label = set_dim.label
                    for scale_name in set_dim:
                        scale = joint_file['scales'][scale_name]
                        joint_dset.dims.create_scale(scale, scale_name)
                        joint_dset.dims[i].attach_scale(scale)

        # Merge sets
        for set_path in set_paths:
            logger.info("Merging data from {}".format(set_path))
            with h5py.File(str(set_path), more='r') as set_file:
                # Get set indices
                trigger = set_file['scales'][trigger_name]
                #cadence = trigger[1] - trigger[0]
                indices = (trigger // cadence).astype(int)
                i0 = min(indices) - min_index
                i1 = max(indices) - min_index + 1
                # Copy scales
                for scale_name in ['sim_time', 'wall_time', 'iteration']:
                    set_dset = set_file['scales'][scale_name]
                    joint_dset = joint_file['scales'][scale_name]
                    joint_dset[i0:i1] = set_dset[:]
                joint_file['scales']['write_number'][i0:i1] = indices
                # Copy tasks
                for task_name in set_file['tasks']:
                    set_dset = set_file['tasks'][task_name]
                    joint_dset = joint_file['tasks'][task_name]
                    joint_dset[i0:i1] = set_dset[:]


if __name__ == "__main__":

    from docopt import docopt
    from mpi4py import MPI
    from dedalus.tools import logging

    args = docopt(__doc__)
    if MPI.COMM_WORLD.rank == 0:
        merge_sets(args['<joint_path>'], *args['<set_paths>'])

