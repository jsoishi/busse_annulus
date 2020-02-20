#!/usr/bin/env bash

# scalars
python3 -m dedalus merge_procs data_scalars
python3 -m dedalus merge_sets data_scalars.h5 data_scalars/*h5

# profiles
python3 -m dedalus merge_procs data_profiles
python3 -m dedalus merge_sets data_profiles.h5 data_profiles/*h5

# snapshots
python3 -m dedalus merge_procs data_snapshots

# checkpoints
python3 -m dedalus merge_procs data_checkpoints


