"""Busse Annulus simulation parameters."""

import numpy as np

# Simulation parameters
Lx = 2 * np.pi
Ly = 1.
beta = 2.8e3
C = 0
Ra = 7.6e4
Pr = 1

# Locally correlated perturbations
# css = A * exp(-(x**2 + (2*sin((y1-y0)/2))**2/2)/δ**2)
pert_amp = 1e-3
pert_width = 0.1

# first cumulant symmetry parameters
cu_ampl = 0
cu_lambda = 1

# Discretization parameters
Nx = 32
Ny = 32
dt = 2.5e-5
stop_sim_time = 1.
stop_wall_time = np.inf
stop_iteration = np.inf
mesh = (8, 8)

# Analysis parameters
checkpoints_sim_dt = 5
snapshots_iter = 100
profiles_iter = 100
scalars_iter = 100


