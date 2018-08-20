"""Pointjet simulation script."""

import numpy as np
np.seterr(all='raise')
import time
import pathlib

from dedalus import public as de
from dedalus.tools.array import reshape_vector
import parameters as param
import diagonal
de.operators.parseables['Diag'] = Diag = diagonal.FourierDiagonal

import logging
logger = logging.getLogger(__name__)


# Bases and domain
x_basis = de.Fourier('x', param.Nx, [-param.Lx/2, param.Lx/2], dealias=3/2)
y0_basis = de.Fourier('y0', param.Ny, [-param.Ly/2, param.Ly/2], dealias=3/2)
y1_basis = de.Fourier('y1', param.Ny, [-param.Ly/2, param.Ly/2], dealias=3/2)
domain = de.Domain([x_basis, y0_basis, y1_basis], grid_dtype=np.float64, mesh=param.mesh)

# Reference jet
cz_ref = domain.new_field()
cz_ref.meta['x']['constant'] = True
x, y0, y1 = domain.grids()
# Build as 1D function of y0
cz_ref['g'] = -param.ref_amp * np.cos(y0*np.pi/param.Ly)**2 * np.tanh(y0/param.ref_width)
# Diagonalize
cz_ref = Diag(cz_ref, 'y0', 'y1').evaluate()

# Problem
problem = de.IVP(domain, variables=['cs','css'])
problem.meta['cs']['x']['constant'] = True
problem.parameters['Lx'] = param.Lx
problem.parameters['Ly'] = param.Ly
problem.parameters['β'] = param.beta
problem.parameters['κ'] = param.kappa
problem.parameters['ν'] = param.nu
problem.parameters['cz_ref'] = cz_ref
problem.substitutions['D(A)'] = "Diag(interp(A, x=0), 'y0', 'y1')"
problem.substitutions['P0(A)'] = "interp(A, y1='left')"
problem.substitutions['P1(A)'] = "interp(A, y0='left')"
problem.substitutions['L0(A)'] = "dx(dx(A)) + dy0(dy0(A))"
problem.substitutions['L1(A)'] = "dx(dx(A)) + dy1(dy1(A))"
problem.substitutions['cz'] = "dy0(dy0(cs))"
problem.substitutions['czs'] = "L0(css)"
problem.substitutions['csz'] = "L1(css)"
problem.substitutions['czz'] = "L1(czs)"
# First cumulant restrictions
problem.add_equation("cs = 0", condition="(nx != 0) or (ny0 != ny1)")
# Stream function gauge
problem.add_equation("cs = 0", condition="(nx == 0) and (ny0 == ny1) and (ny0 == 0)")
# First cumulant evolution
problem.add_equation("dt(cz) + κ*cz - ν*dy0(dy0(cz)) = κ*cz_ref - D(dx(dy0(csz) + dy1(csz)))",
                     condition="(nx == 0) and (ny0 == ny1) and (ny0 != 0)")
# Second cumulant restrictions
problem.add_equation("css = 0", condition="(nx == 0)")
# Second-cumulant evolution (using derived sign for β, opposite Tobias 2013)
problem.add_equation("dt(czz) + β*dx(csz - czs) + 2*κ*czz - ν*L0(czz) - ν*L1(czz) = " +
                     "   dy0(P0(cs))*dx(czz) - dy0(P0(cz))*dx(csz)" +
                     " - dy1(P1(cs))*dx(czz) + dy1(P1(cz))*dx(czs)",
                     condition="(nx != 0)")

# Solver
solver = problem.build_solver(de.timesteppers.RK222)
solver.stop_sim_time = param.stop_sim_time
solver.stop_wall_time = param.stop_wall_time
solver.stop_iteration = param.stop_iteration

# Initial conditions
if pathlib.Path('restart.h5').exists():
    solver.load_state('restart.h5', -1)
else:
    cs = solver.state['cs']
    css = solver.state['css']

    # Invert cz_ref for cs initial condition
    ic_problem = de.LBVP(domain, variables=['cs'])
    ic_problem.parameters['cz_ref'] = cz_ref
    ic_problem.add_equation("cs = 0", condition="(nx != 0) or (ny0 != ny1)")
    ic_problem.add_equation("cs = 0", condition="(nx == 0) and (ny0 == ny1) and (ny0 == 0)")
    ic_problem.add_equation("dy0(dy0(cs)) = cz_ref", condition="(nx == 0) and (ny0 == ny1) and (ny0 != 0)")
    ic_solver = ic_problem.build_solver()
    ic_solver.solve()
    cs['c'] = ic_solver.state['cs']['c']

    # Locally correlated perturbations
    r2 = x**2 + (param.Ly/np.pi*np.sin((y1-y0)*np.pi/param.Ly))**2/2
    czz_ic = domain.new_field()
    czz_ic['g'] = param.pert_amp * np.exp(-r2/2/param.pert_width**2)
    # Invert czz_ic for css initial condition
    ic_problem = de.LBVP(domain, variables=['css'])
    ic_problem.parameters['czz_ic'] = czz_ic
    ic_problem.substitutions['L0(A)'] = "dx(dx(A)) + dy0(dy0(A))"
    ic_problem.substitutions['L1(A)'] = "dx(dx(A)) + dy1(dy1(A))"
    ic_problem.add_equation("css = 0", condition="(nx == 0)")
    ic_problem.add_equation("L1(L0(css)) = czz_ic", condition="(nx != 0)")
    ic_solver = ic_problem.build_solver()
    ic_solver.solve()
    css['c'] = ic_solver.state['css']['c']

# Analysis
an1 = solver.evaluator.add_file_handler('data_checkpoints', sim_dt=param.checkpoints_sim_dt, max_writes=1)
an1.add_system(solver.state)

an2 = solver.evaluator.add_file_handler('data_snapshots', iter=param.snapshots_iter, max_writes=10)
an2.add_task("interp(czz, y1=%.3f)" %(0.0*param.Ly), scales=2)
an2.add_task("interp(czz, y1=%.3f)" %(0.1*param.Ly), scales=2)
an2.add_task("interp(czz, y1=%.3f)" %(0.2*param.Ly), scales=2)

an3 = solver.evaluator.add_file_handler('data_profiles', iter=param.profiles_iter, max_writes=10)
an3.add_task("P1(cz)", name='cz')
an3.add_task("P1(cs)", name='cs')
an3.add_task("-dy1(P1(cs))", name='cu')

an4 = solver.evaluator.add_file_handler('data_scalars', iter=param.scalars_iter, max_writes=10)
an4.add_task("-(Lx/2) * integ(P0(cz)*P0(cs) + P0(D(czs)), 'y0')", name='KE')
an4.add_task(" (Lx/2) * integ(P0(cz)*P0(cz) + P0(D(czz)), 'y0')", name='EN')

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.ok:
        dt = solver.step(param.dt)
        if (solver.iteration-1) % 10 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    run_time = end_time - start_time
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(run_time))
    logger.info('Run cost: %f cpu-hr' %(domain.dist.comm_cart.size*run_time/60/60))

