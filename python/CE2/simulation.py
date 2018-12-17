"""Busse Annulus CE2 script."""
import sys
sys.path = ['.',] + sys.path

import numpy as np
np.seterr(all='raise')
import time
import pathlib

from dedalus import public as de
from dedalus.tools.array import reshape_vector
import parameters as param
import diagonal
import transpose
import reverse

de.operators.parseables['Diag'] = Diag = diagonal.GridDiagonal
de.operators.parseables['Trans'] = Trans = transpose.TransposeOperator
de.operators.parseables['Rev'] = Rev = reverse.ReverseFirst

import logging
logger = logging.getLogger(__name__)

# Bases and domain
x_basis = de.Fourier('x', param.Nx, [-param.Lx/2, param.Lx/2], dealias=3/2)
y0_basis = de.SinCos('y0', param.Ny, [0, param.Ly], dealias=3/2)
y1_basis = de.SinCos('y1', param.Ny, [0, param.Ly], dealias=3/2)
domain = de.Domain([x_basis, y0_basis, y1_basis], grid_dtype=np.float64, mesh=param.mesh)

x, y0, y1 = domain.grids()

# Problem
problem = de.IVP(domain, variables=['cs','css', 'ct', 'cts', 'cst', 'ctt'])
problem.meta['cs']['x']['constant'] = True
problem.meta['ct']['x']['constant'] = True

problem.meta['cs'] ['y0']['parity'] = 1
problem.meta['ct'] ['y0']['parity'] = 1
problem.meta['cs'] ['y1']['parity'] = -1
problem.meta['ct'] ['y1']['parity'] = -1

problem.meta['css']['y0']['parity'] = -1
problem.meta['cts']['y0']['parity'] = -1
problem.meta['cst']['y0']['parity'] = -1
problem.meta['ctt']['y0']['parity'] = -1
problem.meta['css']['y1']['parity'] = -1
problem.meta['cts']['y1']['parity'] = -1
problem.meta['cst']['y1']['parity'] = -1
problem.meta['ctt']['y1']['parity'] = -1

problem.parameters['Lx'] = param.Lx
problem.parameters['Ly'] = param.Ly
problem.parameters['β'] = param.beta
problem.parameters['κ'] = param.C*np.sqrt(np.abs(param.beta))
problem.parameters['Pr'] = param.Pr
problem.parameters['Ra'] = param.Ra
problem.substitutions['D(A)'] = "Diag(interp(A, x=0), 'y0', 'y1')"

problem.substitutions['T(A)'] = "Rev(Trans(A))"

problem.substitutions['P0(A)'] = "Trans(A)"
problem.substitutions['P1(A)'] = "A" # for symmetry
problem.substitutions['L0(A)'] = "dx(dx(A)) + dy0(dy0(A))"
problem.substitutions['L1(A)'] = "dx(dx(A)) + dy1(dy1(A))"
problem.substitutions['cz'] = "dy1(dy1(cs))"
problem.substitutions['czs'] = "L0(css)"
problem.substitutions['csz'] = "L1(css)"
problem.substitutions['czz'] = "L1(czs)"
problem.substitutions['ctz'] = "L1(cts)"
problem.substitutions['czt'] = "L0(cst)"

# First stream function cumulant restrictions
problem.add_equation("cs = 0", condition="(nx != 0) or (ny0 != 0)")
# Stream function gauge
problem.add_equation("cs = 0", condition="(nx == 0) and (ny1 == 0) and (ny0 == 0)")
# First stream function cumulant evolution
problem.add_equation("dt(cz) + κ*cz - dy1(dy1(cz)) = - D(dx(dy0(csz))) - D(dx(dy1(csz)))",
                     condition="(nx == 0) and (ny0 == 0) and (ny1 != 0)")

# Second stream function cumulant restrictions
problem.add_equation("css = 0", condition="(nx == 0)")
# Second stream function cumulant evolution
problem.add_equation("dt(czz) - β*dx(csz - czs) - Ra/Pr * dx(czt - ctz) + 2*κ*czz - L0(czz) - L1(czz) = " +
                     "   dy0(P0(cs))*dx(czz) - dy0(P0(cz))*dx(csz)" +
                     " - dy1(P1(cs))*dx(czz) + dy1(P1(cz))*dx(czs)",
                     condition="(nx != 0)")

# First theta cumulant restrictions
problem.add_equation("ct = 0", condition="(nx != 0) or (ny0 != 0)")
# Theta gauge (THIS MAKES NO SENSE)
problem.add_equation("ct = 0", condition="(nx == 0) and (ny1 == 0) and (ny0 == 0)")
# First theta cumulant evolution
problem.add_equation("dt(ct) - dy1(dy1(ct))/Pr = - D(dx(dy0(cst))) - D(dx(dy1(cst)))",
                     condition="(nx == 0) and (ny0 == 0) and (ny1 != 0)")
# Second theta cumulant restrictions
problem.add_equation("cts = 0", condition="(nx == 0)")
# Second theta cumulant evolution
problem.add_equation("dt(ctz) + κ*ctz - L0(ctz)/Pr - L1(ctz) + β*dx(cts) + dx(csz) - Ra/Pr * dx(ctt) = " +
                     " (dy0(P0(cs)) - dy1(P1(cs)))*dx(ctz) - dy0(P0(ct))*dx(csz) + dy1(P1(cz))*dx(cts)",
                     condition="(nx != 0)")
# Second theta-theta cumulant restrictions
problem.add_equation("ctt = 0", condition="(nx == 0)")
# Second theta-theta cumulant evolution
problem.add_equation("dt(ctt) + dx(cst) - dx(cts) - L0(ctt)/Pr - L1(ctt)/Pr = (dy0(P0(cs)) - dy1(P1(cs))) * dx(ctt) + dy1(P1(ct))*dx(cts) - dy0(P0(ct))*dx(cst)", condition="(nx != 0)")

# symmetry equation for cst
problem.add_equation("cst = T(cts)")

# Solver
solver = problem.build_solver(de.timesteppers.RK222)
solver.stop_sim_time = param.stop_sim_time
solver.stop_wall_time = param.stop_wall_time
solver.stop_iteration = param.stop_iteration

# Initial conditions
if pathlib.Path('restart.h5').exists():
    solver.load_state('restart.h5', -1)
else:
    r2 = x**2 + (param.Ly/np.pi*np.sin((y1-y0)*np.pi/param.Ly))**2/2
    ctt = solver.state['ctt']

    ctt['g'] = param.pert_amp * np.exp(-r2/2/param.pert_width**2) * np.cos(np.pi/param.Ly *y1) * np.cos(np.pi/param.Ly *y0)
# Analysis
an1 = solver.evaluator.add_file_handler('data_checkpoints', sim_dt=param.checkpoints_sim_dt, max_writes=1)
an1.add_system(solver.state)

an2 = solver.evaluator.add_file_handler('data_snapshots', iter=param.snapshots_iter, max_writes=10)
an2.add_task("interp(czz, y1=%.3f)" %(0.0*param.Ly), scales=2)
an2.add_task("interp(czz, y1=%.3f)" %(0.1*param.Ly), scales=2)
an2.add_task("interp(czz, y1=%.3f)" %(0.25*param.Ly), scales=2)
an2.add_task("interp(ctt, y1=%.3f)" %(0.0*param.Ly), scales=2)
an2.add_task("interp(ctt, y1=%.3f)" %(0.1*param.Ly), scales=2)
an2.add_task("interp(ctt, y1=%.3f)" %(0.25*param.Ly), scales=2)

an3 = solver.evaluator.add_file_handler('data_profiles', iter=param.profiles_iter, max_writes=10)
an3.add_task("P1(cz)", name='cz')
an3.add_task("P1(cs)", name='cs')
an3.add_task("-dy1(P1(cs))", name='cu')
an3.add_task("P1(ct)", name='ct')

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

