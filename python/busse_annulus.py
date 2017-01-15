"""
2D Busse Annulus
ref: Brummell & Hart (1993) fig 5a

Usage:
    busse_annulus.py [--Ra=<Ra> --beta=<beta> --C=<C> --Pr=<Pr> --restart=<restart_file> --nx=<nx> --filter=<filter> --seed=<seed>] 

Options:
    --Ra=<Ra>      Rayleigh number [default: 39000]
    --beta=<beta>  beta [default: 2800]
    --Pr=<Pr>      Prandtl number [default: 1]
    --C=<C>        C parameter [default: 0]
    --restart=<restart_file>   Restart from checkpoint
    --nx=<nx>                  x (Fourier) resolution [default: 128]
    --filter=<filter>          fraction of modes to keep in ICs [default: 0.5]
    --seed=<seed>  random seed for ICs [default: None]
"""
import os
import sys
import numpy as np
from mpi4py import MPI
import time

from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools  import post
import logging
logger = logging.getLogger(__name__)

try:
    from dedalus.extras.checkpointing import Checkpoint
except ImportError:
    checkpointing = False
    logger.warn("No Checkpointing module.")

from docopt import docopt

# parse arguments
args = docopt(__doc__)

# Parameters
nx = int(args['--nx']) # resolution
ny = nx
Lx, Ly = (2*np.pi, 1)

Ra = float(args['--Ra'])
Pr = float(args['--Pr'])
beta = float(args['--beta'])
filter_frac = float(args['--filter'])
C = float(args['--C'])

restart = args['--restart']
seed = args['--seed']

if seed == 'None':
    seed = None
else:
    seed = int(seed)

# save data in directory named after script
data_dir = "scratch/" + sys.argv[0].split('.py')[0]
data_dir += "_ra{0:5.02e}_beta{1:5.02e}_C{2:5.02e}_Pr{3:5.02e}_filter{4:5.02e}_nx{5:d}/".format(Ra, beta, C, Pr, filter_frac,nx)

if MPI.COMM_WORLD.rank == 0:
    if not os.path.exists('{:s}/'.format(data_dir)):
        os.mkdir('{:s}/'.format(data_dir))


# Create bases and domain
start_init_time = time.time()

x_basis = de.Fourier('x', nx, interval=(0, Lx), dealias=3/2)
y_basis = de.SinCos('y', ny, interval=(0, Ly), dealias=3/2)
domain = de.Domain([x_basis, y_basis], grid_dtype=np.float64)

# 2D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['psi','theta','zeta'], time='t')
problem.meta['psi','zeta','theta']['y']['parity'] = -1 # sin basis
problem.parameters['Pr'] = Pr
problem.parameters['Ra'] = Ra
problem.parameters['beta'] = beta
problem.parameters['sbeta'] = np.sqrt(np.abs(beta))
problem.parameters['C'] = C

# construct the 2D Jacobian
problem.substitutions['J(A,B)'] = "dx(A) * dy(B) - dy(A) * dx(B)"

problem.add_equation("dt(zeta) - beta*dx(psi) + Ra/Pr * dx(theta) + C * sbeta * zeta - dx(dx(zeta)) - dy(dy(zeta)) = -J(psi,zeta)", condition="ny != 0")
problem.add_equation("dt(theta) + dx(psi) - (dx(dx(theta)) + dy(dy(theta)))/Pr = -J(psi,theta)", condition="ny != 0")
problem.add_equation("dx(dx(psi)) + dy(dy(psi)) - zeta = 0", condition="ny != 0")
problem.add_equation("zeta = 0", condition="ny ==0")
problem.add_equation("theta = 0", condition="ny ==0")
problem.add_equation("psi = 0", condition="ny ==0")

# Build solver
solver = problem.build_solver(de.timesteppers.MCNAB2)
#solver = problem.build_solver(de.timesteppers.RK222)
logger.info('Solver built')

# checkpointing
if checkpointing:
    chk = Checkpoint(data_dir)

# Initial conditions
chk_write = chk_set = 1
if restart:
    chk_write, chk_set, dt = chk.restart(restart, solver)
    counts, sets = chk.find_output_counts()

    snap_count  = counts['snapshots']
    snap_set    = sets['snapshots']
    integ_count = counts['integrals']
    integ_set   = sets['integrals']
    ts_count    = counts['timeseries']
    ts_set      = sets['timeseries']

    logger.info("Restarting from time t = {0:10.5e}".format(solver.sim_time))
else:
    y = domain.grid(1)
    theta = solver.state['theta']

    # Linear background + perturbations damped at walls
    yb, yt = y_basis.interval
    shape = domain.local_grid_shape(scales=1)
    rand = np.random.RandomState(seed)
    pert =  1e-3 * rand.standard_normal(shape) #* (yt - y) * (y - yb)
    theta['g'] = pert

    dt = 1e-4
    snap_count  = 1
    snap_set    = 1
    integ_count = 1
    integ_set   = 1
    ts_count    = 1
    ts_set      = 1

if checkpointing:
    chk.set_checkpoint(solver,wall_dt=15.,write_num=chk_write, set_num=chk_set)

# Integration parameters
solver.stop_sim_time = 20.
solver.stop_wall_time = 60.*60
solver.stop_iteration = 2000 #np.inf

# Analysis
analysis_tasks = []
snap = solver.evaluator.add_file_handler(os.path.join(data_dir,'snapshots'), sim_dt=1e-2, max_writes=10000, write_num=snap_count, set_num=snap_set)
snap.add_task("dx(psi)", name="u_y")
snap.add_task("-dy(psi)", name="u_x")
snap.add_task("zeta")
snap.add_task("theta")
snap.add_task("zeta", name="zeta_kspace", layout='c')
analysis_tasks.append(snap)

integ = solver.evaluator.add_file_handler(os.path.join(data_dir,'integrals'), sim_dt=1e-3, max_writes=10000, write_num=integ_count, set_num=integ_set)
integ.add_task("integ(dx(psi)**2, 'x')", name='<y kin en density>_x', scales=1)
integ.add_task("integ(dy(psi)**2, 'x')", name='<x kin en density>_x', scales=1)
integ.add_task("integ(zeta, 'x')", name='<vorticity>_x', scales=1)
analysis_tasks.append(integ)

timeseries = solver.evaluator.add_file_handler(os.path.join(data_dir,'timeseries'), sim_dt=1e-3, write_num=ts_count, set_num=ts_set)
timeseries.add_task("integ(integ(dx(psi)**2 + dy(psi)**2,'x'),'y')",name='Ekin')
timeseries.add_task("integ(integ(dy(psi)**2,'x'),'y')",name='E_zonal')
timeseries.add_task("2*Ra/Pr * integ(integ(psi*dx(theta),'x'),'y') - 2*integ(integ((dx(dx(zeta)) + dy(dy(zeta)))**2,'x'),'y')", name="dEdt")
analysis_tasks.append(timeseries)

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("dx(psi)**2 + dy(psi)**2", name='Ekin')

try:
    logger.info('Starting loop')
    start_run_time = time.time()

    while solver.ok:
        solver.step(dt)
        if (solver.iteration-1) % 100 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('Max E_kin = %f' %flow.max('Ekin'))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %f' %(end_run_time-start_run_time))


logger.info('beginning join operation')
if checkpointing:
    logger.info(data_dir+'/checkpoint/')
    post.merge_analysis(data_dir+'/checkpoint/')

for task in analysis_tasks:
    logger.info(task.base_path)
    post.merge_analysis(task.base_path)

