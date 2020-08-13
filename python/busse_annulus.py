"""
2D Busse Annulus
ref: Brummell & Hart (1993) fig 5a

Usage:
    busse_annulus.py [--Ra=<Ra> --beta=<beta> --C=<C> --Pr=<Pr> --restart=<restart_file> --nx=<nx> --ny=<ny> --filter=<filter> --seed=<seed> --ICmode=<ICmode> --stop-time=<stop_time> --use-CFL --note=<note> --Jetbias_m=<Jetbias_m> --Jetbias_a=<Jetbias_a> --symmetry --xy] 

Options:
    --Ra=<Ra>                  Rayleigh number [default: 39000]
    --beta=<beta>              beta [default: 2800]
    --Pr=<Pr>                  Prandtl number [default: 1]
    --C=<C>                    C parameter [default: 0]
    --restart=<restart_file>   Restart from checkpoint
    --nx=<nx>                  x (Fourier) resolution [default: 128]
    --ny=<ny>                  y (Sin/Cos) resolution [default: 128]
    --filter=<filter>          fraction of modes to keep in ICs [default: 0.5]
    --seed=<seed>              random seed for ICs [default: None]
    --ICmode=<ICmode>          x mode to initialize [default: None]
    --Jetbias_m=<Jetbias_m>    y mode to initialize jet mode [default: None]
    --Jetbias_a=<Jetbias_a>    amplitude for initialize jet mode [default: 1e-3]
    --stop-time=<stop_time>    simulation time to stop [default: 2.]
    --use-CFL                  use CFL condition
    --symmetry                 enforce reflect-shift symmetry in ICs
    --xy                       order bases [x,y] instead of [y,x]
    --note=<note>              a note to add to directory

"""
import glob
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

from docopt import docopt

# parse arguments
args = docopt(__doc__)

# Parameters
nx = int(args['--nx']) # resolution
ny = int(args['--ny']) # resolution
Lx, Ly = (2*np.pi, 1)

Ra = float(args['--Ra'])
Pr = float(args['--Pr'])
beta = float(args['--beta'])
filter_frac = float(args['--filter'])
C = float(args['--C'])

stop_time = float(args['--stop-time'])
restart = args['--restart']
seed = args['--seed']
ICmode = args['--ICmode']
CFL = args['--use-CFL']
note = args['--note']
Jetbias_m = args['--Jetbias_m']
Jetbias_ampl = float(args['--Jetbias_a'])
symmetry = args['--symmetry']
xy = args['--xy']
if seed == 'None':
    seed = None
else:
    seed = int(seed)

if ICmode == 'None':
    ICmode = None
else:
    ICmode = int(ICmode)

if Jetbias_m == 'None':
    Jetbias_m = None
else:
    Jetbias_m = int(Jetbias_m)

# save data in directory named after script
data_dir = "scratch/" + sys.argv[0].split('.py')[0]
data_dir += "_ra{0:5.02e}_beta{1:5.02e}_C{2:5.02e}_Pr{3:5.02e}_filter{4:5.02e}_nx{5:d}_ny{6:d}".format(Ra, beta, C, Pr, filter_frac,nx,ny)

if ICmode:
    data_dir += "_ICmode{0:d}".format(ICmode)

if CFL:
    data_dir += "_CFL"

if Jetbias_m:
    data_dir += "_Jetbias_m{}_a{}".format(Jetbias_m, Jetbias_ampl)

if symmetry:
    data_dir += "_symmetricIC"

if xy:
    data_dir += "_xy"

if note:
    data_dir += "_{}".format(note)

if restart:
    restart_dirs = glob.glob(data_dir+"restart*")
    if restart_dirs:
        restart_dirs.sort()
        last = int(re.search("_restart(\d+)", restart_dirs[-1]).group(1))
        data_dir += "_restart{}".format(last+1)
    else:
        if os.path.exists(data_dir):
            data_dir += "_restart1"


if MPI.COMM_WORLD.rank == 0:
    if not os.path.exists('{:s}/'.format(data_dir)):
        os.mkdir('{:s}/'.format(data_dir))


# Create bases and domain
start_init_time = time.time()

x_basis = de.Fourier('x', nx, interval=(0, Lx), dealias=3/2)
y_basis = de.SinCos('y', ny, interval=(0, Ly), dealias=3/2)

if xy:
    bases = [x_basis, y_basis]
else:
    bases = [y_basis, x_basis]
domain = de.Domain(bases, grid_dtype=np.float64)

# 2D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['psi','theta'], time='t')
problem.meta['psi','theta']['y']['parity'] = -1 # sin basis
problem.parameters['Pr'] = Pr
problem.parameters['Ra'] = Ra
problem.parameters['beta'] = beta
problem.parameters['sbeta'] = np.sqrt(np.abs(beta))
problem.parameters['C'] = C
problem.parameters['Lx'] = Lx
problem.parameters['Ly'] = Ly

# construct the 2D Jacobian
problem.substitutions['J(A,B)'] = "dx(A) * dy(B) - dy(A) * dx(B)"
problem.substitutions['zeta'] = "dx(dx(psi)) + dy(dy(psi))"
problem.substitutions['Avg_x(A)'] = "integ(A,'x')/Lx"
problem.substitutions['Avg_y(A)'] = "integ(A,'y')/Ly"

problem.add_equation("dt(zeta) - beta*dx(psi) + Ra/Pr * dx(theta) + C * sbeta * zeta - dx(dx(zeta)) - dy(dy(zeta)) = -J(psi,zeta)", condition="(nx != 0) or (ny != 0)")
problem.add_equation("dt(theta) + dx(psi)/Pr - (dx(dx(theta)) + dy(dy(theta)))/Pr = -J(psi,theta)")
problem.add_equation("psi = 0", condition="(nx == 0) and (ny == 0)")

# Build solver
solver = problem.build_solver(de.timesteppers.SBDF2)
#solver = problem.build_solver(de.timesteppers.MCNAB2)
#solver = problem.build_solver(de.timesteppers.RK222)
logger.info('Solver built')

# Initial conditions
if restart:
    logger.info("Restarting from time t = {0:10.5e}".format(solver.sim_time))
    solver.load_state(restart,-1)
else:
    theta = solver.state['theta']
    theta.set_scales(domain.dealias)
    x = domain.grid(axis=1, scales=domain.dealias)
    y = domain.grid(axis=0, scales=domain.dealias)

    if ICmode:
        theta['g'] = 1e-3 * np.sin(np.pi*y)*np.sin(ICmode*2*np.pi/Lx*x)
    else:
        # Linear background + perturbations damped at walls
        yb, yt = y_basis.interval
        # Random perturbations, initialized globally for same results in parallel
        gshape = domain.dist.grid_layout.global_shape(scales=domain.dealias)
        slices = domain.dist.grid_layout.slices(scales=domain.dealias)
        rand = np.random.RandomState(seed=seed)
        noise = rand.standard_normal(gshape)[slices]

        pert =  1e-3 * noise * np.sin(np.pi*y) #* (yt - y) * (y - yb)
        theta['g'] = pert
        theta.set_scales(filter_frac, keep_data=True)
        theta['c']
        theta['g']
        theta.set_scales(domain.dealias, keep_data=True)
    if Jetbias_m:
        psi = solver.state['psi']
        psi['g'] = Jetbias_ampl * np.sin(Jetbias_m*np.pi*y)
    if symmetry:
        # zero out modes with n+m != even 
        cgshape = domain.dist.coeff_layout.global_shape(scales=1)
        cslices = domain.dist.coeff_layout.slices(scales=1)

        # construct mask globally
        pre_mask = np.indices(cgshape)
        mask = (pre_mask[0,...] + pre_mask[1,...]) %2 == 0
        mask[0,:] = 0
        logger.info('mask[cslices] shape = {}'.format(mask[cslices].shape))

        theta['c'] *= mask[cslices]

if CFL:
    CFL = flow_tools.CFL(solver, initial_dt=1e-7, cadence=5, safety=0.3,
                         max_change=1.1, min_change=0.5)
    CFL.add_velocities(('dy(psi)', '-dx(psi)'))
    dt = CFL.compute_dt()
else:
    dt = 1e-4

# Integration parameters
solver.stop_sim_time = stop_time
solver.stop_wall_time = 5*24*60.*60
solver.stop_iteration = np.inf

# Analysis
analysis_tasks = []
check = solver.evaluator.add_file_handler(os.path.join(data_dir,'checkpoints'), wall_dt=3540, max_writes=50)
check.add_system(solver.state)
analysis_tasks.append(check)

snap = solver.evaluator.add_file_handler(os.path.join(data_dir,'snapshots'), sim_dt=1e-3, max_writes=2000)
snap.add_task("dx(psi)", name="u_y")
snap.add_task("-dy(psi)", name="u_x")
snap.add_task("zeta")
snap.add_task("psi")
snap.add_task("theta")
snap.add_task("zeta - Avg_x(zeta)",name="zeta_fluct")
snap.add_task("theta - Avg_x(theta)",name="theta_fluct")
snap.add_task("zeta", name="zeta_kspace", layout='c')
snap.add_task("theta", name="theta_kspace", layout='c')
analysis_tasks.append(snap)

integ = solver.evaluator.add_file_handler(os.path.join(data_dir,'integrals'), sim_dt=1e-3, max_writes=200)
integ.add_task("Avg_x(dx(psi))", name='<v>_x', scales=1)
integ.add_task("Avg_x(-dy(psi))", name='<u>_x', scales=1)
integ.add_task("Avg_x(zeta)", name='<vorticity>_x', scales=1)
integ.add_task("Avg_x(theta)", name='<theta>_x', scales=1)
integ.add_task("Avg_x((theta-Avg_x(theta)) * (zeta - Avg_x(zeta)))", name="<theta_prime zeta_prime>_x",scales=1)
analysis_tasks.append(integ)

timeseries = solver.evaluator.add_file_handler(os.path.join(data_dir,'timeseries'), sim_dt=1e-3)
timeseries.add_task("0.5*integ(dx(psi)**2 + dy(psi)**2)",name='Ekin')
timeseries.add_task("0.5*integ(Avg_x(dy(psi)**2))",name='E_zonal')
timeseries.add_task("2*Ra/Pr * integ(psi*dx(theta))", name="Ekdot_T")
timeseries.add_task("-C*sbeta*integ(dx(psi)**2 + dy(psi)**2)", name="Ekdot_drag")
timeseries.add_task("-integ((dx(dx(psi)) + dy(dy(psi)))**2)", name="Ekdot_visc")
timeseries.add_task("integ(zeta**2)", name="Enstrophy")
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
        if CFL:
            dt = CFL.compute_dt()
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %f' %(end_run_time-start_run_time))


logger.info('beginning join operation')
for task in analysis_tasks:
    logger.info(task.base_path)
    post.merge_analysis(task.base_path)

