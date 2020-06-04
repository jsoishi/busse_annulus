"""Busse Annulus CE2 script."""
import sys
sys.path = ['.',] + sys.path

import numpy as np
#np.seterr(all='raise')
import time
import pathlib

from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools.array import reshape_vector
import parameters as param
import diagonal
import transpose
import reverse
from symmetry import enforce_symmetry
import projection

de.operators.parseables['Diag'] = Diag = diagonal.GridDiagonal
de.operators.parseables['Trans'] = Trans = transpose.TransposeOperator
de.operators.parseables['Rev'] = Rev = reverse.ReverseFirst

import logging
logger = logging.getLogger(__name__)

if not hasattr(param, "force_symmetry"):
    param.force_symmetry = 0

if not hasattr(param, "use_czt"):
    param.use_czt = False

if not hasattr(param, "new_IC"):
    param.new_IC = False

if not hasattr(param, "cu_k"):
    param.cu_k = 1
logger.info("Running with Nx = {:d}, Ny = {:d}".format(param.Nx, param.Ny))
logger.info("Ra = {:e}".format(param.Ra))
logger.info("beta = {:e}".format(param.beta))
logger.info("Pr = {:e}".format(param.Pr))
logger.info("C = {:e}".format(param.C))
logger.info("cu_lambda = {:e}".format(param.cu_lambda))
logger.info("cu_ampl = {:e}".format(param.cu_ampl))
if param.force_symmetry:
    logger.info("enforcing symmetry every {:d} timesteps".format(param.force_symmetry))

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

if param.diagonal_tendancy:
    logger.info("Adding {:e} to dt(czz) nad dt(ctt)".format(param.diagonal_tendancy))
    gamma_field = domain.new_field()
    gamma_field.meta['y0']['parity'] = -1
    gamma_field.meta['y1']['parity'] = -1

    diag_index = (y0 == y1) & (x == 0)
    gamma_field['g'][diag_index] = param.diagonal_tendancy
    problem.parameters['Gamma'] = gamma_field
else:
    gamma_field = domain.new_field()
    gamma_field.meta['y0']['parity'] = -1
    gamma_field.meta['y1']['parity'] = -1
    gamma_field['g'] = 0
    problem.parameters['Gamma'] = gamma_field

# First stream function cumulant restrictions
problem.add_equation("cs = 0", condition="(nx != 0) or (ny0 != 0)")
# Stream function gauge
problem.add_equation("cs = 0", condition="(nx == 0) and (ny1 == 0) and (ny0 == 0)")
# First stream function cumulant evolution
problem.add_equation("dt(cz) + κ*cz - dy1(dy1(cz)) = - D(dx(dy0(csz))) - D(dx(dy1(csz)))",
                     condition="(nx == 0) and (ny0 == 0) and (ny1 != 0)")

# Second stream function cumulant restrictions
problem.add_equation("css = 0", condition="(nx == 0)")
problem.add_equation("css = 0", condition="(nx != 0) and ((ny0 == 0) or (ny1 == 0))")
# Second stream function cumulant evolution
problem.add_equation("dt(czz) + 2*κ*czz - L0(czz) - L1(czz) = " +
                     "   β*dx(csz - czs) + Ra/Pr * dx(czt - ctz)" +
                     " + dy0(P0(cs))*dx(czz) - dy0(P0(cz))*dx(csz)" +
                     " - dy1(P1(cs))*dx(czz) + dy1(P1(cz))*dx(czs) + Gamma",
                     condition="(nx != 0 and ny0 !=0 and ny1 != 0)")

# First theta cumulant restrictions
problem.add_equation("ct = 0", condition="(nx != 0) or (ny0 != 0)")
# Theta gauge (THIS MAKES NO SENSE)
problem.add_equation("ct = 0", condition="(nx == 0) and (ny1 == 0) and (ny0 == 0)")
# First theta cumulant evolution
problem.add_equation("dt(ct) - dy1(dy1(ct))/Pr = - D(dx(dy0(cst))) - D(dx(dy1(cst)))",
                     condition="(nx == 0) and (ny0 == 0) and (ny1 != 0)")
# Second theta cumulant restrictions
problem.add_equation("cts = 0", condition="(nx == 0)")
problem.add_equation("cts = 0", condition="(nx != 0) and ((ny0 == 0) or (ny1 == 0))")
# Second theta cumulant evolution
problem.add_equation("dt(ctz) + κ*ctz - L0(ctz)/Pr - L1(ctz) = " +
                     "- β*dx(cts) - dx(csz) + Ra/Pr * dx(ctt)" +
                     "+ (dy0(P0(cs)) - dy1(P1(cs)))*dx(ctz) - dy0(P0(ct))*dx(csz) + dy1(P1(cz))*dx(cts)",
                     condition="(nx != 0 and ny0 !=0 and ny1 != 0)")
# Second theta-theta cumulant restrictions
problem.add_equation("ctt = 0", condition="(nx == 0)")
problem.add_equation("ctt = 0", condition="(nx != 0) and ((ny0 == 0) or (ny1 == 0))")
# Second theta-theta cumulant evolution
problem.add_equation("dt(ctt) - L0(ctt)/Pr - L1(ctt)/Pr =  - dx(cst) + dx(cts) + (dy0(P0(cs)) - dy1(P1(cs))) * dx(ctt) + dy1(P1(ct))*dx(cts) - dy0(P0(ct))*dx(cst) + Gamma",
                     condition="(nx != 0 and ny0 !=0 and ny1 != 0)")

if param.use_czt:
    logger.info("Using explicit czt equation")
    problem.add_equation("dt(czt) + κ*czt - L0(czt) - L1(czt)/Pr = " +
                         "β*dx(cst) + dx(czs) - Ra/Pr * dx(ctt) " +
                         "+ (dy0(P0(cs)) - dy1(P1(cs)))*dx(czt) + dy1(P1(ct))*dx(czs) - dy0(P0(cz))*dx(cst)",
                     condition="(nx != 0 and ny0 !=0 and ny1 != 0)")
    problem.add_equation("cst = 0", condition="(nx == 0)")
    problem.add_equation("cst = 0", condition="(nx != 0) and ((ny0 == 0) or (ny1 == 0))")
else:
    logger.info("Using symmetry relation for czt.")
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
    cs = solver.state['cs']
    ctt = solver.state['ctt']

    if param.new_IC:
        k = 2
        a = 2
        b = 0
        ctt['g'] = param.pert_amp * np.sin(np.pi*y0/param.Ly)*np.sin(np.pi*y1/param.Ly)/(4*k*param.Lx) * (-2*a*b*np.cos(k*(2*param.Lx - x)) + 2*(a*b+(a**2 + b**2)*k*param.Lx)*np.cos(k*x) + (a**2 - b**2)*(np.sin(k*(2*param.Lx - x)) + np.sin(k*x)))
    else:
        r2 = x**2 + (param.Ly/np.pi*np.sin((y1-y0)*np.pi/param.Ly))**2/2
        ctt['g'] = param.pert_amp * np.exp(-r2/2/param.pert_width**2) * np.sin(np.pi/param.Ly *y1) * np.sin(np.pi/param.Ly *y0)

    # Invert cu_ref for cs initial condition
    # Reference jet: this will have a fractional symmetric component lambda
    if param.cu_ampl != 0:
        cu_ref = domain.new_field()
        cu_ref.meta['x']['constant'] = True
        cu_ref.meta['y0']['parity'] = 1
        cu_ref.meta['y1']['parity'] = 1

        x, y0, y1 = domain.grids()
        # Build as 1D function of y0
        k_cu = param.cu_k
        k_cu_parity = (param.cu_k+1) # opposite parity wrt centerline

        cu_ref['g'] = param.cu_ampl * (param.cu_lambda * np.cos(k_cu_parity*y0*np.pi/param.Ly) + (1 - param.cu_lambda)*np.cos(k_cu*y0*np.pi/param.Ly))
        # Diagonalize
        cu_ref = Diag(cu_ref, 'y0', 'y1').evaluate()
        
        ic_problem = de.LBVP(domain, variables=['cs'])
        ic_problem.meta['cs']['x']['constant'] = True
        ic_problem.meta['cs']['y0']['parity'] = 1
        ic_problem.meta['cs']['y1']['parity'] = -1

        ic_problem.parameters['cu_ref'] = cu_ref
        ic_problem.add_equation("cs = 0", condition="(nx != 0) or (ny0 != 0)")
        ic_problem.add_equation("cs = 0", condition="(nx == 0) and (ny1 == 0) and (ny0 == 0)")
        ic_problem.add_equation("-dy1(cs) = cu_ref", condition="(nx == 0) and (ny0 == 0) and (ny1 != 0)")
        ic_solver = ic_problem.build_solver()
        ic_solver.solve()
        cs['c'] = ic_solver.state['cs']['c']

# Analysis
an1 = solver.evaluator.add_file_handler('data_checkpoints', sim_dt=param.checkpoints_sim_dt, max_writes=1)
an1.add_system(solver.state)

an2 = solver.evaluator.add_file_handler('data_snapshots', iter=param.snapshots_iter, max_writes=10)
an2.add_task("interp(czz, y1=%.3f)" %(0.25*param.Ly), scales=2)
an2.add_task("interp(czz, y1=%.3f)" %(0.5*param.Ly), scales=2)
an2.add_task("interp(czz, y1=%.3f)" %(0.75*param.Ly), scales=2)
an2.add_task("interp(ctt, y1=%.3f)" %(0.25*param.Ly), scales=2)
an2.add_task("interp(ctt, y1=%.3f)" %(0.5*param.Ly), scales=2)
an2.add_task("interp(ctt, y1=%.3f)" %(0.75*param.Ly), scales=2)
an2.add_task("Diag(czz, 'y0', 'y1')", layout='c', name='zeta_power')
an2.add_task("Diag(ctt, 'y0', 'y1')", layout='c', name='theta_power')
an2.add_task("cz", layout='c', name='cz_power')
an2.add_task("ct", layout='c', name='ct_power')

an3 = solver.evaluator.add_file_handler('data_profiles', iter=param.profiles_iter, max_writes=1000)
an3.add_task("P1(cz)", name='cz')
an3.add_task("P1(cs)", name='cs')
an3.add_task("-dy1(P1(cs))", name='cu')
an3.add_task("P1(ct)", name='ct')

an4 = solver.evaluator.add_file_handler('data_scalars', iter=param.scalars_iter, max_writes=10000)
an4.add_task("-(Lx/2) * integ(P0(cz)*P0(cs) + P0(D(czs)))", name='KE')
an4.add_task("-(Lx/2) * integ(P0(cz)*P0(cs))", name='KE_mean')
an4.add_task(" (Lx/2) * integ(P0(cz)*P0(cz) + P0(D(czz)))", name='EN')
an4.add_task("integ((T(css) - css)**2)", name="css_asymm_L2")
an4.add_task("integ((T(ctt) - ctt)**2)", name="ctt_asymm_L2")
an4.add_task("integ((T(cts) - cst)**2)", name="cst_asymm_L2")
an4.add_task("integ(cs)", name="cs_integ")
an4.add_task("integ(ct)", name="ct_integ")

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=1)
flow.add_property("T(css) - css", name='css_sym')
flow.add_property("T(ctt) - ctt", name='ctt_sym')
flow.add_property("T(cts) - cst", name='cst_sym')
flow.add_property("ct", name="ct")
flow.add_property("cs", name="cs")
flow.add_property("css", name="css")
flow.add_property("cst", name="cst")
flow.add_property("cts", name="cts")
flow.add_property("ctt", name="ctt")

# construct projector for eigenvalue projection
if param.project_eigenvalues:
    eval_proj = projection.EigenvalueProjection(domain)

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    projection_time = 0.
    dt = param.dt

    while solver.ok:
        dt = solver.step(param.dt)
        if param.force_symmetry:
            if (solver.iteration-1) % param.force_symmetry == 0:
                logger.info("Iteration: %i, Enforcing symmetry." % (solver.iteration))
                enforce_symmetry(solver.state['css'])
                enforce_symmetry(solver.state['ctt'])
                enforce_symmetry(solver.state['cts'], solver.state['cst'])
        if param.fix_diagonal:
            if (solver.iteration-1) % param.fix_diagonal == 0:
                logger.info("Iteration: %i, Fixing 2nd cumulant diagonals." % (solver.iteration))
                diagonal.fix_diagonal(solver.state['css'])
                diagonal.fix_diagonal(solver.state['ctt'])
        if param.project_eigenvalues:
            if (solver.iteration-1) % param.project_eigenvalues == 0:
                logger.info("Iteration: %i, Projecting off negative eigenvalues." % (solver.iteration))
                start_projection_time = time.time()
                eval_proj.project_all(solver.state)
                end_projection_time = time.time()
                projection_time += end_projection_time - start_projection_time
        # eliminate means in 1st cumulants
        ct_mean = flow.grid_average('ct')
        cs_mean = flow.grid_average('cs')
        css_mean = flow.grid_average('css')
        cst_mean = flow.grid_average('cst')
        cts_mean = flow.grid_average('cts')
        ctt_mean = flow.grid_average('ctt')

        if (solver.iteration-1) % 1 == 0:
            #logger.info("cs_mean = {:e}".format(cs_mean))
            #logger.info("ct_mean = {:e}".format(ct_mean))
            solver.state['ct']['g'] -= ct_mean
            solver.state['cs']['g'] -= cs_mean
            solver.state['css']['g'] -= css_mean
            solver.state['cst']['g'] -= cst_mean
            solver.state['cts']['g'] -= cts_mean
            solver.state['ctt']['g'] -= ctt_mean

        # Hermitian projection
        if (solver.iteration-1) % 100 == 0:
            logger.info("Enforcing Hermitian symmetry")
            for field in solver.state.fields:
                field.require_grid_space()
        if (solver.iteration-1) % 10 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('(min, max) css_sym: (%e, %e), (min, max) ctt_sym: (%e, %e), cst_sym: (%e, %e)' %(flow.min('css_sym'), flow.max('css_sym'), flow.min('ctt_sym'), flow.max('ctt_sym'), flow.min('cst_sym'), flow.max('cst_sym')))
finally:
    try:
        solver.evaluate_handlers_now(dt)
    except:
        logger.error("Failed to write final output.") 
    end_time = time.time()
    run_time = end_time - start_time
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    if param.project_eigenvalues:
        logger.info('Projection time: %f'%projection_time)
        logger.info('Projection time/run_time: %f'%(projection_time/run_time))
    logger.info('Run time: %.2f sec' %(run_time))
    logger.info('Run cost: %f cpu-hr' %(domain.dist.comm_cart.size*run_time/60/60))

