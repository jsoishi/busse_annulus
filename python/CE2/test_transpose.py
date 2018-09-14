import numpy as np

import dedalus.public as de
from transpose import TransposeOperator
import logging
logger = logging.getLogger(__name__)

de.operators.parseables['Trans'] = TransposeOperator

x = de.Fourier('x',16, dealias=3/2)
y0 = de.Fourier('y0',32, dealias=3/2)
y1 = de.Fourier('y1',32, dealias=3/2)

domain = de.Domain([x,y0,y1], grid_dtype='float', mesh=[2,2])

# problem = de.IVP(domain, variables=['f','g'])

# problem.add_equation("dt(f) = 0")
# problem.add_equation("g = Trans(f)")

# solver = problem.build_solver(de.timesteppers.RK111)

# yy0 = domain.grid(1)
# f = solver.state['f']
# f['g'] = yy0

# an = solver.evaluator.add_file_handler('snapshots', iter=1)
# an.add_system(solver.state)

# solver.stop_iteration = 2
# dt = 1
# while solver.ok:
#     solver.step(dt)


f = domain.new_field()
g = domain.new_field()
yy0 = domain.grid(1)
yy1 = domain.grid(2)

f['g'] = yy0
g['g'] = yy1
test = TransposeOperator(f).evaluate()
g.set_scales(domain.dealias)
assert np.allclose(test['g'], g['g'])
