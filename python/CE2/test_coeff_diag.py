import matplotlib.pyplot as plt
import numpy as np

import dedalus.public as de
from dedalus.tools.array import axslice, reshape_vector
from dedalus.core.evaluator import Evaluator
import diagonal

nx = 8
ny = 16

x = de.Fourier('x', nx)
y0 = de.SinCos('y0',ny)
y1 = de.SinCos('y1',ny)

dom = de.Domain([x,y0,y1],grid_dtype='float', mesh=[2,2])

data = dom.new_field()
data.meta['y1']['parity'] = -1
data.meta['y0']['parity'] = -1

xx, yy0, yy1 = dom.grids()
data['g'] = (1 + np.cos(xx) + np.cos(2*xx)) * (np.sin(yy0)*np.sin(2*yy1) + 10*np.sin(yy0)*np.sin(yy1))

op = diagonal.CoefficientDiagonal(data,'y0','y1')
out = op.evaluate()

data = Evaluator(dom, vars={'out':out})

output = data.add_file_handler("test_data",max_writes=1)
output.add_task("out", layout='c')
data.evaluate_handlers(data.handlers,timestep=0,sim_time=0,world_time=0, wall_time=0,iteration=0)


