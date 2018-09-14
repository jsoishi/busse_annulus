import numpy as np

import dedalus.public as de
from reverse import ReverseFirst
import logging
logger = logging.getLogger(__name__)

de.operators.parseables['RFirst'] = ReverseFirst

x = de.Fourier('x',16, dealias=3/2, interval=[-np.pi, np.pi])
y0 = de.SinCos('y0',32, dealias=3/2)
y1 = de.SinCos('y1',32, dealias=3/2)

domain = de.Domain([x,y0,y1], grid_dtype='float', mesh=[2,2])
f = domain.new_field()
g = domain.new_field()
for field in [f, g]:
    for ax in ['y0', 'y1']:
        field.meta[ax]['parity'] = -1
    

xx = domain.grid(0)
yy0 = domain.grid(1)
yy1 = domain.grid(2)

#test_func = lambda xx, yy, zz: np.cos(xx) + np.sin(4*xx)
test_func = lambda xx, yy, zz: np.cos(xx) + np.sin(4*xx)*np.sin(yy0)

f['g'] = test_func(xx, yy0, yy1)
g['g'] = test_func(-xx, yy0, yy1)

g.set_scales(domain.dealias)
test = ReverseFirst(f).evaluate()

assert np.allclose(test['g'], g['g'])

