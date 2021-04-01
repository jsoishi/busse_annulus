"""dns_to_ce2_data.py
Convert DNS output data to CE2 input.

Usage:
    dns_to_ce2_data.py <dns_file_path> <ce2_file_path> [--xres=<xres> --yres=<yres> --snap=<snap>]

Options:
    --xres=<xres>     Output x resolution [default: None] 
    --yres=<yres>     Output y resolution [default: None] 
    --snap=<snap>     snapshot index from DNS [default: -1]
"""
import time
import matplotlib.pyplot as plt
import h5py
import numpy as np
import pathlib
from docopt import docopt

import dedalus.public as de
from dedalus.core.evaluator import Evaluator
from dedalus.core.system import FieldSystem
from dedalus.tools.post import merge_process_files

from second_cumulant import all_second_cumulants, all_second_cumulants_spectral
from file_to_field import field_from_file
import logging
logger = logging.getLogger(__name__)

args = docopt(__doc__)

dns_file_path = pathlib.Path(args['<dns_file_path>'])
ce2_file_path = pathlib.Path(args['<ce2_file_path>'])

snap = int(args['--snap'])

theta = field_from_file('theta', dns_file_path, ['Fourier', 'SinCos'], meta={'y':{'parity':-1,'constant':False}},index=snap)
psi = field_from_file('psi', dns_file_path, ['Fourier', 'SinCos'], meta={'y':{'parity':-1,'constant':False}},index=snap)

# must use DNS with --xy option
x_interval = theta.domain.bases[0].interval
y_interval = theta.domain.bases[1].interval
Lx = x_interval[1] - x_interval[0]
Ly = y_interval[1] - y_interval[0]

# construct CE2 domain
# start with the same resolution as the DNS; during output, scale to the size required.
Nx = theta.domain.bases[0].grid_size(scale=1)
Ny = theta.domain.bases[1].grid_size(scale=1)

print(args['--xres'])
if args['--xres'] == 'None':
    output_Nx = Nx
else:
    output_Nx = int(args['--xres'])
if args['--yres'] == 'None':
    output_Ny = Ny
else:
    output_Ny = int(args['--yres'])


x_basis = de.Fourier('x', Nx, [-Lx/2, Lx/2])
y0_basis = de.SinCos('y0', Ny, [0, Ly])
y1_basis = de.SinCos('y1', Ny, [0, Ly])
domain = de.Domain([x_basis, y0_basis, y1_basis], grid_dtype=np.float64)
xi, y0, y1 = domain.grids()

out_fields = ['cs','ct', 'css', 'ctt', 'cts', 'cst']
y0_parity = [1, 1, -1, -1, -1, -1]

output_vars = {}

for name, parity in zip(out_fields, y0_parity):
    output_vars[name] =domain.new_field(name=name)
    output_vars[name].meta['y0']['parity'] = parity
    output_vars[name].meta['y1']['parity'] = -1

output_vars['ct'].meta['x']['constant'] = True
output_vars['cs'].meta['x']['constant'] = True

# calculate first cumulants
cs = psi.integrate('x')
cs.set_scales(1,keep_data=True)
ct = theta.integrate('x')
ct.set_scales(1,keep_data=True)
print("cs['g'].shape = {}".format(cs['g'].shape))
print("output_vars['cs']['g'].shape = {}".format(output_vars['cs']['g'].shape))
output_vars['cs']['g'] = cs['g'][:,np.newaxis,:]/Lx + 0.*y1
output_vars['ct']['g'] = ct['g'][:,np.newaxis,:]/Lx + 0.*y1
plt.plot(ct['g'][0,:])
plt.plot(output_vars['ct']['g'][0,0,:],'kx')
plt.savefig("output.png",dpi=400)

# calculate second cumulants
ingredients = ((psi, psi),(theta, theta), (theta, psi), (psi, theta))
inputs = dict(zip(out_fields[2:],ingredients))

start_time=time.time()
for k,v in inputs.items():
    v[0].set_scales(1, keep_data=True)
    v[1].set_scales(1, keep_data=True)
    logger.info("calculating {}".format(k))
    all_second_cumulants_spectral(output_vars[k], v[0], g=v[1],layout='xy')

end_time = time.time()
# create FileHandler to output data
field_system = FieldSystem(output_vars.values())
x_scale_factor = output_Nx/Nx
y_scale_factor = output_Ny/Ny
output_evaluator = Evaluator(field_system.domain, out_fields)
output_handler = output_evaluator.add_file_handler(ce2_file_path)
output_handler.add_system(field_system, scales=(x_scale_factor, y_scale_factor, y_scale_factor))

output_evaluator.evaluate_handlers(output_evaluator.handlers, timestep=0,sim_time=0, world_time=0, wall_time=0, iteration=0)

merge_process_files(ce2_file_path, cleanup=True)
logger.info("DNS to CE2 complete in {:f} sec".format(end_time-start_time))
