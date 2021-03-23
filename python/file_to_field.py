"""tools to construct Dedalus fields from Dedalus h5py files.

"""
import h5py
import dedalus.public as de
import numpy as np

def limits_from_grid(basis_type, grid):
    if basis_type == 'Fourier':
        # even grid spacing
        delta = grid[1] - grid[0]
        left_edge = grid[0]
        right_edge = grid[-1] + delta
    elif basis_type == 'SinCos':
        delta = grid[1] - grid[0]
        left_edge = grid[0] - delta/2.
        right_edge = grid[-1] + delta/2.
    elif basis_type == 'Chebyshev':
        # Use Gauss points
        grid_size = len(grid)
        i = np.arange(grid_size)
        theta = pi * (i + 1/2) / grid_size
        native_grid = -np.cos(theta)
        problem_coord = grid

        xr = native_grid[0]/native_grid[1]
        c = (problem_coord[0] - problem_coord[0]*xr)/(1. - xr)
        r = (problem_coord[1] - c)/native_grid[1]
        left_edge = c - r
        right_edge = 2*c - a
    else:
        raise ValueError("Unknown basis type:", basis_type)
    
    return left_edge, right_edge

def domain_from_file(filename, basis_types, dealias=3/2):
    with h5py.File(filename,'r') as dfile:
        testkey = next(iter(dfile['tasks']))
        testdata = dfile['tasks'][testkey]
        dims = testdata.shape[1:] # trim write dimension
        dim_names = [i for i in testdata.attrs['DIMENSION_LABELS'][1:]]
        if not testdata.attrs['grid_space'][-1]:
            dims[-1] *= 2
        bases = []
        for n,b,d in zip(dim_names,basis_types,dims):
            if b == 'Fourier':
                limits = limits_from_grid('Fourier',dfile['/scales/'+n+'/1.0'])
                basis = de.Fourier(n, d, interval=limits, dealias=dealias)
            elif b == 'SinCos':
                limits = limits_from_grid('SinCos',dfile['/scales/'+n+'/1.0'])
                basis = de.SinCos(n, d, interval=limits, dealias=dealias)
            elif b == 'Chebyshev':
                limits = (-1, 1) # limits not implemented for Chebyshev
                basis = de.Chebyshev(n,d, limits=limits, dealias=dealias)
            else:
                raise ValueError("Unknown basis type:", basis_type)
            
            bases.append(basis)
        d = de.Domain(bases,grid_dtype=np.float64) # hardcoded domain dtype for now
        return d

def field_from_file(field, filename, basis_types, meta=None, index=-1, domain=None):
    if not domain:
        domain = domain_from_file(filename, basis_types)
    with h5py.File(filename,'r') as dfile:
        f = domain.new_field()
        dset = dfile['tasks'][field]
        if meta:
            for k,v in meta.items():
                f.meta[k] = v
        for layout in domain.dist.layouts:
            if np.allclose(layout.grid_space, dset.attrs['grid_space']):
                break
        else:
            raise ValueError("No matching layout")
        f[layout] = dset[(index,) + (Ellipsis,)]
        print("DNS sim_time = {:f}".format(dfile['scales/sim_time'][index]))
        f.name = field
    return f

def all_fields_from_file(filename, basis_types, meta=None, index=-1):
    domain = domain_from_file(filename, basis_types)
    field_names = []
    fields = []
    with h5py.File(filename,'r') as dfile:
        for field in dfile['tasks']:
            field_names.append(field)
    
    for fname in field_names:
        fields.append(field_from_file(fname, filename, basis_types,meta=meta, index=index, domain=domain))

    return fields
    
