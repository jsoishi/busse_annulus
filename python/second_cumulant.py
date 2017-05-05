"""second_cumulant.py

Calculate the second cumulant of DNS variables by post processing

"""
import copy
import sys
import h5py
import numpy as np
import dedalus.public as de

def limits_from_grid(basis_type,grid):
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
        raise NotImplementedError("Chebyshev limits not yet supported")
    else:
        raise ValueError("Unknown basis type:", basis_type)
    
    return left_edge, right_edge

def domain_from_file(filename,basis_types,dealias=3/2):
    with h5py.File(filename,'r') as dfile:
        testkey = next(iter(dfile['tasks']))
        testdata = dfile['tasks'][testkey]
        dims = testdata.shape[1:] # trim write dimension
        dim_names = [i.decode('utf-8') for i in testdata.attrs['DIMENSION_LABELS'][1:]]
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

def fields_from_file(filename,basis_types,field,meta=None,index=-1):
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
        f[layout] = dset[(index,) + (slice(None),slice(None))]
    return f

def second_cumulant(y, field1,field2=None):
    """function to calcluate the second cumulant of data. y is the y
    component of the r1 vector.
    """
    if field2 is None:
        field2 = field1
    
    yarr = field1.domain.get_basis_object('y').grid()
    idx = (np.abs(yarr - y)).argmin()
    
    xb_old = field1.domain.get_basis_object('x')
    yb_old = field1.domain.get_basis_object('y')
    x_basis = type(xb_old)(xb_old.name,xb_old.base_grid_size,xb_old.interval,xb_old.dealias)
    y_basis = type(yb_old)(yb_old.name,yb_old.base_grid_size,yb_old.interval,yb_old.dealias)
    xprime_basis = type(xb_old)('xprime',xb_old.base_grid_size,xb_old.interval,xb_old.dealias)
    grid_dtype = field1['g'].dtype
    new_domain = de.Domain([y_basis,x_basis,xprime_basis],grid_dtype=grid_dtype)
    sc_raw = new_domain.new_field()
    if 'parity' in field1.meta['y'].keys():
        sc_raw.meta['y']['parity'] = field1.meta['y']['parity']*field2.meta['y']['parity']
    
    sc_raw['g'] = np.einsum('i,jk->jik',field1['g'][idx,:], field2['g'])
    sc = sc_raw.integrate('xprime')
    sc.set_scales('1')
    
    return sc

if __name__ == "__main__":
    basedir = sys.argv[-1]
    filterfrac = None
    basedir = basedir.rstrip('/')
    print("basedir")
    print(basedir.split('/'))
    run_name = basedir.split('/')[-1]
    basename = 'busse_annulus'
    print(run_name)
    params = parse_params(run_name,basename)
