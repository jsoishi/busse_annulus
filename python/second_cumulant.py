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

def second_cumulant(y,f,g=None,layout='xy'):
    """Computes the second cumulant of two fields, f and g

    c_fg = <f'(x1, y1) g'(x2, y2)>

    where <> is an average in the x1 direction. Because of the
    periodicity in x, x2 = x1 + dx. The y parameter is y2, and is held
    constant. The cumulant is implemented by computing the
    x-correlation of f and g at each y1.

    Parameters
    ----------
    y : float
        y2 in the above equation; this is held constant
    f : dedalus field object
        the first field
    g : dedalus field object, optional
        the second field. if not given, the second cumulant of the f with itself is calculated
    layout : 'xy' or 'yx'
        gives the order of the x and y bases in teh domain

    """
    if g is None:
        g = f
    if layout == 'yx':
        yx = True
    elif layout == 'xy':
        yx = False
    else:
        raise ValueError("Layout must be one of 'yx' or 'xy'.")
    fields = [f,g]
    for f in fields:
        xb = f.domain.get_basis_object('x')
        Lx = xb.interval[1] - xb.interval[0]
        mean = (f.integrate('x')/Lx).evaluate()
        mean.set_scales(1.0)
        f.set_scales(1.0)
        f['g'] = f['g'] - mean['g']
    yarr = f.domain.get_basis_object('y').grid()
    idx = (np.abs(yarr - y)).argmin()

    outdata = np.empty_like(f['g'])
    if yx:
        for i in range(f['g'].shape[0]):
            outdata[i,:] = np.correlate(f['g'][i,:],g['g'][idx,:],mode='same')
    else:
        for i in range(f['g'].shape[1]):
            outdata[:,i] = np.correlate(f['g'][:,i],g['g'][:,idx],mode='same')

    return outdata

def second_cumulant_tavg(datafile, field1, field2, y, start, stop, layout='xy'):
    """compute time averaged second cumulant

    inputs
    ------
    datafile : hdf5 file containing dedalus outputs
    field1 : first field for cumulant
    field2 : second field for cumulant (if the same as field1, will do 2nd cumulant of field with itself)
    y : y point for comparison
    start : starting index in datafile
    start : ending index in datafile
    
    """
    cumulants = []
    for snap in range(start, stop+1):
        f1 = fields_from_file(datafile,['Fourier','SinCos'],field1,meta={'y':{'scale':1.0,'parity':-1,'constant':False}},index=snap)
        if field2 != field1:
            f2 = fields_from_file(datafile,['Fourier','SinCos'],field2,meta={'y':{'scale':1.0,'parity':-1,'constant':False}},index=snap)
        else:
            f2 = None
        c = second_cumulant(y,f1,g=f2,layout='xy')
        cumulants.append(c)
    cumulants = np.array(cumulants)
    
    return cumulants.mean(axis=0)


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
