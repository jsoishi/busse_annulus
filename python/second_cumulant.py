"""second_cumulant.py

Calculate the second cumulant of DNS variables by post processing

"""
import copy
import sys
import h5py
import numpy as np
import dedalus.public as de
from file_to_field import field_from_file
import logging
from numba import jit
logger = logging.getLogger(__name__)

def all_second_cumulants_spectral(output_field, f, g=None, layout='xy'):
    if g is None:
        g = f

    if layout != 'xy':
        raise NotImplementedError("Use xy ordering from DNS.")

    f_3d = output_field.domain.new_field()
    f_3d.meta['y0']['parity'] = f.meta['y']['parity']
    f_3d.meta['y1']['parity'] = 1
    g_3d = output_field.domain.new_field()
    g_3d.meta['y0']['parity'] = 1
    g_3d.meta['y1']['parity'] = g.meta['y']['parity']

    # works in serial, at least
    f_slice = (slice(None), slice(None), 0)
    g_slice = (slice(None), 0, slice(None))
    f_3d['c'][f_slice] = f['c']
    f_3d['c'][0,:,:] = 0
    g_3d['c'][g_slice] = g['c']
    g_3d['c'][0,:,:] = 0
    layout = output_field.domain.dist.layouts[-2]
    for field in [f_3d, g_3d]:
        field.require_layout(layout)

    kx = f_3d.domain.elements(0)
    Lx = f_3d.domain.bases[0].interval[0]
    output_field[layout] = f_3d.data * np.conj(g_3d.data)* np.exp(1j*Lx*kx)

def all_second_cumulants(f, g=None, layout='xy'):
    """Computes a full second cumulants
    
    c_fg(xi, y1, y2) = <f'(x1, y1) g'(x2,y2)>

    assuming periodicity in x, so xi = x2 - x1.

    """
    if layout == 'xy':
        nx,ny = f['g'].shape
        output = np.zeros((nx, ny, ny))
        yx = False
    elif layout == 'yx':
        ny,nx = f['g'].shape
        output = np.zeros((ny, ny, nx))
        yx = True
    else:
        raise ValueError("Layout must be one of 'yx' or 'xy'.")

    n_integ = nx*ny*ny
    n = 0
    for yidx in range(ny):
        logger.info("Integral {}/{}".format(n,n_integ))
        if yx:
            dslice = (yidx, slice(None, None, None), slice(None, None, None))
        else:
            dslice = (slice(None, None, None), slice(None, None, None), yidx)
        output[dslice] = second_cumulant(yidx, f, g=g, layout=layout)
        n += nx*ny
    return output

def second_cumulant(y,f,g=None,layout='xy'):
    """Computes the second cumulant of two fields, f and g

    c_fg = <f'(x1, y1) g'(x2, y2)>

    where <> is an average in the x1 direction. Because of the
    periodicity in x, x2 = x1 + dx. The y parameter is y2, and is held
    constant. 

    WARNING: This is an N**2 operation that uses Dedalus integrate() operators.
    It's slow but accurate.

    Parameters
    ----------
    y : float or int
        y2 in the above equation; this is held constant
        if a float, will calculate corresponding index for g'
        if an int, will use this as the index into the g' array
    f : dedalus field object
        the first field
    g : dedalus field object, optional
        the second field. if not given, the second cumulant of the f with itself is calculated
    layout : 'xy' or 'yx'
        gives the order of the x and y bases in the domain
    """
    if g is None:
        g = f
    if layout == 'yx':
        yx = True
    elif layout == 'xy':
        yx = False
    else:
        print("layout must be xy or yx.")
        raise
    fields = [f,g]
    for field in fields:
        xb = field.domain.get_basis_object('x')
        Lx = xb.interval[1] - xb.interval[0]
        mean = (field.integrate('x')/Lx).evaluate()
        mean.set_scales(1.0)
        field.set_scales(1.0)
        field['g'] = field['g'] - mean['g']

    if type(y) is int:
        idx = y
    elif type(y) is float:
        yarr = f.domain.get_basis_object('y').grid()
        idx = (np.abs(yarr - y)).argmin()
    else:
        print("y must be int or float, not {}".format(type(y)))
        raise 
    

    outdata = np.empty_like(f['g'])
    if yx:
        nx = f['g'].shape[1]
    else:

        nx = f['g'].shape[0]

    tmp_x = de.Fourier('x',nx, interval=xb.interval)
    tmp_d = de.Domain([tmp_x,],grid_dtype=np.float64)
    integrand = tmp_d.new_field()

    if yx:
        ny = f['g'].shape[0]
    else:
        ny = f['g'].shape[1]
    n_int = nx*ny
    n = 0
    for iy in range(ny):
        for xi in range(nx):
            x_roll = xi - nx//2

            if yx:
                f_slice = (iy, slice(None))
                g_slice = (idx, slice(None))
                out_slice = (iy, xi)
            else: 
                f_slice = (slice(None), iy)
                g_slice = (slice(None), idx)
                out_slice = (xi, iy)
            integrand['g'] = f['g'][f_slice] * np.roll(g['g'][g_slice],x_roll)
            outdata[out_slice] = integrand.integrate()['g'][0]/Lx
            n += 1
                                         
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
        f1 = field_from_file(datafile,['Fourier','SinCos'],field1,meta={'y':{'scale':1.0,'parity':-1,'constant':False}},index=snap)
        if field2 != field1:
            f2 = field_from_file(datafile,['Fourier','SinCos'],field2,meta={'y':{'scale':1.0,'parity':-1,'constant':False}},index=snap)
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
