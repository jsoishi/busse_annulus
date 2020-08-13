"""second_cumulant.py

Calculate the second cumulant of DNS variables by post processing

"""
import copy
import sys
import h5py
import numpy as np
import dedalus.public as de
from file_to_field import field_from_file

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

    for yidx in range(ny):
        if yx:
            dslice = slice(yidx, None, None)
        else:
            dslice = slice(None, None, yidx)
        output[dslice] = second_cumulant(yidx, f, g=g, layout=layout)

    return output


def second_cumulant(y,f,g=None,layout='xy'):
    """Computes the second cumulant of two fields, f and g

    c_fg = <f'(x1, y1) g'(x2, y2)>

    where <> is an average in the x1 direction. Because of the
    periodicity in x, x2 = x1 + dx. The y parameter is y2, and is held
    constant. The cumulant is implemented by computing the
    x-correlation of f and g at each y1.

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
        raise ValueError("Layout must be one of 'yx' or 'xy'.")
    fields = [f,g]
    for f in fields:
        xb = f.domain.get_basis_object('x')
        Lx = xb.interval[1] - xb.interval[0]
        mean = (f.integrate('x')/Lx).evaluate()
        mean.set_scales(1.0)
        f.set_scales(1.0)
        f['g'] = f['g'] - mean['g']
    if type(y) is int:
        idx = y
    elif type(y) is float:
        yarr = f.domain.get_basis_object('y').grid()
        idx = (np.abs(yarr - y)).argmin()
    else:
        raise ValueError("y must be int or float, not {}".format(type(y)))
    

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
