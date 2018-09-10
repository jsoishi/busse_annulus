"""Fourier diagonal operator definition."""

import numpy as np

from dedalus.core.field import Operand
from dedalus.core.operators import Operator, FutureField, Interpolate
from dedalus.tools.array import reshape_vector, axslice


class FourierDiagonal(Operator, FutureField):
    """
    Double-Fourier interpolation-on-diagonal operator.

    Parameters
    ----------
    arg : field object
        Field argument
    basis0, basis1 : basis identifiers
        Bases for diagonal interpolation

    Notes
    -----
    The return data is structured such that, for g = Diag(f), and x1,x2 in [a,b],
        g(x1, x2) = f(x1+x2-a, x1+x2-a)
    s.t. interpolation of g at x1=a or x2=a yields the diagonalization of f,
    i.e. f(x, x), arranged as a one-dimensional function of x2 or x1, respectively.

    """

    def __init__(self, arg, basis0, basis1, **kw):
        arg = Operand.cast(arg)
        super().__init__(arg, **kw)
        self.basis0 = self.domain.get_basis_object(basis0)
        self.basis1 = self.domain.get_basis_object(basis1)
        self.axis0 = self.domain.bases.index(self.basis0)
        self.axis1 = self.domain.bases.index(self.basis1)
        if self.axis0 > self.axis1:
            raise ValueError("Cannot evaluate specified axis order.")
        if self.basis0.interval != self.basis1.interval:
            raise ValueError("Bases must occupy same interval.")
        self.name = 'Diag[%s=%s]' %(self.basis0.name, self.basis1.name)
        # Shear array
        k0 = reshape_vector(self.basis0.wavenumbers, dim=self.domain.dim, axis=self.axis0)
        x1 = self.domain.grid(self.axis1, scales=self.domain.dealias)
        dx0 = x1 - self.basis1.interval[0]
        self.shear = np.exp(1j*k0*dx0)
        # Filter mask
        slices = self.domain.dist.coeff_layout.slices(self.domain.dealias)
        k0 = reshape_vector(self.basis0.wavenumbers[slices[self.axis0]], dim=self.domain.dim, axis=self.axis0)
        k1 = reshape_vector(self.basis1.wavenumbers[slices[self.axis1]], dim=self.domain.dim, axis=self.axis1)
        self.filter_mask = (k0 == k1)

    def meta_constant(self, axis):
        # Preserve constancy
        return self.args[0].meta[axis]['constant']

    def check_conditions(self):
        # Shearing layout
        layout = self.args[0].layout
        return ((layout.grid_space[self.axis1]) and
                (not layout.grid_space[self.axis0]) and
                (layout.local[self.axis0]))

    def operate(self, out):
        arg = self.args[0]
        axis0 = self.axis0
        axis1 = self.axis1
        # Enforce conditions for shearing space
        arg.require_grid_space(axis=axis1)
        arg.require_coeff_space(axis=axis0)
        arg.require_local(axis=axis0)
        # Apply Fourier shear to flatten the diagonal
        # s.t out(y0, y1) = arg(y0+y1-a, y1)
        out.layout = arg.layout
        np.multiply(arg.data, self.shear, out=out.data)
        # Interpolate on flattened diagonal
        # s.t. out(y0, y1) = arg(y1, y1)
        self.basis0.Interpolate(out, 'left', out=out).evaluate()
        out.meta['y0']['constant'] = False
        # Broadcast and filter coefficients
        # s.t. out(y0, y1) = arg(y0+y1-a, y0+y1-a)
        out.data[axslice(axis0, None, None)] = out.data[axslice(axis0, 0, 1)]
        out.require_coeff_space(axis=axis1)
        out.data *= self.filter_mask
        # Move back to starting layout
        out.require_layout(arg.layout)

class SinCosDiagonal(Operator, FutureField):
    """SinCos  interpolation-on-diagonal operator.

    Parameters
    ----------
    arg : field object
        Field argument
    basis0, basis1 : basis identifiers
        Bases for diagonal interpolation

    Notes
    -----
    The return data is structured such that, for g = Diag(f), and x1,x2 in [a,b],
        g(x1, x2) = f(x1+x2-a, x1+x2-a)
    s.t. interpolation of g at x1=a or x2=a yields the diagonalization of f,
    i.e. f(x, x), arranged as a one-dimensional function of x2 or x1, respectively.

    """

    def __init__(self, arg, basis0, basis1, **kw):
        arg = Operand.cast(arg)
        super().__init__(arg, **kw)
        self.basis0 = self.domain.get_basis_object(basis0)
        self.basis1 = self.domain.get_basis_object(basis1)
        self.axis0 = self.domain.bases.index(self.basis0)
        self.axis1 = self.domain.bases.index(self.basis1)
        if self.axis0 > self.axis1:
            raise ValueError("Cannot evaluate specified axis order.")
        if self.basis0.interval != self.basis1.interval:
            raise ValueError("Bases must occupy same interval.")
        self.name = 'Diag[%s=%s]' %(self.basis0.name, self.basis1.name)
        # Filter mask
        self.slices = self.domain.dist.coeff_layout.slices(self.domain.dealias)
        k0 = reshape_vector(self.basis0.wavenumbers[self.slices[self.axis0]], dim=self.domain.dim, axis=self.axis0)
        k1 = reshape_vector(self.basis1.wavenumbers[self.slices[self.axis1]], dim=self.domain.dim, axis=self.axis1)
        self.filter_mask = (k0 == k1)

    def meta_constant(self, axis):
        # Preserve constancy
        return self.args[0].meta[axis]['constant']

    def meta_parity(self, axis):
        return self.args[0].meta[axis]['parity']

    def check_conditions(self):
        layout = self.args[0].layout
        return ((layout.grid_space[self.axis1]) and
                (layout.grid_space[self.axis0]) and
                (layout.local[self.axis1]))

    def operate(self, out):
        arg = self.args[0]
        axis0 = self.axis0
        axis1 = self.axis1

        arg.require_grid_space(axis=0)
        arg.require_coeff_space(axis=axis0)
        arg.require_coeff_space(axis=axis1)
        arg.require_local(axis=axis1)

        out.layout = arg.layout
        k0 = self.basis0.elements[self.slices[self.axis0]]
        k1 = self.basis0.elements[self.slices[self.axis1]]
        y = self.basis0.grid(scale=self.domain.dealias[self.axis0])
        phi = np.einsum('ik,jk->ijk',np.sin(k0[:,None]*y),np.sin(k1[:,None]*y))
        di = np.expand_dims(np.einsum('ijk,jkl',arg.data,phi), axis=1)

        out.meta[axis0]['parity'] = 1
        out.meta[axis1]['parity'] = 1
        out.require_grid_space(axis=axis1)
        out.require_coeff_space(axis=axis0)
        out.require_local(axis=axis0)

        out.data[axslice(axis0,0,1)] = di
        out.data[axslice(axis0,1,None)] = 0
        out.data[axslice(axis0,None, None)] = out.data[axslice(axis0,0,1)]
        out.require_coeff_space(axis=axis1)
        out.data *= self.filter_mask

        out.require_layout(arg.layout)


class GridDiagonal(Operator, FutureField):
    """
    Grid-based interpolation-on-diagonal operator.

    Parameters
    ----------
    arg : field object
        Field argument
    basis0, basis1 : basis identifiers
        Bases for diagonal interpolation

    Notes
    -----
    The return data is structured such that, for g = Diag(f), and x1,x2 in [a,b],
        g(x1, x2) = f(x2, x2)

    """

    def __init__(self, arg, basis0, basis1, **kw):
        arg = Operand.cast(arg)
        super().__init__(arg, **kw)
        self.basis0 = self.domain.get_basis_object(basis0)
        self.basis1 = self.domain.get_basis_object(basis1)
        self.axis0 = self.domain.bases.index(self.basis0)
        self.axis1 = self.domain.bases.index(self.basis1)
        if self.axis0 > self.axis1:
            raise ValueError("Cannot evaluate specified axis order.")
        if self.basis0.coeff_size != self.basis1.coeff_size:
            raise ValueError("Bases must be same size.")
        if self.basis0.interval != self.basis1.interval:
            raise ValueError("Bases must occupy same interval.")
        if self.basis0.dealias != self.basis1.dealias:
            raise ValueError("Bases must have same dealiasing.")
        self.name = 'Diag[%s=%s]' %(self.basis0.name, self.basis1.name)

    def meta_constant(self, axis):
        if axis == self.axis0:
            return True
        elif axis == self.axis1:
            c0 = self.args[0].meta[self.axis0]['constant']
            c1 = self.args[0].meta[self.axis1]['constant']
            return (c0 and c1)
        else:
            return self.args[0].meta[axis]['constant']

    def meta_parity(self, axis):
        if axis == self.axis0:
            return 1
        elif axis == self.axis1:
            p0 = self.args[0].meta[self.axis0]['parity']
            p1 = self.args[0].meta[self.axis1]['parity']
            return p0*p1
        else:
            return self.args[0].meta[axis]['parity']

    def check_conditions(self):
        # Shearing layout
        layout = self.args[0].layout
        return ((layout.grid_space[self.axis1]) and
                (layout.grid_space[self.axis0]) and
                (layout.local[self.axis0]))

    def operate(self, out):
        arg = self.args[0]
        axis0 = self.axis0
        axis1 = self.axis1
        dim = self.domain.dim
        # Enforce conditions for shearing space
        arg.require_grid_space(axis=axis1)
        arg.require_grid_space(axis=axis0)
        arg.require_local(axis=axis0)
        # Interpolate along x0
        x1_slices = arg.layout.slices(self.domain.dealias)[axis1]
        x1_start = x1_slices.start
        x1_stop = x1_slices.stop
        out.layout = arg.layout
        # Get same indeces along x0 as x1
        arg_slices = [slice(None) for i in range(dim)]
        arg_slices[axis0] = np.arange(x1_start, x1_stop)
        arg_slices[axis1] = np.arange(0, x1_stop-x1_start)
        # Expand to broadcast over all x0
        exp_slices = [slice(None) for i in range(dim)]
        exp_slices[axis0] = None
        out.data[:] = arg.data[tuple(arg_slices)][tuple(exp_slices)]

