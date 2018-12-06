"""Reverse first axis"""

import numpy as np

from dedalus.core.field import Operand
from dedalus.core.operators import Operator, FutureField, Interpolate
from dedalus.tools.array import reshape_vector, axslice

class ReverseFirst(Operator, FutureField):
    """Reverse the first direction:

    f(x, y, z) -> f(-x, y, z)

    """

    def __init__(self, arg, **kw):
        arg = Operand.cast(arg)
        super().__init__(arg, **kw)
        self.name = 'Reverse[%s]' % (arg.domain.bases[0].name)

    def meta_constant(self, axis):
        # Preserve constancy
        return self.args[0].meta[axis]['constant']

    def meta_parity(self, axis):
        return self.args[0].meta[axis]['parity']

    def check_conditions(self):
        layout = self.args[0].layout
        
        return (not layout.grid_space[0]
                and layout.grid_space[1]
                and layout.grid_space[2])

    def operate(self, out):
        arg = self.args[0]
        arg.require_grid_space()
        arg.require_coeff_space(axis=0)

        out.layout = arg.layout
        out.data = arg.data.conj()
