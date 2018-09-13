"""Transpose operator."""

from mpi4py import MPI
import numpy as np

from dedalus.core.field import Operand
from dedalus.core.operators import Operator, FutureField, Interpolate
from dedalus.tools.array import reshape_vector, axslice

import logging
logger = logging.getLogger(__name__)


class TransposeOperator(Operator, FutureField):
    """Real space transpose operator. Transposes last two dimensions of a
    multidimensional array. Designed to produce second cumulant with
    arguments swapped, e.g.

    c_zt(x, y0, y1) = c_tz(x, y1, y0)

    This is done in real space, with x local. It requires
    communication, so should be carefully monitored for speed.

    """

    def __init__(self, args, **kw):
        super().__init__(args,**kw)
        self.name = 'TransposeLast'

    def meta_constant(self, axis):
        # Preserve constancy
        return self.args[0].meta[axis]['constant']

    def check_conditions(self):
        # Field must be in grid layout
        return (self.args[0].layout is self._grid_layout)

    def operate(self, out):
        arg = self.args[0]
        arg.require_grid_space()
        out.layout = self._grid_layout
        layout = arg.layout
        comm = arg.domain.dist.comm_cart
        
        # exchange data with transposed cartesian communicator rank
        trans_pair = comm.Get_cart_rank(layout.ext_coords[1:][::-1])
        logger.debug("communicating with {:d}".format(trans_pair))

        # transpose local data
        output_buffer = np.transpose(arg.data, axes=[0,2,1]).copy()

        req_send = comm.Isend(output_buffer, dest=trans_pair)
        req_recv = comm.Irecv(out.data, source=trans_pair)

        MPI.Request.Waitall([req_send, req_recv])


