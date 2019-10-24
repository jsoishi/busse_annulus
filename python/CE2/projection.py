import scipy.linalg as linalg
import numpy as np

import logging
logger = logging.getLogger(__name__)

class EigenvalueProjection:
    def __init__(self, domain):
        self.domain = domain
        self._build_subcommunicator()

        dtype = 'complex128'
        self.recbuf = None
        self.workbuf = None
        if self.rank == 0:
            self.local_shape = self.domain.local_coeff_shape
            self.global_shape = self.domain.global_coeff_shape
            self.rec_shape = [self.size,] + list(self.local_shape)
            self.work_shape = [self.local_shape[0],] + list(self.global_shape[1:3])
            self.recbuf = np.empty(self.rec_shape, dtype=dtype)
            self.workbuf = np.empty(self.work_shape, dtype=dtype)
            logger.debug("dtype = {}".format(self.recbuf.dtype))
            logger.debug("local_shape = {}".format(self.local_shape))
            logger.debug("rec_shape = {}".format(self.rec_shape))
            logger.debug("work_shape = {}".format(self.work_shape))

    def _build_subcommunicator(self):
        if self.domain.dist.comm_cart.size == 1:
            self.comm = self.domain.dist.comm_cart
        else:
            self.comm = self.domain.dist.comm_cart.Sub([False,True])
        self.rank = self.comm.rank
        self.size = self.comm.size

    def project(self, data, thresh=1e-12):
        self.gather(data)
        # do the projection
        # A = Q Λ Qinv
        # where Λ is the matrix of eigenvalues and Q is the matrix of eigenvectors,
        # and we delete all negative eigenvalues from Λ
        if self.rank == 0:
            nx = self.local_shape[0]
            for m in range(nx):
                A = self.workbuf[m,:,:]
                H = 0.5*(A + A.conj().T) # hermitian part
                AH = A - H # antihermitian part
                # A is a positive definite matrix iff H, the hermitian part, is positive definite
                evals, Q = linalg.eigh(H)
                index = (evals < -thresh) # find negative eigenvalues
                if np.any(index):
                    logger.warning("{} negative eigenvalues found for {} at x-mode {}".format(np.sum(index),data.name,m))
                evals[index] = 0
                Qinv = Q.conj().T # if input is Normal, Q is unitary, and all H matricies are Normal!
                Lambda = np.diag(evals)
                self.workbuf[m,:,:] = Q @ (Lambda @ Qinv) + AH
        self.scatter(data)

    def gather(self, data):
        self.comm.Gather(data['c'], self.recbuf, root=0)
        if self.rank == 0:
            blocksize = self.rec_shape[2]
            for block in range(self.rec_shape[0]):
                self.workbuf[:,block*blocksize:(block+1)*blocksize,:] = self.recbuf[block]

    def scatter(self, data):
        if self.rank == 0:
            blocksize = self.rec_shape[2]
            for block in range(self.rec_shape[0]):
                self.recbuf[block] = self.workbuf[:,block*blocksize:(block+1)*blocksize,:]
        self.comm.Scatter(self.recbuf,data['c'], root=0)
    

def project_eigenvalues(data, comm, thresh=1e-16):
    rank = comm.rank
    size = comm.size

    # gather (y0, y1) planes on subcommunicator roots
    sendbuf = data['c']
    recbuf = None
    workbuf = None
    if rank == 0:
        local_shape = data.domain.local_coeff_shape
        global_shape = data.domain.global_coeff_shape
        rec_shape = [size,] + list(local_shape)
        work_shape = [local_shape[0],] + list(global_shape[1:3])
        recbuf = np.empty(rec_shape, dtype=data['c'].dtype)

    comm.Gather(sendbuf, recbuf, root=0)
    
    # roots of subcommunicators do projection
    if rank == 0:
        # unpack gathered buffer 
        workbuf = np.empty(work_shape, dtype=data['c'].dtype)
        
        # project eigenvalues
        nx = local_shape[0]
        for m in range(nx):
            sub_data = workbuf[m,:,:]
            # this matrix must be Hermitian
            #evals, evecs = linalg.eigh(sub_data)
            evals, evecs = linalg.eig(sub_data)
            index = (evals < 0) # find negative eigenvalues
            # if np.any(index):
            #     logger.warning("Negative eigenvalue found for {} at x-mode {}".format(data.name,m))
            if np.all(index):
                logger.warning("All eigenvalues negative for {} at x-mode {}".format(data.name,m))
            evals[index] = 0 #thresh
            Pinv = evecs.conj().T # if input is Normal, P is unitary!
            Lambda = np.diag(evals)
            workbuf[m,:,:] = evecs @ (Lambda @ Pinv)

        # pack buffer for scattering
        tmp_shape = list(rec_shape)
        tmp_shape[1], tmp_shape[2] = tmp_shape[2], tmp_shape[1]
        workbuf = workbuf.swapaxes(0,1).reshape(tmp_shape)
        workbuf = workbuf.swapaxes(1,2)
    # scatter data back to all processors
    comm.Scatter(workbuf,data['c'], root=0)

def unpack_buffer(data):
    pass

def pack_buffer(data):
    pass
