import scipy.linalg as linalg
import numpy as np

import logging
logger = logging.getLogger(__name__)

class EigenvalueProjection:
    def __init__(self, domain, variables=['css','cst','cts','ctt']):
        self._variables = variables
        self.domain = domain
        self._build_subcommunicator()

        dtype = 'complex128'
        self.recbuf = None
        self.workbuf = None
        self.rank_recbuf = None
        self.rank_kx = None
        if self.rank == 0:
            self.local_shape = self.domain.local_coeff_shape
            self.global_shape = self.domain.global_coeff_shape
            self.rec_shape = [self.size,] + list(self.local_shape)
            self.work_shape = [self.local_shape[0],] + list(self.global_shape[1:3])
            # next line looks hard coded for 4 cumulants...
            self.total_work_shape = [self.local_shape[0],] + [2*i for i in self.global_shape[1:3]]
            self.recbuf = np.zeros(self.rec_shape, dtype=dtype)
            self.workbuf = np.zeros(self.work_shape, dtype=dtype)
            self.total_workbuf = np.zeros(self.total_work_shape, dtype=dtype)
            logger.debug("dtype = {}".format(self.recbuf.dtype))
            logger.debug("local_shape = {}".format(self.local_shape))
            logger.debug("rec_shape = {}".format(self.rec_shape))
            logger.debug("work_shape = {}".format(self.work_shape))
            logger.debug("total_work_shape = {}".format(self.total_work_shape))
            
        if self.domain.dist.comm.rank == 0:
            self.rank_recbuf = np.zeros(self.global_shape[0], dtype=int)
            self.rank_kx = np.zeros(self.global_shape[0], dtype=np.float64)
        
    def _build_subcommunicator(self):
        comm_cart = self.domain.dist.comm_cart
        comm_w = self.domain.dist.comm
        if comm_cart.size == 1:
            self.comm = comm_cart
            raise ValueError("Not implemented for meshes of dimension other than 2")
        else:
            self.comm = comm_cart.Sub([False,True])

            # build communicator for sending rank at each k_x mode
            group_ranks = []
            for i in range(comm_cart.dims[1]):
                group_ranks.append(comm_cart.Get_cart_rank([i,0]))
            new_group = comm_w.group.Incl(group_ranks)
            self.rank_comm = comm_w.Create_group(new_group)

        self.rank = self.comm.rank
        self.size = self.comm.size

    def project_all(self, state, thresh=1e-12):
        """
        """
        # pack total second cumulant work buffer
        if self.rank == 0:
            nvar = len(self._variables)/2
            Ny = self.work_shape[1]

        for i,field in enumerate(self._variables):
            logger.debug("gathering field {}".format(field))
            self.gather(state[field])
            logger.debug("field {} gathered".format(field))
            if self.rank == 0:
                logger.debug("combining field {}".format(field))
                rstart = int(i/nvar)*Ny
                rend = rstart + Ny
                cstart = int(i%nvar)*Ny
                cend = cstart + Ny
                logger.debug("rstart:rend, cstart:cend = {}:{}, {}:{}".format(rstart,rend,cstart,cend))
                shape = self.total_workbuf[:, rstart:rend, cstart:cend].shape
                logger.debug("total_workbuf shape = {}, {}, {}".format(shape[0],shape[1],shape[2]))
                self.total_workbuf[:, rstart:rend, cstart:cend] = self.workbuf
                logger.debug("field {} combined".format(field))

        # project total second cumulant buffer
        if self.rank == 0:
            logger.debug("projecting all.")
            nx = self.local_shape[0]
            rank2c = np.zeros(nx, dtype=int)
            for m in range(nx):
                A = self.total_workbuf[m,:,:]
                Aproj, mrank = self.projection(A, thresh=thresh)
                rank2c[m] = mrank
                self.total_workbuf[m,:,:] = Aproj
            self.rank_comm.Gather(rank2c,self.rank_recbuf, root=0)
            
            self.rank_comm.Gather(self.domain.elements(0).ravel(),self.rank_kx, root=0)
            if self.rank_comm.rank == 0:
                logger.info("second cumulant kx:   {}".format(self.rank_kx))
                logger.info("second cumulant rank: {}".format(self.rank_recbuf))
            logger.debug("projection complete.")
        # unpack total second cumulant work buffer and scatter
        for i,field in enumerate(self._variables):
            if self.rank == 0:
                rstart = int(i/nvar)*Ny
                rend = rstart + Ny
                cstart = int(i%nvar)*Ny
                cend = cstart + Ny
                self.workbuf = self.total_workbuf[:, rstart:rend, cstart:cend]
            logger.debug("scattering field {}".format(field))            
            self.scatter(state[field])
            logger.debug("field {} scattered".format(field))

    def projection(self, data, thresh=1e-12):
        # do the projection
        # A = Q Λ Qinv
        # where Λ is the matrix of eigenvalues and Q is the matrix of eigenvectors,
        # and we delete all negative eigenvalues from Λ

        H = 0.5*(data + data.conj().T) # hermitian part
        AH = data - H # antihermitian part
        AH_norm = np.linalg.norm(AH)
        H_norm = np.linalg.norm(H)
        if H_norm != 0:
            if H_norm < 1e-12:
                H_norm += 1e-5
            logger.info("||AH||/||H|| : {}".format(AH_norm/H_norm))
        # A is a positive definite matrix iff H, the hermitian part, is positive definite
        evals, Q = linalg.eigh(H)
        index = (evals < -thresh) # find negative eigenvalues
        if np.any(index):
            logger.warning("{} negative eigenvalues found".format(np.sum(index)))

        mrank = (np.abs(evals) > thresh).sum()
        evals[index] = 0
        Qinv = Q.conj().T # if input is Normal, Q is unitary, and all H matricies are Normal!
        Lambda = np.diag(evals)

        return Q @ (Lambda @ Qinv), mrank

    def project(self, data, thresh=1e-12):
        self.gather(data)
        if self.rank == 0:
            nx = self.local_shape[0]
            for m in range(nx):
                A = self.workbuf[m,:,:]
                self.workbuf[m,:,:] = self.projection(A, thresh=thresh)
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
        recbuf = np.zeros(rec_shape, dtype=data['c'].dtype)

    comm.Gather(sendbuf, recbuf, root=0)
    
    # roots of subcommunicators do projection
    if rank == 0:
        # unpack gathered buffer 
        workbuf = np.zeros(work_shape, dtype=data['c'].dtype)
        
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
