#See http://cython.readthedocs.io/en/latest/src/tutorial/numpy.html
import numpy as np
cimport numpy as np
DTYPE = np.float
ctypedef np.float_t DTYPE_t
cdef class MeanFitness(object):
    cdef object d
    def __cinit__(self):
        self.d = []
    def __call__(self,pop):
        cdef size_t i = 0
        cdef size_t N = pop.N
        cdef np.ndarray[DTYPE_t,ndim=1] w = np.zeros(N,dtype=DTYPE)
        while i < N:
            w[i]=pop.diploids[i].w
            i += 1
        self.d.append((pop.generation,w.mean()))
    def get_data(self):
        return self.d
