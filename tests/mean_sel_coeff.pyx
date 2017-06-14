from cython.view cimport array as cvarray

cpdef double mean_s(double[:] arr):
    cdef double rv = 0.0
    cdef size_t i = 0
    cdef len = arr.shape[0]
    while i < len:
        rv += arr[i]
        i+=1
    return rv
    

