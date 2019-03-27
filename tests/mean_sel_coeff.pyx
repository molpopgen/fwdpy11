from cython.view cimport array as cvarray

cpdef double mean_s(const double[:] arr):
    cdef double rv = 0.0
    cdef size_t i = 0
    cdef size_t alen = arr.shape[0]
    while i < alen:
        rv += arr[i]
        i+=1
    return rv/<double>(alen)
    

