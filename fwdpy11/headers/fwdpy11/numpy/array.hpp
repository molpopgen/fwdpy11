#ifndef FWDPY11_NUMPY_ARRAY_HPP
#define FWDPY11_NUMPY_ARRAY_HPP

#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace fwdpy11
{
    template <typename T>
    inline pybind11::array_t<T>
    make_1d_ndarray(const std::vector<T>& v)
    // Returns a 1d numpy array that does not own
    // its data.
    {
        auto rv = pybind11::array_t<T>({ v.size() }, { sizeof(T) }, v.data(),
                                       pybind11::cast(v));
        return rv;
    }

    template <typename T>
    inline pybind11::array_t<T>
    make_1d_ndarray_readonly(const std::vector<T>& v)
    // Returns a readonly 1d numpy array that does not own
    // its data.
    {
        auto rv = pybind11::array_t<T>({ v.size() }, { sizeof(T) }, v.data(),
                                       pybind11::cast(v));
        rv.attr("flags").attr("writeable") = false;
        return rv;
    }

    template <typename T>
    inline pybind11::array_t<T>
    make_2d_ndarray(const std::vector<T>& v, std::size_t dim1,
                    std::size_t dim2)
    // Returns a 2d numpy array that does not own
    // its data.
    {
        auto rv = pybind11::array_t<T>({ dim1, dim2 }, v.data(),
                                       pybind11::cast(v));
        return rv;
    }

    template <typename T>
    inline pybind11::array_t<T>
    make_2d_ndarray_readonly(const std::vector<T>& v, std::size_t dim1,
                             std::size_t dim2)
    // Returns a 2d numpy array that does not own
    // its data.
    {
        auto rv = pybind11::array_t<T>({ dim1, dim2 }, v.data(),
                                       pybind11::cast(v));
        rv.attr("flags").attr("writeable") = false;
        return rv;
    }
} // namespace fwdpy11

#endif
