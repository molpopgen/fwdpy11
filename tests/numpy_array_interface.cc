#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <fwdpy11/numpy/array.hpp>

namespace py = pybind11;

struct VecBackedArray
{
    std::vector<double> x;
    VecBackedArray() : x(std::vector<double>(5, 1.0)) {}
};

struct VecBacked2DArray
{
    std::size_t nrow, ncol;
    std::vector<double> x;
    VecBacked2DArray(std::size_t nr, std::size_t nc)
        : nrow(nr), ncol(nc), x(std::vector<double>(ncol * nrow, 0.0))
    {
    }
};

PYBIND11_MODULE(numpy_array_interface, m)
{
    py::class_<VecBackedArray>(m, "VecBackedArray")
        .def(py::init<>())
        .def("x_readwrite",
             [](const VecBackedArray& self) {
                 return fwdpy11::make_1d_ndarray(self.x);
             })
        .def("x_readonly",
             [](const VecBackedArray& self) {
                 return fwdpy11::make_1d_ndarray_readonly(self.x);
             })
        .def("x_readwrite_via_capsule", [](VecBackedArray& self) {
            return fwdpy11::make_1d_array_with_capsule(std::move(self.x));
        });

    py::class_<VecBacked2DArray>(m, "VecBacked2DArray")
        .def(py::init<std::size_t, std::size_t>())
        .def("x_readwrite",
             [](const VecBacked2DArray& self) {
                 return fwdpy11::make_2d_ndarray(self.x, self.nrow, self.ncol);
             })
        .def("x_readonly",
             [](const VecBacked2DArray& self) {
                 return fwdpy11::make_2d_ndarray_readonly(self.x, self.nrow,
                                                          self.ncol);
             })
        .def_property_readonly("shape", [](const VecBacked2DArray& self) {
            return py::make_tuple(self.nrow, self.ncol);
        });
}

