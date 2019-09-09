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

PYBIND11_MODULE(numpy_array_interface, m)
{
    py::class_<VecBackedArray>(m, "VecBackedArray")
        .def(py::init<>())
        .def("x_readwrite", [](const VecBackedArray& self) {
            return fwdpy11::make_1d_ndarray(self.x);
        })
        .def("x_readonly", [](const VecBackedArray& self) {
            return fwdpy11::make_1d_ndarray_readonly(self.x);
        });

}

