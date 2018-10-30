#include <pybind11/pybind11.h>

#include <fwdpy11/evolvets/samplerecorder.hpp>

namespace py = pybind11;

PYBIND11_MODULE(_tsevolveutils, m)
{
    py::class_<fwdpy11::samplerecorder>(m, "SampleRecorder")
        .def(py::init<>())
        .def_readonly("samples",&fwdpy11::samplerecorder::samples)
        .def("add_sample", &fwdpy11::samplerecorder::add_sample)
        .def("assign", &fwdpy11::samplerecorder::assign);
}
