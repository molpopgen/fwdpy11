#include <pybind11/pybind11.h>

#include <fwdpy11/evolvets/samplerecorder.hpp>

namespace py = pybind11;

PYBIND11_MODULE(_tsevolveutils, m)
{
    m.doc() = "Helper types for simulations with tree sequences.";

    py::class_<fwdpy11::samplerecorder>(
        m, "SampleRecorder",
        "Allow recording of ancient samples during simulations with tree "
        "sequences.")
        .def(py::init<>())
        .def_readonly("samples", &fwdpy11::samplerecorder::samples,
                      "Access to samples. For unit-testing purposes")
        .def("add_sample", &fwdpy11::samplerecorder::add_sample,
             py::arg("individual"),
             "Add the index of an individual to the list of samples")
        .def("assign", &fwdpy11::samplerecorder::assign, py::arg("samples"),
             "Add a list of individuals to the list of samples.  Input is a "
             "numpy array with dtype np.uint32");
}
