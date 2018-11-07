#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

#include <fwdpy11/evolvets/samplerecorder.hpp>
#include <fwdpy11/evolvets/recorders.hpp>

namespace py = pybind11;

PYBIND11_MODULE(tsrecorders, m)
{
    m.doc() = "Classes for recording ancient samples.";

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

    py::class_<fwdpy11::no_ancient_samples>(
        m, "NoAncientSamples",
        "A recorder for tree sequence simulations that does nothing.")
        .def(py::init<>())
        .def("__call__",
             [](fwdpy11::no_ancient_samples& na, const fwdpy11::SlocusPop& pop,
                fwdpy11::samplerecorder& sr) { na(pop, sr); })
        .def("__call__",
             [](fwdpy11::no_ancient_samples& na, const fwdpy11::MlocusPop& pop,
                fwdpy11::samplerecorder& sr) { na(pop, sr); });

    py::class_<fwdpy11::random_ancient_samples>(
        m, "RandomAncientSamples",
        "Preserve random samples of individuals at predetermined time points.")
        .def(py::init<std::uint32_t, fwdpp::uint_t,
                      std::vector<fwdpp::uint_t>>(),
             py::arg("seed"), py::arg("samplesize"), py::arg("timepoints"))
        .def("__call__", [](fwdpy11::random_ancient_samples& na,
                            const fwdpy11::SlocusPop& pop,
                            fwdpy11::samplerecorder& sr) { na(pop, sr); })
        .def("__call__", [](fwdpy11::random_ancient_samples& na,
                            const fwdpy11::MlocusPop& pop,
                            fwdpy11::samplerecorder& sr) { na(pop, sr); });
}

