#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <fwdpy11/evolvets/SampleRecorder.hpp>
#include <fwdpy11/evolvets/recorders.hpp>
#include <fwdpy11/numpy/array.hpp>

namespace py = pybind11;

void
init_tsrecorders(py::module& m)
{
    py::class_<fwdpy11::SampleRecorder>(
        m, "SampleRecorder",
        "Allow recording of ancient samples during simulations with tree "
        "sequences.")
        .def(py::init<>())
        .def_property_readonly("samples",
                               [](const fwdpy11::SampleRecorder& self) {
                                   return fwdpy11::make_1d_ndarray(
                                       self.samples);
                               },
                               "Access to samples. For unit-testing purposes")
        .def("add_sample", &fwdpy11::SampleRecorder::add_sample,
             py::arg("individual"),
             "Add the index of an individual to the list of samples")
        .def("assign", &fwdpy11::SampleRecorder::assign, py::arg("samples"),
             "Add a list of individuals to the list of samples.  Input is a "
             "numpy array with dtype np.uint32");

    py::class_<fwdpy11::no_ancient_samples>(
        m, "NoAncientSamples",
        "A recorder for tree sequence simulations that does nothing.")
        .def(py::init<>())
        .def("__call__", [](fwdpy11::no_ancient_samples& na,
                            const fwdpy11::DiploidPopulation& pop,
                            fwdpy11::SampleRecorder& sr) { na(pop, sr); });

    py::class_<fwdpy11::random_ancient_samples>(
        m, "RandomAncientSamples",
        "Preserve random samples of individuals at predetermined time points.")
        .def(py::init([](std::uint32_t seed, fwdpp::uint_t samplesize,
                         py::array_t<fwdpp::uint_t> timepoints) {
                 auto r = timepoints.unchecked<1>();
                 std::vector<fwdpp::uint_t> tp;
                 for (py::ssize_t i = 0; i < r.shape(0); ++i)
                     {
                         tp.push_back(r(i));
                     };
                 return fwdpy11::random_ancient_samples(seed, samplesize,
                                                        std::move(tp));
             }),
             py::arg("seed"), py::arg("samplesize"), py::arg("timepoints"))
        .def("__call__", [](fwdpy11::random_ancient_samples& na,
                            const fwdpy11::DiploidPopulation& pop,
                            fwdpy11::SampleRecorder& sr) { na(pop, sr); });
}

