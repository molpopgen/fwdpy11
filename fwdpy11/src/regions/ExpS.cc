#include <pybind11/pybind11.h>
#include <fwdpy11/regions/ExpS.hpp>
#include <fwdpy11/mutation_dominance/MutationDominance.hpp>

namespace py = pybind11;

void
init_ExpS(py::module& m)
{
    py::class_<fwdpy11::ExpS, fwdpy11::Sregion>(m, "_ll_ExpS")
        .def(py::init([](double beg, double end, double weight, double mean,
                         const double h, bool coupled, std::uint16_t label,
                         double scaling) {
                 return fwdpy11::ExpS(fwdpy11::Region(beg, end, weight, coupled, label),
                                      scaling, mean, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("mean"),
             py::arg("h"), py::arg("coupled"), py::arg("label"), py::arg("scaling"))
        .def(py::init([](double beg, double end, double weight, double mean,
                         const fwdpy11::MutationDominance& h, bool coupled,
                         std::uint16_t label, double scaling) {
                 return fwdpy11::ExpS(fwdpy11::Region(beg, end, weight, coupled, label),
                                      scaling, mean, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("mean"),
             py::arg("h"), py::arg("coupled"), py::arg("label"), py::arg("scaling"));
}
