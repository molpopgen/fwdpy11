#include <pybind11/pybind11.h>
#include <fwdpy11/regions/GammaS.hpp>
#include <fwdpy11/mutation_dominance/MutationDominance.hpp>

namespace py = pybind11;

void
init_GammaS(py::module& m)
{
    py::class_<fwdpy11::GammaS, fwdpy11::Sregion>(m, "_ll_GammaS")
        .def(
            py::init([](double beg, double end, double weight, double mean, double shape,
                        double h, bool coupled, std::uint16_t label, double scaling) {
                return fwdpy11::GammaS(fwdpy11::Region(beg, end, weight, coupled, label),
                                       scaling, mean, shape, h);
            }),
            py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("mean"),
            py::arg("shape_parameter"), py::arg("h"), py::arg("coupled"),
            py::arg("label"), py::arg("scaling"))
        .def(py::init([](double beg, double end, double weight, double mean,
                         double shape, const fwdpy11::MutationDominance& h, bool coupled,
                         std::uint16_t label, double scaling) {
                 return fwdpy11::GammaS(
                     fwdpy11::Region(beg, end, weight, coupled, label), scaling, mean,
                     shape, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("mean"),
             py::arg("shape_parameter"), py::arg("h"), py::arg("coupled"),
             py::arg("label"), py::arg("scaling"));
}
