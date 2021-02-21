#include <pybind11/pybind11.h>
#include <fwdpy11/regions/LogNormalS.hpp>
#include <fwdpy11/mutation_dominance/MutationDominance.hpp>

namespace py = pybind11;

void
init_LogNormalS(py::module& m)
{
    py::class_<fwdpy11::LogNormalS, fwdpy11::Sregion>(m, "_ll_LogNormalS")
        .def(py::init([](double beg, double end, double weight, py::object zeta,
                         py::object sigma, double h, bool coupled, std::uint16_t label,
                         double scaling) {
                 if (zeta.is_none() == false && sigma.is_none() == false)
                     {
                         return fwdpy11::LogNormalS(
                             fwdpy11::Region(beg, end, weight, coupled, label), scaling,
                             zeta.cast<double>(), sigma.cast<double>(), h);
                     }
                 return fwdpy11 ::LogNormalS(
                     fwdpy11::Region(beg, end, weight, coupled, label), scaling, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("zeta"),
             py::arg("sigma"), py::arg("h"), py::arg("coupled"), py::arg("label"),
             py::arg("scaling"))
        .def(py::init([](double beg, double end, double weight, py::object zeta,
                         py::object sigma, const fwdpy11::MutationDominance& h,
                         bool coupled, std::uint16_t label, double scaling) {
                 if (zeta.is_none() == false && sigma.is_none() == false)
                     {
                         return fwdpy11::LogNormalS(
                             fwdpy11::Region(beg, end, weight, coupled, label), scaling,
                             zeta.cast<double>(), sigma.cast<double>(), h);
                     }
                 return fwdpy11 ::LogNormalS(
                     fwdpy11::Region(beg, end, weight, coupled, label), scaling, h);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("zeta"),
             py::arg("sigma"), py::arg("h"), py::arg("coupled"), py::arg("label"),
             py::arg("scaling"));
}
