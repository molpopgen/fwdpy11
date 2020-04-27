#include <fwdpp/ts/site.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <sstream>

namespace py = pybind11;

void
init_ts_Site(py::module& m)
{
    PYBIND11_NUMPY_DTYPE(fwdpp::ts::site, position, ancestral_state);

    py::class_<fwdpp::ts::site>(m, "Site",
                                R"delim(
        A site in a genome

        .. versionadded:: 0.5.0
        )delim")
        .def_readonly("position", &fwdpp::ts::site::position,
                      "The position of the site in the genome")
        .def_readonly("ancestral_state", &fwdpp::ts::site::ancestral_state,
                      "The ancestral state of the site")
        .def("__repr__",
             [](const fwdpp::ts::site& self) {
                 std::ostringstream out;
                 out << "Site(position = " << self.position
                     << ", ancestral_state = "
                     << static_cast<int>(self.ancestral_state) << ')';
                 return out.str();
             })
        .def(py::pickle(
            [](const fwdpp::ts::site& self) {
                return py::make_tuple(self.position, self.ancestral_state);
            },
            [](py::tuple t) {
                return fwdpp::ts::site{ t[0].cast<double>(),
                                        t[1].cast<std::int8_t>() };
            }));
}

