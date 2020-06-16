#include <pybind11/pybind11.h>
#include <fwdpy11/types/Diploid.hpp>
#include <sstream>
namespace py = pybind11;

void
init_DiploidMetadata(py::module &m)
{
    py::class_<fwdpy11::DiploidMetadata>(m, "DiploidMetadata",
                                         "Diploid meta data.")
        .def_readwrite("g", &fwdpy11::DiploidMetadata::g, "Genetic value.")
        .def_readwrite("e", &fwdpy11::DiploidMetadata::e,
                       "Random component of trait value.")
        .def_readwrite("w", &fwdpy11::DiploidMetadata::w, "Fitness.")
        .def_property(
            "geography",
            [](const fwdpy11::DiploidMetadata &d) {
                return py::make_tuple(d.geography[0], d.geography[1],
                                      d.geography[2]);
            },
            [](fwdpy11::DiploidMetadata &d,
               const std::tuple<double, double, double> &input) {
                d.geography[0] = std::get<0>(input);
                d.geography[1] = std::get<1>(input);
                d.geography[2] = std::get<2>(input);
            },
            "Array containing the geographic location of the individual.")
        .def_property("parents",
                      [](const fwdpy11::DiploidMetadata &d) {
                          return py::make_tuple(d.parents[0], d.parents[1]);
                      },
                      [](fwdpy11::DiploidMetadata &d,
                         const std::pair<std::size_t, std::size_t> &input) {
                          d.parents[0] = input.first;
                          d.parents[1] = input.second;
                      },
                      "Array containing the label fields of the parents.")
        .def_readwrite("sex", &fwdpy11::DiploidMetadata::sex, "Sex.")
        .def_readwrite("deme", &fwdpy11::DiploidMetadata::deme, "Deme.")
        .def_readwrite("label", &fwdpy11::DiploidMetadata::label,
                       "Index of the individual in the population.")
        .def_property_readonly("nodes",
                               [](const fwdpy11::DiploidMetadata &md) {
                                   py::list rv;
                                   rv.append(md.nodes[0]);
                                   rv.append(md.nodes[1]);
                                   return rv;
                               },
                               "Node ids for individual")
        .def("__repr__",
             [](const fwdpy11::DiploidMetadata &self) {
                 std::ostringstream out;
                 out.precision(4);
                 out << "DiploidMetadata("
                     << "g=" << self.g << ',' << "w=" << self.w << ','
                     << "e=" << self.e << ',' << "label=" << self.label << ','
                     << "nodes=[" << self.nodes[0] << ',' << self.nodes[1]
                     << "],"
                     << "parents=[" << self.parents[0] << ','
                     << self.parents[1] << "],"
                     << "sex=" << self.sex << ',' << "deme=" << self.deme
                     << ',' << "geography=[" << self.geography[0] << ','
                     << self.geography[1] << ',' << self.geography[2] << "]"
                     << ')';
                 return out.str();
             })
        .def(py::pickle(
            [](const fwdpy11::DiploidMetadata &md) {
                return py::make_tuple(
                    md.g, md.e, md.w,
                    py::make_tuple(md.geography[0], md.geography[1],
                                   md.geography[2]),
                    md.label, py::make_tuple(md.parents[0], md.parents[1]),
                    md.deme, md.sex, py::make_tuple(md.nodes[0], md.nodes[1]));
            },
            [](py::tuple t) {
                fwdpy11::DiploidMetadata rv;
                rv.g = t[0].cast<double>();
                rv.e = t[1].cast<double>();
                rv.w = t[2].cast<double>();
                auto ttuple = t[3].cast<py::tuple>();
                rv.geography[0] = ttuple[0].cast<double>();
                rv.geography[1] = ttuple[1].cast<double>();
                rv.geography[2] = ttuple[2].cast<double>();
                rv.label = t[4].cast<std::size_t>();
                ttuple = t[5].cast<py::tuple>();
                rv.parents[0] = ttuple[0].cast<std::size_t>();
                rv.parents[1] = ttuple[1].cast<std::size_t>();
                rv.deme = t[6].cast<std::int32_t>();
                rv.sex = t[7].cast<std::int32_t>();
                ttuple = t[8].cast<py::tuple>();
                rv.nodes[0] = ttuple[0].cast<fwdpp::ts::table_index_t>();
                rv.nodes[1] = ttuple[1].cast<fwdpp::ts::table_index_t>();
                return rv;
            }));
}
