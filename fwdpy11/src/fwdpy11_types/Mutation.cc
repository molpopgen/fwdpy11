#include <fwdpy11/types/Mutation.hpp>
#include <pybind11/stl.h>
#include <fwdpy11/numpy/array.hpp>

namespace py = pybind11;

void
init_Mutation(py::module &m)
{
    // Sugar types
    py::class_<fwdpy11::Mutation, fwdpp::mutation_base>(
        m, "Mutation", "Mutation with effect size and dominance")
        .def(py::init<double, double, double, unsigned, std::uint16_t>(),
             py::arg("pos"), py::arg("s"), py::arg("h"), py::arg("g"),
             py::arg("label"),
             R"delim(
                Construct a mutations.

                :param pos: Mutation position (float)
                :param s: Effect size (float)
                :param h: Dominance term (float)
                :param g: Origin time (unsigned integer)
                :param label: Label (16 bit integer)

                .. testcode::

                    import fwdpy11
                    m = fwdpy11.Mutation(1.0, -1.0, 0.25, 0, 0)
                    print(m.pos)
                    print(m.s)
                    print(m.h)
                    print(m.g)
                    print(m.label)

                .. testoutput::
                    
                    1.0
                    -1.0
                    0.25
                    0
                    0
                )delim")
        .def(py::init([](double pos, double s, double h, fwdpp::uint_t g,
                         py::list esizes, py::list heffects,
                         std::uint16_t label) {
                 std::vector<double> esizes_;
                 std::vector<double> heffects_;
                 for (auto i : esizes)
                     {
                         esizes_.push_back(i.cast<double>());
                     }
                 for (auto i : heffects)
                     {
                         heffects_.push_back(i.cast<double>());
                     }
                 return fwdpy11::Mutation(pos, s, h, g, std::move(esizes_),
                                          std::move(heffects_), label);
             }),
             py::arg("pos"), py::arg("s"), py::arg("h"), py::arg("g"),
             py::arg("esizes"), py::arg("heffects"), py::arg("label"),
             R"delim(
			 Construct a mutation with both a constant effect and/or
			 a vector of effects.

			 :param pos: Mutation position (float)
			 :param s: Effect size (float)
			 :param h: Dominance term (float)
			 :param g: Origin time (unsigned integer)
			 :param esizes: List of effect sizes (list of float)
			 :param heffects: List of heterozygous effects (list of float)
			 :param label: Label (16 bit integer)
				
			 .. versionadded:: 0.2.0

			 )delim")
        .def(py::init<fwdpy11::Mutation::constructor_tuple>(),
             R"delim(
                Construct mutation from a tuple.

                The tuple should contain (pos, s, h, g, label)

                .. testcode::

                    import fwdpy11
                    m = fwdpy11.Mutation((1.0, -1.0, 0.25, 0, 0))
                    print(m.pos)
                    print(m.s)
                    print(m.h)
                    print(m.g)
                    print(m.label)

                .. testoutput::
                    
                    1.0
                    -1.0
                    0.25
                    0
                    0
                )delim")
        .def_readonly(
            "g", &fwdpy11::Mutation::g,
            "Generation when mutation arose (origination time). (read-only)")
        .def_readonly("s", &fwdpy11::Mutation::s,
                      "Selection coefficient/effect size. (read-only)")
        .def_readonly("h", &fwdpy11::Mutation::h,
                      "Dominance/effect in heterozygotes. (read-only)")
        .def_property_readonly("heffects",
                               [](const fwdpy11::Mutation &self) {
                                   return fwdpy11::make_1d_ndarray_readonly(
                                       self.heffects);
                               },
                               R"delim(
				Vector of heterozygous effects.

				.. versionadded:: 0.2.0

                .. versionchanged:: 0.4.0

                    Property is now a readonly numpy.ndarray
				)delim")
        .def_property_readonly("esizes",
                               [](const fwdpy11::Mutation &self) {
                                   return fwdpy11::make_1d_ndarray_readonly(
                                       self.esizes);
                               },
                               R"delim(
				Vector of effect sizes.

				.. versionadded:: 0.2.0

                .. versionchanged:: 0.4.0

                    Property is now a readonly numpy.ndarray
				)delim")
        .def_property_readonly(
            "key",
            [](const fwdpy11::Mutation &m) {
                return py::make_tuple(m.pos, m.s, m.g);
            },
            R"delim(It is often useful to have a unique key for
                    tracking mutations.  This property returns 
                    the tuple (pos, esize, origin).

                    .. versionadded:: 0.1.3.a1
                   )delim")
        .def(py::pickle(
            [](const fwdpy11::Mutation &m) {
                return py::make_tuple(m.pos, m.s, m.h, m.g, m.esizes,
                                      m.heffects, m.xtra);
            },
            [](py::tuple p) {
                return std::unique_ptr<fwdpy11::Mutation>(
                    new fwdpy11::Mutation(
                        p[0].cast<double>(), p[1].cast<double>(),
                        p[2].cast<double>(), p[3].cast<unsigned>(),
                        p[4].cast<std::vector<double>>(),
                        p[5].cast<std::vector<double>>(),
                        p[6].cast<std::uint16_t>()));
            }))
        .def("__str__",
             [](const fwdpy11::Mutation &m) {
                 return "Mutation[" + std::to_string(m.pos) + ","
                        + std::to_string(m.s) + "," + std::to_string(m.h) + ","
                        + std::to_string(m.g) + "," + std::to_string(m.xtra)
                        + "]";
             })
        .def("__eq__", [](const fwdpy11::Mutation &a,
                          const fwdpy11::Mutation &b) { return a == b; });
}

