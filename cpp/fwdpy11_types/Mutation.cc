#include <algorithm>
#include <cmath>
#include <fwdpy11/types/Mutation.hpp>
#include <pybind11/stl.h>
#include <fwdpy11/numpy/array.hpp>

namespace py = pybind11;

static const auto INIT_DOCSTRING_1 =
    R"delim(
Construct a mutations.

:param pos: Mutation position (float)
:param s: Effect size (float)
:param h: Dominance term (float)
:param g: Origin time (signed integer)
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
)delim";

static const auto INIT_DOCSTRING_2 =
    R"delim(
Construct a mutation with both a constant effect and/or
a vector of effects.

:param pos: Mutation position (float)
:param s: Effect size (float)
:param h: Dominance term (float)
:param g: Origin time (signed integer)
:param esizes: List of effect sizes (list of float)
:param heffects: List of heterozygous effects (list of float)
:param label: Label (16 bit integer)

.. versionadded:: 0.2.0

)delim";

void
init_Mutation(py::module &m)
{
    // Sugar types
    py::class_<fwdpy11::Mutation, fwdpp::mutation_base>(
        m, "Mutation", "Mutation with effect size and dominance")
        .def(py::init([](double pos, double s, double h, fwdpy11::mutation_origin_time g,
                         std::uint16_t label) {
                 if (!std::isfinite(s))
                     {
                         throw std::invalid_argument("effect size must be finite");
                     }
                 if (!std::isfinite(h))
                     {
                         throw std::invalid_argument("dominance must be finite");
                     }
                 return fwdpy11::Mutation(s == 0., pos, s, h, g, label);
                 return fwdpy11::Mutation(s == 0., pos, s, h, g, label);
             }),
             py::arg("pos"), py::arg("s"), py::arg("h"), py::arg("g"), py::arg("label"),
             INIT_DOCSTRING_1)
        .def(py::init([](double pos, double s, double h, fwdpy11::mutation_origin_time g,
                         py::list esizes, py::list heffects, std::uint16_t label) {
                 if (!std::isfinite(s))
                     {
                         throw std::invalid_argument("effect size must be finite");
                     }
                 if (!std::isfinite(h))
                     {
                         throw std::invalid_argument("dominance must be finite");
                     }
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
                 bool neutral = (s == 0.);
                 if (!esizes_.empty() && neutral == true)
                     {
                         neutral = std::all_of(begin(esizes_), end(esizes_),
                                               [](double d) { return d == 0.; });
                     }
                 return fwdpy11::Mutation(neutral, pos, s, h, g, std::move(esizes_),
                                          std::move(heffects_), label);
             }),
             py::arg("pos"), py::arg("s"), py::arg("h"), py::arg("g"), py::arg("esizes"),
             py::arg("heffects"), py::arg("label"), INIT_DOCSTRING_2)
        .def_readonly("g", &fwdpy11::Mutation::g,
                      "Generation when mutation arose (origination time). (read-only)")
        .def_readonly("s", &fwdpy11::Mutation::s,
                      "Selection coefficient/effect size. (read-only)")
        .def_readonly("h", &fwdpy11::Mutation::h,
                      "Dominance/effect in heterozygotes. (read-only)")
        .def_property_readonly(
            "heffects",
            [](const fwdpy11::Mutation &self) {
                return fwdpy11::make_1d_ndarray_readonly(self.heffects);
            },
            R"delim(
				Vector of heterozygous effects.

				.. versionadded:: 0.2.0

                .. versionchanged:: 0.4.0

                    Property is now a readonly numpy.ndarray
				)delim")
        .def_property_readonly(
            "esizes",
            [](const fwdpy11::Mutation &self) {
                return fwdpy11::make_1d_ndarray_readonly(self.esizes);
            },
            R"delim(
				Vector of effect sizes.

				.. versionadded:: 0.2.0

                .. versionchanged:: 0.4.0

                    Property is now a readonly numpy.ndarray
				)delim")
        .def_property_readonly(
            "key",
            [](const fwdpy11::Mutation &m) { return py::make_tuple(m.pos, m.s, m.g); },
            R"delim(It is often useful to have a unique key for
                    tracking mutations.  This property returns 
                    the tuple (pos, esize, origin).

                    .. versionadded:: 0.1.3.a1
                   )delim")
        .def(py::pickle(
            [](const fwdpy11::Mutation &m) {
                return py::make_tuple(m.neutral, m.pos, m.s, m.h, m.g, m.esizes,
                                      m.heffects, m.xtra);
            },
            [](py::tuple p) {
                return std::unique_ptr<fwdpy11::Mutation>(new fwdpy11::Mutation(
                    p[0].cast<bool>(), p[1].cast<double>(), p[2].cast<double>(),
                    p[3].cast<double>(), p[4].cast<fwdpy11::mutation_origin_time>(),
                    p[5].cast<std::vector<double>>(), p[6].cast<std::vector<double>>(),
                    p[7].cast<std::uint16_t>()));
            }))
        .def("__str__",
             [](const fwdpy11::Mutation &self) {
                 auto rv = "Mutation[position:" + std::to_string(self.pos)
                           + ", effect size:" + std::to_string(self.s)
                           + ", dominance:" + std::to_string(self.h)
                           + ", origin time:" + std::to_string(self.g)
                           + ", label:" + std::to_string(self.xtra);
                 if (!self.esizes.empty())
                     {
                         rv += ", multivariate effect sizes: [";
                         for (auto i : self.esizes)
                             {
                                 rv += std::to_string(i);
                                 rv += ", ";
                             }
                         rv += "], multivariate dominance: [";
                         for (auto i : self.heffects)
                             {
                                 rv += std::to_string(i);
                                 rv += ", ";
                             }
                         rv += "]";
                     }

                 rv += "]";

                 return rv;
             })
        .def("__eq__", [](const fwdpy11::Mutation &a, const fwdpy11::Mutation &b) {
            return a == b;
        });
}
