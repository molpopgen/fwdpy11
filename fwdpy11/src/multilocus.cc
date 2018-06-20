#include <exception>
#include <functional>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/interlocus_recombination.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/multilocus.hpp>
#include <fwdpy11/fitness/fitness.hpp>

namespace py = pybind11;

PYBIND11_MODULE(multilocus, m)
{
    m.doc() = "Types and functions specific to multi-locus/region simulations";

    py::class_<fwdpy11::interlocus_rec>(
        m, "InterlocusRecombination",
        "This class parameterizes interlocus recombination in multilocus "
        "simulations.  You do not make instances of this class directly. "
        " Rather, you make calls to one of "
        ":func:`fwdpy11.multilocus.poisson_rec` or "
        ":func:`fwdpy11.multilocus.binomial_rec`.")
        .def(py::init<double, std::underlying_type<
                                  fwdpy11::interlocus_rec::RECMODEL>::type>())
        .def(py::pickle(
            [](const fwdpy11::interlocus_rec& ir) {
                return py::make_tuple(ir.param, ir.get_model());
            },
            [](py::tuple t) {
                double d = t[0].cast<double>();
                auto m = t[1].cast<fwdpy11::interlocus_rec::mtype>();
                return std::unique_ptr<fwdpy11::interlocus_rec>(
                    new fwdpy11::interlocus_rec(d, m));
                // new (&ir) fwdpy11::interlocus_rec(d, m);
            }))
        .def("__repr__", [](const fwdpy11::interlocus_rec& ir) {
            std::string rv = "multilocus.InterlocusRecombination(";
            rv += std::to_string(ir.param);
            rv += ',';
            rv += std::to_string(ir.m);
            rv += ')';
            return rv;
        });

    m.def("poisson_rec",
          [](const std::vector<double>& means) {
              py::list rv;
              for (auto&& r : means)
                  {
                      rv.append(fwdpy11::interlocus_rec(
                          r, fwdpy11::interlocus_rec::RECMODEL::POISSON));
                  }
              return rv;
          },
          R"delim(
    Parameterize interlocus recomination as a Poisson process.

    :param means: A list of mean values.

    :rtype: list

    :return: A list of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. versionchanged:: 0.1.3
        No longer takes a :class:`fwdpy11.GSLrng` as argument.
    )delim");

    m.def("poisson_rec",
          [](const double mean) {
              return fwdpy11::interlocus_rec(
                  mean, fwdpy11::interlocus_rec::RECMODEL::POISSON);
          },
          R"delim(
    Parameterize interlocus recomination as a Poisson process.

    :param mean: The mean of a Poisson process.

    :rtype: :class:`fwdpy11.multilocus.InterlocusRecombination`

    :return: An instance of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. versionchanged:: 0.1.3
        No longer takes a :class:`fwdpy11.GSLrng` as argument.
    )delim");

    m.def("binomial_rec",
          [](const std::vector<double>& probs) {
              py::list rv;
              for (auto&& r : probs)
                  {
                      rv.append(fwdpy11::interlocus_rec(
                          r, fwdpy11::interlocus_rec::RECMODEL::BINOMIAL));
                  }
              return rv;
          },
          R"delim(
    Parameterize interlocus recomination as a Binomial process.

    :param probs: A list of genetic distance in cM/100.

    :rtype: list

    :return: A list of of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. versionchanged:: 0.1.3
        No longer takes a :class:`fwdpy11.GSLrng` as argument.
    )delim");

    m.def("binomial_rec",
          [](const double prob) {
              return fwdpy11::interlocus_rec(
                  prob, fwdpy11::interlocus_rec::RECMODEL::BINOMIAL);
          },
          R"delim(
    Parameterize interlocus recomination as a Binomial process.

    :param prob: The genetic distance in cM/100.

    :rtype: :class:`fwdpy11.multilocus.InterlocusRecombination`

    :return: An instance of :class:`fwdpy11.multilocus.InterlocusRecombination`.

    .. versionchanged:: 0.1.3
        No longer takes a :class:`fwdpy11.GSLrng` as argument.
    )delim");
}
