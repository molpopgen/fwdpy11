//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/types/diploid.hpp>

namespace py = pybind11;

struct diploid_traits
{
    double g, e, w;
};

struct diploid_gametes
{
    std::size_t locus, first, second;
};

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::diploid_t>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<fwdpy11::diploid_t>>);
PYBIND11_MAKE_OPAQUE(std::vector<diploid_traits>);
PYBIND11_MAKE_OPAQUE(std::vector<diploid_gametes>);

PYBIND11_MODULE(_opaque_diploids, m)
{
    m.doc() = "Expose C++ containers of diploids to Python without copies.";

    py::bind_vector<fwdpy11::dipvector_t>(
        m, "VecDiploid", py::module_local(false),
        "C++ representation of a list of "
        ":class:`fwdpy11."
        "SingleLocusDiploid`.  Typically, access will be read-only.")
        .def("trait_array",
             [](const fwdpy11::dipvector_t& diploids) {
                 std::vector<diploid_traits> rv;
                 rv.reserve(diploids.size());
                 for (auto&& dip : diploids)
                     {
                         rv.push_back(diploid_traits{ dip.g, dip.e, dip.w });
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("trait_array",
             [](const fwdpy11::dipvector_t& diploids,
                py::array_t<std::size_t> individuals) {
                 auto r = individuals.unchecked<1>();
                 std::vector<diploid_traits> rv;
                 rv.reserve(r.shape(0));
                 for (decltype(r.shape(0)) i = 0; i < r.shape(0); ++i)
                     {
                         // range-check here
                         auto&& dip = diploids.at(r(i));
                         rv.push_back(diploid_traits{ dip.g, dip.e, dip.w });
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("trait_array",
             [](const fwdpy11::dipvector_t& diploids, py::slice slice) {
                 size_t start, stop, step, slicelength;

                 if (!slice.compute(diploids.size(), &start, &stop, &step,
                                    &slicelength))
                     throw py::error_already_set();

                 std::vector<diploid_traits> rv;
                 rv.reserve(slicelength);
                 for (size_t i = 0; i < slicelength; ++i)
                     {
                         auto&& dip = diploids.at(start);
                         rv.push_back(diploid_traits{ dip.g, dip.e, dip.w });
                         start += step;
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const fwdpy11::dipvector_t& diploids) {
                 std::vector<diploid_gametes> rv;
                 rv.reserve(diploids.size());
                 for (auto&& dip : diploids)
                     {
                         rv.push_back(
                             diploid_gametes{ 0, dip.first, dip.second });
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const fwdpy11::dipvector_t& diploids,
                py::array_t<std::size_t> individuals) {
                 auto r = individuals.unchecked<1>();
                 std::vector<diploid_gametes> rv;
                 rv.reserve(r.shape(0));
                 for (decltype(r.shape(0)) i = 0; i < r.shape(0); ++i)
                     {
                         rv.push_back(diploid_gametes{
                             0, diploids.at(i).first, diploids.at(i).second });
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const fwdpy11::dipvector_t& diploids, py::slice slice) {
                 size_t start, stop, step, slicelength;

                 if (!slice.compute(diploids.size(), &start, &stop, &step,
                                    &slicelength))
                     throw py::error_already_set();

                 std::vector<diploid_gametes> rv;
                 rv.reserve(slicelength);
                 for (size_t i = 0; i < slicelength; ++i)
                     {
                         rv.push_back(
                             diploid_gametes{ 0, diploids.at(start).first,
                                              diploids.at(start).second });
                         start += step;
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def(py::pickle(
            [](const std::vector<fwdpy11::diploid_t>& v) -> py::list {
                py::list rv;
                for (auto&& vi : v)
                    {
                        rv.append(vi);
                    }
                return rv;
            },
            [](py::list l) {
                std::vector<fwdpy11::diploid_t> rv;
                for (auto&& i : l)
                    {
                        rv.push_back(i.cast<fwdpy11::diploid_t>());
                    }
                return rv;
            }));

    py::bind_vector<std::vector<fwdpy11::dipvector_t>>(
        m, "VecVecDiploid", py::module_local(false),
        "Vector of "
        ":class:`fwdpy11.SingleLocusDiploid`.")
        .def("trait_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids) {
                 std::vector<diploid_traits> rv;
                 rv.reserve(diploids.size());
                 for (auto&& dip : diploids)
                     {
                         rv.push_back(diploid_traits{ dip.at(0).g, dip.at(0).e,
                                                      dip.at(0).w });
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids) {
                 std::vector<diploid_gametes> rv;
                 std::size_t locus;
                 for (auto&& dip : diploids)
                     {
                         locus = 0;
                         for (auto&& di : dip)
                             {
                                 rv.push_back(diploid_gametes{ locus, di.first,
                                                               di.second });
                             }
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("trait_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids,
                py::array_t<std::size_t> individuals) {
                 auto r = individuals.unchecked<1>();
                 std::vector<diploid_traits> rv;
                 rv.reserve(r.shape(0));
                 for (decltype(r.shape(0)) i = 0; i < r.shape(0); ++i)
                     {
                         // range-check here
                         auto&& dip = diploids.at(r(i)).at(0);
                         rv.push_back(diploid_traits{ dip.g, dip.e, dip.w });
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids,
                py::array_t<std::size_t> individuals) {
                 auto r = individuals.unchecked<1>();
                 std::vector<diploid_gametes> rv;
                 rv.reserve(r.shape(0));
                 std::size_t locus;
                 for (decltype(r.shape(0)) i = 0; i < r.shape(0); ++i)
                     {
                         auto&& dip = diploids.at(r(i));
                         locus = 0;
                         for (auto&& di : dip)
                             {
                                 rv.emplace_back(diploid_gametes{
                                     locus++, di.first, di.second });
                             }
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("trait_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids,
                py::slice slice) {
                 size_t start, stop, step, slicelength;

                 if (!slice.compute(diploids.size(), &start, &stop, &step,
                                    &slicelength))
                     throw py::error_already_set();

                 std::vector<diploid_traits> rv;
                 rv.reserve(slicelength);
                 for (size_t i = 0; i < slicelength; ++i)
                     {
                         auto&& dip = diploids.at(start).at(0);
                         rv.push_back(diploid_traits{ dip.g, dip.e, dip.w });
                         start += step;
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidTraits`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def("key_array",
             [](const std::vector<fwdpy11::dipvector_t>& diploids,
                py::slice slice) {
                 size_t start, stop, step, slicelength;

                 if (!slice.compute(diploids.size(), &start, &stop, &step,
                                    &slicelength))
                     throw py::error_already_set();

                 std::vector<diploid_gametes> rv;
                 rv.reserve(slicelength);
                 std::size_t locus;
                 for (size_t i = 0; i < slicelength; ++i)
                     {
                         auto&& dip = diploids.at(start);
                         locus = 0;
                         for (auto&& di : dip)
                             {
                                 rv.push_back(diploid_gametes{
                                     locus++, di.first, di.second });
                             }
                         start += step;
                     }
                 return rv;
             },
             R"delim(
			 :rtype: :class:`fwdpy11.VecDiploidGametes`
			 
			 .. versionadded:: 0.1.2
			 )delim")
        .def(py::pickle(
            [](const std::vector<fwdpy11::dipvector_t>& diploids) {
                py::list rv;
                for (auto&& i : diploids)
                    {
                        rv.append(i);
                    };
                return rv;
            },
            [](py::list l) {
                std::vector<fwdpy11::dipvector_t> rv;
                for (auto&& i : l)
                    {
                        rv.push_back(i.cast<fwdpy11::dipvector_t>());
                    }
                return rv;
            }));

    PYBIND11_NUMPY_DTYPE(diploid_traits, g, e, w);
    PYBIND11_NUMPY_DTYPE(diploid_gametes, locus, first, second);

    py::bind_vector<std::vector<diploid_traits>>(
        m, "VecDiploidTraits", py::buffer_protocol(), py::module_local(false),
        R"delim(
        Vector of the g,e,w data fields in a "
        ":class:`fwdpy11.SingleLocusDiploid`.

        .. versionadded: 0.1.2
        )delim");

    py::bind_vector<std::vector<diploid_gametes>>(
        m, "VecDiploidGametes", py::buffer_protocol(), py::module_local(false),
        R"delim(
        Vector of the locus and gamete index data fields in a "
        ":class:`fwdpy11.SingleLocusDiploid`.

        .. versionadded: 0.1.2
        )delim");
}
