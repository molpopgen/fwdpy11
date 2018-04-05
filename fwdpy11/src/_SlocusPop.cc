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
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/types/create_pops.hpp>
#include <fwdpy11/serialization.hpp>
#include <fwdpy11/serialization/Diploid.hpp>

namespace py = pybind11;

namespace
{
    static const auto DIPLOIDS_DOCSTRING = R"delim(
   A :class:`fwdpy11.VecDiploid`.
   )delim";
}

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Diploid>);

PYBIND11_MODULE(_SlocusPop, m)
{
    py::object base_class_module
        = (pybind11::object)pybind11::module::import("fwdpy11._Population");

    py::class_<fwdpy11::SlocusPop, fwdpy11::Population>(
        m, "_SlocusPop", "Representation of a single-locus population")
        .def(py::init<fwdpp::uint_t>(), "Construct with an unsigned integer "
                                        "representing the initial "
                                        "population size.")
        .def(py::init<const fwdpy11::SlocusPop::dipvector_t&,
                      const fwdpy11::SlocusPop::gcont_t&,
                      const fwdpy11::SlocusPop::mcont_t&>(),
             R"delim(
             Construct with tuple of (diploids, gametes, mutations).
             
             .. versionadded:: 0.1.4
             )delim")
        .def(py::init<const fwdpy11::SlocusPop&>(),
             R"delim(
                Copy constructor

                .. versionadded:: 0.1.4
                )delim")
        .def("clear", &fwdpy11::SlocusPop::clear,
             "Clears all population data.")
        .def("__eq__",
             [](const fwdpy11::SlocusPop& lhs, const fwdpy11::SlocusPop& rhs) {
                 return lhs == rhs;
             })
        .def_readonly("diploids", &fwdpy11::SlocusPop::diploids,
                      DIPLOIDS_DOCSTRING)
        .def_static(
            "create",
            [](fwdpy11::SlocusPop::dipvector_t& diploids,
               fwdpy11::SlocusPop::gcont_t& gametes,
               fwdpy11::SlocusPop::mcont_t& mutations,
               py::tuple args) -> fwdpy11::SlocusPop {
                if (args.size() == 0)
                    {
                        return fwdpy11::create_wrapper<fwdpy11::SlocusPop>()(
                            diploids, gametes, mutations);
                    }
                auto fixations = args[0].cast<fwdpy11::SlocusPop::mcont_t>();
                auto ftimes = args[1].cast<std::vector<fwdpp::uint_t>>();
                auto g = args[2].cast<fwdpp::uint_t>();
                return fwdpy11::create_wrapper<fwdpy11::SlocusPop>()(
                    diploids, gametes, mutations, fixations, ftimes, g);
            })
        .def(py::pickle(
            [](const fwdpy11::SlocusPop& pop) -> py::object {
                auto pb = py::bytes(
                    fwdpy11::serialization::serialize_details(&pop));
                py::list pdata;
                for (auto& d : pop.diploids)
                    {
                        pdata.append(d.parental_data);
                    }
                return py::make_tuple(std::move(pb), std::move(pdata),
                                      pop.popdata, pop.popdata_user);
            },
            [](py::object pickled) -> fwdpy11::SlocusPop {
                try
                    {
                        auto s = pickled.cast<py::bytes>();
                        return fwdpy11::serialization::deserialize_details<
                            fwdpy11::SlocusPop>()(s, 1);
                    }
                catch (std::runtime_error& eas)
                    {
                        PyErr_Clear();
                    }
                auto t = pickled.cast<py::tuple>();
                if (t.size() != 4)
                    {
                        throw std::runtime_error(
                            "expected tuple with 4 elements");
                    }
                auto s = t[0].cast<py::bytes>();
                auto l = t[1].cast<py::list>();
                auto rv = fwdpy11::serialization::deserialize_details<
                    fwdpy11::SlocusPop>()(s, 1);
                for (std::size_t i = 0; i < rv.diploids.size(); ++i)
                    {
                        rv.diploids[i].parental_data = l[i];
                    }
                rv.popdata = t[2];
                rv.popdata_user = t[3];
                return rv;
            }));
}
