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

#include <fwdpy11/discrete_demography/DiscreteDemography.hpp>
#include <fwdpy11/discrete_demography/simulation/demographic_model_state.hpp>

namespace py = pybind11;
namespace ddemog = fwdpy11::discrete_demography;

namespace
{
    ddemog::MigrationMatrix
    decode_migration_matrix_input(py::object o)
    {
        const auto a2v = [](py::array_t<double> m) {
            auto r = m.unchecked<2>();
            if (r.shape(0) != r.shape(1))
                {
                    throw std::invalid_argument("MigrationMatrix must be square");
                }
            auto d = r.data(0, 0);
            std::vector<double> M(d, d + r.shape(0) * r.shape(1));
            return M;
        };

        bool isarray = true;

        try
            {
                py::array_t<double> m = o.cast<py::array_t<double>>();
            }
        catch (...)
            {
                isarray = false;
            }
        if (isarray)
            {
                py::array_t<double> m = o.cast<py::array_t<double>>();
                auto M = a2v(m);
                return ddemog::MigrationMatrix(std::move(M), m.shape(0), false);
            }
        ddemog::MigrationMatrix M = o.cast<ddemog::MigrationMatrix>();
        return M;
    }
} // namespace

void
init_DiscreteDemography(py::module& m)
{
    py::class_<ddemog::DiscreteDemography>(m, "_ll_DiscreteDemography")
        .def(py::init([](py::object mass_migration_events, py::object set_growth_rates,
                         py::object set_deme_sizes, py::object set_selfing_rates,
                         py::object migmatrix, py::object set_migration_rates) {
                 ddemog::DiscreteDemography::mass_migration_vector morc;
                 ddemog::DiscreteDemography::set_growth_rates_vector growth;
                 ddemog::DiscreteDemography::set_deme_sizes_vector change_sizes;
                 ddemog::DiscreteDemography::set_selfing_rates_vector selfing;
                 ddemog::DiscreteDemography::set_migration_rates_vector migrates;

                 if (mass_migration_events.is_none() == false)
                     {
                         morc = mass_migration_events.cast<decltype(morc)>();
                     }
                 if (set_growth_rates.is_none() == false)
                     {
                         growth = set_growth_rates.cast<decltype(growth)>();
                     }
                 if (set_deme_sizes.is_none() == false)
                     {
                         change_sizes = set_deme_sizes.cast<decltype(change_sizes)>();
                     }
                 if (set_selfing_rates.is_none() == false)
                     {
                         selfing = set_selfing_rates.cast<decltype(selfing)>();
                     }
                 if (set_migration_rates.is_none() == false)
                     {
                         migrates = set_migration_rates.cast<decltype(migrates)>();
                     }
                 if (migmatrix.is_none() == true)
                     {
                         return ddemog::DiscreteDemography(
                             std::move(morc), std::move(growth), std::move(change_sizes),
                             std::move(selfing), nullptr, std::move(migrates));
                     }
                 std::unique_ptr<ddemog::MigrationMatrix> M(new ddemog::MigrationMatrix(
                     decode_migration_matrix_input(migmatrix)));
                 return ddemog::DiscreteDemography(
                     std::move(morc), std::move(growth), std::move(change_sizes),
                     std::move(selfing), std::move(M), std::move(migrates));
             }),
             py::arg("mass_migrations") = py::none(),
             py::arg("set_growth_rates") = py::none(),
             py::arg("set_deme_sizes") = py::none(),
             py::arg("set_selfing_rates") = py::none(),
             py::arg("migmatrix") = py::none(),
             py::arg("set_migration_rates") = py::none());
}
