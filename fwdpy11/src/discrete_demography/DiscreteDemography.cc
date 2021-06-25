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

    py::dict
    model_state_as_dict(const ddemog::DiscreteDemographyState& model_state)
    {
        py::dict rv;
        rv["maxdemes"] = model_state.maxdemes;

        // This is the deme_properties stuff
        rv["current_deme_sizes"]
            = model_state.current_deme_parameters.current_deme_sizes.get();
        rv["next_deme_sizes"]
            = model_state.current_deme_parameters.next_deme_sizes.get();
        rv["growth_rate_onset_times"]
            = model_state.current_deme_parameters.growth_rate_onset_times.get();
        rv["growth_initial_sizes"]
            = model_state.current_deme_parameters.growth_initial_sizes.get();
        rv["growth_rates"] = model_state.current_deme_parameters.growth_rates.get();
        rv["selfing_rates"] = model_state.current_deme_parameters.selfing_rates.get();

        rv["migmatrix"]
            = py::make_tuple(model_state.M.M, model_state.M.npops, model_state.M.scaled);

        rv["bookmark"] = py::dict{};
        rv["bookmark"]["starts"] = model_state.fitness_bookmark.starts;
        rv["bookmark"]["stops"] = model_state.fitness_bookmark.stops;
        rv["bookmark"]["offsets"] = model_state.fitness_bookmark.offsets;
        rv["bookmark"]["individuals"] = model_state.fitness_bookmark.individuals;
        rv["bookmark"]["individual_fitness"]
            = model_state.fitness_bookmark.individual_fitness;

        return rv;
    }
} // namespace

void
init_DiscreteDemography(py::module& m)
{
    // TODO: decide if we need to pass fitnesses via asdict?
    // py::class_<ddemog::demographic_model_state>(m, "_ll_DemographicModelState");
    py::class_<ddemog::DiscreteDemography>(m, "_ll_DiscreteDemography")
        .def(py::init([](py::object mass_migration_events, py::object set_growth_rates,
                         py::object set_deme_sizes, py::object set_selfing_rates,
                         py::object migmatrix, py::object set_migration_rates) {
                 std::vector<ddemog::MassMigration> morc;
                 std::vector<ddemog::SetExponentialGrowth> growth;
                 std::vector<ddemog::SetDemeSize> change_sizes;
                 std::vector<ddemog::SetSelfingRate> selfing;
                 std::vector<ddemog::SetMigrationRates> migrates;

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
                             std::move(selfing), ddemog::MigrationMatrix{},
                             std::move(migrates));
                     }
                 ddemog::MigrationMatrix M(decode_migration_matrix_input(migmatrix));
                 return ddemog::DiscreteDemography(
                     std::move(morc), std::move(growth), std::move(change_sizes),
                     std::move(selfing), std::move(M), std::move(migrates));
             }),
             py::arg("mass_migrations") = py::none(),
             py::arg("set_growth_rates") = py::none(),
             py::arg("set_deme_sizes") = py::none(),
             py::arg("set_selfing_rates") = py::none(),
             py::arg("migmatrix") = py::none(),
             py::arg("set_migration_rates") = py::none())
        .def("_clone_state_to", &ddemog::DiscreteDemography::copy_state_to)
        .def("_state_asdict",
             [](ddemog::DiscreteDemography& self) -> py::object {
                 return model_state_as_dict(self.get_model_state());
             })
        .def("_reset_state", [](ddemog::DiscreteDemography& self, py::object o) {
            auto model_state = self.get_model_state();
            auto d = o.cast<py::dict>();
            auto maxdemes = d["maxdemes"].cast<std::int32_t>();

            // Deme properties uses strong types
            ddemog::current_deme_sizes_vector current_deme_sizes(
                d["current_deme_sizes"]
                    .cast<ddemog::current_deme_sizes_vector::value_type>());
            ddemog::next_deme_sizes_vector next_deme_sizes(
                d["next_deme_sizes"].cast<ddemog::next_deme_sizes_vector::value_type>());
            ddemog::growth_rates_onset_times_vector growth_rates_onset_times(
                d["growth_rate_onset_times"]
                    .cast<ddemog::growth_rates_onset_times_vector::value_type>());
            ddemog::growth_initial_size_vector growth_initial_sizes(
                d["growth_initial_sizes"]
                    .cast<ddemog::growth_initial_size_vector::value_type>());
            ddemog::growth_rates_vector growth_rates(
                d["growth_rates"].cast<ddemog::growth_rates_vector::value_type>());
            ddemog::selfing_rates_vector selfing_rates(
                d["selfing_rates"].cast<ddemog::selfing_rates_vector::value_type>());

            ddemog::deme_properties sizes_rates(
                std::move(current_deme_sizes), std::move(next_deme_sizes),
                std::move(growth_rates_onset_times), std::move(growth_initial_sizes),
                std::move(growth_rates), std::move(selfing_rates));

            ddemog::MigrationMatrix M{};
            auto t = d["migmatrix"].cast<py::tuple>();
            auto rates = t[0].cast<std::vector<double>>();
            auto npops = t[1].cast<std::size_t>();
            auto scaled = t[2].cast<bool>();
            M = ddemog::MigrationMatrix(std::move(rates), npops, scaled);

            model_state.maxdemes = maxdemes;
            model_state.current_deme_parameters.current_deme_sizes.get()
                = std::move(current_deme_sizes.get());
            model_state.current_deme_parameters.next_deme_sizes.get()
                = std::move(next_deme_sizes.get());
            model_state.current_deme_parameters.growth_rate_onset_times.get()
                = std::move(growth_rates_onset_times.get());
            model_state.current_deme_parameters.growth_initial_sizes.get()
                = std::move(growth_initial_sizes.get());
            model_state.current_deme_parameters.growth_rates.get()
                = std::move(growth_rates.get());
            model_state.current_deme_parameters.selfing_rates.get()
                = std::move(selfing_rates.get());
            model_state.M = std::move(M);

            model_state.fitness_bookmark.starts
                = d["bookmark"]["starts"].cast<std::vector<std::uint32_t>>();
            model_state.fitness_bookmark.stops
                = d["bookmark"]["stops"].cast<std::vector<std::uint32_t>>();
            model_state.fitness_bookmark.offsets
                = d["bookmark"]["offsets"].cast<std::vector<std::uint32_t>>();
            model_state.fitness_bookmark.individuals
                = d["bookmark"]["individuals"].cast<std::vector<std::uint32_t>>();
            model_state.fitness_bookmark.individual_fitness
                = d["bookmark"]["individual_fitness"].cast<std::vector<double>>();

            self.set_model_state(std::move(model_state));
        });
}
