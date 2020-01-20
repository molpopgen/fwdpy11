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
    ddemog::DiscreteDemography
    init_from_numpy(py::array_t<std::uint32_t> popsizes)
    {
        auto r = popsizes.unchecked<1>();
        auto len = r.shape(0);
        if (len == 0)
            {
                throw std::invalid_argument("empty length of deme sizes");
            }

        std::uint32_t N = r(0);
        ddemog::DiscreteDemography::set_deme_sizes_vector p(
            1, ddemog::SetDemeSize(0, 0, r(0), true));
        for (decltype(len) i = 1; i < len; ++i)
            {
                if (r(i) != N)
                    {
                        p.emplace_back(i, 0, r(i), true);
                        N = r(i);
                    }
            }
        return ddemog::DiscreteDemography({}, {}, p, {}, nullptr, {});
    }

    ddemog::MigrationMatrix
    decode_migration_matrix_input(py::object o)
    {
        const auto a2v = [](py::array_t<double> m) {
            auto r = m.unchecked<2>();
            if (r.shape(0) != r.shape(1))
                {
                    throw std::invalid_argument(
                        "MigrationMatrix must be square");
                }
            auto d = r.data(0, 0);
            std::vector<double> M(d, d + r.shape(0) * r.shape(1));
            return M;
        };

        try
            {
                py::array_t<double> m = o.cast<py::array_t<double>>();
                auto M = a2v(m);
                return ddemog::MigrationMatrix(std::move(M), m.shape(0),
                                                true);
            }
        catch (...)
            {
                try
                    {
                        py::tuple t = o.cast<py::tuple>();
                        if (t.size() != 2)
                            {
                                throw std::invalid_argument(
                                    "MigrationMatrix: tuple contain two "
                                    "elements");
                            }
                        py::array_t<double> m
                            = t[0].cast<py::array_t<double>>();
                        auto M = a2v(m);
                        bool scaled = t[1].cast<bool>();
                        return ddemog::MigrationMatrix(std::move(M),
                                                        m.shape(0), scaled);
                    }
                catch (...)
                    {
                    }
            }
        ddemog::MigrationMatrix M = o.cast<ddemog::MigrationMatrix>();
        return M;
    }
} // namespace

static const auto INIT_DOCSTRING_NUMPY = R"delim(
:param popsizes: A list of deme sizes over time.
:type popsizes: numpy.ndarray with dtype numpy.uint32
)delim";

static const auto INIT_DOCSTRING_LISTS = R"delim(
:param mass_migrations: Instances of :class:`fwdpy11.MassMigration`
:type mass_migrations: list
:param set_growth_rates: Instances of :class:`fwdpy11.SetExponentialGrowth`
:type set_growth_rates: list
:param set_deme_sizes: Instances of :class:`fwdpy11.SetDemeSize`
:type set_deme_sizes: list
:param set_selfing_rates: Instances of :class:`fwdpy11.SetSelfingRate`
:type set_selfing_rates: list
:param migmatrix: A migraton matrix. See :ref:`migration`.
:param set_migration_rates: Instances of :class:`fwdpy11.SetMigrationRates`
:type set_migration_rates: list
)delim";

void
init_DiscreteDemography(py::module& m)
{
    // NOTE: order of __init__ declarations matter
    // so that init from numpy works as expected.
    // TODO: migration matrix & mig matrix changes
    py::class_<ddemog::DiscreteDemography>(m, "DiscreteDemography")
        .def(py::init([](py::array_t<std::uint32_t> popsizes) {
                 return init_from_numpy(popsizes);
             }),
             INIT_DOCSTRING_NUMPY)
        .def(
            py::init([](py::object mass_migration_events,
                        py::object set_growth_rates, py::object set_deme_sizes,
                        py::object set_selfing_rates, py::object migmatrix,
                        py::object set_migration_rates) {
                ddemog::DiscreteDemography::mass_migration_vector morc;
                ddemog::DiscreteDemography::set_growth_rates_vector growth;
                ddemog::DiscreteDemography::set_deme_sizes_vector
                    change_sizes;
                ddemog::DiscreteDemography::set_selfing_rates_vector selfing;
                ddemog::DiscreteDemography::set_migration_rates_vector
                    migrates;

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
                        change_sizes
                            = set_deme_sizes.cast<decltype(change_sizes)>();
                    }
                if (set_selfing_rates.is_none() == false)
                    {
                        selfing = set_selfing_rates.cast<decltype(selfing)>();
                    }
                if (set_migration_rates.is_none() == false)
                    {
                        migrates
                            = set_migration_rates.cast<decltype(migrates)>();
                    }
                if (migmatrix.is_none() == true)
                    {
                        return ddemog::DiscreteDemography(
                            std::move(morc), std::move(growth),
                            std::move(change_sizes), std::move(selfing),
                            nullptr, std::move(migrates));
                    }
                std::unique_ptr<ddemog::MigrationMatrix> M(
                    new ddemog::MigrationMatrix(
                        decode_migration_matrix_input(migmatrix)));
                return ddemog::DiscreteDemography(
                    std::move(morc), std::move(growth),
                    std::move(change_sizes), std::move(selfing), std::move(M),
                    std::move(migrates));
            }),
            py::arg("mass_migrations") = py::none(),
            py::arg("set_growth_rates") = py::none(),
            py::arg("set_deme_sizes") = py::none(),
            py::arg("set_selfing_rates") = py::none(),
            py::arg("migmatrix") = py::none(),
            py::arg("set_migration_rates") = py::none(), INIT_DOCSTRING_LISTS)
        .def_readonly("mass_migrations",
                      &ddemog::DiscreteDemography::mass_migrations)
        .def_readonly("set_growth_rates",
                      &ddemog::DiscreteDemography::set_growth_rates)
        .def_readonly("set_deme_sizes",
                      &ddemog::DiscreteDemography::set_deme_sizes)
        .def_readonly("set_selfing_rates",
                      &ddemog::DiscreteDemography::set_selfing_rates)
        .def_property_readonly("migmatrix",
                               [](const ddemog::DiscreteDemography& self) -> const ddemog::MigrationMatrix*  {
                                   if(self.migmatrix == nullptr) { return nullptr; }
                                   return self.migmatrix.get();
                               })
        .def_readonly("set_migration_rates",
                      &ddemog::DiscreteDemography::set_migration_rates)
        .def(py::pickle(
            [](const ddemog::DiscreteDemography& self) {
                if (self.migmatrix == nullptr)
                    {
                        return py::make_tuple(
                            self.mass_migrations, self.set_growth_rates,
                            self.set_deme_sizes, self.set_selfing_rates,
                            py::none(), self.set_migration_rates);
                    }
                return py::make_tuple(
                    self.mass_migrations, self.set_growth_rates,
                    self.set_deme_sizes, self.set_selfing_rates,
                    self.migmatrix.get(), self.set_migration_rates);
            },
            [](py::tuple t) {
                py::object migobject = t[4].cast<py::object>();
                std::unique_ptr<ddemog::MigrationMatrix> M(nullptr);
                if (migobject.is_none() == false)
                    {
                        M.reset(new ddemog::MigrationMatrix(
                            t[4].cast<ddemog::MigrationMatrix>()));
                    }
                return ddemog::DiscreteDemography(
                    t[0].cast<
                        ddemog::DiscreteDemography::mass_migration_vector>(),
                    t[1].cast<ddemog::DiscreteDemography::
                                  set_growth_rates_vector>(),
                    t[2].cast<
                        ddemog::DiscreteDemography::set_deme_sizes_vector>(),
                    t[3].cast<ddemog::DiscreteDemography::
                                  set_selfing_rates_vector>(),
                    std::move(M),
                    t[5].cast<ddemog::DiscreteDemography::
                                  set_migration_rates_vector>());
            }));
}
