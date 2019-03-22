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
#include <pybind11/functional.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/fitness/fitness.hpp>

namespace py = pybind11;

struct genetic_value : public fwdpy11::single_locus_fitness
{
    fwdpy11::single_locus_fitness_fxn ff;
    std::function<void(const fwdpy11::DiploidPopulation &)> slocuspop_updater;
    std::function<void(const fwdpy11::MlocusPop &)> mlocuspop_updater;

    genetic_value(fwdpy11::single_locus_fitness_fxn ff_)
        : ff{ std::move(ff_) }, slocuspop_updater{}, mlocuspop_updater{}
    {
        slocuspop_updater = [this](const fwdpy11::DiploidPopulation &pop) {};
        mlocuspop_updater = [this](const fwdpy11::MlocusPop &pop) {};
    }

    genetic_value(
        fwdpy11::single_locus_fitness_fxn ff_,
        std::function<void(const fwdpy11::DiploidPopulation &)> slocuspop_updater_)
        : ff{ std::move(ff_) },
          slocuspop_updater{ std::move(slocuspop_updater_) },
          mlocuspop_updater{}
    {
        mlocuspop_updater = [this](const fwdpy11::MlocusPop &pop) {};
    }

    genetic_value(
        fwdpy11::single_locus_fitness_fxn ff_,
        std::function<void(const fwdpy11::DiploidPopulation &)> slocuspop_updater_,
        std::function<void(const fwdpy11::MlocusPop &)> mlocuspop_updater_)
        : ff{ std::move(ff_) },
          slocuspop_updater{ std::move(slocuspop_updater_) },
          mlocuspop_updater{ std::move(mlocuspop_updater_) }
    {
    }

    void
    update(const fwdpy11::DiploidPopulation &pop)
    {
        slocuspop_updater(pop);
    }
    void
    update(const fwdpy11::MlocusPop &pop)
    {
        mlocuspop_updater(pop);
    }
    virtual fwdpy11::single_locus_fitness_fxn
    callback() const
    {
        return std::bind(ff, std::placeholders::_1, std::placeholders::_2,
                         std::placeholders::_3);
    }
    SINGLE_LOCUS_FITNESS_CLONE_SHARED(genetic_value);
    SINGLE_LOCUS_FITNESS_CLONE_UNIQUE(genetic_value);
    SINGLE_LOCUS_FITNESS_CALLBACK_NAME("Custom genetic value object");
};

PYBIND11_MODULE(python_genetic_values, m)
{
    m.doc() = "Enable definition of custom single-locus genetic "
              "value/fitness functions using pure Python functions.";

    FWDPY11_SINGLE_LOCUS_FITNESS()

    py::class_<genetic_value, std::shared_ptr<genetic_value>,
               fwdpy11::single_locus_fitness>(
        m, "GeneticValue", "Allow for custom genetic value calculations in "
                           "Python. See :ref:`customgvalues` for details.")
        .def(py::init<fwdpy11::single_locus_fitness_fxn>())
        .def(py::init<fwdpy11::single_locus_fitness_fxn,
                      std::function<void(const fwdpy11::DiploidPopulation &)>>())
        .def(py::init<fwdpy11::single_locus_fitness_fxn,
                      std::function<void(const fwdpy11::DiploidPopulation &)>,
                      std::function<void(const fwdpy11::MlocusPop &)>>())
        .def("__call__", [](const std::shared_ptr<genetic_value> &aw,
                            const fwdpy11::Diploid &dip,
                            const fwdpy11::DiploidPopulation &pop) {
            return aw->callback()(dip, pop.gametes, pop.mutations);
        });
}
