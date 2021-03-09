//
// Copyright (C) 2019 Kevin Thornton <krthornt@uci.edu>
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
#include <fwdpy11/genetic_value_noise/GeneticValueNoise.hpp>

struct IneritedNoise : public fwdpy11::GeneticValueNoise
{
    std::uint32_t generation;

    IneritedNoise() : GeneticValueNoise{}, generation{} {}

    double
    operator()(const fwdpy11::DiploidGeneticValueNoiseData data) const override
    {
        return data.parent1_metadata.get().e + generation;
    }

    void
    update(const fwdpy11::DiploidPopulation& pop) override
    {
        generation = pop.generation;
    }

    std::shared_ptr<fwdpy11::GeneticValueNoise>
    clone() const override
    {
        return std::make_shared<IneritedNoise>();
    }

    pybind11::object
    pickle() const
    {
        return pybind11::bytes("IneritedNoise");
    }

    static inline IneritedNoise
    unpickle(pybind11::object& o)
    {
        auto s = o.cast<std::string>();
        if (s != "IneritedNoise")
            {
                throw std::runtime_error("invalid object state");
            }
        return IneritedNoise();
    }
};

PYBIND11_MODULE(inherit_noise, m)
{
    pybind11::object imported_base
        = pybind11::module::import("fwdpy11").attr("GeneticValueNoise");

    pybind11::class_<IneritedNoise, fwdpy11::GeneticValueNoise>(m, "IneritedNoise")
        .def(pybind11::init<>())
        .def(pybind11::pickle(
            [](const IneritedNoise& self) { return self.pickle(); },
            [](pybind11::object o) { return IneritedNoise::unpickle(o); }));
}
