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

#ifndef FWDPY11_GENETIC_VALUES_DETAILS_PICKLE_ADDITIVE_HPP
#define FWDPY11_GENETIC_VALUES_DETAILS_PICKLE_ADDITIVE_HPP

#include <pybind11/pybind11.h>

struct pickle_gvalue
{
    template <typename T>
    inline pybind11::object
    operator()(const T& self) const
    {
        return pybind11::make_tuple(self.total_dim, self.scaling(), self.is_fitness());
    }
};

template <typename gvalue_type, typename fitness_type, typename trait_type>
gvalue_type
generate_returned_types(pybind11::tuple t, std::true_type)
{
    if (t.size() != 3)
        {
            throw std::runtime_error("invalid object state");
        }
    auto t0 = t[0].cast<pybind11::tuple>();
    if (t0.size() != 3)
        {
            throw std::runtime_error("invalid tuple size when unpickling Additive");
        }
    auto ndemes = t0[0].cast<std::size_t>();
    auto scaling = t0[1].cast<double>();
    auto is_fitness = t0[2].cast<bool>();
    if (is_fitness == true)
        {
            return gvalue_type(ndemes, scaling, fitness_type());
        }

    auto p = pybind11::module::import("pickle");
    auto t1 = p.attr("loads")(t[1]);
    auto t2 = p.attr("loads")(t[2]);
    //Do the casts in the constructor
    //to avoid any nasty issues w/
    //refs to temp
    return gvalue_type(ndemes, scaling, trait_type(),
                       t1.cast<const fwdpy11::GeneticValueToFitnessMap&>(),
                       t2.cast<const fwdpy11::GeneticValueNoise&>());
}

template <typename gvalue_type, typename fitness_type, typename trait_type>
gvalue_type
generate_returned_types(pybind11::tuple t, std::false_type)
{
    throw std::runtime_error("not implemented");
}

template <typename gvalue_type, typename fitness_type, typename trait_type,
          typename use_scaling>
gvalue_type
unpickle_gvalue(pybind11::tuple t)
{
    return generate_returned_types<gvalue_type, fitness_type, trait_type>(t,
                                                                          use_scaling());
}

#endif
