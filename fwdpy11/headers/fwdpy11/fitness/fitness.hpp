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
#ifndef FWDPY11_FITNESS_HPP__
#define FWDPY11_FITNESS_HPP__

#include <functional>
#include <cmath>
#include <stdexcept>
#include <string>
#include <fwdpp/fitness_models.hpp>

namespace fwdpy11
{
    template <typename... T>
    using singlepop_fitness_signature
        = std::function<double(const fwdpy11::diploid_t &,
                               const fwdpy11::gcont_t &,
                               const fwdpy11::mcont_t &, T...)>;

    /*! Single-deme fitness function signature for standard "popgen"
     *  simulations.
     */
    using singlepop_fitness_fxn = singlepop_fitness_signature<>;

    struct singlepop_fitness
    //! Pure virtual base class for single-deme fitness functions
    {
        virtual ~singlepop_fitness() = default;
        singlepop_fitness() = default;
        virtual singlepop_fitness_fxn callback() const = 0;
        virtual void
        update(const singlepop_t &pop)
        {
        }
    };

    template <typename fitness_model_type>
    struct fwdpp_singlepop_fitness_wrapper : public singlepop_fitness
    {
        using fitness_model = fitness_model_type;
        const double scaling;
        fwdpp_singlepop_fitness_wrapper(const double scaling_ = 2.0)
            : scaling(scaling_)
        {
            if (!std::isfinite(scaling))
                {
                    throw std::runtime_error("non-finite value. "
                                             + std::string(__FILE__) + " line "
                                             + std::to_string(__LINE__));
                }
        }
        inline singlepop_fitness_fxn
        callback() const final
        {
            return std::bind(fitness_model(), std::placeholders::_1,
                             std::placeholders::_2, std::placeholders::_3,
                             scaling);
        }
    };

    using singlepop_mult_wrapper
        = fwdpp_singlepop_fitness_wrapper<KTfwd::multiplicative_diploid>;
    using singlepop_additive_wrapper
        = fwdpp_singlepop_fitness_wrapper<KTfwd::additive_diploid>;

#define FWDPY11_SINGLEPOP_FITNESS()                                           \
    pybind11::object FWDPY11_SINGLEPOP_FITNESS_BASE_IMPORT__                  \
        = (pybind11::object)pybind11::module::import("fwdpy11.fitness")       \
              .attr("SpopFitness");
}
#endif
