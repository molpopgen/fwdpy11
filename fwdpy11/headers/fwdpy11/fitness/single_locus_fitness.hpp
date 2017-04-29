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

#ifndef FWDPY11_SINGLE_LOCUS_FITNESS_HPP__
#define FWDPY11_SINGLE_LOCUS_FITNESS_HPP__

#include <typeinfo>
#include <memory>
#include <fwdpy11/types.hpp>
#include <fwdpp/fitness_models.hpp>

namespace fwdpy11
{
    template <typename... T>
    using single_locus_fitness_signature
        = std::function<double(const fwdpy11::diploid_t &,
                               const fwdpy11::gcont_t &,
                               const fwdpy11::mcont_t &, T...)>;

    /*! Single-deme fitness function signature for standard "popgen"
     *  simulations.
     */
    using single_locus_fitness_fxn = single_locus_fitness_signature<>;

    struct single_locus_fitness
    //! Pure virtual base class for single-deme fitness functions
    {
        virtual ~single_locus_fitness() = default;
        single_locus_fitness() = default;
        virtual void
        update(const singlepop_t &pop)
        {
        }
		virtual void update(const multilocus_t & pop)
		{
		}
        virtual single_locus_fitness_fxn callback() const = 0;
        virtual std::unique_ptr<single_locus_fitness> clone_unique() const = 0;
        virtual std::shared_ptr<single_locus_fitness> clone_shared() const = 0;
        virtual std::string callback_name() const = 0;
    };

#define SINGLE_LOCUS_FITNESS_CLONE_UNIQUE(TYPENAME)                           \
    std::unique_ptr<fwdpy11::single_locus_fitness> clone_unique() const       \
    {                                                                         \
        using ptr = std::unique_ptr<fwdpy11::single_locus_fitness>;           \
        return ptr(new TYPENAME(*this));                                      \
    }

#define SINGLE_LOCUS_FITNESS_CLONE_SHARED(TYPENAME)                           \
    std::shared_ptr<fwdpy11::single_locus_fitness> clone_shared() const       \
    {                                                                         \
        return std::shared_ptr<fwdpy11::single_locus_fitness>(                \
            new TYPENAME(*this));                                             \
    }

#define SINGLE_LOCUS_FITNESS_CALLBACK_NAME(T)                                 \
    std::string callback_name() const { return T; }

    template <typename fitness_model_type>
    struct fwdpp_single_locus_fitness_wrapper : public single_locus_fitness
    {
        using fitness_model = fitness_model_type;
        const double scaling;
        fwdpp_single_locus_fitness_wrapper(const double scaling_ = 2.0)
            : scaling(scaling_)
        {
            if (!std::isfinite(scaling))
                {
                    throw std::runtime_error("non-finite value. "
                                             + std::string(__FILE__) + " line "
                                             + std::to_string(__LINE__));
                }
        }
        inline single_locus_fitness_fxn
        callback() const final
        {
            return std::bind(fitness_model(), std::placeholders::_1,
                             std::placeholders::_2, std::placeholders::_3,
                             scaling);
        }

        SINGLE_LOCUS_FITNESS_CLONE_SHARED(
            fwdpp_single_locus_fitness_wrapper<fitness_model>);
        SINGLE_LOCUS_FITNESS_CLONE_UNIQUE(
            fwdpp_single_locus_fitness_wrapper<fitness_model>);
        SINGLE_LOCUS_FITNESS_CALLBACK_NAME(typeid(fitness_model()).name()
                                           + std::string(" with scaling = ")
                                           + std::to_string(this->scaling));
    };

    using single_locus_mult_wrapper
        = fwdpp_single_locus_fitness_wrapper<KTfwd::multiplicative_diploid>;
    using single_locus_additive_wrapper
        = fwdpp_single_locus_fitness_wrapper<KTfwd::additive_diploid>;

#define FWDPY11_SINGLE_LOCUS_FITNESS()                                        \
    pybind11::object FWDPY11_SINGLE_LOCUS_FITNESS_BASE_IMPORT__               \
        = (pybind11::object)pybind11::module::import("fwdpy11.fitness")       \
              .attr("SlocusFitness");
}

#endif
