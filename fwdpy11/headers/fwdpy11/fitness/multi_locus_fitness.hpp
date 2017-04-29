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

#ifndef FWDPY11_MULTI_LOCUS_FITNESS_HPP__
#define FWDPY11_MULTI_LOCUS_FITNESS_HPP__

#include <algorithm>
#include <memory>
#include <vector>
#include <pybind11/numpy.h>
#include <fwdpy11/fitness/single_locus_fitness.hpp>

namespace fwdpy11
{
    struct multilocus_genetic_value
    {
        /// We store a container of single-locus fitness function wrapper
        /// objects
        std::vector<std::unique_ptr<fwdpy11::single_locus_fitness>>
            fitness_functions;
        /// These are the bound callbacks to the elements in fitness_functions
        std::vector<fwdpy11::single_locus_fitness_fxn> callbacks;
        /// A buffer used to init the NumPy array:
        std::unique_ptr<double> genetic_value_buffer;
        /// The genetic values are accessible to Python directly as a NumPy
        /// array:
        mutable pybind11::array_t<double> genetic_values_np;

        // Constructor takes vector of single region functions
        // and a function mapping individual-locus values -> overall value.
        // It is a bit ugly, but I use x_ for the argument corresponding
        // to class member x.
        multilocus_genetic_value(
            const std::vector<std::shared_ptr<fwdpy11::single_locus_fitness>>&
                fitness_functions_)
            : fitness_functions{}, callbacks{},
              genetic_value_buffer(new double[fitness_functions_.size()]),
              // The numpy array is intialized to share the
              // buffer of the std::vector:
              genetic_values_np{ pybind11::buffer_info(
                  genetic_value_buffer.get(), sizeof(double),
                  pybind11::format_descriptor<double>::format(), 1,
                  { fitness_functions_.size() }, { sizeof(double) }) }
        {
            if (fitness_functions_.empty())
                {
                    throw std::invalid_argument(
                        "empty list of fitness functions not allowed");
                }
            // Clone the input data into unique_ptr.
            // This is for run-time safety.
            for (auto&& ffi : fitness_functions_)
                {
                    // Single-locus ff are required to
                    // provide a function returning a copy
                    // wrapped in std::unique_ptr, which is
                    // what we call here:
                    fitness_functions.emplace_back(ffi->clone_unique());
                }
            // Generate the required callbacks:
            for (auto&& ffi : fitness_functions)
                {
                    callbacks.emplace_back(ffi->callback());
                }
        }

        inline pybind11::array_t<double>
        operator()(const fwdpy11::multilocus_diploid_t& dip,
                   const fwdpy11::gcont_t& gametes,
                   const fwdpy11::mcont_t& mutations) const
        {
            // use std::transform instead of loop b/c
            // genetic_values_np.mutable_at does
            // range-checking, which we do not need
            std::transform(dip.cbegin(), dip.cend(), callbacks.cbegin(),
                           genetic_values_np.mutable_data(),
                           [&gametes, &mutations](
                               const fwdpy11::diploid_t& dip_locus_i,
                               const fwdpy11::single_locus_fitness_fxn& wf) {
                               return wf(dip_locus_i, gametes, mutations);
                           });
            return genetic_values_np;
        }
        inline std::size_t
        size() const
        {
            return callbacks.size();
        }
        inline void
        update(const multilocus_t& pop)
        {
            for (auto&& fi : fitness_functions)
                {
                    fi->update(pop);
                }
        }
    };
}
#endif
