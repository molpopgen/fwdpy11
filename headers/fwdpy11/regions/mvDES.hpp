//
// Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
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

#ifndef FWDPY11_MVDES_HPP
#define FWDPY11_MVDES_HPP

#include <functional>
#include <algorithm>
#include <cmath>
#include "Sregion.hpp"
#include "LogNormalS.hpp"
#include "MultivariateGaussianEffects.hpp"
#include <fwdpy11/policies/mutation.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <limits>
#include <locale>
#include <stdexcept>

namespace fwdpy11
{
    class mvDES : public Sregion
    // NOTE: the class details here are rather complex.
    // In the future, it may be better to define a second
    // class that mvDES encapsulates.  Doing so would solve
    // much of the constructor, pickling, etc., problems,
    // handled here.
    {
      private:
        using matrix_ptr
            = std::unique_ptr<gsl_matrix, std::function<void(gsl_matrix *)>>;
        using vector_ptr
            = std::unique_ptr<gsl_vector, std::function<void(gsl_vector *)>>;
        using callback_type = std::function<std::uint32_t(
            const mvDES *, fwdpp::flagged_mutation_queue &r, std::vector<Mutation> &,
            std::unordered_multimap<double, std::uint32_t> &, const std::uint32_t,
            const GSLrng_t &)>;

        struct default_callback
        {
            inline std::uint32_t
            operator()(const mvDES *outer_this,
                       fwdpp::flagged_mutation_queue &recycling_bin,
                       std::vector<Mutation> &mutations,
                       std::unordered_multimap<double, std::uint32_t> &lookup_table,
                       const std::uint32_t generation, const GSLrng_t &rng) const
            {
                outer_this->generate_deviates(rng);
                for (std::size_t i = 0; i < outer_this->deviates.size(); ++i)
                    {
                        // Subtract means from the deviates so we
                        // can use N(0, sigma[i]) cdf.
                        outer_this->deviates[i]
                            = outer_this->output_distributions[i]->from_mvnorm(
                                outer_this->deviates[i],
                                gsl_cdf_gaussian_P(outer_this->deviates[i]
                                                       - outer_this->means[i],
                                                   outer_this->stddev[i]));
                        outer_this->dominance_values[i]
                            = outer_this->output_distributions[i]->generate_dominance(
                                rng, outer_this->deviates[i]);
                    }
                return outer_this->generate_mutation(recycling_bin, mutations,
                                                     lookup_table, generation, rng);
            }
        };

        struct specialized_callback
        {
            inline std::uint32_t
            operator()(const mvDES *outer_this,
                       fwdpp::flagged_mutation_queue &recycling_bin,
                       std::vector<Mutation> &mutations,
                       std::unordered_multimap<double, std::uint32_t> &lookup_table,
                       const std::uint32_t generation, const GSLrng_t &rng) const
            {
                outer_this->generate_deviates(rng);
                for (std::size_t i = 0; i < outer_this->deviates.size(); ++i)
                    {
                        outer_this->deviates[i]
                            = outer_this->output_distributions[0]->from_mvnorm(
                                outer_this->deviates[i],
                                gsl_cdf_gaussian_P(outer_this->deviates[i],
                                                   outer_this->stddev[i]));
                        outer_this->dominance_values[i]
                            = outer_this->output_distributions[0]->generate_dominance(
                                rng, outer_this->deviates[i]);
                    }
                return outer_this->generate_mutation(recycling_bin, mutations,
                                                     lookup_table, generation, rng);
            }
        };

        void
        generate_deviates(const GSLrng_t &rng) const
        {
            int rv = gsl_ran_multivariate_gaussian(rng.get(), &mu.vector, matrix.get(),
                                                   &res.vector);
            if (rv != GSL_SUCCESS)
                {
                    throw std::runtime_error(
                        "call to gsl_ran_multivariate_gaussian failed");
                }
        }

        std::size_t
        generate_mutation(fwdpp::flagged_mutation_queue &recycling_bin,
                          Population::mutation_container &mutations,
                          Population::lookup_table_t &lookup_table,
                          const fwdpp::uint_t &generation, const GSLrng_t &rng) const
        {
            return infsites_Mutation(
                recycling_bin, mutations, lookup_table, false, generation,
                [this, &rng]() { return this->region(rng); }, []() { return 0.0; },
                [](const double) { return 1.0; }, [this]() { return this->deviates; },
                [this]() { return this->dominance_values; }, this->label());
        }

        Region
        get_region(const std::vector<std::unique_ptr<Sregion>> &odists)
        {
            if (odists.empty())
                {
                    throw std::invalid_argument("empty list of Sregions");
                }
            double beg = odists[0]->beg();
            double end = odists[0]->end();
            double weight = odists[0]->weight();
            double label = odists[0]->label();
            double coupled = odists[0]->region.coupled;
            for (std::size_t i = 0; i < odists.size(); ++i)
                {
                    auto beg_match = beg == odists[i]->beg();
                    auto end_match = end == odists[i]->end();
                    auto weight_match = weight == odists[i]->weight();
                    auto label_match = label == odists[i]->label();
                    auto coupled_match = coupled == odists[i]->region.coupled;
                    if (!beg_match || !end_match || !weight_match || !label_match
                        || !coupled_match)
                        {
                            throw std::invalid_argument("all Region fields must match");
                        }
                }
            return Region(beg, end, weight, coupled, label);
        }

        std::vector<std::unique_ptr<Sregion>>
        fill_output_distributions(const std::vector<std::unique_ptr<Sregion>> &odist,
                                  std::size_t n)
        {
            if (odist.size() > 1 && odist.size() != n)
                {
                    throw std::invalid_argument("invalid number of Sregion objects");
                }
            std::vector<std::unique_ptr<Sregion>> rv;
            for (auto &&o : odist)
                {
                    rv.emplace_back(o->clone());
                }
            if (odist.size() == 1)
                {
                    for (std::size_t i = 1; i < n; ++i)
                        {
                            rv.emplace_back(odist[0]->clone());
                        }
                }
            return rv;
        }

        std::vector<std::unique_ptr<Sregion>>
        clone_into_vector(const Sregion &odist)
        {
            std::vector<std::unique_ptr<Sregion>> temp;
            temp.emplace_back(odist.clone());
            return temp;
        }

        std::vector<double>
        fill_dominance(const std::vector<std::unique_ptr<Sregion>> &odist)
        // NOTE: domimnance needs to be pushed up to Sregion?
        {
            std::vector<double> rv(odist.size(),
                                   std::numeric_limits<double>::quiet_NaN());
            return rv;
        }

        matrix_ptr
        copy_input_matrix(const gsl_matrix &vcov)
        {
            matrix_ptr m(gsl_matrix_alloc(vcov.size1, vcov.size2),
                         [](gsl_matrix *m) { gsl_matrix_free(m); });
            // Assign the matrix and do the Cholesky decomposition
            auto error_handler = gsl_set_error_handler_off();

            int rv = gsl_matrix_memcpy(m.get(), &vcov);
            if (rv != GSL_SUCCESS)
                {
                    // Reset error handler on the way out
                    gsl_set_error_handler(error_handler);
                    throw std::runtime_error("failure copying input matrix");
                }
            // If any values are non-finite, throw an exception
            if (std::any_of(m->data, m->data + m->size1 * m->size2,
                            [](double d) { return !std::isfinite(d); })
                == true)
                {
                    gsl_set_error_handler(error_handler);
                    throw std::invalid_argument(
                        "input matrix contains non-finite values");
                }
            // Reset error handler on the way out
            gsl_set_error_handler(error_handler);
            return m;
        }

        matrix_ptr
        decompose()
        {
            matrix_ptr m(gsl_matrix_alloc(vcov_copy->size1, vcov_copy->size2),
                         [](gsl_matrix *m) { gsl_matrix_free(m); });
            // Assign the matrix and do the Cholesky decomposition
            auto error_handler = gsl_set_error_handler_off();
            int rv = gsl_matrix_memcpy(m.get(), vcov_copy.get());
            if (rv != GSL_SUCCESS)
                {
                    // Reset error handler on the way out
                    gsl_set_error_handler(error_handler);
                    throw std::runtime_error("failure copying input matrix");
                }

            rv = gsl_linalg_cholesky_decomp1(m.get());
            if (rv == GSL_EDOM)
                {
                    // Reset error handler on the way out
                    gsl_set_error_handler(error_handler);
                    throw std::invalid_argument("Cholesky decomposition failed");
                }
            gsl_set_error_handler(error_handler);
            return m;
        }

        void
        finalize_setup()
        {
            if (matrix->size1 != matrix->size2)
                {
                    throw std::invalid_argument("input matrix must be square");
                }
            if (means.size() != matrix->size1)
                {
                    throw std::invalid_argument(
                        "length of means does not match matrix dimensions");
                }
        }

        std::vector<double>
        get_standard_deviations()
        {
            std::vector<double> rv;
            gsl_vector_const_view d = gsl_matrix_const_diagonal(vcov_copy.get());
            for (std::size_t i = 0; i < d.vector.size; ++i)
                {
                    rv.push_back(std::sqrt(gsl_vector_get(&d.vector, i)));
                }
            return rv;
        }

        std::vector<std::unique_ptr<Sregion>> output_distributions;
        matrix_ptr vcov_copy, matrix;
        mutable std::vector<double> deviates, dominance_values, means;
        // Stores the results of gsl_ran_multivariate_gaussian
        mutable gsl_vector_view res;
        // Stores the means of the mvn distribution, which are all zero
        gsl_vector_const_view mu;
        std::vector<double> stddev;
        const bool lognormal_init, mvgaussian_init;
        const callback_type callback;

      public:
        // NOTE: the "scaling" concept is handled by the output_distributions.
        mvDES(const std::vector<std::unique_ptr<Sregion>> &odist,
              std::vector<double> gaussian_means, const gsl_matrix &vcov)
            : Sregion(get_region(odist), 1., odist.size(),
                      process_input_dominance(0.)), // 1. is a dummy param here.
              output_distributions(fill_output_distributions(odist, vcov.size1)),
              vcov_copy(copy_input_matrix(vcov)), matrix(decompose()),
              deviates(vcov.size1), dominance_values(fill_dominance(odist)),
              means(std::move(gaussian_means)),
              res(gsl_vector_view_array(deviates.data(), deviates.size())),
              // NOTE: use of calloc to initialize mu to all zeros
              mu(gsl_vector_const_view_array(means.data(), means.size())),
              stddev(get_standard_deviations()), lognormal_init(false),
              mvgaussian_init(false), callback(default_callback())
        {
            finalize_setup();
        }

        // FIXME: to get the next 2 to work, the callbacks used
        // need some thought...

        mvDES(const LogNormalS &odist, std::vector<double> gaussian_means,
              const gsl_matrix &vcov)
            : Sregion(odist.region, 1., vcov.size1, process_input_dominance(0.)),
              output_distributions(clone_into_vector(odist)),
              vcov_copy(copy_input_matrix(vcov)), matrix(decompose()),
              deviates(vcov.size1),
              // FIXME: this is wrong
              dominance_values(gaussian_means.size(),
                               std::numeric_limits<double>::quiet_NaN()),
              means(std::move(gaussian_means)),
              res(gsl_vector_view_array(deviates.data(), deviates.size())),
              // NOTE: use of calloc to initialize mu to all zeros
              mu(gsl_vector_const_view_array(means.data(), means.size())),
              stddev(get_standard_deviations()), lognormal_init(true),
              mvgaussian_init(false), callback(specialized_callback())
        {
            finalize_setup();
        }

        mvDES(const MultivariateGaussianEffects &odist,
              std::vector<double> gaussian_means)
            : Sregion(odist.region, 1., odist.input_matrix_copy->size1,
                      process_input_dominance(0.)),
              output_distributions(clone_into_vector(odist)),
              vcov_copy(copy_input_matrix(*(odist.input_matrix_copy))),
              matrix(decompose()), deviates(odist.input_matrix_copy->size1),
              dominance_values(odist.dominance_values), means(std::move(gaussian_means)),
              res(gsl_vector_view_array(deviates.data(), deviates.size())),
              // NOTE: use of calloc to initialize mu to all zeros
              mu(gsl_vector_const_view_array(means.data(), means.size())),
              stddev(get_standard_deviations()), lognormal_init(false),
              mvgaussian_init(true), callback(specialized_callback())
        {
            finalize_setup();
        }

        std::uint32_t
        operator()(fwdpp::flagged_mutation_queue &recycling_bin,
                   std::vector<Mutation> &mutations,
                   std::unordered_multimap<double, std::uint32_t> &lookup_table,
                   const std::uint32_t generation, const GSLrng_t &rng) const override
        {
            return callback(this, recycling_bin, mutations, lookup_table, generation,
                            rng);
        }

        std::unique_ptr<Sregion>
        clone() const override
        {
            if (lognormal_init)
                {
                    return std::unique_ptr<mvDES>(new mvDES(
                        *dynamic_cast<LogNormalS *>(output_distributions[0].get()),
                        means, *(this->vcov_copy)));
                }
            else if (mvgaussian_init)
                {
                    return std::unique_ptr<mvDES>(
                        new mvDES(*dynamic_cast<MultivariateGaussianEffects *>(
                                      output_distributions[0].get()),
                                  means));
                }
            return std::unique_ptr<mvDES>(
                new mvDES(output_distributions, means, *(this->vcov_copy)));
        }

        double
        from_mvnorm(const double, const double) const override
        {
            throw std::invalid_argument(
                "mvDES is not allowed to be part of multivariate DES");
        }

        double
        generate_dominance(const GSLrng_t & /*rng*/,
                           const double /*esize*/) const override final
        {
            throw std::runtime_error("mvDES::generate_dominance is not implemented");
            return std::numeric_limits<double>::quiet_NaN();
        }
    };
} // namespace fwdpy11

#endif
