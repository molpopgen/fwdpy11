#ifndef FWDPY11_MULTIVARIATEGAUSSIAN_HPP
#define FWDPY11_MULTIVARIATEGAUSSIAN_HPP

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <functional>
#include <memory>
#include <fwdpy11/policies/mutation.hpp>
#include <fwdpy11/gsl/gsl_error_handler_wrapper.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{
    struct MultivariateGaussianEffects : public Sregion
    {
        using matrix_ptr
            = std::unique_ptr<gsl_matrix, std::function<void(gsl_matrix *)>>;
        using vector_ptr
            = std::unique_ptr<gsl_vector, std::function<void(gsl_vector *)>>;
        matrix_ptr input_matrix_copy;
        // Stores the Cholesky decomposition
        matrix_ptr matrix;
        // Stores the means of the mvn distribution, which are all zero
        vector_ptr mu;
        double fixed_effect;
        mutable std::vector<double> effect_sizes, dominance_values;
        // Stores the results of gsl_ran_multivariate_gaussian
        mutable gsl_vector_view res;

        template <typename Dominance>
        MultivariateGaussianEffects(const Region &r, const double sc,
                                    const gsl_matrix &input_matrix, double s,
                                    Dominance &&h)
            : Sregion(r, sc, input_matrix.size1, std::forward<Dominance>(h)),
              input_matrix_copy(gsl_matrix_alloc(input_matrix.size1, input_matrix.size2),
                                [](gsl_matrix *m) { gsl_matrix_free(m); }),
              matrix(gsl_matrix_alloc(input_matrix.size1, input_matrix.size2),
                     [](gsl_matrix *m) { gsl_matrix_free(m); }),
              // NOTE: use of calloc to initialize mu to all zeros
              mu(gsl_vector_calloc(input_matrix.size1),
                 [](gsl_vector *v) { gsl_vector_free(v); }),
              fixed_effect(s), effect_sizes(input_matrix.size1),
              dominance_values(input_matrix.size1, std::numeric_limits<double>::quiet_NaN()),
              // Holds the results of calls to mvn, and maps the
              // output to effect_sizes
              res(gsl_vector_view_array(effect_sizes.data(), effect_sizes.size()))
        {
            if (!std::isfinite(fixed_effect))
                {
                    throw std::invalid_argument("fixed_effect must be finite");
                }
            if (matrix->size1 != matrix->size2)
                {
                    throw std::invalid_argument("input matrix must be square");
                }
            // If any values are non-finite, throw an exception
            if (std::any_of(input_matrix.data,
                            input_matrix.data + input_matrix.size1 * input_matrix.size2,
                            [](double d) { return !std::isfinite(d); })
                == true)
                {
                    throw std::invalid_argument(
                        "input matrix contains non-finite values");
                }

            gsl_scoped_disable_error_handler_wrapper gsl_error_scope_guard;
            // Assign the matrix and do the Cholesky decomposition

            int rv = gsl_matrix_memcpy(matrix.get(), &input_matrix);
            if (rv != GSL_SUCCESS)
                {
                    throw std::runtime_error("failure copying input matrix");
                }
            rv = gsl_matrix_memcpy(input_matrix_copy.get(), &input_matrix);
            if (rv != GSL_SUCCESS)
                {
                    throw std::runtime_error("failure copying input matrix");
                }
            rv = gsl_linalg_cholesky_decomp1(matrix.get());
            if (rv == GSL_EDOM)
                {
                    throw std::invalid_argument("Cholesky decomposition failed");
                }
        }

        virtual std::unique_ptr<Sregion>
        clone() const override
        {
            return std::unique_ptr<MultivariateGaussianEffects>(
                new MultivariateGaussianEffects(
                    fwdpy11::Region(this->beg(), this->end(), this->weight(),
                                    this->region.coupled, this->label()),
                    1.0, *input_matrix_copy.get(), this->fixed_effect, this->dominance));
        }

        virtual std::uint32_t
        operator()(fwdpp::flagged_mutation_queue &recycling_bin,
                   std::vector<Mutation> &mutations,
                   std::unordered_multimap<double, std::uint32_t> &lookup_table,
                   const std::uint32_t generation, const GSLrng_t &rng) const override
        {
            int rv = gsl_ran_multivariate_gaussian(rng.get(), mu.get(), matrix.get(),
                                                   &res.vector);
            for (std::size_t i = 0; i < effect_sizes.size(); ++i)
                {
                    dominance_values[i]
                        = dominance->generate_dominance(rng, effect_sizes[i]);
                }
            if (rv != GSL_SUCCESS)
                {
                    throw std::runtime_error(
                        "call to gsl_ran_multivariate_gaussian failed");
                }
            return infsites_Mutation(
                recycling_bin, mutations, lookup_table, false, generation,
                [this, &rng]() { return region(rng); },
                [this]() { return fixed_effect; },
                [this, &rng](const double esize) {
                    return dominance->generate_dominance(rng, esize);
                },
                [this]() { return effect_sizes; }, [this]() { return dominance_values; },
                this->label());
        }

        double
        from_mvnorm(const double deviate, const double /*P*/) const override
        {
            return deviate / scaling;
        }
    };
} // namespace fwdpy11

#endif
