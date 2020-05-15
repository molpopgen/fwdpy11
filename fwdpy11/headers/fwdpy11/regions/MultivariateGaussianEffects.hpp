#ifndef FWDPY11_MULTIVARIATEGAUSSIAN_HPP
#define FWDPY11_MULTIVARIATEGAUSSIAN_HPP

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <memory>
#include <fwdpy11/policies/mutation.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{
    struct MultivariateGaussianEffects : public Sregion
    {
        using matrix_ptr
            = std::unique_ptr<gsl_matrix, std::function<void(gsl_matrix *)>>;
        using vector_ptr
            = std::unique_ptr<gsl_vector, std::function<void(gsl_vector *)>>;
        std::vector<double> effect_sizes, dominance_values;
        matrix_ptr input_matrix_copy;
        // Stores the Cholesky decomposition
        matrix_ptr matrix;
        // Stores the results of gsl_ran_multivariate_gaussian
        mutable gsl_vector_view res;
        // Stores the means of the mvn distribution, which are all zero
        vector_ptr mu;
        double fixed_effect, dominance;

        MultivariateGaussianEffects(const Region &r, const double sc,
                                    const gsl_matrix &input_matrix, double s, double h)
            : Sregion(r, sc, input_matrix.size1), effect_sizes(input_matrix.size1),
              dominance_values(input_matrix.size1, h),
              input_matrix_copy(gsl_matrix_alloc(input_matrix.size1, input_matrix.size2),
                                [](gsl_matrix *m) { gsl_matrix_free(m); }),
              matrix(gsl_matrix_alloc(input_matrix.size1, input_matrix.size2),
                     [](gsl_matrix *m) { gsl_matrix_free(m); }),
              // Holds the results of calls to mvn, and maps the
              // output to effect_sizes
              res(gsl_vector_view_array(effect_sizes.data(), effect_sizes.size())),
              // NOTE: use of calloc to initialize mu to all zeros
              mu(gsl_vector_calloc(input_matrix.size1),
                 [](gsl_vector *v) { gsl_vector_free(v); }),
              fixed_effect(s), dominance(h)
        {
            if (!std::isfinite(fixed_effect))
                {
                    throw std::invalid_argument("fixed_effect must be finite");
                }
            if (!std::isfinite(dominance))
                {
                    throw std::invalid_argument("dominance must be finite");
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

            // Assign the matrix and do the Cholesky decomposition
            auto error_handler = gsl_set_error_handler_off();

            int rv = gsl_matrix_memcpy(matrix.get(), &input_matrix);
            if (rv != GSL_SUCCESS)
                {
                    // Reset error handler on the way out
                    gsl_set_error_handler(error_handler);
                    throw std::runtime_error("failure copying input matrix");
                }
            rv = gsl_matrix_memcpy(input_matrix_copy.get(), &input_matrix);
            if (rv != GSL_SUCCESS)
                {
                    // Reset error handler on the way out
                    gsl_set_error_handler(error_handler);
                    throw std::runtime_error("failure copying input matrix");
                }
            rv = gsl_linalg_cholesky_decomp1(matrix.get());
            if (rv == GSL_EDOM)
                {
                    // Reset error handler on the way out
                    gsl_set_error_handler(error_handler);
                    throw std::invalid_argument("Cholesky decomposition failed");
                }
            // Reset error handler on the way out
            gsl_set_error_handler(error_handler);
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
            if (rv != GSL_SUCCESS)
                {
                    throw std::runtime_error(
                        "call to gsl_ran_multivariate_gaussian failed");
                }
            return infsites_Mutation(
                recycling_bin, mutations, lookup_table, false, generation,
                [this, &rng]() { return region(rng); },
                [this]() { return fixed_effect; }, [this]() { return dominance; },
                [this]() { return effect_sizes; }, [this]() { return dominance_values; },
                this->label());
        }

        double
        from_mvnorm(const double deviate, const double /*P*/) const override
        {
            return deviate / scaling;
        }

        std::vector<double>
        get_dominance() const override
        {
            return dominance_values;
        }
    };
} // namespace fwdpy11

#endif
