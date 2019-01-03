#ifndef FWDPY11_MULTIVARIATEGAUSSIAN_HPP
#define FWDPY11_MULTIVARIATEGAUSSIAN_HPP

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <stdexcept>
#include <functional>
#include <memory>
#include "Sregion.hpp"

namespace fwdpy11
{
    struct MultivariateGaussianEffects : public Sregion
    {
        using matrix_ptr
            = std::unique_ptr<gsl_matrix, std::function<void(gsl_matrix *)>>;
        using vector_ptr
            = std::unique_ptr<gsl_vector, std::function<void(gsl_vector *)>>;
        // Stores the Cholesky decomposition
        matrix_ptr matrix;
        // Stores the results of gsl_ran_multivariate_gaussian
        vector_ptr res;
        double fixed_effect, dominance;

        MultivariateGaussianEffects(double beg, double end, double weight,
                                    bool coupled,
                                    const gsl_matrix &input_matrix, double s,
                                    double h, std::uint16_t label,
                                    bool matrix_is_covariance)
            : Sregion(beg, end, weight, coupled, label, 1.0),
              matrix(gsl_matrix_alloc(input_matrix.size1, input_matrix.size2),
                     [](gsl_matrix *m) { gsl_matrix_free(m); }),
              res(gsl_vector_alloc(input_matrix.size1),
                  [](gsl_vector *v) { gsl_vector_free(v); }),
              fixed_effect(s), dominance(h)
        // If matrix_is_covariance is true, then the input_matrix is treated
        // as a covariance matrix, meaning that we copy it and store its
        // Cholesky decomposition.  If matrix_is_covariance is false,
        // then input_matrix is assumed to be a valid Cholesky decomposition.
        {
            if (!std::isfinite(dominance))
                {
                    throw std::invalid_argument("dominance must be finite");
                }
            if (matrix->size1 != matrix->size2)
                {
                    throw std::invalid_argument("input matrix must be square");
                }

            // Assign the matrix and do the Cholesky decomposition
            auto error_handler = gsl_set_error_handler_off();

            int rv = gsl_matrix_memcpy(matrix.get(), &input_matrix);
            if (rv != GSL_SUCCESS)
                {
                    throw std::runtime_error("failure copying input matrix");
                }
            if (matrix_is_covariance)
                {
                    rv = gsl_linalg_cholesky_decomp1(matrix.get());
                    if (rv == GSL_EDOM)
                        {
                            throw std::invalid_argument(
                                "Cholesky decomposition failed");
                        }
                }

            // Reset error handler on the way out
            gsl_set_error_handler(error_handler);
        }

        virtual std::unique_ptr<Sregion>
        clone() const
        {
            return std::unique_ptr<MultivariateGaussianEffects>(
                new MultivariateGaussianEffects(
                    this->beg(), this->end(), this->weight(),
                    this->region.coupled, *matrix.get(), this->fixed_effect,
                    this->dominance, this->label(), false));
        }

        virtual std::uint32_t
        operator()(
            fwdpp::flagged_mutation_queue &recycling_bin,
            std::vector<Mutation> &mutations,
            std::unordered_multimap<double, std::uint32_t> &lookup_table,
            const std::uint32_t generation, const GSLrng_t &rng) const
        {
            return 0;
        }
    };
} // namespace fwdpy11

#endif
