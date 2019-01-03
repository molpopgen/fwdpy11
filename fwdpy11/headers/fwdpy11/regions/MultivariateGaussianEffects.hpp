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
        double dominance;

        MultivariateGaussianEffects(double beg, double end, double weight,
                                    bool coupled, gsl_matrix_const_view &view,
                                    double h, std::uint16_t label)
            : Sregion(beg, end, weight, coupled, label, 1.0),
              matrix(gsl_matrix_alloc(view.matrix.size1, view.matrix.size2),
                     [](gsl_matrix *m) { gsl_matrix_free(m); }),
              res(gsl_vector_alloc(view.matrix.size1),
                  [](gsl_vector *v) { gsl_vector_free(v); }),
              dominance(h)
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

            int rv = gsl_matrix_memcpy(matrix.get(), &view.matrix);
            if (rv != GSL_SUCCESS)
                {
                    throw std::runtime_error("failure copying input matrix");
                }
            rv = gsl_linalg_cholesky_decomp1(matrix.get());
            if (rv == GSL_EDOM)
                {
                    throw std::runtime_error("Cholesky decomposition failed");
                }

            // Reset error handler on the way out
            gsl_set_error_handler(error_handler);
        }

        virtual std::unique_ptr<Sregion>
        clone() const
        {
            gsl_matrix_const_view v = gsl_matrix_const_submatrix(
                matrix.get(), 0, 0, matrix->size1, matrix->size2);
            return std::unique_ptr<MultivariateGaussianEffects>(
                new MultivariateGaussianEffects(
                    this->beg(), this->end(), this->weight(),
                    this->region.coupled, v, this->dominance, this->label()));
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
