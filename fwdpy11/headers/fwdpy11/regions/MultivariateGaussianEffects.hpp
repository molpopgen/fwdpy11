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
        // Stores the Cholesky decomposition
        matrix_ptr matrix;
        // Stores the results of gsl_ran_multivariate_gaussian
        mutable gsl_vector_view res;
        // Stores the means of the mvn distribution, which are all zero
        vector_ptr mu;
        double fixed_effect, dominance;

        MultivariateGaussianEffects(const Region &r, const double sc,
                                    const gsl_matrix &input_matrix, double s,
                                    // NOTE: matrix_is_covariance is
                                    // NOT exposed to Python
                                    double h, bool matrix_is_covariance)
            : Sregion(r, sc), effect_sizes(input_matrix.size1),
              dominance_values(input_matrix.size1, h),
              matrix(gsl_matrix_alloc(input_matrix.size1, input_matrix.size2),
                     [](gsl_matrix *m) { gsl_matrix_free(m); }),
              // Holds the results of calls to mvn, and maps the
              // output to effect_sizes
              res(gsl_vector_view_array(effect_sizes.data(),
                                        effect_sizes.size())),
              // NOTE: use of calloc to initialize mu to all zeros
              mu(gsl_vector_calloc(input_matrix.size1),
                 [](gsl_vector *v) { gsl_vector_free(v); }),
              fixed_effect(s), dominance(h)
        // If matrix_is_covariance is true, then the input_matrix is treated
        // as a covariance matrix, meaning that we copy it and store its
        // Cholesky decomposition.  If matrix_is_covariance is false,
        // then input_matrix is assumed to be a valid Cholesky decomposition.
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
                            input_matrix.data
                                + input_matrix.size1 * input_matrix.size2,
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
            if (matrix_is_covariance)
                {
                    rv = gsl_linalg_cholesky_decomp1(matrix.get());
                    if (rv == GSL_EDOM)
                        {
                            // Reset error handler on the way out
                            gsl_set_error_handler(error_handler);
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
                    fwdpy11::Region(this->beg(), this->end(), this->weight(),
                                    this->region.coupled, this->label()),
                    1.0, *matrix.get(), this->fixed_effect, this->dominance,
                    false));
        }
        std::string
        repr() const
        {
            std::ostringstream out;
            out.precision(4);
            out << "MultivariateGaussianEffects(";
            this->region.region_repr(out);
            out << ", s=" << this->fixed_effect << ", h=" << this->dominance
            << ", matrix at " << matrix.get() << ')';
            return out.str();
        }

        virtual std::uint32_t
        operator()(
            fwdpp::flagged_mutation_queue &recycling_bin,
            std::vector<Mutation> &mutations,
            std::unordered_multimap<double, std::uint32_t> &lookup_table,
            const std::uint32_t generation, const GSLrng_t &rng) const
        {
            int rv = gsl_ran_multivariate_gaussian(rng.get(), mu.get(),
                                                   matrix.get(), &res.vector);
            if (rv != GSL_SUCCESS)
                {
                    throw std::runtime_error(
                        "call to gsl_ran_multivariate_gaussian failed");
                }
            return infsites_Mutation(
                recycling_bin, mutations, lookup_table, generation,
                [this, &rng]() { return region(rng); },
                [this]() { return fixed_effect; },
                [this]() { return dominance; },
                [this]() { return effect_sizes; },
                [this]() { return dominance_values; }, this->label());
        }

        pybind11::tuple
        pickle() const
        {
            pybind11::list matrix_data;
            for (std::size_t i = 0; i < matrix->size1; ++i)
                {
                    for (std::size_t j = 0; j < matrix->size2; ++j)
                        {
                            matrix_data.append(
                                gsl_matrix_get(matrix.get(), i, j));
                        }
                }

            return pybind11::make_tuple(Sregion::pickle_Sregion(), matrix_data,
                                        matrix->size1, matrix->size2,
                                        fixed_effect, dominance);
        }

        static MultivariateGaussianEffects
        unpickle(pybind11::tuple t)
        {
            if (t.size() != 6)
                {
                    throw std::runtime_error("invalid tuple size");
                }
            auto base = t[0].cast<pybind11::tuple>();
            std::vector<double> input_matrix_data;
            pybind11::list input_matrix_list = t[1].cast<pybind11::list>();
            std::size_t size1 = t[2].cast<std::size_t>();
            std::size_t size2 = t[3].cast<std::size_t>();
            for (auto i : input_matrix_list)
                {
                    input_matrix_data.push_back(i.cast<double>());
                }
            auto v = gsl_matrix_const_view_array(input_matrix_data.data(),
                                                 size1, size2);
            return MultivariateGaussianEffects(
                Region::unpickle(base[0]), base[1].cast<double>(), v.matrix,
                t[4].cast<double>(), t[5].cast<double>(), false);
        }
    };

    bool
    operator==(const MultivariateGaussianEffects &lhs,
               const MultivariateGaussianEffects &rhs)
    {
        bool base_equal = lhs.is_equal(rhs);
        if (base_equal == false)
            {
                return false;
            }
        return rhs.fixed_effect == rhs.fixed_effect
               && rhs.dominance == rhs.dominance
               && gsl_matrix_equal(rhs.matrix.get(), rhs.matrix.get());
    }

} // namespace fwdpy11

#endif
