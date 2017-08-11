#ifndef FWDPY11_MULTILOCUS_HPP__
#define FWDPY11_MULTILOCUS_HPP__

#include <functional>
#include <type_traits>
#include <gsl/gsl_randist.h>
#include <fwdpy11/rng.hpp>

namespace fwdpy11
{
    struct interlocus_rec
    {
        enum RECMODEL : int
        {
            BINOMIAL,
            POISSON
        };
        using mtype = std::underlying_type<RECMODEL>::type;
        const double param;
        RECMODEL m;
        interlocus_rec(const double param_, const mtype m_)
            : param(param_), m(static_cast<RECMODEL>(m_))
        {
        }

        std::function<unsigned(void)>
        callback(const fwdpy11::GSLrng_t& rng) const
        {
            if (m == RECMODEL::BINOMIAL)
                {
                    return std::bind(gsl_ran_binomial, rng.get(), param, 1);
                }
            return std::bind(gsl_ran_poisson, rng.get(), param);
        }

        mtype get_model() const
        {
            return static_cast<mtype>(m);
        }
    };
}
#endif
