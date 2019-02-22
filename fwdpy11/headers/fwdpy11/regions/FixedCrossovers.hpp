#ifndef FWDPY11_REGIONS_FIXEDCROSSOVERS_HPP
#define FWDPY11_REGIONS_FIXEDCROSSOVERS_HPP

#include <cmath>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include "GeneticMapUnit.hpp"

namespace fwdpy11
{
    struct FixedCrossovers : public GeneticMapUnit
    {
        double beg, end;
        int nxovers;
        FixedCrossovers(double b, double e, int n) : beg(b), end(e), nxovers(n)
        {
            if (!std::isfinite(b))
                {
                    throw std::invalid_argument("beg must be finite");
                }
            if (!std::isfinite(e))
                {
                    throw std::invalid_argument("end must be finite");
                }
            if (e <= b)
                {
                    throw std::invalid_argument(
                        "end must be greater than beg");
                }
            if (nxovers < 0)
                {
                    throw std::invalid_argument(
                        "number of crossovers must be non-negative");
                }
        }

        void
        operator()(const GSLrng_t& rng,
                   std::vector<double>& breakpoints) const final
        {
            for (int i = 0; i < nxovers; ++i)
                {
                    breakpoints.push_back(gsl_ran_flat(rng.get(), beg, end));
                }
        }

        pybind11::object
        pickle() const final
        {
            return pybind11::make_tuple(beg, end, nxovers);
        }

        std::unique_ptr<GeneticMapUnit>
        clone() const final
        {
            return std::unique_ptr<GeneticMapUnit>(new FixedCrossovers(*this));
        }

        static FixedCrossovers
        unpickle(pybind11::object o)
        {
            auto t = o.cast<pybind11::tuple>();
            if (t.size() != 3)
                {
                    throw std::runtime_error("invalid tuple size");
                }
            return FixedCrossovers(t[0].cast<double>(), t[1].cast<double>(),
                                   t[2].cast<int>());
        }
    };
} // namespace fwdpy11

#endif
