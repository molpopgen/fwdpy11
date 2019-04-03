#ifndef FWDPY11_REGION_HPP
#define FWDPY11_REGION_HPP

#include <gsl/gsl_randist.h>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <cstdint>
#include <fwdpy11/rng.hpp>
#include <pybind11/pybind11.h>

namespace fwdpy11
{
    struct Region
    {
        double beg, end, weight;
        std::uint16_t label;
        bool coupled;
        Region(double b, double e, double w, bool c, std::uint16_t l)
            : beg(b), end(e), weight((c == true) ? (end - beg) * w : w),
              label(l), coupled(c)
        {
            if (!std::isfinite(beg))
                {
                    throw std::invalid_argument("beg must be finite");
                }
            if (!std::isfinite(end))
                {
                    throw std::invalid_argument("end must be finite");
                }
            if (!std::isfinite(weight))
                {
                    throw std::invalid_argument("weight must be finite");
                }
            if (weight < 0.0)
                {
                    throw std::invalid_argument("weight must be non-negative");
                }
            if (!(end > beg))
                {
                    throw std::invalid_argument(
                        "end must be greater than beg");
                }
        }

        inline double
        operator()(const GSLrng_t& rng) const
        {
            return gsl_ran_flat(rng.get(), beg, end);
        }

        pybind11::tuple
        pickle() const
        {
            return pybind11::make_tuple(beg, end, weight, coupled, label);
        }

        static Region
        unpickle(pybind11::tuple t)
        {
            if (t.size() != 5)
                {
                    throw std::runtime_error("invalid tuple size");
                }
            return Region(t[0].cast<double>(), t[1].cast<double>(),
                          t[2].cast<double>(), t[3].cast<bool>(),
                          t[4].cast<std::uint16_t>());
        }

        void
        region_repr(std::ostringstream& o) const
        {
            o << "beg=" << this->beg << ", end=" << this->end
              << ", weight=" << this->weight;
        }

        std::string
        repr() const
        {
            std::ostringstream o;
            o.precision(4);
            o << "Region(";
            this->region_repr(o);
            o << ")";
            return o.str();
        }
    };

    inline bool
    operator==(const Region& lhs, const Region& rhs)
    {
        return lhs.beg == rhs.beg && lhs.end == rhs.end
               && lhs.weight == rhs.weight && lhs.coupled == rhs.coupled
               && lhs.label == rhs.label;
    }
} // namespace fwdpy11

#endif
