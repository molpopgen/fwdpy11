#ifndef FWDPY11_GAMMAS_HPP
#define FWDPY11_GAMMAS_HPP

#include <cmath>
#include <stdexcept>
#include <fwdpy11/policies/mutation.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{

    struct GammaS : public Sregion
    {
        double mean, shape, dominance;

        GammaS(const Region& r, double sc, double m, double s, double h)
            : Sregion(r, sc), mean(m), shape(s), dominance(h)
        {
            if (!std::isfinite(mean))
                {
                    throw std::invalid_argument("mean must be finite");
                }
            if (!std::isfinite(shape))
                {
                    throw std::invalid_argument("shape must be finite");
                }
            if (!std::isfinite(dominance))
                {
                    throw std::invalid_argument("dominance must be finite");
                }
        }

        std::unique_ptr<Sregion>
        clone() const
        {
            return std::unique_ptr<GammaS>(new GammaS(*this));
        }

        std::string
        repr() const
        {
            std::ostringstream out;
            out.precision(4);
            out << "GammaS(";
            this->region.region_repr(out);
            out << ", mean=" << this->mean << ", shape=" << this->shape
                << ", h=" << this->dominance << ", scaling=" << this->scaling
                << ')';
            return out.str();
        }
        std::uint32_t
        operator()(
            fwdpp::flagged_mutation_queue& recycling_bin,
            std::vector<Mutation>& mutations,
            std::unordered_multimap<double, std::uint32_t>& lookup_table,
            const std::uint32_t generation, const GSLrng_t& rng) const
        {
            return infsites_Mutation(
                recycling_bin, mutations, lookup_table, generation,
                [this, &rng]() { return region(rng); },
                [this, &rng]() {
                    return gsl_ran_gamma(rng.get(), shape, mean / shape)
                           / scaling;
                },
                [this]() { return dominance; }, this->label());
        }

        pybind11::tuple
        pickle() const
        {
            return pybind11::make_tuple(Sregion::pickle_Sregion(), mean, shape,
                                        dominance);
        }

        static GammaS
        unpickle(pybind11::tuple t)
        {
            if (t.size() != 4)
                {
                    throw std::runtime_error("invalid tuple size");
                }
            auto base = t[0].cast<pybind11::tuple>();
            return GammaS(Region::unpickle(base[0]), base[1].cast<double>(),
                          t[1].cast<double>(), t[2].cast<double>(),
                          t[3].cast<double>());
        }
    };
} // namespace fwdpy11

#endif
