#ifndef FWDPY11_GAUSSIANS_HPP
#define FWDPY11_GAUSSIANS_HPP

#include <cmath>
#include <stdexcept>
#include <fwdpy11/policies/mutation.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{

    struct GaussianS : public Sregion
    {
        double sd, dominance;

        GaussianS(const Region& r, double sc, double sd_, double h)
            : Sregion(r, sc), sd(sd_), dominance(h)
        {
            if (!std::isfinite(sd))
                {
                    throw std::invalid_argument("sd must be finite");
                }
            if (!(sd > 0))
                {
                    throw std::invalid_argument("sd must be > 0");
                }
            if (!std::isfinite(dominance))
                {
                    throw std::invalid_argument("dominance must be finite");
                }
        }

        std::unique_ptr<Sregion>
        clone() const
        {
            return std::unique_ptr<GaussianS>(new GaussianS(*this));
        }

        std::string
        repr() const
        {
            std::ostringstream out;
            out.precision(4);
            out << "GaussianS(";
            this->region.region_repr(out);
            out << ", sd=" << this->sd << ", h=" << this->dominance
                << ", scaling=" << this->scaling << ')';
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
                    return gsl_ran_gaussian_ziggurat(rng.get(), sd) / scaling;
                },
                [this]() { return dominance; }, this->label());
        }

        pybind11::tuple
        pickle() const
        {
            return pybind11::make_tuple(Sregion::pickle_Sregion(), sd,
                                        dominance);
        }

        static GaussianS
        unpickle(pybind11::tuple t)
        {
            if (t.size() != 3)
                {
                    throw std::runtime_error("invalid tuple size");
                }
            auto base = t[0].cast<pybind11::tuple>();
            return GaussianS(Region::unpickle(base[0]), base[1].cast<double>(),
                             t[1].cast<double>(), t[2].cast<double>());
        }
    };
} // namespace fwdpy11

#endif

