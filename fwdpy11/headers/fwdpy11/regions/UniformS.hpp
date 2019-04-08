#ifndef FWDPY11_UNIFORMS_HPP
#define FWDPY11_UNIFORMS_HPP

#include <cmath>
#include <stdexcept>
#include <fwdpy11/policies/mutation.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{

    struct UniformS : public Sregion
    {
        double lo, hi, dominance;

        UniformS(const Region& r, double sc, double lo_, double hi_, double h)
            : Sregion(r, sc), lo(lo_), hi(hi_), dominance(h)
        {
            if (!std::isfinite(lo))
                {
                    throw std::invalid_argument("lo must be finite");
                }
            if (!std::isfinite(hi))
                {
                    throw std::invalid_argument("hi must be finite");
                }
            if (!std::isfinite(dominance))
                {
                    throw std::invalid_argument("dominance must be finite");
                }
            if (!(hi > lo))
                {
                    throw std::invalid_argument("hi must be > lo");
                }
        }

        std::unique_ptr<Sregion>
        clone() const
        {
            return std::unique_ptr<UniformS>(new UniformS(*this));
        }

        std::string
        repr() const
        {
            std::ostringstream out;
            out.precision(4);
            out << "UniformS(";
            this->region.region_repr(out);
            out << ", lo=" << this->lo << ", hi=" << this->hi
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
                    return gsl_ran_flat(rng.get(), lo, hi) / scaling;
                },
                [this]() { return dominance; }, this->label());
        }
        pybind11::tuple
        pickle() const
        {
            return pybind11::make_tuple(Sregion::pickle_Sregion(), lo, hi,
                                        dominance);
        }

        static UniformS
        unpickle(pybind11::tuple t)
        {
            if (t.size() != 4)
                {
                    throw std::runtime_error("invalid tuple size");
                }
            auto base = t[0].cast<pybind11::tuple>();
            return UniformS(Region::unpickle(base[0]), base[1].cast<double>(),
                            t[1].cast<double>(), t[2].cast<double>(),
                            t[3].cast<double>());
        }
    };
} // namespace fwdpy11

#endif

