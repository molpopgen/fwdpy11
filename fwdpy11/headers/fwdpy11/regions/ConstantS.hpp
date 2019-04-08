#ifndef FWDPY11_CONSTANTS_HPP
#define FWDPY11_CONSTANTS_HPP

#include <cmath>
#include <stdexcept>
#include <fwdpy11/policies/mutation.hpp>
#include "Sregion.hpp"

namespace fwdpy11
{

    struct ConstantS : public Sregion
    {
        double esize, dominance;

        ConstantS(const Region& r, const double s, const double es,
                  const double h)
            : Sregion(r, s), esize(es), dominance(h)
        {
            if (!std::isfinite(esize))
                {
                    throw std::invalid_argument("esize must be finite");
                }
            if (!std::isfinite(dominance))
                {
                    throw std::invalid_argument("dominance must be finite");
                }
        }

        std::unique_ptr<Sregion>
        clone() const
        {
            return std::unique_ptr<ConstantS>(new ConstantS(*this));
        }

        std::string
        repr() const
        {
            std::ostringstream out;
            out.precision(4);
            out << "ConstantS(";
            this->region.region_repr(out);
            out << ", s=" << this->esize << ", h=" << this->dominance
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
                [this]() { return esize / scaling; },
                [this]() { return dominance; }, this->label());
        }

        pybind11::tuple
        pickle() const
        {
            return pybind11::make_tuple(Sregion::pickle_Sregion(), esize,
                                        dominance);
        }

        static ConstantS
        unpickle(pybind11::tuple t)
        {
            if (t.size() != 3)
                {
                    throw std::runtime_error("invalid tuple size");
                }
            auto base = t[0].cast<pybind11::tuple>();
            return ConstantS(Region::unpickle(base[0]), base[1].cast<double>(),
                             t[1].cast<double>(), t[2].cast<double>());
        }
    };
} // namespace fwdpy11

#endif

