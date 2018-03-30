#ifndef FWDPY11_MUTATION_TYPE_HPP__
#define FWDPY11_MUTATION_TYPE_HPP__

/*
 * This file started off via a copy
 * of fwdpp's popgenmut.hpp
 * by Kevin Thornton.
 */

#include <tuple>
#include <cstdint>
#include <vector>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/io/mutation.hpp>
#include <fwdpp/io/scalar_serialization.hpp>

namespace fwdpy11
{
    struct Mutation : public fwdpp::mutation_base
    ///! The fwdpy11 mutation type
    {
        //! The generation when the mutation arose
        fwdpp::uint_t g;
        //! Selection coefficient
        double s;
        //! Dominance of the mutation
        double h;
        std::vector<double> esizes, heffects;
        //! Alias for tuple type that can be used for object construction
        using constructor_tuple
            = std::tuple<double, double, double, unsigned, std::uint16_t>;

        /*!
          \brief Constructor
          \param __pos Mutation position
          \param __s Selection coefficient
          \param __h Dominance coefficient
          \param __g Generation when mutation arose
          \param __x Value to assign to mutation_base::xtra
        */
        Mutation(const double &__pos, const double &__s, const double &__h,
                 const unsigned &__g, const std::uint16_t x = 0) noexcept
            : mutation_base(__pos, (__s == 0.) ? true : false, x), g(__g),
              s(__s), h(__h)
        {
        }

        template <typename vectype>
        Mutation(const double &__pos, const double &__s, const double &__h,
                 const unsigned &__g, vectype &&esizes_, vectype &&heffects_,
                 const std::uint16_t x = 0) noexcept
            : fwdpp::mutation_base(__pos, (__s == 0.) ? true : false, x),
              g(__g), s(__s), h(__h), esizes(std::forward<vectype>(esizes_)),
              heffects(std::forward<vectype>(heffects_))
        {
        }

        Mutation(constructor_tuple t) noexcept
            : mutation_base(std::get<0>(t),
                            (std::get<1>(t) == 0.) ? true : false,
                            std::get<4>(t)),
              g(std::get<3>(t)), s(std::get<1>(t)), h(std::get<2>(t))
        {
        }

        bool
        operator==(const Mutation &rhs) const
        {
            return std::tie(this->g, this->s, this->h, this->esizes,
                            this->heffects)
                       == std::tie(rhs.g, rhs.s, rhs.h, rhs.esizes,
                                   rhs.heffects)
                   && is_equal(rhs);
        }
    };
}

// Make Mutation compatible with fwdpp's
// serialization API:
namespace fwdpp
{
    namespace io
    {
        template <> struct serialize_mutation<Mutation>
        {
            io::scalar_writer writer;
            serialize_mutation<Mutation>() : writer{} {}
            template <typename streamtype>
            inline void
            operator()(streamtype &buffer, const Mutation &m) const
            {
                writer(buffer, &m.g);
                writer(buffer, &m.pos);
                writer(buffer, &m.s);
                writer(buffer, &m.h);
                writer(buffer, &m.xtra);
                std::size_t ns = m.esizes.size(), nh = m.heffects.size();
                writer(buffer, &ns);
                writer(buffer, &nh);
                if (ns)
                    {
                        writer(buffer, m.esizes.data(), ns);
                    }
                if (nh)
                    {
                        writer(buffer, m.heffects.data(), nh);
                    }
            }
        };

        template <> struct deserialize_mutation<Mutation>
        {
            io::scalar_reader reader;
            deserialize_mutation<Mutation>() : reader{} {}
            template <typename streamtype>
            inline Mutation
            operator()(streamtype &buffer) const
            {
                uint_t g;
                double pos, s, h;
                decltype(Mutation::xtra) xtra;
                io::scalar_reader reader;
                reader(buffer, &g);
                reader(buffer, &pos);
                reader(buffer, &s);
                reader(buffer, &h);
                reader(buffer, &xtra);
                std::size_t ns, nh;
                reader(buffer, &ns);
                reader(buffer, &nh);
                std::vector<double> ss, hs;
                if (ns)
                    {
                        ss.resize(ns);
                        reader(buffer, ss.data(), ns);
                    }
                if (nh)
                    {
                        hs.resize(ns);
                        reader(buffer, hs.data(), nh);
                    }
                return Mutation(pos, s, h, g, std::move(ss), std::move(hs),
                                xtra);
            }
        };
    }
}
#endif
