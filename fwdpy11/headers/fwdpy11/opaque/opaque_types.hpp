#ifndef FWDPY11_OPAQUE_TYPES_HPP__
#define FWDPY11_OPAQUE_TYPES_HPP__

#include <pybind11/stl.h>
#include <cstdint>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/sugar/popgenmut.hpp>

namespace fwdpy11
{
    //! Typedef for mutation container
    using mcont_t = std::vector<KTfwd::popgenmut>;
    //! Typedef for gamete type
    using gamete_t = KTfwd::gamete;
    //! Typedef for gamete container
    using gcont_t = std::vector<gamete_t>;

    struct diploid_t
    /*!
      \brief Custom diploid type.
    */
    {
        using first_type = std::size_t;
        using second_type = std::size_t;
        //! First gamete.  A gamete is vector<size_t> where the elements are
        //! indexes to a population's gamete container
        first_type first;
        //! Second gamete. A gamete is vector<size_t> where the elements are
        //! indexes to a population's gamete container
        second_type second;
        //! 64 bits of data to do stuff with.  Initialized to zero upon
        //! construction
        std::size_t label;
        //! Genetic component of trait value.  This is not necessarily written
        //! to by a simulation.
        double g;
        //! Random component of trait value.  This is not necessarily written
        //! to by a simulation.
        double e;
        //! Fitness.  This is not necessarily written to by a simulation.
        double w;
        //! Constructor
        diploid_t() noexcept : first(first_type()),
                               second(second_type()),
                               label(0),
                               g(0.),
                               e(0.),
                               w(1.)
        {
        }
        //! Construct from two indexes to gametes
        diploid_t(first_type g1, first_type g2) noexcept : first(g1),
                                                           second(g2),
                                                           label(0),
                                                           g(0.),
                                                           e(0.),
                                                           w(1.)
        {
        }

        inline bool
        operator==(const diploid_t& dip) const noexcept
        //! Required for py::bind_vector
        {
            return this->first == dip.first && this->second == dip.second
                   && this->w == dip.w && this->g == dip.e
                   && this->label == dip.label;
        }
    };

    //! Typedef for container of diploids
    using dipvector_t = std::vector<diploid_t>;
}

PYBIND11_MAKE_OPAQUE(std::vector<KTfwd::uint_t>);
PYBIND11_MAKE_OPAQUE(fwdpy11::dipvector_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::gcont_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::mcont_t);

#endif
