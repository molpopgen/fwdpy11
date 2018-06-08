#ifndef FWDPY11_TYPES_DIPLOID_HPP__
#define FWDPY11_TYPES_DIPLOID_HPP__

#include <cstdint>
#include <cstddef>
#include <vector>
#include <tuple>

namespace fwdpy11
{
    using Diploid = std::pair<std::size_t, std::size_t>;

    struct diploid_metadata
    /*!
      \brief Bare struct for metadata.  Used as Numpy dtype.
    */
    {
        double g; // Genetic value
        double e; // Random component of trait value
        double w; // Fitness
        std::size_t label;  // Index of individual in pop container
        std::size_t parents[2]; // Indexes of parents
        std::uint32_t deme;
        std::int32_t sex;
    };

    //! Typedef for container of diploids
    using dipvector_t = std::vector<Diploid>;
}

#endif
