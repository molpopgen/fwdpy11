#ifndef FWDPY11_TYPES_DIPLOID_HPP__
#define FWDPY11_TYPES_DIPLOID_HPP__

#include <cstdint>
#include <cstddef>
#include <vector>
#include <tuple>
#include <fwdpp/ts/definitions.hpp>

namespace fwdpy11
{
    struct DiploidGenotype
    /// Kept as "bare"/simple struct so that
    /// we can use this as a Numpy dtype.
    /// Direct initialization of data members must
    /// use construction via brace/"uniform"
    /// initialization.
    {
        using first_type = std::size_t;
        using second_type = std::size_t;
        first_type first;
        second_type second;
    };

    // typedef to ease refactoring in this PR
    // TODO: remove this
    using Diploid = DiploidGenotype;

    inline bool
    operator==(const DiploidGenotype& a, const DiploidGenotype& b)
    {
        return a.first == b.first && a.second == b.second;
    }

    struct DiploidMetadata
    /*!
      \brief Bare struct for metadata.  Used as Numpy dtype.
    */
    {
        double g;               // Genetic value
        double e;               // Random component of trait value
        double w;               // Fitness
        double geography[3];    // Location in geographic space
        std::size_t label;      // Index of individual in pop container
        std::size_t parents[2]; // Indexes of parents
        std::int32_t deme;
        std::int32_t sex;
        std::int32_t nodes[2]; // Nodes in TreeSequence
    };

    struct ancient_sample_record
    /*! When tracking ancient samples, 
     * We want to be able to provide access to
     * the individual's metadata PLUS
     * when they existed and what nodes 
     * on the tree they are.
     */
    {
        double time;
        fwdpp::ts::TS_NODE_INT n1, n2;
    };

    //! Typedef for container of diploids
    using dipvector_t = std::vector<DiploidGenotype>;
} // namespace fwdpy11

#endif
