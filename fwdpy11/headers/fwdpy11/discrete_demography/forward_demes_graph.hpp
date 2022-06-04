#pragma once

#include "forward_demes_graph/demes_model_time.hpp"
#include "forward_demes_graph/asymmetric_migration.hpp"
#include "forward_demes_graph/size_function.hpp"
#include "forward_demes_graph/epoch.hpp"
#include "forward_demes_graph/deme.hpp"
#include "forward_demes_graph/pulse.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        struct ForwardDemesGraph
        /// Limited to:
        /// * construction from a fully-resolved demes.Graph
        ///
        /// Needs:
        ///
        /// * A map of deme id -> name
        /// * Methods to set the current time point
        ///   via the pop.generation value.
        /// * Methods to check compatability of a
        ///   model at a current time point and
        ///   a population. These methods cannot
        ///   depend on knowing about the populaiton
        ///   class! Perhaps they should live elsewhere/
        ///   be uncoupled.
        /// * A method to determine single-generation
        ///   ancestry proportions for any destination
        ///   deme at the current_time.
        {
            std::int64_t current_time; // never public.  just for internal book-keeping.
            std::vector<Deme> demes;
            std::vector<Pulse> pulses;
            std::vector<AsymmetricMigration> migrations;
        };
    }
}
