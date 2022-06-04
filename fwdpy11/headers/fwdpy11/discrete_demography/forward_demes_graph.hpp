#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace fwdpy11
{
    namespace discrete_demography
    {
        /* General comments:
         * * These types should go into separate headers.
         * * And/or...we can probably pimpl them with 
         *   no (noticeable) loss of performance.
         */

        // fwdpy11 time is uint32_t, but that is
        // BAD if you have to do stuff like subtract!
        // So, we use int64_t here and will return casted
        // values as needed (after checks).
        using demes_model_time = std::int64_t;

        struct Pulse
        {
            std::vector<std::int32_t> sources;
            std::vector<double> proportions;
            std::int32_t dest;
            demes_model_time time;
        };

        struct AsymmetricMigration
        {
            std::int32_t source, dest;
            demes_model_time start_time, end_time;
        };

        // Should this be a virtual base class
        // or simply store a std::function?
        struct SizeFunction
        {
            virtual double operator()(demes_model_time epoch_start_time,
                                      demes_model_time epoch_end_time,
                                      demes_model_time current_time) const;
        };

        struct Epoch
        {
            demes_model_time end_time, start_size, end_size;
            double cloning_rate, selfing_rate;

            // Never accept nullptr: "constant" should
            // also be a callback.
            std::unique_ptr<SizeFunction> size_function;
        };

        struct Deme
        {
            std::int32_t id, start_time;
            std::vector<std::int32_t> ancestors;
            std::vector<double> proportions;
            std::vector<Epoch> epochs;
        };

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
