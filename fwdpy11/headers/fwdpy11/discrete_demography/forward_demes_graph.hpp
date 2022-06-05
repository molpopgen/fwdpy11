#pragma once

#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <limits>
#include <unordered_map>
#include <unordered_set>
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
            std::unordered_map<std::string, std::int32_t> deme_name_to_id;

            ForwardDemesGraph()
                : current_time{std::numeric_limits<demes_model_time>::max()}, demes{},
                  pulses{}, migrations{}, deme_name_to_id()
            {
            }

            Deme&
            add_deme(const std::string& name, std::int32_t id,
                     demes_model_time start_time)
            {
                if (deme_name_to_id.find(name) == end(deme_name_to_id))
                    {
                        deme_name_to_id.insert({name, id});
                        demes.push_back(Deme{id, start_time});
                    }
                else
                    {
                        std::ostringstream message;
                        message << "deme name " << name << " already exists";
                        throw std::invalid_argument(message.str());
                    }
                if (demes.empty())
                    {
                        throw std::runtime_error("fatal error: demes is empty");
                    }
                return demes.back();
            }
        };
    }
}
