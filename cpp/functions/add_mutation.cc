//
// Copyright (C) 2021 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//
#include <cmath>
#include <limits>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>
#include <sstream>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/policies/mutation.hpp>
#include <core/internal/gsl_ran_flat.hpp>
#include <fwdpp/ts/marginal_tree_functions/nodes.hpp>
#include <fwdpp/ts/table_collection_functions.hpp>
#include <fwdpp/ts/marginal_tree_functions/samples.hpp>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <gsl/gsl_randist.h>

#include "add_mutation.hpp"

new_mutation_data::new_mutation_data(double e, double h, std::vector<double> esizes,
                                     std::vector<double> heffects,
                                     decltype(fwdpy11::Mutation::xtra) l)
    : effect_size{e}, dominance{h}, esizes{std::move(esizes)},
      heffects{std::move(heffects)}, label{l}
{
    if (e == 0.0
        && std::all_of(begin(esizes), end(esizes), [](double d) { return d == 0.0; }))
        {
            throw std::invalid_argument("new mutation must have non-zero effect size");
        }
    if (!std::isfinite(e) || std::any_of(begin(esizes), end(esizes), [](double d) {
            return !std::isfinite(d);
        }))
        {
            throw std::invalid_argument("all effect size values must be finite");
        }
    if (!std::isfinite(dominance)
        || std::any_of(begin(heffects), end(heffects),
                       [](double d) { return !std::isfinite(d); }))
        {
            throw std::invalid_argument("dominance values must all be finite");
        }
    if (esizes.size() != heffects.size())
        {
            throw std::invalid_argument(
                "effect size and dominance vectors must have the same length");
        }
}

namespace
{
    struct candidate_node_map
    {
        double left, right;
        // The candidate node and its parent id
        fwdpp::ts::table_index_t node, parent;
        // The sample nodes below node
        std::vector<fwdpp::ts::table_index_t> descendants;
        double node_time, parent_time, tree_span;

        candidate_node_map(double l, double r, fwdpp::ts::table_index_t n,
                           fwdpp::ts::table_index_t p,
                           std::vector<fwdpp::ts::table_index_t> d, double a, double b,
                           double c)
            : left{l}, right{r}, node{n}, parent{p},
              descendants{std::move(d)}, node_time{a}, parent_time{b}, tree_span{c}
        {
        }
    };

    bool
    node_has_valid_time(const std::vector<fwdpp::ts::node>& node_table,
                        const std::int32_t mutation_node,
                        const std::int32_t mutation_node_parent)
    {
        if (mutation_node_parent == fwdpp::ts::NULL_INDEX)
            {
                return false;
            }
        auto mutation_node_time = node_table[mutation_node].time;
        // Check for easy case where we can assign mutation time == node time
        if (mutation_node_time < 0.0)
            {
                if (mutation_node_time - std::floor(mutation_node_time) == 0.0)
                    {
                        return true;
                    }
            }
        else if (mutation_node_time - std::ceil(mutation_node_time) == 0.0)
            {
                return true;
            }

        auto mutation_node_parent_time = node_table[mutation_node_parent].time;
        bool valid = false;
        if (mutation_node_parent_time < 0.0)
            {
                if (std::ceil(mutation_node_parent_time) <= mutation_node_time)
                    {
                        valid = true;
                    }
            }
        else
            {
                if (std::floor(mutation_node_parent_time) <= mutation_node_time)
                    {
                        valid = true;
                    }
            }
        return valid;
    }

    std::vector<candidate_node_map>
    generate_canidate_list(const double left, const double right,
                           const unsigned ndescendants,
                           const fwdpp::ts::table_index_t deme,
                           const fwdpy11::DiploidPopulation& pop)
    {
        std::vector<fwdpp::ts::table_index_t> samples;
        for (auto& dip : pop.diploid_metadata)
            {
                for (auto i : dip.nodes)
                    {
                        samples.push_back(i);
                    }
            }

        auto tv = fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection>(
            *pop.tables, samples, fwdpp::ts::update_samples_list(true));

        std::vector<candidate_node_map> candidates;

        while (tv())
            {
                auto& tree = tv.tree();
                if (tree.left < right && tree.right >= left)
                    {
                        fwdpp::ts::node_iterator ni(tree, fwdpp::ts::nodes_preorder());
                        auto n = ni();
                        while (n != fwdpp::ts::NULL_INDEX)
                            {
                                if (tree.leaf_counts[n]
                                    == static_cast<fwdpp::ts::table_index_t>(
                                        ndescendants))
                                    {
                                        // Then this is a node in this tree
                                        // that is a candidate for mutation
                                        // placement
                                        std::vector<fwdpp::ts::table_index_t>
                                            descendants;
                                        fwdpp::ts::samples_iterator si(
                                            tree, n,
                                            fwdpp::ts::convert_sample_index_to_nodes(
                                                true));
                                        auto s = si();
                                        while (s != fwdpp::ts::NULL_INDEX)
                                            {
                                                descendants.push_back(s);
                                                s = si();
                                            }
                                        if (deme < 0
                                            || std::all_of(
                                                begin(descendants), end(descendants),
                                                [&pop,
                                                 deme](fwdpp::ts::table_index_t i) {
                                                    return pop.tables->nodes[i].deme
                                                           == deme;
                                                }))
                                            {
                                                if (node_has_valid_time(
                                                        pop.tables->nodes, n,
                                                        tree.parents[n]))
                                                    {
                                                        candidates.emplace_back(
                                                            std::max(tree.left, left),
                                                            std::min(tree.right, right),
                                                            n, tree.parents[n],
                                                            std::move(descendants),
                                                            pop.tables->nodes[n].time,
                                                            pop.tables
                                                                ->nodes[tree.parents[n]]
                                                                .time,
                                                            tree.right - tree.left);
                                                    }
                                            }
                                    }
                                n = ni();
                            }
                    }
                if (tree.left >= right)
                    {
                        break;
                    }
            }
        return candidates;
    }

    std::int64_t
    generate_mutation_time(const fwdpy11::GSLrng_t& rng,
                           const std::vector<fwdpp::ts::node>& node_table,
                           const std::int32_t mutation_node,
                           const std::int32_t mutation_node_parent)
    {
        auto mutation_node_time = node_table[mutation_node].time;
        // Check for easy case where we can assign mutation time == node time
        if (mutation_node_time < 0.0)
            {
                if (mutation_node_time - std::floor(mutation_node_time) == 0.0)
                    {
                        return static_cast<std::int64_t>(std::floor(mutation_node_time));
                    }
            }
        else if (mutation_node_time - std::ceil(mutation_node_time) == 0.0)
            {
                return static_cast<std::int64_t>(std::ceil(mutation_node_time));
            }

        if (mutation_node_parent == fwdpp::ts::NULL_INDEX)
            {
                // Time is the closest int64_t to the node time
                if (mutation_node_time < 0.0)
                    {
                        return static_cast<std::int64_t>(std::floor(mutation_node_time));
                    }
                return static_cast<std::int64_t>(std::ceil(mutation_node_time));
            }
        // choose a random time
        auto candidate_time = fwdpy11_core::internal::gsl_ran_flat(
            rng, mutation_node_time, node_table[mutation_node_parent].time);
        if (candidate_time < 0.0)
            {
                return static_cast<std::int64_t>(std::floor(mutation_node_time));
            }
        return static_cast<std::int64_t>(std::ceil(mutation_node_time));
    }
}

// Returns size_t max value if no candidates are found.
std::size_t
add_mutation(const fwdpy11::GSLrng_t& rng, const double left, const double right,
             const fwdpp::ts::table_index_t ndescendants,
             const fwdpp::ts::table_index_t deme, const new_mutation_data& data,
             fwdpy11::DiploidPopulation& pop)
{
    // TODO: tables must be indexed!
    if (pop.is_simulating)
        {
            throw std::runtime_error(
                "mutations cannot be added to a population during simulation");
        }
    if (right <= left)
        {
            throw std::invalid_argument("right must be > left");
        }
    if (left < 0.0 || left >= pop.tables->genome_length())
        {
            throw std::invalid_argument("left must be 0 <= x < genome length");
        }
    if (right <= 0. || right > pop.tables->genome_length())
        {
            throw std::invalid_argument("right must be 0 < x <= genome length");
        }
    if (!std::isfinite(right) || !std::isfinite(left))
        {
            throw std::invalid_argument("left and right must both be finite");
        }
    if (pop.tables->edges.empty())
        {
            throw std::invalid_argument(
                "a population must have ancestry in order to add mutations");
        }
    if (ndescendants < 1)
        {
            throw std::invalid_argument("number of descendants must be >= 1");
        }
    if (static_cast<decltype(pop.N)>(ndescendants) >= 2 * pop.N)
        {
            throw std::invalid_argument("ndescendants must be < 2*pop.N");
        }
    if (deme >= 0)
        {
            std::unordered_map<fwdpp::ts::table_index_t, unsigned> deme_sizes;
            bool found = false;
            for (auto& dip : pop.diploid_metadata)
                {
                    auto itr = deme_sizes.find(dip.deme);
                    if (itr == end(deme_sizes))
                        {
                            deme_sizes[dip.deme] = 0;
                        }
                    else
                        {
                            itr->second++;
                        }
                    if (dip.deme == deme)
                        {
                            found = true;
                        }
                }
            if (!found)
                {
                    throw std::invalid_argument(
                        "no alive individuals have the desired deme id");
                }
            if (static_cast<unsigned>(ndescendants) >= 2 * deme_sizes[deme])
                {
                    throw std::invalid_argument("ndescendants must be < 2*(deme size)");
                }
        }

    std::size_t new_mutation_key = std::numeric_limits<std::size_t>::max();
    auto candidates = generate_canidate_list(left, right, ndescendants, deme, pop);

    if (candidates.empty())
        {
            return new_mutation_key;
        }

    // Randomly choose a candidate proportionally to branch length
    std::vector<double> candidate_weights;
    for (auto c : candidates)
        {
            if (c.parent_time >= c.node_time)
                {
                    throw std::runtime_error(
                        "invalid parent/child times for candidate branch");
                }
            if (c.tree_span <= 0.0)
                {
                    throw std::runtime_error("invalid tree span of <= 0.0");
                }
            candidate_weights.push_back((c.node_time - c.parent_time) * c.tree_span);
        }
    auto discrete
        = gsl_ran_discrete_preproc(candidates.size(), candidate_weights.data());
    std::size_t candidate = gsl_ran_discrete(rng.get(), discrete);
    gsl_ran_discrete_free(discrete);

    auto candidate_data = std::move(candidates[candidate]);
    auto mutation_node_time = generate_mutation_time(
        rng, pop.tables->nodes, candidate_data.node, candidate_data.parent);
    if (mutation_node_time > pop.tables->nodes[candidate_data.node].time)
        {
            std::ostringstream o;
            o << "origin time of " << mutation_node_time
              << " for mutation is invalid: it must be <= mutation node time of "
              << pop.tables->nodes[candidate_data.node].time;
            throw std::runtime_error(o.str());
        }
    if (candidate_data.parent != fwdpp::ts::NULL_INDEX)
        {
            if (mutation_node_time <= pop.tables->nodes[candidate_data.parent].time)
                {
                    std::ostringstream o;
                    o << "origin time of " << mutation_node_time
                      << " for mutation is invalid: it must be < parent node time of "
                      << pop.tables->nodes[candidate_data.parent].time;
                    throw std::runtime_error(o.str());
                }
        }

    // NOTE: can we do all of what we need using existing
    // mutate/recombine functions?  Gotta check the back-end.

    // The new variant gets a derived state of 1
    auto empty = fwdpp::empty_mutation_queue();
    new_mutation_key = fwdpy11::infsites_Mutation(
        empty, pop.mutations, pop.mut_lookup, false, mutation_node_time,
        // The lambdas will simply transfer input parameters to the new object.
        [&rng, &candidate_data]() {
            return fwdpy11_core::internal::gsl_ran_flat(rng, candidate_data.left,
                                                        candidate_data.right);
        },
        [&data]() { return data.effect_size; },
        [&data](double) { return data.dominance; }, [&data]() { return data.esizes; },
        [&data]() { return data.heffects; }, data.label);

    pop.mcounts.push_back(ndescendants);
    if (pop.mutations.size() != pop.mcounts.size())
        {
            throw std::runtime_error("mutations.size() != mcounts.size()");
        }

    // Update the tables
    // 0 = ancestral_state
    auto site_id = pop.tables->emplace_back_site(pop.mutations[new_mutation_key].pos,
                                                 static_cast<std::int8_t>(0));

    pop.tables->emplace_back_mutation(candidate_data.node, new_mutation_key, site_id,
                                      static_cast<std::int8_t>(1), false);
    fwdpp::ts::sort_mutation_table_and_rebuild_site_table(*pop.tables);

    // Now, we update individual haploid genomes
    std::vector<std::size_t> nodes_to_haploid_genome(
        pop.tables->nodes.size(), std::numeric_limits<std::size_t>::max());
    std::vector<std::pair<std::size_t, int>> nodes_to_individuals(
        pop.tables->nodes.size(),
        std::make_pair(std::numeric_limits<std::size_t>::max(), -1));

    for (std::size_t i = 0; i < static_cast<std::size_t>(pop.N); ++i)
        {
            auto node0 = pop.diploid_metadata[i].nodes[0];
            auto node1 = pop.diploid_metadata[i].nodes[1];
            nodes_to_haploid_genome[node0] = pop.diploids[i].first;
            nodes_to_individuals[node0] = std::make_pair(i, 0);
            nodes_to_haploid_genome[node1] = pop.diploids[i].second;
            nodes_to_individuals[node1] = std::make_pair(i, 1);
        }

    for (auto n : candidate_data.descendants)
        {
            if (nodes_to_haploid_genome[static_cast<std::size_t>(n)]
                == std::numeric_limits<std::size_t>::max())
                {
                    throw std::runtime_error("bad node to genome mapping");
                }
            if (deme >= 0)
                {
                    if (pop.tables->nodes[n].deme != deme)
                        {
                            throw std::runtime_error("unexpected descendant deme");
                        }
                }
            if (pop.haploid_genomes[nodes_to_haploid_genome[n]].n == 0)
                {
                    throw std::runtime_error(
                        "attempting to mutate an extinct haploid_genome");
                }
            fwdpp::haploid_genome new_genome(1);
            auto itr = std::upper_bound(
                begin(pop.haploid_genomes[nodes_to_haploid_genome[n]].smutations),
                end(pop.haploid_genomes[nodes_to_haploid_genome[n]].smutations),
                new_mutation_key, [&pop](const auto new_key, const auto key) {
                    return pop.mutations[new_key].pos < pop.mutations[key].pos;
                });
            std::copy(begin(pop.haploid_genomes[nodes_to_haploid_genome[n]].smutations),
                      itr, std::back_inserter(new_genome.smutations));
            new_genome.smutations.push_back(new_mutation_key);
            std::copy(itr,
                      end(pop.haploid_genomes[nodes_to_haploid_genome[n]].smutations),
                      std::back_inserter(new_genome.smutations));

            // decrement the original genome count
            pop.haploid_genomes[nodes_to_haploid_genome[n]].n--;

            auto ind = nodes_to_individuals[n];
            if (ind.first == std::numeric_limits<std::size_t>::max())
                {
                    throw std::runtime_error("bad node to individual mapping");
                }
            // Add the new genome to the pop
            if (new_genome.smutations.empty())
                {
                    throw std::runtime_error("new genome is empty");
                }
            pop.haploid_genomes.emplace_back(std::move(new_genome));
            if (ind.second == 0)
                {
                    pop.diploids[ind.first].first = pop.haploid_genomes.size() - 1;
                }
            else if (ind.second == 1)
                {
                    pop.diploids[ind.first].second = pop.haploid_genomes.size() - 1;
                }
            else
                {
                    throw std::runtime_error("bad reference to diploid genome");
                }
        }

    for (auto& dip : pop.diploids)
        {
            if (pop.haploid_genomes[dip.first].n == 0)
                {
                    throw std::runtime_error("diploid refers to extinct genome");
                }
            if (pop.haploid_genomes[dip.second].n == 0)
                {
                    throw std::runtime_error("diploid refers to extinct genome");
                }
        }
    return new_mutation_key;
}
