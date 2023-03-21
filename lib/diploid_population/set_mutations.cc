#include <fwdpy11/types/DiploidPopulation.hpp>

void
set_mutations(const std::vector<fwdpy11::Mutation> &mutations,
              const std::vector<std::int32_t> &mutation_nodes,
              const std::vector<fwdpy11::mutation_origin_time> &origin_times,
              fwdpy11::DiploidPopulation &pop)
{
    if (pop.is_simulating)
        {
            throw std::runtime_error("cannot set mutations for a simulating population");
        }
    if (!pop.mutations.empty())
        {
            throw std::runtime_error("population has existing mutations");
        }
    pop.mut_lookup.clear();
    pop.tables->mutations.clear();
    pop.tables->sites.clear();

    // Here, we cheat a bit
    pop.haploid_genomes.clear();
    for (auto &dip : pop.diploids)
        {
            pop.haploid_genomes.emplace_back(1);
            dip.first = pop.haploid_genomes.size() - 1;
            pop.haploid_genomes.emplace_back(1);
            dip.second = pop.haploid_genomes.size() - 1;
        }

    // Paranoia
    for (auto &g : pop.haploid_genomes)
        {
            if (g.n != 1)
                {
                    throw std::runtime_error(
                        "all haploid_genomes must have a count of 1");
                }
        }

    for (std::size_t i = 0; i < mutations.size(); ++i)
        {
            if (pop.mut_lookup.find(mutations[i].pos) != end(pop.mut_lookup))
                {
                    throw std::invalid_argument("duplicate mutation positions");
                }
            pop.mutations.emplace_back(
                fwdpy11::Mutation(mutations[i].neutral, mutations[i].pos, mutations[i].s,
                                  mutations[i].h, -origin_times[i], mutations[i].esizes,
                                  mutations[i].heffects, mutations[i].xtra));
            pop.mut_lookup.emplace(mutations[i].pos, i);
            pop.tables->emplace_back_site(mutations[i].pos, std::int8_t{0});
            pop.tables->emplace_back_mutation(
                mutation_nodes[i], pop.mutations.size() - 1,
                pop.tables->sites.size() - 1, std::int8_t{1}, mutations[i].neutral);
        }
    pop.mcounts.resize(pop.mutations.size());
    std::fill(begin(pop.mcounts), end(pop.mcounts), 0);
    fwdpp::ts::sort_mutation_table_and_rebuild_site_table(*pop.tables);
    std::vector<fwdpp::ts::table_index_t> samples;
    std::unordered_map<fwdpp::ts::table_index_t, std::size_t> node_to_genome;
    for (std::size_t i = 0; i < pop.diploid_metadata.size(); ++i)
        {
            for (auto n : pop.diploid_metadata[i].nodes)
                {
                    if (node_to_genome.find(n) != end(node_to_genome))
                        {
                            throw std::runtime_error("node id use multiple times");
                        }
                    samples.push_back(n);
                }
            node_to_genome.emplace(pop.diploid_metadata[i].nodes[0],
                                   pop.diploids[i].first);
            node_to_genome.emplace(pop.diploid_metadata[i].nodes[1], pop.diploids[i].second);
        }
    if (samples.empty())
        {
            throw std::runtime_error(
                "samples list for adding tskit mutations to genomes is empty");
        }

    auto sv
        = fwdpp::ts::site_visitor<fwdpp::ts::std_table_collection>(*pop.tables, samples);
    auto site = sv();
    std::size_t nmutations_processed = 0;
    while (site != end(sv))
        {
            auto mutations = sv.get_mutations();
            if (mutations.second - mutations.first != 1)
                {
                    throw std::runtime_error("expected exactly one mutation per site");
                }
            for (auto m = mutations.first; m != mutations.second; ++m)
                {
                    fwdpp::ts::samples_iterator si(
                        sv.current_tree(), m->node,
                        fwdpp::ts::convert_sample_index_to_nodes(true));
                    if (m->key >= pop.mutations.size())
                        {
                            throw std::runtime_error("invalid mutation key");
                        }
                    auto parent_node = sv.current_tree().parents[m->node];
                    auto node_time = pop.tables->nodes[m->node].time;

                    if (pop.mutations[m->key].g > node_time)
                        {
                            std::ostringstream o;
                            o << "invalid mutation origin time";
                            throw std::runtime_error(o.str());
                        }
                    if (parent_node >= 0)
                        {
                            auto parent_time
                                = pop.tables->nodes[sv.current_tree().parents[m->node]]
                                      .time;
                            if (pop.mutations[m->key].g <= parent_time)
                                {
                                    std::ostringstream o;
                                    o << "invalid mutation origin time";
                                    throw std::runtime_error(o.str());
                                }
                        }
                    auto s = si();
                    while (s != fwdpp::ts::NULL_INDEX)
                        {
                            auto lookup = node_to_genome.find(s);
                            if (lookup == end(node_to_genome))
                                {
                                    throw std::runtime_error(
                                        "mutation could not be mapped to a "
                                        "sample node");
                                }
                            if (m->neutral)
                                {
                                    pop.haploid_genomes[lookup->second]
                                        .mutations.emplace_back(m->key);
                                }
                            else
                                {
                                    pop.haploid_genomes[lookup->second]
                                        .smutations.emplace_back(m->key);
                                }
                            pop.mcounts[m->key] += 1;
                            s = si();
                        }
                    nmutations_processed += 1;
                }
            site = sv();
        }
    if (nmutations_processed != pop.mutations.size()
        || nmutations_processed != pop.tables->mutations.size())
        {
            throw std::runtime_error(
                "failed to process the expected number of mutations");
        }
}
