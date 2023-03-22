#include <algorithm>
#include <cstdint>
#include <fwdpp/simparams.hpp>
#include <fwdpp/ts/simplify_tables.hpp>
#include <fwdpp/ts/simplify_tables_output.hpp>
#include <fwdpp/ts/table_collection_functions.hpp>
#include <fwdpp/ts/recording/edge_buffer.hpp>
#include <fwdpp/ts/make_simplifier_state.hpp>
#include <fwdpp/ts/recycling.hpp>
#include <fwdpy11/gsl/gsl_error_handler_wrapper.hpp>

#include <core/demes/forward_graph.hpp>

#include "fwdpy11/discrete_demography/exceptions.hpp"
#include "fwdpy11/types/Diploid.hpp"
#include "util.hpp"
#include "diploid_pop_fitness.hpp"
#include "index_and_count_mutations.hpp"
#include "cleanup_metadata.hpp"
#include "track_mutation_counts.hpp"
#include "remove_extinct_mutations.hpp"
#include "track_ancestral_counts.hpp"
#include "remove_extinct_genomes.hpp"
#include "runtime_checks.hpp"
#include "evolve_generation_ts.hpp"
#include "simplify_tables.hpp"
#include "discrete_demography/simulation/multideme_fitness_bookmark.hpp"

#include <core/evolve_discrete_demes/evolvets.hpp>
#include <sstream>

namespace ddemog = fwdpy11_core::discrete_demography;

void
apply_treseq_resetting_of_ancient_samples(
    const fwdpy11::DiploidPopulation_temporal_sampler &recorder,
    fwdpy11::DiploidPopulation &pop)
{
    recorder(pop);
    if (!pop.ancient_sample_metadata.empty())
        {
            pop.ancient_sample_metadata.clear();
            pop.ancient_sample_genetic_value_matrix.clear();
        }
}

// NOTE: this should be part of fwdpp's edge table functions header!
void
clear_edge_table_indexes(fwdpp::ts::std_table_collection &tables)
{
    tables.input_left.clear();
    tables.output_right.clear();
}

template <typename SimplificationState>
void
simplification(
    bool preserve_selected_fixations, bool suppress_edge_table_indexing,
    bool reset_treeseqs_to_alive_nodes_after_simplification,
    const fwdpy11::DiploidPopulation_temporal_sampler &post_simplification_recorder,
    SimplificationState &simplifier_state,
    fwdpp::ts::simplify_tables_output &simplification_output,
    fwdpp::ts::edge_buffer &new_edge_buffer,
    std::vector<fwdpp::ts::table_index_t> &alive_at_last_simplification,
    fwdpy11::DiploidPopulation &pop)
{
    simplify_tables(pop, pop.mcounts_from_preserved_nodes, alive_at_last_simplification,
                    *pop.tables, simplifier_state, simplification_output,
                    new_edge_buffer, preserve_selected_fixations,
                    suppress_edge_table_indexing);
    if (pop.mcounts.size() != pop.mcounts_from_preserved_nodes.size())
        {
            throw std::runtime_error("evolvets: count vector size mismatch after "
                                     "simplification");
        }
    remap_metadata(pop.ancient_sample_metadata, simplification_output.idmap);
    remap_metadata(pop.diploid_metadata, simplification_output.idmap);
    pop.fill_alive_nodes();
    alive_at_last_simplification.assign(begin(pop.alive_nodes), end(pop.alive_nodes));

    if (reset_treeseqs_to_alive_nodes_after_simplification == true)
        {
            apply_treseq_resetting_of_ancient_samples(post_simplification_recorder, pop);
        }
}

void
final_population_cleanup(
    bool suppress_edge_table_indexing, bool preserve_selected_fixations,
    bool remove_extinct_mutations_at_finish, bool simulating_neutral_variants,
    bool reset_treeseqs_to_alive_nodes_after_simplification,
    std::uint32_t last_preserved_generation,
    const std::vector<std::uint32_t> & /*last_preserved_generation_counts*/,
    fwdpy11::DiploidPopulation &pop)
{
    index_and_count_mutations(suppress_edge_table_indexing, simulating_neutral_variants,
                              reset_treeseqs_to_alive_nodes_after_simplification,
                              last_preserved_generation == pop.generation, pop);
    check_mutation_table_consistency_with_count_vectors(pop, __FILE__, __LINE__);
    // The following block was commented out in GitHub PR 845
    // in order to fix GitHub issue 844
    // if (pop.generation == last_preserved_generation
    //     && !reset_treeseqs_to_alive_nodes_after_simplification
    //     && !simulating_neutral_variants)
    //     {
    //         std::transform(begin(pop.mcounts_from_preserved_nodes),
    //                        end(pop.mcounts_from_preserved_nodes),
    //                        begin(last_preserved_generation_counts),
    //                        begin(pop.mcounts_from_preserved_nodes),
    //                        std::minus<std::uint32_t>());
    //     }
    // check_mutation_table_consistency_with_count_vectors(pop, __FILE__, __LINE__);
    if (!preserve_selected_fixations)
        {
            auto itr = std::remove_if(
                pop.tables->mutations.begin(), pop.tables->mutations.end(),
                [&pop](const fwdpp::ts::mutation_record &mr) {
                    return pop.mcounts[mr.key] == 2 * pop.diploids.size()
                           && pop.mcounts_from_preserved_nodes[mr.key] == 0;
                });
            auto d = std::distance(itr, end(pop.tables->mutations));
            pop.tables->mutations.erase(itr, end(pop.tables->mutations));
            if (d)
                {
                    fwdpp::ts::rebuild_site_table(*pop.tables);
                }
            fwdpp::ts::remove_fixations_from_haploid_genomes(
                pop.haploid_genomes, pop.mutations, pop.mcounts,
                pop.mcounts_from_preserved_nodes, 2 * pop.diploids.size(),
                preserve_selected_fixations);
            // NOTE: this is hacky and should be better-handled upstream
            for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
                {
                    if (pop.mcounts[i] == 2 * pop.N
                        && pop.mcounts_from_preserved_nodes[i] == 0)
                        {
                            pop.mcounts[i] = 0;
                        }
                }
        }
    cleanup_metadata(*pop.tables, pop.generation, pop.ancient_sample_metadata);
    if (remove_extinct_mutations_at_finish)
        {
            remove_extinct_mutations(pop);
        }
    remove_extinct_genomes(pop);
    if (pop.mutations.size() != pop.mcounts.size()
        || pop.mutations.size() != pop.mcounts_from_preserved_nodes.size())
        {
            throw std::runtime_error("mutation container size does not equal mutation "
                                     "count container size(s)");
        }
    pop.is_simulating = false;
}

void
check_initial_deme_sizes(const std::vector<fwdpy11::DiploidMetadata> &metadata,
                         const fwdpy11_core::ForwardDemesGraph &demography)
{
    auto parental_deme_sizes = demography.parental_deme_sizes();
    auto num_extant_parental_demes
        = std::count_if(std::begin(parental_deme_sizes), std::end(parental_deme_sizes),
                        [](auto a) { return a > 0.0; });
    if (num_extant_parental_demes == 0)
        {
            throw fwdpy11::discrete_demography::DemographyError(
                "all parental deme sizes are zero");
        }
    std::vector<std::uint32_t> current_deme_sizes(demography.number_of_demes(), 0);
    for (const auto &md : metadata)
        {
            if (md.deme >= demography.number_of_demes() || md.deme < 0)
                {
                    throw fwdpy11::discrete_demography::DemographyError(
                        "individual has invalid deme");
                }
            current_deme_sizes[md.deme] += 1;
        }
    std::size_t i = 0;
    for (auto psize : parental_deme_sizes)
        {
            auto intsize = static_cast<std::uint32_t>(psize);
            if (intsize != current_deme_sizes[i])
                {
                    std::ostringstream o;
                    o << "initial size of deme " << i << " is " << current_deme_sizes[i]
                      << " but the required size is " << intsize;
                    throw fwdpy11::discrete_demography::DemographyError(o.str());
                }
            ++i;
        }
}

void
evolve_with_tree_sequences(
    const fwdpy11::GSLrng_t &rng, fwdpy11::DiploidPopulation &pop,
    fwdpy11::SampleRecorder &sr, const unsigned simplification_interval,
    fwdpy11_core::ForwardDemesGraph &demography, const std::uint32_t simlen,
    const double mu_neutral, const double mu_selected,
    const fwdpy11::MutationRegions &mmodel, const fwdpy11::GeneticMap &rmodel,
    // NOTE: gvalue_pointers is a change in 0.6.0,
    // and the object holds non-const bare pointers
    // to objects owned by Python.
    fwdpy11::dgvalue_pointer_vector_ &gvalue_pointers,
    fwdpy11::DiploidPopulation_sample_recorder recorder,
    std::function<bool(const fwdpy11::DiploidPopulation &, const bool)>
        &stopping_criteron,
    const fwdpy11::DiploidPopulation_temporal_sampler &post_simplification_recorder,
    evolve_with_tree_sequences_options options)
{
    // FIXME: the pop's state must match what is expected by the ForwardDemesGraph!
    fwdpy11::gsl_scoped_convert_error_to_exception gsl_error_scope_guard;

    if (gvalue_pointers.genetic_values.empty())
        {
            throw std::invalid_argument("empty list of genetic values");
        }
    //validate the input params
    if (pop.tables->genome_length() == std::numeric_limits<double>::max())
        {
            throw std::invalid_argument(
                "Population is not initialized with tree sequence support");
        }
    if (!std::isfinite(mu_selected))
        {
            throw std::invalid_argument("selected mutation rate is not finite");
        }
    if (mu_selected < 0.0)
        {
            throw std::invalid_argument("selected mutation rate must be non-negative");
        }
    if (!std::isfinite(mu_neutral))
        {
            throw std::invalid_argument("neutral mutation rate is not finite");
        }
    if (mu_neutral < 0.0)
        {
            throw std::invalid_argument("neutral mutation rate must be non-negative");
        }
    if (mu_neutral + mu_selected > 0.0 && mmodel.weights.empty())
        {
            throw std::invalid_argument(
                "nonzero mutation rate incompatible with empty regions");
        }
    if (!simlen)
        {
            throw std::invalid_argument("simulation length must be > 0");
        }
    if (pop.tables->nodes.empty())
        {
            throw std::invalid_argument("node table is not initialized");
        }
    const bool simulating_neutral_variants = (mu_neutral > 0.0) ? true : false;
    if (simulating_neutral_variants)
        {
            if (options.track_mutation_counts_during_sim)
                {
                    if (simplification_interval != 1)
                        {
                            throw std::invalid_argument(
                                "when track_mutation_counts is True and simulating "
                                "neutral mutations, the simplification interval must be "
                                "1");
                        }
                }
        }

    for (auto &region : mmodel.regions)
        {
            if (!region->valid(0.0, pop.tables->genome_length()))
                {
                    std::ostringstream o;
                    o << "region contains invalid beg and/or end values: "
                      << region->beg() << ", " << region->end()
                      << ", genome_length = " << pop.tables->genome_length();
                    throw std::invalid_argument(o.str().c_str());
                }
        }

    if (pop.generation >= demography.model_end_time())
        {
            throw fwdpy11::discrete_demography::DemographyError(
                "current time of population is past the end of the model");
        }

    demography.initialize_model(pop.generation);

    if (demography.number_of_demes() <= 0)
        {
            throw std::runtime_error("maxdemes must be > 0");
        }
    check_initial_deme_sizes(pop.diploid_metadata, demography);

    // NOTE: should this be != instead of > ??
    if (gvalue_pointers.genetic_values.size()
        > static_cast<std::size_t>(demography.number_of_demes()))
        {
            throw std::invalid_argument(
                "list of genetic values is longer than maxdemes");
        }
    if (static_cast<std::size_t>(demography.number_of_demes()) > 1
        && gvalue_pointers.genetic_values.size() > 1
        && gvalue_pointers.genetic_values.size()
               < static_cast<std::size_t>(demography.number_of_demes()))
        {
            throw std::invalid_argument("too few genetic value objects");
        }

    double total_mutation_rate = mu_neutral + mu_selected;
    const auto bound_mmodel = [&rng, &mmodel, &pop, total_mutation_rate](
                                  fwdpp::flagged_mutation_queue &recycling_bin,
                                  std::vector<fwdpy11::Mutation> &mutations) {
        std::vector<fwdpp::uint_t> rv;
        unsigned nmuts = gsl_ran_poisson(rng.get(), total_mutation_rate);
        for (unsigned i = 0; i < nmuts; ++i)
            {
                std::size_t x = gsl_ran_discrete(rng.get(), mmodel.lookup.get());
                auto key = mmodel.regions[x]->operator()(
                    recycling_bin, mutations, pop.mut_lookup, pop.generation, rng);
                rv.push_back(key);
            }
        std::sort(begin(rv), end(rv),
                  [&mutations](const fwdpp::uint_t a, const fwdpp::uint_t b) {
                      return mutations[a].pos < mutations[b].pos;
                  });
        return rv;
    };

    const auto bound_rmodel = [&rng, &rmodel]() { return rmodel(rng); };

    auto genetics = fwdpp::make_genetic_parameters(gvalue_pointers.genetic_values,
                                                   std::move(bound_mmodel),
                                                   std::move(bound_rmodel));
    // NOTE: this could be a bit tricky!
    // A demes model gives ids according to order
    // in the graph.
    // We may need an API to allow a Python-side
    // mapping of deme name -> genetic value object.
    std::vector<std::size_t> deme_to_gvalue_map(demography.number_of_demes(), 0);
    if (genetics.gvalue.size() > 1)
        {
            std::iota(begin(deme_to_gvalue_map), end(deme_to_gvalue_map), 0);
        }
    // A stateful fitness model will need its data up-to-date,
    // so we must call update(...) prior to calculating fitness,
    // else bad stuff like segfaults could happen.
    for (auto &i : genetics.gvalue)
        {
            i->update(pop);
            i->gv2w->update(pop);
            i->noise_fxn->update(pop);
        }
    std::vector<fwdpy11::DiploidMetadata> offspring_metadata(pop.diploid_metadata);
    std::vector<fwdpy11::DiploidGenotype> offspring;
    std::vector<double> new_diploid_gvalues;

    calculate_diploid_fitness(rng, pop, genetics.gvalue, deme_to_gvalue_map,
                              offspring_metadata, new_diploid_gvalues,
                              options.record_gvalue_matrix);
    pop.genetic_value_matrix.swap(new_diploid_gvalues);
    pop.diploid_metadata.swap(offspring_metadata);

    ddemog::multideme_fitness_lookups<std::uint32_t> fitness_lookup{
        static_cast<std::int32_t>(demography.number_of_demes())};

    ddemog::multideme_fitness_bookmark fitness_bookmark;

    fitness_bookmark.update(demography.parental_deme_sizes(), pop.diploid_metadata);
    fitness_lookup.update(fitness_bookmark);

    // TODO:  Do we have sufficient test coverage through here?
    if (!pop.mutations.empty())
        {
            // It is possible that pop already has a tree sequence
            // containing neutral variants not in haploid_genome objects.
            // To correctly handle this case, we build the recycling bin
            // from any elements in pop.mutations not in the mutation table.
            if (!pop.tables->nodes.empty())
                {
                    std::vector<std::size_t> indexes(pop.mutations.size());
                    std::iota(begin(indexes), end(indexes), 0);
                    for (auto &m : pop.tables->mutations)
                        {
                            indexes[m.key] = std::numeric_limits<std::size_t>::max();
                        }
                    std::queue<std::size_t> can_recycle;
                    for (auto i : indexes)
                        {
                            if (i != std::numeric_limits<std::size_t>::max())
                                {
                                    can_recycle.push(i);
                                }
                        }
                    genetics.mutation_recycling_bin
                        = fwdpp::flagged_mutation_queue(std::move(can_recycle));
                }
            else // Assume there are no tree sequence data
                {
                    // Then we assume pop exists in an "already simulated"
                    // state and is properly-book-kept
                    genetics.mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                        pop.mcounts, pop.mcounts_from_preserved_nodes);
                }
        }

    fwdpp::ts::table_index_t next_index = pop.tables->nodes.size();
    bool simplified = false;
    auto simplifier_state
        = std::make_unique<decltype(fwdpp::ts::make_simplifier_state(*pop.tables))>(
            fwdpp::ts::make_simplifier_state(*pop.tables));
    auto new_edge_buffer
        = std::make_unique<fwdpp::ts::edge_buffer>(fwdpp::ts::edge_buffer{});
    bool stopping_criteron_met = false;
    std::pair<std::vector<fwdpp::ts::table_index_t>, std::vector<std::size_t>>
        simplification_rv;
    std::uint32_t last_preserved_generation = std::numeric_limits<std::uint32_t>::max();
    decltype(pop.mcounts) last_preserved_generation_counts;
    pop.fill_alive_nodes();
    if (options.preserve_first_generation)
        {
            if (pop.generation != 0)
                {
                    throw std::invalid_argument(
                        "cannot preserve first generation when pop.generation != 0");
                }
            if (pop.tables->edges.empty() == false)
                {
                    throw std::invalid_argument("cannot preserve first generation when "
                                                "the edge table is not empty");
                }
            std::vector<std::uint32_t> individuals;
            for (const auto &md : pop.diploid_metadata)
                {
                    individuals.push_back(md.label);
                }
            pop.record_ancient_samples(individuals);
            pop.update_ancient_sample_genetic_value_matrix(
                individuals, genetics.gvalue[0]->total_dim);
            track_ancestral_counts(individuals, &last_preserved_generation,
                                   last_preserved_generation_counts, pop);
        }

    std::vector<fwdpp::ts::table_index_t> alive_at_last_simplification(pop.alive_nodes);
    new_edge_buffer->reset(alive_at_last_simplification.size());

    clear_edge_table_indexes(*pop.tables);
    fwdpp::ts::simplify_tables_output simplification_output;
    pop.is_simulating = true;
    for (std::uint32_t gen = 0; pop.generation < demography.model_end_time() - 1
                                && gen < simlen && !stopping_criteron_met;
         ++gen)
        {
            if (demography.in_error_state())
                {
                    throw std::runtime_error("forward graph is in an error state");
                }
            ++pop.generation;
            evolve_generation_ts(rng, pop, genetics, demography, fitness_lookup,
                                 fitness_bookmark, // miglookup,
                                 pop.generation, *new_edge_buffer, offspring,
                                 offspring_metadata, next_index,
                                 options.allow_residual_selfing);
            // TODO: abstract out these steps into a "cleanup_pop" function
            // NOTE: by swapping the diploids here, it is not possible
            // for genetics.value to make use of parental genotype information.
            // Although mutations that went extinct this generation still
            // exist in pop, we are limited in our ability to pass down
            // "indirect" genetic effects for more than one generation
            // due to the possibility of recycling.  It is likely
            // that an explicity pedigree structure would be required
            // to make these models "nice".  See GitHub issue 372
            // for a bit more context.
            pop.diploids.swap(offspring);

            // NOTE: the two swaps of the metadata ensure
            // that the update loop below passes the correct
            // metadata on, and that we then have the
            // metadata in the expected places for
            // calculate_diploid_fitness
            pop.diploid_metadata.swap(offspring_metadata);
            // TODO: deal with random effects
            for (auto &i : genetics.gvalue)
                {
                    i->update(pop);
                    i->gv2w->update(pop);
                    i->noise_fxn->update(pop);
                }
            pop.diploid_metadata.swap(offspring_metadata);
            calculate_diploid_fitness(rng, pop, genetics.gvalue, deme_to_gvalue_map,
                                      offspring_metadata, new_diploid_gvalues,
                                      options.record_gvalue_matrix);
            pop.genetic_value_matrix.swap(new_diploid_gvalues);
            // TODO: abstract out these steps into a "cleanup_pop" function
            pop.diploid_metadata.swap(offspring_metadata);
            pop.N = static_cast<std::uint32_t>(pop.diploids.size());

            if (gen % simplification_interval == 0.0)
                {
                    simplification(
                        options.preserve_selected_fixations,
                        options.suppress_edge_table_indexing,
                        options.reset_treeseqs_to_alive_nodes_after_simplification,
                        post_simplification_recorder, *simplifier_state,
                        simplification_output, *new_edge_buffer,
                        alive_at_last_simplification, pop);
                    simplified = true;
                }
            else
                {
                    simplified = false;
                    clear_edge_table_indexes(*pop.tables);
                }
            if (pop.tables->num_nodes()
                >= std::numeric_limits<fwdpp::ts::table_index_t>::max() - 1)
                {
                    throw std::runtime_error("range error for node labels");
                }
            next_index = pop.tables->num_nodes();
            if (options.track_mutation_counts_during_sim)
                {
                    track_mutation_counts(pop, simplified,
                                          options.suppress_edge_table_indexing);
                }

            // The user may now analyze the pop'n and record ancient samples
            recorder(pop, sr);

            if (simplified)
                {
                    if (options.suppress_edge_table_indexing == false)
                        {
                            // Behavior change in 0.5.3: set all fixation counts to 0
                            // to flag for recycling if possible.
                            // NOTE: this may slow things down a touch?
                            if (options.preserve_selected_fixations == false)
                                {
                                    // b/c neutral mutations not in genomes!
                                    fwdpp::ts::remove_fixations_from_haploid_genomes(
                                        pop.haploid_genomes, pop.mutations, pop.mcounts,
                                        pop.mcounts_from_preserved_nodes,
                                        2 * pop.diploids.size(),
                                        options.preserve_selected_fixations);
                                }
                            for (auto &i : simplification_output.preserved_mutations)
                                {
                                    if (pop.mcounts[i] == 2 * pop.diploids.size()
                                        && pop.mcounts_from_preserved_nodes[i] == 0)
                                        {
                                            if (pop.mutations[i].neutral
                                                || !options.preserve_selected_fixations)
                                                {
                                                    // flag variant for recycling
                                                    pop.mcounts[i] = 0;
                                                    // flag item for removal from return value,
                                                    // as mutation is no longer considered "preserved"
                                                    i = std::numeric_limits<
                                                        std::size_t>::max();
                                                }
                                        }
                                }
                            simplification_output.preserved_mutations.erase(
                                std::remove(
                                    begin(simplification_output.preserved_mutations),
                                    end(simplification_output.preserved_mutations),
                                    std::numeric_limits<std::size_t>::max()),
                                end(simplification_output.preserved_mutations));
                            if (simulating_neutral_variants)
                                {
                                    genetics.mutation_recycling_bin
                                        = fwdpp::ts::make_mut_queue(
                                            simplification_output.preserved_mutations,
                                            pop.mutations.size());
                                }
                            else
                                {
                                    genetics.mutation_recycling_bin
                                        = fwdpp::ts::make_mut_queue(
                                            pop.mcounts,
                                            pop.mcounts_from_preserved_nodes);
                                }
                        }
                    else
                        {
                            genetics.mutation_recycling_bin = fwdpp::ts::make_mut_queue(
                                simplification_output.preserved_mutations,
                                pop.mutations.size());
                        }
                }

            demography.iterate_state();

            fitness_bookmark.update(demography.parental_deme_sizes(),
                                    pop.diploid_metadata);
            fitness_lookup.update(fitness_bookmark);
            ddemog::validate_parental_state(pop.generation, fitness_lookup, demography);

            // TODO: deal with the result of the recorder populating sr
            if (!sr.samples.empty())
                {
                    pop.record_ancient_samples(sr.samples);
                    pop.update_ancient_sample_genetic_value_matrix(
                        sr.samples, genetics.gvalue[0]->total_dim);
                    track_ancestral_counts(sr.samples, &last_preserved_generation,
                                           last_preserved_generation_counts, pop);
                    // Finally, clear the input
                    sr.samples.clear();
                }
            stopping_criteron_met = stopping_criteron(pop, simplified);
        }

    // NOTE: if pop.preserved_sample_nodes overlaps with samples,
    // then simplification throws an error. But, since it is annoying
    // for a user to have to remember not to do that, we filter the list
    // here
    pop.fill_alive_nodes();
    std::sort(begin(pop.alive_nodes), end(pop.alive_nodes));
    pop.fill_preserved_nodes();
    auto itr = std::remove_if(
        begin(pop.preserved_sample_nodes), end(pop.preserved_sample_nodes),
        [&pop](const fwdpp::ts::table_index_t l) {
            return std::binary_search(begin(pop.alive_nodes), end(pop.alive_nodes), l);
        });
    pop.preserved_sample_nodes.erase(itr, end(pop.preserved_sample_nodes));
    pop.alive_nodes.clear();
    pop.preserved_sample_nodes.clear();
    pop.ancient_sample_metadata.erase(
        std::remove_if(begin(pop.ancient_sample_metadata),
                       end(pop.ancient_sample_metadata),
                       [&pop](const auto &md) {
                           return pop.tables->nodes[md.nodes[0]].time == pop.generation;
                       }),
        end(pop.ancient_sample_metadata));

    if (!simplified)
        {
            simplification(options.preserve_selected_fixations,
                           options.suppress_edge_table_indexing,
                           options.reset_treeseqs_to_alive_nodes_after_simplification,
                           post_simplification_recorder, *simplifier_state,
                           simplification_output, *new_edge_buffer,
                           alive_at_last_simplification, pop);
            if (!options.preserve_selected_fixations)
                {
                    fwdpp::ts::remove_fixations_from_haploid_genomes(
                        pop.haploid_genomes, pop.mutations, pop.mcounts,
                        pop.mcounts_from_preserved_nodes, 2 * pop.diploids.size(),
                        options.preserve_selected_fixations);
                }
        }
    simplifier_state.reset(nullptr);
    new_edge_buffer.reset(nullptr);
    final_population_cleanup(
        options.suppress_edge_table_indexing, options.preserve_selected_fixations,
        options.remove_extinct_mutations_at_finish, simulating_neutral_variants,
        options.reset_treeseqs_to_alive_nodes_after_simplification,
        last_preserved_generation, last_preserved_generation_counts, pop);
    if (pop.tables->edges.size() != pop.tables->input_left.size()
        || pop.tables->edges.size() != pop.tables->output_right.size())
        {
            std::ostringstream o;
            o << "edge table is not indexed " << strip_unix_path(__FILE__) << ", line "
              << __LINE__;
            throw std::runtime_error(o.str());
        }
}
