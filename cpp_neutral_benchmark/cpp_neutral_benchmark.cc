#include "core/demes/forward_graph.hpp"
#include <iostream>
#include <vector>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/types/Mutation.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <fwdpy11/genetic_values/dgvalue_pointer_vector.hpp>
#include <fwdpy11/genetic_values/fwdpp_wrappers/fwdpp_genetic_value.hpp>
#include <fwdpy11/regions/RecombinationRegions.hpp>
#include <fwdpy11/regions/MutationRegions.hpp>
#include <fwdpy11/evolvets/SampleRecorder.hpp>
#include <fwdpy11/evolvets/recorders.hpp>
#include <fwdpy11/rng.hpp>
#include <core/evolve_discrete_demes/evolvets.hpp>
#include <core/genetic_maps/regions.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

// The next several structs are copied from
// the Python library C++ code. They are needed
// to define how Multiplicative works.

struct single_deme_multiplicative_het
{
    inline void
    operator()(double& d, const fwdpy11::Mutation& m) const
    {
        d *= (1. + m.s * m.h);
    }
};

struct multi_deme_multiplicative_het
{
    inline void
    operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
    {
        d *= (1. + m.esizes[deme] * m.heffects[deme]);
    }
};

struct single_deme_multiplicative_hom
{
    double scaling;
    single_deme_multiplicative_hom(double s) : scaling(s)
    {
    }

    inline void
    operator()(double& d, const fwdpy11::Mutation& m) const
    {
        d *= (1. + scaling * m.s);
    }
};

struct multi_deme_multiplicative_hom
{
    double scaling;
    multi_deme_multiplicative_hom(double s) : scaling(s)
    {
    }

    inline void
    operator()(const std::size_t deme, double& d, const fwdpy11::Mutation& m) const
    {
        d *= (1. + scaling * m.esizes[deme]);
    }
};

struct final_multiplicative_trait
{
    inline double
    operator()(double d) const
    {
        return d - 1.0;
    }
};

struct final_multiplicative_fitness
{
    inline double
    operator()(double d) const
    {
        return std::max(0.0, d);
    }
};

using DiploidMultiplicative = fwdpy11::stateless_site_dependent_genetic_value_wrapper<
    single_deme_multiplicative_het, single_deme_multiplicative_hom,
    multi_deme_multiplicative_het, multi_deme_multiplicative_hom, 1>;

struct RecordNothing
// copied from internal code
{
    inline void
    operator()(const fwdpy11::Population&) const
    {
    }
};

bool
no_stopping(const fwdpy11::Population&, const bool)
// copied from internal code
{
    return false;
}

struct command_line_options
{
    unsigned N;
    unsigned nsteps;
    unsigned simplification_interval;
    double xovers;
    unsigned seed;
    bool index_edge_table_during_sim;

    command_line_options();
};

command_line_options::command_line_options()
    : N{1000}, nsteps{1000}, simplification_interval{100}, xovers{0.}, seed{42},
      index_edge_table_during_sim{true}
{
}

po::options_description
generate_main_options(command_line_options& o)
{
    po::options_description options("Simulation options");
    options.add_options()("help", "Display help");
    options.add_options()("N", po::value<decltype(command_line_options::N)>(&o.N),
                          "Diploid population size. Default = 1000.");

    options.add_options()("nsteps",
                          po::value<decltype(command_line_options::nsteps)>(&o.nsteps),
                          "Number of time steps to evolve. Default = 1000.");
    options.add_options()(
        "simplify",
        po::value<decltype(command_line_options::simplification_interval)>(
            &o.simplification_interval),
        "Time steps between simplifications. Default = 100.");

    options.add_options()("index_edges", po::bool_switch(&o.index_edge_table_during_sim),
                          "If true, index edge tables during simulation (after each "
                          "simplification). Default = false");
    options.add_options()(
        "xovers", po::value<decltype(command_line_options::xovers)>(&o.xovers),
        "Mean number of crossovers (per parent, per mating).  Default=0.");
    options.add_options()("seed",
                          po::value<decltype(command_line_options::seed)>(&o.seed),
                          "Random number seed.  Default = 42.");

    return options;
}

void
simulate(const command_line_options& options)
{
    fwdpy11::DiploidPopulation pop(options.N, 1e7);
    fwdpy11::GSLrng_t rng(options.seed);
    fwdpy11_core::PoissonInterval recombination_region(0., 1e7, options.xovers, true);

    std::vector<std::unique_ptr<fwdpy11::PoissonCrossoverGenerator>> callbacks;
    callbacks.emplace_back(recombination_region.ll_clone());
    fwdpy11::GeneralizedGeneticMap genetic_map(std::move(callbacks), {});

    fwdpy11::MutationRegions mmodel({}, {});

    std::ostringstream o;
    o << "time_units: generations\n";
    o << "demes:\n"
      << " - name: A\n"
      << "   epochs:\n"
      << "    - start_size: " << options.N << '\n';

    fwdpy11_core::ForwardDemesGraph forward_demes_graph(o.str(), options.nsteps);

    DiploidMultiplicative fitness(1, 2., final_multiplicative_fitness(), nullptr,
                                  nullptr);
    fwdpy11::dgvalue_pointer_vector_ gvalue_pointers(fitness);

    fwdpy11::no_ancient_samples no_ancient_samples{};
    RecordNothing record_nothing;
    fwdpy11::SampleRecorder sample_recorder;

    evolve_with_tree_sequences_options tsoptions;

    tsoptions.preserve_selected_fixations = false;
    tsoptions.suppress_edge_table_indexing = !options.index_edge_table_during_sim;
    tsoptions.record_gvalue_matrix = false;
    tsoptions.track_mutation_counts_during_sim = false;
    tsoptions.remove_extinct_mutations_at_finish = true;
    tsoptions.reset_treeseqs_to_alive_nodes_after_simplification = false;
    tsoptions.preserve_first_generation = false;
    std::function<bool(const fwdpy11::DiploidPopulation&, const bool)> stopping_criteron
        = [](const fwdpy11::DiploidPopulation& pop, const bool simplified) {
              return no_stopping(pop, simplified);
          };

    evolve_with_tree_sequences(
        rng, pop, sample_recorder, options.simplification_interval, forward_demes_graph,
        options.nsteps, 0., 0., mmodel, genetic_map, gvalue_pointers,
        [](const fwdpy11::DiploidPopulation& /*pop*/, fwdpy11::SampleRecorder& /*sr*/) {
        },
        stopping_criteron, record_nothing, tsoptions);
    std::cout << pop.generation << ' ' << pop.N << ' ' << pop.tables->edges.size() << ' '
              << pop.tables->nodes.size() << '\n';
}

int
main(int argc, char** argv)
{
    command_line_options options;
    auto cli = generate_main_options(options);
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cli), vm);
    po::notify(vm);

    if (vm.count("help"))
        {
            std::cout << cli << '\n';
            std::exit(1);
        }

    simulate(options);
}
