#include <algorithm>
#include <fwdpp/ts/recycling.hpp>
#include <fwdpp/ts/mutate_tables.hpp>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/policies/mutation.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace
{
    inline std::size_t
    generate_neutral_variants(fwdpp::flagged_mutation_queue& recycling_bin,
                              fwdpy11::Population& pop,
                              const fwdpy11::GSLrng_t& rng, const double left,
                              const double right,
                              const fwdpp::uint_t generation)
    {
        const auto uniform = [left, right, &rng]() {
            return gsl_ran_flat(rng.get(), left, right);
        };
        const auto return_zero = []() { return 0.0; };
        return fwdpy11::infsites_Mutation(recycling_bin, pop.mutations,
                                          pop.mut_lookup, generation, uniform,
                                          return_zero, return_zero, 0);
    }
} // namespace

void
init_infinite_sites(py::module& m)
{
    m.def("infinite_sites", [](const fwdpy11::GSLrng_t& rng,
                               fwdpy11::Population& pop, const double mu) {
        fwdpp::flagged_mutation_queue recycling_bin
            = fwdpp::ts::make_mut_queue(pop.mcounts,
                                        pop.mcounts_from_preserved_nodes);

        std::vector<fwdpp::ts::TS_NODE_INT> samples;
        samples.reserve(pop.N);
        // Assume the pop is simplified
        auto nb = pop.tables.node_table.cbegin();
        auto past_sample = std::find_if(
            nb, pop.tables.node_table.cend(), [&pop](const fwdpp::ts::node& n) {
                return n.time != static_cast<double>(pop.generation);
            });
        for (auto i = nb; i < past_sample; ++i)
            {
                samples.push_back(std::distance(nb, i));
            }

        const auto apply_mutations =
            [&recycling_bin, &rng, &pop](const double left, const double right,
                                         const fwdpp::uint_t generation) {
                return generate_neutral_variants(recycling_bin, pop, rng, left,
                                                 right, generation);
            };
        auto nmuts = fwdpp::ts::mutate_tables(rng, apply_mutations, pop.tables,
                                              samples, mu);
        std::sort(
            pop.tables.mutation_table.begin(), pop.tables.mutation_table.end(),
            [&pop](const fwdpp::ts::mutation_record& a,
                   const fwdpp::ts::mutation_record& b) {
                return pop.mutations[a.key].pos < pop.mutations[b.key].pos;
            });
        fwdpp::ts::count_mutations(pop.tables, pop.mutations, samples,
                                   pop.mcounts,
                                   pop.mcounts_from_preserved_nodes);
        return nmuts;
    },py::arg("rng"),py::arg("pop"),py::arg("mu"));
}
