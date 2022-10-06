#include <cmath>
#include <algorithm>
#include <fwdpp/ts/recycling.hpp>
#include <fwdpp/ts/mutate_tables.hpp>
#include <fwdpp/ts/count_mutations.hpp>
#include <fwdpp/ts/mutation_tools.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/policies/mutation.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace
{
    inline fwdpp::ts::new_variant_record
    generate_neutral_variants(fwdpp::flagged_mutation_queue& recycling_bin,
                              fwdpy11::Population& pop, const fwdpy11::GSLrng_t& rng,
                              const double left, const double right,
                              const fwdpy11::mutation_origin_time generation)
    {
        const auto uniform
            = [left, right, &rng]() { return gsl_ran_flat(rng.get(), left, right); };
        const auto return_zero = []() { return 0.0; };
        const auto return_zero_h = [](const double) { return 0.0; };
        auto key = fwdpy11::infsites_Mutation(recycling_bin, pop.mutations,
                                              pop.mut_lookup, true, generation, uniform,
                                              return_zero, return_zero_h, 0);
        return fwdpp::ts::new_variant_record(
            pop.mutations[key].pos, fwdpp::ts::default_ancestral_state, key,
            pop.mutations[key].neutral, fwdpp::ts::default_derived_state);
    }
} // namespace

void
init_infinite_sites(py::module& m)
{
    m.def(
        "_infinite_sites",
        [](const fwdpy11::GSLrng_t& rng, fwdpy11::DiploidPopulation& pop,
           const double mu) -> unsigned {
            if (pop.is_simulating)
                {
                    throw std::runtime_error(
                        "infinite_sites cannot be called during a simulation");
                }
            if (mu <= 0.0)
                {
                    return 0u;
                }
            fwdpp::flagged_mutation_queue recycling_bin = fwdpp::ts::make_mut_queue(
                pop.mcounts, pop.mcounts_from_preserved_nodes);
            const auto apply_mutations
                = [&recycling_bin, &rng, &pop](const double left, const double right,
                                               const double origin_time) {
                      if (std::isfinite(origin_time))
                          {
                              return generate_neutral_variants(
                                  recycling_bin, pop, rng, left, right,
                                  static_cast<fwdpy11::mutation_origin_time>(
                                      std::ceil(origin_time)));
                          }
                      return generate_neutral_variants(
                          recycling_bin, pop, rng, left, right,
                          std::numeric_limits<fwdpy11::mutation_origin_time>::min());
                  };
            pop.fill_alive_nodes();
            pop.fill_preserved_nodes();
            auto nmuts = fwdpp::ts::mutate_tables(rng, apply_mutations, *pop.tables,
                                                  pop.alive_nodes, mu);
            if (nmuts == 0)
                {
                    return nmuts;
                }
            fwdpp::ts::count_mutations(*pop.tables, pop.mutations, pop.alive_nodes,
                                       pop.preserved_sample_nodes, pop.mcounts,
                                       pop.mcounts_from_preserved_nodes);
            pop.alive_nodes.clear();
            pop.preserved_sample_nodes.clear();
            return nmuts;
        },
        py::arg("rng"), py::arg("pop"), py::arg("mu"));
}
