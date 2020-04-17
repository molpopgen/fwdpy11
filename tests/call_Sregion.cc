#include <pybind11/pybind11.h>
#include <fwdpy11/regions/Sregion.hpp>

fwdpy11::Mutation
call_Sregion(const fwdpy11::Sregion& s, unsigned seed)
{
    fwdpy11::GSLrng_t rng(seed);
    std::vector<fwdpy11::Mutation> mutations;
    std::unordered_multimap<double, std::uint32_t> lookup_table;
    fwdpp::flagged_mutation_queue recycling_bin = fwdpp::empty_mutation_queue();

    auto idx = s(recycling_bin, mutations, lookup_table, 0u, rng);
    return mutations[idx];
}

PYBIND11_MODULE(call_Sregion, m)
{
    m.def("call", &call_Sregion);
}
