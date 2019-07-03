#ifndef FWDPY11_POPULATION_HPP__
#define FWDPY11_POPULATION_HPP__

#include "PyPopulation.hpp"
#include "Mutation.hpp"
#include <unordered_map>
#include <fwdpp/fwd_functional.hpp>

namespace fwdpy11
{
	// Population abstract base class based on
	// fwdpy11's Mutation type.
    using Population
        = PyPopulation<Mutation, std::vector<Mutation>,
                       std::vector<fwdpp::haploid_genome>, std::vector<Mutation>,
                       std::vector<fwdpp::uint_t>,
                       std::unordered_multimap<double,fwdpp::uint_t>>;
}

#endif
