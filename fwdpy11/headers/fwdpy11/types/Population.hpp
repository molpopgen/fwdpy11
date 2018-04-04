#ifndef FWDPY11_POPULATION_HPP__
#define FWDPY11_POPULATION_HPP__

#include "PyPopulation.hpp"
#include <unordered_set>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/sugar/popgenmut.hpp>

namespace fwdpy11
{
	// Population abstract base class based on
	// fwdpy11's Mutation type.
    using Population
        = PyPopulation<fwdpp::popgenmut, std::vector<fwdpp::popgenmut>,
                       std::vector<fwdpp::gamete>, std::vector<fwdpp::popgenmut>,
                       std::vector<fwdpp::uint_t>,
                       std::unordered_set<double, std::hash<double>,
                                          fwdpp::equal_eps>>;
}

#endif
