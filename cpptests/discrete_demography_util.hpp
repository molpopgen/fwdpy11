#pragma once

#include <cstdint>
#include <fwdpy11/discrete_demography/MassMigration.hpp>

fwdpy11::discrete_demography::MassMigration move_individuals(std::uint32_t when,
                                                             int source, int destination,
                                                             double fraction,
                                                             bool resets_growth_rate);

fwdpy11::discrete_demography::MassMigration copy_individuals(std::uint32_t when,
                                                             int source, int destination,
                                                             double fraction,
                                                             bool resets_growth_rate);
