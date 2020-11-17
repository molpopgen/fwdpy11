#include "discrete_demography_util.hpp"
fwdpy11::discrete_demography::MassMigration
move_individuals(std::uint32_t when, int source, int destination, double fraction,
                 bool resets_growth_rate)
{
    return fwdpy11::discrete_demography::MassMigration(
        when, source, destination, 0, -1, fraction, true, false, resets_growth_rate);
}

fwdpy11::discrete_demography::MassMigration
copy_individuals(std::uint32_t when, int source, int destination, double fraction,
                 bool resets_growth_rate)
{
    return fwdpy11::discrete_demography::MassMigration(
        when, source, destination, 0, -1, fraction, false, false, resets_growth_rate);
}
