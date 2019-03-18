// clang-format off
<% 
setup_pybind11(cfg) 
import common_mako
common_mako.setup_mako(cfg)
%>
// clang-format on

#include <pybind11/pybind11.h>
#include <fwdpy11/sim_functions.hpp>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/types/MlocusPop.hpp>

    namespace py = pybind11;

void
update_mutations(fwdpy11::Population& pop, const unsigned generation,
                 const unsigned twoN, const bool remove_selected_fixations)
// Expose some of fwdpy11's internal workings for unit testing
{
    fwdpy11::update_mutations(pop.mutations, pop.fixations, pop.fixation_times,
                              pop.mut_lookup, pop.mcounts, generation, twoN,
                              remove_selected_fixations);
}

PYBIND11_MODULE(fixation_properties, m)
{
    m.def("update_mutations", &update_mutations);
}
