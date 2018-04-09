// clang-format off
<% 
setup_pybind11(cfg) 
#import fwdpy11 so we can find its C++ headers
import fwdpy11 as fp11 
#add fwdpy11 header locations to the include path
cfg['include_dirs'].extend([ fp11.get_includes(), fp11.get_fwdpp_includes() ])
#On OS X using clang, there is more work to do.  Using gcc on OS X
#gets rid of these requirements. The specifics sadly depend on how
#you initially built fwdpy11, and what is below assumes you used
#the provided setup.py + OS X + clang:
#cfg['compiler_args'].extend(['-stdlib=libc++','-mmacosx-version-min=10.7'])
#cfg['linker_args']=['-stdlib=libc++','-mmacosx-version-min=10.7']
#An alternative to the above is to add the first line to CPPFLAGS
#and the second to LDFLAGS when compiling a plugin on OS X using clang.
%>
// clang-format on

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
