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

//This is the only required header
#include <fwdpy11/fitness/single_locus_stateless_fitness.hpp>

//This macro defines the function signature for you,
//but hides the gory details. All you need to know
//is that you have a diploid called "dip", a vector
//of gametes called "gametes", and a vector of mutations
//called "mutations".  These are standard C++ vectors
//and mimic the Python types in fwdpy11 exactly.
//If you want the gory details, see the header 
//included above.
//The macro argument is the name of your function.
STATELESS_SLOCUS_FUNCTION(additive)
{
    double s = 0.0;
    for (auto&& i : gametes[dip.first].smutations)
        s += mutations[i].s;
    for (auto&& i : gametes[dip.second].smutations)
        s += mutations[i].s;
    return std::max(0.0, 1 + s);
}
END_STRUCT()

//Standard pybind11 stuff goes here
PYBIND11_MODULE(custom_additive, m)
{
    //Call this macro so that your custom
    //class is recognizes are part of the 
    //expected Python class hierarchy
    FWDPY11_SINGLE_LOCUS_FITNESS()
    //FWDPY11_SINGLE_LOCUS_STATELESS_FITNESS()
    //This macro creates a Python function
	//that will return an instance of 
	//fwdpy11.fitness.SlocusCustomStatelessGeneticValue.
    //The macro arguments are:
    //1. The name of your function.
    //2. The name of the Python function
    //3. The name of the pybind11::module object
	CREATE_STATELESS_SLOCUS_OBJECT(additive,"additive",m);
}
