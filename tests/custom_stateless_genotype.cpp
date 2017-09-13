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

// This is the only required header
#include <fwdpy11/fitness/single_locus_stateless_fitness.hpp>

    // This is the "1+t" part.
    // The macro defines a function
    // taking fitness (w) as an argument
    // and a mutation, m, as a constant
    // argument.  Here, we update w
    // according to a multiplicative scheme.
    // The macro argument is the name of the function
    // to generate.
    STATELESS_GENOTYPE_POLICY(Aa)
{
    w *= (1. + m.h);
}
END_STRUCT()

// Likewise, this is the 1+s part for an "aa"
//(homozygote) genotype.
STATELESS_GENOTYPE_POLICY(aa) { w *= (1. + m.s); }
END_STRUCT()

// Standard pybind11 stuff goes here
PYBIND11_MODULE(custom_stateless_genotype, m)
{
    // Call this macro so that your custom
    // class is recognizes are part of the
    // expected Python class hierarchy
    FWDPY11_SINGLE_LOCUS_FITNESS()

    // FWDPY11_SINGLE_LOCUS_STATELESS_FITNESS()
    // This macro creates a Python function
    // returning an instance of
    // fwdpy11.fitness.SlocusCustomStatelessGeneticValue
    // The macro arguments are:
    // 1. The name of your function for treating aa genotypes
    // 2. The name of your function for treating Aa genotypes
    // 3. A starting value for fitness/genetic value, which will be one for a
    // multiplicative model.
    // 4. The name of the Python function
    // 5. The name of the pybind11::module object
    CREATE_STATELESS_SLOCUS_GENOTYPE_OBJECT(aa, Aa, 1, "GeneralW", m)
}
