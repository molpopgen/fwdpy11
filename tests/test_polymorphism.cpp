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

#include <pybind11/pybind11.h>
#include <fwdpy11/fitness/fitness.hpp>

    namespace py = pybind11;

PYBIND11_MODULE(test_polymorphism, m)
{
    m.def("test_callback_names",
          [](const std::vector<std::shared_ptr<fwdpy11::single_locus_fitness>>&
                 f) {
              std::vector<std::string> rv;
              for (auto&& i : f)
                  {
                      rv.emplace_back(i->callback_name());
                  }
              return rv;
          });
}
