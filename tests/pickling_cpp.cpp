// clang-format off
<% 
setup_pybind11(cfg) 
#import fwdpy11 so we can find its C++ headers
import fwdpy11 as fp11 
#add fwdpy11 header locations to the include path
cfg['include_dirs'].extend([ fp11.get_includes(), fp11.get_fwdpp_includes()])
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
#include <fwdpy11/types.hpp>

namespace py = pybind11;

//Example of pickling a specific C++ type
py::bytes
pickle_mutation(const KTfwd::popgenmut& p)
{
    py::object m = py::cast<decltype(p)>(p);
    return py::module::import("pickle").attr("dumps")(m);
}

//General pickler for any Python type.
//Also shows how to save pickle.dumps
//as a callable on the C++ side
py::bytes
general_pickler(py::object p)
{
    auto f = py::module::import("pickle").attr("dumps");
    return f(p);
}

PYBIND11_MODULE(pickling_cpp, m)
{
    m.def("pickle_mutation", &pickle_mutation);
    m.def("general_pickler", &general_pickler);
}
