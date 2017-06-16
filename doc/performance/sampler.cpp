// clang-format off
<% 
setup_pybind11(cfg) 
#import fwdpy11 so we can find its C++ headers
import fwdpy11 as fp11 
#add fwdpy11 header locations to the include path
cfg['include_dirs'] = [ fp11.get_includes(), fp11.get_fwdpp_includes() ] 
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
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/samplers.hpp>
#include <fwdpy11/types.hpp>
#include <cstdint>
#include <type_traits>

struct recorder_data
// This is our struct to hold
// our data
{
    std::uint32_t generation;
    double wbar, s1, s2;
};

struct recorder
// This is the C++ class representing our recorder.
// In terms of mechanics, it works identically to our
// Python implementation.
{
    mutable std::vector<recorder_data> data;
    mutable recorder_data d;

    recorder(const std::uint32_t simlen) : data(std::vector<recorder_data>())
    // Our constructor initializes our vector and reserves memory.
    {
        data.reserve(simlen);
    }

    void
    operator()(const fwdpy11::singlepop_t& pop) const
    // This is the __call__ body
    {
        double w = 0., s1 = 0., s2 = 0.;
        for (auto&& dip : pop.diploids)
            // Note that the syntax for this
            // for look is very similar to how
            // we would write it in Python!
            {
                w += dip.w;
                for (auto&& key : pop.gametes[dip.first].smutations)
                    {
                        s1 += pop.mutations[key].s;
                    }
                for (auto&& key : pop.gametes[dip.second].smutations)
                    {
                        s2 += pop.mutations[key].s;
                    }
            }
        d.generation = pop.generation;
        d.wbar = w / static_cast<double>(pop.N);
        d.s1 = s1 / static_cast<double>(pop.N);
        d.s2 = s2 / static_cast<double>(pop.N);
        data.push_back(d);
    }
};

// Now, we get to pybind11-specific stuff
namespace py = pybind11;

// This tells pybind11 to disable vector <-> list
// conversion for vectors of recorder data. This
// keeps the access "C++ only", as far as we are
// concerned.
PYBIND11_MAKE_OPAQUE(std::vector<recorder_data>);

// Now, we define a Python module called "sampler":
PYBIND11_PLUGIN(sampler)
{
    pybind11::module m("sampler", "Example of temporal sampler in C++");

    // Register recorder_data's members as a Numpy dtype
    PYBIND11_NUMPY_DTYPE(recorder_data, generation, wbar, s1, s2);

    // Allow our vector of recorder_data to be accessible in Python.
    // The py::buffer_protocol is a signal to wrap this type to be
    // compatible with Python's buffer protocol.
    py::bind_vector<std::vector<recorder_data>>(m, "VecRecorderData",
                                                py::buffer_protocol());

    // Finally, we make a Python class based on our C++ class.
    // This is pretty easy.
    py::class_<recorder>(m, "cppRecorder")
        // Create an __init__ function:
        .def(py::init<std::uint32_t>())
        // Give access to our buffer:
        .def_readonly("data", &recorder::data)
        // Define our __call__.  We must use a C++11
        // lambda here:
        .def("__call__",
             [](const recorder& r, const fwdpy11::singlepop_t& pop) -> void {
                 r(pop);
             });

    return m.ptr();
}
