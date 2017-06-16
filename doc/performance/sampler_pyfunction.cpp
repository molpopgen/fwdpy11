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

struct genetic_values_per_diploid
{
    std::uint32_t generation;
    std::uint32_t individual;
    double g, s1, s2;
};

namespace py = pybind11;

struct recorder_with_func
{
    mutable std::vector<genetic_values_per_diploid> data;
    //Our C++ recorder stores a python function
    py::function f;
    mutable genetic_values_per_diploid d;
    //We need to construct it with a function/callable:
    recorder_with_func(const std::uint32_t simlen, py::function f_)
        : data(std::vector<genetic_values_per_diploid>()), f(f_)
    {
        data.reserve(simlen);
    }

    void
    operator()(const fwdpy11::singlepop_t& pop) const
    {
        std::uint32_t i = 0;
        d.generation = pop.generation;
        for (auto&& dip : pop.diploids)
            {
                d.individual = i++;
                d.g = dip.g;
                d.s1=d.s2=0.;
                for (auto&& key : pop.gametes[dip.first].smutations)
                    {
                        d.s1 += pop.mutations[key].s;
                    }
                for (auto&& key : pop.gametes[dip.second].smutations)
                    {
                        d.s2 += pop.mutations[key].s;
                    }
                data.push_back(d);
            }
        //After we do our analysis, we send the data
        //to the Python function and clear out our 
        //C++ vector.
        f(&data);
        data.clear();
    }
};


PYBIND11_MAKE_OPAQUE(std::vector<genetic_values_per_diploid>);

PYBIND11_PLUGIN(sampler_pyfunction)
{
    pybind11::module m(
        "sampler_pyfunction",
        "Example of temporal sampler in C++ that holds a Python function.");

    PYBIND11_NUMPY_DTYPE(genetic_values_per_diploid, generation, individual,
                         g,s1,s2);

    py::bind_vector<std::vector<genetic_values_per_diploid>>(
        m, "VecGeneticValues", py::buffer_protocol());
    py::class_<recorder_with_func>(m, "cppRecorderFunc")
        //The constructor for our Python type needs
        //to have a function passed to it:
        .def(py::init<std::uint32_t, py::function>())
        .def_readonly("data", &recorder_with_func::data)
        .def("__call__",
             [](const recorder_with_func& r,
                const fwdpy11::singlepop_t& pop) -> void { r(pop); });

    return m.ptr();
}
