#include <pybind11/pybind11.h>
#include <fwdpy11/types/Mutation.hpp>

namespace py = pybind11;

//Example of pickling a specific C++ type
py::bytes
pickle_mutation(const fwdpy11::Mutation& p)
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
    return f(p,-1);
}

PYBIND11_MODULE(pickling_cpp, m)
{
    m.def("pickle_mutation", &pickle_mutation);
    m.def("general_pickler", &general_pickler);
}
