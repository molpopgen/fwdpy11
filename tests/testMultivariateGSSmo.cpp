// Low-level tests of the C++ API for this
// class that is not exposed to Python
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/genetic_values/MultivariateGSSmo.hpp>

void
update(fwdpy11::DiploidPopulation& pop, fwdpy11::MultivariateGSSmo& gssmo)
{
    gssmo.update(pop);
}

std::vector<double>
get_optima(const fwdpy11::MultivariateGSSmo& gssmo)
{
    return std::vector<double>(begin(gssmo.optima) + gssmo.optima_offset,
                               begin(gssmo.optima) + gssmo.optima_offset
                                   + gssmo.ndim);
}

PYBIND11_MODULE(testMultivariateGSSmo, m)
{
    m.def("update", &update);
    m.def("get_optima", &get_optima);
}
