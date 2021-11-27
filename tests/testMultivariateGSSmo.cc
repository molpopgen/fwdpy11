// Low-level tests of the C++ API for this
// class that is not exposed to Python
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/genetic_value_to_fitness/MultivariateGSSmo.hpp>

void
update(fwdpy11::DiploidPopulation& pop, fwdpy11::MultivariateGSSmo& gssmo)
{
    gssmo.update(pop);
}

std::vector<double>
get_optima(const fwdpy11::MultivariateGSSmo& gssmo)
{
    return gssmo.current_timepoint_optima;
}

PYBIND11_MODULE(testMultivariateGSSmo, m)
{
    m.def("update", &update);
    m.def("get_optima", &get_optima);
}
