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
    if (gssmo.current_timepoint >= gssmo.optima.size())
        {
            throw std::runtime_error("current_timepoint out of range");
        }
    return gssmo.optima[gssmo.current_timepoint].optima;
}

PYBIND11_MODULE(testMultivariateGSSmo, m)
{
    m.def("update", &update);
    m.def("get_optima", &get_optima);
}
