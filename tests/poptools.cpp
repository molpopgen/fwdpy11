// Tools for interacting with population types
// in ways not allowed by the Python API (for
// reasons of safety).
#include <pybind11/pybind11.h>
#include <fwdpy11/types/DiploidPopulation.hpp>

void
change_generation(fwdpy11::Population& pop, std::uint32_t g)
{
    pop.generation = g;
}

PYBIND11_MODULE(poptools, m)
{
    m.def("change_generation", &change_generation);
}
