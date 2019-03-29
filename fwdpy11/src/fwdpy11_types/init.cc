#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_Mutation(py::module & m);
void init_MutationVector(py::module & m);
void init_DiploidGenotype(py::module &m);
void init_DiploidMetadata(py::module &m);
void init_DiploidVector(py::module & m);
void init_HaploidGenomeVector(py::module & m);
void init_rng(py::module &);
void init_PopulationBase(py::module & m);
void init_DiploidPopulation(py::module & m);
void init_tsrecorders(py::module & m);
void
init_RecordNothing(pybind11::module &);

void initialize_fwdpy11_types(py::module & m)
{
    init_Mutation(m);
    init_MutationVector(m);
    init_DiploidGenotype(m);
    init_DiploidMetadata(m);
    init_DiploidVector(m);
    init_HaploidGenomeVector(m);
    init_rng(m);
    init_PopulationBase(m);
    init_DiploidPopulation(m);
    init_RecordNothing(m);
    init_tsrecorders(m);
}
