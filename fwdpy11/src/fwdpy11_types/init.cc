#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_Mutation(py::module & m);
void init_MutationVector(py::module & m);
void init_DiploidGenotype(py::module &m);
void init_DiploidMetadata(py::module &m);
void init_DiploidVector(py::module & m);
void init_GameteVector(py::module & m);
void init_rng(py::module &);

void initialize_fwdpy11_types(py::module & m)
{
    init_Mutation(m);
    init_MutationVector(m);
    init_DiploidGenotype(m);
    init_DiploidMetadata(m);
    init_DiploidVector(m);
    init_GameteVector(m);
    init_rng(m);
}
