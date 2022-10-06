#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <fwdpy11/types/Population.hpp>

struct RecordNothing
{
    inline void
    operator()(const fwdpy11::Population &) const
    {
    }
};

void
init_RecordNothing(pybind11::module &m)
{
    pybind11::class_<RecordNothing>(m, "RecordNothing")
        .def(pybind11::init<>())
        .def("__call__", &RecordNothing::operator());
}
