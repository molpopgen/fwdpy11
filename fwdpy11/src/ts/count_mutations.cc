#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/Population.hpp>
#include <fwdpp/ts/count_mutations.hpp>

namespace py = pybind11;

void
init_count_mutations(py::module& m)
{
    m.def("count_mutations",
          [](const fwdpy11::Population& pop,
             const std::vector<fwdpp::ts::TS_NODE_INT>& samples) {
              decltype(pop.mcounts) mc;
              mc.resize(pop.mutations.size());
              fwdpp::ts::count_mutations(pop.tables, pop.mutations, samples,
                                         mc);
              return mc;
          });
}

