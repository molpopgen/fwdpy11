#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/regions/RecombinationRegions.hpp>

namespace py = pybind11;

void
init_RecombinationRegions(py::module& m)
{
    py::class_<fwdpy11::GeneticMap>(m, "GeneticMap", "ABC for genetic maps");

    py::class_<fwdpy11::RecombinationRegions, fwdpy11::GeneticMap>(
        m, "RecombinationRegions")
        .def(py::init<double, std::vector<fwdpy11::Region>>())
        .def_readonly("weights", &fwdpy11::RecombinationRegions::weights);

    py::class_<fwdpy11::GeneralizedGeneticMap, fwdpy11::GeneticMap>(
        m, "GeneralizedGeneticMap")
        .def("_num_poisson_callbacks",
             [](const fwdpy11::GeneralizedGeneticMap& self) {
                 return self.poisson_callbacks.size();
             })
        .def("_num_non_poisson_callbacks",
             [](const fwdpy11::GeneralizedGeneticMap& self) {
                 return self.non_poisson_callbacks.size();
             });

    m.def("dispatch_create_GeneticMap",
          [](py::object o, std::vector<fwdpy11::Region>& regions) {
              if (regions.empty() && o.is_none())
                  {
                      return fwdpy11::RecombinationRegions(0.0, regions);
                  }
              return fwdpy11::RecombinationRegions(o.cast<double>(), regions);
          });

    m.def("dispatch_create_GeneticMap_non_Region",
          [](py::list poisson, py::list non_poisson) {
              std::vector<std::unique_ptr<fwdpy11::PoissonCrossoverGenerator>>
                  poisson_callbacks;
              std::vector<std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator>>
                  non_poisson_callbacks;
              for (auto& i : poisson)
                  {
                      auto& ref = i.cast<fwdpy11::PoissonCrossoverGenerator&>();
                      poisson_callbacks.emplace_back(ref.ll_clone());
                  }
              for (auto& i : non_poisson)
                  {
                      auto& ref = i.cast<fwdpy11::NonPoissonCrossoverGenerator&>();
                      non_poisson_callbacks.emplace_back(ref.ll_clone());
                  }
              std::unique_ptr<fwdpy11::GeneralizedGeneticMap> rv(
                  new fwdpy11::GeneralizedGeneticMap(std::move(poisson_callbacks),
                                                     std::move(non_poisson_callbacks)));
              return rv;
          });
}
