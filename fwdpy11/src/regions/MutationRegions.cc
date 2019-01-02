#include <pybind11/pybind11.h>
#include <fwdpy11/regions/MutationRegions.hpp>
#include <fwdpy11/regions/ConstantS.hpp>

namespace py = pybind11;

void
init_MutationRegions(py::module& m)
{
    py::class_<fwdpy11::MutationRegions>(m, "MutationRegions")
        .def_static("create",
                    [](double pneutral,
                       const std::vector<fwdpy11::Region>& neutral,
                       py::list selected) -> fwdpy11::MutationRegions {
                        std::vector<std::unique_ptr<fwdpy11::Sregion>> regions;

                        std::vector<double> nweights, sweights;

                        for (auto& n : neutral)
                            {
                                nweights.push_back(n.weight);
                                regions.emplace_back(new fwdpy11::ConstantS(
                                    n.beg, n.end, n.weight, 0.0, 0.0,
                                    n.coupled, n.label, 1.0));
                            }

                        std::vector<double> weights;

                        return fwdpy11::MutationRegions(std::move(regions),
                                                        std::move(weights));
                    });
}
