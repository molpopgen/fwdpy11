#include <algorithm>
#include <numeric>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpy11/regions/MutationRegions.hpp>
#include <fwdpy11/regions/ConstantS.hpp>

namespace py = pybind11;

void
init_MutationRegions(py::module& m)
{
    py::class_<fwdpy11::MutationRegions>(m, "MutationRegions")
        .def_readonly("weights", &fwdpy11::MutationRegions::weights)
        .def_static(
            "create",
            [](double pneutral, const std::vector<fwdpy11::Region>& neutral,
               py::list selected) -> fwdpy11::MutationRegions {
                std::vector<std::unique_ptr<fwdpy11::Sregion>> nregions, sregions;

                std::vector<double> nweights, sweights;

                for (auto& n : neutral)
                    {
                        nweights.push_back(n.weight);
                        nregions.emplace_back(new fwdpy11::ConstantS(fwdpy11::Region(
                            n.beg, n.end, n.weight, n.coupled, n.label)));
                    }

                for (auto& s : selected)
                    {
                        auto& ref = s.cast<fwdpy11::Sregion&>();
                        sweights.push_back(ref.weight());
                        sregions.emplace_back(ref.clone());
                    }

                return fwdpy11::MutationRegions::create(pneutral, nweights, sweights,
                                                        nregions, sregions);
            });
}
