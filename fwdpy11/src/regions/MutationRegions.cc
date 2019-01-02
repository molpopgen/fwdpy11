#include <algorithm>
#include <numeric>
#include <pybind11/pybind11.h>
#include <fwdpy11/regions/MutationRegions.hpp>
#include <fwdpy11/regions/ConstantS.hpp>

namespace py = pybind11;

void
init_MutationRegions(py::module& m)
{
    py::class_<fwdpy11::MutationRegions>(m, "MutationRegions")
        .def_static(
            "create",
            [](double pneutral, const std::vector<fwdpy11::Region>& neutral,
               py::list selected) -> fwdpy11::MutationRegions {
                std::vector<std::unique_ptr<fwdpy11::Sregion>> regions;

                std::vector<double> nweights, sweights;

                for (auto& n : neutral)
                    {
                        nweights.push_back(n.weight);
                        regions.emplace_back(new fwdpy11::ConstantS(
                            n.beg, n.end, n.weight, 0.0, 0.0, n.coupled,
                            n.label, 1.0));
                    }

                for (auto& s : selected)
                    {
                        auto& ref = s.cast<fwdpy11::Sregion&>();
                        sweights.push_back(ref.weight());
                        regions.emplace_back(ref.clone());
                    }

                // Have to reweight input weights
                double sum_neutral
                    = std::accumulate(begin(nweights), end(nweights), 0.0);
                double sum_selected
                    = std::accumulate(begin(sweights), end(sweights), 0.0);

                std::transform(
                    begin(nweights), end(nweights), begin(nweights),
                    [sum_neutral](double d) { return d / sum_neutral; });
                std::transform(begin(nweights), end(nweights), begin(nweights),
                               [pneutral](double d) { return d * pneutral; });
                std::transform(
                    begin(sweights), end(sweights), begin(sweights),
                    [sum_selected](double d) { return d / sum_selected; });
                std::transform(
                    begin(sweights), end(sweights), begin(sweights),
                    [pneutral](double d) { return d * (1. - pneutral); });

                std::vector<double> weights(begin(nweights), end(nweights));
                weights.insert(end(weights), begin(sweights), end(sweights));

                return fwdpy11::MutationRegions(std::move(regions),
                                                std::move(weights));
            });
}
