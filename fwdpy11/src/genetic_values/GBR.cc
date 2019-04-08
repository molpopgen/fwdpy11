#include <fwdpy11/genetic_values/DiploidGBR.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

static const auto GBR_CONSTRUCTOR1 =
    R"delim(
 Construct object with specific genetic value to fitness map.
 
 :param gv2w: Genetic value to fitness map
 :type gv2w: :class:`fwdpy11.GeneticValueIsTrait`
 )delim";

static const auto GBR_CONSTRUCTOR2 =
    R"delim(
Construct object with specific genetic value to fitness map 
and random effects on trait value.

:param gv2w: Genetic value to fitness map
:type gv2w: :class:`fwdpy11.GeneticValueIsTrait`
:param noise: Model of random effects on trait value.
:type noise: :class:`fwdpy11.GeneticValueNoise`
)delim";

void
init_GBR(py::module& m)
{
    py::class_<fwdpy11::DiploidGBR,
               fwdpy11::DiploidPopulationGeneticValueWithMapping>(m,
                                                                  "GBR",
                                                                  R"delim(
        The "gene-based recessive" trait model described in Thornton et al.
        2013 http://dx.doi.org/10.1371/journal.pgen.1003258 and Sanjak et al. 2017
        http://dx.doi.org/10.1371/journal.pgen.1006573.

        The trait value is the geometric mean of the sum of effect sizes on each haplotype.
        It is undefined for the case where these sums are negative.
        )delim")
        .def(py::init([](const fwdpy11::GeneticValueIsTrait& gv2w) {
                 return fwdpy11::DiploidGBR(fwdpy11::GBR{}, gv2w);
             }),
             py::arg("gv2w"), GBR_CONSTRUCTOR1)
        .def(py::init([](const fwdpy11::GeneticValueIsTrait& gv2w,
                         const fwdpy11::GeneticValueNoise& noise) {
                 return fwdpy11::DiploidGBR(fwdpy11::GBR{}, gv2w, noise);
             }),
             py::arg("gv2w"), py::arg("noise"), GBR_CONSTRUCTOR2)
        .def(py::pickle(
            [](const fwdpy11::DiploidGBR& g) {
                auto p = py::module::import("pickle");
                return py::make_tuple(
                    g.pickle(), p.attr("dumps")(g.gv2w->clone(), -1),
                    p.attr("dumps")(g.noise_fxn->clone(), -1));
            },
            [](py::tuple t) {
                if (t.size() != 3)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                std::string s = t[0].cast<std::string>();
                if (s != "GBR")
                    {
                        throw std::runtime_error("invalid object state");
                    }
                auto p = py::module::import("pickle");
                auto t1 = p.attr("loads")(t[1]);
                auto t2 = p.attr("loads")(t[2]);
                //Do the casts in the constructor
                //to avoid any nasty issues w/
                //refs to temp
                return fwdpy11::DiploidGBR(
                    fwdpy11::GBR{},
                    t1.cast<const fwdpy11::GeneticValueIsTrait&>(),
                    t2.cast<const fwdpy11::GeneticValueNoise&>());
            }));
}
