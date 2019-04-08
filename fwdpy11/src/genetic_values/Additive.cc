#include <fwdpy11/genetic_values/DiploidAdditive.hpp>
#include <pybind11/pybind11.h>

namespace
{
    static const auto ADDITIVE_CONSTRUCTOR_1 =
        R"delim(
Additive effects on fitness.

:param scaling: How to treat mutant homozygotes.
:type scaling: float

For a model of fitness, the genetic value is 1, 1+e*h,
1+scaling*e for genotypes AA, Aa, and aa, respectively.
)delim";

    static const auto ADDITIVE_CONSTRUCTOR_2 =
        R"delim(
Construct an object of additive effects on a trait with a specific
functional mapping from genetic value to fitness.

:param scaling: How to treat mutant homozygotes.
:type scaling: float
:param gv2w: Map from genetic value to fitness.
:type gv2w: :class:`fwdpy11.GeneticValueIsTrait`
)delim";

    static const auto ADDITIVE_CONSTRUCTOR_3 =
        R"delim(
Additive effects on a trait with a specific mapping from 
genetic value to fitness and random effects ("noise").

:param scaling: How to treat mutant homozygotes.
:type scaling: float
:param gv2w: Map from genetic value to fitness.
:type gv2w: :class:`fwdpy11.GeneticValueIsTrait`
:param noise: Function to generate random effects on trait value.
:type noise: :class:`fwdpy11.GeneticValueNoise`
)delim";
} // namespace

namespace py = pybind11;

void
init_Additive(py::module& m)
{
    py::class_<fwdpy11::DiploidAdditive,
               fwdpy11::DiploidPopulationGeneticValueWithMapping>(
        m, "Additive", "Additive genetic values.")
        .def(py::init([](const double scaling) {
                 return fwdpy11::DiploidAdditive(
                     fwdpp::additive_diploid(fwdpp::fitness(scaling)));
             }),
             py::arg("scaling"), ADDITIVE_CONSTRUCTOR_1)
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& g) {
                 return fwdpy11::DiploidAdditive(
                     fwdpp::additive_diploid(fwdpp::trait(scaling)), g);
             }),
             py::arg("scaling"), py::arg("gv2w"), ADDITIVE_CONSTRUCTOR_2)
        .def(py::init([](const double scaling,
                         const fwdpy11::GeneticValueIsTrait& g,
                         const fwdpy11::GeneticValueNoise& n) {
                 return fwdpy11::DiploidAdditive(
                     fwdpp::additive_diploid(fwdpp::trait(scaling)), g, n);
             }),
             py::arg("scaling"), py::arg("gv2w"), py::arg("noise"),
             ADDITIVE_CONSTRUCTOR_3)
        .def_property_readonly(
            "scaling",
            [](const fwdpy11::DiploidAdditive& wa) { return wa.gv.scaling; },
            "Access to the scaling parameter.")
        .def_property_readonly(
            "is_fitness",
            [](const fwdpy11::DiploidAdditive& wa) {
                return wa.gv.gvalue_is_fitness;
            },
            "Returns True if instance calculates fitness as the genetic value "
            "and False if the genetic value is a trait value.")
        .def(py::pickle(
            [](const fwdpy11::DiploidAdditive& a) {
                auto p = py::module::import("pickle");
                return py::make_tuple(
                    a.pickle(), p.attr("dumps")(a.gv2w->clone(), -1),
                    p.attr("dumps")(a.noise_fxn->clone(), -1));
            },
            [](py::tuple t) {
                if (t.size() != 3)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                auto t0 = t[0].cast<py::tuple>();
                int pol = t0[0].cast<int>();

                double scaling = t0[1].cast<double>();
                auto p = py::module::import("pickle");
                auto a
                    = (pol == 1)
                          ? fwdpp::additive_diploid(fwdpp::trait(scaling))
                          : fwdpp::additive_diploid(fwdpp::fitness(scaling));
                auto t1 = p.attr("loads")(t[1]);
                auto t2 = p.attr("loads")(t[2]);
                //Do the casts in the constructor
                //to avoid any nasty issues w/
                //refs to temp
                return fwdpy11::DiploidAdditive(
                    std::move(a),
                    t1.cast<const fwdpy11::GeneticValueToFitnessMap&>(),
                    t2.cast<const fwdpy11::GeneticValueNoise&>());
            }));
}
