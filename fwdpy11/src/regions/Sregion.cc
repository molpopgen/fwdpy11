#include <pybind11/pybind11.h>
#include <fwdpy11/regions/Sregion.hpp>

namespace py = pybind11;

void
init_Sregion(py::module& m)
{
    py::class_<fwdpy11::Sregion>(m, "Sregion",
                                 R"delim(
        Representation of a "region" in a simulation with a dominance term.

        This class is the base class for a general set
        of objects representing distributions of fitness effects.

        Attributes:
            * b: the beginning of the region
            * e: the end of the region
            * w: the "weight" assigned to the region
            * h: the dominance term
            * l: A label assigned to the region.
                Labels must be integers, and can be
                used to 'tag' mutations arising in different regions.
            * scaling: The scaling of the distrubution.  See note below.

        See :func:`fwdpy11.wright_fisher.evolve` for how this class
        may be used to parameterize a simulation.

        .. note:: This class cannot be used directly to parameterize a simulation.
            Rather, you must use a derived type that specifies a
            distribution of fitness effects.  These types include:
           :class:`fwdpy11.ConstantS`,
           :class:`fwdpy11.UniformS`,
           :class:`fwdpy11.ExpS`,
           :class:`fwdpy11.GammaS`, and
           :class:`fwdpy11.GaussianS`

        .. note:: The scaling of a distribution refers to the distribution of effect sizes.
            For example, if scaling = 1.0, then the distribution is on the effect size itself.  If
            scaling = 2N (where N is the population size), then the DFE is on 2Ns.  If N
            is not constant during a simulation, then the scaling is with respect to some
            "reference" population size.

        .. versionchanged:: 0.13.a2
            Added "scaling" attribute.

        .. versionchanged:: 0.3.0
            Refactored from a pure Python class to a C++/pybind11 class

        )delim")
        .def_property_readonly(
            "b", [](const fwdpy11::Sregion& s) { return s.beg(); },
            "Beginning of region")
        .def_property_readonly(
            "e", [](const fwdpy11::Sregion& s) { return s.end(); },
            "End of region")
        .def_property_readonly(
            "w", [](const fwdpy11::Sregion& s) { return s.weight(); },
            "Weight")
        .def_property_readonly(
            "l", [](const fwdpy11::Sregion& s) { return s.label(); }, "Label")
        .def_property_readonly(
            "c", [](const fwdpy11::Sregion& s) { return s.region.coupled; },
            "Coupling parameter")
        .def_readonly("scaling", &fwdpy11::Sregion::scaling,
                      "Scaling parameter");
}

