#include <fwdpy11/types.hpp>
#include <fwdpy11/fitness/fitness.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

struct additive_diploid_trait_fxn
{
    inline double
    operator()(const fwdpy11::diploid_t &dip, const fwdpy11::gcont_t &gametes,
               const fwdpy11::mcont_t &mutations, const double scaling) const
    {
        return KTfwd::site_dependent_fitness()(
            dip, gametes, mutations,
            [scaling](double &fitness,
                      const fwdpy11::mcont_t::value_type &mut) noexcept {
                fitness += (scaling * mut.s);
            },
            [](double &fitness,
               const fwdpy11::mcont_t::value_type &mut) noexcept {
                fitness += (mut.h * mut.s);
            },
            0.);
    }
};

struct multiplicative_diploid_trait_fxn
{
    inline double
    operator()(const fwdpy11::diploid_t &dip, const fwdpy11::gcont_t &gametes,
               const fwdpy11::mcont_t &mutations, const double scaling) const
    {
        return KTfwd::site_dependent_fitness()(
                   dip, gametes, mutations,
                   [scaling](
                       double &fitness,
                       const fwdpy11::mcont_t::value_type &mut) noexcept {
                       fitness *= (1. + scaling * mut.s);
                   },
                   [](double &fitness,
                      const fwdpy11::mcont_t::value_type &mut) noexcept {
                       fitness *= (1. + mut.h * mut.s);
                   },
                   1.)
               - 1.0;
    }
};

using singlepop_multiplicative_trait_wrapper = fwdpy11::
    fwdpp_singlepop_fitness_wrapper<multiplicative_diploid_trait_fxn>;
using singlepop_additive_trait_wrapper
    = fwdpy11::fwdpp_singlepop_fitness_wrapper<additive_diploid_trait_fxn>;

PYBIND11_PLUGIN(trait_values)
{
    py::module m("trait_values", "Trait values.");

    FWDPY11_SINGLEPOP_FITNESS()

    py::class_<singlepop_additive_trait_wrapper, fwdpy11::singlepop_fitness>(
        m, "SpopAdditiveTrait",
        R"delim(
                Additive trait value, centered on zero.

                Trait value is :math:`\sum_{i} x_i`), 
                where :math:`x_i = 0, sh,\ \mathrm{or\ }scaling \times s`
                for AA, Aa, and aa, respectively.
                )delim")
        .def(py::init<double>(), py::arg("scaling"));

    py::class_<singlepop_multiplicative_trait_wrapper,
               fwdpy11::singlepop_fitness>(m, "SpopMultTrait",
                                           R"delim(
                Multiplicative trait value, centered on zero.
                )delim")
        .def(py::init<double>(), py::arg("scaling"));
    return m.ptr();
}
