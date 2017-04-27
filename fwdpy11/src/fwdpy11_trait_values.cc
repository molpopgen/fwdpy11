#include <fwdpy11/types.hpp>
#include <fwdpy11/fitness/fitness.hpp>
#include <pybind11/pybind11.h>
#include <numeric>
#include <cmath>

namespace py = pybind11;

struct additive_diploid_trait_fxn
{
    const KTfwd::additive_diploid w;
    additive_diploid_trait_fxn()
        : w{ KTfwd::additive_diploid(KTfwd::atrait()) }
    {
    }
    inline double
    operator()(const fwdpy11::diploid_t &dip, const fwdpy11::gcont_t &gametes,
               const fwdpy11::mcont_t &mutations, const double scaling) const
    {
        return w(dip, gametes, mutations, scaling);
    }
};

struct multiplicative_diploid_trait_fxn
{
    const KTfwd::multiplicative_diploid w;
    multiplicative_diploid_trait_fxn()
        : w{ KTfwd::multiplicative_diploid(KTfwd::mtrait()) }
    {
    }
    inline double
    operator()(const fwdpy11::diploid_t &dip, const fwdpy11::gcont_t &gametes,
               const fwdpy11::mcont_t &mutations, const double scaling) const
    {
        return w(dip, gametes, mutations, scaling);
    }
};

struct gbr_diploid_trait_fxn
{
    inline double
    operator()(const fwdpy11::diploid_t &dip, const fwdpy11::gcont_t &gametes,
               const fwdpy11::mcont_t &mutations, const double) const
    {
        auto sum1 = std::accumulate(
            gametes[dip.first].smutations.cbegin(),
            gametes[dip.first].smutations.cend(), 0.,
            [&mutations](const double s, const KTfwd::uint_t key) {
                return s + mutations[key].s;
            });
        auto sum2 = std::accumulate(
            gametes[dip.second].smutations.cbegin(),
            gametes[dip.second].smutations.cend(), 0.,
            [&mutations](const double s, const KTfwd::uint_t key) {
                return s + mutations[key].s;
            });
        return std::sqrt(sum1 * sum2);
    }
};

using single_locus_multiplicative_trait_wrapper = fwdpy11::
    fwdpp_single_locus_fitness_wrapper<multiplicative_diploid_trait_fxn>;
using single_locus_additive_trait_wrapper
    = fwdpy11::fwdpp_single_locus_fitness_wrapper<additive_diploid_trait_fxn>;
using gbr_trait_wrapper
    = fwdpy11::fwdpp_single_locus_fitness_wrapper<gbr_diploid_trait_fxn>;

PYBIND11_PLUGIN(trait_values)
{
    py::module m("trait_values", "Trait values.");

    FWDPY11_SINGLE_LOCUS_FITNESS()

    py::class_<single_locus_additive_trait_wrapper,
               std::shared_ptr<single_locus_additive_trait_wrapper>,
               fwdpy11::single_locus_fitness>(m, "SlocusAdditiveTrait",
                                           R"delim(
                Additive trait value, centered on zero.

                Trait value is :math:`\sum_{i} x_i`, 
                where :math:`x_i = 0, sh,\ \mathrm{or\ }scaling \times s`
                for AA, Aa, and aa, respectively.
                )delim")
        .def(py::init<double>(), py::arg("scaling"));

    py::class_<single_locus_multiplicative_trait_wrapper,
               std::shared_ptr<single_locus_multiplicative_trait_wrapper>,
               fwdpy11::single_locus_fitness>(m, "SlocusMultTrait",
                                           R"delim(
                Multiplicative trait value, centered on zero.
                )delim")
        .def(py::init<double>(), py::arg("scaling"));

    py::class_<gbr_trait_wrapper, std::shared_ptr<gbr_trait_wrapper>,
               fwdpy11::single_locus_fitness>(m, "SlocusGBRTrait",
                                           R"delim(
            The "gene-based recessive" model from Thornton et al.
            2013 http://dx.doi.org/10.1371/journal.pgen.1003258 
            and Sanjak et al. 2017 http://dx.doi.org/10.1371/journal.pgen.1006573.

            The trait value is the geometric mean of the sum of effect sizes on
            each haplotype.  It is undefined for the case where these sums are negative.
            )delim")
        .def(py::init<>());

    return m.ptr();
}
