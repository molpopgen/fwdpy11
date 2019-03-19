#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/fitness/fitness.hpp>
#include <pybind11/pybind11.h>
#include <numeric>
#include <cmath>

namespace py = pybind11;

struct additive_diploid_trait_fxn
{
    const fwdpp::additive_diploid w;
    additive_diploid_trait_fxn(const double scaling)
        : w{ fwdpp::additive_diploid(scaling,
                                     fwdpp::additive_diploid::policy::atrait) }
    {
    }
    inline double
    operator()(const fwdpy11::Diploid &dip,
               const fwdpy11::Population::gcont_t &gametes,
               const fwdpy11::Population::mcont_t &mutations) const
    {
        return w(dip, gametes, mutations);
    }
};

struct multiplicative_diploid_trait_fxn
{
    const fwdpp::multiplicative_diploid w;
    multiplicative_diploid_trait_fxn(const double scaling)
        : w{ fwdpp::multiplicative_diploid(
              scaling, fwdpp::multiplicative_diploid::policy::mtrait) }
    {
    }
    inline double
    operator()(const fwdpy11::Diploid &dip,
               const fwdpy11::Population::gcont_t &gametes,
               const fwdpy11::Population::mcont_t &mutations) const
    {
        return w(dip, gametes, mutations);
    }
};

struct gbr_diploid_trait_fxn
{
    gbr_diploid_trait_fxn() {}
    gbr_diploid_trait_fxn(double) {} // hack to make API happy
    inline double
    operator()(const fwdpy11::Diploid &dip,
               const fwdpy11::Population::gcont_t &gametes,
               const fwdpy11::Population::mcont_t &mutations) const
    {
        auto sum1 = std::accumulate(
            gametes[dip.first].smutations.cbegin(),
            gametes[dip.first].smutations.cend(), 0.,
            [&mutations](const double s, const fwdpp::uint_t key) {
                return s + mutations[key].s;
            });
        auto sum2 = std::accumulate(
            gametes[dip.second].smutations.cbegin(),
            gametes[dip.second].smutations.cend(), 0.,
            [&mutations](const double s, const fwdpp::uint_t key) {
                return s + mutations[key].s;
            });
        return std::sqrt(sum1 * sum2);
    }
};

using single_locus_multiplicative_trait_wrapper
    = fwdpy11::fwdpp_single_locus_fitness_wrapper<
        multiplicative_diploid_trait_fxn>;
using single_locus_additive_trait_wrapper
    = fwdpy11::fwdpp_single_locus_fitness_wrapper<additive_diploid_trait_fxn>;
using gbr_trait_wrapper
    = fwdpy11::fwdpp_single_locus_fitness_wrapper<gbr_diploid_trait_fxn>;

PYBIND11_MODULE(trait_values, m)
{
    m.doc() = "Trait values.";

    FWDPY11_SINGLE_LOCUS_FITNESS()

    py::class_<single_locus_additive_trait_wrapper,
               std::shared_ptr<single_locus_additive_trait_wrapper>,
               fwdpy11::single_locus_fitness>(m, "DiploidAdditiveTrait",
                                              R"delim(
                Additive trait value, centered on zero.

                Trait value is :math:`\sum_{i} x_i`, 
                where :math:`x_i = 0, sh,\ \mathrm{or\ }scaling \times s`
                for AA, Aa, and aa, respectively.
                )delim")
        .def(py::init<double>(), py::arg("scaling"))
        .def_readonly("scaling", &single_locus_additive_trait_wrapper::scaling,
                      "Get the scaling attribute.")
        .def("__repr__",
             [](const single_locus_additive_trait_wrapper &a) {
                 std::string rv = "trait_values.DiploidAdditiveTrait(";
                 rv += std::to_string(a.scaling) + ")";
                 return rv;
             })
        .def(py::pickle(
            [](const single_locus_additive_trait_wrapper &w) {
                return py::make_tuple(w.scaling);
            },
            [](py::tuple t) {
                return std::make_shared<single_locus_additive_trait_wrapper>(
                    t[0].cast<double>());
            }));

    py::class_<single_locus_multiplicative_trait_wrapper,
               std::shared_ptr<single_locus_multiplicative_trait_wrapper>,
               fwdpy11::single_locus_fitness>(m, "DiploidMultTrait",
                                              R"delim(
                Multiplicative trait value, centered on zero.
                )delim")
        .def(py::init<double>(), py::arg("scaling"))
        .def_readonly("scaling",
                      &single_locus_multiplicative_trait_wrapper::scaling,
                      "Get the scaling attribute.")
        .def("__repr__",
             [](const single_locus_multiplicative_trait_wrapper &a) {
                 std::string rv = "trait_values.DiploidMultTrait(";
                 rv += std::to_string(a.scaling) + ")";
                 return rv;
             })
        .def(py::pickle(
            [](const single_locus_multiplicative_trait_wrapper &w) {
                return py::make_tuple(w.scaling);
            },
            [](py::tuple t) {
                return std::make_shared<
                    single_locus_multiplicative_trait_wrapper>(
                    t[0].cast<double>());
            }));

    py::class_<gbr_trait_wrapper, std::shared_ptr<gbr_trait_wrapper>,
               fwdpy11::single_locus_fitness>(m, "DiploidGBRTrait",
                                              R"delim(
            The "gene-based recessive" model from Thornton et al.
            2013 http://dx.doi.org/10.1371/journal.pgen.1003258 
            and Sanjak et al. 2017 http://dx.doi.org/10.1371/journal.pgen.1006573.

            The trait value is the geometric mean of the sum of effect sizes on
            each haplotype.  It is undefined for the case where these sums are negative.
            )delim")
        .def(py::init<>())
        .def("__repr__",
             [](const gbr_trait_wrapper &a) {
                 std::string rv = "trait_values.DiploidGBRTrait()";
                 return rv;
             })
        .def(py::pickle(
            [](const gbr_trait_wrapper &g) {
                return py::make_tuple("gbr_trait_wrapper");
            },
            [](py::tuple t) {
                std::string n = t[0].cast<std::string>();
                if (n == "gbr_trait_wrapper")
                    {
                        return std::make_shared<gbr_trait_wrapper>();
                    }
                else
                    {
                        throw std::invalid_argument(
                            "incorrect type name found when unpickling");
                    }
            }));
}
