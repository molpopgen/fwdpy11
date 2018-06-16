#ifndef FWDPY11_GENETIC_VALUES_SLOCUSMULT_HPP__
#define FWDPY11_GENETIC_VALUES_SLOCUSMULT_HPP__

#include <fwdpp/fitness_models.hpp>
#include "fwdpp_wrappers/fwdpp_slocus_gvalue.hpp"

namespace fwdpy11
{
    template <>
    fwdpp_slocus_gvalue<fwdpp::multiplicative_diploid>::fwdpp_slocus_gvalue(
        const double scaling)
        : fwdpy11::SlocusPopGeneticValueWithMapping{ gvalue_map_ptr(
              new fwdpy11::GeneticValueIsFitness()) },
          gv{ scaling, fwdpp::multiplicative_diploid::policy::mw }
    {
        if (!std::isfinite(scaling))
            {
                throw std::invalid_argument("scaling must be finite");
            }
    }

    template <>
    fwdpp_slocus_gvalue<fwdpp::multiplicative_diploid>::fwdpp_slocus_gvalue(
        const double scaling, const fwdpy11::GeneticValueIsTrait& g2w)
        : fwdpy11::SlocusPopGeneticValueWithMapping{ g2w.clone() }, gv{
              scaling, fwdpp::multiplicative_diploid::policy::mtrait
          }
    {
    }

    template <>
    fwdpp_slocus_gvalue<fwdpp::multiplicative_diploid>::fwdpp_slocus_gvalue(
        const double scaling, const fwdpy11::GeneticValueIsTrait& g2w,
        const fwdpy11::GeneticValueNoise& noise_fxn)
        : fwdpy11::SlocusPopGeneticValueWithMapping{ g2w.clone(),
                                                     noise_fxn.clone() },
          gv{ scaling, fwdpp::multiplicative_diploid::policy::mtrait }
    {
    }

    using SlocusMult = fwdpp_slocus_gvalue<fwdpp::multiplicative_diploid>;
} // namespace fwdpy11
#endif
