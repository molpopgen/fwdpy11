#ifndef FWDPY11_GENETIC_VALUES_SLOCUSGBR_HPP__
#define FWDPY11_GENETIC_VALUES_SLOCUSGBR_HPP__

#include "fwdpp_wrappers/fwdpp_slocus_gvalue.hpp"
#include "details/GBR.hpp"

namespace fwdpy11
{
    template <>
    fwdpp_slocus_gvalue<GBR>::fwdpp_slocus_gvalue(
        const fwdpy11::GeneticValueIsTrait& g2w)
        : fwdpy11::SlocusPopGeneticValueWithMapping{ g2w.clone() }, gv{}
    {
    }

    template <>
    fwdpp_slocus_gvalue<GBR>::fwdpp_slocus_gvalue(
        const fwdpy11::GeneticValueIsTrait& g2w,
        const fwdpy11::GeneticValueNoise& noise_fxn)
        : fwdpy11::SlocusPopGeneticValueWithMapping{ g2w.clone(),
                                                     noise_fxn.clone() },
          gv{}
    {
    }

    using SlocusGBR = fwdpp_slocus_gvalue<GBR>;
} // namespace fwdpy11

#endif
