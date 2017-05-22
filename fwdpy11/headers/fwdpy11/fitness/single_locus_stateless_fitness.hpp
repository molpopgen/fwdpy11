#ifndef FWDPY11_SINGLE_LOCUS_STATELESS_FITNESS_HPP_
#define FWDPY11_SINGLE_LOCUS_STATELESS_FITNESS_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include "single_locus_fitness.hpp"

namespace fwdpy11
{
    struct single_locus_stateless_fitness : public single_locus_fitness
    {
        const single_locus_fitness_fxn ff;
        single_locus_stateless_fitness(single_locus_fitness_fxn ff_)
            : ff{ std::move(ff_) }
        {
        }
        inline single_locus_fitness_fxn
        callback() const final
        {
            return ff;
        }

        SINGLE_LOCUS_FITNESS_CLONE_SHARED(single_locus_stateless_fitness);
        SINGLE_LOCUS_FITNESS_CLONE_UNIQUE(single_locus_stateless_fitness);
        SINGLE_LOCUS_FITNESS_CALLBACK_NAME(typeid(*this).name());
    };

#define STATELESS_SLOCUS_FUNCTION(FUNCTION)                                   \
    double FUNCTION(const fwdpy11::diploid_t& dip,                            \
                    const fwdpy11::gcont_t& gametes,                          \
                    const fwdpy11::mcont_t& mutations)

#define FWDPY11_SINGLE_LOCUS_STATELESS_FITNESS()                              \
    FWDPY11_SINGLE_LOCUS_FITNESS()                                            \
    pybind11::object                                                          \
        FWDPY11_SINGLE_LOCUS_CUSTOM_STATELESS_FITNESS_BASE_IMPORT__           \
        = (pybind11::object)pybind11::module::import("fwdpy11.fitness")       \
              .attr("SlocusCustomStatelessGeneticValue");

#define STATELESS_GENOTYPE_POLICY(FUNCTION)                                   \
    void FUNCTION(double& w, const KTfwd::popgenmut& m)

#define CREATE_STATELESS_SLOCUS_OBJECT(FUNCTION, FUNCNAME, MODOBJ)            \
    MODOBJ.def(FUNCNAME, []() {                                               \
        return fwdpy11::single_locus_stateless_fitness(FUNCTION);             \
    });

#define CREATE_STATELESS_SLOCUS_GENOTYPE_OBJECT(HOM, HET, SVALUE, FUNCNAME,   \
                                                MODOBJ)                       \
    MODOBJ.def(FUNCNAME, []() {                                               \
        return fwdpy11::single_locus_stateless_fitness(std::bind(             \
            KTfwd::site_dependent_genetic_value(), std::placeholders::_1,     \
            std::placeholders::_2, std::placeholders::_3, HOM, HET, SVALUE)); \
    });

}
#endif
