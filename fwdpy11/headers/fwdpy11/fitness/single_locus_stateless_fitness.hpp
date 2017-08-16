#ifndef FWDPY11_SINGLE_LOCUS_STATELESS_FITNESS_HPP_
#define FWDPY11_SINGLE_LOCUS_STATELESS_FITNESS_HPP_

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <memory>
#include <typeinfo>
#include "single_locus_fitness.hpp"

namespace fwdpy11
{
    template <typename callback_type>
    struct single_locus_stateless_fitness : public single_locus_fitness
    {
        inline single_locus_fitness_fxn
        callback() const final
        {
            return callback_type();
        }

        SINGLE_LOCUS_FITNESS_CLONE_SHARED(
            single_locus_stateless_fitness<callback_type>);
        SINGLE_LOCUS_FITNESS_CLONE_UNIQUE(
            single_locus_stateless_fitness<callback_type>);
        SINGLE_LOCUS_FITNESS_CALLBACK_NAME(typeid(*this).name());
    };

    template <typename AA, typename Aa, int starting_value>
    struct single_locus_stateless_fitness_genotype
        : public single_locus_fitness
    {
        inline single_locus_fitness_fxn
        callback() const final
        {
            return std::bind(KTfwd::site_dependent_genetic_value(),
                             std::placeholders::_1, std::placeholders::_2,
                             std::placeholders::_3, AA(), Aa(),
                             static_cast<double>(starting_value));
        }

        // We need this typedef so that the commas
        // don't confuse the macro calls that
        // come next.
        using this_type
            = single_locus_stateless_fitness_genotype<AA, Aa, starting_value>;

        SINGLE_LOCUS_FITNESS_CLONE_SHARED(this_type);
        SINGLE_LOCUS_FITNESS_CLONE_UNIQUE(this_type);
        SINGLE_LOCUS_FITNESS_CALLBACK_NAME(typeid(*this).name());
    };

#define STATELESS_SLOCUS_FUNCTION(FUNCTION)                                   \
    struct FUNCTION                                                           \
    {                                                                         \
        inline double operator()(const fwdpy11::diploid_t& dip,               \
                                 const fwdpy11::gcont_t& gametes,             \
                                 const fwdpy11::mcont_t& mutations) const

#define END_STRUCT()                                                          \
    }                                                                         \
    ;

#define FWDPY11_SINGLE_LOCUS_STATELESS_FITNESS()                              \
    FWDPY11_SINGLE_LOCUS_FITNESS()                                            \
    pybind11::object                                                          \
        FWDPY11_SINGLE_LOCUS_CUSTOM_STATELESS_FITNESS_BASE_IMPORT__           \
        = (pybind11::object)pybind11::module::import("fwdpy11.fitness")       \
              .attr("SlocusCustomStatelessGeneticValue");

#define STATELESS_GENOTYPE_POLICY(FUNCTION)                                   \
    struct FUNCTION                                                           \
    {                                                                         \
        inline void operator()(double& w, const KTfwd::popgenmut& m) const

#define CREATE_STATELESS_SLOCUS_OBJECT(FUNCTION, FUNCNAME, MODOBJ)            \
    pybind11::                                                                \
        class_<fwdpy11::single_locus_stateless_fitness<FUNCTION>,             \
               std::shared_ptr<fwdpy11::                                      \
                                   single_locus_stateless_fitness<FUNCTION>>, \
               fwdpy11::single_locus_fitness>(MODOBJ, FUNCNAME)               \
            .def(pybind11::init<>())                                          \
            .def("__getstate__",                                              \
                 [](const fwdpy11::single_locus_stateless_fitness<FUNCTION>&  \
                        ff) { return pybind11::make_tuple(FUNCNAME); })       \
            .def("__setstate__",                                              \
                 [](fwdpy11::single_locus_stateless_fitness<FUNCTION>& ff,    \
                    pybind11::tuple t) {                                      \
                     auto s = t[0].cast<std::string>();                       \
                     if (s != FUNCNAME)                                       \
                         {                                                    \
                             throw std::invalid_argument(                     \
                                 "incorrect type name found");                \
                         }                                                    \
                     new (&ff)                                                \
                         fwdpy11::single_locus_stateless_fitness<FUNCTION>(); \
                 });

#define CREATE_STATELESS_SLOCUS_GENOTYPE_OBJECT(HOM, HET, SVALUE, FUNCNAME,         \
                                                MODOBJ)                             \
    pybind11::                                                                      \
        class_<fwdpy11::single_locus_stateless_fitness_genotype<HOM, HET,           \
                                                                SVALUE>,            \
               std::                                                                \
                   shared_ptr<fwdpy11::                                             \
                                  single_locus_stateless_fitness_genotype<HOM,      \
                                                                          HET,      \
                                                                          SVALUE>>, \
               fwdpy11::single_locus_fitness>(MODOBJ, FUNCNAME)                     \
            .def(pybind11::init<>())                                                \
            .def(                                                                   \
                "__getstate__",                                                     \
                [](const fwdpy11::                                                  \
                       single_locus_stateless_fitness_genotype<HOM, HET,            \
                                                               SVALUE>& ff) {       \
                    return pybind11::make_tuple(FUNCNAME);                          \
                })                                                                  \
            .def(                                                                   \
                "__setstate__",                                                     \
                [](fwdpy11::single_locus_stateless_fitness_genotype<HOM, HET,       \
                                                                    SVALUE>&        \
                       ff,                                                          \
                   pybind11::tuple t) {                                             \
                    auto s = t[0].cast<std::string>();                              \
                    if (s != FUNCNAME)                                              \
                        {                                                           \
                            throw std::invalid_argument(                            \
                                "incorrect type name found");                       \
                        }                                                           \
                    new (&ff) fwdpy11::                                             \
                        single_locus_stateless_fitness_genotype<HOM, HET,           \
                                                                SVALUE>();          \
                });
}
#endif
