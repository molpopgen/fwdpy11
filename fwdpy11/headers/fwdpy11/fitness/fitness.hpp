#ifndef FWDPY11_FITNESS_HPP__
#define FWDPY11_FITNESS_HPP__

#include <functional>

namespace fwdpy11
{
    template <typename... T>
    using singlepop_fitness_signature
        = std::function<double(const fwdpy11::diploid_t &,
                               const fwdpy11::gcont_t &,
                               const fwdpy11::mcont_t &, T...)>;

    /*! Single-deme fitness function signature for standard "popgen"
     *  simulations.
     */
    using singlepop_fitness_fxn = singlepop_fitness_signature<>;
    //! For quantitative traits, optimum and VS are additional arguments
    using singlepop_fitness_fxn_qtrait
        = singlepop_fitness_signature<double, double>;

    struct singlepop_fitness
    //! Pure virtual base class for single-deme fitness functions
    {
        virtual ~singlepop_fitness() = default;
        singlepop_fitness() = default;
        virtual singlepop_fitness_fxn callback() const = 0;
    };

    struct singlepop_fitness_qtrait
        /*! Pure virtual base class for single-deme fitness functions
         *  for simulations of gaussian stabilizing selection.
         */
    {
        virtual ~singlepop_fitness_qtrait() = default;
        singlepop_fitness_qtrait() = default;
        virtual singlepop_fitness_fxn_qtrait callback() const = 0;
    };

    template <typename fitness_model_type>
    struct fwdpp_singlepop_fitness_wrapper : public singlepop_fitness
    {
        using fitness_model = fitness_model_type;
        const double scaling;
        fwdpp_singlepop_fitness_wrapper(const double scaling_)
            : scaling(scaling_)
        {
        }
        inline singlepop_fitness_fxn
        callback() const final
        {
            return std::bind(fitness_model(), std::placeholders::_1,
                             std::placeholders::_2, std::placeholders::_3,
                             scaling);
        }
    };

    using singlepop_mult_wrapper
        = fwdpp_singlepop_fitness_wrapper<KTfwd::multiplicative_diploid>;
    using singlepop_additive_wrapper
        = fwdpp_singlepop_fitness_wrapper<KTfwd::additive_diploid>;
}
#endif
