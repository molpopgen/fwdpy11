#ifndef FWDPY11_FITNESS_HPP__
#define FWDPY11_FITNESS_HPP__

#include <functional>

namespace fwdpy11
{
    using singlepop_fitness = std::function<double(const fwdpy11::diploid_t &,
                                                   const fwdpy11::gcont_t &,
                                                   const fwdpy11::mcont_t &)>;

    template <typename fitness_model_type>
    struct fwdpp_singlepop_fitness_wrapper
    {
        using fitness_model = fitness_model_type;
        const double scaling;
        fwdpp_singlepop_fitness_wrapper(const double scaling_)
            : scaling(scaling_)
        {
        }
        inline singlepop_fitness
        callback() const
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
