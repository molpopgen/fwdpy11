#ifndef FWDPY11_CREATE_POPS_HPP__
#define FWDPY11_CREATE_POPS_HPP__

#include <limits>
#include <vector>
#include <utility>

namespace fwdpy11
{
    template <typename poptype> struct create_wrapper
    {
        template <typename diploids_input, typename gametes_input,
                  typename mutations_input>
        inline poptype
        operator()(diploids_input &&diploids, gametes_input &&gametes,
                   mutations_input &&mutations) const
        {
            return poptype(std::forward<diploids_input>(diploids),
                           std::forward<gametes_input>(gametes),
                           std::forward<mutations_input>(mutations));
        }

        template <typename diploids_input, typename gametes_input,
                  typename mutations_input>
        inline poptype
        operator()(
            diploids_input &&diploids, gametes_input &&gametes,
            mutations_input &&mutations,
            std::vector<std::pair<double, double>> &&locus_boundaries) const
        {
            return poptype(std::forward<diploids_input>(diploids),
                           std::forward<gametes_input>(gametes),
                           std::forward<mutations_input>(mutations),
                           std::move(locus_boundaries));
        }

        template <typename diploids_input, typename gametes_input,
                  typename mutations_input>
        inline poptype
        operator()(diploids_input &&diploids, gametes_input &&gametes,
                   mutations_input &&mutations, mutations_input &&fixations,
                   std::vector<fwdpp::uint_t> &&fixation_times,
                   fwdpp::uint_t generation) const
        {
            auto rv
                = this->operator()(std::forward<diploids_input>(diploids),
                                   std::forward<gametes_input>(gametes),
                                   std::forward<mutations_input>(mutations));

            rv.fixations.swap(fixations);
            rv.fixation_times.swap(fixation_times);
            rv.generation = generation;
            return rv;
            return rv;
        }
        template <typename diploids_input, typename gametes_input,
                  typename mutations_input>
        inline poptype
        operator()(diploids_input &&diploids, gametes_input &&gametes,
                   mutations_input &&mutations,
                   std::vector<std::pair<double, double>> &&locus_boundaries,
                   mutations_input &&fixations,
                   std::vector<fwdpp::uint_t> &&fixation_times,
                   fwdpp::uint_t generation) const
        {
            auto rv
                = this->operator()(std::forward<diploids_input>(diploids),
                                   std::forward<gametes_input>(gametes),
                                   std::forward<mutations_input>(mutations),
                                   std::move(locus_boundaries));

            rv.fixations.swap(fixations);
            rv.fixation_times.swap(fixation_times);
            rv.generation = generation;
            return rv;
            return rv;
        }
    };
}

#endif
