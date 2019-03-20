#ifndef FWDPY11_GET_INDIVIDUALS_HPP__
#define FWDPY11_GET_INDIVIDUALS_HPP__

#include <cstdint>
#include <stdexcept>
#include <vector>
#include <fwdpy11/rng.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


std::vector<std::size_t>
get_individuals(const fwdpp::uint_t popsize, pybind11::kwargs kwargs)
{
    std::vector<std::size_t> ind;
    bool has_ind = kwargs.contains("individuals");
    bool has_nsam = kwargs.contains("nsam");
    bool has_rng = kwargs.contains("rng");

    if (has_ind && !(has_nsam || has_rng))
        {
            ind = kwargs["individuals"].cast<decltype(ind)>();
        }
    else if (has_rng && has_nsam && !has_ind)
        {
            const auto& rng = kwargs["rng"].cast<const fwdpy11::GSLrng_t&>();
            const auto nsam = kwargs["nsam"].cast<fwdpp::uint_t>();
            for (fwdpp::uint_t i = 0; i < nsam; ++i)
                {
                    ind.push_back(static_cast<std::size_t>(
                        gsl_rng_uniform_int(rng.get(), popsize)));
                }
        }
    else
        {
            throw std::invalid_argument("invalid kwargs");
        }

    return ind;
}

#endif
