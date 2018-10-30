#ifndef FWDPY11_EVOLVETS_SAMPLE_RECORDER_HPP
#define FWDPY11_EVOLVETS_SAMPLE_RECORDER_HPP

#include <stdexcept>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <fwdpp/forward_types.hpp>

namespace fwdpy11
{
    struct samplerecorder
    {
        std::vector<fwdpp::uint_t> samples;

        samplerecorder() : samples{} {}

        void
        add_sample(const fwdpp::uint_t i)
        {
            samples.push_back(i);
        }

        void
        assign(pybind11::array_t<fwdpp::uint_t> a)
        {
            pybind11::buffer_info info = a.request();
            if (info.ndim != 1)
                {
                    throw std::invalid_argument(
                        "preserved sample list must have ndim == 1");
                }
            std::size_t size = info.shape[0];
            if (size == 0)
                {
                    throw std::invalid_argument(
                        "empty list of samples to preserve");
                }

            fwdpp::uint_t* ptr = static_cast<fwdpp::uint_t*>(info.ptr);
            samples.assign(ptr, ptr + size);
        }
    };
} // namespace fwdpy11
#endif
