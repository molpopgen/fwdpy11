#pragma once

#include <cstdint>
#include <sstream>
#include <cmath>
#include <memory>
#include <stdexcept>
#include "demes_model_time.hpp"
#include "size_function.hpp"

namespace fwdpy11
{
    namespace discrete_demography
    {
        enum class SelfingPolicy
        {
            WrightFisher,
            Strict
        };

        struct Selfing
        {
            double rate;
            SelfingPolicy policy;

            Selfing(double rate, SelfingPolicy policy) : rate{rate}, policy{policy}
            {
                if (!std::isfinite(rate) || rate < 0. || rate > 1.0)
                    {
                        std::ostringstream message;
                        message << "invalid selfing rate: " << rate;
                        throw std::invalid_argument(message.str());
                    }
                if (policy == SelfingPolicy::WrightFisher && rate != 0.0)
                    {
                        throw std::invalid_argument(
                            "Wright-Fisher selfing must set a rate of 0.0");
                    }
            }

            static Selfing
            wright_fisher(double rate)
            {
                return Selfing(rate, SelfingPolicy::WrightFisher);
            }

            static Selfing
            strict(double rate)
            {
                return Selfing(rate, SelfingPolicy::Strict);
            }

            inline bool
            is_wright_fisher() const
            {
                return policy == SelfingPolicy::WrightFisher;
            }
        };

        // Selfing rate ideas:
        // NaN will be interpreted as W-F selfing.
        // Any other value 0. <= s <= 1. will be
        // taken as the probabilty that an individual selfs
        // or out-crosses.
        //
        // OR:
        //
        // We use some type of struct/enum to distinguish
        // W-F selfing from a "strict" selfing rate.
        //
        // Ideally, Epoch only stores the minimal info
        // and doesn't have to know about/interact with
        // the rng type, etc., used in fwdpy11.
        struct Epoch
        {
            demes_model_time end_time;
            std::uint32_t start_size, end_size;
            double cloning_rate;
            Selfing selfing;
            std::unique_ptr<SizeFunction> size_function;

            Epoch(demes_model_time end_time, std::uint32_t start_size,
                  std::uint32_t end_size, double cloning_rate, Selfing selfing,
                  std::unique_ptr<SizeFunction> size_function)
                : end_time{end_time}, start_size{start_size}, end_size{end_size},
                  cloning_rate{cloning_rate}, selfing{std::move(selfing)},
                  size_function{std::move(size_function)}
            {
                // NOTE: fail early!
                size_function->validate(start_size, end_size);
            }
        };
    }
}
