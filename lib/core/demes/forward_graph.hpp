#pragma once

#include <memory>
#include <cstdint>
#include <string>
#include <vector>

namespace fwdpy11_core
{
    template <typename T> struct ForwardDemesGraphDataIterator
    {
        const T *first;
        const T *last;

        const T *
        begin() const
        {
            return first;
        }

        const T *
        end() const
        {
            return last;
        }
    };

    class ForwardDemesGraph
    {
      private:
        struct forward_graph_implementation;
        std::unique_ptr<forward_graph_implementation> pimpl;

      public:
        ForwardDemesGraph(const std::string &yaml, std::uint32_t burnin,
                          bool round = false);
        // The constructor has default behavior.
        // We must declare it as non-default
        // so that this class looks like a complete type
        // at compile time. This is a C++ "gotcha".
        ~ForwardDemesGraph();

        void initialize_model(std::uint32_t);
        bool iterating_model() const;
        bool in_error_state() const;
        void iterate_state();
        std::uint32_t model_end_time() const;
        std::ptrdiff_t number_of_demes() const;
        ForwardDemesGraphDataIterator<double> parental_deme_sizes() const;
        ForwardDemesGraphDataIterator<double> offspring_deme_sizes() const;
        ForwardDemesGraphDataIterator<double> offspring_selfing_rates() const;
        ForwardDemesGraphDataIterator<double> offspring_cloning_rates() const;
        ForwardDemesGraphDataIterator<double>
        offspring_ancestry_proportions(std::size_t offspring_deme) const;

        // TODO: this can become const once upstream
        // gets a fn to handle this
        std::uint32_t sum_deme_sizes_at_time_zero();
        std::vector<std::uint32_t> parental_deme_sizes_at_time_zero() const;
    };
}
