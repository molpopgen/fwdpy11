#pragma once

#include <memory>
#include <cstdint>
#include <string>

namespace fwdpy11_core
{
    struct DemeSizeIterator
    {
        const double *first;
        const double *last;

        const double *
        begin() const
        {
            return first;
        }

        const double *
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
        ForwardDemesGraph(const std::string &yaml, std::uint32_t burnin);
        // The constructor has default behavior.
        // We must declare it as non-default
        // so that this class looks like a complete type
        // at compile time. This is a C++ "gotcha".
        ~ForwardDemesGraph();

        void initialize_model(std::uint32_t);
        bool iterating_model() const;
        void iterate_state();
        std::uint32_t model_end_time() const;
        std::ptrdiff_t number_of_demes() const;
        DemeSizeIterator parental_deme_sizes() const;
        DemeSizeIterator offspring_deme_sizes() const;
    };
}
