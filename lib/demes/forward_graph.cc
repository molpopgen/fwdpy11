#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <stdexcept>
#include <core/fp11rust.h>
#include <core/demes/forward_graph.hpp>
#include <fwdpy11/discrete_demography/exceptions.hpp>

namespace fwdpy11_core
{
    using deleter = decltype(&demes_forward_graph_deallocate);
    using graph_ptr_t = std::unique_ptr<OpaqueForwardGraph, deleter>;

    struct ForwardDemesGraph::forward_graph_implementation
    {
        graph_ptr_t graph;

        forward_graph_implementation(const std::string &yaml, std::uint32_t burnin)
            : graph(demes_forward_graph_allocate(), &demes_forward_graph_deallocate)
        {
            auto code = demes_forward_graph_initialize_from_yaml(
                yaml.c_str(), static_cast<double>(burnin), graph.get());
            handle_error_code(code);
        }

        void handle_error_code(std::int32_t code) const;
        void update_model_state_to_time(std::uint32_t);
        void initialize_time_iteration();
        const double *iterate();
        bool offspring_demes_exist() const;
        void update_internal_state(double time);
    };

    void
    ForwardDemesGraph::forward_graph_implementation::handle_error_code(
        std::int32_t code) const
    {
        if (code < 0)
            {
                // NOTE: the rust lib currently doesn't
                // set the status code.  Rather, it returns
                // nullptr if the model is not, in fact,
                // in an error state.
                std::int32_t status;
                auto message
                    = demes_forward_graph_get_error_message(graph.get(), &status);
                if (message == nullptr)
                    {
                        throw std::runtime_error(
                            "graph in error state but message is nullptr");
                    }
                throw fwdpy11::discrete_demography::DemographyError(message);
            }
    }

    void
    ForwardDemesGraph::forward_graph_implementation::update_model_state_to_time(
        std::uint32_t time)
    {
        auto code
            = demes_forward_graph_update_state(static_cast<double>(time), graph.get());
        handle_error_code(code);
    }

    void
    ForwardDemesGraph::forward_graph_implementation::initialize_time_iteration()
    {
        auto code = demes_forward_graph_initialize_time_iteration(graph.get());
        handle_error_code(code);
    }

    const double *
    ForwardDemesGraph::forward_graph_implementation::iterate()
    {
        std::int32_t code;
        auto rv = demes_forward_graph_iterate_time(graph.get(), &code);
        handle_error_code(code);
        return rv;
    }

    bool
    ForwardDemesGraph::forward_graph_implementation::offspring_demes_exist() const
    {
        std::int32_t code;
        auto rv = demes_forward_graph_any_extant_offspring_demes(graph.get(), &code);
        handle_error_code(code);
        return rv;
    }

    void
    ForwardDemesGraph::forward_graph_implementation::update_internal_state(double time)
    {
        std::int32_t code = demes_forward_graph_update_state(time, graph.get());
        handle_error_code(code);
    }

    ForwardDemesGraph::~ForwardDemesGraph()
    {
    }

    ForwardDemesGraph::ForwardDemesGraph(const std::string &yaml, std::uint32_t burnin)
        : pimpl{new forward_graph_implementation(yaml, burnin)}
    {
    }

    std::uint32_t
    ForwardDemesGraph::model_end_time() const
    {
        std::int32_t code;
        auto end_time = demes_forward_graph_model_end_time(&code, pimpl->graph.get());
        pimpl->handle_error_code(code);
        return static_cast<std::uint32_t>(end_time);
    }

    void
    ForwardDemesGraph::initialize_model(std::uint32_t time)
    {
        pimpl->update_model_state_to_time(time);
        pimpl->initialize_time_iteration();
        pimpl->iterate();
    }

    bool
    ForwardDemesGraph::iterating_model() const
    {
        return pimpl->offspring_demes_exist();
    };

    void
    ForwardDemesGraph::iterate_state()
    {
        auto t = pimpl->iterate();
        if (t != nullptr)
            {
                pimpl->update_internal_state(*t);
            }
    };
}
