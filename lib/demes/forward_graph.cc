#include <cstdint>
#include <stdexcept>
#include <core/fp11rust.h>
#include <core/demes/forward_graph.hpp>

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
                    throw std::invalid_argument(message);
                }
        }
    };

    ForwardDemesGraph::~ForwardDemesGraph()
    {
    }

    ForwardDemesGraph::ForwardDemesGraph(const std::string &yaml, std::uint32_t burnin)
        : pimpl{new forward_graph_implementation(yaml, burnin)}
    {
    }
}
