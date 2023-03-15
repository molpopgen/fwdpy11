#include <limits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <stdexcept>
#include <sstream>
#include <core/fp11rust.h>
#include <core/demes/forward_graph.hpp>
#include <fwdpy11/discrete_demography/exceptions.hpp>

namespace
{
    template <typename T>
    void
    throw_if_null(const T *t, std::string filename, int line)
    {
        if (t == nullptr)
            {
                std::ostringstream o;
                o << "unexpected NULL pointer: " << filename << ", " << line;
                throw std::runtime_error(o.str());
            }
    }
}

namespace fwdpy11_core
{

    using deleter = decltype(&demes_forward_graph_deallocate);
    using graph_ptr_t = std::unique_ptr<OpaqueForwardGraph, deleter>;

    struct ForwardDemesGraph::forward_graph_implementation
    {
        graph_ptr_t graph;
        std::ptrdiff_t number_of_demes;

        forward_graph_implementation(const std::string &yaml, std::uint32_t burnin,
                                     bool round)
            : graph(demes_forward_graph_allocate(), &demes_forward_graph_deallocate),
              number_of_demes{-1}
        {
            auto code = 0;
            if (round == true)
                {
                    code = demes_forward_graph_initialize_from_yaml_round_epoch_sizes(
                        yaml.c_str(), static_cast<double>(burnin), graph.get());
                }
            else
                {
                    code = demes_forward_graph_initialize_from_yaml(
                        yaml.c_str(), static_cast<double>(burnin), graph.get());
                }

            // FIXME: dup of handle_error_code
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
            if (demes_forward_graph_is_error_state(graph.get()))
                {
                    int status = 0;
                    auto message
                        = demes_forward_graph_get_error_message(graph.get(), &status);
                    if (message == nullptr)
                        {
                            throw std::runtime_error(
                                "graph in error state but message is nullptr");
                        }
                    throw std::invalid_argument(message);
                }
            number_of_demes = demes_forward_graph_number_of_demes(graph.get());
            if (number_of_demes < 1)
                {
                    throw std::invalid_argument("number of demes must be >= 1");
                }
        }

        void handle_error_code(std::int32_t code) const;
        void update_model_state_to_time(std::uint32_t);
        void initialize_time_iteration();
        const double *iterate();
        bool offspring_demes_exist() const;
        bool in_error_state() const;
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

    bool
    ForwardDemesGraph::forward_graph_implementation::in_error_state() const
    {
        auto rv = demes_forward_graph_is_error_state(graph.get());
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

    ForwardDemesGraph::ForwardDemesGraph(const std::string &yaml, std::uint32_t burnin,
                                         bool round)
        : pimpl{new forward_graph_implementation(yaml, burnin, round)}
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

    bool
    ForwardDemesGraph::in_error_state() const
    {
        return pimpl->in_error_state();
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

    std::ptrdiff_t
    ForwardDemesGraph::number_of_demes() const
    {
        return pimpl->number_of_demes;
    }
    ForwardDemesGraphDataIterator<double>
    ForwardDemesGraph::parental_deme_sizes() const
    {
        std::int32_t status;
        auto begin
            = demes_forward_graph_parental_deme_sizes(pimpl->graph.get(), &status);
        pimpl->handle_error_code(status);
        // NOTE: this might be overly strict,
        // but we keep it for now.
        throw_if_null(begin, __FILE__, __LINE__);
        return ForwardDemesGraphDataIterator<double>{begin, begin + number_of_demes()};
    }

    ForwardDemesGraphDataIterator<double>
    ForwardDemesGraph::offspring_deme_sizes() const
    {
        std::int32_t status;
        auto begin
            = demes_forward_graph_offspring_deme_sizes(pimpl->graph.get(), &status);
        pimpl->handle_error_code(status);
        // NOTE: this might be overly strict,
        // but we keep it for now.
        throw_if_null(begin, __FILE__, __LINE__);
        return ForwardDemesGraphDataIterator<double>{begin, begin + number_of_demes()};
    }

    ForwardDemesGraphDataIterator<double>
    ForwardDemesGraph::offspring_selfing_rates() const
    {
        std::int32_t status;
        auto begin = demes_forward_graph_selfing_rates(pimpl->graph.get(), &status);
        pimpl->handle_error_code(status);
        // NOTE: this might be overly strict,
        // but we keep it for now.
        throw_if_null(begin, __FILE__, __LINE__);
        return ForwardDemesGraphDataIterator<double>{begin, begin + number_of_demes()};
    }

    ForwardDemesGraphDataIterator<double>
    ForwardDemesGraph::offspring_cloning_rates() const
    {
        std::int32_t status;
        auto begin = demes_forward_graph_cloning_rates(pimpl->graph.get(), &status);
        pimpl->handle_error_code(status);
        // NOTE: this might be overly strict,
        // but we keep it for now.
        throw_if_null(begin, __FILE__, __LINE__);
        return ForwardDemesGraphDataIterator<double>{begin, begin + number_of_demes()};
    }

    ForwardDemesGraphDataIterator<double>
    ForwardDemesGraph::offspring_ancestry_proportions(std::size_t offspring_deme) const
    {
        std::int32_t status;
        auto begin = demes_forward_graph_ancestry_proportions(offspring_deme, &status,
                                                              pimpl->graph.get());
        pimpl->handle_error_code(status);
        // NOTE: this might be overly strict,
        // but we keep it for now.
        throw_if_null(begin, __FILE__, __LINE__);
        return ForwardDemesGraphDataIterator<double>{begin, begin + number_of_demes()};
    }

    std::uint32_t
    ForwardDemesGraph::sum_deme_sizes_at_time_zero()
    {
        std::int32_t status;
        auto rv
            = demes_forward_graph_sum_sizes_at_time_zero(&status, pimpl->graph.get());
        pimpl->handle_error_code(status);
        if (rv >= static_cast<double>(std::numeric_limits<std::uint32_t>::max()))
            {
                throw std::runtime_error(
                    "sum of sizes at time zero is too large for integer type");
            }
        return rv;
    }

    std::vector<std::uint32_t>
    ForwardDemesGraph::parental_deme_sizes_at_time_zero() const
    {
        pimpl->update_model_state_to_time(0.0);
        std::vector<std::uint32_t> rv;
        auto iter = this->parental_deme_sizes();
        for (auto i = std::begin(iter); i != std::end(iter); ++i)
            {
                rv.push_back(static_cast<std::uint32_t>(*i));
            }

        return rv;
    }
}
