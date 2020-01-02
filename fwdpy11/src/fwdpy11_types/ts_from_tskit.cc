#include <vector>
#include <utility>
#include <cstdint>
#include <cmath>
#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <fwdpp/ts/definitions.hpp>
#include <fwdpp/ts/edge.hpp>
#include <fwdpp/ts/node.hpp>
#include <fwdpy11/types/DiploidPopulation.hpp>

namespace py = pybind11;

// TODO put this in a header in case anyone else finds it useful?
std::tuple<fwdpp::ts::node_vector, fwdpp::ts::edge_vector, int, double>
convert_tables_from_tskit(py::object ts)
{
    py::object tstables = ts.attr("tables");
    py::object tsnodes = tstables.attr("nodes");
    py::object tsedges = tstables.attr("edges");

    auto left = tsedges.attr("left").cast<py::array_t<double>>();
    auto right = tsedges.attr("right").cast<py::array_t<double>>();
    auto parent = tsedges.attr("parent").cast<py::array_t<std::int32_t>>();
    auto child = tsedges.attr("child").cast<py::array_t<std::int32_t>>();

    // Use unchecked access for speed
    auto left_ = left.unchecked<1>();
    auto right_ = right.unchecked<1>();
    auto parent_ = parent.unchecked<1>();
    auto child_ = child.unchecked<1>();

    std::vector<fwdpp::ts::edge> edges;
    edges.reserve(left_.shape(0));
    for (py::ssize_t i = 0; i < left_.shape(0); ++i)
        {
            edges.push_back(
                fwdpp::ts::edge{ left_(i), right_(i), parent_(i), child_(i) });
        }

    // same treatment for nodes
    auto time = tsnodes.attr("time").cast<py::array_t<double>>();
    auto population
        = tsnodes.attr("population").cast<py::array_t<std::int32_t>>();
    auto time_ = time.unchecked<1>();
    auto population_ = population.unchecked<1>();

    std::vector<fwdpp::ts::node> nodes;
    nodes.reserve(time_.shape(0));
    int ntips = 0;
    for (py::ssize_t i = 0; i < time_.shape(0); ++i)
        {
            // reverse direction of time
            auto t = time_(i);
            if (t != 0.0)
                {
                    t *= -1.0;
                }
            else
                {
                    ++ntips;
                }
            nodes.push_back(fwdpp::ts::node{ population_(i), t });
        }

    double l = ts.attr("get_sequence_length")().cast<double>();
    return std::make_tuple(std::move(nodes), std::move(edges), ntips, l);
}

fwdpy11::DiploidPopulation
create_DiploidPopulation_from_tree_sequence(py::object ts)
{
    auto t = convert_tables_from_tskit(ts);
    auto nodes(std::move(std::get<0>(t)));
    auto edges(std::move(std::get<1>(t)));
    auto twoN = std::get<2>(t);
    auto l = std::get<3>(t);
    if (twoN % 2 != 0.0)
        {
            throw std::invalid_argument(
                "tree sequence has odd number of tips");
        }

    fwdpy11::DiploidPopulation pop(twoN / 2, l);
    pop.tables.node_table.swap(nodes);
    pop.tables.edge_table.swap(edges);
    // NOTE: this may be an issue when we allow variable survival probabilitites!
    pop.tables.update_offset();
    if (!pop.tables.edges_are_sorted())
        {
            throw std::runtime_error("edge table is not sorted");
        }
    pop.tables.build_indexes();

    // Fixed in 0.5.3: metadata nodes now correct
    // This was GitHub issue 333
    for (std::uint32_t i = 0; i < pop.N; ++i)
        {
            pop.diploid_metadata[i].nodes[0] = 2 * i;
            pop.diploid_metadata[i].nodes[1] = 2 * i + 1;
        }
    return pop;
}
