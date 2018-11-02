#include <vector>
#include <utility>
#include <cstdint>
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <fwdpp/ts/definitions.hpp>
#include <fwdpp/ts/edge.hpp>
#include <fwdpp/ts/node.hpp>
#include <fwdpy11/types/SlocusPop.hpp>

namespace py = pybind11;

py::tuple
convert_tables_from_msprime(py::object ts, const bool discretize_time)
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
            if (discretize_time)
                {
                    t = std::ceil(t);
                }
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
    return py::make_tuple(std::move(nodes), std::move(edges), ntips, l);
}

fwdpy11::SlocusPop
create_SlocusPop_from_tree_sequence(py::object ts, const bool discretize_time)
{
    auto t = convert_tables_from_msprime(ts, discretize_time);
    auto nodes = t[0].cast<std::vector<fwdpp::ts::node>>();
    auto edges = t[1].cast<std::vector<fwdpp::ts::edge>>();
    auto twoN = t[2].cast<int>();
    if (twoN % 2 != 0.0)
        {
            throw std::invalid_argument(
                "tree sequence has odd number of tips");
        }
    auto l = t[3].cast<double>();
    t.release();

    fwdpy11::SlocusPop pop(twoN / 2, l);
    pop.tables.node_table.swap(nodes);
    pop.tables.edge_table.swap(edges);
    pop.tables.update_offset();
    if (!pop.tables.edges_are_sorted())
        {
            throw std::runtime_error("edge table is not sorted");
        }
    pop.tables.build_indexes();
    return pop;
}

PYBIND11_MODULE(ts_from_msprime, m)
{
    m.doc() = "Converts node and edge data from tree sequences generated in "
              "msprime";

    // Import fwdpy11.ts so that relevant C++ types
    // have their Python wrappers visible here
    auto imported_ts = static_cast<pybind11::object>(
        pybind11::module::import("fwdpy11.ts"));

    m.def("convert_tables", &convert_tables_from_msprime, py::arg("ts"),
          py::arg("discretize_time"),
          R"delim(
        Get node and edge tables from msprime.

        This function primarily exists to help construct populations
        from the output of a coalescent simulation.  You probably won't call
        it directly.

        :param ts: A tree sequence from msprime
        :type ts: msprime.TreeSequence
        :param discretize_time: Whether to convert time into integer values
        :type discretize_time: boolean

        :rtype: tuple
        :returns: :class:`fwdpy11.ts.NodeTable`, :class:`fwdpy11.ts.EdgeTable`, sequence length

        The node table has had time values reversed so that tip times are zero
        and ancestral node times are negative.

        When `discretize_time` is `True`, all time values are converted to their floors.

        .. versionadded:: 0.2.0

        .. note::

            It is not always straightforward to start a forward simulation
            with a valid tree from a coalescent simulation! Pay very careful
            attention to parameter scaling.  Also, be very careful when 
            mixing trees from coalescent simulations with forward simulations 
            that are not doing Wright-Fisher dynamics.
        )delim");

    m.def("_create_SlocusPop", &create_SlocusPop_from_tree_sequence);
}
