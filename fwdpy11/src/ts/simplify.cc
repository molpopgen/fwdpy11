#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpp/ts/table_simplifier.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Mutation>);

py::tuple
simplify(const fwdpy11::Population& pop,
         const std::vector<fwdpp::ts::TS_NODE_INT>& samples)
{
    if (pop.tables.genome_length() == std::numeric_limits<double>::max())
        {
            throw std::invalid_argument(
                "population is not using tree sequences");
        }
    if (pop.tables.num_nodes() == 0)
        {
            throw std::invalid_argument(
                "population has empty TableCollection");
        }
    if (samples.empty())
        {
            throw std::invalid_argument("empty sample list");
        }
    if (std::any_of(samples.begin(), samples.end(),
                    [&pop](const fwdpp::ts::TS_NODE_INT s) {
                        return s == fwdpp::ts::TS_NULL_NODE
                               || static_cast<std::size_t>(s)
                                      >= pop.tables.num_nodes();
                    }))
        {
            throw std::invalid_argument("invalid sample list");
        }
    auto t(pop.tables);
    // NOTE: If the user wants to keep ancient samples,
    // they must pass them in to the samples list
    t.preserved_nodes.clear();
    fwdpp::ts::table_simplifier simplifier(pop.tables.genome_length());
    auto rv = simplifier.simplify(t, samples);
    t.build_indexes();
    return py::make_tuple(std::move(t), fwdpy11::make_1d_array_with_capsule(
                                            std::move(rv.first)));
}

void
init_simplify_functions(py::module& m)
{
    m.def("simplify", &simplify, py::arg("pop"), py::arg("samples"),
          R"delim(
            Simplify a TableCollection stored in a Population.

            :param pop: A :class:`fwdpy11.PopulationBase`
            :param samples: A list of samples (node indexes).
                
            :return: The simplified tables and array mapping input sample IDs to output IDS

            :rtype: tuple

            Note that the samples argument is agnostic with respect to the time of
            the nodes in the input tables. Thus, you may do things like simplify
            to a set of "currently-alive" nodes plus some or all ancient samples by
            including some node IDs from :attr:`fwdpy11.TableCollection.preserved_nodes`.
            
            If the input contains ancient samples, and you wish to include them in the output,
            then you need to include their IDs in the samples argument.

            .. note::

                Due to node ID remapping, the metadata corresponding to nodes becomes a bit more
                difficult to look up.  You need to use the output ID map, the original IDs, and 
                the population's metadata containers.

            .. deprecated:: 0.3.0

                Prefer :func:`fwdpp.simplify_tables`

            
            .. versionchanged:: 0.3.0

                Ancient samples are no longer kept by default

            .. versionchanged:: 0.5.0

                No longer requires a :class:`MutationVector` argument.

            )delim");

    m.def(
        "simplify_tables",
        [](const fwdpp::ts::table_collection& tables,
           const std::vector<fwdpp::ts::TS_NODE_INT>& samples) -> py::tuple {
            auto t(tables);
            t.preserved_nodes.clear();
            fwdpp::ts::table_simplifier simplifier(tables.genome_length());
            auto rv = simplifier.simplify(t, samples);
            t.build_indexes();
            return py::make_tuple(std::move(t),fwdpy11::make_1d_array_with_capsule(std::move(rv.first)));
        },
        py::arg("tables"), py::arg("samples"),
        R"delim(
          Simplify a TableCollection.
          
          :param pop: A table collection.
          :type pop: :class:`fwdpy11.TableCollection`
          :param samples: list of samples
          :type list: list-like or array-like

          :returns: A simplified TableCollection and an array containing remapped sample ids.
          :rtype: tuple

          .. versionadded:: 0.3.0
          )delim");
}

