#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <fwdpp/ts/site_visitor.hpp>
#include <fwdpp/ts/marginal_tree_functions/samples.hpp>

namespace py = pybind11;

class VariantIterator
{
  private:
    using site_table_itr = fwdpp::ts::site_vector::const_iterator;
    fwdpp::ts::site_visitor sv;
    site_table_itr scurrent, send;
    std::vector<std::int8_t> genotype_data;
    bool include_neutral, include_selected;
    double from, to;
    fwdpp::ts::convert_sample_index_to_nodes convert;
    py::list mutation_records;

  public:
    py::array_t<std::int8_t> genotypes;
    double current_position;
    VariantIterator(const fwdpp::ts::table_collection& tc,
                    const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                    const double beg, const double end,
                    const bool include_neutral_variant,
                    const bool include_selected_variants)
        : sv(tc, samples), scurrent(begin(tc.site_table)),
          send(std::end(tc.site_table)), genotype_data(samples.size(), 0),
          include_neutral(include_neutral_variant),
          include_selected(include_selected_variants), from(beg), to(end),
          convert(false), mutation_records{},
          genotypes(fwdpy11::make_1d_ndarray(genotype_data)),
          current_position(std::numeric_limits<double>::quiet_NaN())
    {
        if (!include_selected && !include_neutral)
            {
                throw std::invalid_argument(
                    "excluding neutral and selected variants is invalid");
            }
        if (!std::isnan(from) && !std::isnan(to))
            {
                if (!(to > from))
                    {
                        throw std::invalid_argument(
                            "invalid position interval");
                    }
            }
    }

    VariantIterator&
    next_variant()
    {
        if (!(scurrent < send) || scurrent->position >= to)
            {
                throw py::stop_iteration();
            }
        while ((scurrent = sv()) != end(sv) && scurrent->position < to)
            {
                if (scurrent->position >= from && scurrent->position < to)
                    {
                        mutation_records.attr("clear")();
                        current_position = scurrent->position;
                        std::fill(begin(genotype_data), end(genotype_data),
                                  scurrent->ancestral_state);
                        auto m = sv.get_mutations();
                        unsigned n = 0;
                        int neutral = -1, selected = -1;
                        for (auto i = m.first; i < m.second; ++i)
                            {
                                bool is_neutral = (i->neutral);
                                bool is_selected = (!is_neutral);
                                neutral += is_neutral;
                                selected += is_selected;
                                int ni = 0;
                                fwdpp::ts::process_samples(
                                    sv.current_tree(), convert, i->node,
                                    [this, i, &ni](fwdpp::ts::TS_NODE_INT u) {
                                        ++ni;
                                        genotype_data[u] = i->derived_state;
                                    });
                                if (ni)
                                    {
                                        fwdpp::ts::mutation_record temp(*i);
                                        auto o = py::cast(std::move(temp));
                                        mutation_records.append(std::move(o));
                                    }
                                n += ni;
                            }
                        if (neutral != -1 && selected != -1)
                            {
                                throw fwdpp::ts::tables_error(
                                    "invalid mutation data");
                            }
                        if (n)
                            {
                                if ((neutral > -1 && include_neutral)
                                    || include_selected)
                                    {
                                        return *this;
                                    }
                            }
                    }
            }
        if (!(scurrent < send) || scurrent->position >= to)
            {
                throw py::stop_iteration();
            }
        return *this;
    }

    fwdpp::ts::site
    current_site() const
    {
        if (!(scurrent < send))
            {
                throw std::runtime_error("end of sites");
            }
        return *scurrent;
    }

    py::list
    records() const
    {
        return mutation_records;
    }
};

void
init_variant_iterator(py::module& m)
{
    py::class_<VariantIterator>(
        m, "VariantIterator",
        "An iterable class for traversing genotypes in a tree sequence.")
        .def(py::init([](const fwdpp::ts::table_collection& tables,
                         const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                         double begin, double end,
                         bool include_neutral_variants,
                         bool include_selected_variants) {
                 return VariantIterator(tables, samples, begin, end,
                                        include_neutral_variants,
                                        include_selected_variants);
             }),
             py::arg("tables"), py::arg("samples"), py::arg("begin") = 0.0,
             py::arg("end") = std::numeric_limits<double>::max(),
             py::arg("include_neutral_variants") = true,
             py::arg("include_selected_variants") = true,
             R"delim(
             :param tables: The table collection
             :type tables: :class:`fwdpy11.TableCollection`
             :param samples: Samples list
             :type samples: list
             :param begin: (0.0) First position, inclusive.
             :param end: (max float) Last position, exclusive.
             :param include_neutral_variants: (True) Include neutral variants during traversal
             :type include_neutral_variants: boolean
             :param include_selected_variants: (True) Include selected variants during traversal
             :type include_selected_variants: boolean

             .. versionchanged:: 0.4.1
        
                 Add begin, end options as floats

            .. versionchanged:: 0.4.2

                 Add include_neutral_variants and include_selected_variants

            .. versionchanged:: 0.5.0

                 No longer requires a :class:`fwdpy11.MutationVector`.
            )delim")
        .def(py::init([](const fwdpy11::Population& pop,
                         const bool include_preserved, double begin,
                         double end, bool include_neutral_variants,
                         bool include_selected_variants) {
                 std::vector<fwdpp::ts::TS_NODE_INT> samples(2 * pop.N, 0);
                 std::iota(samples.begin(), samples.end(), 0);
                 if (include_preserved)
                     {
                         samples.insert(samples.end(),
                                        pop.tables.preserved_nodes.begin(),
                                        pop.tables.preserved_nodes.end());
                     }
                 return VariantIterator(pop.tables, samples, begin, end,
                                        include_neutral_variants,
                                        include_selected_variants);
             }),
             py::arg("pop"), py::arg("include_preserved_nodes") = false,
             py::arg("begin") = 0.0,
             py::arg("end") = std::numeric_limits<double>::max(),
             py::arg("include_selected_variants") = true,
             py::arg("include_selected_variants") = true,
             R"delim(
             :param pop: The population 
             :type pop: :class:`fwdpy11.DiploidPopulation`
             :param include_preserved_nodes: (False) Whether to include preserved samples during traversal
             :type include_preserved_nodes: boolean
             :param begin: (0.0) First position, inclusive.
             :param end: (max float) Last position, exclusive.
             :param include_neutral_variants: (True) Include neutral variants during traversal
             :type include_neutral_variants: boolean
             :param include_selected_variants: (True) Include selected variants during traversal
             :type include_selected_variants: boolean

             .. versionchanged:: 0.4.1
        
                 Add begin, end options as floats

            .. versionchanged:: 0.4.2

                 Add include_neutral_variants and include_selected_variants

            )delim")
        .def("__iter__",
             [](VariantIterator& v) -> VariantIterator& { return v; })
        .def("__next__", &VariantIterator::next_variant)
        .def_readonly("genotypes", &VariantIterator::genotypes,
                      "Genotype array.  Index order is same as sample input")
        .def_property_readonly("site", &VariantIterator::current_site,
                               "Current :class:`fwdpy11.Site`")
        .def_property_readonly(
            "records", &VariantIterator::records,
            "Returns a copy of the :class:`fwdpy11.MutationRecord` objects "
            "corresponding to the current site and sample")
        .def_readonly("position", &VariantIterator::current_position,
                      "Current mutation position");
}

