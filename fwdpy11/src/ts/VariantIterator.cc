#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <fwdpp/ts/tree_visitor.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Mutation>);

class VariantIterator
{
  private:
    using mut_table_itr
        = std::vector<fwdpp::ts::mutation_record>::const_iterator;

    mut_table_itr
    advance_trees_and_mutations()
    {
        while (mbeg < mend)
            {
                const auto& m = tv.tree();
                while (pos[mbeg->key] < m.left || pos[mbeg->key] >= m.right)
                    {
                        auto flag = tv(std::true_type(), std::true_type());
                        if (flag == false)
                            {
                                throw std::runtime_error(
                                    "VariantIterator: tree traversal "
                                    "error");
                            }
                    }
                if (m.leaf_counts[mbeg->node]
                        + m.preserved_leaf_counts[mbeg->node]
                    != 0)
                    {
                        if ((neutral[mbeg->key] && include_neutral)
                            || (!neutral[mbeg->key] && include_selected))
                            {
                                return mbeg;
                            }
                    }
                ++mbeg;
            }
        return mbeg;
    }

    mut_table_itr
    set_mbeg(mut_table_itr mtbeg, mut_table_itr mtend, const double start,
             const std::vector<fwdpy11::Mutation>& mutations)
    {
        if (std::isnan(start))
            {
                return mtbeg;
            }
        return std::lower_bound(
            mtbeg, mtend, start,
            [&mutations](const fwdpp::ts::mutation_record& mr,
                         const double v) {
                return mutations[mr.key].pos < v;
            });
    }

    mut_table_itr
    set_mend(mut_table_itr mtbeg, mut_table_itr mtend, const double end,
             const std::vector<fwdpy11::Mutation>& mutations)
    {
        if (std::isnan(end))
            {
                return mtend;
            }
        return std::upper_bound(
            mtbeg, mtend, end,
            [&mutations](const double v,
                         const fwdpp::ts::mutation_record& mr) {
                return v < mutations[mr.key].pos;
            });
    }
    std::vector<double> pos;
    std::vector<std::int8_t> neutral;
    const bool include_neutral, include_selected;

  public:
    std::vector<fwdpp::ts::mutation_record>::const_iterator mbeg, mend;
    fwdpp::ts::tree_visitor tv;
    std::vector<std::int8_t> genotype_data;
    py::array_t<std::int8_t> genotypes;
    double current_position;
    fwdpp::ts::mutation_record current_record;
    VariantIterator(const fwdpp::ts::table_collection& tc,
                    const std::vector<fwdpy11::Mutation>& mutations,
                    const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                    const double beg, const double end,
                    const bool include_neutral_variant,
                    const bool include_selected_variants)
        : pos(), neutral(), include_neutral(include_neutral_variant),
          include_selected(include_selected_variants),
          mbeg(set_mbeg(tc.mutation_table.begin(), tc.mutation_table.end(),
                        beg, mutations)),
          mend(set_mend(mbeg, tc.mutation_table.end(), end, mutations)),
          tv(tc, samples), genotype_data(samples.size(), 0),
          genotypes(fwdpy11::make_1d_ndarray(genotype_data)),
          current_position(std::numeric_limits<double>::quiet_NaN()),
          current_record{ fwdpp::ts::TS_NULL_NODE,
                          std::numeric_limits<std::size_t>::max() }
    {
        if (!include_selected && !include_neutral)
            {
                throw std::invalid_argument(
                    "excluding neutral and selected variants is invalid");
            }
        if (!std::isnan(beg) && !std::isnan(end))
            {
                if (!(end > beg))
                    {
                        throw std::invalid_argument(
                            "invalid position interval");
                    }
            }
        // Advance to first tree
        auto flag = tv(std::true_type(), std::true_type());
        if (flag == false)
            {
                throw std::invalid_argument(
                    "TableCollection contains no trees");
            }
        for (auto& m : mutations)
            {
                pos.push_back(m.pos);
                neutral.push_back(m.neutral);
            }
        mbeg = advance_trees_and_mutations();
    }

    VariantIterator&
    next_variant()
    {
        if (!(mbeg < mend))
            {
                throw py::stop_iteration();
            }
        const auto& m = tv.tree();
        std::fill(genotype_data.begin(), genotype_data.end(), 0);
        auto ls = m.left_sample[mbeg->node];
        current_position = std::numeric_limits<double>::quiet_NaN();
        current_record.key = std::numeric_limits<std::size_t>::max();
        current_record.node = fwdpp::ts::TS_NULL_NODE;
        if (ls != fwdpp::ts::TS_NULL_NODE)
            {
                current_position = pos[mbeg->key];
                current_record = *mbeg;
                auto rs = m.right_sample[mbeg->node];
                int nsteps = 1;
                while (true)
                    {
                        if (genotype_data[ls] == 1)
                            {
                                throw std::runtime_error(
                                    "VariantIterator error");
                            }
                        genotype_data[ls] = 1;
                        if (ls == rs)
                            {
                                break;
                            }
                        ls = m.next_sample[ls];
                        ++nsteps;
                    }
                if (nsteps != m.leaf_counts[mbeg->node])
                    {
                        throw std::runtime_error(
                            "VariantIterator: sample traversal error");
                    }
            }
        mbeg++;
        mbeg = advance_trees_and_mutations();
        return *this;
    }
};

void
init_variant_iterator(py::module& m)
{
    py::class_<VariantIterator>(
        m, "VariantIterator",
        "An iterable class for traversing genotypes in a tree sequence.")
        .def(py::init([](const fwdpp::ts::table_collection& tables,
                         const std::vector<fwdpy11::Mutation>& mutations,
                         const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                         double begin, double end,
                         bool include_neutral_variants,
                         bool include_selected_variants) {
                 return VariantIterator(tables, mutations, samples, begin, end,
                                        include_neutral_variants,
                                        include_selected_variants);
             }),
             py::arg("tables"), py::arg("mutations"), py::arg("samples"),
             py::arg("begin") = 0.0,
             py::arg("end") = std::numeric_limits<double>::max(),
             py::arg("include_neutral_variants") = true,
             py::arg("include_selected_variants") = true,
             R"delim(
             :param tables: The table collection
             :type tables: :class:`fwdpy11.TableCollection`
             :param mutations: Mutation container
             :type mutations: :class:`fwdpy11.MutationVector`
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
                 return VariantIterator(pop.tables, pop.mutations, samples,
                                        begin, end, include_neutral_variants,
                                        include_selected_variants);
             }),
             py::arg("pop"), py::arg("include_preserved_nodes") = false,
             py::arg("begin") = 0.0,
             py::arg("end") = std::numeric_limits<double>::max(),
             py::arg("include_selected_variants") = true,
             py::arg("include_selected_variants") = true,
             R"delim(
             :param pop: The table collection
             :type pop: :class:`fwdpy11.TableCollection`
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
        .def_readonly("record", &VariantIterator::current_record,
                      "Current mutation record")
        .def_readonly("position", &VariantIterator::current_position,
                      "Current mutation position");
}

