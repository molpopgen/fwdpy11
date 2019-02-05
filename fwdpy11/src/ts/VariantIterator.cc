#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <fwdpp/ts/tree_visitor.hpp>

namespace py = pybind11;

class VariantIterator
{
  private:
    std::vector<fwdpp::ts::mutation_record>::const_iterator
    advance_trees_and_mutations()
    {
        ++mbeg;
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
                        return mbeg;
                    }
                ++mbeg;
            }
        return mbeg;
    }

  public:
    std::vector<fwdpp::ts::mutation_record>::const_iterator mbeg, mend;
    std::vector<double> pos;
    fwdpp::ts::tree_visitor tv;
    std::vector<std::int8_t> genotype_data;
    py::array_t<std::int8_t> genotypes;
    double current_position;
    VariantIterator(const fwdpp::ts::table_collection& tc,
                    const std::vector<fwdpy11::Mutation>& mutations,
                    const std::vector<fwdpp::ts::TS_NODE_INT>& samples)
        : mbeg(tc.mutation_table.begin()), mend(tc.mutation_table.end()),
          pos(), tv(tc, samples), genotype_data(samples.size(), 0),
          genotypes(fwdpy11::make_1d_ndarray(genotype_data)),
          current_position(std::numeric_limits<double>::quiet_NaN())
    {
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
            }
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
        if (ls != fwdpp::ts::TS_NULL_NODE)
            {
                current_position = pos[mbeg->key];
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
        .def(py::init<const fwdpp::ts::table_collection&,
                      const std::vector<fwdpy11::Mutation>&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&>(),
             py::keep_alive<1, 2>(), py::arg("tables"), py::arg("mutations"),
             py::arg("samples"),
             R"delim(
             :param tables: The table collection
             :type tables: :class:`fwdpy11.ts.TableCollection`
             :param mutations: Mutation container
             :type mutations: :class:`fwdpy11.VecMutation`
             :param samples: Samples list
             :type samples: list
            )delim")
        .def(py::init([](const fwdpy11::Population& pop,
                         const bool include_preserved) {
                 std::vector<fwdpp::ts::TS_NODE_INT> samples(2 * pop.N, 0);
                 std::iota(samples.begin(), samples.end(), 0);
                 if (include_preserved)
                     {
                         samples.insert(samples.end(),
                                        pop.tables.preserved_nodes.begin(),
                                        pop.tables.preserved_nodes.end());
                     }
                 return VariantIterator(pop.tables, pop.mutations, samples);
             }),
             py::arg("pop"), py::arg("include_preserved_nodes") = false)
        .def("__iter__",
             [](VariantIterator& v) -> VariantIterator& { return v; },
             py::keep_alive<0, 1>())
        .def("__next__", &VariantIterator::next_variant)
        .def_readonly("genotypes", &VariantIterator::genotypes,
                      "Genotype array.  Index order is same as sample input")
        .def_readonly("position", &VariantIterator::current_position,
                      "Current mutation position");
}

