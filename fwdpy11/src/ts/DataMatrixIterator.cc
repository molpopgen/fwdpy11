#include <memory>
#include <cmath>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/detail/generate_data_matrix_details.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Mutation>);

class DataMatrixIterator
// Encapsulate fwdpp::ts::tree_visitor and fwdpp::data_matrix
// to provide very fast interation over multiple genomic intervals.
// We accomplish this by copying the state of the current tree_visitor
// whenever we recognize that the current tree overlaps with the next interval.
// The saved state can be restored in O(1) time, "snapping" us to the start of
// next interval without having to initiate a new tree_visitor.
{
  private:
    using mut_table_itr
        = std::vector<fwdpp::ts::mutation_record>::const_iterator;
    std::unique_ptr<fwdpp::ts::tree_visitor> current_tree, next_tree;
    const std::vector<std::pair<double, double>> position_ranges;
    std::vector<std::int8_t> genotypes, is_neutral;
    std::vector<double> mutation_positions;
    std::unique_ptr<fwdpp::data_matrix> dmatrix;
    const mut_table_itr mbeg, mend;
    mut_table_itr mcurrent;
    std::size_t current_range;
    const bool include_neutral_variants, include_selected_variants,
        include_fixations;
    bool matrix_requires_clearing;

    std::unique_ptr<fwdpp::ts::tree_visitor>
    initialize_current_tree(const fwdpp::ts::table_collection& tables,
                            const std::vector<fwdpp::ts::TS_NODE_INT>& samples)
    {
        std::unique_ptr<fwdpp::ts::tree_visitor> rv(
            new fwdpp::ts::tree_visitor(tables, samples));
        auto flag = rv->operator()(std::true_type(), std::true_type());
        if (flag == false)
            {
                throw std::invalid_argument(
                    "TableCollection contains no trees");
            }
        return rv;
    }

    std::vector<std::pair<double, double>>
    init_intervals(
        const std::vector<std::pair<double, double>>& input_intervals)
    {
        if (input_intervals.empty())
            {
                throw std::invalid_argument("empty interval list");
            }
        for (auto& i : input_intervals)
            {
                if (!std::isfinite(i.first) || !std::isfinite(i.second))
                    {
                        throw std::invalid_argument(
                            "invalid interval: all values must be finite");
                    }
                if (i.second < 0.0 || i.first < 0.0)
                    {
                        throw std::invalid_argument(
                            "invalid interval: all positions must be >= "
                            "0.0");
                    }
                if (!(i.second > i.first))
                    {
                        throw std::invalid_argument(
                            "invalid interval: end <= beg");
                    }
            }
        for (std::size_t i = 1; i < input_intervals.size(); ++i)
            {
                if (!(input_intervals[i].first > input_intervals[i - 1].first))
                    {
                        throw std::invalid_argument(
                            "invalid interval start positions");
                    }
            }
        return input_intervals;
    }

    std::vector<std::int8_t>
    set_neutral(const std::vector<fwdpy11::Mutation>& mutations)
    {
        std::vector<std::int8_t> n;
        n.reserve(mutations.size());
        for (auto& m : mutations)
            {
                n.push_back(m.neutral);
            }
        return n;
    }

    std::vector<double>
    set_positions(const std::vector<fwdpy11::Mutation>& mutations)
    {
        std::vector<double> p;
        p.reserve(mutations.size());
        for (auto& m : mutations)
            {
                p.push_back(m.pos);
            }
        return p;
    }

    mut_table_itr
    find_first_mutation_record(mut_table_itr b, mut_table_itr e,
                               const double start)
    {
        return std::lower_bound(
            b, e, start,
            [this](const fwdpp::ts::mutation_record& mr, const double v) {
                return mutation_positions[mr.key] < v;
            });
    }

    mut_table_itr
    init_trees_and_mutations()
    {
        if (current_tree == nullptr || position_ranges.empty())
            {
                throw std::runtime_error(
                    "DataMatrixIterator __init__ failure");
            }
        double first_left = position_ranges.front().first;
        double current_right = current_tree->tree().right;

        while (current_right < first_left)
            {
                current_tree->operator()(std::true_type(), std::true_type());
                current_right = current_tree->tree().right;
            }
        double current_tree_left = current_tree->tree().left;
        mut_table_itr firstmut = mbeg;
        while (firstmut < mend
               && mutation_positions[firstmut->key] < current_tree_left)
            {
                ++firstmut;
            }
        return firstmut;
    }

    mut_table_itr
    find_first_mutation_record_after_current_max()
    // After a call to cleanup_matrix, we may need to reset mcurrent
    // to the first mutation AFTER dmatrix's current data.  We do so
    // here.
    {
        if (dmatrix == nullptr)
            {
                throw std::runtime_error("DataMatrix is nullptr");
            }
        double m = std::numeric_limits<double>::max();
        if (!dmatrix->neutral.positions.empty())
            {
                m = dmatrix->neutral.positions.back();
            }
        if (!dmatrix->selected.positions.empty())
            {
                if (m == std::numeric_limits<double>::max())
                    {
                        m = dmatrix->selected.positions.back();
                    }
                else
                    {
                        m = std::max(m, dmatrix->selected.positions.back());
                    }
            }
        if (m == std::numeric_limits<double>::max())
            {
                return mcurrent;
            }
        return std::upper_bound(
            mcurrent, mend, m,
            [this](double v, const fwdpp::ts::mutation_record& mr) {
                return v < mutation_positions[mr.key];
            });
    }

    void
    advance_trees()
    // Advance the tree visitor to the left edge of current interval.
    {
        if (current_range >= position_ranges.size())
            {
                throw std::runtime_error("DataMatrixIterator fatal error");
            }
        auto l = position_ranges[current_range].first;
        while (current_tree->tree().right < l)
            {
                auto flag = current_tree->operator()(std::true_type(),
                                                     std::true_type());
                if (!flag)
                    {
                        break;
                    }
            }
    }

    mut_table_itr
    advance_mutations()
    {
        matrix_requires_clearing = false;
        if (next_tree != nullptr)
            {
                current_tree.swap(next_tree);
                next_tree.reset(nullptr);
                double left = position_ranges[current_range].first;
                double right = position_ranges[current_range].second;
                mcurrent = find_first_mutation_record(mbeg, mend, left);
                cleanup_matrix(left, right);
                mcurrent = find_first_mutation_record_after_current_max();
            }
        else
            {
                matrix_requires_clearing = true;
                const auto& m = current_tree->tree();
                while (mcurrent < mend
                       && mutation_positions[mcurrent->key] < m.left)
                    {
                        ++mcurrent;
                    }
            }
        return mcurrent;
    }

    void
    release_memory()
    // Called when iteration stops.
    // Frees potentially-large data
    // structures on the C++ side
    // in case the Python object isn't
    // GC'd anytime soon.
    {
        current_tree.reset(nullptr);
        next_tree.reset(nullptr);
        dmatrix.reset(nullptr);
    }

    void
    clear_matrix()
    // Clear out member data of dmatrix,
    // but keep the memory allocated for
    // reuse
    {
        if (dmatrix == nullptr)
            {
                throw std::runtime_error("DataMatrix is nullptr");
            }
        dmatrix->neutral_keys.clear();
        dmatrix->selected_keys.clear();
        dmatrix->neutral.data.clear();
        dmatrix->neutral.positions.clear();
        dmatrix->selected.data.clear();
        dmatrix->selected.positions.clear();
    }

    void
    cleanup_matrix_details(fwdpp::state_matrix& sm,
                           std::vector<std::size_t>& keys, double l, double r)
    // When genomic intervals overlap, they have mutations in common.
    // This function removes all mutations from the previous window,
    // keeping any mutations shared by both windows.
    {
        // find first key corresponding to position >= l
        auto itr = std::lower_bound(begin(keys), end(keys), l,
                                    [this](std::size_t k, double v) {
                                        return mutation_positions[k] < v;
                                    });
        // This is the number of mutations
        // with positions < l
        auto d = std::distance(begin(keys), itr);
        if (d > 0)
            {
                // Sanity check
                auto pitr = std::lower_bound(begin(sm.positions),
                                             end(sm.positions), l);
                if (d != std::distance(begin(sm.positions), pitr))
                    {
                        throw std::runtime_error(
                            "DataMatrix internal state inconsistent");
                    }
                // erase all keys where position < p...
                keys.erase(begin(keys), itr);

                // ...and positions...
                sm.positions.erase(begin(sm.positions), pitr);

                // ...and genotypes.
                sm.data.erase(begin(sm.data),
                              begin(sm.data) + d * dmatrix->ncol);
            }

        // Give all mutations >= r the same treatment
        itr = std::lower_bound(begin(keys), end(keys), r,
                               [this](std::size_t k, double v) {
                                   return mutation_positions[k] < v;
                               });
        d = std::distance(itr, end(keys));
        if (d > 0)
            {
                auto pitr = std::lower_bound(begin(sm.positions),
                                             end(sm.positions), r);
                if (d != std::distance(pitr, end(sm.positions)))
                    {
                        throw std::runtime_error(
                            "DataMatrix internal state inconsistent");
                    }
                keys.erase(itr, end(keys));
                sm.positions.erase(pitr, end(sm.positions));
                sm.data.erase(end(sm.data) - (d * dmatrix->ncol),
                              end(sm.data));
            }
    }

    void
    cleanup_matrix(double l, double r)
    // Removes all data in dmatrix
    // corresponding to position < p
    {
        cleanup_matrix_details(dmatrix->neutral, dmatrix->neutral_keys, l, r);
        cleanup_matrix_details(dmatrix->selected, dmatrix->selected_keys, l,
                               r);
    }

    void
    check_if_still_iterating()
    {
        if (!(mcurrent < mend) || current_range >= position_ranges.size())
            {
                release_memory();
                throw py::stop_iteration();
            }
    }

    void
    update_data_matrix(const bool mut_is_neutral, const std::size_t key)
    // NOTE: this is a re-implementation
    // of fwdpp::ts::detail::update_data_matrix
    {
        auto& sm = (mut_is_neutral) ? dmatrix->neutral : dmatrix->selected;
        auto& k = (mut_is_neutral) ? dmatrix->neutral_keys
                                   : dmatrix->selected_keys;
        k.push_back(key);
        sm.positions.push_back(mutation_positions[key]);
        sm.data.insert(sm.data.end(), begin(genotypes), end(genotypes));
    }

    bool
    mutation_in_current_range(const mut_table_itr mitr)
    {
        if (!(mitr < mend))
            {
                return false;
            }
        if (current_range >= position_ranges.size())
            {
                throw std::runtime_error("DataMatrixIterator fatal error");
            }
        auto pos = mutation_positions[mitr->key];
        return pos >= position_ranges[current_range].first
               && pos < position_ranges[current_range].second;
    }

    void
    process_current_mutation(const fwdpp::ts::marginal_tree& tree,
                             const mut_table_itr mitr)
    {
        // Make sure we skip mutants on the
        // current tree but not in the current
        // genomic interval.
        if (!mutation_in_current_range(mitr))
            {
                return;
            }
        auto lc = tree.leaf_counts[mitr->node];
        bool fixed = (lc == tree.sample_size);
        if (lc > 0 && (!fixed || (fixed && include_fixations)))
            {
                bool mut_is_neutral = is_neutral[mitr->key];
                bool tracking_n = (mut_is_neutral && include_neutral_variants);
                bool tracking_s
                    = (!mut_is_neutral && include_selected_variants);
                if (tracking_n || tracking_s)
                    {
                        auto index = tree.left_sample[mitr->node];
                        if (index != fwdpp::ts::TS_NULL_NODE)
                            {
                                fwdpp::ts::detail::process_samples(
                                    tree, mitr->node, index, genotypes);
                                update_data_matrix(mut_is_neutral,
                                                   mcurrent->key);
                            }
                    }
            }
    }

    py::array_t<std::int8_t>
    genotype_matrix_wrapper(const fwdpp::state_matrix& sm,
                            const std::vector<std::size_t>& keys) const
    // Implementation details behing public getter for genotype matrix
    {
        std::size_t ncol = 0;
        if (!keys.empty())
            {
                ncol = sm.data.size() / keys.size();
            }
        return fwdpy11::make_2d_ndarray_readonly(sm.data, keys.size(), ncol);
    }

  public:
    DataMatrixIterator(const fwdpp::ts::table_collection& tables,
                       const std::vector<fwdpy11::Mutation>& mutations,
                       const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                       const std::vector<std::pair<double, double>>& intervals,
                       bool neutral, bool selected, bool fixations)
        : current_tree(initialize_current_tree(tables, samples)),
          next_tree(nullptr), position_ranges(init_intervals(intervals)),
          genotypes(samples.size(), 0), is_neutral(set_neutral(mutations)),
          mutation_positions(set_positions(mutations)),
          dmatrix(new fwdpp::data_matrix(samples.size())),
          mbeg(tables.mutation_table.begin()),
          mend(tables.mutation_table.end()),
          mcurrent(init_trees_and_mutations()), current_range(0),
          include_neutral_variants(neutral),
          include_selected_variants(selected), include_fixations(fixations),
          matrix_requires_clearing(false)
    {
    }

    DataMatrixIterator&
    next_data_matrix()
    // This is the back end for the Python
    // class's __next__ function
    {
        check_if_still_iterating();
        advance_trees();
        mcurrent = advance_mutations();
        next_tree.reset(nullptr);
        if (matrix_requires_clearing)
            {
                clear_matrix();
            }

        bool iteration_flag = true;
        bool tree_in_current_range = true;
        do
            {
                const auto& tree = current_tree->tree();
                if (next_tree == nullptr
                    && current_range + 1 < position_ranges.size())
                    {
                        double right = tree.right;
                        if (right > position_ranges[current_range + 1].first)
                            {
                                next_tree.reset(new fwdpp::ts::tree_visitor(
                                    *current_tree));
                            }
                    }
                for (; mcurrent < mend
                       && mutation_positions[mcurrent->key] < tree.right;
                     ++mcurrent)
                    {
                        process_current_mutation(tree, mcurrent);
                    }
                iteration_flag = current_tree->operator()(std::true_type(),
                                                          std::true_type());
                tree_in_current_range
                    = (current_tree->tree().left
                       < position_ranges[current_range].second);
            }
        while (iteration_flag == true && tree_in_current_range);
        ++current_range;
        return *this;
    }

    // Getter functions to establish Python properties
    py::array_t<std::int8_t>
    neutral() const
    {
        return genotype_matrix_wrapper(dmatrix->neutral,
                                       dmatrix->neutral_keys);
    }

    py::array_t<double>
    neutral_positions() const
    {
        return fwdpy11::make_1d_ndarray_readonly(dmatrix->neutral.positions);
    }

    py::array_t<std::size_t>
    neutral_keys() const
    {
        return fwdpy11::make_1d_ndarray_readonly(dmatrix->neutral_keys);
    }

    py::array_t<std::int8_t>
    selected() const
    {
        return genotype_matrix_wrapper(dmatrix->selected,
                                       dmatrix->selected_keys);
    }

    py::array_t<double>
    selected_positions() const
    {
        return fwdpy11::make_1d_ndarray_readonly(dmatrix->selected.positions);
    }

    py::array_t<std::size_t>
    selected_keys() const
    {
        return fwdpy11::make_1d_ndarray_readonly(dmatrix->selected_keys);
    }
};

void
init_DataMatrixIterator(py::module& m)
{
    py::class_<DataMatrixIterator>(m, "DataMatrixIterator",
                                   R"delim(
        Efficient iteration over genomic windows.

        This class allows efficient traversal across multiple
        genomic intervals.  The class is iterable,
        and encapsulates the data fields of 
        :class:`fwdpy11.DataMatrix`, which are accessible
        as properties.

        .. versionadded:: 0.4.4
        )delim")
        .def(py::init<const fwdpp::ts::table_collection&,
                      const std::vector<fwdpy11::Mutation>&,
                      const std::vector<fwdpp::ts::TS_NODE_INT>&,
                      const std::vector<std::pair<double, double>>&, bool,
                      bool, bool>(),
             py::arg("tables"), py::arg("mutations"), py::arg("samples"),
             py::arg("intervals"), py::arg("neutral"), py::arg("selected"),
             py::arg("fixations") = false,
             R"delim(
        :param tables: A table collection
        :type tables: :class:`fwdpy11.TableCollection`
        :param mutations: The mutations from the simulation
        :type mutations: :class:`fwdpy11.MutationVector`
        :param samples: A list of samples
        :type samples: List-like
        :param intervals: The :math:`[start, stop)` positions of each interval
        :type intervals: list of tuples
        :param neutral: If True, include neutral variants
        :type neutral: boolean
        :param selected: If True, include selected variants
        :type selected: boolean
        :param fixations: (False) If True, include fixations in the sample
        :type fixations: boolean
        )delim")
        .def("__iter__",
             [](DataMatrixIterator& v) -> DataMatrixIterator& { return v; })
        .def("__next__", &DataMatrixIterator::next_data_matrix)
        .def_property_readonly("neutral", &DataMatrixIterator::neutral,
                               "Genotypes at neutral variants")
        .def_property_readonly("neutral_keys",
                               &DataMatrixIterator::neutral_keys,
                               "Indexes of the neutral variants.")
        .def_property_readonly("neutral_positions",
                               &DataMatrixIterator::neutral_positions,
                               "Positions of neutral variants.")
        .def_property_readonly("selected", &DataMatrixIterator::selected,
                               "Genotypes at selected variants.")
        .def_property_readonly("selected_keys",
                               &DataMatrixIterator::selected_keys,
                               "Indexes of selected variants.")
        .def_property_readonly("selected_positions",
                               &DataMatrixIterator::selected_positions,
                               "Positions of selected variants.");
}
