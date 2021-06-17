#include <memory>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <fwdpp/ts/tree_visitor.hpp>
#include <fwdpp/ts/detail/generate_data_matrix_details.hpp>

namespace py = pybind11;

class DataMatrixIterator
// Encapsulate fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection> and fwdpp::data_matrix
// to provide very fast interation over multiple genomic intervals.
// We accomplish this by copying the state of the current tree_visitor
// whenever we recognize that the current tree overlaps with the next interval.
// The saved state can be restored in O(1) time, "snapping" us to the start of
// next interval without having to initiate a new tree_visitor.
{
  private:
    using site_table_itr = fwdpp::ts::std_table_collection::site_table::const_iterator;
    using mut_table_itr
        = fwdpp::ts::std_table_collection::mutation_table::const_iterator;
    std::unique_ptr<fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection>>
        current_tree, next_tree;
    const std::vector<std::pair<double, double>> position_ranges;
    std::vector<std::int8_t> genotypes;
    std::unordered_map<std::size_t, double> mutation_positions;
    std::unique_ptr<fwdpp::data_matrix> dmatrix;
    const site_table_itr sbeg, send;
    site_table_itr scurrent;
    const mut_table_itr mbeg, mend;
    mut_table_itr mcurrent;
    std::size_t current_range;
    const bool include_neutral_variants, include_selected_variants, include_fixations;
    bool matrix_requires_clearing;

    std::unique_ptr<fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection>>
    initialize_current_tree(const fwdpp::ts::std_table_collection& tables,
                            const std::vector<fwdpp::ts::table_index_t>& samples)
    {
        std::unique_ptr<fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection>> rv(
            new fwdpp::ts::tree_visitor<fwdpp::ts::std_table_collection>(
                tables, samples, fwdpp::ts::update_samples_list(true)));
        auto flag = rv->operator()();
        if (flag == false)
            {
                throw std::invalid_argument("TableCollection contains no trees");
            }
        return rv;
    }

    std::vector<std::pair<double, double>>
    init_intervals(const std::vector<std::pair<double, double>>& input_intervals)
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
                        throw std::invalid_argument("invalid interval: end <= beg");
                    }
            }
        for (std::size_t i = 1; i < input_intervals.size(); ++i)
            {
                if (!(input_intervals[i].first > input_intervals[i - 1].first))
                    {
                        throw std::invalid_argument("invalid interval start positions");
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

    std::unordered_map<std::size_t, double>
    set_positions(const fwdpp::ts::std_table_collection& tables)
    {
        std::unordered_map<std::size_t, double> rv;
        for (auto& m : tables.mutations)
            {
                if (rv.find(m.key) != end(rv))
                    {
                        throw fwdpp::ts::tables_error(
                            "mutation key present more than once in "
                            "MutationTable");
                    }
                if (m.site >= tables.sites.size())
                    {
                        throw fwdpp::ts::tables_error("invalid site id");
                    }
                rv[m.key] = tables.sites[m.site].position;
            }
        return rv;
    }

    site_table_itr
    find_first_site(site_table_itr b, site_table_itr e, const double start)
    {
        return std::lower_bound(
            b, e, start,
            [](const fwdpp::ts::site& s, const double v) { return s.position < v; });
    }

    void
    set_current_mutation_to_current_site()
    {
        mcurrent = std::lower_bound(
            mbeg, mend, scurrent->position,
            [this](const fwdpp::ts::mutation_record& mr, const double p) {
                return (sbeg + mr.site)->position < p;
            });
    }

    site_table_itr
    init_trees_and_sites()
    {
        if (current_tree == nullptr || position_ranges.empty())
            {
                throw std::runtime_error("DataMatrixIterator __init__ failure");
            }
        double first_left = position_ranges.front().first;
        double current_right = current_tree->tree().right;

        while (current_right < first_left)
            {
                current_tree->operator()();
                current_right = current_tree->tree().right;
            }
        site_table_itr s = sbeg;
        while (s < send && s->position < first_left)
            {
                ++s;
            }
        return s;
    }

    site_table_itr
    find_first_site_after_current_max()
    // After a call to cleanup_matrix, we may need to reset scurrent
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
                return scurrent;
            }
        return std::upper_bound(
            scurrent, send, m,
            [](double v, const fwdpp::ts::site& s) { return v < s.position; });
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
                auto flag = current_tree->operator()();
                if (!flag)
                    {
                        break;
                    }
            }
    }

    site_table_itr
    advance_sites()
    {
        matrix_requires_clearing = false;
        if (next_tree != nullptr)
            {
                current_tree.swap(next_tree);
                next_tree.reset(nullptr);
                double left = position_ranges[current_range].first;
                double right = position_ranges[current_range].second;
                scurrent = find_first_site(sbeg, send, left);
                cleanup_matrix(left, right);
                scurrent = find_first_site_after_current_max();
            }
        else
            {
                matrix_requires_clearing = true;
                double left = position_ranges[current_range].first;
                while (scurrent < send && scurrent->position < left)
                    {
                        ++scurrent;
                    }
            }
        set_current_mutation_to_current_site();
        return scurrent;
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
    cleanup_matrix_details(fwdpp::state_matrix& sm, std::vector<std::size_t>& keys,
                           double l, double r)
    // When genomic intervals overlap, they have mutations in common.
    // This function removes all mutations from the previous window,
    // keeping any mutations shared by both windows.
    {
        // find first key corresponding to position >= l
        auto itr = std::lower_bound(
            begin(keys), end(keys), l,
            [this](std::size_t k, double v) { return mutation_positions[k] < v; });
        // This is the number of mutations
        // with positions < l
        auto d = std::distance(begin(keys), itr);
        if (d > 0)
            {
                // Sanity check
                auto pitr = std::lower_bound(begin(sm.positions), end(sm.positions), l);
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
                sm.data.erase(begin(sm.data), begin(sm.data) + d * dmatrix->ncol);
            }

        // Give all mutations >= r the same treatment
        itr = std::lower_bound(
            begin(keys), end(keys), r,
            [this](std::size_t k, double v) { return mutation_positions[k] < v; });
        d = std::distance(itr, end(keys));
        if (d > 0)
            {
                auto pitr = std::lower_bound(begin(sm.positions), end(sm.positions), r);
                if (d != std::distance(pitr, end(sm.positions)))
                    {
                        throw std::runtime_error(
                            "DataMatrix internal state inconsistent");
                    }
                keys.erase(itr, end(keys));
                sm.positions.erase(pitr, end(sm.positions));
                sm.data.erase(end(sm.data) - (d * dmatrix->ncol), end(sm.data));
            }
    }

    void
    cleanup_matrix(double l, double r)
    // Removes all data in dmatrix
    // corresponding to position < p
    {
        cleanup_matrix_details(dmatrix->neutral, dmatrix->neutral_keys, l, r);
        cleanup_matrix_details(dmatrix->selected, dmatrix->selected_keys, l, r);
    }

    void
    check_if_still_iterating()
    {
        if (!(scurrent < send) || current_range >= position_ranges.size())
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
        auto& k = (mut_is_neutral) ? dmatrix->neutral_keys : dmatrix->selected_keys;
        k.push_back(key);
        sm.positions.push_back(mutation_positions[key]);
        sm.data.insert(sm.data.end(), begin(genotypes), end(genotypes));
    }

    bool
    site_in_current_range(const site_table_itr itr)
    {
        if (!(itr < send))
            {
                return false;
            }
        if (current_range >= position_ranges.size())
            {
                throw std::runtime_error("DataMatrixIterator fatal error");
            }
        auto pos = itr->position;
        return pos >= position_ranges[current_range].first
               && pos < position_ranges[current_range].second;
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
    DataMatrixIterator(const fwdpp::ts::std_table_collection& tables,
                       const std::vector<fwdpp::ts::table_index_t>& samples,
                       const std::vector<std::pair<double, double>>& intervals,
                       bool neutral, bool selected, bool fixations)
        : current_tree(initialize_current_tree(tables, samples)), next_tree(nullptr),
          position_ranges(init_intervals(intervals)), genotypes(samples.size(), 0),
          mutation_positions(set_positions(tables)),
          dmatrix(new fwdpp::data_matrix(samples.size())), sbeg(begin(tables.sites)),
          send(end(tables.sites)), scurrent(init_trees_and_sites()),
          mbeg(begin(tables.mutations)), mend(end(tables.mutations)),
          mcurrent(begin(tables.mutations)), current_range(0),
          include_neutral_variants(neutral), include_selected_variants(selected),
          include_fixations(fixations), matrix_requires_clearing(false)
    {
    }

    DataMatrixIterator&
    next_data_matrix()
    // This is the back end for the Python
    // class's __next__ function
    {
        check_if_still_iterating();
        advance_trees();
        scurrent = advance_sites();
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
                if (next_tree == nullptr && current_range + 1 < position_ranges.size())
                    {
                        double right = tree.right;
                        if (right > position_ranges[current_range + 1].first)
                            {
                                next_tree.reset(
                                    new fwdpp::ts::tree_visitor<
                                        fwdpp::ts::std_table_collection>(*current_tree));
                            }
                    }
                for (; scurrent < send && scurrent->position < tree.right
                       && site_in_current_range(scurrent);
                     ++scurrent)
                    {
                        while (mcurrent < mend
                               && (sbeg + mcurrent->site)->position < scurrent->position)
                            {
                                ++mcurrent;
                            }
                        auto m = mcurrent;
                        m++;
                        while (m < mend
                               && (sbeg + m->site)->position == scurrent->position)
                            {
                                ++m;
                            }
                        fwdpp::ts::detail::process_site_range(
                            tree, scurrent, std::make_pair(mcurrent, m),
                            include_neutral_variants, include_selected_variants,
                            !include_fixations, genotypes, *dmatrix);
                        mcurrent = m;
                    }
                iteration_flag = current_tree->operator()();
                tree_in_current_range = (current_tree->tree().left
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
        return genotype_matrix_wrapper(dmatrix->neutral, dmatrix->neutral_keys);
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
        return genotype_matrix_wrapper(dmatrix->selected, dmatrix->selected_keys);
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
    py::class_<DataMatrixIterator>(m, "ll_DataMatrixIterator")
        .def(py::init<const fwdpp::ts::std_table_collection&,
                      const std::vector<fwdpp::ts::table_index_t>&,
                      const std::vector<std::pair<double, double>>&, bool, bool, bool>(),
             py::arg("tables"), py::arg("samples"), py::arg("intervals"),
             py::arg("neutral"), py::arg("selected"), py::arg("fixations") = false)
        .def("_ll_iter", [](DataMatrixIterator& v) -> DataMatrixIterator& { return v; })
        .def("_ll_next", &DataMatrixIterator::next_data_matrix)
        .def_property_readonly("_neutral", &DataMatrixIterator::neutral)
        .def_property_readonly("_neutral_keys", &DataMatrixIterator::neutral_keys)
        .def_property_readonly("_neutral_positions",
                               &DataMatrixIterator::neutral_positions)
        .def_property_readonly("_selected", &DataMatrixIterator::selected)
        .def_property_readonly("_selected_keys", &DataMatrixIterator::selected_keys)
        .def_property_readonly("_selected_positions",
                               &DataMatrixIterator::selected_positions);
}
