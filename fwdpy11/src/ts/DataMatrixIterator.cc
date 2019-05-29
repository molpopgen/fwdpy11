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
                            "invalid interval: all positions must be >= 0.0");
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
    advance_trees_and_mutations()
    {
        if (next_tree != nullptr)
            {
                current_tree.swap(next_tree);
                next_tree.reset(nullptr);
                double left = current_tree->tree().left;
                mcurrent = find_first_mutation_record(mbeg, mend, left);
            }
        else
            {
                while (mcurrent < mend)
                    {
                        const auto& m = current_tree->tree();
                        while (mutation_positions[mcurrent->key] < m.left
                               || mutation_positions[mcurrent->key] >= m.right)
                            {
                                auto flag = current_tree->operator()(
                                    std::true_type(), std::true_type());
                                if (flag == false)
                                    {
                                        throw std::runtime_error(
                                            "DataMatrixIterator: tree "
                                            "traversal "
                                            "error");
                                    }
                            }
                        if (m.leaf_counts[mcurrent->node] != 0)
                            {
                                // Skip over mutations in this tree
                                // that don't lead to samples
                                if ((is_neutral[mcurrent->key]
                                     && include_neutral_variants)
                                    || (!is_neutral[mcurrent->key]
                                        && include_selected_variants))
                                    {
                                        return mcurrent;
                                    }
                            }
                        ++mcurrent;
                    }
            }
        return mcurrent;
    }

    void
    release_memory()
    {
        current_tree.reset(nullptr);
        next_tree.reset(nullptr);
        dmatrix.reset(nullptr);
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

    void
    process_current_mutation(const fwdpp::ts::marginal_tree& tree,
                             const mut_table_itr mitr)
    {
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

  public:
    DataMatrixIterator(const fwdpp::ts::table_collection& tables,
                       const std::vector<fwdpy11::Mutation>& mutations,
                       const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                       const std::vector<std::pair<double, double>>& intervals,
                       bool neutral, bool selected, bool fixations)
        : current_tree(new fwdpp::ts::tree_visitor(tables, samples)),
          next_tree(nullptr), position_ranges(init_intervals(intervals)),
          genotypes(samples.size(), 0), is_neutral(set_neutral(mutations)),
          mutation_positions(set_positions(mutations)),
          dmatrix(new fwdpp::data_matrix(samples.size())),
          mbeg(find_first_mutation_record(tables.mutation_table.begin(),
                                          tables.mutation_table.end(),
                                          intervals[0].first)),
          mend(tables.mutation_table.end()), mcurrent(mbeg), current_range(0),
          include_neutral_variants(neutral),
          include_selected_variants(selected), include_fixations(fixations)
    {
        mcurrent = advance_trees_and_mutations();
    }

    DataMatrixIterator&
    next_data_matrix()
    {
        check_if_still_iterating();
        next_tree.reset(nullptr);

        bool iteration_flag = true;
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
            }
        while (iteration_flag == true);
        ++current_range;
        mcurrent = advance_trees_and_mutations();
        return *this;
    }
};

void
init_DataMatrixIterator(py::module& m)
{
    py::class_<DataMatrixIterator>(m, "DataMatrixIterator");
}
