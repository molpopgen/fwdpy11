#include <memory>
#include <cmath>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/Population.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <fwdpp/ts/tree_visitor.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpy11::Mutation>);

class DataMatrixIterator
{
  private:
    std::unique_ptr<fwdpp::ts::tree_visitor> current_tree, next_tree;
    const std::vector<std::pair<double, double>> position_ranges;
    const bool include_neutral_variants, include_selected_variants,
        include_fixations;

    std::vector<std::pair<double, double>>
    init_intervals(
        const std::vector<std::pair<double, double>>& input_intervals)
    {
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

  public:
    DataMatrixIterator(const fwdpp::ts::table_collection& tables,
                       const std::vector<fwdpy11::Mutation>& mutations,
                       const std::vector<fwdpp::ts::TS_NODE_INT>& samples,
                       const std::vector<std::pair<double, double>>& intervals,
                       bool neutral, bool selected, bool fixations)
        : current_tree(new fwdpp::ts::tree_visitor(tables, samples)),
          next_tree(nullptr), position_ranges(init_intervals(intervals)),
          include_neutral_variants(neutral),
          include_selected_variants(selected), include_fixations(fixations)
    {
    }
};

void
init_DataMatrixIterator(py::module& m)
{
    py::class_<DataMatrixIterator>(m, "DataMatrixIterator");
}
