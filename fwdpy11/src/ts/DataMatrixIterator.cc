#include <memory>
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
        // TODO: validate ranges, sorting, and overlaps
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
          include_neutral_variants(neutral), include_selected_variants(selected),
          include_fixations(fixations)
    {
    }
};

void
init_DataMatrixIterator(py::module& m)
{
    py::class_<DataMatrixIterator>(m, "DataMatrixIterator");
}
