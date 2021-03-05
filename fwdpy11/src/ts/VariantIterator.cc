#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <fwdpp/ts/std_table_collection.hpp>
#include <fwdpp/ts/site_visitor.hpp>
#include <fwdpp/ts/marginal_tree_functions/samples.hpp>

namespace py = pybind11;

class VariantIterator
{
  private:
    using site_table_itr = fwdpp::ts::std_table_collection::site_table::const_iterator;
    fwdpp::ts::site_visitor<fwdpp::ts::std_table_collection> sv;
    site_table_itr scurrent, send;
    std::vector<std::int8_t> genotype_data;
    bool include_neutral, include_selected;
    double from, to;
    fwdpp::ts::convert_sample_index_to_nodes convert;
    py::list mutation_records;

  public:
    py::array_t<std::int8_t> genotypes;
    double current_position;
    VariantIterator(const fwdpp::ts::std_table_collection& tc,
                    const std::vector<fwdpp::ts::table_index_t>& samples,
                    const double beg, const double end,
                    const bool include_neutral_variant,
                    const bool include_selected_variants)
        : sv(tc, samples), scurrent(begin(tc.sites)), send(std::end(tc.sites)),
          genotype_data(samples.size(), 0), include_neutral(include_neutral_variant),
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
                        throw std::invalid_argument("invalid position interval");
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
                                    [this, i, &ni](fwdpp::ts::table_index_t u) {
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
                                throw fwdpp::ts::tables_error("invalid mutation data");
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
    py::class_<VariantIterator>(m, "ll_VariantIterator")
        .def(py::init([](const fwdpp::ts::std_table_collection& tables,
                         const std::vector<fwdpp::ts::table_index_t>& samples,
                         double begin, double end, bool include_neutral_variants,
                         bool include_selected_variants) {
                 return VariantIterator(tables, samples, begin, end,
                                        include_neutral_variants,
                                        include_selected_variants);
             }),
             py::arg("tables"), py::arg("samples"), py::arg("begin") = 0.0,
             py::arg("end") = std::numeric_limits<double>::max(),
             py::arg("include_neutral_variants") = true,
             py::arg("include_selected_variants") = true)
        .def("__iter__", [](VariantIterator& v) -> VariantIterator& { return v; })
        .def("__next__", &VariantIterator::next_variant)
        .def_readonly("_genotypes", &VariantIterator::genotypes)
        .def_property_readonly("site", &VariantIterator::current_site)
        .def_property_readonly("_records", &VariantIterator::records)
        .def_readonly("_position", &VariantIterator::current_position);
}

