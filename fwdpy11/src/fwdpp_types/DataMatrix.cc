#include <fwdpp/data_matrix.hpp>
#include <fwdpp/io/scalar_serialization.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace
{
    static const auto STATE_MATRIX_DOCSTRING = R"delim(
Simple matrix representation of variation data.

These are not constructed directly.
Rather, they are generated when a :class:`fwdpy11.DataMatrix` is generated .

This object supports the buffer protocol .

    .. versionadded:: 0.2.0
)delim";

    static const auto STATE_MATRIX_POSITIONS = R"delim(
The mutation positions.

.. versionchanged:: 0.6.1

    Type changed to :class:`numpy.ndarray`
)delim";

}

void
init_data_matrix(py::module &m)
{
    py::class_<fwdpp::state_matrix>(
        m, "StateMatrix", py::buffer_protocol(), STATE_MATRIX_DOCSTRING)
        .def_property_readonly(
            "shape",
            [](const fwdpp::state_matrix &sm) {
                if (sm.positions.empty())
                    {
                        return py::make_tuple(std::size_t(0), std::size_t(0));
                    }
                if (sm.data.empty())
                    {
                        throw std::runtime_error("StatMatrix data are empty");
                    }
                return py::make_tuple(sm.positions.size(),
                                      sm.data.size() / sm.positions.size());
            },
            "Shape of the matrix.")
        .def_property_readonly(
            "positions",
            [](const fwdpp::state_matrix &self) {
                return fwdpy11::make_1d_ndarray_readonly(self.positions);
            },
            STATE_MATRIX_POSITIONS)
        .def_buffer([](const fwdpp::state_matrix &sm) -> py::buffer_info {
            using value_type = std::int8_t;
            auto nrow = sm.positions.size();
            auto ncol = (nrow > 0) ? sm.data.size() / nrow : 0;
            return py::buffer_info(
                const_cast<value_type *>(sm.data.data()), sizeof(value_type),
                py::format_descriptor<value_type>::format(), 2, {nrow, ncol},
                {sizeof(value_type) * ncol, sizeof(value_type)});
        });

    py::class_<fwdpp::data_matrix, std::shared_ptr<fwdpp::data_matrix>>(
        m, "ll_DataMatrix")
        .def(py::init([](std::shared_ptr<fwdpp::data_matrix> &p) {
            if (p == nullptr)
                {
                    throw std::invalid_argument("input pointer is nullptr");
                }
            return p;
        }))
        .def_readwrite("_neutral", &fwdpp::data_matrix::neutral)
        .def_readwrite("_selected", &fwdpp::data_matrix::selected)
        .def_readonly("_ncol", &fwdpp::data_matrix::ncol)
        .def_property_readonly(
            "_neutral_keys",
            [](const fwdpp::data_matrix &self) {
                return fwdpy11::make_1d_ndarray_readonly(self.neutral_keys);
            })
        .def_property_readonly(
            "_selected_keys",
            [](const fwdpp::data_matrix &self) {
                return fwdpy11::make_1d_ndarray_readonly(self.selected_keys);
            })
        .def(py::pickle(
            [](const fwdpp::data_matrix &d) {
                std::ostringstream o;
                fwdpp::io::scalar_writer w;
                auto nsites = d.neutral.positions.size();
                auto dsize = d.neutral.data.size();
                w(o, &nsites, 1);
                w(o, &dsize, 1);
                if (nsites)
                    {
                        w(o, d.neutral.data.data(), d.neutral.data.size());
                        w(o, d.neutral.positions.data(), nsites);
                        w(o, d.neutral_keys.data(), nsites);
                    }
                nsites = d.selected.positions.size();
                dsize = d.selected.data.size();
                w(o, &nsites, 1);
                w(o, &dsize, 1);

                {
                    w(o, d.selected.data.data(), d.selected.data.size());
                    w(o, d.selected.positions.data(), nsites);
                    w(o, d.selected_keys.data(), nsites);
                }
                return py::bytes(o.str());
            },
            [](py::bytes b) {
                std::istringstream data(b);
                fwdpp::io::scalar_reader r;
                std::size_t nsites, dsize;
                r(data, &nsites);
                r(data, &dsize);
                fwdpp::data_matrix d(dsize);
                if (nsites)
                    {
                        d.neutral.data.resize(dsize);
                        r(data, d.neutral.data.data(), dsize);
                        d.neutral.positions.resize(nsites);
                        r(data, d.neutral.positions.data(), nsites);
                        d.neutral_keys.resize(nsites);
                        r(data, d.neutral_keys.data(), nsites);
                    }
                r(data, &nsites);
                r(data, &dsize);
                if (nsites)
                    {
                        d.selected.data.resize(dsize);
                        r(data, d.selected.data.data(), dsize);
                        d.selected.positions.resize(nsites);
                        r(data, d.selected.positions.data(), nsites);
                        d.selected_keys.resize(nsites);
                        r(data, d.selected_keys.data(), nsites);
                    }
                return std::unique_ptr<fwdpp::data_matrix>(
                    new fwdpp::data_matrix(std::move(d)));
            }));
}

