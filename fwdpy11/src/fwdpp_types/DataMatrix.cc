#include <fwdpp/data_matrix.hpp>
#include <fwdpp/io/scalar_serialization.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void
init_data_matrix(py::module &m)
{
    py::class_<fwdpp::state_matrix>(m, "StateMatrix", py::buffer_protocol(),
                                    R"delim(
            Simple matrix representation of variation data.

            These are not constructed directly.  Rather,
            they are generated when a 
            :class:`fwdpy11.DataMatrix` is generated.

            This object supports the buffer protocol.

            .. versionadded:: 0.2.0
            )delim")
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
        .def_readonly("positions", &fwdpp::state_matrix::positions,
                      "The mutation positions")
        .def_buffer([](const fwdpp::state_matrix &sm) -> py::buffer_info {
            using value_type = std::int8_t;
            auto nrow = sm.positions.size();
            auto ncol = (nrow > 0) ? sm.data.size() / nrow : 0;
            return py::buffer_info(
                const_cast<value_type *>(sm.data.data()), sizeof(value_type),
                py::format_descriptor<value_type>::format(), 2, { nrow, ncol },
                { sizeof(value_type) * ncol, sizeof(value_type) });
        });

    py::class_<fwdpp::data_matrix>(m, "DataMatrix",
                                   R"delim(
		Represent a sample from a population in a matrix format.

		There are two possible representations of the data:

		1. As a genotype matrix, where individuals are encoded a 0,1, or 2
		copies of the derived mutation. There is one column per diploid here,
        and one row per variable site.

		2. As a haplotype matrix, with two columns per diploid, and each
		column containing a 0 (ancestral) or 1 (derived) label. Each row
        represents a variable site.

        .. versionchanged:: 0.2.0

            Changed layout to row = variable site. 
            Changed to match fwdpp 0.7.0 layout where the neutral
            and selected data are represented as a 
            :class:`fwdpy11.StateMatrix`
		)delim")
        .def_readwrite("neutral", &fwdpp::data_matrix::neutral,
                       R"delim(
                Return a buffer representing neutral variants.
                This buffer may be used to create a NumPy
                ndarray object.

                .. versionchanged:: 0.1.2
                    Return a buffer instead of 1d numpy.array

                .. versionchanged:: 0.1.4
                    Allow read/write access instead of readonly

                .. versionchanged:: 0.2.0
                    Type is :class:`fwdpy11.StateMatrix`
                )delim")
        .def_readwrite("selected", &fwdpp::data_matrix::selected,
                       R"delim(
                Return a buffer representing neutral variants.
                This buffer may be used to create a NumPy
                ndarray object.

                .. versionchanged:: 0.1.2
                    Return a buffer instead of 1d numpy.array

                .. versionchanged:: 0.1.4
                    Allow read/write access instead of readonly

                .. versionchanged:: 0.2.0
                    Type is :class:`fwdpy11.StateMatrix`
                )delim")
        .def_readonly("ncol", &fwdpp::data_matrix::ncol,
                      "Sample size of the matrix")
        .def_readonly("neutral_keys", &fwdpp::data_matrix::neutral_keys,
                      "Keys for neutral mutations used to generate matrix")
        .def_readonly("selected_keys", &fwdpp::data_matrix::selected_keys,
                      "Keys for selected mutations used to generate matrix")
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
                if (nsites)
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

