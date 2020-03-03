#include <memory>
#include <functional>
#include <gsl/gsl_matrix.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fwdpy11/regions/mvDES.hpp>
#include <fwdpy11/numpy/array.hpp>

namespace py = pybind11;

namespace
{
    using matrix_ptr = std::unique_ptr<gsl_matrix, std::function<void(gsl_matrix *)>>;
    fwdpy11::mvDES
    create_from_python(py::list sregions, py::array_t<double> means,
                       py::array_t<double> vcov)
    {
        auto means_buffer = means.request();
        if (means_buffer.ndim != 1)
            {
                throw std ::invalid_argument("means must be a 1d ndarray");
            }
        std::vector<double> mu(static_cast<double *>(means_buffer.ptr),
                               static_cast<double *>(means_buffer.ptr) + means.shape(0));
        auto vcov_unchecked = vcov.unchecked<2>();
        if (vcov_unchecked.shape(0) != vcov_unchecked.shape(1))
            {
                throw std::invalid_argument("input matrix must be square");
            }
        matrix_ptr m(gsl_matrix_alloc(vcov_unchecked.shape(0), vcov_unchecked.shape(1)),
                     [](gsl_matrix *m) { gsl_matrix_free(m); });
        for (py::ssize_t i = 0; i < vcov_unchecked.shape(0); ++i)
            {
                for (py::ssize_t j = 0; j < vcov_unchecked.shape(1); ++j)
                    {
                        gsl_matrix_set(m.get(), i, j, vcov_unchecked(i, j));
                    }
            }

        std::vector<std::unique_ptr<fwdpy11::Sregion>> output_distributions;
        for (auto s : sregions)
            {
                output_distributions.emplace_back(
                    s.cast<const fwdpy11::Sregion &>().clone());
            }
        return fwdpy11::mvDES(output_distributions, std::move(mu), *m);
    }
}

static const auto mvDES_CLASS_DOCSTRING = R"delim(
General multivariate distribution of effect sizes.

.. versionadded:: 0.7.0
)delim";

static const auto mvDES_INIT_DOCSTRING = R"delim(
:param des: Distributions of effect sizes
:type des: list
:param means: means marginal gaussian Distributions
:type means: numpy.ndarray
:param matrix: Variance/covariance matrix
:type matrix: numpy.ndarray
)delim";

void
init_mvDES(py::module &m)
{
    py::class_<fwdpy11::mvDES, fwdpy11::Sregion>(m, "mvDES", mvDES_CLASS_DOCSTRING)
        .def(py::init([](py::list sregions, py::array_t<double> means,
                         py::array_t<double> vcov) {
                 return create_from_python(sregions, means, vcov);
             }),
             py::arg("des"), py::arg("means"), py::arg("matrix"), mvDES_INIT_DOCSTRING)
        .def("__repr__", &fwdpy11::mvDES::repr)
        .def_property_readonly("means", &fwdpy11::mvDES::get_means)
        .def_property_readonly("matrix", &fwdpy11::mvDES::get_matrix)
        .def(py::pickle([](const fwdpy11::mvDES &self) { return self.pickle(); },
                        &fwdpy11::mvDES::unpickle));
}
