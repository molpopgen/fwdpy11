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
        pybind11::handle lognormal
            = pybind11::module::import("fwdpy11").attr("LogNormalS");
        pybind11::handle mvgaussian
            = pybind11::module::import("fwdpy11").attr("MultivariateGaussianEffects");
        for (auto s : sregions)
            {
                if (pybind11::isinstance(s, lognormal))
                    {
                        throw std::invalid_argument(
                            "incorrect init method using LogNormalS");
                    }
                if (pybind11::isinstance(s, mvgaussian))
                    {
                        throw std::invalid_argument(
                            "incorrect init method using MultivariateGaussianEffects");
                    }
                output_distributions.emplace_back(
                    s.cast<const fwdpy11::Sregion &>().clone());
            }
        return fwdpy11::mvDES(output_distributions, std::move(mu), *m);
    }

    std::vector<double>
    convert_means(py::array_t<double> means)
    {
        auto means_buffer = means.request();
        if (means_buffer.ndim != 1)
            {
                throw std ::invalid_argument("means must be a 1d ndarray");
            }
        std::vector<double> mu(static_cast<double *>(means_buffer.ptr),
                               static_cast<double *>(means_buffer.ptr) + means.shape(0));
        return mu;
    }

    matrix_ptr
    convert_matrix(py::array_t<double> vcov)
    {
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
        return m;
    }
}

void
init_mvDES(py::module &m)
{
    py::class_<fwdpy11::mvDES, fwdpy11::Sregion>(m, "_ll_mvDES")
        .def(py::init([](py::list sregions, py::array_t<double> means,
                         py::array_t<double> vcov) {
                 return create_from_python(sregions, means, vcov);
             }),
             py::arg("des"), py::arg("means"), py::arg("matrix"))
        .def(py::init([](const fwdpy11::LogNormalS &mvln, py::array_t<double> means,
                         py::array_t<double> vcov) {
                 if (mvln.univariate == true)
                     {
                         throw std::invalid_argument(
                             "LogNormalS instance must be created via "
                             "fwdpy11.LogNormalS.mv");
                     }
                 auto mu = convert_means(means);
                 auto matrix = convert_matrix(vcov);
                 return fwdpy11::mvDES(mvln, std::move(mu), *matrix);
             }),
             py::arg("des"), py::arg("means"), py::arg("matrix"))
        .def(py::init([](const fwdpy11::MultivariateGaussianEffects &mvgauss,
                         py::array_t<double> means, py::object) {
                 auto mu = convert_means(means);
                 if (mu.size() != mvgauss.input_matrix_copy->size1)
                     {
                         throw std::invalid_argument("incorrect number of mean values");
                     }
                 return fwdpy11::mvDES(mvgauss, std::move(mu));
             }),
             py::arg("des"), py::arg("means"), py::arg("matrix"));
}
