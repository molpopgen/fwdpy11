#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <fwdpy11/regions/MultivariateGaussianEffects.hpp>

namespace py = pybind11;

void
init_MultivariateGaussianEffects(py::module& m)
{
    py::class_<fwdpy11::MultivariateGaussianEffects, fwdpy11::Sregion>(
        m, "MultivariateGaussianEffects")
        .def(py::init([](double beg, double end, double weight,
                         py::array_t<double> cov_matrix, double fixed_effect,
                         double h, bool coupled, std::uint16_t label) {
                 auto r = cov_matrix.unchecked<2>();
                 if (r.shape(0) != r.shape(1))
                     {
                         throw std::invalid_argument(
                             "input matrix is not square");
                     }
                 gsl_matrix_const_view v = gsl_matrix_const_view_array(
                     r.data(0, 0), r.shape(0), r.shape(1));
                 return fwdpy11::MultivariateGaussianEffects(
                     beg, end, weight, coupled, v.matrix, fixed_effect, h,
                     label, true);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"),
             py::arg("covariance_matrix"), py::arg("fixed_effect") = 0.0,
             py::arg("h") = 1.0, py::arg("coupled") = true,
             py::arg("label") = 0);
}
