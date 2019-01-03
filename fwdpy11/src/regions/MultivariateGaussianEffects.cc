#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <fwdpy11/regions/MultivariateGaussianEffects.hpp>

namespace py = pybind11;

void
init_MultivariateGaussianEffects(py::module& m)
{
    py::class_<fwdpy11::MultivariateGaussianEffects, fwdpy11::Sregion>(
        m, "MultivariateGaussianEffects")
        .def(py::init([](double beg, double end, double weight, bool coupled,
                         py::array_t<double> cov_matrix, double h,
                         std::uint16_t label) {
            auto r = cov_matrix.unchecked<2>();
            if (r.shape(0) != r.shape(1))
                {
                    throw std::invalid_argument("input matrix is not square");
                }
            gsl_matrix_const_view v = gsl_matrix_const_view_array(
                r.data(0,0), r.shape(0), r.shape(1));
            return fwdpy11::MultivariateGaussianEffects(beg, end, weight,
                                                        coupled, v, h, label);
        }));
}
