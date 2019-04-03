#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <fwdpy11/regions/MultivariateGaussianEffects.hpp>

namespace py = pybind11;

void
init_MultivariateGaussianEffects(py::module& m)
{
    py::class_<fwdpy11::MultivariateGaussianEffects, fwdpy11::Sregion>(
        m, "MultivariateGaussianEffects",
        R"delim(
        Pleiotropic effects via a multivariate Gaussian distribution.

        This class can be used to generate mutations with both vectors
        of effect sizes as well as a separate fixed effect.

        .. versionadded:: 0.3.0
        )delim")
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
                     fwdpy11::Region(beg, end, weight, coupled, label), 1.0,
                     v.matrix, fixed_effect, h, true);
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"),
             py::arg("matrix"), py::arg("fixed_effect") = 0.0,
             py::arg("h") = 1.0, py::arg("coupled") = true,
             py::arg("label") = 0,
             R"delim(
             Constructor

             :param beg: Beginning of the region
             :type beg: float
             :param end: End of the region
             :type end: float
             :param weight: Weight on the region
             :type weight: float
             :param matrix: Variance-covariance matrix
             :type matrix: numpy ndarray
             :param fixed_effect: Fixed effect size. Defaults to 0.0.
             :type fixed_effect: float
             :param h: Dominance. Defaults to 1.0
             :type h: float
             :param coupled: Specify if weight is function of end-beg or not. Defaults to True
             :type coupled: bool
             :param label: Label for mutations from this region. Defaults to 0.
             :type label: np.uint16

             The input matrix must be square and semi-positive definite.   If either
             of these conditions are not met, ValueError will be raised. ValueError
             will also be raised if the input matrix contains any non-finite values.
            
             .. note::

                The dominance parameter (`h`) applies to both the fixed effect and those
                drawn from a multivariate normal.
             )delim")
        .def(py::pickle(
            [](const fwdpy11::MultivariateGaussianEffects& self) {
                return self.pickle();
            },
            [](py::tuple t) {
                return fwdpy11::MultivariateGaussianEffects::unpickle(t);
            }))
        .def("__repr__",&fwdpy11::MultivariateGaussianEffects::repr)
        .def("__eq__", [](const fwdpy11::MultivariateGaussianEffects& lhs,
                          const fwdpy11::MultivariateGaussianEffects& rhs) {
            return lhs == rhs;
        });
}
