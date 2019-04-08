#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <fwdpy11/genetic_values/MultivariateGSS.hpp>

namespace py = pybind11;

void
init_MultivariateGSS(py::module& m)
{
    py::class_<fwdpy11::MultivariateGSS,
               fwdpy11::MultivariateGeneticValueToFitnessMap>(
        m, "MultivariateGSS",
        R"delim(
        Multivariate gaussian stablizing selection.

        Maps a multidimensional trait to fitness using the Euclidian 
        distance of a vector of trait values to a vector of optima.

        Essentially, this is Equation 1 of

        Simons, Yuval B., Kevin Bullaughey, Richard R. Hudson, and Guy Sella. 2018.
        "A Population Genetic Interpretation of GWAS Findings for Human Quantitative Traits."
        PLoS Biology 16 (3): e2002985.

        For the case of moving optima, see :class:`fwdpy11.MultivariateGSSmo`.
        )delim")
        .def(py::init([](py::array_t<double> optima, double VS) {
            auto r = optima.unchecked<1>();
            std::vector<double> voptima(r.data(0), r.data(0) + r.shape(0));
            return fwdpy11::MultivariateGSS(std::move(voptima), VS);
        }),
        R"delim(
        :param optima: The optimum value for each trait
        :type optima: np.array
        :param VS: Strength of stablizing selection
        :type VS: float

        .. note::

            VS is :math:`\omega^2` in the Simons et al. notation
        )delim")
        .def(py::pickle(
            [](const fwdpy11::MultivariateGSS& self) { return self.pickle(); },
            [](py::object o) {
                py::tuple t = o.cast<py::tuple>();
                if (t.size() != 2)
                    {
                        throw std::invalid_argument("incorrect tuple size");
                    }
                std::vector<double> optima = t[0].cast<std::vector<double>>();
                double VS = t[1].cast<double>();
                return fwdpy11::MultivariateGSS(std::move(optima), VS);
            }));
}
