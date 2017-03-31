<% 
setup_pybind11(cfg) 
import fwdpy11 as fp11 
cfg['include_dirs'] = [ fp11.get_includes(), fp11.get_fwdpp_includes() ] 
%>

#include <pybind11/pybind11.h>
#include <fwdpy11/types.hpp>
#include <fwdpy11/fitness/fitness.hpp>

namespace py = pybind11;

struct snowdrift_diploid
{
    using result_type = double;
    inline result_type
    operator()(const fwdpy11::diploid_t &dip, const fwdpy11::gcont_t &,
               const fwdpy11::mcont_t &, const std::vector<double> &phenotypes,
               const double b1, const double b2, const double c1,
               const double c2) const
    {
        unsigned N = phenotypes.size();
        auto i = dip.label;
        double zself = phenotypes[i];
        result_type fitness = 0;
        for (unsigned j = 0; j < N; ++j)
            {
                if (i != j)
                    {
                        double zpair = zself + phenotypes[j];
                        // Payoff function from Fig 1
                        double a = b1 * zpair + b2 * zpair * zpair - c1 * zself
                                   - c2 * zself * zself;
                        fitness += 1 + std::max(a, 0.0);
                    }
            }
        return fitness / double(N - 1);
    }
};

struct snowdrift : public fwdpy11::singlepop_fitness
{
    double b1, b2, c1, c2;
    std::vector<double> phenotypes;

    snowdrift(double b1_, double b2_, double c1_, double c2_)
        : b1(b1_), b2(b2_), c1(c1_), c2(c2_), phenotypes(std::vector<double>())
    {
    }

    inline fwdpy11::singlepop_fitness_fxn
    callback() const final
    {
        return std::bind(snowdrift_diploid(), std::placeholders::_1,
                         std::placeholders::_2, std::placeholders::_3,
                         std::cref(phenotypes), b1, b2, c1, c2);
    }
    void
    update(fwdpy11::singlepop_t &pop) final
    {
        phenotypes.resize(pop.N);
        std::size_t d = 0;
        for (auto &&dip : pop.diploids)
            {
                dip.label = d;
                phenotypes[d++] = KTfwd::additive_diploid()(
                    dip, pop.gametes, pop.mutations, 2.0);
            }
    };
};

PYBIND11_PLUGIN(snowdrift)
{
    pybind11::module m("snowdrift", "Example of custom stateful fitness model.");

    py::object base = (py::object) py::module::import("fwdpy11.fitness").attr("SpopFitness");


    py::class_<snowdrift, fwdpy11::singlepop_fitness>(m, "SpopSnowdrift")
        .def(py::init<double, double, double, double>(), py::arg("b1"),
             py::arg("b2"), py::arg("c1"), py::arg("c2"));

    return m.ptr();
}
