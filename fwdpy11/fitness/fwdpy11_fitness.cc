#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/types.hpp>

namespace py = pybind11;

using singlepop_fitness = std::function<double(const fwdpy11::diploid_t &,
                                               const fwdpy11::gcont_t &,
                                               const fwdpy11::mcont_t &)>;

template <typename fitness_model> struct fwdpp_singlepop_fitness_wrapper
{
    fitness_model f;
    const double scaling;
    fwdpp_singlepop_fitness_wrapper(const double scaling_)
        : f(fitness_model()), scaling(scaling_)
    {
    }
    inline double
    operator()(const fwdpy11::diploid_t &dip, const fwdpy11::gcont_t &gametes,
               const fwdpy11::mcont_t &mutations)
    {
        return f(dip, gametes, mutations, scaling);
    }
};

PYBIND11_PLUGIN(fwdpy11_fitness)
{
    py::module m("fwdpy11_fitness", "Standard fitness models.");

    py::class_<fwdpp_singlepop_fitness_wrapper<KTfwd::multiplicative_diploid>>(
        m, "SpopMult")
        .def(py::init<double>())
        .def_property_readonly(
            "callback",
            [](const fwdpp_singlepop_fitness_wrapper<KTfwd::
                                                         multiplicative_diploid>
                   &m) -> singlepop_fitness { return m; });

    py::class_<fwdpp_singlepop_fitness_wrapper<KTfwd::additive_diploid>>(
        m, "SpopAdditive")
        .def(py::init<double>())
        .def_property_readonly(
            "callback",
            [](const fwdpp_singlepop_fitness_wrapper<KTfwd::additive_diploid> &m)
                -> singlepop_fitness { return m; });

    return m.ptr();
}
