#include <fwdpy11/genetic_values/StrictAdditive.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void
init_StrictAdditive(py::module& m)
{
    py::class_<fwdpy11::StrictAdditive>(m, "StrictAdditive")
        .def(py::init<>())
        .def(py::init<const fwdpy11::GeneticValueToFitnessMap&>())
        .def(py::init<const fwdpy11::GeneticValueToFitnessMap&,
                      const fwdpy11::GeneticValueNoise&>())
        .def(py::pickle(
            [](const fwdpy11::StrictAdditive& self) {
                auto p = py::module::import("pickle");
                return py::make_tuple(
                    self.pickle(), p.attr("dumps")(self.gv2w->clone(), -1),
                    p.attr("dumps")(self.noise_fxn->clone(), -1));
            },
            [](py::tuple t) {
                if (t.size() != 3)
                    {
                        throw std::runtime_error("invalid object state");
                    }

                auto p = py::module::import("pickle");
                auto t1 = p.attr("loads")(t[1]);
                auto t2 = p.attr("loads")(t[2]);
                //Do the casts in the constructor
                //to avoid any nasty issues w/
                //refs to temp
                return fwdpy11::StrictAdditive(
                    t1.cast<const fwdpy11::GeneticValueToFitnessMap&>(),
                    t2.cast<const fwdpy11::GeneticValueNoise&>());
            }));
}
