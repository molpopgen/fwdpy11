#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/evolvets/SampleRecorder.hpp>

PYBIND11_MODULE(gvalue_recorder, m)
{
    m.def("record_gvalue", []() {
        return pybind11::cpp_function(
            [](const fwdpy11::DiploidPopulation &pop,
                pybind11::list l) {
                double mean_trait_value = 0.0;
                for (auto &md : pop.diploid_metadata)
                    {
                        mean_trait_value += md.g;
                    }
                l.append(mean_trait_value / static_cast<double>(pop.N));
            });
    });
}
