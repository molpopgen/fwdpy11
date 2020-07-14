#include <fwdpy11/genetic_value_to_fitness/GeneticValueIsTrait.hpp>
#include <fwdpy11/numpy/array.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

namespace py = pybind11;

struct genetic_values_buffer_proxy
{
    double* data;
    std::size_t size;
    genetic_values_buffer_proxy() : data{nullptr}, size{0}
    {
    }
};

struct GeneticValueIsTraitData
{
    fwdpy11::DiploidMetadata offspring_metadata_copy;
    py::object offspring_metadata;
    genetic_values_buffer_proxy buffer;
    py::object genetic_values;

    GeneticValueIsTraitData()
        : offspring_metadata_copy{},
          offspring_metadata{
              py::cast<fwdpy11::DiploidMetadata*>(&offspring_metadata_copy)},
          buffer{}, genetic_values{py::cast<genetic_values_buffer_proxy*>(&buffer)}
    {
    }
};

class GeneticValueIsTraitTrampoline : public fwdpy11::GeneticValueIsTrait
// Trampoline class allowing custom
// GeneticValueIsTrait to be written
// in Python
{
  private:
    mutable GeneticValueIsTraitData data;
    py::object pydata;

  public:
    GeneticValueIsTraitTrampoline(std::size_t ndim)
        : fwdpy11::GeneticValueIsTrait(ndim), data{},
          pydata{py::cast<GeneticValueIsTraitData*>(&data)}
    {
    }

    double
    operator()(const fwdpy11::DiploidMetadata& metadata,
               const std::vector<double>& genetic_values) const override
    {
        data.offspring_metadata_copy = metadata;
        data.buffer.data = const_cast<double*>(genetic_values.data());
        data.buffer.size = genetic_values.size();
        PYBIND11_OVERLOAD_PURE_NAME(double, fwdpy11::GeneticValueIsTrait,
                                    "__call__", operator(), pydata);
    }

    void
    update(const fwdpy11::DiploidPopulation& pop) override
    {
        PYBIND11_OVERLOAD_PURE(void, fwdpy11::GeneticValueIsTrait, update, pop);
    }

    std::shared_ptr<fwdpy11::GeneticValueToFitnessMap>
    clone() const override
    // Implementation details from pybind11 issue 1049
    {
        auto self = py::cast(this);
        auto cloned = self.attr("clone")();

        auto keep_python_state_alive = std::make_shared<py::object>(cloned);
        auto ptr = cloned.cast<GeneticValueIsTraitTrampoline*>();
        return std::shared_ptr<GeneticValueToFitnessMap>(keep_python_state_alive, ptr);
    }
};

void
init_GeneticValueIsTrait(py::module& m)
{
    py::class_<fwdpy11::GeneticValueIsTrait, fwdpy11::GeneticValueToFitnessMap,
               GeneticValueIsTraitTrampoline>(
        m, "GeneticValueIsTrait",
        "ABC for functions mapping genetic values representing traits to "
        "fitness.")
        .def(py::init<std::size_t>(), py::arg("ndim") = 1);

    py::class_<genetic_values_buffer_proxy>(m, "_GeneticValuesBufferProxy",
                                            py::buffer_protocol())
        .def_buffer([](const genetic_values_buffer_proxy& self) {
            return py::buffer_info(self.data, sizeof(double),
                                   py::format_descriptor<double>::format(), 1,
                                   {self.size}, {sizeof(double)}
                                   );
        });

    py::class_<GeneticValueIsTraitData>(m, "PyGeneticValueIsTraitData")
        .def_readonly("offspring_metadata", &GeneticValueIsTraitData::offspring_metadata)
        .def_readonly("genetic_values", &GeneticValueIsTraitData::genetic_values);
}
