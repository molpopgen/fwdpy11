// clang-format off
<% 
setup_pybind11(cfg) 
#import fwdpy11 so we can find its C++ headers
import fwdpy11 as fp11 
#add fwdpy11 header locations to the include path
cfg['include_dirs'].extend([ fp11.get_includes(), fp11.get_fwdpp_includes() ])
#On OS X using clang, there is more work to do.  Using gcc on OS X
#gets rid of these requirements. The specifics sadly depend on how
#you initially built fwdpy11, and what is below assumes you used
#the provided setup.py + OS X + clang:
#cfg['compiler_args'].extend(['-stdlib=libc++','-mmacosx-version-min=10.7'])
#cfg['linker_args']=['-stdlib=libc++','-mmacosx-version-min=10.7']
#An alternative to the above is to add the first line to CPPFLAGS
#and the second to LDFLAGS when compiling a plugin on OS X using clang.
%>
// clang-format on

#include <algorithm> //for std::max
#include <pybind11/pybind11.h>
#include <fwdpp/fitness_models.hpp>
#include <fwdpy11/genetic_values/SlocusPopGeneticValue.hpp>
#include <fwdpy11/genetic_values/default_update.hpp>

    struct GeneralW : public fwdpy11::SlocusPopGeneticValue
{
    fwdpp::site_dependent_genetic_value w;

    GeneralW() : fwdpy11::SlocusPopGeneticValue{}, w{} {}

    inline double
    operator()(const std::size_t diploid_index,
               const fwdpy11::SlocusPop& pop) const
    {
        return std::max(
            0.0,
            w(pop.diploids[diploid_index], pop.gametes, pop.mutations,
              [](double& g, const fwdpy11::Mutation& m) { g *= (1.0 + m.s); },
              [](double& g, const fwdpy11::Mutation& m) {
                  g *= (1.0 + m.h);
              }));
    }

    double
    genetic_value_to_fitness(const fwdpy11::DiploidMetadata& metadata) const
    {
        return metadata.g;
    }

    double
    noise(const fwdpy11::GSLrng_t& /*rng*/,
          const fwdpy11::DiploidMetadata& /*offspring_metadata*/,
          const std::size_t /*parent1*/, const std::size_t /*parent2*/,
          const fwdpy11::SlocusPop& /*pop*/) const
    {
        return 0.0;
    }

    pybind11::object
    pickle() const
    {
        return pybind11::bytes("custom_stateless_genotype");
    }

    DEFAULT_SLOCUSPOP_UPDATE();
};

//Standard pybind11 stuff goes here
PYBIND11_MODULE(custom_stateless_genotype, m)
{
    pybind11::object imported_custom_stateless_genotype_base_class_type
        = pybind11::module::import("fwdpy11.genetic_values")
              .attr("SlocusPopGeneticValue");
    pybind11::class_<GeneralW, fwdpy11::SlocusPopGeneticValue>(m, "GeneralW")
        .def(pybind11::init<>())
        .def(pybind11::pickle([](const GeneralW& g) { return g.pickle(); },
                              [](pybind11::object o) {
                                  pybind11::bytes b(o);
                                  auto s = b.cast<std::string>();
                                  if (s != "custom_stateless_genotype")
                                      {
                                          throw std::runtime_error(
                                              "invalid object state");
                                      }
                                  return GeneralW();
                              }));
}
