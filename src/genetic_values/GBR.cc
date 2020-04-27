#include <fwdpy11/genetic_values/DiploidGeneticValue.hpp>
#include <fwdpy11/genetic_value_to_fitness/GeneticValueToFitnessMap.hpp>
#include <fwdpy11/genetic_value_noise/GeneticValueNoise.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace
{
    class GBR : public fwdpy11::DiploidGeneticValue
    {
      private:
        mutable std::size_t deme;
        template <typename mutation_key_cont_t, typename mcont_t>
        inline double
        sum_haplotype_effect_sizes(const mutation_key_cont_t& keys,
                                   const mcont_t& mutations) const
        {
            double rv = 0.0;
            for (auto& k : keys)
                {
                    rv += mutations[k].s;
                }
            return rv;
        }

        template <typename mutation_key_cont_t, typename mcont_t>
        inline double
        sum_haplotype_effect_sizes(std::size_t deme, const mutation_key_cont_t& keys,
                                   const mcont_t& mutations) const
        {
            double rv = 0.0;
            for (auto& k : keys)
                {
                    rv += mutations[k].esizes[deme];
                }
            return rv;
        }

        std::function<double(const std::size_t diploid_index,
                             const fwdpy11::DiploidPopulation& pop)>
        generate_backend_function(std::size_t ndemes, std::size_t& deme)
        {
            if (ndemes == 1)
                {
                    return [this](const std::size_t diploid_index,
                                  const fwdpy11::DiploidPopulation& pop) {
                        double h1 = sum_haplotype_effect_sizes(
                            pop.haploid_genomes[pop.diploids[diploid_index].first]
                                .smutations,
                            pop.mutations);
                        double h2 = sum_haplotype_effect_sizes(
                            pop.haploid_genomes[pop.diploids[diploid_index].second]
                                .smutations,
                            pop.mutations);
                        return sqrt(h1 * h2);
                    };
                }
            return [&deme, this](const std::size_t diploid_index,
                                 const fwdpy11::DiploidPopulation& pop) {
                double h1 = sum_haplotype_effect_sizes(
                    deme,
                    pop.haploid_genomes[pop.diploids[diploid_index].first].smutations,
                    pop.mutations);
                double h2 = sum_haplotype_effect_sizes(
                    deme,
                    pop.haploid_genomes[pop.diploids[diploid_index].second].smutations,
                    pop.mutations);
                return sqrt(h1 * h2);
            };
        }

        const std::function<double(const std::size_t diploid_index,
                                   const fwdpy11::DiploidPopulation& pop)>
            f;

      public:
        GBR(std::size_t ndim, const fwdpy11::GeneticValueToFitnessMap& gv2w_)
            : fwdpy11::DiploidGeneticValue{ndim, gv2w_},
              deme(0), f{generate_backend_function(total_dim, deme)}
        {
        }

        GBR(std::size_t ndim, const fwdpy11::GeneticValueToFitnessMap& gv2w_,
            const fwdpy11::GeneticValueNoise& noise_)
            : fwdpy11::DiploidGeneticValue{ndim, gv2w_, noise_},
              deme(0), f{generate_backend_function(total_dim, deme)}
        {
        }

        double
        calculate_gvalue(const std::size_t diploid_index,
                         const fwdpy11::DiploidMetadata& /*metadata*/,
                         const fwdpy11::DiploidPopulation& pop) const override
        {
            deme = pop.diploid_metadata[diploid_index].deme;
            gvalues[0] = f(diploid_index, pop);
            return gvalues[0];
        }

        pybind11::object
        pickle() const override
        {
            return py::int_(total_dim);
        }
    };
}

static const auto GBR_CONSTRUCTOR1 =
    R"delim(
 Construct object with specific genetic value to fitness map.
 
 :param gv2w: Genetic value to fitness map
 :type gv2w: :class:`fwdpy11.GeneticValueIsTrait`
 )delim";

static const auto GBR_CONSTRUCTOR2 =
    R"delim(
Construct object with specific genetic value to fitness map 
and random effects on trait value.

:param gv2w: Genetic value to fitness map
:type gv2w: :class:`fwdpy11.GeneticValueIsTrait`
:param noise: Model of random effects on trait value.
:type noise: :class:`fwdpy11.GeneticValueNoise`
)delim";

void
init_GBR(py::module& m)
{
    py::class_<GBR, fwdpy11::DiploidGeneticValue>(m, "GBR",
                                                  R"delim(
        The "gene-based recessive" trait model described in Thornton et al.
        2013 http://dx.doi.org/10.1371/journal.pgen.1003258 and Sanjak et al. 2017
        http://dx.doi.org/10.1371/journal.pgen.1006573.

        The trait value is the geometric mean of the sum of effect sizes on each haplotype.
        It is undefined for the case where these sums are negative.
        )delim")
        .def(py::init(
                 [](const fwdpy11::GeneticValueIsTrait& gv2w) { return GBR(1, gv2w); }),
             py::arg("gv2w"), GBR_CONSTRUCTOR1)
        .def(py::init([](const fwdpy11::GeneticValueIsTrait& gv2w,
                         const fwdpy11::GeneticValueNoise& noise) {
                 return GBR(1, gv2w, noise);
             }),
             py::arg("gv2w"), py::arg("noise"), GBR_CONSTRUCTOR2)
        .def(py::pickle(
            [](const GBR& g) {
                auto p = py::module::import("pickle");
                return py::make_tuple(g.pickle(), p.attr("dumps")(g.gv2w->clone(), -1),
                                      p.attr("dumps")(g.noise_fxn->clone(), -1));
            },
            [](py::tuple t) {
                if (t.size() != 3)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                auto ndim = t[0].cast<std::size_t>();
                auto p = py::module::import("pickle");
                auto t1 = p.attr("loads")(t[1]);
                auto t2 = p.attr("loads")(t[2]);
                //Do the casts in the constructor
                //to avoid any nasty issues w/
                //refs to temp
                return GBR(ndim, t1.cast<const fwdpy11::GeneticValueIsTrait&>(),
                           t2.cast<const fwdpy11::GeneticValueNoise&>());
            }));
}
