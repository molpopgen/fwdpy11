#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <fwdpp/gsl_discrete.hpp>
#include <fwdpy11/regions/Sregion.hpp>
#include <fwdpy11/policies/mutation.hpp>

namespace py = pybind11;

class DiscreteDESD : public fwdpy11::Sregion
{
  private:
    std::vector<double> esize, h, weight;
    fwdpp::gsl_ran_discrete_t_ptr sh_lookup;

    fwdpp::gsl_ran_discrete_t_ptr
    init_lookup()
    {
        for (auto wi : weight)
            {
                if (wi < 0)
                    {
                        throw std::invalid_argument("all weights must be >= 0");
                    }
                if (!std::isfinite(wi))
                    {
                        throw std::invalid_argument("all weights must be finite");
                    }
            }
        return fwdpp::gsl_ran_discrete_t_ptr(
            gsl_ran_discrete_preproc(weight.size(), weight.data()));
    }

  public:
    DiscreteDESD(const fwdpy11::Region& r, const double sc, std::vector<double> esize_,
                 std::vector<double> h_, std::vector<double> w_)
        : fwdpy11::Sregion(r, sc, 1, fwdpy11::process_input_dominance(0.)),
          esize{std::move(esize_)}, h{std::move(h_)}, weight{std::move(w_)},
          sh_lookup{init_lookup()}
    {
        if (esize.size() != h.size() || esize.size() != weight.size()
            || h.size() != weight.size())
            {
                throw std::invalid_argument("all arrays must be equal-length");
            }
    }

    std::unique_ptr<fwdpy11::Sregion>
    clone() const override
    {
        return std::make_unique<DiscreteDESD>(this->region, this->scaling, this->esize,
                                              this->h, this->weight);
    }

    std::uint32_t
    operator()(fwdpp::flagged_mutation_queue& recycling_bin,
               std::vector<fwdpy11::Mutation>& mutations,
               std::unordered_multimap<double, std::uint32_t>& lookup_table,
               const std::uint32_t generation,
               const fwdpy11::GSLrng_t& rng) const override
    {
        auto idx = gsl_ran_discrete(rng.get(), sh_lookup.get());
        return fwdpy11::infsites_Mutation(
            recycling_bin, mutations, lookup_table, false, generation,
            [this, &rng]() { return region(rng); },
            [this, idx]() { return esize[idx] / scaling; },
            [this, idx](const double /*esize*/) { return h[idx]; }, this->label());
    }

    double
    from_mvnorm(const double /*deviate*/, const double /*P*/) const override
    {
        throw std::runtime_error("not implemented yet");
        return 1.;
    }

    double
    generate_dominance(const fwdpy11::GSLrng_t& /*rng*/,
                       const double /*esize*/) const override
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
};

void
init_DiscreteDESD(py::module& m)
{
    py::class_<DiscreteDESD, fwdpy11::Sregion>(m, "_ll_DiscreteDESD")
        .def(py::init([](double beg, double end, double weight,
                         std::vector<std::tuple<double, double, double>> joint_dist,
                         bool coupled, std::uint16_t label, double scaling) {
                 if (joint_dist.empty())
                     {
                         throw std::invalid_argument("empty input list");
                     }
                 std::vector<double> esizes, dominance, weights;
                 for (auto&& t : joint_dist)
                     {
                         auto e = std::get<0>(t);
                         auto h = std::get<1>(t);
                         auto w = std::get<2>(t);
                         if (!std::isfinite(e))
                             {
                                 throw std::invalid_argument(
                                     "non-finite effect size input");
                             }
                         if (!std::isfinite(h))
                             {
                                 throw std::invalid_argument(
                                     "non-finite dominance input");
                             }
                         if (!std::isfinite(w))
                             {
                                 throw std::invalid_argument("non-finite weight input");
                             }
                         if (w < 0.0)
                             {
                                 throw std::invalid_argument("weights must be >= 0.0");
                             }
                         esizes.push_back(e);
                         dominance.push_back(h);
                         weights.push_back(w);
                     }
                 return DiscreteDESD(fwdpy11::Region(beg, end, weight, coupled, label),
                                     scaling, std::move(esizes), std::move(dominance),
                                     std::move(weights));
             }),
             py::arg("beg"), py::arg("end"), py::arg("weight"), py::arg("joint_dist"),
             py::arg("coupled"), py::arg("label"), py::arg("scaling"));
}

