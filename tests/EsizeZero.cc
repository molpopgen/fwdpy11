#include <pybind11/pybind11.h>
#include <fwdpy11/regions/Sregion.hpp>
#include <fwdpy11/policies/mutation.hpp>

using namespace pybind11::literals;

struct EsizeZero : public fwdpy11::Sregion
{
    EsizeZero(const fwdpy11::Region& r) : fwdpy11::Sregion(r, 1., 1)
    {
    }

    std::unique_ptr<fwdpy11::Sregion>
    clone() const override
    {
        return std::unique_ptr<EsizeZero>(new EsizeZero(*this));
    }

    std::string
    repr() const override
    {
        return "EsizeZero()"_s;
    }

    pybind11::tuple
    pickle() const override
    {
        return pybind11::make_tuple(Sregion::pickle_Sregion());
    }

    static EsizeZero
    unpickle(pybind11::tuple t)
    {
        return EsizeZero(fwdpy11::Region::unpickle(t[0]));
    }

    double
    from_mvnorm(const double, const double) const override
    {
        return 0.0;
    }

    std::vector<double>
    get_dominance() const override
    {
        return {1.};
    }

    std::uint32_t
    operator()(fwdpp::flagged_mutation_queue& recycling_bin,
               std::vector<fwdpy11::Mutation>& mutations,
               std::unordered_multimap<double, std::uint32_t>& lookup_table,
               const std::uint32_t generation,
               const fwdpy11::GSLrng_t& rng) const override
    {
        return fwdpy11::infsites_Mutation(
            recycling_bin, mutations, lookup_table, false, generation,
            [this, &rng]() { return region(rng); }, []() { return 0.; },
            []() { return 1.; }, this->label());
    }
};

PYBIND11_MODULE(EsizeZero, m)
{
    pybind11::object base = pybind11::module::import("fwdpy11").attr("Sregion");
    pybind11::class_<EsizeZero, fwdpy11::Sregion>(m, "EsizeZero")
        .def(pybind11::init([](double beg, double end, double weight, bool coupled,
                               std::uint16_t label) {
            return EsizeZero(fwdpy11::Region(beg, end, weight, coupled, label));
        }));
}

