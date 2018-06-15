template <typename fwdppT>
struct wrap_fwdpp_genetic_value
    : public fwdpy11::SlocusPopGeneticValueWithMapping
{
    using gvalue_map_ptr = std::unique_ptr<fwdpy11::GeneticValueToFitnessMap>;
    const fwdppT gv;

    wrap_fwdpp_genetic_value(const double);

    wrap_fwdpp_genetic_value(const double scaling,
                             const fwdpy11::GeneticValueIsTrait& g2w);

    wrap_fwdpp_genetic_value(const double scaling,
                             const fwdpy11::GeneticValueIsTrait& g2w,
                             const fwdpy11::GeneticValueNoise& noise_fxn);

    inline double
    operator()(const std::size_t diploid_index,
               const fwdpy11::SlocusPop& pop) const
    {
        return gv(pop.diploids[diploid_index], pop.gametes, pop.mutations);
    }

    inline void
    update(const fwdpy11::SlocusPop& pop)
    {
        gv2w->update(pop);
        noise_fxn->update(pop);
    }
};

template <>
wrap_fwdpp_genetic_value<fwdpp::additive_diploid>::wrap_fwdpp_genetic_value(
    const double scaling)
    : fwdpy11::SlocusPopGeneticValueWithMapping{ gvalue_map_ptr(
          new fwdpy11::GeneticValueIsFitness()) },
      gv{ scaling, fwdpp::additive_diploid::policy::aw }
{
    if (!std::isfinite(scaling))
        {
            throw std::invalid_argument("scaling must be finite");
        }
}

template <>
wrap_fwdpp_genetic_value<fwdpp::multiplicative_diploid>::
    wrap_fwdpp_genetic_value(const double scaling)
    : fwdpy11::SlocusPopGeneticValueWithMapping{ gvalue_map_ptr(
          new fwdpy11::GeneticValueIsFitness()) },
      gv{ scaling, fwdpp::multiplicative_diploid::policy::mw }
{
    if (!std::isfinite(scaling))
        {
            throw std::invalid_argument("scaling must be finite");
        }
}

template <>
wrap_fwdpp_genetic_value<fwdpp::additive_diploid>::wrap_fwdpp_genetic_value(
    const double scaling, const fwdpy11::GeneticValueIsTrait& g2w)
    : fwdpy11::SlocusPopGeneticValueWithMapping{ g2w.clone() }, gv{
          scaling, fwdpp::additive_diploid::policy::atrait
      }
{
}

template <>
wrap_fwdpp_genetic_value<fwdpp::additive_diploid>::wrap_fwdpp_genetic_value(
    const double scaling, const fwdpy11::GeneticValueIsTrait& g2w,
    const fwdpy11::GeneticValueNoise& noise_fxn)
    : fwdpy11::SlocusPopGeneticValueWithMapping{ g2w.clone(),
                                                 noise_fxn.clone() },
      gv{ scaling, fwdpp::additive_diploid::policy::atrait }
{
}

template <>
wrap_fwdpp_genetic_value<fwdpp::multiplicative_diploid>::
    wrap_fwdpp_genetic_value(const double scaling,
                             const fwdpy11::GeneticValueIsTrait& g2w)
    : fwdpy11::SlocusPopGeneticValueWithMapping{ g2w.clone() }, gv{
          scaling, fwdpp::multiplicative_diploid::policy::mtrait
      }
{
}

template <>
wrap_fwdpp_genetic_value<fwdpp::multiplicative_diploid>::
    wrap_fwdpp_genetic_value(const double scaling,
                             const fwdpy11::GeneticValueIsTrait& g2w,
                             const fwdpy11::GeneticValueNoise& noise_fxn)
    : fwdpy11::SlocusPopGeneticValueWithMapping{ g2w.clone(),
                                                 noise_fxn.clone() },
      gv{ scaling, fwdpp::multiplicative_diploid::policy::mtrait }
{
}
