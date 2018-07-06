//Challenge here is in the aggregation:
//Addtive trait
//Addtive fitness
//Mult trait
//Mult fitness
//The above 4 need to be auto-detected at construction time
template <typename fwdppT, typename aggregator>
struct wrap_fwdpp_genetic_value_mlocus
    : public fwdpy11::MlocusPopGeneticValueWithMapping
{
    using gvalue_map_ptr = std::unique_ptr<fwdpy11::GeneticValueToFitnessMap>;
    const fwdppT gv;
    aggregator agg;
    const double starting_value;

    wrap_fwdpp_genetic_value_mlocus(const double);

    wrap_fwdpp_genetic_value_mlocus(const double scaling,
                                    const fwdpy11::GeneticValueIsTrait& g2w);

    wrap_fwdpp_genetic_value_mlocus(
        const double scaling, const fwdpy11::GeneticValueIsTrait& g2w,
        const fwdpy11::GeneticValueNoise& noise_fxn);

    inline double
    operator()(const std::size_t diploid_index,
               const fwdpy11::MlocusPop& pop) const
    {
        return gv(pop.diploids[diploid_index], pop.gametes, pop.mutations);
    }

    inline void
    update(const fwdpy11::MlocusPop& pop)
    {
        gv2w->update(pop);
        noise_fxn->update(pop);
    }
};

template <>
wrap_fwdpp_genetic_value_mlocus<fwdpp::additive_diploid,std::plus<double>>::wrap_fwdpp_genetic_value_mlocus(
    const double scaling)
    : fwdpy11::MlocusPopGeneticValueWithMapping{ gvalue_map_ptr(
          new fwdpy11::GeneticValueIsFitness()) },
      gv{ scaling, fwdpp::additive_diploid::policy::aw },agg{}
{
    if (!std::isfinite(scaling))
        {
            throw std::invalid_argument("scaling must be finite");
        }
}
