#include "fwdpy11/regions/RecombinationRegions.hpp"
#include <cmath>
#include <gsl/gsl_randist.h>
#include <limits>
#include <memory>
#include <stdexcept>
#include <core/genetic_maps/regions.hpp>
#include <type_traits>

static void
validate_position(double pos)
{
    if (pos < 0.0)
        {
            throw std::invalid_argument("positions must be >= 0.0");
        }
    if (!std::isfinite(pos))
        {
            throw std::invalid_argument("positions must be finite");
        }
}

static void
validate_interval(double left, double right, bool discrete)
{
    if (left >= right)
        {
            throw std::invalid_argument("left must be < right");
        }
    if (discrete == true && right - left <= 1.0)
        {
            throw std::invalid_argument(
                "interval length must be > 1 when discrete == True");
        }
}

static void
validate_parameter(double param, double minimum, double maximum)
{
    if (!std::isfinite(param))
        {
            throw std::invalid_argument("parameter values must be finite");
        }
    if (param < minimum || param > maximum)
        {
            throw std::invalid_argument("parameter value is invalid");
        }
}

namespace fwdpy11_core
{
    PoissonInterval::PoissonInterval(double left, double right, double mean,
                                     bool discrete)
        : left_boundary(left), right_boundary(right), mean(mean), discrete(discrete)
    {
        validate_position(left);
        validate_position(right);
        validate_interval(left, right, discrete);
        validate_parameter(mean, 0.0, std::numeric_limits<double>::max());
    }

    std::unique_ptr<fwdpy11::PoissonCrossoverGenerator>
    PoissonInterval::ll_clone()
    {
        return std::unique_ptr<fwdpy11::PoissonCrossoverGenerator>(
            new PoissonInterval(*this));
    }

    PoissonPoint::PoissonPoint(double position, double mean, bool discrete)
        : position(position), mean(mean), discrete(discrete)
    {
        validate_position(position);
        validate_parameter(mean, 0.0, std::numeric_limits<double>::max());
    }

    std::unique_ptr<fwdpy11::PoissonCrossoverGenerator>
    PoissonPoint::ll_clone()
    {
        return std::unique_ptr<fwdpy11::PoissonCrossoverGenerator>(
            new PoissonPoint(*this));
    }

    BinomialInterval::BinomialInterval(double left, double right, double probability,
                                       bool discrete)
        : left_boundary(left), right_boundary(right), probability(probability),
          discrete(discrete)
    {
        validate_position(left);
        validate_position(right);
        validate_interval(left, right, discrete);
        validate_parameter(probability, 0.0, 1.0);
    }

    std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator>
    BinomialInterval::ll_clone()
    {
        return std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator>(
            new BinomialInterval(*this));
    }

    BinomialPoint::BinomialPoint(double position, double probability, bool discrete)
        : position(position), probability(probability), discrete(discrete)
    {
        validate_position(position);
        validate_parameter(probability, 0.0, 1.0);
    }

    std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator>
    BinomialPoint::ll_clone()
    {
        return std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator>(
            new BinomialPoint(*this));
    }

    FixedCrossovers::FixedCrossovers(double left, double right, int number,
                                     bool discrete)
        : left_boundary(left), right_boundary(right), number(number), discrete(discrete)
    {
        validate_position(left);
        validate_position(right);
        validate_interval(left, right, discrete);
        if (number <= 0)
            {
                throw std::invalid_argument("nummber of crossovers must be >= 0");
            }
    }

    std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator>
    FixedCrossovers::ll_clone()
    {
        return std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator>(
            new FixedCrossovers(*this));
    }

    BinomialIntervalMap::BinomialIntervalMap(double probability, const bool discrete,
                                             const std::vector<fwdpy11::Region>& regions)
        : probability(probability), discrete{discrete}, lookup{nullptr},
          regions{regions}, segments{}
    {
        validate_parameter(probability, 0.0, 1.0);
        std::vector<double> weights;
        weights.reserve(regions.size());
        segments.reserve(regions.size());
        for (const auto& r : this->regions)
            {
                validate_interval(r.beg, r.end, discrete);
                auto w = r.weight;
                if (r.coupled)
                    {
                        w += (r.end - r.beg);
                    }
                weights.push_back(w);
                segments.push_back(Segment{r.beg, r.end});
            }
        lookup.reset(gsl_ran_discrete_preproc(segments.size(), weights.data()));
    }

    double
    BinomialIntervalMap::left()
    {
        auto rv = std::numeric_limits<double>::max();
        for (const auto& s : this->segments)
            {
                rv = std::min(s.left, rv);
            }
        return rv;
    }

    double
    BinomialIntervalMap::right()
    {
        auto rv = std::numeric_limits<double>::min();
        for (const auto& s : this->segments)
            {
                rv = std::max(s.right, rv);
            }
        return rv;
    }

    void
    BinomialIntervalMap::breakpoint(const fwdpy11::GSLrng_t& rng,
                                    std::vector<double>& breakpoints)
    {
        if (gsl_rng_uniform(rng.get()) <= probability)
            {
                auto index = gsl_ran_discrete(rng.get(), this->lookup.get());
                auto bp = gsl_ran_flat(rng.get(), this->segments[index].left,
                                       this->segments[index].right);
                if (discrete)
                    {
                        bp = std::floor(bp);
                    }
                breakpoints.push_back(bp);
            }
    }

    std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator>
    BinomialIntervalMap::ll_clone()
    {
        return std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator>(
            new BinomialIntervalMap(this->probability, this->discrete, this->regions));
    }
}
