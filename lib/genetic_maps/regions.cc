#include "fwdpy11/regions/RecombinationRegions.hpp"
#include <cmath>
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
}
