#pragma once

#include "fwdpp/gsl_discrete.hpp"
#include "fwdpy11/regions/Region.hpp"
#include "fwdpy11/rng.hpp"
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <fwdpy11/regions/RecombinationRegions.hpp>
#include <gsl/gsl_rng.h>

namespace fwdpy11_core
{
    struct PoissonInterval : public fwdpy11::PoissonCrossoverGenerator
    {
        double left_boundary;
        double right_boundary;
        double mean;
        bool discrete;
        PoissonInterval(double, double, double, bool);
        virtual double
        left()
        {
            return left_boundary;
        }
        virtual double
        right()
        {
            return right_boundary;
        }
        virtual double
        mean_number_xovers()
        {
            return mean;
        }
        virtual void
        breakpoint(const fwdpy11::GSLrng_t& rng, std::vector<double>& breakpoints)
        {
            auto bp = gsl_ran_flat(rng.get(), left_boundary, right_boundary);
            if (discrete)
                {
                    bp = std::floor(bp);
                }
            breakpoints.push_back(bp);
        }
        virtual std::unique_ptr<fwdpy11::PoissonCrossoverGenerator> ll_clone();
    };

    struct PoissonPoint : public fwdpy11::PoissonCrossoverGenerator
    {
        double position;
        double mean;
        bool discrete;
        PoissonPoint(double, double, bool);
        virtual double
        left()
        {
            return position;
        }
        virtual double
        right()
        {
            return position + 1.;
        }
        virtual double
        mean_number_xovers()
        {
            return mean;
        }
        virtual void
        breakpoint(const fwdpy11::GSLrng_t& /*rng*/, std::vector<double>& breakpoints)
        {
            breakpoints.push_back(position);
        }
        virtual std::unique_ptr<fwdpy11::PoissonCrossoverGenerator> ll_clone();
    };

    struct BinomialInterval : public fwdpy11::NonPoissonCrossoverGenerator
    {
        double left_boundary;
        double right_boundary;
        double probability;
        bool discrete;
        BinomialInterval(double, double, double, bool);
        virtual double
        left()
        {
            return left_boundary;
        }
        virtual double
        right()
        {
            return right_boundary;
        }
        virtual void
        breakpoint(const fwdpy11::GSLrng_t& rng, std::vector<double>& breakpoints)
        {
            if (gsl_rng_uniform(rng.get()) <= probability)
                {
                    auto bp = gsl_ran_flat(rng.get(), left_boundary, right_boundary);
                    if (discrete)
                        {
                            bp = std::floor(bp);
                        }
                    breakpoints.push_back(bp);
                }
        }
        virtual std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator> ll_clone();
    };

    struct BinomialIntervalMap : public fwdpy11::NonPoissonCrossoverGenerator
    {
        struct Segment
        {
            double left;
            double right;
        };
        double probability;
        bool discrete;
        fwdpp::gsl_ran_discrete_t_ptr lookup;
        std::vector<fwdpy11::Region> regions;
        std::vector<Segment> segments;
        BinomialIntervalMap(double probability, bool discrete,
                            const std::vector<fwdpy11::Region>&);

        virtual double left();
        virtual double right();
        virtual void breakpoint(const fwdpy11::GSLrng_t& rng,
                                std::vector<double>& breakpoints);
        virtual std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator> ll_clone();
    };

    struct BinomialPoint : public fwdpy11::NonPoissonCrossoverGenerator
    {
        double position;
        double probability;
        bool discrete;
        BinomialPoint(double, double, bool);
        virtual double
        left()
        {
            return position;
        }
        virtual double
        right()
        {
            return position + 1.;
        }
        virtual void
        breakpoint(const fwdpy11::GSLrng_t& rng, std::vector<double>& breakpoints)
        {
            if (gsl_rng_uniform(rng.get()) <= probability)
                {
                    breakpoints.push_back(position);
                }
        }
        virtual std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator> ll_clone();
    };

    struct FixedCrossovers : public fwdpy11::NonPoissonCrossoverGenerator
    {
        double left_boundary;
        double right_boundary;
        int number;
        bool discrete;
        FixedCrossovers(double, double, int, bool);
        virtual double
        left()
        {
            return left_boundary;
        }
        virtual double
        right()
        {
            return right_boundary;
        }
        virtual void
        breakpoint(const fwdpy11::GSLrng_t& rng, std::vector<double>& breakpoints)
        {
            for (int i = 0; i < number; ++i)
                {
                    auto bp = gsl_ran_flat(rng.get(), left_boundary, right_boundary);
                    if (discrete)
                        {
                            bp = std::floor(bp);
                        }
                    breakpoints.push_back(bp);
                }
        }
        virtual std::unique_ptr<fwdpy11::NonPoissonCrossoverGenerator> ll_clone();
    };
}
