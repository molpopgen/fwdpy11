#ifndef FWDPY11_REGIONS_GENETICMAPUNIT_HPP
#define FWDPY11_REGIONS_GENETICMAPUNIT_HPP

#include <vector>
#include <memory>
#include <fwdpy11/rng.hpp>
#include <pybind11/pybind11.h>

namespace fwdpy11
{
    struct GeneticMapUnit
    {
        virtual ~GeneticMapUnit() = default;
        virtual void operator()(const GSLrng_t&, std::vector<double>&) const = 0;
        virtual pybind11::object pickle()  const= 0;
        virtual std::unique_ptr<GeneticMapUnit> clone()  const = 0;
    };
} // namespace fwdpy11

#endif

