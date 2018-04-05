//
// Copyright (C) 2017 Kevin Thornton <krthornt@uci.edu>
//
// This file is part of fwdpy11.
//
// fwdpy11 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fwdpy11 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
//
#include <cmath>
#include <cstdint>
#include <exception>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <fwdpy11/types/SlocusPop.hpp>
#include <fwdpy11/types/MlocusPop.hpp>
#include <fwdpy11/rng.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/extensions/regions.hpp>

namespace py = pybind11;

using dfe_callback_type = std::function<double(const gsl_rng *)>;

template <typename S> struct make_sh_model_fixed_dom
{
    fwdpp::extensions::shmodel
    operator()(const double h, const std::int64_t scaling,
               const double param) const
    {
        S S_(param);
        auto h_ = fwdpp::extensions::constant(h);
        if (scaling == 1)
            {
                return fwdpp::extensions::shmodel{ std::move(S_),
                                                   std::move(h_) };
            }
        if (scaling < 1)
            {
                throw std::invalid_argument("DES scaling must be >= 1");
            }
        auto bound_S = [scaling, S_](const gsl_rng *r) {
            return S_(r) / static_cast<double>(scaling);
        };
        return fwdpp::extensions::shmodel{ std::move(bound_S), std::move(h_) };
    }

    fwdpp::extensions::shmodel
    operator()(const double h, const std::int64_t scaling, const double param1,
               const double param2) const
    {
        S S_(param1, param2);
        auto h_ = fwdpp::extensions::constant(h);
        if (scaling == 1)
            {
                return fwdpp::extensions::shmodel{ std::move(S_),
                                                   std::move(h_) };
            }
        if (scaling < 1)
            {
                throw std::invalid_argument("DES scaling must be >= 1");
            }
        auto bound_S = [scaling, S_](const gsl_rng *r) {
            return S_(r) / static_cast<double>(scaling);
        };
        return fwdpp::extensions::shmodel{ std::move(bound_S), std::move(h_) };
    }
};

#define RETURN_DFE_FIXEDH(DIST, SCALING, S, H)                                \
    return make_sh_model_fixed_dom<fwdpp::extensions::DIST>()(H, SCALING, S);
#define RETURN_DFE2_FIXEDH(DIST, SCALING, A, B, H)                            \
    return make_sh_model_fixed_dom<fwdpp::extensions::DIST>()(H, SCALING, A,  \
                                                              B);

PYBIND11_MODULE(fwdpp_extensions, m)
{
    m.doc() = "Expose fwdpp's extensions library.";

    py::class_<fwdpp::extensions::shmodel>(m, "DFEFixedDominance")
        .def(py::init<>())
        .def(py::init<dfe_callback_type, dfe_callback_type>())
        .def("__call__", [](const fwdpp::extensions::shmodel &sh,
                            const fwdpy11::GSLrng_t &rng) {
            return py::make_tuple(sh.s(rng.get()), sh.h(rng.get()));
        });

    m.def("makeConstantSH",
          ([](const double s, const double h, const std::int64_t scaling) {
              RETURN_DFE_FIXEDH(constant, scaling, s, h);
          }));

    m.def("makeExpSH",
          ([](const double mean, const double h, const std::int64_t scaling) {
              RETURN_DFE_FIXEDH(exponential, scaling, mean, h);
          }));

    m.def("makeGaussianSH",
          ([](const double sd, const double h, const std::int64_t scaling) {
              RETURN_DFE_FIXEDH(gaussian, scaling, sd, h);
          }));

    m.def("makeUniformSH", ([](const double lo, const double hi,
                               const double h, const std::int64_t scaling) {
              RETURN_DFE2_FIXEDH(uniform, scaling, lo, hi, h);
          }));

    m.def("makeGammaSH", ([](const double mean, const double shape,
                             const double h, const std::int64_t scaling) {
              RETURN_DFE2_FIXEDH(gamma, scaling, mean, shape, h);
          }));

    // start, stop, weight, label
    using mutregion_tuple = std::tuple<double, double, double,
                                       decltype(fwdpp::mutation_base::xtra)>;
    using mutfunction = fwdpp::extensions::discrete_mut_model<
        std::vector<fwdpp::popgenmut>>::function_type;

    //This is a hack-ish means of doing things.
    //Once we have a proper base class in place, 
    //we can reduce what is below to a single function
    //taking const ref to base.
    py::class_<
        fwdpp::extensions::discrete_mut_model<std::vector<fwdpp::popgenmut>>>(
        m, "MutationRegions")
        .def(py::init(
            [](const fwdpy11::GSLrng_t &r, fwdpy11::SlocusPop &pop,
               const std::vector<mutregion_tuple> &neutral_mut_regions,
               const std::vector<mutregion_tuple> &selected_mut_regions,
               std::vector<fwdpp::extensions::shmodel> &dist_effect_sizes) {
                std::vector<mutfunction> functions;
                std::vector<double> weights;
                for (auto &&i : neutral_mut_regions)
                    {
                        auto start = std::get<0>(i);
                        auto stop = std::get<1>(i);
                        auto label = std::get<3>(i);
                        weights.push_back(std::get<2>(i));
                        auto mf
                            = [&r, &pop, start, stop, label](
                                  std::queue<std::size_t> &recbin,
                                  fwdpy11::SlocusPop::mcont_t &mutations) {
                                  return fwdpp::infsites_popgenmut(
                                      recbin, mutations, r.get(),
                                      pop.mut_lookup, pop.generation, 0,
                                      [&r, start, stop]() {
                                          return gsl_ran_flat(r.get(), start,
                                                              stop);
                                      },
                                      []() { return 0.; }, []() { return 1.; },
                                      label);
                              };
                        functions.push_back(std::move(mf));
                    }
                for (std::size_t i = 0; i < selected_mut_regions.size(); ++i)
                    {
                        auto start = std::get<0>(selected_mut_regions.at(i));
                        auto stop = std::get<1>(selected_mut_regions.at(i));
                        auto label = std::get<3>(selected_mut_regions.at(i));
                        auto sh = dist_effect_sizes.at(i);
                        weights.push_back(
                            std::get<2>(selected_mut_regions.at(i)));
                        auto mf
                            = [&r, &pop, start, stop, label,
                               sh](std::queue<std::size_t> &recbin,
                                   fwdpy11::SlocusPop::mcont_t &mutations) {
                                  return fwdpp::infsites_popgenmut(
                                      recbin, mutations, r.get(),
                                      pop.mut_lookup, pop.generation, 1.,
                                      [&r, start, stop]() {
                                          return gsl_ran_flat(r.get(), start,
                                                              stop);
                                      },
                                      [&r, sh]() { return sh.s(r.get()); },
                                      [&r, sh]() { return sh.h(r.get()); },
                                      label);
                              };
                        functions.push_back(std::move(mf));
                    }
                return fwdpp::extensions::discrete_mut_model<
                    fwdpy11::SlocusPop::mcont_t>(std::move(functions),
                                                   std::move(weights));
            }))
        .def(py::init(
            [](const fwdpy11::GSLrng_t &r, fwdpy11::MlocusPop &pop,
               const std::vector<mutregion_tuple> &neutral_mut_regions,
               const std::vector<mutregion_tuple> &selected_mut_regions,
               std::vector<fwdpp::extensions::shmodel> &dist_effect_sizes) {
                std::vector<mutfunction> functions;
                std::vector<double> weights;
                for (auto &&i : neutral_mut_regions)
                    {
                        auto start = std::get<0>(i);
                        auto stop = std::get<1>(i);
                        auto label = std::get<3>(i);
                        weights.push_back(std::get<2>(i));
                        auto mf
                            = [&r, &pop, start, stop, label](
                                  std::queue<std::size_t> &recbin,
                                  fwdpy11::MlocusPop::mcont_t &mutations) {
                                  return fwdpp::infsites_popgenmut(
                                      recbin, mutations, r.get(),
                                      pop.mut_lookup, pop.generation, 0,
                                      [&r, start, stop]() {
                                          return gsl_ran_flat(r.get(), start,
                                                              stop);
                                      },
                                      []() { return 0.; }, []() { return 1.; },
                                      label);
                              };
                        functions.push_back(std::move(mf));
                    }
                for (std::size_t i = 0; i < selected_mut_regions.size(); ++i)
                    {
                        auto start = std::get<0>(selected_mut_regions.at(i));
                        auto stop = std::get<1>(selected_mut_regions.at(i));
                        auto label = std::get<2>(selected_mut_regions.at(i));
                        auto sh = dist_effect_sizes.at(i);
                        weights.push_back(
                            std::get<2>(selected_mut_regions.at(i)));
                        auto mf
                            = [&r, &pop, start, stop, label,
                               sh](std::queue<std::size_t> &recbin,
                                   fwdpy11::MlocusPop::mcont_t &mutations) {
                                  return fwdpp::infsites_popgenmut(
                                      recbin, mutations, r.get(),
                                      pop.mut_lookup, pop.generation, 1.,
                                      [&r, start, stop]() {
                                          return gsl_ran_flat(r.get(), start,
                                                              stop);
                                      },
                                      [&r, sh]() { return sh.s(r.get()); },
                                      [&r, sh]() { return sh.h(r.get()); },
                                      label);
                              };
                        functions.push_back(std::move(mf));
                    }
                return fwdpp::extensions::discrete_mut_model<
                    fwdpy11::MlocusPop::mcont_t>(std::move(functions),
                                                    std::move(weights));
            }));

    py::class_<fwdpp::extensions::discrete_rec_model>(m,
                                                      "RecombinationRegions")
        .def(py::init([](const fwdpy11::GSLrng_t &r, const double recrate,
                         std::vector<double> beg, std::vector<double> end,
                         std::vector<double> weights) {
            std::vector<fwdpp::extensions::discrete_rec_model::function_type>
                functions;
            for (std::size_t i = 0; i < beg.size(); ++i)
                {
                    auto start = beg.at(i);
                    auto stop = end.at(i);
                    auto f = [&r, start, stop](std::vector<double> &b) {
                        b.push_back(gsl_ran_flat(r.get(), start, stop));
                    };
                    functions.push_back(std::move(f));
                }
            return fwdpp::extensions::discrete_rec_model(
                recrate, std::move(functions), std::move(weights));
        }));
}
