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
#include <exception>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <fwdpy11/rng.hpp>
#include <fwdpp/extensions/callbacks.hpp>
#include <fwdpp/extensions/regions.hpp>

namespace py = pybind11;

using dfe_callback_type = std::function<double(const gsl_rng *)>;

template <typename S> struct make_sh_model_fixed_dom
{
    fwdpp::extensions::shmodel
    operator()(const double h, const fwdpp::uint_t scaling,
               const double param) const
    {
        S S_(param);
        auto h_ = fwdpp::extensions::constant(h);
        if (scaling == 1)
            {
                return fwdpp::extensions::shmodel{ std::move(S_),
                                                   std::move(h_) };
            }
        auto bound_S = [scaling, S_](const gsl_rng *r) {
            return S_(r) / static_cast<double>(scaling);
        };
        return fwdpp::extensions::shmodel{ std::move(bound_S),
                                           std::move(h_) };
    }

    fwdpp::extensions::shmodel
    operator()(const double h, const fwdpp::uint_t scaling,
               const double param1, const double param2) const
    {
        S S_(param1, param2);
        auto h_ = fwdpp::extensions::constant(h);
        if (scaling == 1)
            {
                return fwdpp::extensions::shmodel{ std::move(S_),
                                                   std::move(h_) };
            }
        auto bound_S = [scaling, S_](const gsl_rng *r) {
            return S_(r) / static_cast<double>(scaling);
        };
        return fwdpp::extensions::shmodel{ std::move(bound_S),
                                           std::move(h_) };
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
          ([](const double s, const double h, const fwdpp::uint_t scaling) {
              RETURN_DFE_FIXEDH(constant, scaling, s, h);
          }));

    m.def("makeExpSH",
          ([](const double mean, const double h, const fwdpp::uint_t scaling) {
              RETURN_DFE_FIXEDH(exponential, scaling, mean, h);
          }));

    m.def("makeGaussianSH",
          ([](const double sd, const double h, const fwdpp::uint_t scaling) {
              RETURN_DFE_FIXEDH(gaussian, scaling, sd, h);
          }));

    m.def("makeUniformSH", ([](const double lo, const double hi,
                               const double h, const fwdpp::uint_t scaling) {
              RETURN_DFE2_FIXEDH(uniform, scaling, lo, hi, h);
          }));

    m.def("makeGammaSH", ([](const double mean, const double shape,
                             const double h, const fwdpp::uint_t scaling) {
              RETURN_DFE2_FIXEDH(gamma, scaling, mean, shape, h);
          }));

    py::class_<fwdpp::extensions::discrete_mut_model>(m, "MutationRegions")
        .def(py::init<std::vector<double>, std::vector<double>,
                      std::vector<double>, std::vector<double>,
                      std::vector<double>, std::vector<double>,
                      std::vector<fwdpp::extensions::shmodel>,
                      std::vector<decltype(fwdpp::mutation::xtra)>,
                      std::vector<decltype(fwdpp::mutation::xtra)>>());

    py::class_<fwdpp::extensions::discrete_rec_model>(m,
                                                      "RecombinationRegions")
        .def(py::init([](const fwdpy11::GSLrng_t &r, const double recrate,
                         std::vector<double> beg, std::vector<double> end,
                         std::vector<double> weights) {
            return std::unique_ptr<fwdpp::extensions::discrete_rec_model>(
                new fwdpp::extensions::discrete_rec_model(
                    r.get(), recrate, std::move(beg), std::move(end),
                    std::move(weights)));
        }));
}
