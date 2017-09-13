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
    KTfwd::extensions::shmodel
    operator()(const double h, const KTfwd::uint_t scaling,
               const double param) const
    {
        if (scaling == 1)
            {
                return KTfwd::extensions::shmodel(
                    std::bind(S(param), std::placeholders::_1),
                    std::bind(KTfwd::extensions::constant(h),
                              std::placeholders::_1));
            }
        auto scaling_fxn = [scaling, param](const gsl_rng *r) {
            return S(param)(r) / static_cast<double>(scaling);
        };
        return KTfwd::extensions::shmodel(
            std::move(scaling_fxn),
            std::bind(KTfwd::extensions::constant(h), std::placeholders::_1));
    }

    KTfwd::extensions::shmodel
    operator()(const double h, const KTfwd::uint_t scaling,
               const double param1, const double param2) const
    {
        if (scaling == 1)
            {
                return KTfwd::extensions::shmodel(
                    std::bind(S(param1, param2), std::placeholders::_1),
                    std::bind(KTfwd::extensions::constant(h),
                              std::placeholders::_1));
            }
        auto scaling_fxn = [scaling, param1, param2](const gsl_rng *r) {
            return S(param1, param2)(r) / static_cast<double>(scaling);
        };
        return KTfwd::extensions::shmodel(
            std::move(scaling_fxn),
            std::bind(KTfwd::extensions::constant(h), std::placeholders::_1));
    }
};

#define RETURN_DFE_FIXEDH(DIST, SCALING, S, H)                                \
    return make_sh_model_fixed_dom<KTfwd::extensions::DIST>()(H, SCALING, S);
#define RETURN_DFE2_FIXEDH(DIST, SCALING, A, B, H)                            \
    return make_sh_model_fixed_dom<KTfwd::extensions::DIST>()(H, SCALING, A,  \
                                                              B);

PYBIND11_MODULE(fwdpp_extensions,m)
{
    m.doc() = "Expose fwdpp's extensions library.";

    py::class_<KTfwd::extensions::shmodel>(m, "DFEFixedDominance")
        .def(py::init<>())
        .def(py::init<dfe_callback_type, dfe_callback_type>())
        .def("__call__", [](const KTfwd::extensions::shmodel &sh,
                            const fwdpy11::GSLrng_t &rng) {
            return py::make_tuple(sh.s(rng.get()), sh.h(rng.get()));
        });

    m.def("makeConstantSH",
          ([](const double s, const double h, const KTfwd::uint_t scaling) {
              RETURN_DFE_FIXEDH(constant, scaling, s, h);
          }));

    m.def("makeExpSH",
          ([](const double mean, const double h, const KTfwd::uint_t scaling) {
              RETURN_DFE_FIXEDH(exponential, scaling, mean, h);
          }));

    m.def("makeGaussianSH",
          ([](const double sd, const double h, const KTfwd::uint_t scaling) {
              RETURN_DFE_FIXEDH(gaussian, scaling, sd, h);
          }));

    m.def("makeUniformSH", ([](const double lo, const double hi,
                               const double h, const KTfwd::uint_t scaling) {
              RETURN_DFE2_FIXEDH(uniform, scaling, lo, hi, h);
          }));

    m.def("makeGammaSH", ([](const double mean, const double shape,
                             const double h, const KTfwd::uint_t scaling) {
              RETURN_DFE2_FIXEDH(gamma, scaling, mean, shape, h);
          }));

    py::class_<KTfwd::extensions::discrete_mut_model>(m, "MutationRegions")
        .def(py::init<std::vector<double>, std::vector<double>,
                      std::vector<double>, std::vector<double>,
                      std::vector<double>, std::vector<double>,
                      std::vector<KTfwd::extensions::shmodel>>());

    py::class_<KTfwd::extensions::discrete_rec_model>(m,
                                                      "RecombinationRegions")
        .def(py::init<std::vector<double>, std::vector<double>,
                      std::vector<double>>());
}
