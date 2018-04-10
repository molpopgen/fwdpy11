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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/types/Population.hpp>
#include <fwdpp/sugar/change_neutral.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<fwdpp::uint_t>);
PYBIND11_MAKE_OPAQUE(fwdpy11::Population::mcont_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::Population::gcont_t);

void
check_finite(const double d, const std::string& error)
{
    if (!std::isfinite(d))
        {
            throw std::invalid_argument(error);
        }
}

template <typename mkeys, typename mutation_container>
void
sort_keys(mkeys& k, const mutation_container& mutations)
{
    std::sort(std::begin(k), std::end(k),
              [&mutations](const typename mkeys::value_type a,
                           const typename mkeys::value_type b) {
                  return mutations[a].pos < mutations[b].pos;
              });
}

PYBIND11_MODULE(util, m)
{
    m.doc() = "Miscellaneous utilities for simulations.";

    m.def("change_effect_size",
          [](fwdpy11::Population& pop, const std::size_t index,
             const double new_esize, const double new_dominance) {
              if (index >= pop.mutations.size())
                  {
                      throw std::range_error("mutation index out of range");
                  }
              check_finite(new_esize, "new effect size is not finite");
              check_finite(new_dominance, "new dominance is not finite");
              // Check if we'll need to call fwdpp's internals
              bool need_to_update_storage = false;
              if (pop.mutations[index].neutral && new_esize != 0.0)
                  {
                      need_to_update_storage = true;
                  }
              else if (!pop.mutations[index].neutral && new_esize == 0.0)
                  {
                      need_to_update_storage = true;
                  }
              pop.mutations[index].s = new_esize;
              pop.mutations[index].h = new_dominance;
              // Update the storage of the mutation,
              // which requires a call into fwdpp
              if (need_to_update_storage)
                  {
                      fwdpp::change_neutral(pop, index);
                  }
          },
          py::arg("pop"), py::arg("index"), py::arg("new_esize"),
          py::arg("new_dominance") = 1.0,
          R"delim(
        Change effect sizes and/or dominance of mutations.
        
        From the Python size, the population objects are immutable.
        This function allows you to change the effect of a mutation
        on genetic value.
        
        :param pop: A :class:`fwdpy11.Population`
        :param index: The index of the mutation to change
        :param new_esize: The new value for the `s` field.
        :param new_dominance: (1.0) The new value for the `h` field.

        :versionadded: 0.13.0
        )delim");

    m.def("sort_gamete_keys",
          [](fwdpy11::Population::gcont_t& gametes,
             const fwdpy11::Population::mcont_t& mutations) {
              for (auto& g : gametes)
                  {
                      sort_keys(g.mutations, mutations);
                      sort_keys(g.smutations, mutations);
                  }
          },
          R"delim(
          Sorts mutation keys in gametes.  Useful when constructing
          population objects directly.  See :ref:`popobjects` for example
          use.
          
          .. versionadded:: 0.1.4
          )delim");
}
