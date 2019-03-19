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

#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpy11/types/Population.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(fwdpy11::Population::mcont_t);
PYBIND11_MAKE_OPAQUE(fwdpy11::Population::gcont_t);

namespace
{
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
} // namespace

void init_sort_gamete_keys(py::module & m)
{
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
