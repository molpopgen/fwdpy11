#include <cmath>
#include <stdexcept>
#include <fwdpp/sugar/change_neutral.hpp>
#include <fwdpy11/types/Population.hpp>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace
{
    void
    check_finite(const double d, const std::string& error)
    {
        if (!std::isfinite(d))
            {
                throw std::invalid_argument(error);
            }
    }
} // namespace

void
init_change_effect_size(py::module& m)
{
    m.def("change_effect_size",
          [](fwdpy11::Population& pop, const std::size_t index,
             const double new_esize, const double new_dominance,
             py::list new_esizes, py::list new_heffects) {
              if (index >= pop.mutations.size())
                  {
                      throw std::range_error("mutation index out of range");
                  }
              check_finite(new_esize, "new effect size is not finite");
              check_finite(new_dominance, "new dominance is not finite");
              std::vector<double> esizes, heffects;
              for (auto i : new_esizes)
                  {
                      auto d = i.cast<double>();
                      check_finite(d,
                                   "new effect size in new_esizes not finite");
                      esizes.push_back(d);
                  }
              for (auto i : new_heffects)
                  {
                      auto d = i.cast<double>();
                      check_finite(
                          d, "new effect size in new_heffects not finite");
                      heffects.push_back(d);
                  }
              // Check if we'll need to call fwdpp's internals
              bool need_to_update_storage = false;
              if (pop.mutations[index].neutral
                  && (new_esize != 0.0
                      || std::any_of(std::begin(esizes), std::end(esizes),
                                     [](double d) { return d != 0.0; })))
                  {
                      need_to_update_storage = true;
                  }
              else if (!pop.mutations[index].neutral && new_esize == 0.0
                       && std::all_of(std::begin(esizes), std::end(esizes),
                                      [](double d) { return d == 0.0; }))
                  {
                      need_to_update_storage = true;
                  }
              pop.mutations[index].s = new_esize;
              pop.mutations[index].h = new_dominance;
              pop.mutations[index].esizes.swap(esizes);
              pop.mutations[index].heffects.swap(heffects);
              // Update the storage of the mutation,
              // which requires a call into fwdpp
              if (need_to_update_storage)
                  {
                      fwdpp::change_neutral(pop, index);
                  }
          },
          py::arg("pop"), py::arg("index"), py::arg("new_esize") = 0,
          py::arg("new_dominance") = 1.0, py::arg("new_esizes") = py::list{},
          py::arg("new_heffects") = py::list{},
          R"delim(
        Change effect sizes and/or dominance of mutations.
        
        From the Python size, the population objects are immutable.
        This function allows you to change the effect of a mutation
        on genetic value.
        
        :param pop: A :class:`fwdpy11.Population`
        :param index: The index of the mutation to change
        :param new_esize: (0.0) The new value for the `s` field.
        :param new_dominance: (1.0) The new value for the `h` field.
        :param new_esizes: (empty list) New values for :attr:`fwdpy11.Mutation.esizes`
        :param new_heffects: (empty list) New values for :attr:`fwdpy11.Mutation.heffects`

        .. versionadded:: 0.1.3

        .. versionchanged:: 0.2.0

            Modified to act on base class and handle vectors of effect sizes. Default values also updated.
        )delim");
}
