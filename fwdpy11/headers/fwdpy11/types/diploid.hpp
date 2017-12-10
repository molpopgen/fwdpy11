#ifndef FWDPY11_TYPES_DIPLOID_HPP__
#define FWDPY11_TYPES_DIPLOID_HPP__

#include <cstddef>
#include <vector>
#include <pybind11/pybind11.h>

namespace fwdpy11
{
    struct diploid_t
    /*!
      \brief Custom diploid type.
    */
    {
        using first_type = std::size_t;
        using second_type = std::size_t;
        //! First gamete.  A gamete is vector<size_t> where the elements are
        //! indexes to a population's gamete container
        first_type first;
        //! Second gamete. A gamete is vector<size_t> where the elements are
        //! indexes to a population's gamete container
        second_type second;
        //! 64 bits of data to do stuff with.  Initialized to zero upon
        //! construction
        std::size_t label;
        //! Genetic component of trait value.  This is not necessarily written
        //! to by a simulation.
        double g;
        //! Random component of trait value.  This is not necessarily written
        //! to by a simulation.
        double e;
        //! Fitness.  This is not necessarily written to by a simulation.
        double w;
        //! IDs of parents.  NB: this will be changed in future releases
        pybind11::object parental_data;
        //! Constructor
        diploid_t() noexcept
            : first(first_type()), second(second_type()), label(0), g(0.),
              e(0.), w(1.), parental_data(pybind11::none())
        {
        }
        //! Construct from two indexes to gametes
        diploid_t(first_type g1, first_type g2) noexcept
            : first(g1), second(g2), label(0), g(0.), e(0.), w(1.),
              parental_data(pybind11::none())
        {
        }

        diploid_t(first_type g1, first_type g2, std::size_t label_, double g_,
                  double e_, double w_)
            : first(g1), second(g2), label(label_), g(g_), e(e_), w(w_),
              parental_data(pybind11::none())
        {
        }

        static inline diploid_t
        create(first_type g1, first_type g2, std::size_t label_, double g_,
               double e_, double w_)
        {
            return diploid_t(g1, g2, label_, g_, e_, w_);
        }

        inline bool
        operator==(const diploid_t& dip) const noexcept
        //! Required for py::bind_vector
        {
            auto cpp_data_comparison
                = this->first == dip.first && this->second == dip.second
                  && this->w == dip.w && this->g == dip.g && this->e == dip.e
                  && this->label == dip.label;
            // We now attempt to compare the parental data.
            // pybind11 doesn't have a rich comparison support,
            // so we rely on the __eq__ attribute.  If no
            // such attribute exists, we simply clear the
            // error and move on.
            try
                {
                    bool parental_data_comp
                        = this->parental_data.attr("__eq__")(dip.parental_data)
                              .cast<bool>();
                    return cpp_data_comparison && parental_data_comp;
                }
            catch (...)
                {
                    PyErr_Clear();
                }
            return cpp_data_comparison;
        }
    };
    
    //! Typedef for container of diploids
    using dipvector_t = std::vector<diploid_t>;
}

#endif
