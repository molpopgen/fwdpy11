#ifndef FWDPY11_TYPES_DIPLOID_HPP__
#define FWDPY11_TYPES_DIPLOID_HPP__

#include <cstdint>
#include <cstddef>
#include <vector>
#include <pybind11/pybind11.h>
#include <fwdpp/io/scalar_serialization.hpp>
#include <fwdpp/io/diploid.hpp>

namespace fwdpy11
{
    struct Diploid
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
        std::uint32_t deme;
        std::int32_t sex;
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
        Diploid() noexcept
            : first(first_type()), second(second_type()), label(0), deme(0),
              sex(-1), g(0.), e(0.), w(1.), parental_data(pybind11::none())
        {
        }
        //! Construct from two indexes to gametes
        Diploid(first_type g1, first_type g2) noexcept
            : first(g1), second(g2), label(0), deme(0), sex(-1), g(0.), e(0.),
              w(1.), parental_data(pybind11::none())
        {
        }

        Diploid(first_type g1, first_type g2, std::size_t label_,
                  double deme_, double sex_, double g_, double e_, double w_)
            : first(g1), second(g2), label(label_), deme(deme_), sex(sex_),
              g(g_), e(e_), w(w_), parental_data(pybind11::none())
        {
        }

        static inline Diploid
        create(first_type g1, first_type g2, std::size_t label_, double deme_,
               double sex_, double g_, double e_, double w_)
        {
            return Diploid(g1, g2, label_, deme_, sex_, g_, e_, w_);
        }

        inline bool
        operator==(const Diploid& dip) const noexcept
        //! Required for py::bind_vector
        {
            auto cpp_data_comparison
                = this->first == dip.first && this->second == dip.second
                  && this->w == dip.w && this->g == dip.g && this->e == dip.e
                  && this->label == dip.label && this->deme == dip.deme
                  && this->sex == dip.sex;
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
    using dipvector_t = std::vector<Diploid>;
}

namespace fwdpp
{
    namespace io
    {
        template <> struct serialize_diploid<fwdpy11::Diploid>
        {
            template <typename streamtype>
            inline void
            operator()(streamtype& buffer, const fwdpy11::Diploid& dip) const
            {
                fwdpp::io::scalar_writer w;
                w(buffer, &dip.first);
                w(buffer, &dip.second);
                w(buffer, &dip.g);
                w(buffer, &dip.e);
                w(buffer, &dip.w);
                w(buffer, &dip.label);
                w(buffer, &dip.deme);
                w(buffer, &dip.sex);
            }
        };

        template <> struct deserialize_diploid<fwdpy11::Diploid>
        {
            template <typename streamtype>
            inline void
            operator()(streamtype& buffer, fwdpy11::Diploid& dip) const
            {
                fwdpp::io::scalar_reader r;
                r(buffer, &dip.first);
                r(buffer, &dip.second);
                r(buffer, &dip.g);
                r(buffer, &dip.e);
                r(buffer, &dip.w);
                r(buffer, &dip.label);
                r(buffer, &dip.deme);
                r(buffer, &dip.sex);
            }
        };
    }
}

#endif
