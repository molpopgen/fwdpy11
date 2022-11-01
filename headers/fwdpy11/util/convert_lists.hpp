#ifndef FWDPY11_UTIL_CONVERT_LISTS_HPP
#define FWDPY11_UTIL_CONVERT_LISTS_HPP

#include <pybind11/pybind11.h>

namespace fwdpy11
{
    template <typename T>
    inline pybind11::list
    vector_to_list(const T& t)
    {
        pybind11::list rv;
        for (auto& i : t)
            {
                rv.append(i);
            }
        return rv;
    }

    template <typename T>
    inline T
    list_to_vector(pybind11::list l)
    {
        T rv;
        rv.reserve(l.size());
        for (auto& i : l)
            {
                rv.push_back(i.cast<typename T::value_type>());
            }
        return rv;
    }
} // namespace fwdpy11

#endif
