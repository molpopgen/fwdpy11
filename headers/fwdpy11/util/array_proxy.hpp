//
// Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
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

#ifndef FWPY11_UTIL_ARRAY_PROXY_HPP
#define FWPY11_UTIL_ARRAY_PROXY_HPP

#include <cstdint>
#include <vector>
#include <pybind11/pybind11.h>

namespace fwdpy11
{
    template <typename T> struct array_proxy
    // NOTE: not the most const-correct type
    {
        T* data;
        std::size_t size;
        array_proxy() : data{nullptr}, size{0}
        {
        }

        template <typename VT>
        explicit array_proxy(const VT& v) : data{v.data()}, size{v.size()}
        {
        }

        template <typename VT>
        explicit array_proxy(VT& v) : data{v.data()}, size{v.size()}
        {
        }

        void
        set(const std::vector<T>& v)
        {
            data = const_cast<T*>(v.data());
            size = v.size();
        }

        void
        set(std::vector<T>& v)
        {
            data = v.data();
            size = v.size();
        }
    };

    template <typename T>
    inline pybind11::buffer_info
    as_buffer(const array_proxy<T>& self)
    {
        return pybind11::buffer_info(self.data, sizeof(T),
                                     pybind11::format_descriptor<T>::format(), 1,
                                     {self.size}, {sizeof(T)});
    }

    // fwdpy11._Uint32ArrayProxy
    using uint32_array_proxy = array_proxy<std::uint32_t>;
    // fwdpy11._FloatArrayProxy
    using double_array_proxy = array_proxy<double>;
}

#endif
