//
//Copyright(C) 2017 Kevin Thornton < krthornt @uci.edu >
//
//This file is part of fwdpy11.
//
//fwdpy11 is free software : you can redistribute it and / or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//fwdpy11 is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with fwdpy11.If not, see < http: //www.gnu.org/licenses/>.
//

#include <pybind11/pybind11.h>
#include <fwdpy11/types/TemporalSampler.hpp>

namespace py = pybind11;

struct RecordNothing : public fwdpy11::TemporalSampler
{
    std::string serialize() const final{ return std::string(); }
    void deserialize(std::string) final{};
    std::unique_ptr<fwdpy11::TemporalSampler>
    clone() const final
    {
        return std::unique_ptr<fwdpy11::TemporalSampler>(new RecordNothing());
    }
    void
    operator()(const fwdpy11::SlocusPop &) final
    {
    }
    void
    operator()(const fwdpy11::MlocusPop &) final
    {
    }
};

PYBIND11_MODULE(temporal_samplers, m)
{
    m.doc() = "Support for temporal samplers implemented in C++";

    py::class_<fwdpy11::TemporalSampler>(m, "TemporalSampler",
                                         R"delim(
            ABC for temporal samplers.

            This type implements nothing but a name. Derived 
            types are required to implement serialization
            methods on the C++ side in order to implement
            pickling methods on the Python side.
            )delim");

    py::class_<RecordNothing, fwdpy11::TemporalSampler>(m, "RecordNothing")
        .def(py::init<>())
        .def(py::pickle(
             [](const RecordNothing &r) { return py::bytes(r.serialize()); },
             [](py::bytes) { return RecordNothing(); }));
}
