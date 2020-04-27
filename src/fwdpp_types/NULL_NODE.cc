#include <fwdpp/ts/definitions.hpp>
#include <pybind11/pybind11.h>

void init_NULL_NODE(pybind11::module & m)
{
    // TODO: how do I docstring this?
    m.attr("NULL_NODE") = pybind11::int_(fwdpp::ts::TS_NULL_NODE);
}
