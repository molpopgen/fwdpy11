#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpp/sugar/generalmut.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<KTfwd::generalmut_vec>);

PYBIND11_MODULE(_opaque_generalmutvecs,m)
{
    py::bind_vector<std::vector<KTfwd::generalmut_vec>>(
        m, "VecGeneralMutVec", py::module_local(false),
        "A list of :class:`fwdpy11.GeneralMutVec`.")
        .def(py::pickle(
            [](const std::vector<KTfwd::generalmut_vec> &mutations) {
                py::list rv;
                for (auto &&i : mutations)
                    {
                        rv.append(i);
                    }
                return rv;
            },
            [](py::list l) {
                std::vector<KTfwd::generalmut_vec> rv;
                for (auto &&i : l)
                    {
                        rv.push_back(i.cast<KTfwd::generalmut_vec>());
                    }
                return rv;
            }));
}
