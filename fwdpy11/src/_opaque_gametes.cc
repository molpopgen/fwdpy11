#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <fwdpp/forward_types.hpp>

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<KTfwd::gamete>);

PYBIND11_MODULE(_opaque_gametes, m)
{
    m.doc() = "Expose C++ vectors of gametes to Python without copies.";

    py::bind_vector<std::vector<KTfwd::gamete>>(
        m, "VecGamete", py::module_local(false),
        "C++ representations of a list of "
        ":class:`fwdpy11.fwdpp_types.Gamete`.  "
        "Typically, access will be read-only.")
        .def(py::pickle(
            [](const std::vector<KTfwd::gamete>& gametes) {
                py::list rv;
                for (auto&& g : gametes)
                    {
                        rv.append(g);
                    }
                return rv;
            },
            [](const py::list l) {
                std::vector<KTfwd::gamete> rv;
                for (auto&& li : l)
                    {
                        rv.push_back(li.cast<KTfwd::gamete>());
                    }
                return rv;
            }));
}
