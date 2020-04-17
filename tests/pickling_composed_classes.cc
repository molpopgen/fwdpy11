// Mocking what it takes to get pickling support for a class structure
// and then an object that is composed of unique_ptr to base classes

#include <pybind11/pybind11.h>
#include <string>
#include <memory>

    namespace py = pybind11;

struct Base
{
    virtual ~Base() = default;
    virtual py::object repr() const = 0;
    virtual std::unique_ptr<Base> clone() const = 0;
};

struct Derived : public Base
{
    virtual py::object
    repr() const
    {
        return py::bytes("Derived");
    }

    virtual std::unique_ptr<Base>
    clone() const
    {
        return std::unique_ptr<Base>(new Derived(*this));
    }
};

struct Composed
// This type is a composition of a double
// and an object in the Base class hierarchy.
// This roughly mocks how things like genetic
// value objects are implemented on the C++ side.
{
    double x;
    std::unique_ptr<Base> b;
    template <typename X> Composed(double a, X&& b_) : x(a), b{ b_.clone() } {}
};

PYBIND11_MODULE(pickling_composed_classes, m)
{
    // This is an ABC
    py::class_<Base>(m, "Base");

    // This is a concrete class that can be pickled.
    py::class_<Derived, Base>(m, "Derived")
        .def(py::pickle(
            [](const Derived& o) -> py::object { return py::bytes(o.repr()); },
            [](py::object o) -> Derived {
                auto s = o.cast<std::string>();
                if (s.find("Derived") == std::string::npos)
                    {
                        throw std::runtime_error("invalid obect state");
                    }
                return Derived();
            }))
        .def(py::init<>());

    py::class_<Composed>(m, "Composed")
        .def(py::init(
            [](const double d, const Base& b) { return Composed(d, b); }))
        .def_property_readonly("b",
                               [](const Composed& c) { return c.b->repr(); })
        .def(py::pickle(
            // In order to pickle,
            // we return a tuple of c's data plus
            // the result of pickling c.b
            [](const Composed& c) {
                auto x = py::module::import("pickle");
                auto y = x.attr("dumps")(c.b->clone(), -1);
                return py::make_tuple(c.x, y);
            },
            // Unpickling requires some casting magic
            [](py::tuple t) {
                if (t.size() != 2)
                    {
                        throw std::runtime_error("invalid object state");
                    }
                auto x = py::module::import("pickle");
                auto b = x.attr("loads")(t[1]);
                // We take a const reference to our abstract
                // base in order to compose a new object
                // from the pickled member data.
                // This ONLY works if pickling of classes in
                // our ABC hierarchy is working.
                const Base& br = b.cast<const Base&>();
                return Composed(t[0].cast<double>(), br);
            }));
}
