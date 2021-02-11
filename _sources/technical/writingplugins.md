---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(writingplugins)=

# Writing "plugins" with C++

New functionality may be added through new Python code and/or new C++ code.
Further, you may use the existing C++ types in Python extensions depending on
`fwdpy11`.  For example, you could write a custom "evolve" function for
non-Wright-Fisher models.  Or, you could write custom genetic value objects.
There are several examples of custom genetic value objects in the unit tests.

## Finding the headers

You can find the location of the installed header files using Python:

```{code-cell} python
import fwdpy11

print(fwdpy11.get_includes())
print(fwdpy11.get_fwdpp_includes())
```

The above is useful for generating a functioning `setup.py` file.  Note that you will have
to join the output with the proper compiler flags indicating include paths (typically `-I`).

If using a `Makefile`, it is handy to get the above info via the shell, which is done as follows:

```{code-block} bash

python3 -m fwdpy11 --includes

```

The above command prepends the paths with `-I`.  If that is not desired,
you may get the paths separately for `fwdpp` and `fwdpy11`:

```{code-block} bash

python3 -m fwdpy11 --fwdpp_headers
python3 -m fwdpy11 --fwdpy11_headers

```

## Building with `cmake`

The above two commands are useful when using tools like `cmake` to configure build systems.  Here is an example
from one of the examples that comes with `fwdpy11`:

```{literalinclude} ../../examples/plugin/CMakeLists.txt

```

## Mako headers for cppimport

Extensions using [cppimport][cppimport] require `mako` headers to guide compilation.  You make get a minimal header via the shell:

```{code-block} bash

python3 -m fwdpy11 --mako

```

[cppimport]: https://github.com/tbenthompson/cppimport

## Dealing with `GSL` errors in custom C++ code

The `GSL` uses C-like error handling.  Here, this means that there is a global error handling function
that will print the error to screen and then abort the running process.  The behavior of abort-on-error is not
acceptable here, as the Python session itself will abort!

When `fwdpy11` is imported, the default `GSL` behavior changes.  Instead of aborting, a `RuntimeError` will be raised.
This exception will contain the entire string of text from the `GSL` error message.

However, if you write and C++ code using the `GSL`, you may wish to handle errors locally and maybe return `None` or throw
an exception with your own message in it.  To do so, you must do the following in your C++ code:

1. Temporarily disable `GSL` error handling.
2. Restore the `fwdpy11` error handler at all possible exit paths from your code.

On paper, one could do all of the above using the C API found in the `GSL`.  However, `fwdpy11` provides a class that
provides an idiomatic C++ approach to managing the error handler.  The C++ class `fwdpy11::gsl_error_handler_wrapper`,
defined in `<fwdpy11/gsl/gsl_error_handler_wrapper.hpp>` provide a "smart pointer"-like wrapper around the error
handler. The constructor disables the error handler, storing the value of the disabled handler.  The destructor restores
the handler.

To see this in action, check out the unit test file `tests/test_GSLerror.py` and its associated C++ file
`tests/gsl_error.cc`.


