.. _developers:

Information for developers
====================================================================================

New functionality may be added through new Python code and/or new C++ code.  Further, you may use the existing C++ types
in Python extensions depending on fwdpy11.  For example, you could write a custom "evolve" function for
non-Wright-Fisher models.  Or, you could write custom fitness functions.  See :ref:`customgvalues` for an example of the
latter, which uses cppimport_ to create an "on the fly" module that is compiled upon import.

Finding the headers
---------------------------------------

You can find the location of the installed header files programatically within Python:

.. ipython:: python

    import fwdpy11
    print(fwdpy11.get_includes())
    print(fwdpy11.get_fwdpp_includes())

The above is useful for generating a functioning `setup.py` file.  Note that you will have to join the output with the
proper compiler flags indicating include paths (typically `-I`).

If using a `Makefile`, it is handy to get the above info via the shell, which is done as follows:

.. code-block:: bash

    python3 -m fwdpy11 --includes

Mako headers for cppimport
------------------------------------------

Extensions using cppimport_ require "mako" headers to guide compilation.  You make get a minimal header via the shell:

.. code-block:: bash

    python3 -m fwdpy11 --mako

.. _cppimport: https://github.com/tbenthompson/cppimport

