.. _developersguide:

Developer's guide
====================================================================================

C++ code standards
--------------------------------------------

Code formatting style
+++++++++++++++++++++++++++++++++++++++++

Python code standards
--------------------------------------------

Code formatting style
+++++++++++++++++++++++++++++++++++++++++

The formatting style is black_ with a line length of 89, which is specified
in the `.flake8` file that is part of the source repository.

The beauty of formatting with `black` is that you can just write all your Python
code and then format at the end. For example, with `vim` or `neovim`, you may use the
plugin from the black_ repository and the following `normal` mode command to
auto-format your code:

.. code-block:: vim

    :Black

.. note::

    At the time of this writing (April 11, 2020), `vim/nvim` integration
    with `black` requires installing `black` from the GitHub repo
    and not from PyPi or conda.

    .. code-block:: bash

       pip3 install git+git://github.com/psf/black

Placing the following in your `vim` or `neovim` configuration file will
map this command to the `F6` button:

.. code-block:: vim

    autocmd FileType python nnoremap <buffer> <f6> :Black<CR>

The black_ repository describes methods for integration with other editors.

.. note::

   Other than line length, all `black` parameters used are the defaults!

In addition to `black`, `import` statements should be sorted using isort_:

.. code-block:: bash

    isort file.py

Or, in `vim/neovim`:

.. code-block:: vim

    :!isort %

.. note::

    In the future, we may rely entirely on `black` to sort includes as
    `black/isort` compatibility evolves.


Writing documentation
--------------------------------------------

Formatting the Python code in rst files
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Code blocks in documentation should be formatted using blacken-docs_.  In general,
you can write your code blocks and format them after the fact using the following
command:

.. code-block:: bash

    blacken-docs -E -l 89 -t py36 file.rst

Using an editor like `vim` or `neovim`, you can format from within the editor using
the following `normal` mode command:

.. code-block:: vim

    :!blacken-docs -E -l 89 -t py36 %

Code blocks written using the `iPython` directive may fail to execute after formatting.
Usually, this is due to blank likes being added in order to generate `PEP8`-compliant
Python code.  However, this procedure often generates too many blank lines and the `iPython`
parser raises an exception.  Sigh.  To work around this:

* Add more `.. ipython:: python` directives to split the blocks up
* Wrap class definitions in `# fmt: off` and `# fmt: on` comments to disable formatting using `black`.
  You can find examples by grepping for `fmt` within the `.rst` files.

The latter recommendation is needed because `PEP8` wants blank lines between the definitions of 
class functions, yet the `iPython` parser will fail to properly parse such a class.  I tend to
run `blacken-docs` after defining a class, manually delete the blank lines, and then wrap in the `fmt`
comments.  This procedure is manual but not too burdensome.

Docstrings in C++ code
++++++++++++++++++++++++++++++++++++++++++++

.. _blacken-docs: https://github.com/asottile/blacken-docs
.. _black: https://black.readthedocs.io/en/stable/
.. _isort: https://github.com/timothycrosley/isort
