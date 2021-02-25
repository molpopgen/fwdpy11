(developersguide)=

# Developer's guide

## Enabling compiler command support for developing the C++ code

Many language servers support processing a local file called `compile_commands.json`.
See [here](https://clang.llvm.org/docs/JSONCompilationDatabase.html) for background.

To generate this file, execute the following steps from the root of the source code repository:

```sh
mkdir ccomands
cd ccomands
cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=1
mv compile_commands.json ..
cd ..
rm -rf ccomands
```

Now, language servers supporting `clangd` will have nice error checking and code completion for the C++ code!

## C++ code standards

## Enabling code profiling

By default, fwdpy11 is compiled with aggressive optimizations to help reduce the library size. One side effect
is that it becomes impossible to accurately profile the code.  To override these defaults:

```{code-block} bash

python setup.py build_ext -i --enable-profiling

```

:::{note}

The package should not be installed with profiling enabled. This method of building
is for developers who need to accurately profile the C++ back-end.  Also note that
only the main package is affected.  Building the unit test modules is not affected.

:::

## Disabling link-time optimization (LTO)

LTO is enabled by default and reduced the final library size substantially. However, it takes a
long time and is therefore a drag during development.  To disable it:

```{code-block} bash

python setup.py build_ext -i --disable_lto

```

:::{note}

This option only affects the main package and not the unit tests.

:::

## Enabling debugging symbols in the C++ code

```{code-block} bash

python setup.py build_ext -i --debug

```

Debug mode disables all compiler optimizations, allows C-like assertions, and generated debug symbols.

:::{note}

Never install the package compiled in debug mode!  First, things will run much more slowly.
Second, triggering an assertion will cause the Python interpreter to crash.  These assertions
exist as a brute-force method to help developers quickly identify bugs.

:::

## Enabling assertions in the C++ code

The fwdpp library code uses C's assert macros in several places.  These are disabled by default.  However, it can be useful to
enable them when hacking the code.  To do so, you must manually set your compiler flags with cmake:

```{code-block} bash

cmake . -DCMAKE_CXX_FLAGS="-UNDEBUG -O2 -g"

```

When compiling this way, fwdpy11 makes some extra checks that will throw `RuntimeError` if they fail.  The fwdpp back
end also makes extra checks.  If those fail, `abort` will be called, which will crash the Python interpreter.  Thus,
compiling with this option is a "serious debugging mode only" option.

## Enabling aggressive debugging of C++ STL templates using GCC

Use the following flags to enable an "extreme" debugging mode of the C++ standard template library:

```{code-block} bash

CXXFLAGS="-D_GLIBCXX_CONCEPT_CHECKS -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" \
   CPPFLAGS="-D_GLIBCXX_CONCEPT_CHECKS -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" python3 setup.py build_ext -i

```

## Static analysis using clang-tidy

It is sometimes useful to go through the code and to a "static" analysis to look for problems. The clang-tidy
tool is especially useful.  For example:

```{code-block} bash

find fwdpy11/src/ -name '*.cc' | xargs -I  clang-tidy -checks='performance-*'   -- -I/usr/include/python3.7m -I./fwdpy11/headers/fwdpp -I./fwdpy11/headers

```

### Code formatting style

## Python code standards

### Code formatting style

The formatting style is [black][black] with a line length of 89, which is specified
in the `.flake8` file that is part of the source repository.

The beauty of formatting with `black` is that you can just write all your Python
code and then format at the end. For example, with `vim` or `neovim`, you may use the
plugin from the [black][black] repository and the following `normal` mode command to
auto-format your code:

```{code-block} vim

:Black

```

::::{note}

At the time of this writing (April 11, 2020), `vim/nvim` integration
with `black` requires installing `black` from the GitHub repo
and not from PyPi or conda.

```{code-block} bash

pip3 install git+git://github.com/psf/black

```

::::

Placing the following in your `vim` or `neovim` configuration file will
map this command to the `F6` button:

```{code-block} vim

autocmd FileType python nnoremap <buffer> <f6> :Black<CR>

```

The [black][black] repository describes methods for integration with other editors.

:::{note}

Other than line length, all `black` parameters used are the defaults!

:::

In addition to `black`, `import` statements should be sorted using [isort][isort]:

```{code-block} bash

isort file.py

```

Or, in `vim/neovim`:

```{code-block} vim

:!isort %

```

:::{note}

In the future, we may rely entirely on `black` to sort includes as
`black/isort` compatibility evolves.

:::

## Writing documentation

### Formatting the Python code in rst files

Code blocks in documentation should be formatted using [blacken-docs][blacken-docs].  In general,
you can write your code blocks and format them after the fact using the following
command:

```{code-block} bash

blacken-docs -E -l 89 -t py36 file.rst

```

Using an editor like `vim` or `neovim`, you can format from within the editor using
the following `normal` mode command:

```{code-block} vim

:!blacken-docs -E -l 89 -t py36 %

```

### Docstrings in C++ code

[blacken-docs]: https://github.com/asottile/blacken-docs

[black]: https://black.readthedocs.io/en/stable/

[isort]: https://github.com/timothycrosley/isort


