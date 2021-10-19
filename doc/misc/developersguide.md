(developersguide)=

# Developer's guide

## Working with the source code repository

### Dependencies

#### non-Python

* A C++ compiler that supports C++14 or later.
  On Linux, `gcc` is preferred, although `clang` is very useful to compare against.
  On macOS, use `clang`.
* `cmake`
* The GNU Scientific Library, version 2.3 or later.

:::{todo}
Describe dependencies for the C++ test suite
:::

:::{warning}
`conda` users **must** install all of these dependencies using that system, including the compilers!
See [here](https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html) for info about `conda` compilers for Linux and for macOS.

The `requirements/` folder of the source code repository contains a list of minimal `conda` dependencies that we use for CI.
:::

#### Python

To reproduce the environments that we use for CI and deployment:

```sh
python3 -m pip install -r requirements.txt
# If working on the manual:
python3 -m pip install -r doc/requirements.txt
```

To use the latest version of the dependencies:

```sh
python3 -m pip install -r requirements.in
# If working on the manual:
python3 -m pip install -r doc/requirements.in
```

(sec_requirements_files)=

#### Comments on "requirements" files.

Python uses files with names like `requirements.txt` to specify the package environment that is expected to go along with a certain tool such as this one.
Unfortunately, the contents of such files vary quite a bit across projects.

Here, we take the files `requirements.txt` and `doc/requirements.txt` to specify completely reproducible development and production environments.
Within these files, all dependencies are pinned to specific versions, as are the dependencies of those dependencies.
These files specify the exact versions used for testing the package on Ubuntu environments and for building Docker containers.

The files `requirements.in` and `doc/requirements.in` are looser specifications of these environments.
Here, the dependency pinnings are minimal, trying to keep the build environment compatible with the lowest Python version number that is still supported.
Many users working in extant Python environments may want to use the `.in` versions for installing packages.
Some may also want to remove the installation modifications to various packages in those files.

These pinnings will be updated periodically.

People using `conda` to install `fwdpy11` may consider using the fully-pinned files to install dependencies using `pip`.

#### Updating the CI/deployment dependencies

Use `pip-compile` from [pip-tools](https://pypi.org/project/pip-tools/).

### Building the code

```sh
git submodule init
git submodule update
python3 setup.py build_ext -i
```

### Installing from a clone of the source repository

```{code-block} bash
# from the repository root
pip install .
```

### Installing from the source distribution at PyPi

With all of the dependencies in place, tell `pip` to ignore the binary wheels:

```{code-block} bash
# from the repository root
pip install fwdpy11 --no-binary
```

### Running the Python test suite

```sh
python3 -p ptest tests -n 6
```

Replace `6` with a number of threads appropriate for your machine.

### Building and running the C++ test suite.

:::{todo}
Describe building and executing the C++ test suite.
:::

### Enabling compiler command support for developing the C++ code

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

### C++ code standards

This is a large topic.
Briefly,

* Functions should not be concerned with ownership!
  If an input argument cannot cannot be `nullptr`, then pass it by reference.
  (Or by value if it will be consumed.)
  Likewise, smart pointers are not function arguments unless they will be consumed.
  If the calling environment manages an object with a smart pointer, and `nullptr` is not a valid value to pass to a given function, then the function should be written to take a reference.
* If a pointer type is an argument, `nullptr` must be checked for and handled.
* Always use smart pointers instead of raw allocations.
* Classes should respect the "rule of 5".
  For abstract base classes, this usually means that the copy constructors and assignment operators must be marked `delete` to avoid slicing.
  (In practice, we don't have any code that could lead to slicing, but we design the classes with this in mind.)
* Virtual functions in derived classes are marked `override`.
  Where appropriate, they are also marked `final`.

### Enabling code profiling

By default, fwdpy11 is compiled with aggressive optimizations to help reduce the library size.
One side effect is that it becomes impossible to accurately profile the code.
To override these defaults:

```{code-block} bash

python setup.py build_ext -i --enable-profiling

```

:::{note}

The package should not be installed with profiling enabled.
This method of building is for developers who need to accurately profile the C++ back-end.
Also note that only the main package is affected.
Building the unit test modules is not affected.

:::

### Disabling link-time optimization (LTO)

LTO is enabled by default and reduced the final library size substantially.
However, it takes a long time and is therefore a drag during development.
To disable it:

```{code-block} bash

python setup.py build_ext -i --disable_lto

```

:::{note}

This option only affects the main package and not the unit tests.

:::

### Enabling debugging symbols in the C++ code

```{code-block} bash

python setup.py build_ext -i --debug

```

Debug mode disables all compiler optimizations, allows C-like assertions, and generated debug symbols.

:::{note}

Never install the package compiled in debug mode!
First, things will run much more slowly.
Second, triggering an assertion will cause the Python interpreter to crash.
These assertions exist as a brute-force method to help developers quickly identify bugs.

:::

### Enabling assertions in the C++ code

The fwdpp library code uses C's assert macros in several places.
These are disabled by default.
However, it can be useful to enable them when hacking the code.
To do so, you must manually set your compiler flags with cmake:

```{code-block} bash

cmake . -DCMAKE_CXX_FLAGS="-UNDEBUG -O2 -g"

```

When compiling this way, fwdpy11 makes some extra checks that will throw `RuntimeError` if they fail.
The fwdpp back end also makes extra checks.
If those fail, `abort` will be called, which will crash the Python interpreter.
Thus, compiling with this option is a "serious debugging mode only" option.

### Enabling aggressive debugging of C++ STL templates using GCC

Use the following flags to enable an "extreme" debugging mode of the C++ standard template library:

```{code-block} bash

CXXFLAGS="-D_GLIBCXX_CONCEPT_CHECKS -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" \
   CPPFLAGS="-D_GLIBCXX_CONCEPT_CHECKS -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" python3 setup.py build_ext -i

```

### Static analysis using clang-tidy

It is sometimes useful to go through the code and to a "static" analysis to look for problems.
The clang-tidy tool is especially useful.

For example:

```{code-block} bash
#!/bin/bash

for i in $(find fwdpy11/src/ -name '*.cc')
do
    n=`basename $i .cc`
    echo $i
    clang-tidy $i -checks="cppcoreguidelines-special-member-functions" -header-filter=fwdpy11 -export-fixes=fixes/$n.yaml -- -I/usr/include/python3.8d -I./fwdpy11/headers/fwdpp -I./fwdpy11/headers -I$HOME/.local/include
done
```

### Code formatting style

### Python code standards

#### Code formatting style

The formatting style is [black][black] with a line length of 89, which is specified in the `.flake8` file that is part of the source repository.

The beauty of formatting with `black` is that you can just write all your Python code and then format at the end.
For example, with `vim` or `neovim`, you may use the plugin from the [black][black] repository and the following `normal` mode command to auto-format your code:

```{code-block} vim

:Black

```

::::{note}

At the time of this writing (April 11, 2020), `vim/nvim` integration with `black` requires installing `black` from the GitHub repo and not from PyPi or conda.

```{code-block} bash

pip3 install git+git://github.com/psf/black

```

::::

Placing the following in your `vim` or `neovim` configuration file will map this command to the `F6` button:

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

In the future, we may rely entirely on `black` to sort includes as `black/isort` compatibility evolves.

:::

## Writing documentation

We use [jupyter book](https://www.jupyterbook.org) to write the manual.
This tool uses a flavor of markdown called `MYsT`.
At the moment, there aren't many tools available to do things like auto-format code, etc..

The dependencies for the documentation are described above.

### Building the documentation

```sh
cd doc
make 
```

Then, point your browser to `_build/html/index.html`.

To remove the manual:

```sh
make clean
```

Remember to `cd ..` to get out of the `doc/` folder and get back to the repository root when you're done.

### Docstrings in C++ code

We are moving away from docstrings in C++ code.
Where possible, we prefer to either have a Python class inherit from a `pybind11` class.
Alternately, encapsulation of the `pybind11` class in a Python class is a possibility.
Then, the Python classes contain the docstrings.

For functions, we apply the same idea, defining a Python function that calls the low-level C++ function.

Although this is a bit more work, it provides better support for tools like `mypy`, `jedi`, etc..

[blacken-docs]: https://github.com/asottile/blacken-docs

[black]: https://black.readthedocs.io/en/stable/

[isort]: https://github.com/timothycrosley/isort


