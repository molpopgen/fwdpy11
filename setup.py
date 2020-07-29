import glob
import os
import platform
import re
import subprocess
import sys
from distutils.version import LooseVersion

import pybind11
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

PYBIND11_MIN_VERSION = "2.4.3"
TSKIT_MIN_VERSION = "0.2.3"

if sys.version_info[0] < 3:
    raise ValueError("Python 3 is required!")

if pybind11.__version__ < PYBIND11_MIN_VERSION:
    raise RuntimeError(f"pybind11 >= {PYBIND11_MIN_VERSION} required")


if "--weffcpp" in sys.argv:
    USE_WEFFCPP = True
    sys.argv.remove("--weffcpp")
else:
    USE_WEFFCPP = False

if "--debug" in sys.argv:
    DEBUG_MODE = True
    sys.argv.remove("--debug")
else:
    DEBUG_MODE = False

if "--enable-profiling" in sys.argv:
    ENABLE_PROFILING = True
    sys.argv.remove("--enable-profiling")
else:
    ENABLE_PROFILING = False

if "--skip_tests" in sys.argv:
    SKIP_BUILDING_TESTS = True
    sys.argv.remove("--skip_tests")
else:
    SKIP_BUILDING_TESTS = False

if "--disable_lto" in sys.argv:
    SKIP_LTO = True
    sys.argv.remove("--disable_lto")
else:
    SKIP_LTO = False

if "--cpp17" in sys.argv:
    USE_CPP17 = True
    sys.argv.remove("--cpp17")
else:
    USE_CPP17 = False

# if '--assert' in sys.argv:
#     ASSERT_MODE = True
#     sys.argv.remove('--assert')
# else:
#     ASSERT_MODE = False


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r"version\s*([\d.]+)", out.decode()).group(1)
            )
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]

        cfg = "Debug" if DEBUG_MODE is True else "Release"

        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
            ]
            if sys.maxsize > 2 ** 32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j2"]

        if ENABLE_PROFILING is True:
            cmake_args += ["-DENABLE_PROFILING=ON"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )

        if USE_WEFFCPP is True:
            cmake_args.append("-DUSE_WEFFCPP=ON")
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        if SKIP_BUILDING_TESTS is True:
            cmake_args.append("-DBUILD_UNIT_TESTS=OFF")

        if SKIP_LTO is True:
            cmake_args.append("-DDISABLE_LTO=ON")

        if USE_CPP17 is True:
            cmake_args.append("-DUSECPP17=ON")

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11

        return pybind11.get_include(self.user)


PKGS = [
    "fwdpy11",
    "fwdpy11.demographic_models",
    "fwdpy11.tskit_tools",
    "fwdpy11._types",
    "fwdpy11._monkeypatch",
]

INCLUDES = [
    "fwdpy11/headers",
    "fwdpy11/headers/fwdpp",
    # Path to pybind11 headers
    get_pybind_include(),
    get_pybind_include(user=True),
    os.path.join(sys.prefix, "include"),
]

LIBRARY_DIRS = [os.path.join(sys.prefix, "lib")]

ext_modules = [
    CMakeExtension("fwdpy11._fwdpy11"),
]


# Figure out the headers we need to install:
generated_package_data = {}
for root, dirnames, filenames in os.walk("fwdpy11/headers"):
    if (
        "testsuite" not in root
        and "examples" not in root
        and "python_examples" not in root
    ):
        g = glob.glob(root + "/*.hh")
        if len(g) > 0:
            replace = root.replace("/", ".")
            # If there's a header file, we add the directory as a package
            if replace not in PKGS:
                PKGS.append(replace)
            generated_package_data[replace] = ["*.hh"]
        g = glob.glob(root + "/*.hpp")
        if len(g) > 0:
            replace = root.replace("/", ".")
            # If there's a header file, we add the directory as a package
            if replace not in PKGS:
                PKGS.append(replace)
            try:
                if "*.hpp" not in generated_package_data[replace]:
                    generated_package_data[replace].append("*.hpp")
            except:  # NOQA
                generated_package_data[replace] = ["*.hpp"]
        g = glob.glob(root + "/*.tcc")
        if len(g) > 0:
            replace = root.replace("/", ".")
            # If there's a template implementation file,
            # we add the directory as a package
            if replace not in PKGS:
                PKGS.append(replace)
            try:
                if "*.tcc" not in generated_package_data[replace]:
                    generated_package_data[replace].append("*.tcc")
            except:  # NOQA
                generated_package_data[replace] = ["*.tcc"]

long_desc = open("README.rst").read()

setup(
    name="fwdpy11",
    author="Kevin Thornton",
    author_email="krthornt@uci.edu",
    url="http://molpopgen.github.io/fwdpy11",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],
    description="Forward-time population genetic simulation in Python",
    license="GNU GPLv3+",
    provides=["fwdpy11"],
    obsoletes=["none"],
    data_files=[("fwdpy11", ["COPYING", "README.rst"])],
    long_description=long_desc,
    long_description_content_type="text/x-rst",
    ext_modules=ext_modules,
    python_requires=">=3.6",
    install_requires=[
        f"pybind11>={PYBIND11_MIN_VERSION}",
        "numpy",
        f"tskit>={TSKIT_MIN_VERSION}",
        "sparse",
        "attrs",
    ],
    setup_requires=["setuptools_scm"],
    use_scm_version={"write_to": "fwdpy11/_version.py"},
    cmdclass={"build_ext": CMakeBuild},
    packages=PKGS,
    package_data=generated_package_data,
    zip_safe=False,
)
