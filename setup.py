from setuptools import setup
from setuptools import Extension
from setuptools.command.build_ext import build_ext
import sys
import pybind11
import setuptools
import os
import platform
import glob
import subprocess
from distutils.version import LooseVersion


if sys.version_info[0] < 3:
    raise ValueError("Python 3 is required!")

__version__ = '0.2.0a0'

if sys.version_info < (3, 3):
    raise RuntimeError("Python >= 3.3 required")

if pybind11.__version__ < '2.2.3':
    raise RuntimeError("pybind11 >= " + '2.2.3' + " required")

if sys.version_info >= (3, 7):
    if pybind11.__version__ < '2.3.0':
        raise RuntimeError(
            "Python 3.7 and newer required pybind11 2.3 or greater")


# clang/llvm is default for OS X builds.
# can over-ride darwin-specific options
# with setup.py --gcc install
if '--gcc' in sys.argv:
    USE_GCC = True
    sys.argv.remove('--gcc')
else:
    USE_GCC = False

if '--debug' in sys.argv:
    DEBUG_MODE = True
else:
    DEBUG_MODE = False

if '--assert' in sys.argv:
    ASSERT_MODE = True
    sys.argv.remove('--assert')
else:
    ASSERT_MODE = False


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(
            self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += [
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j4']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] +
                              cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] +
                              build_args, cwd=self.build_temp)


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


PKGS = ['fwdpy11']

INCLUDES = [
    'fwdpy11/headers',
    'fwdpy11/headers/fwdpp',
    # Path to pybind11 headers
    get_pybind_include(),
    get_pybind_include(user=True),
    os.path.join(sys.prefix, 'include')
]

LIBRARY_DIRS = [
    os.path.join(sys.prefix, 'lib')
]

ext_modules = [
    CMakeExtension('fwdpy11._init'),
    CMakeExtension('fwdpy11.fwdpp_types'),
    CMakeExtension('fwdpy11._opaque_gametes'),
    CMakeExtension('fwdpy11._opaque_mutations'),
    CMakeExtension('fwdpy11._opaque_diploids'),
    CMakeExtension('fwdpy11.fwdpp_extensions'),
    CMakeExtension('fwdpy11._Population'),
    CMakeExtension('fwdpy11._Populations'),
    CMakeExtension('fwdpy11.fwdpy11_types'),
    CMakeExtension('fwdpy11.genetic_value_noise',),
    CMakeExtension('fwdpy11.wright_fisher_slocus',),
    CMakeExtension('fwdpy11.wright_fisher_mlocus',),
    CMakeExtension('fwdpy11.gsl_random',),
    CMakeExtension('fwdpy11.multilocus'),
    CMakeExtension('fwdpy11.util'),
    CMakeExtension('fwdpy11.ts',),
    CMakeExtension('fwdpy11.tsrecorders',),
    CMakeExtension('fwdpy11.ts_from_msprime',),
    CMakeExtension('fwdpy11.wright_fisher_slocus_ts',),
]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin' and USE_GCC is False:
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' %
                        self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden') and \
                       (sys.platform != 'darwin' or USE_GCC is True):
                opts.append('-fvisibility=hidden')
            if has_flag(self.compiler, '-g0') and DEBUG_MODE is False:
                opts.append('-g0')
            if ASSERT_MODE is True:
                opts.append('-UNDEBUG')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' %
                        self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
            if sys.platform == 'darwin' and USE_GCC is False:
                ext.extra_link_args = ['-stdlib=libc++',
                                       '-mmacosx-version-min=10.7']
        build_ext.build_extensions(self)


# Figure out the headers we need to install:
generated_package_data = {}
for root, dirnames, filenames in os.walk('fwdpy11/headers'):
    if 'testsuite' not in root and 'examples' not in root and \
       'python_examples' not in root:
        g = glob.glob(root + '/*.hh')
        if len(g) > 0:
            replace = root.replace('/', '.')
            # If there's a header file, we add the directory as a package
            if replace not in PKGS:
                PKGS.append(replace)
            generated_package_data[replace] = ['*.hh']
        g = glob.glob(root + '/*.hpp')
        if len(g) > 0:
            replace = root.replace('/',  '.')
            # If there's a header file, we add the directory as a package
            if replace not in PKGS:
                PKGS.append(replace)
            try:
                if '*.hpp' not in generated_package_data[replace]:
                    generated_package_data[replace].append('*.hpp')
            except:
                generated_package_data[replace] = ['*.hpp']
        g = glob.glob(root + '/*.tcc')
        if len(g) > 0:
            replace = root.replace('/', '.')
            # If there's a template implementation file,
            # we add the directory as a package
            if replace not in PKGS:
                PKGS.append(replace)
            try:
                if '*.tcc' not in generated_package_data[replace]:
                    generated_package_data[replace].append('*.tcc')
            except:
                generated_package_data[replace] = ['*.tcc']

long_desc = open("README.rst").read()

setup(
    name='fwdpy11',
    version=__version__,
    author='Kevin Thornton',
    author_email='krthornt@uci.edu',
    url='http://molpopgen.github.io/fwdpy11',
    classifiers=['Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'],
    description='Forward-time population genetic simulation in Python',
    license='GPL >= 3',
    requires=['pybind11', 'numpy'],
    provides=['fwdpy11'],
    obsoletes=['none'],
    data_files=[('fwdpy11', ['COPYING', 'README.rst'])],
    long_description=long_desc,
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.2.3', 'numpy'],
    cmdclass={'build_ext': CMakeBuild},
    packages=PKGS,
    package_data=generated_package_data,
    zip_safe=False,
)
