PYTHON=$1

# SETUP
set -e -x
rm -rf venv

# INSTALL SYSTEM DEPENDENCIES

yum -y install curl gsl-devel

# INSTALL RUST

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y 
source "$HOME/.cargo/env" 
# Pin the rustc toolchain to a specific version.
# Rust 1.64.0 will change the minimum glibc ABI
# to a version incompatible with manylinux_2014,
# so we need to be careful in general.
rustup override set 1.62.1 
# Pin cbindgen 
cargo install --locked cbindgen@0.24.3 

# Taken from msprime/#2043
# We're running as root in the docker container so git commands issued by
# setuptools_scm will fail without this:
git config --global --add safe.directory /project

# BUILD WHEEL INSIDE A VENV

$(which $PYTHON) -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip setuptools build oldest-supported-numpy
python -m build .
deactivate
rm -rf build venv

cd dist
for whl in *.whl; do
    # Need to set this so that the core library
    # can be found
    LD_LIBRARY_PATH=fwdpy11 auditwheel repair "$whl"
    rm "$whl"
done
cd ..

# INSTALL WHEEL INTO VENV FOR TESTING

$(which $PYTHON) -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip
python -m pip install wheel
python -m pip install fwdpy11 --pre --no-cache-dir --only-binary fwdpy11 --find-links dist/wheelhouse
# cd to avoid having the fwdpy11/ directory get mistaken for a package
TESTDIR=wheel_tests_$PYTHON
if [ -e  $TESTDIR ]
then
    rm -rf $TESTDIR
fi
mkdir $TESTDIR
cd $TESTDIR
python -c "import fwdpy11;print(fwdpy11.__version__)"
python -m pip install pytest pytest-xdist hypothesis msprime
# Copy the tests b/c some of what they do depends on paths
# and will fail don't run them from right outside that directory.
cp -r ../tests .
python -m pytest tests -n 4
rm -rf $TESTDIR
cd ..
deactivate
rm -rf venv
