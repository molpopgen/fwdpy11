#!/bin/bash

set -e -x

yum update -y
yum -y install cmake curl gsl-devel

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y 
source "$HOME/.cargo/env" 
# Pin the rustc toolchain to a specific version.
# Rust 1.64.0 will change the minimum glibc ABI
# to a version incompatibly with manylinux_2014,
# so we need to be careful in general.
rustup override set 1.62.1 
# Pin cbindgen 
cargo install cbindgen@0.24.3 

# GSL_VERSION=2.5
# curl -o gsl-${GSL_VERSION}.tar.gz "ftp://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz" 
# tar -zxf gsl-${GSL_VERSION}.tar.gz 
# cd gsl-${GSL_VERSION} 
# ./configure --prefix=/usr/local
# make -j 2 
# make install 
# cd .. 

# Taken from msprime/#2043
# We're running as root in the docker container so git commands issued by
# setuptools_scm will fail without this:
git config --global --add safe.directory /project
# Fetch the full history as we'll be missing tags otherwise.
# git fetch --unshallow
  
for py in cp38-cp38 cp39-cp39 cp310-cp310 cp311-cp311
do
    PYPATH=/opt/python/${py}
    PYBIN=${PYPATH}/bin/python 
    rm -rf build/ 
    # Instead of letting setup.py install a newer numpy we install it here
    # using the oldest supported version for ABI compatibility
    ${PYBIN} -m pip install --no-cache oldest-supported-numpy build
    PATH=${PYPATH}:$PATH ${PYBIN} -m build
done

cd dist
for whl in *.whl; do
    # Need to set this so that the core library
    # can be found
    LD_LIBRARY_PATH=fwdpy11 auditwheel repair "$whl"
    rm "$whl"
done
cd ..
