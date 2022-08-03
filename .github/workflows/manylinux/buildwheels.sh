#!/bin/bash

set -e -x

yum update -y
yum -y install cmake 

GSL_VERSION=2.5
curl -o gsl-${GSL_VERSION}.tar.gz "ftp://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz" 
tar -zxf gsl-${GSL_VERSION}.tar.gz 
cd gsl-${GSL_VERSION} 
./configure --prefix=/usr/local
make -j 2 
make install 
cd .. 

# Taken from msprime/#2043
# We're running as root in the docker container so git commands issued by
# setuptools_scm will fail without this:
git config --global --add safe.directory /project
# Fetch the full history as we'll be missing tags otherwise.
# git fetch --unshallow
  
for py in cp37-cp37m cp38-cp38 cp39-cp39 cp310-cp310
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
    auditwheel repair "$whl"
    rm "$whl"
done
cd ..
