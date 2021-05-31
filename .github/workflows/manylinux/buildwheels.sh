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
  
for py in cp37-cp37m cp38-cp38 cp39-cp39
do
    PYPATH=/opt/python/${py}
    PYBIN=${PYPATH}/bin/python 
    rm -rf build/ 
    # Instead of letting setup.py install a newer numpy we install it here
    # using the oldest supported version for ABI compatibility
    ${PYBIN} -m pip install --no-cache oldest-supported-numpy pybind11-global
    PATH=${PYPATH}:$PATH ${PYBIN} setup.py build_ext -i
    PATH=${PYPATH}:$PATH ${PYBIN} setup.py bdist_wheel
done

cd dist
for whl in *.whl; do
    auditwheel repair "$whl"
    rm "$whl"
done
cd ..
