#!/usr/bin/env bash

    # conda install -c conda-forge sphinx nbsphinx ipython matplotlib msprime
if [ "$USECONDA" == "1" ];
then
    echo "HERE"
    if [ "$TRAVIS_OS_NAME" == "linux" ]; 
    then 
        wget https://repo.continuum.io/miniconda/Miniconda3-4.3.21-Linux-x86_64.sh -O miniconda.sh;
    fi
    if [ "$TRAVIS_OS_NAME" == "osx" ];
    then 
        wget https://repo.continuum.io/miniconda/Miniconda3-4.3.21-MacOSX-x86_64.sh -O miniconda.sh;
    fi
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    conda config --set always_yes yes --set changeps1 no
    #conda update -q conda
    # Useful for debugging any issues with conda
    conda info -a
    if [ "$TRAVIS_OS_NAME" == "linux" ];
    then conda install gcc;
    fi
    if [ "$TRAVIS_OS_NAME" == "osx" -a "$OSXGCC" == "1" ];
    then
        conda install gcc;
    fi
    conda install cython numpy python==3.6 gsl
    conda install -c conda-forge pybind11==2.2.4 numpy==1.16.2 msprime==0.7.0 openblas cmake
    # Weird stuff is happening right now, so we force numpy and numpy base:
    # conda install -c conda-forge numpy==1.16.2 numpy-base==1.16.2
    # conda install -c conda-forge sphinx nbsphinx ipython matplotlib msprime
    pip install cppimport
    echo `which python`
    echo `which python3`
else
    sudo apt-get update -qq
    sudo apt-get -f install python-dev libffi-dev libssl-dev libhdf5-dev cmake
    ./travis_scripts/gsl2.sh
    pip install -r requirements.txt

    # We need to get pybind11 NOT from pip
    curl -L https://github.com/pybind/pybind11/archive/v2.2.4.tar.gz > pybind11-2.2.4.tar.gz
    tar xzf pybind11-2.2.4.tar.gz
    echo "python version is $TRAVIS_PYTHON_VERSION"
    # Skip running the pybind11 unit suite--not needed...
    cd pybind11-2.2.4 && cmake . -DPYBIND11_PYTHON_VERSION=$TRAVIS_PYTHON_VERSION -DPYBIND11_TEST=0 && sudo make install && cd ..
fi
