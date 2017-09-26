#!/usr/bin/env bash

if [ "$USECONDA" == "1" ];
then
    if [ "$TRAVIS_OS_NAME" == "linux" ]; 
    then 
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
    fi
    if [ "$TRAVIS_OS_NAME" == "osx" ];
    then 
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
    fi
    bash miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    conda config --set always_yes yes --set changeps1 no
    conda update -q conda
    Useful for debugging any issues with conda
    conda info -a
    if [ "$TRAVIS_OS_NAME" == "linux" ];
    then conda install gcc;
    fi
    if [ "$TRAVIS_OS_NAME" == "osx" -a "$OSXGCC" == "1" ];
    then
        conda install gcc;
    fi
    conda install cython numpy python==3.5 gsl
    conda install -c conda-forge pybind11==2.2.0
    conda install -c conda-forge sphinx nbsphinx ipython matplotlib
    pip install cppimport
else
    sudo apt-get update -qq
    sudo apt-get -f install python-dev libffi-dev libssl-dev
    pip install -r requirements.txt
fi
