#!/bin/bash

curl https://ftp.gnu.org/gnu/gsl/gsl-2.7.tar.gz -o gsl-2.7.tar.gz
tar xzf gsl-2.7.tar.gz
cd gsl-2.7
./configure
make -j 6
sudo make install
