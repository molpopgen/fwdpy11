#!/bin/bash

wget http://mirror.us-midwest-1.nexcess.net/gnu/gsl/gsl-2.2.tar.gz
tar xzf gsl-2.2.tar.gz
cd gsl-2.2
./configure
make
sudo make install
cd ..
