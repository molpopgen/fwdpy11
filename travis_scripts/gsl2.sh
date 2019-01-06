#!/bin/bash

wget http://mirror.us-midwest-1.nexcess.net/gnu/gsl/gsl-2.3.tar.gz
tar xzf gsl-2.3.tar.gz
cd gsl-2.3
./configure
make -j 2 > /dev/null 
sudo make install > /dev/null 2> /dev/null
cd ..
