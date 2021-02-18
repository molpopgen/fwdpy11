GSL_VERSION=2.5
curl -o gsl-${GSL_VERSION}.tar.gz "ftp://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz" 
tar -zxf gsl-${GSL_VERSION}.tar.gz 
cd gsl-${GSL_VERSION} 
./configure --prefix=/usr/local
make -j 2 
make install 
cd .. 
  
