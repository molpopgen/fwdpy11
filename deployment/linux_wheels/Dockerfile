# This file builds an image that does a complete
# test of building and installing wheels and
# then running the test suite for multiple 
# Python versions
FROM quay.io/pypa/manylinux2014_x86_64

WORKDIR /app

COPY . /app

RUN yum update -y && yum install -y cmake \
  # The GSL version available from yum install is too old so we manually install.
  && bash deployment/linux_wheels/install_gsl.sh \
  && bash deployment/linux_wheels/build_and_audit.sh \
  && rm -rf build

RUN bash deployment/linux_wheels/install_wheels_run_tests.sh

