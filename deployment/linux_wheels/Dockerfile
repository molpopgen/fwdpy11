# This file builds an image that does a complete
# test of building and installing wheels and
# then running the test suite for multiple 
# Python versions
FROM quay.io/pypa/manylinux_2_28_x86_64:2022-10-02-69a0972

WORKDIR /app

COPY . /app

RUN yum update -y && yum install -y cmake curl gsl-devel \
  && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y \
  && source "$HOME/.cargo/env" \
  # Pin the rustc toolchain to a specific version.
  # Rust 1.64.0 will change the minimum glibc ABI
  # to a version incompatibly with manylinux_2014,
  # so we need to be careful in general.
  && rustup override set 1.62.1 \
  # Pin cbindgen 
  && cargo install cbindgen@0.24.3 \
  && rustc --version \
  && cbindgen --version \
  # The GSL version available from yum install is too old so we manually install.
  # && bash deployment/linux_wheels/install_gsl.sh \
  && bash deployment/linux_wheels/build_and_audit.sh \
  && rm -rf build

RUN bash deployment/linux_wheels/install_wheels_run_tests.sh

